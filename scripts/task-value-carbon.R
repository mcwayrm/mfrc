required_packages <- c("terra", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
script_argument <- grep("^--file=", commandArgs(), value = TRUE)
if (length(script_argument) != 1L) {
  stop("Run with Rscript scripts/task-value-carbon.R.")
}
root_dir <- dirname(dirname(normalizePath(
  sub("^--file=", "", script_argument), mustWork = TRUE
)))
input_dir <- file.path(root_dir, "data", "forest-layer")
output_dir <- file.path(root_dir, "data", "GEP")
years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
prices <- c(21.40, 21.40, 36, 42, 46)
epa_table <- paste0(
  "https://19january2017snapshot.epa.gov/climatechange/",
  "social-cost-carbon_.html"
)
early_source <- "https://www.epa.gov/sites/default/files/2016-12/documents/scc_tsd_2010.pdf"
price_table <- data.frame(
  year = years, scc_2007_usd_per_metric_tonne_co2 = prices,
  scenario = "3% discount rate, average",
  source = c(rep(early_source, 2), rep(epa_table, 3)),
  price_basis = c(
    "User-specified 21.40; applied to 2005, not a published 2005 table entry",
    "User-specified 21.40 from 2010 TSD",
    rep("EPA snapshot table: emission-year-specific price", 3)
  )
)

input_paths <- file.path(input_dir, sprintf("forest-carbon-stock-%d.tif", years))
metadata_paths <- sub("\\.tif$", ".json", input_paths)
output_paths <- file.path(output_dir, sprintf("gep-carbon-stock-%d.tif", years))
output_metadata <- sub("\\.tif$", ".json", output_paths)
missing <- c(input_paths, metadata_paths)[!file.exists(c(input_paths, metadata_paths))]
if (length(missing)) {
  stop("Missing carbon-stock inputs (basal area is not carbon): ", paste(missing, collapse = ", "))
}
price_path <- file.path(output_dir, "carbon-stock-scc-prices.csv")
existing <- c(output_paths, output_metadata, price_path)
if (any(file.exists(existing))) {
  stop("Refusing to overwrite GEP outputs: ", paste(existing[file.exists(existing)], collapse = ", "))
}

validate_input <- function(index) {
  metadata <- jsonlite::fromJSON(metadata_paths[index])
  if (!identical(metadata$attribute, "Live aboveground tree carbon") ||
      !identical(metadata$units, "metric tonnes C per hectare") ||
      !identical(as.integer(metadata$target_year), years[index])) {
    stop("Unexpected carbon quantity, units, or year: ", metadata_paths[index])
  }
  carbon <- terra::rast(input_paths[index])
  if (terra::nlyr(carbon) != 2L || !nzchar(terra::crs(carbon))) {
    stop("Expected georeferenced two-band carbon-stock raster: ", input_paths[index])
  }
  carbon <- carbon[[1]]
  minimum <- terra::global(carbon, "min", na.rm = TRUE)[1, 1]
  maximum <- terra::global(carbon, "max", na.rm = TRUE)[1, 1]
  if (!is.finite(minimum) || !is.finite(maximum) || minimum < 0) {
    stop("Carbon raster must contain finite nonnegative quantities: ", input_paths[index])
  }
  list(carbon = carbon, metadata = metadata)
}

# Preflight all years before writing any valuation outputs.
reference <- NULL
for (index in seq_along(years)) {
  input <- validate_input(index)
  if (is.null(reference)) {
    reference <- input$carbon
  } else {
    terra::compareGeom(reference, input$carbon, stopOnError = TRUE)
  }
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(output_dir)) {
  stop("Could not create output directory: ", output_dir)
}
warning(
  "SCC-valued live aboveground carbon is a storage-value proxy, not annual GEP ",
  "sequestration or total ecosystem carbon. 2005 uses the user-specified 2010 price."
)
utils::write.csv(price_table, price_path, row.names = FALSE)
for (index in seq_along(years)) {
  input <- validate_input(index)
  carbon <- input$carbon
  # SCC prices CO2, whereas the stock layer measures elemental carbon.
  multiplier <- (44 / 12) * prices[index]
  value <- carbon * multiplier
  names(value) <- "carbon_stock_value_2007_USD_per_hectare"
  terra::writeRaster(
    value, output_paths[index], overwrite = FALSE,
    datatype = "FLT8S", NAflag = -9999,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES", "BIGTIFF=IF_SAFER")
  )
  saved <- terra::rast(output_paths[index])
  terra::compareGeom(carbon, saved, stopOnError = TRUE)
  difference <- terra::global(abs(saved - value), "max", na.rm = TRUE)[1, 1]
  mask_difference <- terra::global(
    is.na(saved) != is.na(carbon), "sum", na.rm = TRUE
  )[1, 1]
  if (!is.finite(difference) || difference > 1e-8 || mask_difference != 0) {
    stop("Written valuation raster failed value/NoData verification: ", output_paths[index])
  }
  jsonlite::write_json(
    list(
      target_year = years[index],
      quantity = "Live aboveground tree carbon",
      input_raster = input_paths[index],
      input_units = "metric tonnes C per hectare",
      output_units = "2007 USD per hectare",
      formula = "carbon_tonnes_C_per_ha * (44/12) * SCC_2007_USD_per_tonne_CO2",
      scc = price_table[index, ],
      forest_layer_metadata = input$metadata,
      interpretation = paste(
        "Stock-value proxy assuming CO2-equivalent avoided emissions;",
        "not annual sequestration GEP; excludes other carbon pools;",
        "not inflation-adjusted or summed across years."
      )
    ),
    output_metadata[index], auto_unbox = TRUE, pretty = TRUE, na = "null"
  )
  message("Saved ", output_paths[index], " (2007 USD/hectare).")
}
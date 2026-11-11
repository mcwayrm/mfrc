required_packages <- c("terra", "sf", "tmap", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") {
  stop("Carbon maps require tmap >= 4.0.")
}
source_paths <- vapply(sys.frames(), function(frame) {
  if (is.character(frame$ofile) && length(frame$ofile) == 1L) {
    frame$ofile
  } else {
    NA_character_
  }
}, character(1))
source_paths <- source_paths[!is.na(source_paths)]
if (length(source_paths)) {
  script_path <- tail(source_paths, 1L)
} else {
  script_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_argument) != 1L) {
    stop("Run this script with Rscript or source() using its file path.")
  }
  script_path <- sub("^--file=", "", script_argument)
}
root_dir <- dirname(dirname(normalizePath(script_path, mustWork = TRUE)))
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
boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
inputs <- c(input_paths, metadata_paths, boundary_path)
missing <- inputs[!file.exists(inputs)]
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
  if (terra::nlyr(carbon) != 2L || !nzchar(terra::crs(carbon)) ||
      !terra::same.crs(carbon, "EPSG:26915") || any(terra::res(carbon) != 300)) {
    stop("Expected two-band, 300-metre EPSG:26915 carbon-stock raster: ", input_paths[index])
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
boundary <- sf::st_read(boundary_path, quiet = TRUE)
if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) {
  stop("Minnesota boundary must have a CRS and valid geometry.")
}
boundary <- sf::st_sf(geometry = sf::st_union(sf::st_transform(boundary, 26915)))
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

visualize_carbon <- function() {
  map_dir <- file.path(root_dir, "outputs", "maps")
  figure_dir <- file.path(root_dir, "outputs", "figures")
  table_dir <- file.path(root_dir, "outputs", "tables")
  for (directory in c(map_dir, figure_dir, table_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create output directory: ", directory)
  }
  carbon_layers <- lapply(seq_along(years), function(index) validate_input(index)$carbon)
  value_layers <- lapply(output_paths, terra::rast)
  cell_area_ha <- prod(terra::res(carbon_layers[[1]])) / 10000
  summary <- do.call(rbind, lapply(seq_along(years), function(index) {
    carbon <- terra::values(carbon_layers[[index]], mat = FALSE)
    value <- terra::values(value_layers[[index]], mat = FALSE)
    valid <- !is.na(carbon)
    expected <- carbon[valid] * (44 / 12) * prices[index]
    if (!identical(is.na(carbon), is.na(value)) ||
        any(!is.finite(value[valid])) ||
        any(abs(value[valid] - expected) > 1e-8)) {
      stop("Carbon/value inputs disagree for visualization year ", years[index])
    }
    metadata <- jsonlite::fromJSON(metadata_paths[index])
    total_carbon <- sum(carbon[valid]) * cell_area_ha
    total_value <- sum(value[valid]) * cell_area_ha
    if (!isTRUE(all.equal(total_value, total_carbon * (44 / 12) * prices[index],
                          tolerance = 1e-12))) {
      stop("Area-weighted total value failed conversion verification for ", years[index])
    }
    data.frame(
      target_year = years[index],
      measurement_start_year = metadata$measurement_start_year,
      measurement_end_year = metadata$measurement_end_year,
      esa_mask_year = metadata$mask_year,
      forest_cells = sum(valid), mapped_forest_area_ha = sum(valid) * cell_area_ha,
      total_live_aboveground_carbon_metric_tonnes_C = total_carbon,
      mean_carbon_metric_tonnes_C_per_ha = mean(carbon[valid]),
      scc_2007_USD_per_metric_tonne_CO2 = prices[index],
      total_stock_value_2007_USD = total_value,
      total_stock_value_at_fixed_2005_scc_2007_USD = total_carbon * (44 / 12) * prices[1],
      scope = "Statewide mapped forest; not restricted to family ownership",
      interpretation = "SCC-valued live aboveground carbon stock proxy; not annual GEP"
    )
  }))
  first <- value_layers[[1]]
  last <- value_layers[[length(years)]]
  change <- terra::ifel(
    !is.na(first) & !is.na(last) & first > 0, 100 * (last - first) / first, NA
  )
  first_values <- terra::values(first, mat = FALSE)
  last_values <- terra::values(last, mat = FALSE)
  change_values <- terra::values(change, mat = FALSE)
  shared <- !is.na(first_values) & !is.na(last_values)
  eligible <- shared & first_values > 0
  if (!any(eligible) || any(!is.finite(change_values[eligible])) ||
      any(!is.na(change_values[!eligible])) ||
      !isTRUE(all.equal(
        change_values[eligible],
        100 * (last_values[eligible] - first_values[eligible]) / first_values[eligible],
        tolerance = 1e-10
      ))) {
    stop("Carbon-value percentage-change raster failed shared-mask/formula verification.")
  }
  change_summary <- data.frame(
    first_year = years[1], last_year = tail(years, 1),
    shared_forest_cells = sum(shared), positive_baseline_cells = sum(eligible),
    zero_baseline_cells_excluded = sum(shared & first_values == 0),
    forest_only_in_2005_cells = sum(!is.na(first_values) & is.na(last_values)),
    forest_only_in_2025_cells = sum(is.na(first_values) & !is.na(last_values)),
    mean_pixel_percent_value_change = mean(change_values[eligible]),
    median_pixel_percent_value_change = stats::median(change_values[eligible]),
    increasing_cells = sum(change_values[eligible] > 0),
    decreasing_cells = sum(change_values[eligible] < 0),
    unchanged_cells = sum(change_values[eligible] == 0)
  )
  utils::write.csv(summary, file.path(table_dir, "carbon-stock-value-time-series.csv"),
                   row.names = FALSE)
  utils::write.csv(change_summary, file.path(table_dir, "carbon-stock-value-change-summary.csv"),
                   row.names = FALSE)
  names(last) <- "stock_value"
  names(change) <- "percent_value_change"
  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  map_base <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
    tmap::tm_polygons(fill = "#f0f0f0", col = NA)
  finish_map <- function(map, title, note) {
    map + tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
      tmap::tm_title(title, size = 0.9, fontface = "bold",
                     position = tmap::tm_pos_out("center", "top")) +
      tmap::tm_credits(note, size = 0.65,
                       position = tmap::tm_pos_out("left", "bottom")) +
      tmap::tm_scalebar(breaks = c(0, 100, 200), text.size = 0.6) +
      tmap::tm_layout(
        frame = FALSE, bg.color = "white", meta.margins = c(0.12, 0, 0.07, 0.25),
        legend.position = tmap::tm_pos_out("right", "center"),
        legend.title.size = 0.8, legend.text.size = 0.7
      )
  }
  maps <- list(
    "carbon-stock-value-hotspots-2025" = finish_map(
      map_base + tmap::tm_shape(last) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_continuous(
            values = rev(grDevices::hcl.colors(100, "YlOrRd")), value.na = NA,
            limits = c(0, max(last_values, na.rm = TRUE)), midpoint = NA
          ),
          col.legend = tmap::tm_legend(title = "Stock-value proxy\n2007 USD / ha", na.show = FALSE)
        ),
      "A. Carbon stock value: 2025",
      "Live aboveground carbon; SCC = $46/t CO2; ESA mask: 2022.\nHeat map of value density, not a statistical hotspot test or annual GEP."
    ),
    "carbon-stock-value-percent-change-2005-2025" = finish_map(
      map_base + tmap::tm_shape(change) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_intervals(
            breaks = c(-Inf, -50, -10, 10, 50, 100, 200, Inf), midpoint = NA,
            labels = c("< -50%", "-50 to -10%", "-10 to +10%", "+10 to +50%",
                       "+50 to +100%", "+100 to +200%", ">= +200%"),
            values = c("#b2182b", "#ef8a62", "#f7f7f7", "#d1e5f0",
                       "#92c5de", "#4393c3", "#2166ac"), value.na = NA
          ),
          col.legend = tmap::tm_legend(title = "Stock-value change\n% of 2005 value", na.show = FALSE)
        ),
      "B. Carbon stock value change: 2005-2025",
      "Shared forest pixels, positive 2005 value; 2025 mask: 2022.\nIncludes carbon AND SCC price change ($21.40 to $46/t CO2)."
    )
  )
  check_export <- function(path) {
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Visualization export failed: ", path)
    }
  }
  for (name in names(maps)) {
    for (extension in c("png", "pdf")) {
      path <- file.path(map_dir, paste0(name, ".", extension))
      tmap::tmap_save(maps[[name]], path, width = 7, height = 7, units = "in", dpi = 300)
      check_export(path)
    }
  }
  draw_series <- function() {
    old_par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old_par), add = TRUE)
    graphics::par(mfrow = c(2, 1), mar = c(4.5, 5, 2.5, 1),
                  oma = c(4, 0, 2, 0), las = 1)
    graphics::plot(
      years, summary$total_live_aboveground_carbon_metric_tonnes_C / 1e6,
      type = "o", pch = 16, col = "#238443", lwd = 2, xaxt = "n",
      xlab = "", ylab = "Million metric tonnes C",
      main = "Live aboveground carbon stock"
    )
    graphics::axis(1, at = years)
    values <- summary$total_stock_value_2007_USD / 1e9
    fixed <- summary$total_stock_value_at_fixed_2005_scc_2007_USD / 1e9
    graphics::plot(
      years, values, type = "o", pch = 16, col = "#2166ac", lwd = 2, xaxt = "n",
      ylim = range(c(values, fixed)), xlab = "Composite target year",
      ylab = "Billion 2007 USD", main = "Gross-ecosystem Product (GEP)"
    )
    graphics::axis(1, at = years)
    graphics::lines(years, fixed, type = "o", pch = 1, col = "#666666", lty = 2, lwd = 2)
    graphics::legend(
      "topleft", c("Year-specific SCC", "Fixed $21.40/t CO2 (quantity/area only)"),
      col = c("#2166ac", "#666666"), lty = c(1, 2), pch = c(16, 1), bty = "n", cex = 0.8
    )
    graphics::mtext("C. Statewide carbon stock and value", outer = TRUE,
                    side = 3, line = 0.3, font = 2)
    graphics::mtext(
      "Each point uses a five-year FIA window and its own forest mask; 2025 uses ESA 2022.\nExploratory totals; not family-ownership estimates. Do not sum stocks across years.",
      outer = TRUE, side = 1, line = 0.5, cex = 0.75
    )
  }
  for (extension in c("png", "pdf")) {
    path <- file.path(figure_dir, paste0("carbon-stock-value-time-series.", extension))
    if (extension == "png") {
      grDevices::png(path, width = 8, height = 8, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::pdf(path, width = 8, height = 8)
    }
    tryCatch(draw_series(), finally = grDevices::dev.off())
    check_export(path)
  }
  warning(
    "Carbon value change includes changing SCC prices; annual totals use different forest masks. ",
    "These are unvalidated statewide stock proxies, not annual GEP or family-forest totals."
  )
  message("Saved two carbon maps (PNG/PDF), the two-panel time series (PNG/PDF), ",
          "and two summary CSVs under outputs/. Visualization files replace prior versions.")
  invisible(summary)
}

visualize_carbon()
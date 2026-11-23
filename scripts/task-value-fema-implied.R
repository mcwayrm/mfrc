# FEMA implied ecosystem-service values for Minnesota forests.
#
# Quantity: forest acreage from the 300 m FIA/ESA forest layers in
# data/forest-layer (each forest pixel covers 9 hectares = 22.24 acres).
#
# Price: FEMA's 2022 update of its per-acre ecosystem-service values for the
# "Forest" land-cover category (Table 5, page 20, 2021 USD/acre/year):
# https://www.fema.gov/sites/default/files/documents/fema_ecosystem-service-value-updates_2022.pdf
# Forest-column marginal values: Aesthetic Value 1,477; Air Quality 711;
# Climate Regulation 199; Erosion Control 1,672; Existence Value 7,531;
# Flood and Storm Hazard Reduction 368; Recreation/Tourism 94; Water
# Filtration 435; Water Supply 103. FEMA assigns no forest value for
# Biological Control, Food Provisioning, Habitat, Pollination, or Research
# and Education. Components sum to $12,590/acre/year; FEMA's published
# rounded total is $12,589.
#
# The uniform per-acre bundle is applied to every forest pixel (national
# benefit transfer), so the implied GEP layer varies over time only through
# forest extent. Values are 2021 nominal USD with no inflation adjustment.

required_packages <- c("terra", "sf", "tmap", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") stop("FEMA maps require tmap >= 4.0.")
source_paths <- vapply(sys.frames(), function(frame) {
  if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else NA_character_
}, character(1))
source_paths <- source_paths[!is.na(source_paths)]
if (length(source_paths)) {
  script_path <- tail(source_paths, 1L)
} else {
  script_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_argument) != 1L) stop("Run with Rscript or source() using the script path.")
  script_path <- sub("^--file=", "", script_argument)
}
root_dir <- dirname(dirname(normalizePath(script_path, mustWork = TRUE)))

fema_source_url <- paste0("https://www.fema.gov/sites/default/files/documents/",
                          "fema_ecosystem-service-value-updates_2022.pdf")

# Forest column of FEMA Table 5 (2021 USD per acre per year); NA marks
# services for which FEMA assigns no forest value.
fema_forest_values <- data.frame(
  ecosystem_service = c(
    "Aesthetic Value", "Air Quality", "Biological Control", "Climate Regulation",
    "Erosion Control", "Existence Value", "Flood and Storm Hazard Reduction",
    "Food Provisioning", "Habitat", "Pollination", "Recreation/Tourism",
    "Research and Education", "Water Filtration", "Water Supply"
  ),
  usd_per_acre_per_year_2021 = c(
    1477, 711, NA, 199,
    1672, 7531, 368,
    NA, NA, NA, 94,
    NA, 435, 103
  )
)

run_fema_implied <- function(root_dir) {
  years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
  acres_per_hectare <- 2.4710538147
  total_per_acre <- sum(fema_forest_values$usd_per_acre_per_year_2021, na.rm = TRUE)
  if (abs(total_per_acre - 12590) > 1e-9) {
    stop("FEMA forest per-acre components must sum to $12,590 (published total $12,589).")
  }
  value_per_ha <- total_per_acre * acres_per_hectare
  input_dir <- file.path(root_dir, "data", "forest-layer")
  output_dir <- file.path(root_dir, "data", "GEP")
  map_dir <- file.path(root_dir, "outputs", "maps")
  figure_dir <- file.path(root_dir, "outputs", "figures")
  table_dir <- file.path(root_dir, "outputs", "tables")
  template_paths <- file.path(input_dir, sprintf("forest-layer-%d.tif", years))
  metadata_paths <- sub("\\.tif$", ".json", template_paths)
  value_paths <- file.path(output_dir, sprintf("gep-fema-implied-%d.tif", years))
  price_path <- file.path(output_dir, "fema-implied-forest-prices-2021.csv")
  service_table_path <- file.path(table_dir, "fema-implied-service-values.csv")
  series_path <- file.path(table_dir, "fema-implied-time-series.csv")
  boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
  inputs <- c(template_paths, metadata_paths, boundary_path)
  if (any(!file.exists(inputs))) {
    stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  }
  outputs <- c(value_paths, sub("\\.tif$", ".json", value_paths), price_path,
               service_table_path, series_path)
  if (any(file.exists(outputs))) {
    stop("Refusing to overwrite FEMA implied-value outputs: ",
         paste(outputs[file.exists(outputs)], collapse = ", "))
  }
  metadata <- lapply(metadata_paths, jsonlite::fromJSON)
  templates <- lapply(seq_along(years), function(index) {
    raster <- terra::rast(template_paths[index])
    info <- metadata[[index]]
    if (!terra::same.crs(raster, "EPSG:26915") || any(terra::res(raster) != 300) ||
        info$target_year != years[index] ||
        info$mask_year != if (years[index] == 2025L) 2022L else years[index]) {
      stop("Unexpected forest template or metadata: ", template_paths[index])
    }
    raster[[1]]
  })
  for (index in seq_along(years)[-1]) {
    terra::compareGeom(templates[[1]], templates[[index]], stopOnError = TRUE)
  }
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) {
    stop("Invalid Minnesota boundary.")
  }
  boundary <- sf::st_sf(geometry = sf::st_union(sf::st_transform(boundary, 26915)))
  for (directory in c(output_dir, map_dir, figure_dir, table_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create directory: ", directory)
  }
  warning("FEMA per-acre values are national benefit-transfer averages (2021 USD) ",
          "applied uniformly to every Minnesota forest pixel. The bundle overlaps ",
          "services valued elsewhere in this project (for example climate regulation ",
          "and recreation), so do not add this layer to the other GEP estimates.")
  price_table <- cbind(
    fema_forest_values,
    price_year = 2021L,
    land_cover_category = "Forest",
    source = paste0(fema_source_url, " (Table 5, page 20)"),
    assumption = paste("National per-acre benefit-transfer values applied uniformly",
                       "to all Minnesota forest pixels; no inflation adjustment")
  )
  utils::write.csv(price_table, price_path, row.names = FALSE)
  cell_area_ha <- prod(terra::res(templates[[1]])) / 10000
  summaries <- vector("list", length(years))
  values_by_year <- vector("list", length(years))
  for (index in seq_along(years)) {
    year <- years[index]
    template <- templates[[index]]
    cells <- which(!is.na(terra::values(template, mat = FALSE)))
    if (length(cells) != metadata[[index]]$forest_pixels) {
      stop("Template forest count mismatch: ", year)
    }
    forest_area_ha <- length(cells) * cell_area_ha
    forest_area_acres <- forest_area_ha * acres_per_hectare
    value <- terra::rast(template)
    terra::values(value) <- NA_real_
    value[cells] <- value_per_ha
    names(value) <- "annual_fema_implied_value_USD_per_ha_per_year"
    terra::writeRaster(value, value_paths[index], overwrite = FALSE, datatype = "FLT8S",
                       NAflag = -9999, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3",
                                                "TILED=YES"))
    saved <- terra::rast(value_paths[index])
    terra::compareGeom(value, saved, stopOnError = TRUE)
    actual <- terra::values(saved, mat = FALSE)
    expected <- terra::values(value, mat = FALSE)
    if (!identical(is.na(actual), is.na(expected)) ||
        any(abs(actual[!is.na(actual)] - expected[!is.na(expected)]) > 1e-8)) {
      stop("Written FEMA value raster failed verification: ", value_paths[index])
    }
    service_values <- fema_forest_values$usd_per_acre_per_year_2021
    names(service_values) <- fema_forest_values$ecosystem_service
    info <- list(
      target_year = year,
      price_year = 2021L,
      currency = "2021 nominal USD",
      value_units = "USD per forest hectare per year",
      quantity = "Forest acreage from the 300 m FIA/ESA forest layer",
      formula = paste("uniform $", format(total_per_acre, big.mark = ","),
                      "/acre/year FEMA forest bundle *", round(acres_per_hectare, 7),
                      "acres/ha on every forest pixel"),
      fema_usd_per_acre_per_year = as.list(service_values[!is.na(service_values)]),
      total_usd_per_acre_per_year = total_per_acre,
      published_rounded_total_usd_per_acre_per_year = 12589,
      price_source = paste0(fema_source_url, " (Table 5, page 20)"),
      mask_template = template_paths[index],
      esa_mask_year = metadata[[index]]$mask_year,
      forest_pixels = length(cells),
      forest_area_ha = forest_area_ha,
      forest_area_acres = forest_area_acres,
      total_annual_value_USD = forest_area_acres * total_per_acre,
      interpretation = paste(
        "National FEMA benefit-transfer averages applied uniformly to Minnesota",
        "forest; spatial and temporal variation comes only from forest extent.",
        "The bundle includes services valued separately in this project (for",
        "example climate regulation, recreation, and water services), so this",
        "layer must not be added to the other GEP layers. Values are 2021",
        "nominal USD without Minnesota-specific calibration or ownership",
        "attribution; riparian-specific FEMA values are not used."
      )
    )
    jsonlite::write_json(info, sub("\\.tif$", ".json", value_paths[index]),
                         auto_unbox = TRUE, pretty = TRUE, digits = NA)
    summaries[[index]] <- data.frame(
      target_year = year,
      esa_mask_year = metadata[[index]]$mask_year,
      forest_pixels = length(cells),
      forest_area_ha = forest_area_ha,
      forest_area_acres = forest_area_acres,
      usd_per_acre_per_year_2021 = total_per_acre,
      total_annual_value_USD = forest_area_acres * total_per_acre
    )
    values_by_year[[index]] <- value
    message("Saved ", year, " FEMA implied-value raster: ",
            format(round(forest_area_acres), big.mark = ","), " forest acres, $",
            format(round(forest_area_acres * total_per_acre / 1e9, 2)), "B/year.")
  }
  summary <- do.call(rbind, summaries)
  summary$pct_change_from_2005 <-
    100 * (summary$total_annual_value_USD / summary$total_annual_value_USD[1] - 1)
  utils::write.csv(summary, series_path, row.names = FALSE)
  acres_2025 <- summary$forest_area_acres[summary$target_year == 2025L]
  service_table <- cbind(
    fema_forest_values,
    statewide_2025_million_usd_per_year =
      fema_forest_values$usd_per_acre_per_year_2021 * acres_2025 / 1e6
  )
  service_table <- rbind(service_table, data.frame(
    ecosystem_service = "Total",
    usd_per_acre_per_year_2021 = total_per_acre,
    statewide_2025_million_usd_per_year = total_per_acre * acres_2025 / 1e6
  ))
  utils::write.csv(service_table, service_table_path, row.names = FALSE)
  visualize_fema_implied(values_by_year, summary, value_per_ha, boundary,
                         map_dir, figure_dir)
  invisible(summary)
}

visualize_fema_implied <- function(values, summary, value_per_ha, boundary,
                                   map_dir, figure_dir) {
  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  last_year <- tail(summary$target_year, 1)
  # The FEMA bundle is uniform per acre, so the 2025 layer is mapped as a
  # single valued-forest class rather than a degenerate continuous scale.
  display <- terra::classify(values[[length(values)]], cbind(-Inf, Inf, 1))
  levels(display) <- data.frame(
    id = 1,
    valued = sprintf("Forest: $%s / ha / year",
                     format(round(value_per_ha), big.mark = ","))
  )
  map <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
    tmap::tm_polygons(fill = "#eeeeee", col = NA) +
    tmap::tm_shape(display) +
    tmap::tm_raster(
      col.scale = tmap::tm_scale_categorical(values = "#1b7837"),
      col.legend = tmap::tm_legend(title = "FEMA implied value", na.show = FALSE)
    ) +
    tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
    tmap::tm_title(paste0("FEMA implied forest ecosystem-service value: ", last_year),
                   size = 0.85, position = tmap::tm_pos_out("center", "top")) +
    tmap::tm_credits(
      paste0("Uniform FEMA 2022 forest bundle ($12,590/acre/year, 2021 USD) on the ",
             "ESA ", tail(summary$esa_mask_year, 1), " forest mask.\nNational ",
             "benefit-transfer averages; not additive with the other GEP layers."),
      size = 0.65, position = tmap::tm_pos_out("left", "bottom")
    ) +
    tmap::tm_layout(frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
                    legend.position = tmap::tm_pos_out("right", "center"),
                    legend.title.size = 0.8, legend.text.size = 0.7)
  for (extension in c("png", "pdf")) {
    path <- file.path(map_dir, paste0("fema-implied-value-", last_year, ".", extension))
    tmap::tmap_save(map, path, width = 7, height = 7, units = "in", dpi = 300)
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Map export failed: ", path)
    }
  }
  for (extension in c("png", "pdf")) {
    path <- file.path(figure_dir, paste0("fema-implied-time-series.", extension))
    if (extension == "png") {
      grDevices::png(path, width = 8, height = 8, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::pdf(path, width = 8, height = 8)
    }
    tryCatch({
      graphics::par(mfrow = c(2, 1), mar = c(4.5, 5, 2.5, 1), oma = c(4, 0, 2, 0), las = 1)
      acres <- summary$forest_area_acres / 1e6
      graphics::plot(summary$target_year, acres, type = "o", pch = 16, col = "#1b7837",
                     lwd = 2, xaxt = "n", xlab = "", ylab = "Million acres",
                     main = "A. Mapped Minnesota forest extent")
      graphics::axis(1, at = summary$target_year)
      totals <- summary$total_annual_value_USD / 1e9
      graphics::plot(summary$target_year, totals, type = "o", pch = 16, col = "#2166ac",
                     lwd = 2, xaxt = "n", xlab = "Target year",
                     ylab = "Billion USD / year (2021 USD)",
                     main = "B. FEMA implied annual forest ecosystem-service value")
      graphics::axis(1, at = summary$target_year)
      graphics::text(summary$target_year, totals,
                     sprintf("%+.1f%%", summary$pct_change_from_2005),
                     pos = c(3, 1, 3, 3, 1), cex = 0.75, col = "#444444")
      graphics::mtext("FEMA implied forest ecosystem-service value, 2005-2025",
                      side = 3, outer = TRUE, font = 2)
      graphics::mtext(
        paste("Uniform FEMA 2022 forest bundle ($12,590/acre/year, 2021 USD) times",
              "mapped forest acreage;\nvalue changes track forest extent only.",
              "Percent labels are changes from 2005."),
        side = 1, outer = TRUE, line = 0.5, cex = 0.75
      )
    }, finally = grDevices::dev.off())
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Time-series export failed: ", path)
    }
  }
}

run_fema_implied(root_dir)

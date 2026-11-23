# FEMA implied ecosystem-service values for Minnesota forests.
#
# Quantity: forest acreage from the 300 m FIA/ESA forest layers in
# data/forest-layer (each forest pixel covers 9 hectares = 22.24 acres).
#
# Price: FEMA's 2022 update of its per-acre ecosystem-service values for the
# Forest and Riparian land-cover categories (Table 5, page 20, 2021 USD/acre/year):
# https://www.fema.gov/sites/default/files/documents/fema_ecosystem-service-value-updates_2022.pdf
# A 50 m buffer around selected DNR streams, lakes, and wetlands identifies
# riparian area. Water polygon interiors are excluded. FEMA's riparian bundle
# replaces, rather than adds to, the forest bundle in this area. Fractions of
# each 300 m cell are estimated on a 30 m grid, assuming forest is distributed
# uniformly within each cell. Values are in 2021 USD, without inflation adjustment.

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

# FEMA Table 5 values (2021 USD per acre per year); NA marks an unreported value.
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
  ),
  riparian_usd_per_acre_per_year_2021 = c(
    767, 254, 199, 96,
    13823, NA, 6052,
    736, 2547, NA, 6215,
    NA, 6239, 272
  )
)

read_hydrography <- function(zip_path) {
  files <- utils::unzip(zip_path, list = TRUE)$Name
  entries <- grep(
    "^water_dnr_hydrography/fgdb/water_dnr_hydrography_uncompressed[.]gdb/[^/]+$",
    files, value = TRUE
  )
  if (!length(entries)) stop("No uncompressed DNR geodatabase in: ", zip_path)
  hydro_dir <- file.path(tempdir(), "fema-hydrography")
  dir.create(hydro_dir, recursive = TRUE, showWarnings = FALSE)
  utils::unzip(zip_path, files = entries, exdir = hydro_dir)
  gdb <- file.path(hydro_dir, "water_dnr_hydrography", "fgdb",
                   "water_dnr_hydrography_uncompressed.gdb")
  if (!dir.exists(gdb)) stop("Could not extract DNR geodatabase: ", gdb)
  streams <- sf::st_read(
    gdb, query = paste("SELECT Strm_type FROM dnr_rivers_and_streams",
                       "WHERE Strm_type IN (20, 21, 22, 62)"), quiet = TRUE
  )
  waters <- sf::st_read(
    gdb, query = paste("SELECT wb_class FROM dnr_hydro_features_all WHERE wb_class IN",
                       "('Lake or Pond', 'Reservoir', 'Riverine polygon',",
                       "'Innundation Area', 'Intermittent Water', 'Wetland')"),
    quiet = TRUE
  )
  if (!nrow(streams) || !nrow(waters) || is.na(sf::st_crs(streams)) ||
      is.na(sf::st_crs(waters))) {
    stop("Missing selected hydrography features or CRS.")
  }
  streams <- sf::st_transform(sf::st_zm(streams), 26915)
  waters <- sf::st_transform(sf::st_zm(waters), 26915)
  invalid <- which(!sf::st_is_valid(waters))
  if (length(invalid)) {
    message("Repairing ", length(invalid), " invalid water polygons.")
    sf::st_geometry(waters)[invalid] <- sf::st_geometry(sf::st_make_valid(waters[invalid, ]))
  }
  if (any(sf::st_geometry_type(waters) == "GEOMETRYCOLLECTION")) {
    waters <- sf::st_collection_extract(waters, "POLYGON", warn = FALSE)
  }
  if (any(!sf::st_is_valid(waters)) || any(sf::st_is_empty(streams)) ||
      any(sf::st_is_empty(waters))) {
    stop("Invalid or empty selected hydrography geometries.")
  }
  list(streams = streams, waters = waters)
}

riparian_fraction <- function(template, streams, waters, buffer_m = 50, factor = 10L) {
  fine <- terra::disagg(terra::rast(template), fact = factor)
  options <- list(datatype = "INT1U", NAflag = 255,
                  gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  waters <- terra::vect(waters)
  water_cells <- terra::rasterize(waters, fine, field = 1, background = 0,
                                  filename = tempfile(fileext = ".tif"), wopt = options)
  water_buffer <- terra::buffer(waters, width = buffer_m)
  water_buffer_cells <- terra::rasterize(water_buffer, fine, field = 1, background = 0,
                                         filename = tempfile(fileext = ".tif"), wopt = options)
  rm(waters, water_buffer)
  invisible(gc())
  stream_buffer <- terra::buffer(terra::vect(streams), width = buffer_m)
  stream_buffer_cells <- terra::rasterize(stream_buffer, fine, field = 1, background = 0,
                                          filename = tempfile(fileext = ".tif"), wopt = options)
  rm(stream_buffer)
  eligible <- terra::ifel(water_cells == 1, 0,
                         (water_buffer_cells == 1) | (stream_buffer_cells == 1),
                         filename = tempfile(fileext = ".tif"), wopt = options)
  fraction <- terra::aggregate(eligible, fact = factor, fun = "mean",
                               filename = tempfile(fileext = ".tif"),
                               wopt = list(datatype = "FLT4S",
                                           gdal = c("COMPRESS=DEFLATE", "TILED=YES")))
  terra::compareGeom(fraction, template, stopOnError = TRUE)
  names(fraction) <- "riparian_fraction"
  fraction
}

run_fema_implied <- function(root_dir) {
  years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
  acres_per_hectare <- 2.4710538147
  forest_rates <- fema_forest_values$usd_per_acre_per_year_2021
  riparian_rates <- fema_forest_values$riparian_usd_per_acre_per_year_2021
  forest_per_acre <- sum(forest_rates, na.rm = TRUE)
  riparian_per_acre <- sum(riparian_rates, na.rm = TRUE)
  if (forest_per_acre != 12590 || riparian_per_acre != 37200) {
    stop("Unexpected FEMA component sums for Forest or Riparian.")
  }
  forest_per_ha <- forest_per_acre * acres_per_hectare
  riparian_per_ha <- riparian_per_acre * acres_per_hectare
  input_dir <- file.path(root_dir, "data", "forest-layer")
  output_dir <- file.path(root_dir, "data", "GEP")
  map_dir <- file.path(root_dir, "outputs", "maps")
  figure_dir <- file.path(root_dir, "outputs", "figures")
  table_dir <- file.path(root_dir, "outputs", "tables")
  template_paths <- file.path(input_dir, sprintf("forest-layer-%d.tif", years))
  metadata_paths <- sub("\\.tif$", ".json", template_paths)
  value_paths <- file.path(output_dir, sprintf("gep-fema-implied-riparian-%d.tif", years))
  fraction_path <- file.path(root_dir, "data", "forest-riparian", "fema-riparian-fraction-50m.tif")
  price_path <- file.path(output_dir, "fema-implied-riparian-prices-2021.csv")
  service_table_path <- file.path(table_dir, "fema-implied-riparian-service-values.csv")
  series_path <- file.path(table_dir, "fema-implied-riparian-time-series.csv")
  map_paths <- file.path(map_dir, paste0("fema-implied-riparian-value-2025.", c("png", "pdf")))
  figure_paths <- file.path(figure_dir, paste0("fema-implied-riparian-time-series.", c("png", "pdf")))
  boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
  hydro_path <- file.path(root_dir, "data", "DNR-hydropgraphy", "water_dnr_hydrography.zip")
  inputs <- c(template_paths, metadata_paths, boundary_path, hydro_path)
  if (any(!file.exists(inputs))) {
    stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  }
  outputs <- c(value_paths, sub("\\.tif$", ".json", value_paths), fraction_path,
               sub("\\.tif$", ".json", fraction_path), price_path,
               service_table_path, series_path, map_paths, figure_paths)
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
  for (directory in unique(dirname(outputs))) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create directory: ", directory)
  }
  warning("FEMA per-acre values are national benefit-transfer averages (2021 USD). ",
          "The forest and riparian bundles replace one another and overlap ",
          "services valued elsewhere in this project (for example climate regulation ",
          "and recreation), so do not add this layer to the other GEP estimates.")
  hydro <- read_hydrography(hydro_path)
  hydro_counts <- c(streams = nrow(hydro$streams), water_polygons = nrow(hydro$waters))
  fraction <- riparian_fraction(templates[[1]], hydro$streams, hydro$waters)
  rm(hydro)
  invisible(gc())
  fraction <- round(fraction, digits = 2)
  terra::writeRaster(fraction, fraction_path, overwrite = FALSE, datatype = "FLT8S",
                     NAflag = -9999, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
  fraction <- terra::rast(fraction_path)
  fraction_values <- terra::values(fraction, mat = FALSE)
  if (any(!is.finite(fraction_values)) ||
      any(fraction_values < 0 | fraction_values > 1)) {
    stop("Invalid riparian fractions.")
  }
  jsonlite::write_json(list(
    source = hydro_path, buffer_m = 50, sample_resolution_m = 30,
    stream_types = c(20, 21, 22, 62),
    water_classes = c("Lake or Pond", "Reservoir", "Riverine polygon",
                      "Innundation Area", "Intermittent Water", "Wetland"),
    stream_features = unname(hydro_counts["streams"]),
    water_features = unname(hydro_counts["water_polygons"]),
    method = paste("50 m buffers unioned by 30 m cell-center rasterization;",
                   "water polygon interiors excluded; mean of 10 by 10 cells",
                   "gives the 300 m fraction. Assumes forest is uniformly",
                   "distributed within each forest pixel."),
    temporal_scope = "One current DNR hydrography snapshot is used for all target years."
  ), sub("\\.tif$", ".json", fraction_path), auto_unbox = TRUE, pretty = TRUE, digits = NA)
  price_table <- cbind(
    fema_forest_values,
    price_year = 2021L,
    source = paste0(fema_source_url, " (Table 5, page 20)"),
    assumption = paste("The Riparian bundle replaces Forest within the buffer;",
                       "unreported components are omitted, not assumed physically absent.")
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
    riparian_area_ha <- sum(fraction_values[cells]) * cell_area_ha
    nonriparian_area_ha <- forest_area_ha - riparian_area_ha
    forest_area_acres <- forest_area_ha * acres_per_hectare
    riparian_area_acres <- riparian_area_ha * acres_per_hectare
    nonriparian_area_acres <- nonriparian_area_ha * acres_per_hectare
    total_value <- nonriparian_area_acres * forest_per_acre +
      riparian_area_acres * riparian_per_acre
    value <- terra::rast(template)
    terra::values(value) <- NA_real_
    value[cells] <- forest_per_ha + fraction_values[cells] * (riparian_per_ha - forest_per_ha)
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
    if (!isTRUE(all.equal(sum(actual, na.rm = TRUE) * cell_area_ha, total_value,
                          tolerance = 1e-10))) {
      stop("Saved FEMA raster does not reproduce the statewide value: ", year)
    }
    forest_components <- stats::setNames(forest_rates, fema_forest_values$ecosystem_service)
    riparian_components <- stats::setNames(riparian_rates, fema_forest_values$ecosystem_service)
    info <- list(
      target_year = year,
      price_year = 2021L,
      currency = "2021 nominal USD",
      value_units = "USD per forest hectare per year",
      quantity = "FIA/ESA forest pixels partitioned by a DNR hydrography buffer",
      formula = paste("2.4710538147 * ((1 - riparian_fraction) * 12590 +",
                      "riparian_fraction * 37200) USD/ha/year"),
      forest_usd_per_acre_per_year = as.list(forest_components[!is.na(forest_components)]),
      riparian_usd_per_acre_per_year = as.list(riparian_components[!is.na(riparian_components)]),
      forest_total_usd_per_acre_per_year = forest_per_acre,
      riparian_total_usd_per_acre_per_year = riparian_per_acre,
      published_forest_total_usd_per_acre_per_year = 12589,
      published_riparian_total_usd_per_acre_per_year = 37199,
      price_source = paste0(fema_source_url, " (Table 5, page 20)"),
      mask_template = template_paths[index],
      riparian_fraction_raster = fraction_path,
      esa_mask_year = metadata[[index]]$mask_year,
      forest_pixels = length(cells),
      forest_area_ha = forest_area_ha,
      forest_area_acres = forest_area_acres,
      riparian_forest_area_ha = riparian_area_ha,
      nonriparian_forest_area_ha = nonriparian_area_ha,
      total_annual_value_USD = total_value,
      interpretation = paste(
        "National FEMA benefit-transfer averages applied to Minnesota forest;",
        "the Riparian bundle replaces Forest within 50 m of selected water features.",
        "Forest is assumed uniformly distributed within each 300 m cell;",
        "a fixed DNR hydrography snapshot is used for all years.",
        "Unreported service values are omitted, not assumed physically absent.",
        "The bundle includes services valued separately in this project (for",
        "example climate regulation, recreation, and water services), so this",
        "layer must not be added to the other GEP layers. Values are 2021",
        "nominal USD without Minnesota-specific calibration or ownership",
        "attribution."
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
      riparian_forest_area_ha = riparian_area_ha,
      nonriparian_forest_area_ha = nonriparian_area_ha,
      riparian_forest_area_acres = riparian_area_acres,
      nonriparian_forest_area_acres = nonriparian_area_acres,
      forest_usd_per_acre_per_year_2021 = forest_per_acre,
      riparian_usd_per_acre_per_year_2021 = riparian_per_acre,
      total_annual_value_USD = total_value
    )
    values_by_year[[index]] <- value
    message("Saved ", year, " FEMA implied-value raster: ",
            format(round(riparian_area_acres), big.mark = ","), " riparian forest acres, $",
            format(round(total_value / 1e9, 2)), "B/year.")
  }
  summary <- do.call(rbind, summaries)
  summary$pct_change_from_2005 <-
    100 * (summary$total_annual_value_USD / summary$total_annual_value_USD[1] - 1)
  utils::write.csv(summary, series_path, row.names = FALSE)
  result_2025 <- summary[summary$target_year == 2025L, ]
  forest_contribution <- forest_rates * result_2025$nonriparian_forest_area_acres / 1e6
  riparian_contribution <- riparian_rates * result_2025$riparian_forest_area_acres / 1e6
  combined <- ifelse(is.na(forest_contribution), 0, forest_contribution) +
    ifelse(is.na(riparian_contribution), 0, riparian_contribution)
  combined[is.na(forest_contribution) & is.na(riparian_contribution)] <- NA_real_
  service_table <- cbind(
    fema_forest_values,
    nonriparian_2025_million_usd_per_year = forest_contribution,
    riparian_2025_million_usd_per_year = riparian_contribution,
    combined_2025_million_usd_per_year = combined
  )
  service_table <- rbind(service_table, data.frame(
    ecosystem_service = "Total",
    usd_per_acre_per_year_2021 = forest_per_acre,
    riparian_usd_per_acre_per_year_2021 = riparian_per_acre,
    nonriparian_2025_million_usd_per_year = sum(forest_contribution, na.rm = TRUE),
    riparian_2025_million_usd_per_year = sum(riparian_contribution, na.rm = TRUE),
    combined_2025_million_usd_per_year = sum(combined, na.rm = TRUE)
  ))
  if (!isTRUE(all.equal(tail(service_table$combined_2025_million_usd_per_year, 1) * 1e6,
                        result_2025$total_annual_value_USD, tolerance = 1e-10))) {
    stop("Service table does not reproduce the statewide value.")
  }
  utils::write.csv(service_table, service_table_path, row.names = FALSE)
  visualize_fema_implied(values_by_year, summary, c(forest_per_ha, riparian_per_ha), boundary,
                         map_dir, figure_dir)
  invisible(summary)
}

visualize_fema_implied <- function(values, summary, value_per_ha, boundary,
                                   map_dir, figure_dir) {
  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  last_year <- tail(summary$target_year, 1)
  map <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
    tmap::tm_polygons(fill = "#eeeeee", col = NA) +
    tmap::tm_shape(values[[length(values)]]) +
    tmap::tm_raster(
      col.scale = tmap::tm_scale_continuous(
        limits = value_per_ha, values = c("#92bf65", "#27955e", "#173e91"), value.na = NA
      ),
      col.legend = tmap::tm_legend(title = "Annual value\n2021 USD / ha", na.show = FALSE)
    ) +
    tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
    tmap::tm_title(paste0("FEMA implied forest ecosystem-service value: ", last_year),
                   size = 0.85, position = tmap::tm_pos_out("center", "top")) +
    tmap::tm_credits(
      paste0("Forest and riparian FEMA bundles are weighted by 50 m buffer fraction; ",
             "ESA ", tail(summary$esa_mask_year, 1), " forest mask.\nNational ",
             "benefit-transfer averages; not additive with other GEP layers."),
      size = 0.65, position = tmap::tm_pos_out("left", "bottom")
    ) +
    tmap::tm_layout(frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
                    legend.position = tmap::tm_pos_out("right", "center"),
                    legend.title.size = 0.8, legend.text.size = 0.7)
  for (extension in c("png", "pdf")) {
    path <- file.path(map_dir, paste0("fema-implied-riparian-value-", last_year, ".", extension))
    tmap::tmap_save(map, path, width = 7, height = 7, units = "in", dpi = 300)
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Map export failed: ", path)
    }
  }
  for (extension in c("png", "pdf")) {
    path <- file.path(figure_dir, paste0("fema-implied-riparian-time-series.", extension))
    if (extension == "png") {
      grDevices::png(path, width = 8, height = 8, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::pdf(path, width = 8, height = 8)
    }
    tryCatch({
      graphics::par(mfrow = c(2, 1), mar = c(4.5, 5, 2.5, 1), oma = c(4, 0, 2, 0), las = 1)
      acres <- summary$forest_area_acres / 1e6
      riparian_acres <- summary$riparian_forest_area_acres / 1e6
      graphics::plot(summary$target_year, acres, type = "o", pch = 16, col = "#1b7837",
                     lwd = 2, xaxt = "n", ylim = c(0, max(acres) * 1.15),
                     xlab = "", ylab = "Million acres", main = "A. Mapped forest area")
      graphics::axis(1, at = summary$target_year)
      graphics::lines(summary$target_year, riparian_acres, type = "o", pch = 17,
                      col = "#2166ac", lwd = 2)
      graphics::legend("right", c("All forest", "Riparian forest"),
                       col = c("#1b7837", "#2166ac"), pch = c(16, 17), lty = 1, bty = "n")
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
        paste("FEMA Forest and Riparian bundles applied to mutually exclusive areas",
              "using one fixed 50 m hydrography buffer.\nPercent labels are changes",
              "from 2005; values are 2021 USD."),
        side = 1, outer = TRUE, line = 0.5, cex = 0.75
      )
    }, finally = grDevices::dev.off())
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Time-series export failed: ", path)
    }
  }
}

run_fema_implied(root_dir)

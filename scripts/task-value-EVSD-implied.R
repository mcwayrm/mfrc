# Exploratory ESVD implied values for Minnesota forests.
# ESVD standardized values are 2020 international dollars per hectare per year.
# Forest and selected all-biome observations are pooled only within each service.

required_packages <- c("terra", "sf", "tmap", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") stop("ESVD maps require tmap >= 4.0.")

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

read_esvd <- function(path) {
  data <- utils::read.csv(path, check.names = FALSE, colClasses = "character",
                          fileEncoding = "UTF-8", na.strings = "")
  required <- c("ValueId", "StudyId", "ESVD2.0_Biome_Codes", "ESVD2.0_Biome",
                "TEEB_ES", "ES_1", "Value_Change", "Incl_Excl",
                "Int$ Per Hectare Per Year", "Reference")
  if (!all(required %in% names(data))) {
    stop("Missing ESVD columns in ", path, ": ",
         paste(setdiff(required, names(data)), collapse = ", "))
  }
  if (anyNA(data$ValueId) || any(!nzchar(trimws(data$ValueId)))) {
    stop("ESVD value identifiers must be present: ", path)
  }
  duplicate_rows <- sum(duplicated(data))
  data <- data[!duplicated(data), , drop = FALSE]
  if (anyDuplicated(data$ValueId)) stop("Conflicting ESVD value identifiers in: ", path)
  message("Read ", nrow(data), " unique ESVD values from ", basename(path),
          "; removed ", duplicate_rows, " identical duplicate rows.")
  data
}

select_esvd_values <- function(forest, all_biomes) {
  forest_codes <- c("5", "6", "7")
  fallback_services <- c("4", "15", "16", "20", "21")
  prepare <- function(data, forest_only) {
    raw_value <- data[["Int$ Per Hectare Per Year"]]
    value <- suppressWarnings(as.numeric(raw_value))
    invalid <- !is.na(raw_value) & !is.finite(value)
    if (any(invalid)) stop("Non-numeric or nonfinite standardized ESVD value.")
    if (any(value < 0, na.rm = TRUE)) stop("Negative standardized ESVD values require review.")
    biome_codes <- strsplit(ifelse(is.na(data[["ESVD2.0_Biome_Codes"]]), "",
                                   data[["ESVD2.0_Biome_Codes"]]), ";", fixed = TRUE)
    biome_codes <- lapply(biome_codes, trimws)
    biome_ok <- vapply(biome_codes, function(codes) {
      length(codes) > 0L && all(nzchar(codes)) &&
        if (forest_only) all(codes %in% forest_codes) else length(codes) == 1L
    }, logical(1))
    service_code <- trimws(data$TEEB_ES)
    keep <- !is.na(value) & !is.na(service_code) & grepl("^[0-9]+$", service_code) &
      !is.na(data$ES_1) & nzchar(trimws(data$ES_1)) &
      !is.na(data$Value_Change) & data$Value_Change == "Total" &
      (is.na(data$Incl_Excl) | data$Incl_Excl != "0") & biome_ok
    data.frame(
      value_id = data$ValueId[keep], study_id = data$StudyId[keep],
      service_code = service_code[keep], service = trimws(data$ES_1[keep]),
      int2020_per_ha_per_year = value[keep],
      biome = data[["ESVD2.0_Biome"]][keep],
      reference = data$Reference[keep], stringsAsFactors = FALSE
    )
  }
  forest_values <- prepare(forest, TRUE)
  all_values <- prepare(all_biomes, FALSE)
  missing_services <- setdiff(fallback_services, unique(forest_values$service_code))
  fallback <- all_values[all_values$service_code %in% missing_services, , drop = FALSE]
  if (!nrow(forest_values) || !all(missing_services %in% fallback$service_code)) {
    stop("No forest values or missing selected all-biome fallback service values.")
  }
  forest_values$source <- "forest biomes"
  fallback$source <- rep("all biomes (fallback)", nrow(fallback))
  result <- rbind(forest_values, fallback)
  if (anyDuplicated(result$value_id)) stop("Selected ESVD value identifiers overlap.")
  labels <- split(result$service, result$service_code)
  if (any(vapply(labels, function(x) length(unique(x)) != 1L, logical(1)))) {
    stop("ESVD service codes do not have unique labels.")
  }
  result[order(as.integer(result$service_code), result$value_id), , drop = FALSE]
}

summarize_esvd <- function(values) {
  groups <- split(values, factor(values$service_code,
                                 levels = unique(values$service_code)))
  do.call(rbind, lapply(groups, function(group) {
    x <- group$int2020_per_ha_per_year
    q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
    data.frame(
      service_code = group$service_code[1], service = group$service[1],
      source = group$source[1], observations = length(x),
      studies = length(unique(group$study_id[!is.na(group$study_id)])),
      minimum = min(x), q1 = q[1], median = q[2], q3 = q[3], maximum = max(x),
      mean_int2020_per_ha_per_year = mean(x), stringsAsFactors = FALSE
    )
  }))
}

run_esvd_implied <- function(root_dir) {
  years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
  forest_path <- file.path(root_dir, "data", "EVSD", "EVSD-forests-ecoregions.csv")
  all_path <- file.path(root_dir, "data", "EVSD", "EVDS-all-ecoregions.csv")
  input_dir <- file.path(root_dir, "data", "forest-layer")
  output_dir <- file.path(root_dir, "data", "GEP")
  table_dir <- file.path(root_dir, "outputs", "tables")
  figure_dir <- file.path(root_dir, "outputs", "figures")
  map_dir <- file.path(root_dir, "outputs", "maps")
  boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
  template_paths <- file.path(input_dir, sprintf("forest-layer-%d.tif", years))
  metadata_paths <- sub("\\.tif$", ".json", template_paths)
  value_paths <- file.path(output_dir, sprintf("gep-evsd-implied-%d.tif", years))
  price_path <- file.path(output_dir, "evsd-implied-prices-2020.csv")
  records_path <- file.path(table_dir, "evsd-implied-selected-records.csv")
  service_path <- file.path(table_dir, "evsd-implied-service-values.csv")
  series_path <- file.path(table_dir, "evsd-implied-time-series.csv")
  box_paths <- file.path(figure_dir, paste0("evsd-implied-service-boxplots.", c("png", "pdf")))
  series_paths <- file.path(figure_dir, paste0("evsd-implied-time-series.", c("png", "pdf")))
  map_paths <- file.path(map_dir, paste0("evsd-implied-value-2025.", c("png", "pdf")))
  inputs <- c(forest_path, all_path, boundary_path, template_paths, metadata_paths)
  if (any(!file.exists(inputs))) {
    stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  }
  outputs <- c(value_paths, sub("\\.tif$", ".json", value_paths), price_path,
               records_path, service_path, series_path, box_paths, series_paths, map_paths)
  if (any(file.exists(outputs))) {
    stop("Refusing to overwrite ESVD implied-value outputs: ",
         paste(outputs[file.exists(outputs)], collapse = ", "))
  }
  values <- select_esvd_values(read_esvd(forest_path), read_esvd(all_path))
  prices <- summarize_esvd(values)
  bundle_rate <- sum(prices$mean_int2020_per_ha_per_year)
  if (!is.finite(bundle_rate) || bundle_rate <= 0) stop("Invalid ESVD bundle rate.")
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
  warning("ESVD global means are illustrative benefit transfers, not Minnesota-specific ",
          "prices. Service estimates can overlap; the bundle is not additive with FEMA ",
          "or other GEP estimates. International dollars are not nominal US dollars.")
  utils::write.csv(values, records_path, row.names = FALSE, na = "")
  utils::write.csv(prices, price_path, row.names = FALSE)
  cell_area_ha <- prod(terra::res(templates[[1]])) / 10000
  summaries <- vector("list", length(years))
  last_value <- NULL
  for (index in seq_along(years)) {
    year <- years[index]
    template <- templates[[index]]
    cells <- which(!is.na(terra::values(template, mat = FALSE)))
    if (length(cells) != metadata[[index]]$forest_pixels) {
      stop("Template forest count mismatch: ", year)
    }
    area <- length(cells) * cell_area_ha
    total <- area * bundle_rate
    value <- terra::rast(template)
    terra::values(value) <- NA_real_
    value[cells] <- bundle_rate
    names(value) <- "annual_esvd_implied_value_Int2020_per_ha_per_year"
    terra::writeRaster(value, value_paths[index], overwrite = FALSE, datatype = "FLT8S",
                       NAflag = -9999, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
    saved <- terra::rast(value_paths[index])
    terra::compareGeom(value, saved, stopOnError = TRUE)
    actual <- terra::values(saved, mat = FALSE)
    expected <- terra::values(value, mat = FALSE)
    if (!identical(is.na(actual), is.na(expected)) ||
        any(abs(actual[!is.na(actual)] - expected[!is.na(expected)]) > 1e-8) ||
        !isTRUE(all.equal(sum(actual, na.rm = TRUE) * cell_area_ha, total, tolerance = 1e-10))) {
      stop("Saved ESVD raster failed value, mask or total verification: ", value_paths[index])
    }
    info <- list(
      target_year = year, price_year = 2020L, currency = "2020 international dollars (PPP)",
      value_units = "Int$2020 per forest hectare per year",
      formula = "Sum of selected service-level arithmetic means, applied uniformly to forest pixels",
      bundle_int2020_per_ha_per_year = bundle_rate,
      source_version = "ESVD SEP2026V1.0", forest_source = forest_path,
      all_biomes_source = all_path, price_table = price_path, selected_records = records_path,
      selection = paste("Single service; total value; standardized value present;",
                        "Incl_Excl != 0; forest-only biome combinations; identical",
                        "ValueId duplicates removed; selected fallback services use",
                        "single-biome all-ecoregion observations only when forest values are absent."),
      mask_template = template_paths[index], esa_mask_year = metadata[[index]]$mask_year,
      forest_pixels = length(cells), forest_area_ha = area, total_annual_value_Int2020 = total,
      interpretation = paste("Illustrative global benefit transfer, not calibrated for Minnesota;",
                             "some service means may overlap and are not demonstrably additive.",
                             "Do not add to FEMA or other GEP layers. Rates are constant in time",
                             "and space, forest ownership is not attributed, and Int$ are not",
                             "nominal USD. ESVD advises against general summary-value transfer.")
    )
    jsonlite::write_json(info, sub("\\.tif$", ".json", value_paths[index]),
                         auto_unbox = TRUE, pretty = TRUE, digits = NA)
    summaries[[index]] <- data.frame(
      target_year = year, esa_mask_year = metadata[[index]]$mask_year,
      forest_pixels = length(cells), forest_area_ha = area,
      bundle_int2020_per_ha_per_year = bundle_rate, total_annual_value_Int2020 = total
    )
    if (year == 2025L) last_value <- saved
    message("Saved ", year, " ESVD implied-value raster.")
  }
  summary <- do.call(rbind, summaries)
  summary$pct_change_from_2005 <- 100 *
    (summary$total_annual_value_Int2020 / summary$total_annual_value_Int2020[1] - 1)
  prices$estimated_2025_million_int2020_per_year <-
    prices$mean_int2020_per_ha_per_year * tail(summary$forest_area_ha, 1) / 1e6
  if (!isTRUE(all.equal(sum(prices$estimated_2025_million_int2020_per_year) * 1e6,
                        tail(summary$total_annual_value_Int2020, 1), tolerance = 1e-10))) {
    stop("ESVD service values do not reproduce the statewide total.")
  }
  utils::write.csv(prices, service_path, row.names = FALSE)
  utils::write.csv(summary, series_path, row.names = FALSE)
  visualize_esvd(values, prices, summary, last_value, boundary, box_paths, series_paths, map_paths)
  invisible(summary)
}

visualize_esvd <- function(values, prices, summary, value, boundary,
                           box_paths, series_paths, map_paths) {
  labels <- paste0(prices$service, ifelse(prices$source == "forest biomes", "", " *"))
  groups <- lapply(prices$service_code, function(code) {
    log10(1 + values$int2020_per_ha_per_year[values$service_code == code])
  })
  names(groups) <- labels
  for (path in box_paths) {
    if (grepl("\\.png$", path)) {
      grDevices::png(path, width = 12, height = 10, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::pdf(path, width = 12, height = 10)
    }
    tryCatch({
      graphics::par(mar = c(5, 19, 3, 2), las = 1)
      graphics::boxplot(groups, horizontal = TRUE, xaxt = "n", col = "#a6dba0",
                        cex = 0.4, xlab = "Int$2020 / ha / year (log10(1 + value) scale)",
                        main = "Distribution of selected ESVD service values")
      ticks <- seq.int(0, ceiling(max(unlist(groups))))
      graphics::axis(1, at = ticks, labels = format(10^ticks - 1, big.mark = ",", scientific = FALSE))
      graphics::points(log10(1 + prices$mean_int2020_per_ha_per_year),
                       seq_along(groups), pch = 18, col = "#b2182b")
      graphics::mtext("* All-biome fallback; red diamonds are arithmetic means.",
                      side = 3, line = 0.2, cex = 0.8)
    }, finally = grDevices::dev.off())
    if (!file.exists(path) || file.info(path)$size == 0) stop("Boxplot export failed: ", path)
  }
  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  rate <- tail(summary$bundle_int2020_per_ha_per_year, 1)
  map <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
    tmap::tm_polygons(fill = "#eeeeee", col = NA) +
    tmap::tm_shape(value) +
    tmap::tm_raster(
      col.scale = tmap::tm_scale_intervals(
        breaks = c(0, rate * 1.001), values = "#238443",
        labels = format(round(rate), big.mark = ",", scientific = FALSE), value.na = NA
      ),
      col.legend = tmap::tm_legend(title = "Int$2020 / ha / year", na.show = FALSE)
    ) +
    tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
    tmap::tm_title("ESVD implied forest value: 2025", size = 0.85,
                   position = tmap::tm_pos_out("center", "top")) +
    tmap::tm_credits(
      "Uniform global service means on ESA 2022 forest mask; not locally calibrated.\nNot additive with FEMA or other GEP layers.",
      size = 0.65, position = tmap::tm_pos_out("left", "bottom")
    ) +
    tmap::tm_layout(frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
                    legend.position = tmap::tm_pos_out("right", "center"))
  for (path in map_paths) {
    tmap::tmap_save(map, path, width = 7, height = 7, units = "in", dpi = 300)
    if (!file.exists(path) || file.info(path)$size == 0) stop("Map export failed: ", path)
  }
  for (path in series_paths) {
    if (grepl("\\.png$", path)) {
      grDevices::png(path, width = 8, height = 5, units = "in", res = 300, type = "cairo")
    } else {
      grDevices::pdf(path, width = 8, height = 5)
    }
    tryCatch({
      graphics::par(mar = c(5, 5, 3, 1), las = 1)
      totals <- summary$total_annual_value_Int2020 / 1e9
      graphics::plot(summary$target_year, totals, type = "o", pch = 16, col = "#2166ac",
                     lwd = 2, xaxt = "n", xlab = "Target year",
                     ylab = "Billion Int$2020 / year",
                     main = "ESVD implied forest value")
      graphics::axis(1, at = summary$target_year)
      graphics::mtext("Fixed rates; changes reflect forest area, not changing service delivery.",
                      side = 3, line = 0.2, cex = 0.8)
    }, finally = grDevices::dev.off())
    if (!file.exists(path) || file.info(path)$size == 0) stop("Series export failed: ", path)
  }
}

if (sys.nframe() == 0L) run_esvd_implied(root_dir)

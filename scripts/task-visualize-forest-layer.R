required_packages <- c("terra", "sf", "tmap", "data.table", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") {
  stop("These maps require tmap >= 4.0.")
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

percent_change <- function(first, last) {
  terra::ifel(!is.na(first) & !is.na(last) & first > 0,
              100 * (last - first) / first, NA)
}

create_forest_maps <- function(root_dir) {
  years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
  forest_classes <- c(50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100)
  input_dir <- file.path(root_dir, "data", "forest-layer")
  map_dir <- file.path(root_dir, "outputs", "maps")
  table_dir <- file.path(root_dir, "outputs", "tables")
  boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
  plot_path <- file.path(root_dir, "data", "FIA", "MN_PLOT.csv")
  layer_paths <- file.path(input_dir, sprintf("forest-layer-%d.tif", years))
  metadata_paths <- sub("\\.tif$", ".json", layer_paths)
  inputs <- c(boundary_path, plot_path, layer_paths, metadata_paths)
  if (any(!file.exists(inputs))) {
    stop("Missing map inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  }

  metadata <- lapply(metadata_paths, jsonlite::fromJSON)
  layers <- lapply(seq_along(years), function(index) {
    info <- metadata[[index]]
    if (!identical(as.integer(info$target_year), years[index]) ||
        !identical(info$attribute, "Forest-condition-weighted COND.BALIVE") ||
        !identical(info$units, "ft2 per forest acre") ||
        !identical(as.integer(info$forest_classes), as.integer(forest_classes)) ||
        info$measurement_start_year != years[index] - 4L ||
        info$measurement_end_year != years[index] ||
        info$mask_year != if (years[index] == 2025L) 2022L else years[index]) {
      stop("Unexpected basal-area quantity, sampling window, or mask: ", metadata_paths[index])
    }
    raster <- terra::rast(layer_paths[index])
    if (terra::nlyr(raster) != 2L || !nzchar(terra::crs(raster)) ||
        !terra::same.crs(raster, "EPSG:26915") ||
        any(terra::res(raster) != 300)) {
      stop("Expected two-band, 300-metre EPSG:26915 forest layer: ", layer_paths[index])
    }
    raster <- raster[[1]]
    values <- terra::values(raster, mat = FALSE)
    valid <- values[!is.na(values)]
    if (!length(valid) || any(!is.finite(valid)) || any(valid < 0) ||
        length(valid) != info$forest_pixels) {
      stop("Invalid basal area or forest-pixel count: ", layer_paths[index])
    }
    raster
  })
  for (index in seq_along(layers)[-1]) {
    terra::compareGeom(layers[[1]], layers[[index]], stopOnError = TRUE)
  }
  message("Validated five annual basal-area rasters and their metadata.")
  template <- layers[[length(years)]]
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) {
    stop("Minnesota boundary must have a CRS and valid geometry.")
  }
  boundary <- sf::st_sf(geometry = sf::st_union(sf::st_transform(boundary, 26915)))
  # The producer predicts every forest-mask pixel; NoData also includes unavailable ESA.
  forest_mask <- terra::ifel(is.na(template), NA, 1)
  names(forest_mask) <- "esa_forest"

  keys <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT")
  plots <- data.table::fread(
    plot_path,
    select = c("CN", keys, "MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "LAT", "LON"),
    colClasses = list(character = c("CN", keys)), showProgress = FALSE
  )
  if (anyNA(plots$CN) || anyDuplicated(plots$CN) ||
      any(is.na(plots$STATECD) | plots$STATECD != "27")) {
    stop("FIA PLOT must contain unique CNs and Minnesota state codes only.")
  }
  for (column in c("MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "LAT", "LON")) {
    if (!is.numeric(plots[[column]])) {
      stop("FIA PLOT has a nonnumeric column: ", column)
    }
  }
  if (anyNA(plots$MEASYEAR)) {
    warning("Excluding ", sum(is.na(plots$MEASYEAR)), " FIA visits without MEASYEAR.")
  }
  plots <- plots[!is.na(MEASYEAR) & MEASYEAR >= 2000 & MEASYEAR <= 2025]
  if (!nrow(plots) || anyNA(plots[, c(keys, "INVYR"), with = FALSE])) {
    stop("No usable 2000-2025 FIA visits, or missing physical plot identifiers/INVYR.")
  }
  for (column in c("MEASMON", "MEASDAY")) {
    maximum <- if (column == "MEASMON") 12 else 31
    if (any(!is.na(plots[[column]]) &
            (plots[[column]] < 0 | plots[[column]] > maximum))) {
      stop("Invalid FIA measurement date component: ", column)
    }
    data.table::set(plots, which(is.na(plots[[column]])), column, 0L)
  }
  frequencies <- plots[, .(measurement_years = data.table::uniqueN(MEASYEAR)), by = keys]
  data.table::setorderv(plots, c(keys, "MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "CN"))
  latest <- plots[, .SD[.N], by = keys]
  latest <- merge(latest, frequencies, by = keys, sort = FALSE)
  valid_coordinates <- is.finite(latest$LAT) & is.finite(latest$LON) &
    latest$LAT >= -90 & latest$LAT <= 90 & latest$LON >= -180 & latest$LON <= 180
  if (any(!valid_coordinates)) {
    warning("Excluding ", sum(!valid_coordinates), " FIA sites with invalid latest coordinates.")
  }
  sites <- sf::st_transform(
    sf::st_as_sf(latest[valid_coordinates], coords = c("LON", "LAT"), crs = 4326), 26915
  )
  in_state <- lengths(sf::st_intersects(sites, boundary)) > 0L
  if (any(!in_state)) {
    warning("Excluding ", sum(!in_state), " public FIA locations outside the Minnesota boundary.")
  }
  mapped_keys <- data.table::as.data.table(sf::st_drop_geometry(sites[in_state, ]))[, ..keys]
  mapped_frequency <- merge(frequencies, mapped_keys, by = keys)
  frequency_table <- merge(
    frequencies[, .(physical_sites = .N), by = measurement_years],
    mapped_frequency[, .(mapped_sites = .N), by = measurement_years],
    by = "measurement_years", all.x = TRUE
  )
  frequency_table[is.na(mapped_sites), mapped_sites := 0L]
  frequency_table[, excluded_locations := physical_sites - mapped_sites]
  data.table::setorder(frequency_table, measurement_years)
  sites <- sites[in_state, ]
  if (!nrow(sites)) {
    stop("No FIA sites have valid Minnesota public coordinates.")
  }
  frequency_levels <- sort(unique(frequencies$measurement_years))
  sites$frequency <- factor(sites$measurement_years, levels = frequency_levels)
  sites <- sites[order(sites$measurement_years, decreasing = TRUE), ]
  message("Prepared ", nrow(sites), " physical FIA plot locations.")

  change <- percent_change(layers[[1]], template)
  names(change) <- "percent_change"
  first_values <- terra::values(layers[[1]], mat = FALSE)
  last_values <- terra::values(template, mat = FALSE)
  change_values <- terra::values(change, mat = FALSE)
  comparable <- !is.na(first_values) & !is.na(last_values)
  eligible <- comparable & first_values > 0
  if (!any(eligible) || any(!is.finite(change_values[eligible])) ||
      any(!is.na(change_values[!eligible]))) {
    stop("Percentage change must be finite only on shared forest pixels with positive 2005 density.")
  }
  if (!isTRUE(all.equal(
    change_values[eligible],
    100 * (last_values[eligible] - first_values[eligible]) / first_values[eligible],
    tolerance = 1e-10
  ))) {
    stop("Percentage-change raster failed formula verification.")
  }
  cell_area_ha <- prod(terra::res(template)) / 10000
  annual_table <- data.table::rbindlist(lapply(seq_along(years), function(index) {
    values <- terra::values(layers[[index]], mat = FALSE)
    values <- values[!is.na(values)]
    data.table::data.table(
      target_year = years[index],
      measurement_start_year = metadata[[index]]$measurement_start_year,
      measurement_end_year = metadata[[index]]$measurement_end_year,
      esa_mask_year = metadata[[index]]$mask_year,
      mapped_forest_cells = length(values),
      mapped_forest_area_ha = length(values) * cell_area_ha,
      mean_basal_area_ft2_per_forest_acre = mean(values),
      median_basal_area_ft2_per_forest_acre = stats::median(values),
      p10_basal_area_ft2_per_forest_acre = unname(stats::quantile(values, 0.1)),
      p90_basal_area_ft2_per_forest_acre = unname(stats::quantile(values, 0.9)),
      fia_sites_used = metadata[[index]]$sample_statistics$unique_sites
    )
  }))
  change_table <- data.table::data.table(
    first_target_year = 2005L, last_target_year = 2025L, last_esa_mask_year = 2022L,
    shared_forest_cells = sum(comparable),
    positive_baseline_cells = sum(eligible),
    zero_baseline_cells_excluded = sum(comparable & first_values == 0),
    forest_only_in_2005_cells = sum(!is.na(first_values) & is.na(last_values)),
    forest_only_in_2025_cells = sum(is.na(first_values) & !is.na(last_values)),
    comparable_area_ha = sum(eligible) * cell_area_ha,
    mean_pixel_percent_change = mean(change_values[eligible]),
    median_pixel_percent_change = stats::median(change_values[eligible]),
    p10_pixel_percent_change = unname(stats::quantile(change_values[eligible], 0.1)),
    p90_pixel_percent_change = unname(stats::quantile(change_values[eligible], 0.9)),
    increasing_cells = sum(change_values[eligible] > 0),
    decreasing_cells = sum(change_values[eligible] < 0)
  )

  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  map_base <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
    tmap::tm_polygons(fill = "#f0f0f0", col = NA)
  finish_map <- function(map, title, note) {
    map + tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
      tmap::tm_title(title, size = 1.1, fontface = "bold",
                     position = tmap::tm_pos_out("center", "top")) +
      tmap::tm_credits(note, size = 0.65,
                       position = tmap::tm_pos_out("left", "bottom")) +
      tmap::tm_scalebar(breaks = c(0, 100, 200), text.size = 0.6) +
      tmap::tm_layout(
        frame = FALSE, bg.color = "white",
        meta.margins = c(0.10, 0, 0.07, 0.23),
        legend.position = tmap::tm_pos_out("right", "center"),
        legend.title.size = 0.8, legend.text.size = 0.7
      )
  }
  maps <- list(
    "a-fia-measurement-frequency" = finish_map(
      map_base + tmap::tm_shape(sites) +
        tmap::tm_dots(
          fill = "frequency", size = 0.025,
          fill.scale = tmap::tm_scale_categorical(
            values = grDevices::hcl.colors(length(frequency_levels), "Viridis"),
            labels = paste(frequency_levels, ifelse(frequency_levels == 1, "year", "years"))
          ),
          fill.legend = tmap::tm_legend(title = "Measurement years\n2000-2025")
        ),
      "A. FIA measurement frequency",
      "One point per physical plot; latest public location, no jitter.\nAll plot statuses. FIA coordinates are fuzzed/swapped."
    ),
    "b-kriged-basal-area-2025" = finish_map(
      map_base + tmap::tm_shape(template) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_intervals(
            breaks = c(0, 50, 100, 150, 200, 300, Inf),
            values = rev(grDevices::hcl.colors(6, "YlGn")), value.na = NA
          ),
          col.legend = tmap::tm_legend(title = "Basal area\nft2 / forest acre", na.show = FALSE)
        ),
      "B. Kriged basal-area density: 2025",
      "FIA measurements: 2021-2025; ESA forest mask: 2022.\nExploratory interpolation, not design-based FIA estimates."
    ),
    "c-esa-forest-mask-2022" = finish_map(
      map_base + tmap::tm_shape(forest_mask) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_categorical(
            levels = 1, labels = "Forest",
            values = "#238443", value.na = NA
          ),
          col.legend = tmap::tm_legend(title = "ESA land cover\n2022 on 300-m grid", na.show = FALSE)
        ),
      "C. ESA forest mask used for 2025",
      "Exact ESA 2022 forest footprint preserved in the saved layer.\nGray: outside forest mask (including unavailable land cover)."
    ),
    "d-basal-area-percent-change-2005-2025" = finish_map(
      map_base + tmap::tm_shape(change) +
        tmap::tm_raster(
          col.scale = tmap::tm_scale_intervals(
            breaks = c(-Inf, -50, -25, -10, 10, 25, 50, Inf),
            labels = c("< -50%", "-50 to -25%", "-25 to -10%", "-10 to +10%",
                       "+10 to +25%", "+25 to +50%", ">= +50%"),
            midpoint = NA,
            values = c("#b2182b", "#d6604d", "#f4a582", "#f7f7f7",
                       "#92c5de", "#4393c3", "#2166ac"), value.na = NA
          ),
          col.legend = tmap::tm_legend(title = "Basal-area change\n% of 2005 density", na.show = FALSE)
        ),
      "D. Basal-area density change: 2005-2025",
      "100 x (2025 - 2005) / 2005; forest in both masks, 2005 > 0.\nFive-year composites; 2025 mask is 2022. Not forest-cover change."
    )
  )
  for (directory in c(map_dir, table_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create output directory: ", directory)
  }
  save_map <- function(map, stem, width, height) {
    for (extension in c("png", "pdf")) {
      path <- file.path(map_dir, paste0(stem, ".", extension))
      message("Exporting ", path)
      tmap::tmap_save(map, path, width = width, height = height, units = "in", dpi = 300)
      if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
        stop("Map export failed: ", path)
      }
    }
  }
  for (name in names(maps)) save_map(maps[[name]], name, 7, 7)
  panel <- tmap::tmap_arrange(
    maps[[1]], maps[[2]], maps[[3]], maps[[4]], ncol = 2, nrow = 2, sync = TRUE
  )
  save_map(panel, "forest-layer-four-panel", 14, 14)
  data.table::fwrite(frequency_table, file.path(table_dir, "fia-measurement-frequency.csv"))
  data.table::fwrite(annual_table, file.path(table_dir, "forest-layer-annual-summary.csv"))
  data.table::fwrite(change_table, file.path(table_dir, "forest-layer-change-summary.csv"))
  warning(
    "Maps and equal-cell summaries are exploratory, not design-based FIA estimates. ",
    "Basal area is a density proxy, not carbon or annual ecosystem-service supply. ",
    "Public coordinates and unvalidated kriging limit inference; small baselines amplify percent change."
  )
  message("Saved four individual maps and a four-panel figure (PNG/PDF) to ", map_dir,
          "; three summary CSVs to ", table_dir, ". Rerunning replaces these visualization outputs.")
  invisible(maps)
}

create_forest_maps(root_dir)

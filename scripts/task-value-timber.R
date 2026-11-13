# First-pass spatial timber provisioning valuation for Minnesota forests.
# The specification and limitations are documented in papers/report.tex.

required_packages <- c("terra", "sf", "tmap", "jsonlite", "data.table", "gstat", "rFIA")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") stop("Timber maps require tmap >= 4.0.")

source_paths <- vapply(sys.frames(), function(frame) {
  if (is.character(frame$ofile) && length(frame$ofile) == 1L) frame$ofile else NA_character_
}, character(1))
source_paths <- source_paths[!is.na(source_paths)]
if (length(source_paths)) {
  script_path <- tail(source_paths, 1L)
} else {
  script_argument <- grep("^--file=", commandArgs(), value = TRUE)
  if (length(script_argument) != 1L) {
    stop("Run with Rscript or source() using the script path.")
  }
  script_path <- sub("^--file=", "", script_argument)
}
root_dir <- dirname(dirname(normalizePath(script_path, mustWork = TRUE)))

format_fia_id <- function(x) {
  if (is.numeric(x)) sprintf("%.0f", x) else as.character(x)
}

stumpage_prices_2021 <- function() {
  data.table::data.table(
    product = c(
      rep("cord", 18L),
      rep("sawtimber", 12L)
    ),
    species_group = c(
      "Ash", "Aspen", "Balm", "Birch", "Oak", "Basswood", "Maple",
      "Mixed/Other Hdwds", "Balsam fir", "W Spruce", "B Spruce",
      "Tamarack", "Cedar", "Jack Pine", "Red Pine", "White Pine",
      "Hdwd Cd Composite", "Sftwd Cd Composite",
      "Ash", "Birch", "Elm", "Oak", "Basswood", "Maple",
      "Mixed/Other Hdwds", "W Spruce", "Jack Pine", "Red Pine",
      "White Pine", "MBF Composite"
    ),
    price_per_unit = c(
      8.69, 30.73, 26.26, 12.57, 27.88, 12.19, 12.12, 13.43,
      8.13, 17.54, 17.60, 6.17, 6.18, 23.78, 36.67, 33.10,
      26.73, 21.64,
      89.96, 80.54, 43.99, 108.64, 76.32, 110.28, 72.59,
      96.89, 103.91, 149.11, 109.09, 140.73
    ),
    unit = c(rep("USD per cord", 18L), rep("USD per Mbf", 12L)),
    price_year = 2021L,
    price_source = paste0(
      "Minnesota DNR, 2021 Minnesota Public Agencies Stumpage Price Review, ",
      "Tables 2 and 3"
    ),
    source_url = paste0(
      "https://files.dnr.state.mn.us/forestry/timber_sales/stumpage/",
      "stumpage-review-report-2021.pdf"
    )
  )
}

species_price_groups <- function(species_code, common_name) {
  common <- tolower(trimws(common_name))
  softwood <- species_code < 300L
  species_group <- data.table::fcase(
    grepl("balsam poplar", common), "Balm",
    grepl("aspen", common), "Aspen",
    grepl("ash", common), "Ash",
    grepl("birch", common), "Birch",
    grepl("oak", common), "Oak",
    grepl("basswood", common), "Basswood",
    grepl("maple", common), "Maple",
    grepl("elm", common), "Elm",
    grepl("balsam fir", common), "Balsam fir",
    grepl("white spruce", common), "W Spruce",
    grepl("black spruce", common), "B Spruce",
    grepl("tamarack", common), "Tamarack",
    grepl("cedar|redcedar", common), "Cedar",
    grepl("jack pine", common), "Jack Pine",
    grepl("red pine", common), "Red Pine",
    grepl("white pine", common), "White Pine",
    softwood, "Sftwd Cd Composite",
    default = "Mixed/Other Hdwds"
  )
  list(
    species_group = species_group,
    cord_group = data.table::fifelse(
      species_group %in% c(
        "Ash", "Aspen", "Balm", "Birch", "Oak", "Basswood", "Maple",
        "Mixed/Other Hdwds", "Balsam fir", "W Spruce", "B Spruce",
        "Tamarack", "Cedar", "Jack Pine", "Red Pine", "White Pine",
        "Sftwd Cd Composite"
      ),
      species_group,
      "Mixed/Other Hdwds"
    ),
    saw_group = data.table::fifelse(
      species_group %in% c(
        "Ash", "Birch", "Elm", "Oak", "Basswood", "Maple",
        "Mixed/Other Hdwds", "W Spruce", "Jack Pine", "Red Pine",
        "White Pine"
      ),
      species_group,
      "MBF Composite"
    )
  )
}

fit_timber_variogram <- function(sites) {
  if (nrow(sites) < 32L || stats::var(sites$annual_rate) <= 0) {
    stop("Insufficient spatial sites or variation for timber kriging.")
  }
  set.seed(42)
  fitting <- sites[sample.int(nrow(sites), min(1000L, nrow(sites))), ]
  coordinates <- sf::st_coordinates(fitting)
  cutoff <- max(stats::dist(coordinates))
  if (!is.finite(cutoff) || cutoff <= 0) stop("Invalid sample-coordinate spread.")
  empirical <- gstat::variogram(
    annual_rate ~ 1, fitting, cutoff = cutoff, width = cutoff / 20
  )
  empirical <- empirical[empirical$np >= 30, ]
  if (nrow(empirical) < 4L) stop("Insufficient populated timber variogram bins.")
  variance <- stats::var(fitting$annual_rate)
  distances <- as.numeric(stats::dist(coordinates))
  minimum_scale <- min(distances[distances > 0]) / (10 * cutoff)
  objective <- function(parameters) {
    modeled <- parameters[3] + parameters[1] *
      (1 - exp(-empirical$dist / (parameters[2] * cutoff)))
    sum((modeled - empirical$gamma / variance)^2 * empirical$np / sum(empirical$np))
  }
  fit <- stats::optim(
    c(0.9, 1 / 3, 0.1), objective, method = "L-BFGS-B",
    lower = c(1e-8, minimum_scale, 1e-8),
    upper = c(10, 10, 10),
    control = list(maxit = 2000, factr = 1e7)
  )
  if (fit$convergence != 0 || any(!is.finite(fit$par)) || !is.finite(fit$value)) {
    stop("Bounded timber variogram failed to converge: ", fit$message)
  }
  list(
    model = gstat::vgm(
      fit$par[1] * variance, "Exp", fit$par[2] * cutoff,
      fit$par[3] * variance
    ),
    fitting_sites = nrow(fitting),
    fitting_lags = nrow(empirical),
    normalized_fit_cost = fit$value
  )
}

krige_timber_rate <- function(sites, template, cells, model) {
  xy <- terra::xyFromCell(template, cells)
  prediction <- variance <- numeric(length(cells))
  for (start in seq.int(1L, length(cells), by = 20000L)) {
    rows <- start:min(start + 19999L, length(cells))
    queries <- sf::st_as_sf(
      data.frame(x = xy[rows, 1], y = xy[rows, 2]),
      coords = c("x", "y"), crs = 26915
    )
    estimate <- gstat::krige(
      annual_rate ~ 1, sites, queries, model = model,
      nmax = 32, debug.level = 0
    )
    prediction[rows] <- estimate$var1.pred
    variance[rows] <- estimate$var1.var
  }
  if (any(!is.finite(prediction)) || any(!is.finite(variance)) ||
      any(variance < -1e-8)) {
    stop("Nonfinite timber kriging estimates or negative prediction variance.")
  }
  clipped <- sum(prediction < 0)
  list(rate = pmax(prediction, 0), variance = pmax(variance, 0), clipped = clipped)
}

write_verified_raster <- function(raster, path) {
  partial <- sub("\\.tif$", ".partial.tif", path)
  if (file.exists(partial)) stop("Refusing to overwrite partial raster: ", partial)
  terra::writeRaster(
    raster, partial, overwrite = FALSE, datatype = "FLT8S", NAflag = -9999,
    gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES")
  )
  saved <- terra::rast(partial)
  terra::compareGeom(raster, saved, stopOnError = TRUE)
  for (band in seq_len(terra::nlyr(saved))) {
    actual <- terra::values(saved[[band]], mat = FALSE)
    expected <- terra::values(raster[[band]], mat = FALSE)
    if (!identical(is.na(actual), is.na(expected)) ||
        any(abs(actual[!is.na(actual)] - expected[!is.na(expected)]) > 1e-7)) {
      stop("Written timber raster failed verification: ", path)
    }
  }
  if (!file.rename(partial, path)) stop("Could not finalize timber raster: ", path)
  invisible(path)
}

save_timber_map <- function(raster, boundary, path, title, legend, note) {
  values <- terra::values(raster, mat = FALSE)
  maximum <- max(values, na.rm = TRUE)
  if (!is.finite(maximum) || maximum <= 0) stop("No positive values to map: ", path)
  map <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
    tmap::tm_polygons(fill = "#eeeeee", col = NA) +
    tmap::tm_shape(raster) +
    tmap::tm_raster(
      col.scale = tmap::tm_scale_continuous(
        limits = c(0, maximum),
        values = c("#ffffcc", "#78c679", "#006837"),
        value.na = NA
      ),
      col.legend = tmap::tm_legend(title = legend, na.show = FALSE)
    ) +
    tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
    tmap::tm_title(title, size = 0.85, position = tmap::tm_pos_out("center", "top")) +
    tmap::tm_credits(note, size = 0.65, position = tmap::tm_pos_out("left", "bottom")) +
    tmap::tm_layout(
      frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
      legend.position = tmap::tm_pos_out("right", "center"),
      legend.title.size = 0.8, legend.text.size = 0.7
    )
  for (extension in c("png", "pdf")) {
    map_path <- sub("\\.png$", paste0(".", extension), path)
    tmap::tmap_save(map, map_path, width = 7, height = 7, units = "in", dpi = 300)
    if (!file.exists(map_path) || is.na(file.info(map_path)$size) ||
        file.info(map_path)$size == 0) {
      stop("Timber map export failed: ", map_path)
    }
  }
}

run_timber_valuation <- function(root_dir) {
  target_year <- 2025L
  start_year <- target_year - 4L
  price_year <- 2021L
  conversion_cubic_feet_per_cord <- 80
  dnr_cords_per_mbf <- 2
  acre_to_hectare <- 2.4710538147
  forest_dir <- file.path(root_dir, "data", "forest-layer")
  fia_dir <- file.path(root_dir, "data", "FIA")
  timber_dir <- file.path(root_dir, "data", "forest-timber")
  gep_dir <- file.path(root_dir, "data", "GEP")
  map_dir <- file.path(root_dir, "outputs", "maps")
  table_dir <- file.path(root_dir, "outputs", "tables")
  forest_path <- file.path(forest_dir, "forest-layer-2025.tif")
  forest_metadata_path <- file.path(forest_dir, "forest-layer-2025.json")
  boundary_path <- file.path(
    root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp"
  )
  quantity_path <- file.path(timber_dir, "forest-timber-yield-2025.tif")
  value_path <- file.path(gep_dir, "gep-timber-2025.tif")
  plot_path <- file.path(timber_dir, "fia-timber-yield-plots-2025.csv")
  prices_path <- file.path(gep_dir, "timber-stumpage-prices-2021.csv")
  summary_path <- file.path(table_dir, "timber-valuation-summary-2025.csv")
  diagnostics_path <- file.path(table_dir, "timber-fia-sustainability-check-2025.csv")
  metadata_paths <- c(
    sub("\\.tif$", ".json", quantity_path),
    sub("\\.tif$", ".json", value_path)
  )
  map_paths <- file.path(
    map_dir,
    c(
      "timber-value-2025.png", "timber-value-2025.pdf",
      "timber-cord-yield-2025.png", "timber-cord-yield-2025.pdf",
      "timber-sawtimber-yield-2025.png", "timber-sawtimber-yield-2025.pdf"
    )
  )
  required_inputs <- c(
    forest_path, forest_metadata_path, boundary_path,
    file.path(fia_dir, c(
      "MN_PLOT.csv", "MN_COND.csv", "MN_POP_EVAL.csv", "MN_POP_EVAL_TYP.csv",
      "MN_TREE.csv", "MN_TREE_GRM_COMPONENT.csv", "MN_TREE_GRM_BEGIN.csv",
      "MN_TREE_GRM_MIDPT.csv", "MN_TREE_GRM_THRESHOLD.csv"
    ))
  )
  if (any(!file.exists(required_inputs))) {
    stop("Missing timber inputs: ", paste(required_inputs[!file.exists(required_inputs)], collapse = ", "))
  }
  outputs <- c(
    quantity_path, value_path, sub("\\.tif$", ".partial.tif", quantity_path),
    sub("\\.tif$", ".partial.tif", value_path), plot_path, prices_path,
    summary_path, diagnostics_path, metadata_paths, map_paths
  )
  if (any(file.exists(outputs))) {
    stop("Refusing to overwrite timber outputs: ",
         paste(outputs[file.exists(outputs)], collapse = ", "))
  }

  forest_metadata <- jsonlite::fromJSON(forest_metadata_path)
  template_stack <- terra::rast(forest_path)
  if (terra::nlyr(template_stack) != 2L ||
      !terra::same.crs(template_stack, "EPSG:26915") ||
      any(terra::res(template_stack) != 300) ||
      forest_metadata$target_year != target_year ||
      forest_metadata$mask_year != 2022L) {
    stop("Unexpected 2025 forest raster geometry or metadata.")
  }
  template <- template_stack[[1]]
  cells <- which(!is.na(terra::values(template, mat = FALSE)))
  if (!length(cells) || length(cells) != forest_metadata$forest_pixels) {
    stop("Forest raster mask does not match its metadata.")
  }

  evaluations <- data.table::fread(
    file.path(fia_dir, "MN_POP_EVAL.csv"),
    select = c("CN", "STATECD", "EVALID", "START_INVYR", "END_INVYR"),
    showProgress = FALSE
  )
  evaluation_types <- data.table::fread(
    file.path(fia_dir, "MN_POP_EVAL_TYP.csv"),
    select = c("EVAL_CN", "EVAL_TYP"), showProgress = FALSE
  )
  current_evaluations <- merge(
    evaluations[STATECD == 27 & END_INVYR == target_year],
    unique(evaluation_types[, .(EVAL_CN, EVAL_TYP)]),
    by.x = "CN", by.y = "EVAL_CN", all.x = TRUE, sort = FALSE
  )
  required_types <- c("EXPGROW", "EXPMORT", "EXPREMV")
  growth_evaluation <- unique(
    current_evaluations[EVAL_TYP %in% required_types, EVALID]
  )
  if (length(growth_evaluation) != 1L ||
      !all(required_types %in% current_evaluations[EVALID %in% growth_evaluation, EVAL_TYP])) {
    stop("Could not uniquely verify the 2021-2025 FIA growth/removal/mortality evaluation.")
  }

  plots <- data.table::fread(
    file.path(fia_dir, "MN_PLOT.csv"),
    select = c(
      "CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "MEASYEAR", "MEASMON",
      "MEASDAY", "INVYR", "PLOT_STATUS_CD", "LAT", "LON"
    ),
    colClasses = list(character = c("CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT")),
    showProgress = FALSE
  )
  conditions <- data.table::fread(
    file.path(fia_dir, "MN_COND.csv"),
    select = c("PLT_CN", "COND_STATUS_CD", "CONDPROP_UNADJ"),
    colClasses = list(character = "PLT_CN"), showProgress = FALSE
  )
  if (anyNA(plots$CN) || anyDuplicated(plots$CN) ||
      anyNA(conditions$PLT_CN) ||
      any(!is.na(conditions$CONDPROP_UNADJ) &
          (!is.finite(conditions$CONDPROP_UNADJ) |
             conditions$CONDPROP_UNADJ < 0 | conditions$CONDPROP_UNADJ > 1))) {
    stop("Invalid or duplicate FIA plot/condition records.")
  }
  if (anyNA(conditions$COND_STATUS_CD) || any(!is.finite(conditions$CONDPROP_UNADJ))) {
    stop("FIA conditions have missing status or forest-proportion values.")
  }
  plot_forest_shares <- conditions[COND_STATUS_CD == 1L, .(
    forest_share = sum(CONDPROP_UNADJ, na.rm = TRUE),
    valid_share = all(is.finite(CONDPROP_UNADJ) & CONDPROP_UNADJ >= 0 &
                        CONDPROP_UNADJ <= 1)
  ), by = PLT_CN]
  if (any(plot_forest_shares$forest_share > 1.001)) {
    stop("FIA forest-condition proportions exceed one.")
  }
  plots[, forest_share := plot_forest_shares$forest_share[
    match(CN, plot_forest_shares$PLT_CN)
  ]]
  plots[, valid_forest_share := plot_forest_shares$valid_share[
    match(CN, plot_forest_shares$PLT_CN)
  ]]
  plots[, `:=`(
    MEASYEAR = as.integer(MEASYEAR),
    MEASMON = as.integer(MEASMON),
    MEASDAY = as.integer(MEASDAY),
    INVYR = as.integer(INVYR),
    PLOT_STATUS_CD = as.integer(PLOT_STATUS_CD),
    LAT = as.numeric(LAT),
    LON = as.numeric(LON)
  )]
  visit_window <- plots[
    MEASYEAR >= start_year & MEASYEAR <= target_year
  ]
  if (!nrow(visit_window)) stop("No FIA plot visits in the 2021-2025 window.")
  date_columns <- c("MEASMON", "MEASDAY")
  for (column in date_columns) {
    x <- visit_window[[column]]
    maximum <- if (column == "MEASMON") 12L else 31L
    if (any(!is.na(x) & (!is.finite(x) | x < 0 | x > maximum | x %% 1 != 0))) {
      stop("Invalid FIA measurement date component: ", column)
    }
    data.table::set(visit_window, which(is.na(x)), column, 0L)
  }

  message("Loading Minnesota FIA tables for rFIA plot-level annual growth estimates.")
  fia_database <- rFIA::readFIA(dir = fia_dir, states = "MN", inMemory = FALSE)
  net_volume <- rFIA::growMort(
    fia_database, byPlot = TRUE, bySpecies = TRUE, landType = "forest",
    treeType = "gs", method = "annual", stateVar = "NETVOL"
  )
  saw_volume <- rFIA::growMort(
    fia_database, byPlot = TRUE, bySpecies = TRUE, landType = "forest",
    treeType = "gs", method = "annual", stateVar = "SAWVOL_BF"
  )
  net_volume <- data.table::as.data.table(net_volume)
  saw_volume <- data.table::as.data.table(saw_volume)
  net_required <- c(
    "PLT_CN", "YEAR", "SPCD", "COMMON_NAME", "RECR_NETVOL_ACRE",
    "MORT_NETVOL_ACRE", "REMV_NETVOL_ACRE", "GROW_NETVOL_ACRE"
  )
  saw_required <- c(
    "PLT_CN", "YEAR", "SPCD", "COMMON_NAME", "RECR_SAWVOL_BF_ACRE",
    "MORT_SAWVOL_BF_ACRE", "REMV_SAWVOL_BF_ACRE", "GROW_SAWVOL_BF_ACRE"
  )
  if (!all(net_required %in% names(net_volume)) ||
      !all(saw_required %in% names(saw_volume))) {
    stop("rFIA returned unexpected plot-level timber volume fields.")
  }
  net_volume[, `:=`(
    plot_id = format_fia_id(PLT_CN),
    YEAR = as.integer(YEAR),
    SPCD = as.integer(SPCD)
  )]
  saw_volume[, `:=`(
    plot_id = format_fia_id(PLT_CN),
    YEAR = as.integer(YEAR),
    SPCD = as.integer(SPCD)
  )]
  net_volume <- net_volume[YEAR >= start_year & YEAR <= target_year]
  saw_volume <- saw_volume[YEAR >= start_year & YEAR <= target_year]
  if (!nrow(net_volume) || !nrow(saw_volume)) {
    stop("rFIA returned no annualized merchantable or sawtimber records for 2021-2025.")
  }
  net_columns <- c(
    "GROW_NETVOL_ACRE", "RECR_NETVOL_ACRE", "MORT_NETVOL_ACRE",
    "REMV_NETVOL_ACRE"
  )
  saw_columns <- c(
    "GROW_SAWVOL_BF_ACRE", "RECR_SAWVOL_BF_ACRE", "MORT_SAWVOL_BF_ACRE",
    "REMV_SAWVOL_BF_ACRE"
  )
  for (column in net_columns) {
    if (any(!is.finite(net_volume[[column]]))) {
      stop("Nonfinite rFIA plot-level rate: ", column)
    }
  }
  for (column in saw_columns) {
    if (any(!is.finite(saw_volume[[column]]))) {
      stop("Nonfinite rFIA plot-level rate: ", column)
    }
  }
  if (anyNA(net_volume$SPCD) || anyNA(net_volume$COMMON_NAME) ||
      anyNA(saw_volume$SPCD) || anyNA(saw_volume$COMMON_NAME)) {
    stop("rFIA timber records have missing species identifiers.")
  }
  if (anyDuplicated(net_volume[, .(plot_id, YEAR, SPCD)]) ||
      anyDuplicated(saw_volume[, .(plot_id, YEAR, SPCD)])) {
    stop("rFIA returned duplicate plot-year-species records.")
  }

  plot_ids <- unique(c(net_volume$plot_id, saw_volume$plot_id))
  visit_window <- visit_window[CN %in% plot_ids]
  measurement_years <- visit_window[, .(plot_id = CN, MEASYEAR)]
  if (anyDuplicated(measurement_years$plot_id)) {
    stop("FIA plot identifiers map to multiple measurement years.")
  }
  for (records in list(net_volume, saw_volume)) {
    year_check <- merge(
      unique(records[, .(plot_id, YEAR)]), measurement_years,
      by = "plot_id", all.x = TRUE, sort = FALSE
    )
    if (anyNA(year_check$MEASYEAR) || any(year_check$YEAR != year_check$MEASYEAR)) {
      stop("rFIA annual growth records do not match FIA plot measurement years.")
    }
  }
  keys <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT")
  data.table::setorderv(
    visit_window,
    c(keys, "MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "CN")
  )
  selected_visits <- visit_window[, .SD[.N], by = keys]
  selected_visits[, plot_id := CN]
  selected_visits[, exclusion := ""]
  selected_visits[
    is.na(PLOT_STATUS_CD) | PLOT_STATUS_CD != 1L,
    exclusion := "visit_not_sampled"
  ]
  selected_visits[
    exclusion == "" & (is.na(valid_forest_share) | !valid_forest_share |
                         !is.finite(forest_share) | forest_share <= 0),
    exclusion := "no_valid_forest_condition"
  ]
  selected_visits[
    exclusion == "" & (!is.finite(LAT) | !is.finite(LON) |
                         LAT < -90 | LAT > 90 | LON < -180 | LON > 180),
    exclusion := "invalid_public_coordinates"
  ]
  selected_visits <- selected_visits[exclusion == ""]
  if (nrow(selected_visits) < 32L) {
    stop("Fewer than 32 valid sampled forest plots remain in the target window.")
  }
  selected_ids <- selected_visits$plot_id
  net_volume <- net_volume[plot_id %in% selected_ids]
  saw_volume <- saw_volume[plot_id %in% selected_ids]
  data.table::setnames(
    net_volume, net_columns,
    c("growth_cf_acre", "recruitment_cf_acre", "mortality_cf_acre", "removals_cf_acre")
  )
  data.table::setnames(
    saw_volume, saw_columns,
    c("growth_bf_acre", "recruitment_bf_acre", "mortality_bf_acre", "removals_bf_acre")
  )
  species_rates <- merge(
    net_volume[, .(
      plot_id, YEAR, SPCD, COMMON_NAME, growth_cf_acre,
      recruitment_cf_acre, mortality_cf_acre, removals_cf_acre
    )],
    saw_volume[, .(
      plot_id, YEAR, SPCD, saw_common_name = COMMON_NAME, growth_bf_acre,
      recruitment_bf_acre, mortality_bf_acre, removals_bf_acre
    )],
    by = c("plot_id", "YEAR", "SPCD"), all = TRUE, sort = FALSE
  )
  species_rates[, COMMON_NAME := data.table::fifelse(
    is.na(COMMON_NAME) | COMMON_NAME == "", saw_common_name, COMMON_NAME
  )]
  for (column in c(
    "growth_cf_acre", "recruitment_cf_acre", "mortality_cf_acre",
    "removals_cf_acre", "growth_bf_acre", "recruitment_bf_acre",
    "mortality_bf_acre", "removals_bf_acre"
  )) {
    data.table::set(species_rates, which(is.na(species_rates[[column]])), column, 0)
  }
  if (any(!is.finite(as.matrix(species_rates[, c(
    "growth_cf_acre", "recruitment_cf_acre", "mortality_cf_acre",
    "removals_cf_acre", "growth_bf_acre", "recruitment_bf_acre",
    "mortality_bf_acre", "removals_bf_acre"
  ), with = FALSE])))) {
    stop("Invalid joined FIA timber rates.")
  }
  groups <- species_price_groups(species_rates$SPCD, species_rates$COMMON_NAME)
  species_rates[, `:=`(
    species_group = groups$species_group,
    cord_price_group = groups$cord_group,
    saw_price_group = groups$saw_group
  )]
  prices <- stumpage_prices_2021()
  cord_prices <- prices[product == "cord", .(
    cord_price_group = species_group, cord_price = price_per_unit
  )]
  saw_prices <- prices[product == "sawtimber", .(
    saw_price_group = species_group, saw_price = price_per_unit
  )]
  species_rates <- merge(species_rates, cord_prices, by = "cord_price_group",
                         all.x = TRUE, sort = FALSE)
  species_rates <- merge(species_rates, saw_prices, by = "saw_price_group",
                         all.x = TRUE, sort = FALSE)
  if (anyNA(species_rates$cord_price) || anyNA(species_rates$saw_price)) {
    stop("A FIA species group has no matched 2021 DNR stumpage price.")
  }
  species_rates[, `:=`(
    sustainable_cf_acre = pmax(
      growth_cf_acre + recruitment_cf_acre - mortality_cf_acre, 0
    ),
    sustainable_bf_acre = pmax(
      growth_bf_acre + recruitment_bf_acre - mortality_bf_acre, 0
    ),
    observed_removals_cf_acre = pmax(removals_cf_acre, 0),
    observed_removals_bf_acre = pmax(removals_bf_acre, 0)
  )]
  species_rates[, `:=`(
    sustainable_mbf_acre = sustainable_bf_acre / 1000,
    sustainable_cords_acre = pmax(
      sustainable_cf_acre / conversion_cubic_feet_per_cord -
        sustainable_bf_acre / 1000 * dnr_cords_per_mbf,
      0
    ),
    observed_removals_mbf_acre = observed_removals_bf_acre / 1000,
    observed_removals_cords_acre = pmax(
      observed_removals_cf_acre / conversion_cubic_feet_per_cord -
        observed_removals_bf_acre / 1000 * dnr_cords_per_mbf,
      0
    )
  )]
  species_rates[, `:=`(
    sustainable_cords_ha = sustainable_cords_acre * acre_to_hectare,
    sustainable_mbf_ha = sustainable_mbf_acre * acre_to_hectare,
    observed_removals_cords_ha = observed_removals_cords_acre * acre_to_hectare,
    observed_removals_mbf_ha = observed_removals_mbf_acre * acre_to_hectare
  )]
  species_rates[, value_2021_usd_ha := (
    sustainable_cords_ha * cord_price +
      sustainable_mbf_ha * saw_price
  )]
  plot_rates <- species_rates[, .(
    sustainable_cords_ha = sum(sustainable_cords_ha),
    sustainable_mbf_ha = sum(sustainable_mbf_ha),
    value_2021_usd_ha = sum(value_2021_usd_ha),
    observed_removals_cords_ha = sum(observed_removals_cords_ha),
    observed_removals_mbf_ha = sum(observed_removals_mbf_ha)
  ), by = plot_id]
  plot_rates <- merge(
    selected_visits[, .(
      plot_id, STATECD, UNITCD, COUNTYCD, PLOT, MEASYEAR, INVYR,
      LAT, LON, forest_share
    )],
    plot_rates, by = "plot_id", all = FALSE, sort = FALSE
  )
  coordinates <- sf::st_transform(
    sf::st_as_sf(plot_rates, coords = c("LON", "LAT"), crs = 4326, remove = FALSE),
    26915
  )
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) {
    stop("Invalid Minnesota boundary geometry or CRS.")
  }
  boundary <- sf::st_sf(
    geometry = sf::st_union(sf::st_transform(boundary, 26915))
  )
  inside <- lengths(sf::st_intersects(coordinates, boundary)) > 0L
  outside_count <- sum(!inside)
  if (sum(inside) < 32L) stop("Fewer than 32 FIA sample plots fall inside Minnesota.")
  plot_rates <- data.table::as.data.table(sf::st_drop_geometry(coordinates[inside, ]))
  xy <- sf::st_coordinates(coordinates[inside, ])
  plot_rates[, `:=`(x = xy[, 1], y = xy[, 2])]
  site_rates <- plot_rates[, .(
    sustainable_cords_ha = mean(sustainable_cords_ha),
    sustainable_mbf_ha = mean(sustainable_mbf_ha),
    value_2021_usd_ha = mean(value_2021_usd_ha)
  ), by = .(x, y)]
  sites <- sf::st_as_sf(site_rates, coords = c("x", "y"), crs = 26915)
  if (qr(scale(sf::st_coordinates(sites), center = TRUE, scale = FALSE))$rank < 2) {
    stop("FIA timber sample locations are spatially collinear.")
  }
  fits <- lapply(
    c("sustainable_cords_ha", "sustainable_mbf_ha", "value_2021_usd_ha"),
    function(field) {
      rates <- sites
      rates$annual_rate <- rates[[field]]
      fit_timber_variogram(rates)
    }
  )
  names(fits) <- c("cords", "sawtimber", "value")
  effective_ranges_m <- vapply(fits, function(item) item$model$range[2], numeric(1))
  short_range_warning <- any(effective_ranges_m < min(terra::res(template)))
  if (short_range_warning) {
    warning(
      "At least one fitted timber variogram range is below the 300-m forest-grid ",
      "resolution; do not interpret pixel-scale value hotspots as validated spatial patterns."
    )
  }
  message(
    "Prepared ", nrow(plot_rates), " FIA plot visits at ", nrow(sites),
    " unique coordinates for the 2021-2025 interval."
  )

  for (directory in c(timber_dir, gep_dir, map_dir, table_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create output directory: ", directory)
  }
  price_table_path <- "https://www.dnr.state.mn.us/forestry/timbersales/stumpage.html"
  prices[, price_page := price_table_path]
  data.table::fwrite(prices, prices_path)
  data.table::fwrite(plot_rates, plot_path)

  fields <- c("sustainable_cords_ha", "sustainable_mbf_ha", "value_2021_usd_ha")
  estimates <- lapply(seq_along(fields), function(index) {
    sites_with_rate <- sites
    sites_with_rate$annual_rate <- sites_with_rate[[fields[index]]]
    krige_timber_rate(
      sites_with_rate, template, cells, fits[[index]]$model
    )
  })
  names(estimates) <- c("cords", "sawtimber", "value")
  rate_layers <- lapply(estimates, function(item) {
    raster <- terra::rast(template)
    terra::values(raster) <- NA_real_
    raster[cells] <- item$rate
    raster
  })
  names(rate_layers) <- names(estimates)
  variance_layer <- terra::rast(template)
  terra::values(variance_layer) <- NA_real_
  variance_layer[cells] <- estimates$value$variance
  names(rate_layers$cords) <- "annual_sustainable_cord_products_cords_per_ha"
  names(rate_layers$sawtimber) <- "annual_sustainable_sawtimber_Mbf_per_ha"
  names(rate_layers$value) <- "annual_sustainable_stumpage_value_2021_USD_per_ha"
  names(variance_layer) <- "kriging_prediction_variance_USD_per_ha_squared"
  quantity_raster <- c(rate_layers$cords, rate_layers$sawtimber)
  names(quantity_raster) <- c(
    "annual_sustainable_cord_products_cords_per_ha",
    "annual_sustainable_sawtimber_Mbf_per_ha"
  )
  value_raster <- c(rate_layers$value, variance_layer)
  write_verified_raster(quantity_raster, quantity_path)
  write_verified_raster(value_raster, value_path)

  pixel_area_ha <- prod(terra::res(template)) / 10000
  mapped_area_ha <- length(cells) * pixel_area_ha
  summary <- data.frame(
    target_year = target_year,
    FIA_measurement_window_start = start_year,
    FIA_measurement_window_end = target_year,
    FIA_growth_evalid = growth_evaluation,
    price_year = price_year,
    price_currency = "2021 nominal USD",
    fia_sample_visits = nrow(plot_rates),
    unique_kriging_sites = nrow(sites),
    plots_outside_minnesota = outside_count,
    forest_pixels = length(cells),
    mapped_area_ha = mapped_area_ha,
    pixel_area_ha = pixel_area_ha,
    mean_annual_cords_per_ha = mean(estimates$cords$rate),
    total_annual_sustainable_cords = sum(estimates$cords$rate) * pixel_area_ha,
    mean_annual_sawtimber_Mbf_per_ha = mean(estimates$sawtimber$rate),
    total_annual_sustainable_sawtimber_Mbf = sum(estimates$sawtimber$rate) * pixel_area_ha,
    mean_annual_stumpage_value_2021_USD_per_ha = mean(estimates$value$rate),
    total_annual_stumpage_value_2021_USD = sum(estimates$value$rate) * pixel_area_ha,
    sample_mean_observed_removals_cords_per_ha = mean(plot_rates$observed_removals_cords_ha),
    sample_mean_observed_removals_Mbf_per_ha = mean(plot_rates$observed_removals_mbf_ha),
    sample_share_cord_removals_above_sustainable = mean(
      plot_rates$observed_removals_cords_ha > plot_rates$sustainable_cords_ha
    ),
    sample_share_sawtimber_removals_above_sustainable = mean(
      plot_rates$observed_removals_mbf_ha > plot_rates$sustainable_mbf_ha
    ),
    negative_kriging_predictions_clipped = sum(vapply(estimates, `[[`, integer(1), "clipped")),
    scope = "Statewide forest pixels; all ownerships; exploratory spatial interpolation"
  )
  utils::write.csv(summary, summary_path, row.names = FALSE)
  diagnostics <- plot_rates[, .(
    plot_id, MEASYEAR, forest_share,
    sustainable_cords_ha, observed_removals_cords_ha,
    sustainable_mbf_ha, observed_removals_mbf_ha
  )]
  diagnostics[, `:=`(
    cord_removals_above_sustainable =
      observed_removals_cords_ha > sustainable_cords_ha,
    sawtimber_removals_above_sustainable =
      observed_removals_mbf_ha > sustainable_mbf_ha
  )]
  data.table::fwrite(diagnostics, diagnostics_path)

  limitations <- c(
    "FIA plot coordinates are public, fuzzed/swapped coordinates; results are exploratory.",
    "Plot-level rFIA annual growth, recruitment, mortality, and removal rates are kriged without FIA population/stratum expansion weights or cross-validation.",
    "Annual sustainable harvest ceiling is max(growth on survivors + recruitment - natural mortality, 0); observed removals are a diagnostic and are not subtracted from the ceiling.",
    "Cord products are estimated as residual merchantable cubic volume after sawtimber deduction, using 80 solid cubic feet per cord and the DNR report's 2 cords per Mbf equivalence.",
    "Species-specific DNR price groups distinguish cord products and sawtimber, not detailed log grades; unlisted sawtimber species use the Mbf composite price.",
    "FIA board-foot volume is treated as compatible with the DNR Scribner Mbf price scale; confirm this scale conversion before final inference.",
    "2021 public-agency stumpage prices are applied to the latest 2021-2025 FIA annualized rates; no post-2021 prices or variable harvest/haul costs are available here.",
    "The result is a gross stumpage-value proxy, not net timber rent; costs are not subtracted from stumpage prices.",
    "No ownership overlay is applied. The map covers all masked forest pixels and must not be interpreted as family-forest-only value.",
    "The forest mask is the existing 2025 layer using the 2022 ESA land-cover mask; all forest pixels are predicted, including beyond local FIA sample support."
  )
  if (short_range_warning) {
    limitations <- c(
      limitations,
      paste0(
        "At least one fitted exponential variogram range (minimum ",
        signif(min(effective_ranges_m), 4),
        " m) is shorter than the 300-m forest grid; spatial hotspots are not validated."
      )
    )
  }
  metadata <- list(
    target_year = target_year,
    measurement_window = c(start_year, target_year),
    price_year = price_year,
    FIA_growth_evalid = growth_evaluation,
    FIA_EVAL_TYP = required_types,
    quantity_method = paste(
      "rFIA::growMort(method='annual', byPlot=TRUE, bySpecies=TRUE,",
      "treeType='gs', landType='forest'); annual ceiling=max(GROW+RECR-MORT,0)"
    ),
    products = list(
      cord_products = "Residual net cubic volume / 80 cf per cord, minus sawtimber at 2 cords per Mbf",
      sawtimber = "FIA SAWVOL_BF converted from board feet to Mbf by dividing by 1000"
    ),
    pricing = "2021 Minnesota DNR public-agency average stumpage prices by species group and product",
    price_source_page = price_table_path,
    price_source_report = unique(prices$source_url),
    currency = "2021 nominal USD",
    annual_value_units = "USD per forest hectare per year",
    quantity_raster = quantity_path,
    value_raster = value_path,
    mapped_pixels = length(cells),
    mapped_area_ha = mapped_area_ha,
    sample_visits = nrow(plot_rates),
    unique_kriging_sites = nrow(sites),
    outside_boundary = outside_count,
    effective_variogram_ranges_m = as.list(effective_ranges_m),
    variogram_range_below_pixel_resolution = short_range_warning,
    variograms = lapply(fits, function(x) list(
      model = as.data.frame(x$model),
      fitting_sites = x$fitting_sites,
      fitting_lags = x$fitting_lags,
      normalized_fit_cost = x$normalized_fit_cost
    )),
    kriging = "Local ordinary kriging with 32 nearest unique FIA sites; negative rates clipped to zero",
    mask_source = forest_metadata$mask_source,
    limitations = limitations
  )
  jsonlite::write_json(
    metadata, metadata_paths[1], auto_unbox = TRUE, pretty = TRUE, digits = NA
  )
  jsonlite::write_json(
    metadata, metadata_paths[2], auto_unbox = TRUE, pretty = TRUE, digits = NA
  )

  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  spatial_note <- if (short_range_warning) {
    paste0(
      "\nShortest fitted variogram range: ", signif(min(effective_ranges_m), 3),
      " m (<300-m grid); hotspots are not validated."
    )
  } else {
    ""
  }
  save_timber_map(
    rate_layers$value, boundary,
    file.path(map_dir, "timber-value-2025.png"),
    "Annual sustainable stumpage value: 2025",
    "2021 USD / ha / year",
    paste0(
      "FIA annualized growth ceiling; 2021 DNR prices.\n",
      "All forest ownerships; gross stumpage-value proxy.", spatial_note
    )
  )
  save_timber_map(
    rate_layers$cords, boundary,
    file.path(map_dir, "timber-cord-yield-2025.png"),
    "Annual sustainable cord-product yield: 2025",
    "Cords / ha / year",
    paste0(
      "Cord products are a residual-volume proxy; 80 cf/cord and 2 cords/Mbf.",
      spatial_note
    )
  )
  save_timber_map(
    rate_layers$sawtimber, boundary,
    file.path(map_dir, "timber-sawtimber-yield-2025.png"),
    "Annual sustainable sawtimber yield: 2025",
    "Mbf / ha / year",
    paste0(
      "FIA sawtimber board-foot volume; local ordinary kriging; not detailed grade.",
      spatial_note
    )
  )
  message(
    "Saved 2025 timber quantity/value rasters, maps, FIA plot diagnostics, and summaries."
  )
  invisible(summary)
}

if (sys.nframe() == 0L) run_timber_valuation(root_dir)

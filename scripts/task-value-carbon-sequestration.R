required_packages <- c("terra", "sf", "tmap", "jsonlite", "data.table", "gstat")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") stop("Carbon maps require tmap >= 4.0.")
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

annual_carbon_change <- function(first, last, interval) {
  if (any(!is.finite(first)) || any(!is.finite(last)) ||
      any(first < 0 | last < 0) || any(!is.finite(interval) | interval <= 0) ||
      length(first) != length(last) || length(first) != length(interval)) {
    stop("Annual stock change requires matching finite nonnegative stocks and positive intervals.")
  }
  (last - first) / interval
}

read_matched_carbon <- function(fia_dir) {
  keys <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT")
  plots <- data.table::fread(
    file.path(fia_dir, "MN_PLOT.csv"),
    select = c("CN", "PREV_PLT_CN", keys, "MEASYEAR", "MEASMON", "MEASDAY",
               "INVYR", "REMPER", "PLOT_STATUS_CD", "DESIGNCD", "LAT", "LON"),
    colClasses = list(character = c("CN", "PREV_PLT_CN", keys)), showProgress = FALSE
  )
  conditions <- data.table::fread(
    file.path(fia_dir, "MN_COND.csv"),
    select = c("PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ"),
    colClasses = list(character = c("PLT_CN", "CONDID")), showProgress = FALSE
  )
  if (anyNA(plots$CN) || anyDuplicated(plots$CN) ||
      any(is.na(plots$STATECD) | plots$STATECD != "27") ||
      anyNA(conditions[, c("PLT_CN", "CONDID"), with = FALSE]) ||
      anyDuplicated(conditions[, c("PLT_CN", "CONDID"), with = FALSE])) {
    stop("FIA plot/condition identifiers must be unique and present, with Minnesota plots only.")
  }
  endpoints <- plots[!is.na(MEASYEAR) & MEASYEAR >= 2001 & MEASYEAR <= 2025]
  if (!nrow(endpoints) || anyNA(endpoints[, c(keys, "INVYR"), with = FALSE])) {
    stop("No endpoint visits or missing physical identifiers/inventory year.")
  }
  for (column in c("MEASMON", "MEASDAY")) {
    maximum <- if (column == "MEASMON") 12 else 31
    x <- endpoints[[column]]
    if (any(!is.na(x) & (!is.finite(x) | x < 0 | x > maximum | x %% 1 != 0))) {
      stop("Invalid measurement date component: ", column)
    }
    data.table::set(endpoints, which(is.na(x)), column, 0L)
  }
  endpoint_ids <- endpoints$CN
  needed_ids <- unique(c(endpoint_ids, endpoints$PREV_PLT_CN))
  needed_ids <- needed_ids[!is.na(needed_ids)]
  conditions <- conditions[PLT_CN %in% needed_ids]
  condition_stats <- conditions[, .(
    valid_conditions = all(!is.na(COND_STATUS_CD)) &&
      all(is.finite(CONDPROP_UNADJ) & CONDPROP_UNADJ >= 0 & CONDPROP_UNADJ <= 1),
    all_forest = all(!is.na(COND_STATUS_CD) & COND_STATUS_CD == 1),
    forest_share = sum(CONDPROP_UNADJ[COND_STATUS_CD == 1], na.rm = TRUE)
  ), by = PLT_CN]
  trees <- data.table::fread(
    file.path(fia_dir, "MN_TREE.csv"),
    select = c("PLT_CN", "CONDID", "STATUSCD", "CARBON_AG", "TPA_UNADJ"),
    colClasses = list(character = c("PLT_CN", "CONDID")), showProgress = FALSE
  )
  if (!nrow(trees)) stop("FIA TREE contains no records.")
  trees <- trees[PLT_CN %in% needed_ids]
  if (anyNA(trees$STATUSCD)) stop("Relevant FIA TREE records have missing live/dead status.")
  trees <- trees[STATUSCD == 1]
  trees <- merge(trees, conditions[, .(PLT_CN, CONDID, COND_STATUS_CD)],
                 by = c("PLT_CN", "CONDID"), all.x = TRUE, sort = FALSE)
  if (anyNA(trees$COND_STATUS_CD)) stop("Live-tree records have missing matching condition status.")
  trees <- trees[COND_STATUS_CD == 1]
  totals <- trees[, .(
    valid_carbon = all(is.finite(CARBON_AG) & CARBON_AG >= 0 &
                       is.finite(TPA_UNADJ) & TPA_UNADJ > 0),
    carbon_lb_per_plot_acre = sum(CARBON_AG * TPA_UNADJ)
  ), by = PLT_CN]
  condition_stats <- merge(condition_stats, totals, by = "PLT_CN", all.x = TRUE, sort = FALSE)
  # A sampled, fully forested visit without live trees has zero live-tree carbon.
  condition_stats[is.na(valid_carbon), `:=`(valid_carbon = TRUE, carbon_lb_per_plot_acre = 0)]
  condition_stats[, eligible_stock := valid_conditions & all_forest & valid_carbon &
                    abs(forest_share - 1) <= 1e-6]
  condition_stats[, carbon_tonnes_per_ha := carbon_lb_per_plot_acre / forest_share *
                    (0.00045359237 / 0.40468564224)]
  stocks <- merge(plots, condition_stats, by.x = "CN", by.y = "PLT_CN", all.x = TRUE, sort = FALSE)
  previous <- stocks[, c("CN", keys, "MEASYEAR", "PLOT_STATUS_CD", "DESIGNCD",
                        "eligible_stock", "carbon_tonnes_per_ha"), with = FALSE]
  data.table::setnames(previous, names(previous), paste0("previous_", names(previous)))
  endpoints <- stocks[CN %in% endpoint_ids]
  for (column in c("MEASMON", "MEASDAY")) {
    data.table::set(endpoints, which(is.na(endpoints[[column]])), column, 0L)
  }
  pairs <- merge(endpoints, previous, by.x = "PREV_PLT_CN", by.y = "previous_CN",
                 all.x = TRUE, sort = FALSE)
  pairs[, exclusion := ""]
  exclude <- function(test, reason) {
    rows <- which(pairs$exclusion == "" & (is.na(test) | test))
    data.table::set(pairs, rows, "exclusion", reason)
  }
  exclude(is.na(pairs$previous_MEASYEAR), "missing_previous_visit_or_measurement_year")
  same_site <- Reduce(`&`, lapply(keys, function(key) {
    pairs[[key]] == pairs[[paste0("previous_", key)]]
  }))
  exclude(!same_site, "physical_plot_key_changed")
  exclude(pairs$PLOT_STATUS_CD != 1 |
            pairs$previous_PLOT_STATUS_CD != 1, "visit_not_sampled")
  exclude(pairs$DESIGNCD != pairs$previous_DESIGNCD, "sampling_design_changed")
  exclude(!pairs$eligible_stock | !pairs$previous_eligible_stock,
          "not_fully_forested_or_invalid_carbon")
  exclude(!is.finite(pairs$REMPER) | pairs$REMPER <= 0 |
            pairs$MEASYEAR <= pairs$previous_MEASYEAR, "invalid_remeasurement_interval")
  exclude(!is.finite(pairs$LAT) | !is.finite(pairs$LON) |
            pairs$LAT < -90 | pairs$LAT > 90 | pairs$LON < -180 | pairs$LON > 180,
          "invalid_public_coordinates")
  pairs[exclusion == "", annual_change := annual_carbon_change(
    previous_carbon_tonnes_per_ha, carbon_tonnes_per_ha, REMPER
  )]
  message("Prepared FIA endpoint records and linked live-carbon stocks.")
  pairs
}

fit_change_variogram <- function(sites) {
  set.seed(42)
  fitting <- sites[sample.int(nrow(sites), min(1000L, nrow(sites))), ]
  coordinates <- sf::st_coordinates(fitting)
  cutoff <- max(stats::dist(coordinates))
  if (!is.finite(cutoff) || cutoff <= 0 || stats::var(fitting$annual_change) <= 0) {
    stop("Insufficient spatial/quantity variation for a sequestration variogram.")
  }
  empirical <- gstat::variogram(annual_change ~ 1, fitting, cutoff = cutoff, width = cutoff / 20)
  empirical <- empirical[empirical$np >= 30, ]
  if (nrow(empirical) < 4L) stop("Insufficient populated variogram bins.")
  variance <- stats::var(fitting$annual_change)
  distances <- as.numeric(stats::dist(coordinates))
  minimum_scale <- min(distances[distances > 0]) / (10 * cutoff)
  objective <- function(parameters) {
    modeled <- parameters[3] + parameters[1] *
      (1 - exp(-empirical$dist / (parameters[2] * cutoff)))
    sum((modeled - empirical$gamma / variance)^2 * empirical$np / sum(empirical$np))
  }
  fit <- stats::optim(
    c(0.9, 1 / 3, 0.1), objective, method = "L-BFGS-B",
    lower = c(1e-8, minimum_scale, 0), upper = c(10, 10, 10),
    control = list(maxit = 2000, factr = 1e7)
  )
  if (fit$convergence != 0 || any(!is.finite(fit$par)) || !is.finite(fit$value)) {
    stop("Bounded exponential annual-change variogram failed to converge: ", fit$message)
  }
  model <- gstat::vgm(fit$par[1] * variance, "Exp", fit$par[2] * cutoff, fit$par[3] * variance)
  list(model = model, fitting_sites = nrow(fitting), fitting_lags = nrow(empirical),
       normalized_fit_cost = fit$value)
}

run_sequestration <- function(root_dir) {
  years <- c(2005L, 2010L, 2015L, 2020L, 2025L)
  prices <- c(21.40, 21.40, 36, 42, 46)
  input_dir <- file.path(root_dir, "data", "forest-layer")
  output_dir <- file.path(root_dir, "data", "GEP")
  quantity_dir <- file.path(root_dir, "data", "forest-sequestration")
  map_dir <- file.path(root_dir, "outputs", "maps")
  figure_dir <- file.path(root_dir, "outputs", "figures")
  table_dir <- file.path(root_dir, "outputs", "tables")
  template_paths <- file.path(input_dir, sprintf("forest-carbon-stock-%d.tif", years))
  metadata_paths <- sub("\\.tif$", ".json", template_paths)
  quantity_paths <- file.path(quantity_dir, sprintf("forest-carbon-sequestration-%d.tif", years))
  value_paths <- file.path(output_dir, sprintf("gep-carbon-sequestration-%d.tif", years))
  sample_paths <- file.path(quantity_dir, sprintf("matched-fia-carbon-%d.csv", years))
  price_path <- file.path(output_dir, "carbon-sequestration-scc-prices.csv")
  summary_path <- file.path(table_dir, "carbon-sequestration-time-series.csv")
  coverage_path <- file.path(table_dir, "carbon-sequestration-year-coverage.csv")
  exclusions_path <- file.path(table_dir, "carbon-sequestration-pair-exclusions.csv")
  boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
  inputs <- c(template_paths, metadata_paths, boundary_path,
              file.path(root_dir, "data", "FIA", c("MN_PLOT.csv", "MN_COND.csv", "MN_TREE.csv")))
  if (any(!file.exists(inputs))) stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  outputs <- c(quantity_paths, sub("\\.tif$", ".json", quantity_paths), value_paths,
               sub("\\.tif$", ".json", value_paths), sample_paths, price_path, summary_path,
               coverage_path, exclusions_path)
  if (any(file.exists(outputs))) {
    stop("Refusing to overwrite sequestration outputs: ",
         paste(outputs[file.exists(outputs)], collapse = ", "))
  }
  metadata <- lapply(metadata_paths, jsonlite::fromJSON)
  templates <- lapply(seq_along(years), function(index) {
    raster <- terra::rast(template_paths[index])
    info <- metadata[[index]]
    if (terra::nlyr(raster) != 2L || !terra::same.crs(raster, "EPSG:26915") ||
        any(terra::res(raster) != 300) || info$target_year != years[index] ||
        info$attribute != "Live aboveground tree carbon" ||
        info$units != "metric tonnes C per hectare" ||
        info$mask_year != if (years[index] == 2025L) 2022L else years[index]) {
      stop("Unexpected forest template or metadata: ", template_paths[index])
    }
    raster[[1]]
  })
  for (index in seq_along(years)[-1]) terra::compareGeom(templates[[1]], templates[[index]], stopOnError = TRUE)
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) stop("Invalid Minnesota boundary.")
  boundary <- sf::st_sf(geometry = sf::st_union(sf::st_transform(boundary, 26915)))
  pairs <- read_matched_carbon(file.path(root_dir, "data", "FIA"))
  site_keys <- c("STATECD", "UNITCD", "COUNTYCD", "PLOT")
  prepared <- lapply(seq_along(years), function(index) {
    year <- years[index]
    selected <- data.table::copy(pairs[MEASYEAR >= year - 4 & MEASYEAR <= year])
    data.table::setorderv(selected, c(site_keys, "MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "CN"))
    selected <- selected[, .SD[.N], by = site_keys]
    exclusions <- selected[exclusion != "", .N, by = exclusion]
    selected <- selected[exclusion == ""]
    if (nrow(selected) < 32L) {
      warning(year, ": unavailable; fewer than 32 valid matched forest sites. No rate was fabricated.")
      return(list(available = FALSE, exclusions = exclusions, retained_pairs = nrow(selected)))
    }
    sites <- sf::st_transform(sf::st_as_sf(
      as.data.frame(selected), coords = c("LON", "LAT"), crs = 4326, remove = FALSE
    ), 26915)
    inside <- lengths(sf::st_intersects(sites, boundary)) > 0L
    outside <- sum(!inside)
    selected <- selected[inside]
    coordinates <- sf::st_coordinates(sites[inside, ])
    selected[, `:=`(x = coordinates[, 1], y = coordinates[, 2])]
    unique_sites <- selected[, .(annual_change = mean(annual_change)), by = .(x, y)]
    sites <- sf::st_as_sf(as.data.frame(unique_sites), coords = c("x", "y"), crs = 26915)
    if (nrow(sites) < 32L ||
        qr(scale(sf::st_coordinates(sites), center = TRUE, scale = FALSE))$rank < 2) {
      stop(year, ": insufficient noncollinear matched sites.")
    }
    fitted <- fit_change_variogram(sites)
    message(year, ": ", nrow(sites), " matched sites; ",
            sum(selected$annual_change < 0), " observed live-carbon losses.")
    list(available = TRUE, samples = selected, sites = sites, fitted = fitted,
         exclusions = exclusions, outside = outside)
  })
  available <- which(vapply(prepared, function(item) item$available, logical(1)))
  if (!length(available)) stop("No target years have sufficient matched-plot support.")
  for (directory in c(output_dir, quantity_dir, map_dir, figure_dir, table_dir)) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create directory: ", directory)
  }
  warning("Annualized live-aboveground stock changes retain losses but are not a full atmospheric ",
          "CO2 balance. Fully forested matched plots are interpolated statewide without FIA ",
          "population weights or family-ownership attribution; validate before GEP inference.")
  price_table <- data.frame(
    target_year = years, scc_2007_USD_per_tonne_CO2 = prices,
    scenario = "3% discount rate; average; endpoint price applied to annualized interval rate",
    source = c(rep("https://www.epa.gov/sites/default/files/2016-12/documents/scc_tsd_2010.pdf", 2),
               rep("https://19january2017snapshot.epa.gov/climatechange/social-cost-carbon_.html", 3)),
    assumption = c("User-specified 2010 price applied to 2005", "User-specified 2010 price",
                   rep("EPA emission-year-specific price", 3))
  )
  utils::write.csv(price_table, price_path, row.names = FALSE)
  summaries <- vector("list", length(years))
  values_by_year <- vector("list", length(years))
  for (index in seq_along(years)) {
    year <- years[index]
    support <- prepared[[index]]
    if (!support$available) next
    template <- templates[[index]]
    cells <- which(!is.na(terra::values(template, mat = FALSE)))
    if (length(cells) != metadata[[index]]$forest_pixels) stop("Template forest count mismatch: ", year)
    xy <- terra::xyFromCell(template, cells)
    estimates <- uncertainty <- numeric(length(cells))
    for (start in seq.int(1L, length(cells), by = 20000L)) {
      rows <- start:min(start + 19999L, length(cells))
      queries <- sf::st_as_sf(data.frame(x = xy[rows, 1], y = xy[rows, 2]),
                              coords = c("x", "y"), crs = 26915)
      prediction <- gstat::krige(
        annual_change ~ 1, support$sites, queries, model = support$fitted$model,
        nmax = 32, debug.level = 0
      )
      estimates[rows] <- prediction$var1.pred
      uncertainty[rows] <- prediction$var1.var
    }
    if (any(!is.finite(estimates)) || any(!is.finite(uncertainty)) ||
        any(uncertainty < -1e-8)) stop("Nonfinite kriging results or negative variance: ", year)
    quantity <- terra::rast(template)
    variance <- terra::rast(template)
    terra::values(quantity) <- NA_real_
    terra::values(variance) <- NA_real_
    quantity[cells] <- estimates
    variance[cells] <- pmax(uncertainty, 0)
    names(quantity) <- "annual_live_AG_carbon_change_tonnes_C_per_ha_per_year"
    names(variance) <- "kriging_prediction_variance"
    value <- quantity * (44 / 12) * prices[index]
    names(value) <- "annual_live_AG_carbon_value_2007_USD_per_ha_per_year"
    for (item in list(list(raster = c(quantity, variance), path = quantity_paths[index]),
                      list(raster = value, path = value_paths[index]))) {
      terra::writeRaster(item$raster, item$path, overwrite = FALSE, datatype = "FLT8S",
                          NAflag = -9999, gdal = c("COMPRESS=DEFLATE", "PREDICTOR=3", "TILED=YES"))
      saved <- terra::rast(item$path)
      terra::compareGeom(item$raster, saved, stopOnError = TRUE)
      for (band in seq_len(terra::nlyr(saved))) {
        actual <- terra::values(saved[[band]], mat = FALSE)
        expected <- terra::values(item$raster[[band]], mat = FALSE)
        if (!identical(is.na(actual), is.na(expected)) ||
            any(abs(actual[!is.na(actual)] - expected[!is.na(expected)]) > 1e-8)) {
          stop("Written quantity/value raster failed verification: ", item$path)
        }
      }
    }
    info <- list(
      target_year = year, endpoint_measurement_window = c(year - 4, year),
      previous_measurement_year_range = range(support$samples$previous_MEASYEAR),
      remeasurement_interval_years = range(support$samples$REMPER),
      quantity = "Signed annualized change in live aboveground tree carbon",
      quantity_units = "metric tonnes C per hectare per year",
      value_units = "2007 USD per hectare per year",
      formula = "(current_stock_C - previous_stock_C) / FIA_PLOT.REMPER * (44/12) * SCC",
      matching = "PREV_PLT_CN; identical physical keys and DESIGNCD; sampled at both visits",
      forest_rule = "All conditions forest at both visits; forest shares within 1e-6 of one",
      selection = "Newest endpoint per physical plot in window, before exclusions",
      retained_pair_count = nrow(support$samples), unique_sites = nrow(support$sites),
      excluded_pairs = as.data.frame(support$exclusions), outside_boundary = support$outside,
      observed_negative_pairs = sum(support$samples$annual_change < 0),
      variogram = as.data.frame(support$fitted$model),
      fitting_sites = support$fitted$fitting_sites, fitting_lags = support$fitted$fitting_lags,
      variogram_fitting = "Bounded exponential; pair-count weighted normalized least squares; seed 42",
      normalized_fit_cost = support$fitted$normalized_fit_cost,
      kriging = "Local ordinary kriging; 32 nearest sites; signed predictions; point prediction variance",
      mask_template = template_paths[index], esa_mask_year = metadata[[index]]$mask_year,
      scc = price_table[index, ],
      interpretation = paste(
        "Annualized live-pool change proxy, not full net ecosystem/atmospheric sequestration.",
        "Mortality/removals leave the live pool but may enter deadwood/products rather than atmosphere.",
        "Excludes other pools, land-use transitions, mixed-forest plots, and ownership attribution.",
        "Fully forested matched-plot selection is extrapolated to all masked forest pixels.",
        "No design-based expansion or spatial cross-validation; coordinate and temporal uncertainty omitted.",
        "REMPER may exceed five years; endpoint composites are not observations of a single annual flow."
      )
    )
    for (path in c(sub("\\.tif$", ".json", quantity_paths[index]),
                    sub("\\.tif$", ".json", value_paths[index]))) {
      jsonlite::write_json(info, path, auto_unbox = TRUE, pretty = TRUE, digits = NA)
    }
    data.table::fwrite(support$samples, sample_paths[index])
    area <- prod(terra::res(template)) / 10000
    multiplier <- (44 / 12) * prices[index]
    summaries[[index]] <- data.frame(
      target_year = year, endpoint_start_year = year - 4, esa_mask_year = metadata[[index]]$mask_year,
      matched_pairs = nrow(support$samples), unique_sites = nrow(support$sites),
      minimum_REMPER = min(support$samples$REMPER), maximum_REMPER = max(support$samples$REMPER),
      forest_cells = length(cells), mapped_area_ha = length(cells) * area,
      mean_annual_change_tonnes_C_per_ha = mean(estimates),
      total_annual_change_tonnes_C = sum(estimates) * area,
      total_annual_CO2_equivalent_tonnes = sum(estimates) * area * (44 / 12),
      scc_2007_USD_per_tonne_CO2 = prices[index],
      mean_annual_value_2007_USD_per_ha = mean(estimates) * multiplier,
      total_annual_value_2007_USD = sum(estimates) * area * multiplier,
      total_annual_value_fixed_2005_SCC_2007_USD = sum(estimates) * area * (44 / 12) * prices[1],
      positive_change_tonnes_C_per_year = sum(pmax(estimates, 0)) * area,
      negative_change_tonnes_C_per_year = sum(pmin(estimates, 0)) * area,
      scope = "Statewide interpolation of fully forested matched plots; live-pool change proxy"
    )
    values_by_year[[index]] <- value
    message("Saved ", year, " annual-change quantity and value rasters; negative values retained.")
  }
  summary <- do.call(rbind, summaries)
  utils::write.csv(summary, summary_path, row.names = FALSE)
  coverage <- data.frame(
    target_year = years, available = seq_along(years) %in% available,
    reason = vapply(prepared, function(item) {
      if (item$available) "Supported matched-plot estimate" else "Fewer than 32 valid matched sites"
    }, character(1))
  )
  utils::write.csv(coverage, coverage_path, row.names = FALSE)
  exclusion_table <- data.table::rbindlist(lapply(seq_along(years), function(index) {
    data <- data.table::copy(prepared[[index]]$exclusions)
    data[, target_year := years[index]]
    data
  }))
  data.table::fwrite(exclusion_table, exclusions_path)
  visualize_sequestration(values_by_year[available], summary, boundary, map_dir, figure_dir)
  invisible(summary)
}

visualize_sequestration <- function(values, summary, boundary, map_dir, figure_dir) {
  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  first <- values[[1]]
  last <- values[[length(values)]]
  first_year <- summary$target_year[1]
  last_year <- tail(summary$target_year, 1)
  # Absolute change remains interpretable when the baseline rate is zero or negative.
  change <- last - first
  maps <- list(
    "carbon-sequestration-value-2025" = list(
      raster = last, title = "A. Annualized live-carbon value: 2025",
      legend = "Annualized value\n2007 USD / ha / year",
      note = "Matched FIA visits; ESA mask: 2022; SCC: $46/t CO2.\nSigned live-pool change proxy, not a complete atmospheric balance."
    ),
    "carbon-sequestration-value-change" = list(
      raster = change, title = paste0("B. Annualized value difference: ", first_year, "-", last_year),
      legend = "Value difference\n2007 USD / ha / year",
      note = "Shared forest pixels only; absolute difference, not percentage.\nIncludes changing SCC prices and matched-plot sampling support."
    )
  )
  names(maps)[2] <- paste0("carbon-sequestration-value-change-", first_year, "-", last_year)
  names(maps)[1] <- paste0("carbon-sequestration-value-", last_year)
  maps[[1]]$title <- paste0("A. Annualized live-carbon value: ", last_year)
  maps[[1]]$note <- paste0(
    "Matched FIA visits; ESA mask: ", tail(summary$esa_mask_year, 1),
    "; SCC: $", tail(summary$scc_2007_USD_per_tonne_CO2, 1),
    "/t CO2.\nSigned live-pool change proxy, not a complete atmospheric balance."
  )
  for (name in names(maps)) {
    item <- maps[[name]]
    vals <- terra::values(item$raster, mat = FALSE)
    bound <- max(abs(vals), na.rm = TRUE)
    if (!is.finite(bound) || bound <= 0) stop("No nonzero values to visualize: ", name)
    map <- tmap::tm_shape(boundary, crs = 26915, is.main = TRUE) +
      tmap::tm_polygons(fill = "#eeeeee", col = NA) +
      tmap::tm_shape(item$raster) +
      tmap::tm_raster(
        col.scale = tmap::tm_scale_continuous(
          limits = c(-bound, bound), midpoint = 0,
          values = c("#b2182b", "#f7f7f7", "#2166ac"), value.na = NA
        ),
        col.legend = tmap::tm_legend(title = item$legend, na.show = FALSE)
      ) +
      tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.7) +
      tmap::tm_title(item$title, size = 0.85, position = tmap::tm_pos_out("center", "top")) +
      tmap::tm_credits(item$note, size = 0.65, position = tmap::tm_pos_out("left", "bottom")) +
      tmap::tm_layout(frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
                      legend.position = tmap::tm_pos_out("right", "center"),
                      legend.title.size = 0.8, legend.text.size = 0.7)
    for (extension in c("png", "pdf")) {
                        path <- file.path(map_dir, paste0(name, ".", extension))
                        tmap::tmap_save(map, path,
                                        width = 7, height = 7, units = "in", dpi = 300)
                        if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
                          stop("Map export failed: ", path)
                        }
    }
  }
  for (extension in c("png", "pdf")) {
    path <- file.path(figure_dir, paste0("carbon-sequestration-time-series.", extension))
    if (extension == "png") grDevices::png(path, width = 8, height = 8, units = "in", res = 300, type = "cairo")
    else grDevices::pdf(path, width = 8, height = 8)
    tryCatch({
      graphics::par(mfrow = c(2, 1), mar = c(4.5, 5, 2.5, 1), oma = c(4, 0, 2, 0), las = 1)
      graphics::plot(summary$target_year, summary$total_annual_change_tonnes_C / 1e6,
                      type = "o", pch = 16, col = "#238443", lwd = 2, xaxt = "n",
                      xlab = "", ylab = "Million tonnes C / year",
                      main = "Annualized net live-aboveground carbon change")
      graphics::axis(1, at = summary$target_year)
      graphics::abline(h = 0, col = "#aaaaaa", lty = 3)
      annual <- summary$total_annual_value_2007_USD / 1e9
      fixed <- summary$total_annual_value_fixed_2005_SCC_2007_USD / 1e9
      graphics::plot(summary$target_year, annual, type = "o", pch = 16, col = "#2166ac",
                      lwd = 2, xaxt = "n", ylim = range(c(annual, fixed)),
                      xlab = "Endpoint composite target year", ylab = "Billion 2007 USD / year",
                      main = "SCC-valued annualized live-pool change")
      graphics::axis(1, at = summary$target_year)
      graphics::lines(summary$target_year, fixed, type = "o", pch = 1, lty = 2, col = "#666666")
      graphics::legend("topleft", c("Endpoint SCC", "Fixed $21.40/t CO2"),
                        lty = c(1, 2), col = c("#2166ac", "#666666"), bty = "n", cex = 0.8)
      graphics::mtext("C. Exploratory statewide annual-change valuation", side = 3, outer = TRUE, font = 2)
      graphics::mtext(
        "Matched-plot intervals may exceed five years; each endpoint composite uses its own forest mask.\nLive-pool change, not full atmospheric sequestration or family-forest-only GEP.",
        side = 1, outer = TRUE, line = 0.5, cex = 0.75
      )
    }, finally = grDevices::dev.off())
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Time-series export failed: ", path)
    }
  }
}

run_sequestration(root_dir)
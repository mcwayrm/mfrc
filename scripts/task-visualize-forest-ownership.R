# Map the circa-2022 forest ownership dataset and summarize Minnesota forest area.
# Source: Harris, Caputo, and Butler (2025), https://doi.org/10.2737/RDS-2025-0045.
# Class 4 includes corporate and other private owners, not just forest industry.

required_packages <- c("terra", "sf", "tmap", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
if (utils::packageVersion("tmap") < "4.0") stop("Ownership maps require tmap >= 4.0.")

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

create_ownership_map <- function(root_dir) {
  source_path <- file.path(root_dir, "data", "david-wilson", "NationalWoodlandOwnerSurvey",
                           "Data", "US_forest_ownership.tif")
  boundary_path <- file.path(root_dir, "data", "david-wilson", "GIS", "bdry_counpy2.shp")
  raster_path <- file.path(root_dir, "data", "forest-ownership", "forest-ownership-2022-100m.tif")
  metadata_path <- sub("\\.tif$", ".json", raster_path)
  table_path <- file.path(root_dir, "outputs", "tables", "forest-ownership-summary.csv")
  map_paths <- file.path(root_dir, "outputs", "maps", paste0("forest-ownership-2022.", c("png", "pdf")))
  inputs <- c(source_path, boundary_path)
  if (any(!file.exists(inputs))) {
    stop("Missing inputs: ", paste(inputs[!file.exists(inputs)], collapse = ", "))
  }
  outputs <- c(raster_path, metadata_path, table_path, map_paths)
  if (any(file.exists(outputs))) {
    stop("Refusing to overwrite ownership outputs: ",
         paste(outputs[file.exists(outputs)], collapse = ", "))
  }
  for (directory in unique(dirname(outputs))) {
    dir.create(directory, recursive = TRUE, showWarnings = FALSE)
    if (!dir.exists(directory)) stop("Could not create directory: ", directory)
  }

  source <- terra::rast(source_path)
  if (terra::nlyr(source) != 1L || !terra::same.crs(source, "EPSG:4269")) {
    stop("Expected a single-band ownership raster in NAD83: ", source_path)
  }
  boundary <- sf::st_read(boundary_path, quiet = TRUE)
  if (is.na(sf::st_crs(boundary)) || any(!sf::st_is_valid(boundary))) {
    stop("Minnesota boundary must have a CRS and valid geometry.")
  }
  boundary <- sf::st_sf(geometry = sf::st_union(sf::st_transform(boundary, 5070)))
  bbox <- sf::st_bbox(boundary)
  resolution_m <- 100
  grid <- terra::rast(
    xmin = floor(bbox["xmin"] / resolution_m) * resolution_m,
    xmax = ceiling(bbox["xmax"] / resolution_m) * resolution_m,
    ymin = floor(bbox["ymin"] / resolution_m) * resolution_m,
    ymax = ceiling(bbox["ymax"] / resolution_m) * resolution_m,
    resolution = resolution_m, crs = "EPSG:5070"
  )
  temp_path <- tempfile(fileext = ".tif")
  on.exit(unlink(temp_path), add = TRUE)
  message("Reprojecting categorical ownership to a 100 m equal-area Minnesota grid.")
  projected <- terra::project(
    source, grid, method = "near", filename = temp_path,
    wopt = list(datatype = "INT1U", NAflag = 255,
                gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  )
  ownership <- terra::mask(
    projected, terra::vect(boundary), touches = FALSE, filename = raster_path,
    wopt = list(datatype = "INT1U", NAflag = 255,
                gdal = c("COMPRESS=DEFLATE", "TILED=YES"))
  )
  terra::compareGeom(ownership, grid, stopOnError = TRUE)
  frequencies <- terra::freq(ownership)
  if (is.null(frequencies) || any(!frequencies$value %in% 0:8)) {
    stop("The Minnesota ownership raster contains no data or unexpected class codes.")
  }
  counts <- numeric(9)
  counts[frequencies$value + 1L] <- frequencies$count
  state_mask <- terra::rasterize(terra::vect(boundary), grid, field = 1, touches = FALSE)
  state_cells <- terra::global(state_mask, "sum", na.rm = TRUE)[1, 1]
  missing_cells <- state_cells - sum(counts)
  if (!is.finite(missing_cells) || missing_cells < 0) stop("Invalid boundary coverage count.")
  if (missing_cells > 0) {
    warning(missing_cells, " Minnesota grid cells have no ownership source data.")
  }
  forest_codes <- c(0L, 3L:8L)
  forest_cells <- sum(counts[forest_codes + 1L])
  if (forest_cells <= 0 || any(counts[4:9] == 0)) {
    stop("The ownership raster has no forest or is missing an expected ownership class.")
  }
  labels <- c("Unknown forest", "Non-forest", "Water", "Family forest",
              "Corporate/other private forest", "Tribal forest",
              "Federal forest", "State forest", "Local forest")
  forest_counts <- counts[forest_codes + 1L]
  summary <- data.frame(
    class_code = forest_codes,
    ownership_class = labels[forest_codes + 1L],
    cells = forest_counts,
    area_hectares = forest_counts * resolution_m^2 / 10000,
    percent_of_all_forest = 100 * forest_counts / forest_cells
  )
  utils::write.csv(summary, table_path, row.names = FALSE)
  jsonlite::write_json(list(
    source = "Harris, Caputo, and Butler (2025), https://doi.org/10.2737/RDS-2025-0045",
    source_raster = source_path,
    source_period = "Forest ownership circa 2022; forest cover based on 2016 canopy",
    boundary = boundary_path,
    crs = "EPSG:5070",
    resolution_m = resolution_m,
    resampling = "Nearest neighbor, with Minnesota boundary tested at pixel center",
    class_codes = as.list(stats::setNames(labels, 0:8)),
    forest_codes = forest_codes,
    forest_hectares = sum(summary$area_hectares),
    state_cells_without_source = missing_cells,
    interpretation = paste(
      "Area shares are approximate equal-area grid-cell counts, including unknown",
      "forest in the denominator; they are not design-based FIA estimates.",
      "Corporate/other private includes organizations other than forest industry.",
      "The 100 m sample does not retain all 10 m source detail. Ownership is",
      "modeled and may not coincide with ESA forest or current ownership."
    )
  ), metadata_path, auto_unbox = TRUE, pretty = TRUE, digits = NA)

  display <- terra::as.factor(ownership)
  levels(display) <- data.frame(value = 0:8, ownership = labels)
  old_mode <- tmap::tmap_mode("plot")
  on.exit(tmap::tmap_mode(old_mode), add = TRUE)
  map <- tmap::tm_shape(boundary, crs = 5070, is.main = TRUE) +
    tmap::tm_polygons(fill = "white", col = NA) +
    tmap::tm_shape(display) +
    tmap::tm_raster(
      col.scale = tmap::tm_scale_categorical(
        values = c("#8c8c8c", "#eeeeee", "#d9eff5", "#2a9d55", "#d98016",
                   "#9253a1", "#306bb1", "#c5a500", "#b93333"),
        value.na = NA
      ),
      col.legend = tmap::tm_legend(title = "Ownership / land cover", na.show = FALSE)
    ) +
    tmap::tm_shape(boundary) + tmap::tm_borders(col = "#555555", lwd = 0.6) +
    tmap::tm_title("Minnesota forest ownership, circa 2022", size = 0.95,
                   position = tmap::tm_pos_out("center", "top")) +
    tmap::tm_credits(
      paste("Harris, Caputo and Butler (2025); ownership circa 2022, forest cover based on 2016 canopy.",
            "100 m nearest-neighbor grid; classifications are modeled, not parcel boundaries.", sep = "\n"),
      size = 0.65, position = tmap::tm_pos_out("left", "bottom")
    ) +
    tmap::tm_layout(frame = FALSE, meta.margins = c(0.12, 0, 0.07, 0.25),
                    legend.position = tmap::tm_pos_out("right", "center"),
                    legend.title.size = 0.8, legend.text.size = 0.7)
  for (path in map_paths) {
    tmap::tmap_save(map, path, width = 8, height = 7, units = "in", dpi = 300)
    if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size == 0) {
      stop("Map export failed: ", path)
    }
  }
  print(summary)
  message("Saved ownership raster, area summary, and map.")
  invisible(summary)
}

create_ownership_map(root_dir)

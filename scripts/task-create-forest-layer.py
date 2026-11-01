"""Create exploratory Minnesota basal-area surfaces; run only after review."""

from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path

import fiona
import numpy as np
import pandas as pd
import rasterio
from numpy.typing import NDArray
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.features import geometry_mask
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.spatial.distance import pdist
from shapely import contains_xy
from shapely.geometry import mapping, shape
from shapely.geometry.base import BaseGeometry
from shapely.ops import unary_union
from rasterio.warp import transform_geom

ROOT = Path(__file__).resolve().parents[1]
YEARS = (2005, 2010, 2015, 2020, 2025)
MASK_YEARS = {2005: 2005, 2010: 2010, 2015: 2015, 2020: 2020, 2025: 2022}
FOREST_CLASSES = (50, 60, 61, 62, 70, 71, 72, 80, 81, 82, 90, 100)
CRS = "EPSG:26915"
RESOLUTION = 300
NODATA = -9999.0
NEIGHBORS = 32
PREDICTION_BATCH = 256
LOG = logging.getLogger(__name__)


def read_fia(fia_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    plot_columns = [
        "CN", "STATECD", "UNITCD", "COUNTYCD", "PLOT", "MEASYEAR", "MEASMON", "MEASDAY",
        "INVYR", "PLOT_STATUS_CD", "LAT", "LON",
    ]
    condition_columns = [
        "PLT_CN", "CONDID", "COND_STATUS_CD", "CONDPROP_UNADJ", "BALIVE",
    ]
    plots = pd.read_csv(
        fia_dir / "MN_PLOT.csv", usecols=plot_columns, dtype=str
    )
    conditions = pd.read_csv(
        fia_dir / "MN_COND.csv", usecols=condition_columns, dtype=str
    )
    if plots["CN"].isna().any() or plots["CN"].duplicated().any():
        raise ValueError("PLOT.CN must be present and unique.")
    if conditions[["PLT_CN", "CONDID"]].isna().any().any():
        raise ValueError("Condition plot and condition identifiers must be present.")
    if conditions.duplicated(["PLT_CN", "CONDID"]).any():
        raise ValueError("Duplicate FIA plot-condition records.")
    for column in ("MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "PLOT_STATUS_CD", "LAT", "LON"):
        plots[column] = pd.to_numeric(plots[column], errors="raise")
    for column in ("COND_STATUS_CD", "CONDPROP_UNADJ", "BALIVE"):
        conditions[column] = pd.to_numeric(conditions[column], errors="raise")
    if not plots["STATECD"].eq("27").all():
        raise ValueError("The plot database must contain Minnesota (STATECD=27) only.")
    return plots, conditions


def samples_for_year(
    plots: pd.DataFrame,
    conditions: pd.DataFrame,
    year: int,
    boundary: BaseGeometry,
) -> tuple[NDArray[np.float64], NDArray[np.float64], dict[str, int]]:
    selected = plots.loc[plots["MEASYEAR"].between(year - 4, year)].copy()
    site_columns = ["STATECD", "UNITCD", "COUNTYCD", "PLOT"]
    if selected[site_columns + ["INVYR"]].isna().any().any():
        raise ValueError(f"{year}: missing physical plot identifiers or inventory year.")
    for column, maximum in (("MEASMON", 12), ("MEASDAY", 31)):
        unknown = selected[column].isna() | selected[column].eq(0)
        LOG.warning(
            "%s: %s visits have unknown %s; order these before known dates",
            year, int(unknown.sum()), column,
        )
        selected[column] = selected[column].fillna(0)
        if (
            ~selected[column].between(0, maximum)
            | selected[column].mod(1).ne(0)
        ).any():
            raise ValueError(f"{year}: invalid measurement date component {column}.")
    # CN distinguishes visits; the physical plot key prevents repeated-site weighting.
    selected = selected.sort_values(
        ["MEASYEAR", "MEASMON", "MEASDAY", "INVYR", "CN"]
    ).drop_duplicates(
        site_columns, keep="last"
    )
    visits = len(selected)
    selected = selected.loc[selected["PLOT_STATUS_CD"].eq(1)]
    forest = conditions.loc[
        conditions["PLT_CN"].isin(selected["CN"])
        & conditions["COND_STATUS_CD"].eq(1)
    ].copy()
    invalid = (
        ~np.isfinite(forest["BALIVE"])
        | forest["BALIVE"].lt(0)
        | ~np.isfinite(forest["CONDPROP_UNADJ"])
        | ~forest["CONDPROP_UNADJ"].between(0, 1, inclusive="right")
    )
    bad_plots = forest.loc[invalid, "PLT_CN"].unique()
    LOG.warning(
        "%s: excluding %s plots with missing/invalid forest basal area or proportions",
        year, len(bad_plots),
    )
    forest = forest.loc[~forest["PLT_CN"].isin(bad_plots)].copy()
    forest["weighted_ba"] = forest["BALIVE"] * forest["CONDPROP_UNADJ"]
    aggregated = forest.groupby("PLT_CN").agg(
        weighted_ba=("weighted_ba", "sum"),
        forest_share=("CONDPROP_UNADJ", "sum"),
    )
    if aggregated["forest_share"].gt(1.001).any():
        raise ValueError(f"{year}: forest condition proportions exceed one.")
    aggregated["basal_area"] = aggregated["weighted_ba"] / aggregated["forest_share"]
    selected = selected.merge(
        aggregated[["basal_area"]], left_on="CN", right_index=True, validate="one_to_one"
    )
    valid_coordinates = (
        np.isfinite(selected["LAT"]) & np.isfinite(selected["LON"])
        & selected["LAT"].between(-90, 90)
        & selected["LON"].between(-180, 180)
    )
    LOG.warning("%s: excluding %s invalid coordinates", year, int((~valid_coordinates).sum()))
    selected = selected.loc[valid_coordinates]
    transformer = Transformer.from_crs("EPSG:4326", CRS, always_xy=True)
    x, y = transformer.transform(
        selected["LON"].to_numpy(), selected["LAT"].to_numpy(), errcheck=True
    )
    selected["x"], selected["y"] = x, y
    inside = contains_xy(boundary, x, y)
    LOG.warning("%s: excluding %s points outside Minnesota", year, int((~inside).sum()))
    selected = selected.loc[inside]
    # Public coordinates can coincide after fuzzing; average rather than invert
    # a singular kriging matrix or treat coincident observations as extra sites.
    sites = selected.groupby(["x", "y"], as_index=False)["basal_area"].mean()
    stats = {
        "newest_visits_in_window": visits,
        "invalid_forest_attribute_plots": len(bad_plots),
        "usable_forest_visits": len(selected),
        "coincident_visits_collapsed": len(selected) - len(sites),
        "unique_sites": len(sites),
    }
    coordinates = sites[["x", "y"]].to_numpy(dtype=np.float64)
    values = sites["basal_area"].to_numpy(dtype=np.float64)
    if len(sites) < NEIGHBORS or not np.isfinite(values).all():
        raise ValueError(f"{year}: need at least {NEIGHBORS} valid, unique forest sites.")
    if np.linalg.matrix_rank(coordinates - coordinates.mean(axis=0)) < 2:
        raise ValueError(f"{year}: sample locations are collinear.")
    LOG.info("%s: using %s forest sites from %s-%s", year, len(sites), year - 4, year)
    return coordinates, values, stats


def fit_variogram(
    coordinates: NDArray[np.float64], values: NDArray[np.float64]
) -> dict[str, float]:
    # Bound pairwise fitting memory; prediction still uses every eligible site.
    indices = np.random.default_rng(42).choice(
        len(values), size=min(1000, len(values)), replace=False
    )
    distance = pdist(coordinates[indices])
    semivariance = 0.5 * pdist(values[indices, None], metric="sqeuclidean")
    edges = np.linspace(0, distance.max(), 21)
    lag, empirical, counts = [], [], []
    for lower, upper in zip(edges[:-1], edges[1:]):
        chosen = (distance > lower) & (distance <= upper)
        if chosen.sum() >= 30:
            lag.append(float(distance[chosen].mean()))
            empirical.append(float(semivariance[chosen].mean()))
            counts.append(int(chosen.sum()))
    if len(lag) < 4 or np.var(values) <= 0:
        raise ValueError("Insufficient spatial variation to fit a variogram.")
    h = np.asarray(lag)
    gamma = np.asarray(empirical)
    weights = np.sqrt(np.asarray(counts) / np.sum(counts))
    variance = float(np.var(values))
    minimum_scale = float(distance[distance > 0].min())
    maximum_scale = float(distance.max())

    def residual(parameters: NDArray[np.float64]) -> NDArray[np.float64]:
        partial_sill, scale, nugget = parameters
        modeled = nugget + partial_sill * (1 - np.exp(-h / scale))
        return (modeled - gamma) * weights / variance

    fitted = least_squares(
        residual,
        x0=[variance * 0.9, maximum_scale / 3, variance * 0.1],
        bounds=(
            [variance * 1e-8, minimum_scale / 10, 0],
            [variance * 10, maximum_scale * 10, variance * 10],
        ),
        x_scale="jac",
    )
    if not fitted.success or not np.isfinite(fitted.x).all():
        raise ValueError(f"Variogram fit failed: {fitted.message}")
    return {
        "partial_sill": float(fitted.x[0]),
        "scale_m": float(fitted.x[1]),
        "nugget": float(fitted.x[2]),
        "fitting_sites": len(indices),
        "fitting_lags": len(lag),
        "normalized_fit_cost": float(fitted.cost),
    }


def predict(
    queries: NDArray[np.float64],
    coordinates: NDArray[np.float64],
    values: NDArray[np.float64],
    tree: cKDTree,
    model: dict[str, float],
) -> tuple[NDArray[np.float64], NDArray[np.float64], float]:
    distances, indices = tree.query(queries, k=NEIGHBORS)
    neighbors = coordinates[indices]
    pair_distance = np.linalg.norm(
        neighbors[:, :, None, :] - neighbors[:, None, :, :], axis=-1
    )
    sill = model["partial_sill"] + model["nugget"]
    covariance = model["partial_sill"] * np.exp(-pair_distance / model["scale_m"])
    diagonal = np.arange(NEIGHBORS)
    covariance[:, diagonal, diagonal] = sill
    cross_covariance = model["partial_sill"] * np.exp(-distances / model["scale_m"])
    matrix = np.ones((len(queries), NEIGHBORS + 1, NEIGHBORS + 1))
    matrix[:, :-1, :-1] = covariance
    matrix[:, -1, -1] = 0
    rhs = np.ones((len(queries), NEIGHBORS + 1))
    rhs[:, :-1] = cross_covariance
    solution = np.linalg.solve(matrix, rhs[..., None])[..., 0]
    estimates = np.sum(solution[:, :-1] * values[indices], axis=1)
    # Predict latent basal area: nugget represents noisy/fuzzed observations,
    # not extra variance in the underlying spatial surface.
    variance = (
        model["partial_sill"]
        - np.sum(solution[:, :-1] * cross_covariance, axis=1)
        - solution[:, -1]
    )
    if not np.isfinite(estimates).all() or not np.isfinite(variance).all():
        raise ValueError("Kriging returned nonfinite predictions or variances.")
    if (variance < -1e-6 * sill).any():
        raise ValueError("Kriging returned materially negative variances.")
    return estimates, np.maximum(variance, 0), float(distances[:, 0].max())


def create_layer(
    year: int, mask_path: Path, output_dir: Path,
    plots: pd.DataFrame, conditions: pd.DataFrame, boundary: BaseGeometry,
) -> None:
    coordinates, values, sample_stats = samples_for_year(plots, conditions, year, boundary)
    model = fit_variogram(coordinates, values)
    LOG.info("%s: exponential variogram %s", year, model)
    tree = cKDTree(coordinates)
    left, bottom, right, top = boundary.bounds
    left = math.floor(left / RESOLUTION) * RESOLUTION
    top = math.ceil(top / RESOLUTION) * RESOLUTION
    width = math.ceil((right - left) / RESOLUTION)
    height = math.ceil((top - bottom) / RESOLUTION)
    transform = from_origin(left, top, RESOLUTION, RESOLUTION)
    output = output_dir / f"forest-layer-{year}.tif"
    partial = output.with_suffix(".partial.tif")
    if partial.exists():
        raise FileExistsError(f"Remove or review previous incomplete output: {partial}")
    predicted_pixels, clipped_pixels, maximum_distance = 0, 0, 0.0
    profile = {
        "driver": "GTiff", "height": height, "width": width, "count": 2,
        "dtype": "float32", "crs": CRS, "transform": transform, "nodata": NODATA,
        "compress": "deflate", "predictor": 3, "tiled": True,
        "blockxsize": 256, "blockysize": 256, "BIGTIFF": "IF_SAFER",
    }
    try:
        with rasterio.open(mask_path) as source:
            if source.crs is None or source.count != 1:
                raise ValueError(f"ESA mask must be a single-band georeferenced raster: {mask_path}")
            with WarpedVRT(
                source, crs=CRS, transform=transform, width=width, height=height,
                resampling=Resampling.nearest,
                src_nodata=source.nodata if source.nodata is not None else 0,
                nodata=0,
            ) as land_cover, rasterio.open(partial, "w", **profile) as destination:
                destination.set_band_description(1, "Forest basal area (ft2/forest acre)")
                destination.set_band_description(2, "Latent kriging variance ((ft2/forest acre)^2)")
                destination.update_tags(
                    target_year=year, measurement_years=f"{year - 4}-{year}",
                    mask_year=MASK_YEARS[year], mask_source=str(mask_path),
                    forest_classes=",".join(map(str, FOREST_CLASSES)),
                    method="Local ordinary kriging; exponential covariance; 32 neighbors",
                    coordinate_caution="Public FIA coordinates are fuzzed/swapped",
                    negative_prediction_policy="Clipped to zero; count in JSON sidecar",
                )
                for _, window in destination.block_windows(1):
                    classes = land_cover.read(1, window=window, masked=True)
                    block_transform = destination.window_transform(window)
                    in_state = geometry_mask(
                        [mapping(boundary)], out_shape=classes.shape,
                        transform=block_transform, invert=True,
                    )
                    forest = (
                        in_state & ~np.ma.getmaskarray(classes)
                        & np.isin(classes.data, FOREST_CLASSES)
                    )
                    rows, columns = np.nonzero(forest)
                    predictions = np.full(classes.shape, NODATA, dtype=np.float32)
                    variances = np.full(classes.shape, NODATA, dtype=np.float32)
                    for start in range(0, len(rows), PREDICTION_BATCH):
                        r = rows[start:start + PREDICTION_BATCH]
                        c = columns[start:start + PREDICTION_BATCH]
                        x, y = rasterio.transform.xy(block_transform, r, c)
                        query = np.column_stack((x, y))
                        estimates, uncertainty, nearest_distance = predict(
                            query, coordinates, values, tree, model
                        )
                        clipped_pixels += int((estimates < 0).sum())
                        maximum_distance = max(maximum_distance, nearest_distance)
                        predictions[r, c] = np.maximum(estimates, 0).astype(np.float32)
                        variances[r, c] = uncertainty.astype(np.float32)
                    predicted_pixels += len(rows)
                    destination.write(predictions, 1, window=window)
                    destination.write(variances, 2, window=window)
        if predicted_pixels == 0:
            raise ValueError(f"{year}: mask has no valid Minnesota forest pixels.")
        metadata = {
            "target_year": year, "measurement_start_year": year - 4,
            "measurement_end_year": year, "mask_year": MASK_YEARS[year],
            "mask_source": str(mask_path), "forest_classes": FOREST_CLASSES,
            "crs": CRS, "resolution_m": RESOLUTION,
            "attribute": "Forest-condition-weighted COND.BALIVE",
            "units": "ft2 per forest acre",
            "condition_weight": "CONDPROP_UNADJ normalized over forest conditions",
            "sample_statistics": sample_stats, "variogram": model,
            "method": "Local ordinary kriging", "neighbors": NEIGHBORS,
            "forest_pixels": predicted_pixels, "negative_predictions_clipped": clipped_pixels,
            "maximum_nearest_sample_distance_m": maximum_distance,
            "limitations": [
                "Public FIA coordinates are fuzzed/swapped; exploratory use only.",
                "Five-year composites are not exact single-year measurements.",
                "Not design-based rFIA estimates; no population expansion weights.",
                "No spatial cross-validation; review variogram and validate before inference.",
                "All forest pixels predicted, including extrapolation beyond sample coverage.",
                "Band 2 excludes land-cover, temporal, and coordinate uncertainty.",
            ],
        }
        output.with_suffix(".json").write_text(
            json.dumps(metadata, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
        partial.replace(output)
    finally:
        if partial.exists():
            partial.unlink()
    LOG.info("%s: saved %s (%s forest pixels)", year, output, predicted_pixels)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fia-dir", type=Path, default=ROOT / "data" / "FIA")
    parser.add_argument(
        "--esa-dir", type=Path,
        default=Path(r"G:\Shared drives\NatCapTEEMs\Files\base_data\lulc\esa"),
    )
    parser.add_argument(
        "--boundary", type=Path,
        default=ROOT / "data" / "david-wilson" / "GIS" / "bdry_counpy2.shp",
        help="Minnesota polygon boundary; county polygons will be dissolved.",
    )
    parser.add_argument("--output-dir", type=Path, default=ROOT / "data" / "forest-layer")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    masks = {
        year: args.esa_dir / f"lulc_esa_{mask_year}.tif"
        for year, mask_year in MASK_YEARS.items()
    }
    for path in [
        args.fia_dir / "MN_PLOT.csv", args.fia_dir / "MN_COND.csv",
        args.boundary, *masks.values(),
    ]:
        if not path.is_file():
            raise FileNotFoundError(f"Required input does not exist: {path}")
    for year in YEARS:
        for suffix in (".tif", ".json", ".partial.tif"):
            path = args.output_dir / f"forest-layer-{year}{suffix}"
            if path.exists():
                raise FileExistsError(f"Refusing to overwrite existing output: {path}")
    with fiona.open(args.boundary) as source:
        if not source.crs:
            raise ValueError("Minnesota boundary has no CRS.")
        polygons = [
            shape(transform_geom(source.crs, CRS, feature["geometry"]))
            for feature in source if feature["geometry"] is not None
        ]
    boundary = unary_union(polygons)
    if boundary.is_empty or not boundary.is_valid or boundary.geom_type not in (
        "Polygon", "MultiPolygon",
    ):
        raise ValueError("Minnesota boundary must contain valid polygon geometry.")
    plots, conditions = read_fia(args.fia_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    LOG.warning("Exploratory surfaces: public FIA coordinates are fuzzed/swapped.")
    LOG.warning("2025 uses the explicitly agreed 2022 ESA forest mask.")
    for year in YEARS:
        create_layer(year, masks[year], args.output_dir, plots, conditions, boundary)


if __name__ == "__main__":
    main()
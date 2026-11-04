"""Public single-volume API for OCT-derived CVI phenotypes.

The image processing and numerical reductions are delegated unchanged to
``processing.py``. A volume here means five grayscale central B-scans with
five matching, parsed SVG layer segmentations.
"""

from __future__ import annotations

from dataclasses import dataclass
import inspect
import json
from pathlib import Path
from typing import Sequence

import numpy as np
from lxml import etree

from .processing import (
    calculate_cvi,
    calculate_pachy_area,
    segment_oct,
    select_central_oct_slices,
)


# Taken from the submitted main.py, not verified as final cohort-specific values.
DEFAULT_PARAMETERS = {
    "kernel_size_o_pachy": 17,
    "kernel_size_c_pachy": 5,
    "block_size_pachy": 91,
    "cst_sub_mean_pachy": 13,
    "kernel_size_o": 5,
    "block_size": 41,
    "cst_sub_mean": 16,
    "kernel_size_o_small": 5,
    "pachy_only": False,
    "medium_only": False,
}


@dataclass(frozen=True)
class CVIResult:
    """Median-over-B-scan phenotypes from the supplied algorithm."""

    cvi_mean: float
    cvi_spatial_sd: float
    pachy_area_mm2: float | None = None


def resolve_parameters(overrides: dict | None = None) -> dict:
    """Validate settings; only combined vessel segmentation defines CVI."""
    parameters = DEFAULT_PARAMETERS.copy()
    if overrides is not None:
        if not isinstance(overrides, dict):
            raise ValueError("Segmentation parameters must be a JSON object/dict")
        unknown = set(overrides) - set(parameters)
        if unknown:
            raise ValueError(f"Unknown segmentation parameters: {sorted(unknown)}")
        parameters.update(overrides)
    for key in (
        "kernel_size_o_pachy", "kernel_size_c_pachy", "block_size_pachy",
        "kernel_size_o", "block_size", "kernel_size_o_small",
    ):
        value = parameters[key]
        if type(value) is not int or value < 3 or value % 2 == 0:
            raise ValueError(f"{key} must be an odd integer >= 3")
    for key in ("cst_sub_mean_pachy", "cst_sub_mean"):
        if type(parameters[key]) not in (int, float):
            raise ValueError(f"{key} must be numeric")
    if parameters["pachy_only"] is not False or parameters["medium_only"] is not False:
        raise ValueError("Both *_only options must be false for combined CVI")
    if set(parameters) - set(inspect.signature(segment_oct).parameters):
        raise ValueError("Parameters do not match segment_oct signature")
    return parameters


def extract_cvi(
    oct_bscans: Sequence[np.ndarray],
    choroid_segmentations: Sequence[etree._Element],
    *,
    spacings: Sequence[float] | None = None,
    parameters: dict | None = None,
) -> CVIResult:
    """Extract CVI from exactly five *central* grayscale B-scans and SVGs.

    ``oct_bscans`` and ``choroid_segmentations`` must already be aligned and
    ordered. The provided processing algorithm is not changed. The SD is the
    median of the five within-B-scan SDs of column-wise CVI, not a pooled SD.
    If spacings are supplied, indices 2 and 3 must be positive mm/pixel.
    """
    if len(oct_bscans) != 5 or len(choroid_segmentations) != 5:
        raise ValueError("Five aligned central B-scans and five SVGs are required")
    shapes = [np.asarray(scan).shape for scan in oct_bscans]
    if any(np.asarray(scan).ndim != 2 for scan in oct_bscans) or len(set(shapes)) != 1:
        raise ValueError("All five OCT B-scans must be grayscale and have the same shape")
    if spacings is not None:
        if len(spacings) < 4:
            raise ValueError("Spacings need at least four entries")
        if not all(np.isfinite(float(x)) and float(x) > 0 for x in spacings[2:4]):
            raise ValueError("Spacing indices 2 and 3 must be positive finite mm/pixel")

    vessels, masks = segment_oct(
        oct_cube=oct_bscans,
        oct_seg=choroid_segmentations,
        **resolve_parameters(parameters),
    )
    mean, sd = calculate_cvi(vessels, masks)
    if mean is None or sd is None or not np.isfinite([mean, sd]).all():
        raise ValueError("CVI cannot be calculated from these segmentations")
    area = calculate_pachy_area(vessels, spacings=spacings) if spacings is not None else None
    return CVIResult(float(mean), float(sd), None if area is None else float(area))


def load_volume_folders(cube_dir: Path | str, seg_dir: Path | str):
    """Read a vendor-exported OCT folder and corresponding SVG segmentations.

    This adapter is separate from the array-based extraction API. Both folders
    must use consecutive, zero-based four-digit stems, and JPEG input.
    """
    cube_dir, seg_dir = Path(cube_dir), Path(seg_dir)
    if not cube_dir.is_dir() or not seg_dir.is_dir():
        raise ValueError("OCT or segmentation directory does not exist")
    images = sorted(cube_dir.glob("*.jpg"))
    if len(images) < 5:
        raise ValueError(f"Expected at least 5 JPEG B-scans, found {len(images)}")
    mid_start = (len(images) - 5) // 2
    missing = [f"{i:04d}" for i in range(mid_start, mid_start + 5)
               if not ((cube_dir / f"{i:04d}.jpg").is_file() and
                       (seg_dir / f"{i:04d}.svg").is_file())]
    if missing:
        raise ValueError(f"Missing central JPEG/SVG pairs: {', '.join(missing)}")
    bscans, svgs = select_central_oct_slices(cube_dir, seg_dir, nb_slices=5)
    return bscans, svgs


def extract_cvi_from_folders(
    cube_dir: Path | str,
    seg_dir: Path | str,
    *,
    parameters: dict | None = None,
) -> CVIResult:
    """Convenience adapter for the original directory export convention."""
    bscans, svgs = load_volume_folders(cube_dir, seg_dir)
    info_file = Path(cube_dir) / "info.json"
    if not info_file.is_file():
        raise ValueError("Missing OCT metadata info.json (needed for pachy area)")
    info = json.loads(info_file.read_text(encoding="utf-8"))
    if "spacing" not in info:
        raise ValueError("Missing spacing in info.json (needed for pachy area)")
    return extract_cvi(bscans, svgs, spacings=info["spacing"], parameters=parameters)

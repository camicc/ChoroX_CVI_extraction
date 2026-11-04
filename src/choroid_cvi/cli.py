"""Command-line entry point for one OCT study or a manifest of studies."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import inspect
import json
from pathlib import Path
import platform
import sys

import numpy as np

from . import __version__
from .processing import (
    calculate_cvi,
    calculate_pachy_area,
    segment_oct,
    select_central_oct_slices,
)


# These are the parameters in the supplied main.py, not independently
# verified as the final UK Biobank or OphtalmoLaus production settings.
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

OUTPUT_COLUMNS = (
    "sample_id", "n_bscans", "cvi_mean", "cvi_spatial_sd",
    "pachy_area_mm2", "status", "error",
)


def load_parameters(path: Path | None) -> dict:
    params = DEFAULT_PARAMETERS.copy()
    if path is not None:
        overrides = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(overrides, dict):
            raise ValueError("Parameter JSON must contain an object")
        unknown = set(overrides) - set(params)
        if unknown:
            raise ValueError(f"Unknown segmentation parameters: {sorted(unknown)}")
        params.update(overrides)

    for key in (
        "kernel_size_o_pachy", "kernel_size_c_pachy", "block_size_pachy",
        "kernel_size_o", "block_size", "kernel_size_o_small",
    ):
        value = params[key]
        if type(value) is not int or value < 3 or value % 2 == 0:
            raise ValueError(f"{key} must be an odd integer >= 3")
    for key in ("cst_sub_mean_pachy", "cst_sub_mean"):
        if type(params[key]) not in (int, float):
            raise ValueError(f"{key} must be numeric")
    if params["pachy_only"] is not False or params["medium_only"] is not False:
        raise ValueError("Both *_only options must be false for the combined CVI phenotype")

    # Protect against changing the algorithm via a future signature mismatch.
    accepted = set(inspect.signature(segment_oct).parameters)
    if set(params) - accepted:
        raise ValueError("Parameters do not match segment_oct signature")
    return params


def extract_one(sample_id: str, cube_dir: Path, seg_dir: Path, params: dict) -> dict:
    """Preserve the supplied segmentation and summary calculations."""
    row = dict.fromkeys(OUTPUT_COLUMNS)
    row.update(sample_id=sample_id, n_bscans=0, status="error", error="")
    try:
        if not cube_dir.is_dir() or not seg_dir.is_dir():
            raise ValueError("OCT or segmentation directory does not exist")

        images = sorted(cube_dir.glob("*.jpg"))
        if len(images) < 5:
            raise ValueError(f"Expected at least 5 JPEG B-scans, found {len(images)}")
        mid_start = (len(images) - 5) // 2
        expected = [f"{i:04d}" for i in range(mid_start, mid_start + 5)]
        missing = [stem for stem in expected if not (
            (cube_dir / f"{stem}.jpg").is_file()
            and (seg_dir / f"{stem}.svg").is_file()
        )]
        if missing:
            raise ValueError(f"Missing central JPEG/SVG pairs: {', '.join(missing)}")

        octs, svgs = select_central_oct_slices(cube_dir, seg_dir, nb_slices=5)
        row["n_bscans"] = len(octs)
        if len(octs) != 5 or any(np.asarray(scan).ndim != 2 for scan in octs):
            raise ValueError("Expected exactly 5 grayscale OCT B-scans")

        info_file = cube_dir / "info.json"
        if not info_file.is_file():
            raise ValueError("Missing OCT metadata info.json (needed for pachy area)")
        info = json.loads(info_file.read_text(encoding="utf-8"))
        spacings = info.get("spacing")
        if not isinstance(spacings, list) or len(spacings) < 4:
            raise ValueError("info.json requires spacing array of at least four elements")
        if not all(np.isfinite(float(v)) and float(v) > 0 for v in (spacings[2], spacings[3])):
            raise ValueError("Spacing indices 2 and 3 must be positive finite values")

        vessels, masks = segment_oct(oct_cube=octs, oct_seg=svgs, **params)
        cvi_mean, cvi_sd = calculate_cvi(vessels, masks)
        pachy_area = calculate_pachy_area(vessels, spacings=spacings)
        if cvi_mean is None or cvi_sd is None or not np.isfinite([cvi_mean, cvi_sd]).all():
            raise ValueError("CVI is undefined or non-finite")
        row.update(cvi_mean=float(cvi_mean), cvi_spatial_sd=float(cvi_sd),
                   pachy_area_mm2=None if pachy_area is None else float(pachy_area),
                   status="ok")
    except Exception as exc:
        row["error"] = f"{type(exc).__name__}: {exc}"
    return row


def _resolve_path(value: str, base: Path) -> Path:
    path = Path(value)
    return path if path.is_absolute() else base / path


def _manifest_rows(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        expected = {"sample_id", "cube_dir", "seg_dir"}
        if not reader.fieldnames or not expected.issubset(reader.fieldnames):
            raise ValueError("Manifest columns must include sample_id,cube_dir,seg_dir")
        seen = set()
        for record in reader:
            sample_id = record["sample_id"].strip()
            if not sample_id or sample_id in seen:
                raise ValueError(f"Missing or duplicate sample_id: {sample_id!r}")
            seen.add(sample_id)
            yield (sample_id, _resolve_path(record["cube_dir"], path.parent),
                   _resolve_path(record["seg_dir"], path.parent))


def _package_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not installed"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Extract CVI phenotypes from five central OCT B-scans")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--manifest", type=Path, help="CSV with sample_id,cube_dir,seg_dir")
    group.add_argument("--cube-dir", type=Path, help="OCT directory for one sample")
    parser.add_argument("--seg-dir", type=Path, help="SVG segmentation directory for one sample")
    parser.add_argument("--sample-id", help="Pseudonymous ID for one sample")
    parser.add_argument("--parameters", type=Path, help="JSON overrides of vessel segmentation settings")
    parser.add_argument("--output", type=Path, required=True, help="CSV output (excluded from git)")
    args = parser.parse_args(argv)
    if args.cube_dir and (not args.seg_dir or not args.sample_id):
        parser.error("Single-sample mode requires --seg-dir and --sample-id")
    if args.manifest and (args.seg_dir or args.sample_id):
        parser.error("--seg-dir and --sample-id only apply to single-sample mode")
    try:
        params = load_parameters(args.parameters)
        samples = list(_manifest_rows(args.manifest)) if args.manifest else [
            (args.sample_id, args.cube_dir, args.seg_dir)
        ]
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        parser.error(str(exc))

    args.output.parent.mkdir(parents=True, exist_ok=True)
    rows = [extract_one(sid, cube, seg, params) for sid, cube, seg in samples]
    with args.output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    provenance = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "pipeline_version": __version__,
        "python": platform.python_version(),
        "dependencies": {name: _package_version(name) for name in (
            "numpy", "scipy", "opencv-python-headless", "Pillow", "lxml",
            "svglib", "reportlab", "rlPyCairo",
        )},
        "central_bscans": 5,
        "segmentation_parameters": params,
        "n_samples": len(rows),
        "n_success": sum(row["status"] == "ok" for row in rows),
        "phenotype_definition": "median of per-B-scan column-wise CVI mean and SD",
    }
    metadata_path = args.output.with_suffix(".run.json")
    metadata_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {len(rows)} rows to {args.output} ({provenance['n_success']} successful)")
    if provenance["n_success"] != len(rows):
        print("Some samples failed; see status/error columns", file=sys.stderr)
        return 1
    return 0

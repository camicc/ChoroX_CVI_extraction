"""Command-line entry point for one OCT study or a manifest of studies."""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
from importlib.metadata import PackageNotFoundError, version
import json
from pathlib import Path
import platform
import sys

import cv2

from . import __version__
from .api import extract_cvi_from_folders, resolve_parameters

OUTPUT_COLUMNS = (
    "sample_id", "n_bscans", "cvi_mean", "cvi_spatial_sd",
    "pachy_area_mm2", "status", "error",
)


def load_parameters(path: Path | None) -> dict:
    overrides = json.loads(path.read_text(encoding="utf-8")) if path is not None else None
    return resolve_parameters(overrides)


def extract_one(sample_id: str, cube_dir: Path, seg_dir: Path, params: dict) -> dict:
    """Run one volume and record success or an actionable error."""
    row = dict.fromkeys(OUTPUT_COLUMNS)
    row.update(sample_id=sample_id, n_bscans=0, status="error", error="")
    try:
        result = extract_cvi_from_folders(cube_dir, seg_dir, parameters=params)
        row.update(n_bscans=5, cvi_mean=result.cvi_mean,
                   cvi_spatial_sd=result.cvi_spatial_sd,
                   pachy_area_mm2=result.pachy_area_mm2,
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
        "opencv_cv2": cv2.__version__,
        "dependencies": {name: _package_version(name) for name in (
            "numpy", "scipy", "Pillow", "lxml",
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

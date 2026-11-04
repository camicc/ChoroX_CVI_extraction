"""Minimal single-volume API example using artificial B-scans and SVGs."""

import json
from pathlib import Path

from choroid_cvi import extract_cvi
from choroid_cvi.processing import read_img


root = Path(__file__).resolve().parent / "synthetic_case"
bscans = [read_img(root / "oct_volume" / f"{i:04d}.jpg") for i in range(5)]
segmentations = [read_img(root / "svg_segmentation" / f"{i:04d}.svg") for i in range(5)]
info = json.loads((root / "oct_volume" / "info.json").read_text(encoding="utf-8"))

result = extract_cvi(bscans, segmentations, spacings=info["spacing"])
print(f"CVI mean: {result.cvi_mean:.6f}")
print(f"CVI spatial variability: {result.cvi_spatial_sd:.6f}")
print(f"Largest vessel area (mm2): {result.pachy_area_mm2:.6f}")

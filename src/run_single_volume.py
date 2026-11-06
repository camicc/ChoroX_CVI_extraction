"""Extract CVI and save vessel segmentations for one OCT volume.

Usage (Windows):
    python reproduce_one_volume.py --root "V:\\...\\patient_id"

Input directory:
    patient_id/
        oct/info.json
        oct/0000.jpg ... oct/0127.jpg
        segmentation/0000.svg ... segmentation/0127.svg

The five central B-scans are 0062–0066, as in the original workflow.
"""

import argparse
import json
from pathlib import Path
import numpy as np

from processing import (
    read_img,
    extract_cvi,
    segment_oct,
    preprocess_oct,
    seg2mask,
    detect_bm_from_seg,
    calculate_cvi,
    plot_enface_oct,
    save_segmentation_overlays,
)

CENTRAL_INDICES = range(62, 67)
FULL_INDICES = range(128)


def load_slices(root: Path, indices):
    """Load matching OCT JPEGs and segmentation SVGs."""
    bscans = [read_img(root / "oct" / f"{i:04d}.jpg") for i in indices]
    segmentations = [read_img(root / "segmentation" / f"{i:04d}.svg") for i in indices]
    return bscans, segmentations


def display(value):
    return "N/A" if value is None else f"{value:.6f}"


def main(root: Path, *, skip_full_volume: bool = False):
    root = root.expanduser()
    with (root / "oct" / "info.json").open(encoding="utf-8") as file:
        spacings = json.load(file)["spacing"]

    # 1. Read the five original central B-scans.
    bscans, segmentations = load_slices(root, CENTRAL_INDICES)

    # 2. Original CVI calculation (median of per-scan column statistics).
    result = extract_cvi(bscans, segmentations, spacings=spacings)
    print(f"CVI mean: {display(result.cvi_mean)}")
    print(f"CVI spatial variability: {display(result.cvi_spatial_sd)}")
    print(f"Legacy largest-component area (mm²): {result.pachy_area_mm2:.6f}")

    # 3. Segment once more to obtain debug images from the same processing run.
    segmented, cropped_masks, processed, _ = segment_oct(
        bscans, segmentations, return_intermediates=True
    )

    for i, svg, mask in zip(CENTRAL_INDICES, segmentations, cropped_masks):
        raw_mask = seg2mask(svg, layers="CC_CS")
        bm, _ = detect_bm_from_seg(svg)
        print(
            f"Slice {i}: raw mask={raw_mask.mean():.1%}, "
            f"cropped mask={mask.mean():.1%}, BM median={np.median(bm):.1f}"
        )

    # 4. Save choroid/vessel overlays. These do not influence CVI values.
    overlay_dir = root / "segmentation_debug"
    save_segmentation_overlays(
        processed, segmented, cropped_masks, overlay_dir, scan_indices=CENTRAL_INDICES
    )

    # 5. Individual scan statistics (SD is across columns, not across scans).
    print("\nIndividual B-scan CVI:")
    for i, (vessels, mask) in zip(CENTRAL_INDICES, zip(segmented, cropped_masks)):
        mean, sd = calculate_cvi([vessels], [mask])
        print(f"Slice {i}: mean={display(mean)}, std={display(sd)}")

    # 6. Optionally generate the original full-volume en face projection.
    if not skip_full_volume:
        all_bscans, all_svgs = load_slices(root, FULL_INDICES)
        full_vessels, _ = segment_oct(all_bscans, all_svgs)
        enface_path = root / "enface_projection.png"
        plot_enface_oct(full_vessels, root.name, save_path=enface_path, spacings=spacings)

    # 7. Summary.
    print("\nProcessing completed.")
    print("Central B-scans:", len(bscans))
    print("Full-volume B-scans:", 0 if skip_full_volume else len(all_bscans))
    if not skip_full_volume:
        print("Vessel volume shape:", np.asarray(full_vessels).shape)
        print("Vessel values:", np.unique(np.asarray(full_vessels)))
    print("Spacings:", spacings)
    print("Segmentation images:", overlay_dir)
    if not skip_full_volume:
        print("En face projection:", enface_path)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", required=True, type=Path, help="OCT volume directory")
    parser.add_argument(
        "--skip-full-volume", action="store_true",
        help="Run only the five central B-scans (faster for CVI checks)",
    )
    args = parser.parse_args()
    main(args.root, skip_full_volume=args.skip_full_volume)

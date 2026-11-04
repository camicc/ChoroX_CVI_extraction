# OCT-derived choroidal vascularity: phenotype extraction

Code accompanying work on the genetic architecture of choroidal vascularity using UK Biobank macular OCT and external evaluation in OphtalmoLaus. This repository covers **OCT-to-phenotype extraction only**. It does not yet reproduce the GWAS, clumping, heritability, gene-based tests, eQTL/colocalisation, MR, or replication analyses.

> **Pre-submission status.** The algorithm was reorganized from three supplied scripts (`main.py`, `helper.py`, `svg_util.py`). The numerical extraction steps have been retained, but the defaults supplied with `main.py` have **not** been verified as the final production parameters for either cohort. Do not cite this as a complete reproduction of the paper until the checks in [docs/reproducibility_checklist.md](docs/reproducibility_checklist.md) are completed.

## Core API: one OCT volume

The intended public interface is `extract_cvi`: one call takes five central grayscale B-scans and their five matching parsed SVG segmentations. It does not depend on participant identifiers or dataset-specific directory structures.

```python
from choroid_cvi import extract_cvi

result = extract_cvi(
    oct_bscans,                 # five 2-D grayscale NumPy arrays
    choroid_segmentations,     # five corresponding SVG XML elements
    spacings=[1, 1, 0.01, 0.01],  # optional; last two are mm/pixel
)
print(result.cvi_mean, result.cvi_spatial_sd)
```

`spacings` is optional if only CVI is needed; including it also produces `result.pachy_area_mm2`. The example spacing values above are illustrative **synthetic** values, not calibrated scan metadata. In real analyses, use the scanner-exported spacing. Use `parameters={...}` to override the documented vessel segmentation defaults where justified.

The dataset-level interface `extract_cvi_from_folders(cube_dir, seg_dir)` simply loads the original JPEG/SVG export convention and calls `extract_cvi`. The CLI supports single samples and cohorts by calling the same API.

### Run the synthetic example

The repository includes **five artificial JPEG B-scans, five matching SVG masks, and `info.json` metadata** under `examples/synthetic_case/`. These are programmatically generated illustrations, not clinical images.

```bash
python examples/run_single_volume.py
```

The script prints both CVI phenotypes and the area estimate. You may regenerate the example deterministically using `python examples/make_synthetic_example.py`. Numerical reference tests for the phenotype reduction logic are in `tests/`; compare the *full* synthetic example output after validating the rendering environment and final method parameters.

## What it computes

1. Select five central OCT B-scans and their corresponding SVG choroid segmentations.
2. Detect Bruch's membrane from the choroidal mask; flatten, crop, denoise, depth-equalize, and enhance contrast.
3. Detect large and medium/small vascular regions with adaptive Gaussian thresholds and morphology.
4. For each B-scan and image column `x`, compute `CVI(x) = vessel_pixels(x) / mask_pixels(x)`.
5. Calculate the **mean** and **population standard deviation** of the column-wise CVI values within each B-scan, then take the **median across the five B-scans** of each statistic.
6. Calculate the median area of the largest connected vascular region in each B-scan, using `info.json` spacing entries 2 and 3.

These definitions reflect the supplied implementation. The standard deviation reported as `cvi_spatial_sd` is not the standard deviation of all voxels or all pooled columns across the cube. Missing-mask columns contribute zero to the per-B-scan summaries in the current code; this is an important methodological choice to validate.

## Install

Use Python 3.10–3.12 in a clean environment:

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -e .
python -m unittest discover -s tests -v
```

Dependencies are specified in `pyproject.toml`. On Linux, the ReportLab SVG rasterization backend also requires system Cairo development libraries (`libcairo2-dev` and `pkg-config` on Ubuntu); the GitHub Actions workflow installs them. Record the exact resolved environment used for published results (for example `python -m pip freeze > requirements-lock.txt`) **after validating it**. Current dependency ranges are an installation specification, not an archival lockfile.

## Expected input

```text
data/
  example_001/
    oct_volume/
      info.json
      0000.jpg
      0001.jpg
      ...
    svg_segmentation/
      0000.svg
      0001.svg
      ...
```

- OCT scans must be grayscale JPEGs with consecutive, zero-based four-digit filenames.
- SVGs must correspond by B-scan stem and contain the choroid region tagged `layer layer-6` (`CC_CS`).
- `info.json` must contain `"spacing": [...]` with positive entries at indices 2 and 3 representing **millimetres per pixel** for the relevant B-scan axes. Confirm that convention against your source export.
- The current image-preprocessing crop uses fixed vertical rows 280–380 after BM alignment; input geometry must match that choice.
- Five central image/segmentation pairs are required. Input files are not distributed or uploaded.

## Run

For a single pseudonymized study:

```bash
choroid-cvi --cube-dir data/example_001/oct_volume \
  --seg-dir data/example_001/svg_segmentation \
  --sample-id example_001 --output results/phenotypes.csv
```

For a cohort, create a CSV with `sample_id,cube_dir,seg_dir` (see [example manifest](examples/manifest.example.csv), which points to the included synthetic case). Relative paths in a manifest are relative to its directory:

```bash
choroid-cvi --manifest examples/manifest.example.csv \
  --parameters examples/segmentation_parameters.example.json \
  --output results/phenotypes.csv
```

The results include `sample_id`, `n_bscans`, `cvi_mean`, `cvi_spatial_sd`, `pachy_area_mm2`, `status`, and `error`. A companion `results/phenotypes.run.json` captures versions and segmentation parameters. If a study fails, it receives a descriptive error row; the process exits with status 1 after writing all rows.

### Original code mapping

| Original | Here | Change |
| --- | --- | --- |
| `main.py` | `src/choroid_cvi/api.py` and `cli.py` | Added a reusable single-volume API and made dataset batching a convenience wrapper |
| `helper.py` | `src/choroid_cvi/processing.py` | Retained extraction calculations; removed unused imports and an unused asynchronous JSON reader |
| `svg_util.py` | `svglib==1.4.1` dependency | Avoids distributing an incomplete third-party module with unresolved relative imports; retains its renderer API |

## Reproducibility and data availability

The original imaging and participant-level genetic/clinical data are subject to UK Biobank and local cohort access conditions. Do not push raw images, segmentation files, participant identifiers, phenotypes, or credentials to public GitHub. Release eligible aggregate results through the destinations described in the manuscript's Data Availability statement.

See [methods and assumptions](docs/methods_and_assumptions.md) and the [release checklist](docs/reproducibility_checklist.md). **License, final author/citation metadata, and a permanent release DOI remain to be determined with the research team.** No licence is granted by this working copy.


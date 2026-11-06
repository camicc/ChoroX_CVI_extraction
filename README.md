# ChoroX — Choroidal vascularity from OCT

Automated extraction of choroidal vascularity index (CVI) phenotypes from macular optical coherence tomography (OCT) volumes. This repository contains the image-processing component of the workflow described in:

**Pittet C, et al.** *Genetic Architecture of Choroidal Vascularity Revealed by Population-scale Retinal Imaging*. Research Square, 2026. [https://doi.org/10.21203/rs.3.rs-8842359/v1](https://doi.org/10.21203/rs.3.rs-8842359/v1)

![ChoroX CVI extraction workflow](src/workflow.png)

## Overview

The pipeline uses OCT B-scans and their retinal/choroidal layer segmentations to derive two complementary phenotypes:

- **CVI mean:** average choroidal vascular proportion across the B-scan width.
- **CVI spatial variability:** variation in vascular proportion across the B-scan width.

The extraction procedure comprises:

1. Select the five central B-scans from a macular OCT volume.
2. Identify the choroid using the segmentation of the choriocapillaris and choroidal stroma (`CC_CS`) and align B-scans to Bruch's membrane.
3. Apply denoising, depth-dependent intensity equalization, and contrast enhancement.
4. Segment larger and smaller choroidal vessels using adaptive thresholds and morphological filtering.
5. Calculate column-wise CVI and aggregate the resulting statistics across B-scans.

The code in this repository focuses on **image processing and phenotype extraction**. The downstream genome-wide association and other genetic analyses reported in the paper are not part of this extraction script.

## Requirements

- Python 3.11 (tested on Windows with Conda)
- NumPy, SciPy, Pillow, Matplotlib, lxml, and svglib
- OpenCV
- ReportLab and a compatible bitmap-rendering backend for rasterizing SVG segmentations

One compatible installation, based on the Windows `cvi-repro` environment, is:

```powershell
conda create -n cvi-repro python=3.11 -y
conda activate cvi-repro
python -m pip install numpy scipy Pillow matplotlib lxml "svglib>=1.5,<2" "opencv-python-headless==5.0.0.93" "reportlab==4.5.1" "rl-renderPM==4.0.3"
```

**Rendering note:** ReportLab 5 no longer supports the legacy `_renderPM` backend. With the versions above, select that backend in `src/processing.py` before rasterizing an SVG:

```python
from reportlab import rl_config
rl_config.renderPMBackend = "_renderPM"
```

For reproducibility, record the complete package versions used for the final analysis. In particular, changes to the SVG renderer or OpenCV may affect pixel-level segmentations.

## Input data

The command-line script processes one OCT volume at a time. The directory passed to `--root` must contain:

```text
volume/
├── oct/
│   ├── info.json
│   ├── 0000.jpg
│   ├── ...
│   └── 0127.jpg
└── segmentation/
    ├── 0000.svg
    ├── ...
    └── 0127.svg
```

`info.json` provides image spacing under the `spacing` key. The SVGs must contain the segmentation layers expected by the extraction code, including `CC_CS`. The example follows the 128-B-scan OCT layout used for the original workflow; other acquisition geometries require validation of slice selection, cropping, and spacing.

Raw OCT images and their segmentations are **not distributed** with this repository. Use data for which you have the necessary access permissions.

## Run the extraction

From the repository root, with `cvi-repro` activated:

```powershell
python src/run_single_volume.py --root "data/your_volume_directory"
```

To calculate the CVI phenotypes from the five central B-scans without processing the entire volume for the en face projection:

```powershell
python src/run_single_volume.py --root "data/your_volume_directory" --skip-full-volume
```

On a 128-B-scan volume, the default central indices are **62–66** (files `0062.jpg` through `0066.jpg` and matching SVGs). Paths containing spaces must be quoted.

### Outputs

The script prints the volume-level **CVI mean** and **CVI spatial variability**, as well as per-B-scan diagnostics. It also saves:

| Output | Description |
| --- | --- |
| `segmentation_debug/` | Vessel masks and OCT overlays for the five central B-scans |
| `enface_projection.png` | Full-volume en face vessel projection (omitted with `--skip-full-volume`) |

Outputs are written under the input volume directory. The script may also report a legacy connected-component area for compatibility; **it should not be interpreted as a validated vessel-area biomarker**, because the original area calculation includes the background component.

## CVI definitions

For B-scan $i$ and image column $x$, the column-wise vascularity index is

$$
\mathrm{CVI}_i(x) = \frac{N_{\mathrm{vessel},i}(x)}{N_{\mathrm{choroid},i}(x)},
$$

where the numerator and denominator are the vessel and choroidal-mask pixel counts within the column. The implementation handles columns without valid choroidal pixels explicitly.

For each B-scan, the mean and standard deviation are computed **across columns**. The final volume-level phenotypes are the medians of those five per-scan statistics:

$$
\mathrm{CVI}_{\mathrm{mean}} = \operatorname{median}_i\!\left[\operatorname{mean}_x\mathrm{CVI}_i(x)\right],
\qquad
\mathrm{CVI}_{\mathrm{spatial\ SD}} = \operatorname{median}_i\!\left[\operatorname{SD}_x\mathrm{CVI}_i(x)\right].
$$

**CVI spatial variability is not the standard deviation of CVI means between B-scans.** It measures the horizontal variation within each B-scan, followed by median aggregation across the five central scans.

## Reproducibility and scope

The pipeline retains the original image-processing operations and their order. Its settings were designed for the OCT geometry used in the study; applying them to other devices, datasets, or segmentation formats requires technical validation. Before claiming numerical equivalence with the original analysis, compare the resulting masks and CVI values on the same reference OCT volume.

The linked article is a **Research Square preprint (version 1)**. Refer to that version when citing the methods associated with this code.

## Citation

If you use the extraction pipeline, please cite:

> Pittet C, Böttger L, Eandi CM, et al. Genetic Architecture of Choroidal Vascularity Revealed by Population-scale Retinal Imaging. *Research Square*. 2026. https://doi.org/10.21203/rs.3.rs-8842359/v1

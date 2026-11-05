# CVI extraction: methods and assumptions

## Scientific definition

The workflow estimates the relative area occupied by choroidal vessels within each lateral OCT column. For each volume, it computes two complementary phenotypes from the five central B-scans:

1. **CVI mean**: median across the five B-scans of the per-column CVI means.
2. **CVI spatial variability**: median across the five B-scans of the per-column CVI population standard deviations.

The column-wise CVI denominator is the **total choroidal-mask area within the column**, not luminal area. This definition deliberately differs from taking one global lumen-to-choroid area ratio, because the lateral variation is itself a phenotype.

## Processing steps

| Stage | Implementation | Interpretation |
| --- | --- | --- |
| OCT input | An exported macular OCT volume, matched to its SVG segmentation for each B-scan | Raw acquisition and AI choroid-layer segmentation are upstream |
| Central region | Select five consecutive central B-scans | UK Biobank: 128 B-scans per volume; OphtalmoLaus: 256 |
| Choroidal mask | Extract SVG `CC_CS` / `layer layer-6` region | Depends on anatomical accuracy of the upstream segmentation |
| Alignment | Estimate Bruch's membrane as the smoothed upper choroidal-mask boundary (Gaussian σ = 4); align it to row 300 | The image geometry must support the fixed alignment/crop |
| Crop | Use aligned vertical rows 280–380 | Fixed 100-pixel window; not dynamically resized for a different scanner |
| 3-D preprocessing | Median and anisotropic Gaussian filters, depth-dependent intensity equalization, contrast scaling | Uses information across the five selected B-scans |
| Vessel segmentation | Separate adaptive Gaussian thresholds for large and smaller vessels, morphological opening/closing, exclusion of mask boundaries | Output is a binary vessel representation within the choroid |
| Column-wise CVI | Vessel pixel count divided by choroidal-mask pixel count for each lateral column | Zero-mask columns are assigned 0 by the current implementation |
| Phenotype aggregation | Per-B-scan mean and population SD of column-wise CVI; median of each statistic across five B-scans | No pooling of the full volume before computing SD |
| Largest connected vessel area | Median largest connected-region area, converted with pixel spacings | Optional output, `pachy_area_mm2` |

## Segmentation parameters

The [default parameters](../src/choroid_cvi/api.py) and an [editable JSON example](../examples/segmentation_parameters.example.json) cover:

| Role | Parameter | Default |
| --- | --- | ---: |
| Large-vessel threshold window | `block_size_pachy` | 91 |
| Large-vessel threshold subtraction | `cst_sub_mean_pachy` | 13 |
| Large-vessel morphological closing | `kernel_size_c_pachy` | 5 |
| Large-vessel morphological opening | `kernel_size_o_pachy` | 17 |
| Smaller-vessel threshold window | `block_size` | 41 |
| Smaller-vessel threshold subtraction | `cst_sub_mean` | 16 |
| Smaller-vessel morphological opening | `kernel_size_o` / `kernel_size_o_small` | 5 / 5 |
| Exclusion near choroidal boundary | distance transform `> 8` pixels | fixed |

In the manuscript, OphtalmoLaus processing uses a device-calibrated vessel-segmentation threshold to account for acquisition and luminosity differences, while retaining the same phenotype definitions. Record the chosen parameter file and image-processing environment for **each** cohort; do not assume an uncalibrated parameter set transfers between OCT platforms.

## Image quality and cohort selection

The extraction API analyzes one supplied volume; it does **not** implement cohort selection, image-quality grading, or participant/genotype QC. The manuscript describes upstream UK Biobank image QC using PLUME-OCT, with SNR ≥ 10 and spatially weighted displacement ≤ 1 among the scan-quality conditions, and exclusion of missing/empty choroid segmentations. Exact sample selection and genetic analysis exclusions belong to the study-level analysis workflow.

The method requires choroid masks from a separate AI-based retinal layer segmentation system; those model weights and inference steps are not supplied with these OCT/SVG processing functions. Both the mask generation method and image-quality assessment should be documented separately if the entire raw-OCT-to-GWAS process is to be reproduced.

## Units and edge cases

- CVI mean and CVI spatial variability are dimensionless proportions. `cvi_spatial_sd` is not a standard deviation across participants.
- `spacings[2]` and `spacings[3]` are assumed to be in mm/pixel and correspond to the two within-B-scan axes for area calculation.
- The folder adapter expects `info.json` with those spacings. The array-based `extract_cvi` call does not require spacings when only CVI is needed.
- If an individual B-scan contains no choroid-mask pixels, its contribution is skipped by the numerical reduction. If all masks are empty, phenotype computation fails. This condition should be filtered by image QC.
- If a column contains no choroid-mask pixels, its CVI value is set to 0 and included in within-B-scan mean and SD. Because this can affect both outcomes, segmentation coverage must be inspected, especially across OCT devices.
- Fixed crop coordinates can truncate choroidal structures outside the expected depth window. Verify representative outputs before large-scale extraction.
- ReportLab/SVG rasterization may depend on graphics-library versions. Use `environment.yml`, compare masks visually, and validate reference phenotypes before publication.

## What a downstream user should do

Check one volume from the intended OCT device, inspect its SVG masks and vessel segmentation, verify pixel spacing and cropping, and confirm that the chosen threshold parameters produce anatomically plausible choroidal vessel masks. Then process the cohort and retain per-run parameter and version metadata.

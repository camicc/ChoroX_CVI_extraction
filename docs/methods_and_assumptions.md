# Methods mapping and implementation assumptions

The code was reorganized without intentionally changing the scientific calculations in the supplied `helper.py`. Its numerical outputs have **not** been benchmarked against the original production outputs.

| Stage | Implementation | Item requiring paper/production verification |
| --- | --- | --- |
| B-scan selection | Five central consecutive B-scans by image count | Final image count, orientation, and numbering for UKB/OphtalmoLaus |
| Segmentation | Existing SVG `CC_CS` masks define the choroid | How SVG masks were produced (AI model/version), validation, and availability |
| BM alignment | Upper boundary of choroid mask, Gaussian smoothing (sigma 4); BM target row 300 | Correct scanner geometry and BM definition |
| Cropping | Rows 280–380 after flattening | Is the same crop correct across imaging platforms? |
| Denoising | 3-D median filtering, intensity depth equalization, Gaussian filtering | Whether cross-B-scan filtering matches the archived pipeline |
| Vessel segmentation | Two Gaussian adaptive thresholds (block sizes 91 and 41), ellipse morphology, mask edge exclusion (distance >8 pixels) | Parameters tuned by cohort, final defaults, and visual QC |
| CVI per column | Count foreground pixels / count choroid-mask pixels; missing-mask columns set to 0 | Whether invalid columns should be excluded, and whether masks retain exact scale |
| Phenotypes | Median over B-scans of per-column mean and SD | Confirm definition of *spatial variability* in Methods and supplied UKB GWAS phenotypes |
| Pachy area | Median over B-scans of largest connected-component pixel area × `spacing[2]` × `spacing[3]` | Spacing axes and units in the exported JSON |

## Known caveats in the supplied implementation

- Incomplete choroidal masks may cause errors; current batch wrapper records a failed subject instead of hiding it.
- Rescaling a slice and whole cube to 8-bit can change small intensity differences; platform-specific acquisition settings matter.
- The code processes JPEG B-scans and SVG masks, not native OCT/DICOM. It cannot reproduce the upstream segmentation/model inference stage without that separate code and weights.
- The image-level `cvi_spatial_sd` reported here is the **median of within-B-scan standard deviations**. This differs mathematically from the SD of all columns pooled across five B-scans.
- For zero-mask columns, `CVI(x)` is assigned zero and included in the mean and SD. This can lower the measured CVI and alter its spatial variability. Reassess only as a separately versioned method change after comparing the original and revised phenotypes.
- The `pachy_only`/`medium_only` modes produce images on potentially different numeric scales; they are not exposed by the combined-phenotype CLI.
- Fixed crop indices can truncate a choroid that lies outside the expected depth window.
- `svg_util.py` was a copy of a third-party `svglib` module containing unresolved `from .utils` and `from .fonts` imports, so it is not redistributed. Check the SVG rendering against the environment used in production.

Any scientifically motivated fix to the above should be a distinct version accompanied by a comparison table and correlation/shift plots versus the published phenotypes. Do not silently update the code underlying the reported GWAS.


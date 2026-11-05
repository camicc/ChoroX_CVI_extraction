# Reproducibility and public-release checklist

## Extraction and environment

- [x] Provide a single-volume Python API and a minimal OCT/SVG example.
- [x] Include five synthetic central B-scans, corresponding SVG masks, and example input metadata.
- [x] Document input folder format, phenotype definitions, device-specific parameters, and batch execution.
- [x] Supply a Conda/Mamba environment file and a cross-platform CI test workflow.
- [ ] Confirm the five-B-scan central selection and per-column CVI calculation against archived output for authorized OCT volumes.
- [ ] Record validated UK Biobank and OphtalmoLaus segmentation parameter files (and their scan geometry assumptions).
- [ ] Verify SVG rasterization, selected masks, and resulting vessel maps visually on authorized samples.
- [ ] Run all tests and the synthetic example successfully in a fresh Conda/Mamba environment on Windows and Linux; record expected numeric output.
- [ ] Archive a resolved dependency lock or explicit working environment export.
- [ ] Document upstream choroid segmentation model/weights and PLUME-OCT image-quality processing, or supply controlled-access references.

## Genetic analysis

- [ ] Add the exact commands/configurations for GWAS input QC, phenotype transformations, covariate construction, and BGENIE association tests.
- [ ] Add LD-clumping and genetic signal reporting configuration, reference genome/build and reference data versions.
- [ ] Add LDSC, PascalX, functional/regulatory annotation, GWAS Catalog overlap, MR, and genetic replication scripts or link to a versioned companion repository.
- [ ] Map manuscript GWAS figures/tables to reproducible commands and shareable, nonrestricted inputs/summary statistics.
- [ ] Review the [GWAS inventory](gwas_reproducibility.md) with the genetic analysis team.

## Publication

- [ ] Obtain permission to release the code and ensure no participant data, restricted OCT material, patient IDs, tokens, or institutional paths are tracked.
- [ ] Choose an appropriate software license with the research team.
- [ ] Complete author/citation metadata and manuscript citation.
- [ ] Verify the final manuscript Methods and repository parameter documentation agree.
- [ ] Publish a tagged code release and archive its DOI when appropriate.

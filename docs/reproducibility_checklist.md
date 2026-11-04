# Pre-release checklist for a Nature Communications code repository

- [ ] Confirm the exact commit of the UK Biobank extraction pipeline used to generate the final GWAS phenotypes.
- [ ] Confirm all segmentation and preprocessing parameters, including any separate OphtalmoLaus settings.
- [ ] Confirm the SVG mask provenance, upstream choroid model, software/weights availability, and relevant restrictions.
- [ ] Verify per-column CVI formula, handling of missing-mask columns, B-scan aggregation, units, and correspondence with manuscript equations.
- [ ] Compare output against archived production phenotypes for at least a small, authorized test subset; document numeric tolerance.
- [x] Add a distributable synthetic five-B-scan OCT+SVG example and a single-volume usage script.
- [ ] Validate the complete synthetic example and record trusted reference numbers after SVG rasterization is verified in the final environment.
- [ ] Freeze the exact dependency environment; run tests in a fresh environment and record a successful end-to-end example.
- [ ] Add the actual GWAS, QC, LD clumping, heritability, gene-based analysis, functional annotation, eQTL/coloc, MR and external replication workflows—or point to versioned, accessible companion repositories.
- [ ] Publish sample inclusion/QC logic, input data schemas and source-access instructions without disclosing sensitive data.
- [ ] Include a script/table mapping paper figures and tables to the precise analysis commands and aggregated output files.
- [ ] Confirm approval for code release, choose the repository license, and complete author/citation metadata.
- [ ] Check repository for personal data, institution-specific paths, tokens, licensed assets, large files, and generated outputs.
- [ ] Tag the accepted manuscript code version; archive that release in Zenodo (or another persistent repository) and cite its DOI when available.


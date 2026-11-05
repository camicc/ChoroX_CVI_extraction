# GWAS and downstream analyses: reproducibility inventory

This repository contains the OCT-to-phenotype extraction module. Reproducing the genetic results additionally requires documented analysis of the resulting CVI mean and CVI spatial variability phenotypes. The genetic analyses are not run by the `choroid-cvi` command.

## What should be released?

The standard scientific approach is to **publish the study-specific scripts/configuration and exact calls to established software**, without redistributing BGENIE, PLINK2, LDSC, PascalX, TwoSampleMR, or access-controlled UK Biobank individual-level data. For each analysis step, a reproducible package should specify input field names, inclusion criteria, parameter values, software versions, reference-build coordinates, and expected output columns.

| Analysis stage in manuscript | Details stated in the draft | Minimum accompanying code or record |
| --- | --- | --- |
| Participant and image QC | Right-eye macular OCT, most recent UK Biobank scan where repeated; PLUME-OCT QC; nonempty choroid mask | Filtering logic, sample flow, OCT-to-participant linkage schema, excluded reason codes |
| Phenotype processing | CVI mean and CVI spatial variability; normal-score transformation and z-score sensitivity analysis | Exact transformation, treatment of missing values/outliers, order of covariate adjustment, input/output schema |
| GWAS covariates | Sex, age, age², refractive power, first 20 genetic principal components | Covariate construction code, encoding, units, exclusion logic |
| Variant QC and association | MAF ≥ 0.01, imputation quality ≥ 0.8; BGENIE two-step linear modelling; `P < 5 × 10⁻⁸` | BGENIE version, executable commands, filtering commands, genotype build, summary-statistics schema |
| Lead-variant selection | PLINK2 clumping; UK10K hg19 reference; `r² < 0.1` and 500 kb separation in the draft | Exact PLINK2 invocations, reference panel/variant IDs, ancestry, algorithm and selection settings; lead-SNP table |
| Heritability and cross-trait correlation | Univariate and bivariate LDSC | Versions, munging commands, LD score source, intercept options and output files |
| Gene/pathway/cell-type enrichment | PascalX; gene window ±50 kb; UK10K LD reference; reference pathways/MSigDB and cell-type sets | Exact gene coordinate/build mapping, PascalX configuration, set versions, multiple-testing implementation |
| eQTL, sQTL, and colocalisation | EyeGEx retina and GTEx vascular/blood sources | QTL release/version, harmonization, allele and reference matching, priors, statistical models, significance decisions |
| GWAS Catalog overlap | Comparison of lead variants/LD loci to prior traits | GWAS Catalog snapshot date, LD criteria, filtering and annotation scripts |
| Bidirectional Mendelian randomization | TwoSampleMR v0.5.7; AMD FinnGen/OpenGWAS `finn-b-H7_AMD`; European 1000 Genomes LD pruning at `r² < 0.001` | Instrument selection, harmonization, palindromic SNP handling, IVW and sensitivity calls, multiple testing |
| OphtalmoLaus genetic replication | Device-calibrated vessel threshold, covariate correction, direct/proxy SNP matching, UKB White British reference, `r² ≥ 0.8` | Cohort-specific phenotype settings, direct/proxy matching, LD sign/orientation, harmonization, FDR script |
| Tables and figures | GWAS Manhattan/QQ plots, genetic results and replication tables | Figure-generating scripts, plot input schemas, exact numeric output provenance |

This table records the methods currently described in the manuscript draft, **not a validated run configuration**. Check the final manuscript, input files, clumping settings, and reference data before publishing any analysis scripts.

## Recommended organization if the GWAS scripts are shared

One companion repository, or a separate `analysis/` directory here, could contain:

```text
analysis/
  README.md               # execution order, dependencies, reference datasets
  config/                 # analysis settings without participant paths/credentials
  phenotype_qc/           # inclusion, transformation, covariate preparation
  association/            # BGENIE invocation and summary-stat harmonization
  clumping_ldsc/          # lead signals and genetic correlations
  enrichment_qtl/         # PascalX and expression-based analyses
  mendelian_randomization/
  replication/            # independent cohort and LD proxy handling
  figures_tables/         # manuscript output assembly
```

This is a **suggested structure**, not a claim that these directories or executable scripts already exist. A single workflow script is equally acceptable if its configuration and outputs are documented.

### Minimal release versus full reproduction

- **Required for the imaging extraction repository:** the executable CVI code, environment, input schema, parameter files, and a synthetic example. Those are provided here.
- **Strongly recommended for the genetic results reported in the paper:** exact analysis commands and custom scripts for the main GWAS, clumping, MR, replication and figures, plus software versions. If the code resides in another repository, link to its stable tag.
- **Not necessary to distribute:** copies of standard statistical packages, individual-level UKB/OphtalmoLaus files, or licensed LD panels. Document their source and authorized access route instead.

If the full analysis code cannot be released yet, keep the manuscript's Code Availability statement precise: distinguish the public OCT extraction module from the downstream genetic analysis material, and avoid implying that the entire paper can be regenerated from this repository alone.

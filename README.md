# Plasma Proteome and Pancreatic Ductal Adenocarcinoma Risk: A Mendelian Randomization and Observational Analysis​

Mendelian randomization (MR) and colocalization code for UKB-PPP and deCODE proteins, using PanScan1234 + FinnGen R13 pancreatic cancer GWAS.

## Analysis files

- [`analysis/run_mr.R`](analysis/run_mr.R): IVW, weighted median, MR-Egger, PIVW and MR-RAPS.
- [`analysis/run_colocalization.R`](analysis/run_colocalization.R): regional Bayesian colocalization.

Both scripts are self-contained and support UKB and deCODE.

## Inputs and scope

These scripts use the final prepared inputs; raw-data cleaning, coordinate conversion, harmonization and LD clumping are separate upstream steps.

Place these four files in one input directory:

- `MR_instruments.tsv.gz`: harmonized, LD-clumped IV-level exposure/outcome statistics.
- `Assay_manifest.tsv`: assay identifiers, target labels and bootstrap seeds.
- `Coloc_inputs.tsv.gz`: full quality-filtered regional statistics, not just IVs.
- `Coloc_regions.tsv`: lead variants, region boundaries and SNP counts.

Required columns are checked in each script. `cohort` is `ukb` or `decode`; `protein` identifies an assay. Variants use normalized GRCh37 `chr:position:REF:ALT` keys, with exposure/outcome effects already aligned. Outcome beta is log(OR) for ALT.

Data and reference panels are not included. 

## Run

Requires R and `data.table`, `R.utils`, `MendelianRandomization`, `mr.pivw`, `mr.raps` and `coloc`.

```bash
Rscript analysis/run_mr.R /path/to/prepared_inputs mr_results all
Rscript analysis/run_colocalization.R /path/to/prepared_inputs coloc_results all
```

Replace `all` with `ukb` or `decode` to run one cohort. Output directories must be new or empty. 

## Analysis settings

MR preserves per-assay seeds; PIVW uses 10,000 bootstraps. Unavailable estimates remain `NA`. Bonferroni and BH-FDR corrections use all positive-IV assays within each cohort and method: 1,985 UKB and 1,789 deCODE in the complete primary inputs. Complete cohort inputs are required to reproduce study-adjusted P-values.

Colocalization uses ±150 kb IV-lead regions for IVW-FDR-selected assays, quantitative exposures (`sdY = 1`), a binary outcome, and priors p1 = p2 = 10⁻⁴ and p12 = 10⁻⁵. PPH4 ≥ 0.8 indicates support. The beta/variance model does not require outcome sample-size assumptions.

Study versions: R 4.4.3; data.table 1.17.0; MendelianRandomization 0.10.0; mr.pivw 0.1.3; mr.raps 0.4.2; coloc 5.2.3.

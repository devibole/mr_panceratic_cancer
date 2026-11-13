# Plasma Proteome and Pancreatic Ductal Adenocarcinoma Risk: A Mendelian Randomization and Observational Analysis​
MR and colocalization code for UKB-PPP and deCODE proteins using a meta-analyzed PanScan1234 + FinnGen R13 pancreatic cancer GWAS for the outcome.

## Scripts

- [`analysis/run_ukb_mr.R`](analysis/run_ukb_mr.R): raw UKB `.tar` → numeric cis filtering → quality control → allele-aware GWAS matching → LD clumping → MR.
- [`analysis/run_decode_mr.R`](analysis/run_decode_mr.R): the same analysis from raw deCODE `.txt.gz`, using its native column names and allele annotations.
- [`analysis/run_colocalization.R`](analysis/run_colocalization.R): regional colocalization for either cohort, from full regional summary statistics.

## Inputs and running

```bash
Rscript analysis/run_ukb_mr.R /path/to/protein.tar ukb_output
Rscript analysis/run_decode_mr.R /path/to/protein.txt.gz decode_output
Rscript analysis/run_colocalization.R /path/to/regional_inputs coloc_output all
```

Colocalization requires `Coloc_inputs.tsv.gz` and `Coloc_regions.tsv`; these contain full regional variants, not just MR IVs. Replace `all` with `ukb` or `decode` for one cohort. 

## Settings

Cis windows are gene body ±1 Mb, assessed in UKB GRCh38 or deCODE GRCh37 coordinates. Analysis variants use normalized GRCh37 `chr:position:REF:ALT`. Exposure P ≤ 5×10⁻⁸ and the original quality filters are retained; LD clumping uses p1 = 0.5, r² = 0.1 and 2,000 kb. Single-assay runs report unadjusted P-values; cohort-wide Bonferroni/FDR correction requires the complete positive-IV assay family.

Colocalization uses IV-lead ±150 kb regions for IVW-FDR-selected assays, `coloc.abf`, exposure `sdY = 1`, a binary outcome and priors p1 = p2 = 10⁻⁴, p12 = 10⁻⁵.

Requires R, PLINK 1.9, gzip, awk, and the R packages named in each script. Study versions: R 4.4.3, data.table 1.17.0, MendelianRandomization 0.10.0, mr.pivw 0.1.3, mr.raps 0.4.2 and coloc 5.2.3. 

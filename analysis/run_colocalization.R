#!/usr/bin/env Rscript
# Basic regional colocalization for UKB and/or deCODE using the study's
# cleaned, allele-aligned GRCh37 regional inputs. This does not rebuild sources.
# Usage: Rscript analysis/run_colocalization.R INPUTS_DIR NEW_OUTPUT_DIR [all|ukb|decode]

required_packages <- c("data.table", "R.utils", "coloc")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Install required R packages: ", paste(missing_packages, collapse = ", "))
}
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: Rscript analysis/run_colocalization.R INPUTS_DIR NEW_OUTPUT_DIR [all|ukb|decode]")
}
input_dir <- normalizePath(args[1L], mustWork = TRUE)
output_dir <- args[2L]
cohort_choice <- if (length(args) == 3L) args[3L] else "all"
if (!cohort_choice %in% c("all", "ukb", "decode")) {
  stop("Cohort must be all, ukb, or decode")
}
check_output_dir <- function() {
  if ((file.exists(output_dir) && !dir.exists(output_dir)) ||
      (dir.exists(output_dir) && length(list.files(output_dir, all.files = TRUE,
                                                 no.. = TRUE)))) {
    stop("Use a new or empty output directory; existing files will not be overwritten")
  }
}
check_output_dir()

input_file <- file.path(input_dir, "Coloc_inputs.tsv.gz")
region_file <- file.path(input_dir, "Coloc_regions.tsv")
if (!all(file.exists(c(input_file, region_file)))) {
  stop("INPUTS_DIR must contain Coloc_inputs.tsv.gz and Coloc_regions.tsv")
}
inputs <- fread(input_file)
regions <- fread(region_file)
input_columns <- c("cohort", "protein", "normalized_key", "chr", "pos_bp37",
                   "beta_exposure", "se_exposure", "beta_outcome", "se_outcome")
region_columns <- c("cohort", "protein", "region_id", "lead_variant", "chr",
                    "lead_pos_bp37", "start_bp37", "end_bp37", "n_snps")
if (!all(input_columns %in% names(inputs)) ||
    !all(region_columns %in% names(regions))) {
  stop("Missing required input columns; see README.md")
}
# Ignore any supplied result/posterior columns: every posterior is fitted anew.
inputs <- inputs[, ..input_columns]
regions <- regions[, ..region_columns]
if (anyNA(inputs$cohort) || anyNA(regions$cohort) ||
    any(!inputs$cohort %in% c("ukb", "decode")) ||
    any(!regions$cohort %in% c("ukb", "decode"))) {
  stop("Input cohort labels must be ukb or decode")
}
if (cohort_choice != "all") {
  inputs <- inputs[cohort == cohort_choice]
  regions <- regions[cohort == cohort_choice]
}
if (!nrow(inputs) || !nrow(regions)) stop("No regional inputs for the requested cohort")
nonempty <- function(x) !anyNA(x) && all(nzchar(as.character(x)))
finite_numeric <- function(x) is.numeric(x) && all(is.finite(x))
if (!all(vapply(inputs[, .(protein, normalized_key)], nonempty, logical(1))) ||
    !all(vapply(regions[, .(protein, region_id, lead_variant)], nonempty, logical(1))) ||
    anyDuplicated(inputs[, .(cohort, protein, normalized_key)]) ||
    anyDuplicated(regions$region_id)) {
  stop("Missing identifiers or duplicate variant/region identifiers")
}
input_numeric <- setdiff(input_columns, c("cohort", "protein", "normalized_key"))
region_numeric <- setdiff(region_columns, c("cohort", "protein", "region_id", "lead_variant"))
if (!all(vapply(inputs[, ..input_numeric], finite_numeric, logical(1))) ||
    !all(vapply(regions[, ..region_numeric], finite_numeric, logical(1))) ||
    any(inputs$se_exposure <= 0 | inputs$se_outcome <= 0) ||
    any(!inputs$chr %in% 1:22) || any(!regions$chr %in% 1:22) ||
    any(inputs$pos_bp37 < 1 | inputs$pos_bp37 != floor(inputs$pos_bp37)) ||
    any(regions$lead_pos_bp37 < 1 |
          regions$lead_pos_bp37 != floor(regions$lead_pos_bp37)) ||
    any(regions$n_snps < 1 | regions$n_snps != floor(regions$n_snps)) ||
    any(regions$start_bp37 != pmax(1, regions$lead_pos_bp37 - 150000)) ||
    any(regions$end_bp37 != regions$lead_pos_bp37 + 150000)) {
  stop("Invalid effects, SEs, coordinates, SNP counts, or lead +/-150 kb region bounds")
}
key_prefix <- paste0(inputs$chr, ":", format(inputs$pos_bp37,
                                           scientific = FALSE, trim = TRUE), ":")
if (any(!grepl("^[^:]+:[0-9]+:[^:]+:[^:]+$", inputs$normalized_key)) ||
    any(!startsWith(inputs$normalized_key, key_prefix))) {
  stop("Variant keys must be GRCh37 chr:position:REF:ALT matching the coordinate columns")
}

# Regions are the production FDR-selected lead windows. Use ALL available
# quality-filtered regional variants, not only significant/clumped MR IVs.
results <- vector("list", nrow(regions))
for (i in seq_len(nrow(regions))) {
  region <- regions[i]
  x <- inputs[cohort == region$cohort & protein == region$protein &
                chr == region$chr & pos_bp37 >= region$start_bp37 &
                pos_bp37 <= region$end_bp37]
  lead <- x[normalized_key == region$lead_variant]
  if (nrow(x) != region$n_snps || nrow(lead) != 1L ||
      lead$pos_bp37 != region$lead_pos_bp37 ||
      anyDuplicated(x$normalized_key)) {
    stop("Incomplete or inconsistent regional data: ", region$region_id)
  }
  # With beta/varbeta supplied, the ABF model does not need N or case fraction.
  # FinnGen-only counts must not be presented as whole-meta-GWAS sample sizes.
  fit <- coloc::coloc.abf(
    list(beta = x$beta_exposure, varbeta = x$se_exposure^2,
         snp = x$normalized_key, type = "quant", sdY = 1),
    list(beta = x$beta_outcome, varbeta = x$se_outcome^2,
         snp = x$normalized_key, type = "cc"),
    p1 = 1e-4, p2 = 1e-4, p12 = 1e-5
  )
  pp <- unname(fit$summary[paste0("PP.H", 0:4, ".abf")])
  if (length(pp) != 5L || any(!is.finite(pp)) || abs(sum(pp) - 1) > 1e-9) {
    stop("Invalid posterior probabilities: ", region$region_id)
  }
  results[[i]] <- cbind(region, as.data.table(as.list(setNames(pp, paste0("PPH", 0:4)))))
}

check_output_dir()
if (!dir.exists(output_dir) && !dir.create(output_dir, recursive = TRUE)) {
  stop("Could not create output directory")
}
fwrite(rbindlist(results), file.path(output_dir, "Coloc_regions.tsv"),
       sep = "\t", na = "NA")
capture.output(sessionInfo(), file = file.path(output_dir, "R_session.txt"))
message("Completed ", nrow(regions), " regions for ", cohort_choice, ".")

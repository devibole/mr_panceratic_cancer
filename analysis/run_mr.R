#!/usr/bin/env Rscript
# Primary MR for UKB and/or deCODE from the study's exact, already harmonized
# and LD-clumped IV inputs. This script does not rebuild or reselect IVs.
# Run computational analyses on an allocated compute node, not a login node.

args <- commandArgs(trailingOnly=TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop('Usage: Rscript analysis/run_mr.R INPUTS_DIR NEW_OUTPUT_DIR [all|ukb|decode]')
}
cohort_choice <- if (length(args) == 3L) args[3] else 'all'
if (!cohort_choice %in% c('all', 'ukb', 'decode')) {
  stop('Cohort must be all, ukb or decode')
}
inputs <- normalizePath(args[1], mustWork=TRUE)
out <- args[2]
if (file.exists(out) && !dir.exists(out)) stop('Output path is an existing file')
if (dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE))) {
  stop('Use a new or empty output directory; existing results will not be overwritten')
}
required_packages <- c('data.table', 'R.utils', 'MendelianRandomization', 'mr.pivw', 'mr.raps')
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                             logical(1), quietly=TRUE)]
if (length(missing_packages)) stop('Missing R packages: ', paste(missing_packages, collapse=', '))
suppressPackageStartupMessages(library(data.table))
setDTthreads(1L)
Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')

# All five original calls, with the original package defaults. Errors are kept
# as unavailable results, not converted to P=1 or zero. The Egger normal-test
# P value uses the upper tail to prevent 1-CDF cancellation to a false zero.
fit_original_method <- function(d,method,seed) {
  if(!nrow(d))return(data.table::data.table(method,estimate=NA_real_,se=NA_real_,
    lower=NA_real_,upper=NA_real_,p=NA_real_,egger_intercept_p=NA_real_,
    status='No eligible instruments',error=NA_character_))
  set.seed(seed)
  z<-MendelianRandomization::mr_input(bx=d$beta_exposure,bxse=d$se_exposure,
    by=d$beta_outcome,byse=d$se_outcome,snps=d$SNPID)
  error<-NA_character_
  fit<-tryCatch(suppressWarnings(switch(method,
    IVW=MendelianRandomization::mr_ivw(z),
    Median=MendelianRandomization::mr_median(z),
    Egger=MendelianRandomization::mr_egger(z),
    PIVW=mr.pivw::mr_pivw(Bx=d$beta_exposure,Bxse=d$se_exposure,
      By=d$beta_outcome,Byse=d$se_outcome,n.boot=10000),
    RAPS=mr.raps::mr.raps(data=data.frame(beta.exposure=d$beta_exposure,
      beta.outcome=d$beta_outcome,se.exposure=d$se_exposure,se.outcome=d$se_outcome),diagnostics=FALSE))),
    error=function(e){error<<-conditionMessage(e);NULL})
  getslot<-function(n)if(isS4(fit)&&n %in% slotNames(fit))as.numeric(slot(fit,n)) else NA_real_
  if(method=='RAPS') {
    b<-if(is.null(fit))NA_real_ else fit$beta.hat;s<-if(is.null(fit))NA_real_ else fit$beta.se
    lo<-b-1.96*s;hi<-b+1.96*s;p<-2*pnorm(abs(b/s),lower.tail=FALSE);ip<-NA_real_
  } else {
    eg<-method=='Egger';b<-getslot('Estimate');s<-getslot(if(eg)'StdError.Est' else 'StdError')
    lo<-getslot(if(eg)'CILower.Est' else 'CILower');hi<-getslot(if(eg)'CIUpper.Est' else 'CIUpper')
    p<-getslot(if(eg)'Pvalue.Est' else 'Pvalue');ip<-if(eg)getslot('Pvalue.Int') else NA_real_
    if(eg&&is.finite(b)&&is.finite(s)&&s>0)p<-2*pnorm(abs(b/s),lower.tail=FALSE)
  }
  data.table::data.table(method,estimate=b,se=s,lower=lo,upper=hi,p,
    egger_intercept_p=ip,status=if(is.finite(b)&&is.finite(s)&&s>0)'Estimated' else 'Method unavailable',error)
}

input_files <- file.path(inputs, c('MR_instruments.tsv.gz', 'Assay_manifest.tsv'))
if (!all(file.exists(input_files))) stop('INPUTS_DIR must contain MR_instruments.tsv.gz and Assay_manifest.tsv')
d <- fread(input_files[1], colClasses=c(normalized_key='character'))
m <- fread(input_files[2])
need_d <- c('cohort', 'protein', 'normalized_key', 'beta_exposure', 'se_exposure',
            'beta_outcome', 'se_outcome')
need_m <- c('cohort', 'protein', 'protein_symbol', 'job_id', 'bootstrap_seed')
if (!all(need_d %in% names(d))) stop('Missing MR input columns: ', paste(setdiff(need_d, names(d)), collapse=', '))
if (!all(need_m %in% names(m))) stop('Missing manifest columns: ', paste(setdiff(need_m, names(m)), collapse=', '))
nonempty <- function(x) !is.na(x) & nzchar(trimws(as.character(x)))
valid_integer <- function(x) is.numeric(x) && all(is.finite(x) & x > 0 &
                                                  x == floor(x) & x <= .Machine$integer.max)
if (!all(m$cohort %in% c('ukb', 'decode')) || !all(d$cohort %in% c('ukb', 'decode')) ||
    !all(nonempty(m$protein)) || !all(nonempty(m$protein_symbol)) ||
    !all(nonempty(d$protein)) || anyDuplicated(m[, .(cohort, protein)]) ||
    !valid_integer(m$job_id) || anyDuplicated(m$job_id) ||
    !valid_integer(m$bootstrap_seed) || anyDuplicated(m$bootstrap_seed)) {
  stop('Manifest must have unique cohort/assay IDs, unique positive job IDs and seeds, and nonempty target labels')
}
if (anyNA(d$normalized_key) ||
    !all(grepl('^([1-9]|1[0-9]|2[0-2]):[1-9][0-9]*:[ACGT]+:[ACGT]+$', d$normalized_key)) ||
    anyDuplicated(d[, .(cohort, protein, normalized_key)])) {
  stop('IVs require unique normalized GRCh37 chr:position:REF:ALT keys within each assay')
}
alleles <- tstrsplit(d$normalized_key, ':', fixed=TRUE, keep=3:4)
if (length(alleles) != 2L || any(alleles[[1]] == alleles[[2]])) stop('REF and ALT must differ')
for (field in c('beta_exposure', 'beta_outcome', 'se_exposure', 'se_outcome')) {
  if (!is.numeric(d[[field]]) || any(!is.finite(d[[field]]))) stop('Invalid values in ', field)
}
if (any(d$se_exposure <= 0 | d$se_outcome <= 0)) stop('All standard errors must be positive')
if ('distance_to_gene_bp' %in% names(d) &&
    (!is.numeric(d$distance_to_gene_bp) || any(!is.finite(d$distance_to_gene_bp) |
         d$distance_to_gene_bp < 0 | d$distance_to_gene_bp > 1000000))) {
  stop('An instrument is outside the numeric gene-body +/-1 Mb cis window')
}
if (nrow(d[!m, on=.(cohort, protein)])) stop('Every IV assay must match exactly one manifest entry')
if (cohort_choice != 'all') {
  d <- d[cohort == cohort_choice]
  m <- m[cohort == cohort_choice]
}
# Only assays with an eligible IV enter the reported primary testing family.
m <- m[unique(d[, .(cohort, protein)]), on=.(cohort, protein), nomatch=0L]
setorder(m, job_id)
if (!nrow(m)) stop('No positive-IV assays are available for the requested cohort')
d[, SNPID := normalized_key]
family <- m[, .(testing_family=.N), by=cohort]
rows <- vector('list', nrow(m))
methods <- c('IVW', 'Median', 'Egger', 'PIVW', 'RAPS')
for (j in seq_len(nrow(m))) {
  assay <- m[j]
  x <- d[cohort == assay$cohort & protein == assay$protein]
  stopifnot(nrow(x) > 0L)
  fit <- rbindlist(lapply(methods, function(method)
    fit_original_method(x, method, assay$bootstrap_seed)))
  fit[, `:=`(job_id=assay$job_id, cohort=assay$cohort, protein=assay$protein,
             protein_symbol=assay$protein_symbol, n_IVs=nrow(x),
             mean_F=mean((x$beta_exposure/x$se_exposure)^2),
             OR=exp(estimate), OR_lower=exp(lower), OR_upper=exp(upper))]
  rows[[j]] <- fit
  if (j %% 50L == 0L) message(j, '/', nrow(m), ' assays')
}
r <- merge(rbindlist(rows), family, by='cohort', sort=FALSE)
# The denominator is all positive-IV assays, including assays for which an
# individual secondary method is unavailable, not just rows with finite P.
r[, `:=`(FDR=p.adjust(p, method='BH', n=unique(testing_family)),
          Bonferroni=p.adjust(p, method='bonferroni', n=unique(testing_family))),
  by=.(cohort, method)]
setorder(r, job_id, method)
setcolorder(r, c('cohort', 'protein', 'protein_symbol', 'job_id', 'n_IVs', 'mean_F',
                'method', 'estimate', 'se', 'lower', 'upper', 'OR', 'OR_lower',
                'OR_upper', 'p', 'FDR', 'Bonferroni', 'testing_family',
                'egger_intercept_p', 'status', 'error'))
if ((file.exists(out) && !dir.exists(out)) ||
    (dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE)))) {
  stop('Output path changed during analysis; existing files will not be overwritten')
}
if (!dir.create(out, recursive=TRUE, showWarnings=FALSE) && !dir.exists(out)) {
  stop('Could not create output directory')
}
fwrite(r, file.path(out, 'MR_results.tsv'), sep='\t', na='NA')
capture.output(sessionInfo(), file=file.path(out, 'R_session.txt'))
message('Saved ', nrow(r), ' method results for ', nrow(m), ' positive-IV assays.')

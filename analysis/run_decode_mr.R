#!/usr/bin/env Rscript
# Raw deCODE protein pQTL file -> numeric cis selection -> allele-key matching
# -> PLINK clumping -> MR. Run on allocated compute, never a login node.
# Set these reference paths for your environment; no MR input cache is used.
# maps: decode_chrN.rds plus Map_chromosome_index.tsv (source -> GRCh37 chr).
# ld: chrN.{bed,bim,fam,exclude}, chrN_clump_lookup.tsv.gz and
# chrN_CLUMP_LABELS_COMPLETE from the verified reference preparation.
config <- list(
  gwas = Sys.getenv('PAN_MR_GWAS', '/path/to/PanScan1234_FinnGenR13.GRCh37.tsv.gz'),
  maps = Sys.getenv('PAN_MR_VARIANT_MAP_DIR', '/path/to/variant_maps'),
  ld = Sys.getenv('PAN_MR_LD_DIR', '/path/to/ld_reference'),
  manifest = Sys.getenv('PAN_MR_MANIFEST', '/path/to/Assay_manifest.tsv'),
  genes = Sys.getenv('DECODE_GENE_COORDS', '/path/to/Gene_coordinates_GRCh37_full_panel.tsv'),
  plink = Sys.getenv('PLINK_BIN', Sys.which('plink'))
)

# Original study calls/defaults; unavailable methods remain NA, never P=1.
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

main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) != 2L) stop('Usage: Rscript analysis/run_decode_mr.R RAW_DECODE_TXT_GZ NEW_OUTPUT_DIR')
  packages <- c('data.table', 'R.utils', 'MendelianRandomization', 'mr.pivw', 'mr.raps')
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly=TRUE)]
  if (length(missing)) stop('Missing R packages: ', paste(missing, collapse=', '))
  suppressPackageStartupMessages(library(data.table))
  setDTthreads(1L)
  Sys.setenv(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
  raw_file <- normalizePath(args[1], mustWork=TRUE)
  out <- args[2]
  if ((file.exists(out) && !dir.exists(out)) ||
      (dir.exists(out) && length(list.files(out, all.files=TRUE, no..=TRUE)))) {
    stop('Output directory must be new or empty; no existing files will be overwritten')
  }
  needed <- c(config$gwas, config$manifest, config$genes,
              file.path(config$maps, 'Map_chromosome_index.tsv'))
  if (!all(file.exists(needed)) || !dir.exists(config$ld)) stop('Set the reference paths in the config header or environment')
  if (!nzchar(config$plink) || !file.exists(config$plink)) stop('Set PLINK_BIN to the study PLINK 1.9 executable')
  for (command in c('gzip', 'awk')) if (!nzchar(Sys.which(command))) stop('Missing executable: ', command)
  assay_id <- sub('\\.txt\\.gz$', '', basename(raw_file))
  if (identical(assay_id, basename(raw_file))) stop('Expected the original deCODE .txt.gz filename')
  manifest <- fread(config$manifest)
  required_manifest <- c('cohort', 'protein', 'protein_symbol', 'job_id', 'bootstrap_seed')
  if (!all(required_manifest %in% names(manifest))) stop('Incomplete assay manifest')
  assay <- manifest[cohort == 'decode' & protein == assay_id]
  if (nrow(assay) != 1L) stop('Raw filename must match exactly one deCODE manifest assay')
  seed <- assay$bootstrap_seed
  if (length(seed) != 1L || !is.finite(seed) || seed < 1 || seed != floor(seed) ||
      seed > .Machine$integer.max) stop('Invalid original bootstrap seed')
  work <- tempfile('decode_mr_'); if (!dir.create(work)) stop('Cannot create temporary directory')
  on.exit(unlink(work, recursive=TRUE), add=TRUE)
  # Reject damaged/truncated sources. No historical partial-gzip recovery is used.
  if (system2('gzip', c('-t', shQuote(raw_file))) != 0L) stop('Raw deCODE gzip integrity check failed')
  raw_info <- file.info(raw_file)[, c('size', 'mtime')]

  # Stream the WHOLE gzip, selecting wanted IDs only after checking every row.
  # pipefail propagates decompression errors rather than accepting a prefix.
  stream_subset <- function(path, key_column, wanted, label) {
    keys <- file.path(work, paste0(label, '_keys.tsv'))
    plain <- file.path(work, paste0(label, '.tsv'))
    fwrite(data.table(key=unique(wanted)), keys, col.names=FALSE, quote=FALSE, sep='\t')
    program <- paste0('NR==FNR { wanted[$1]=1; next } FNR==1 { ',
      'for(i=1;i<=NF;i++) if($i==key_column) k=i; ',
      'if(!k) exit 2; width=NF; print; next } ',
      'NF!=width { exit 3 } ($k in wanted) { print }')
    script <- paste('set -o pipefail; export LC_ALL=C;',
      'gzip -cd --', shQuote(path), '| awk -F', shQuote('\t'),
      '-v', shQuote(paste0('key_column=', key_column)), shQuote(program),
      shQuote(keys), '- >', shQuote(plain))
    before <- file.info(path)[, c('size', 'mtime')]
    status <- system2('/bin/bash', c('-c', shQuote(script)), stderr=file.path(work, paste0(label, '.stderr')))
    if (status != 0L || !identical(before, file.info(path)[, c('size', 'mtime')])) {
      stop('Incomplete, malformed, or changing gzip input: ', path)
    }
    fread(plain)
  }
  empty_iv <- function() data.table(cohort=character(), protein=character(), protein_symbol=character(),
    normalized_key=character(), SNPID=character(), chr=integer(), pos_bp37=integer(),
    effect_allele=character(), other_allele=character(), beta_exposure=numeric(), se_exposure=numeric(),
    beta_outcome=numeric(), se_outcome=numeric(), freq_exposure=numeric(), freq_outcome=numeric(),
    source_minus_log10_pvalue_exposure=numeric(), pvalue_exposure_numeric_used=numeric(),
    pvalue_outcome=numeric(), source_variant_id=character(), source_pos_bp38=integer(),
    source_chr38=integer(), N_exposure=numeric(), Fstat=numeric(), distance_to_gene_bp=numeric())
  finish <- function(iv, reason, conflicting=0L) {
    fits <- rbindlist(lapply(c('IVW', 'Median', 'Egger', 'PIVW', 'RAPS'), function(method)
      fit_original_method(iv, method, seed)))
    fits[, `:=`(cohort='decode', protein=assay$protein, protein_symbol=assay$protein_symbol,
      job_id=assay$job_id, n_IVs=nrow(iv), mean_F=if(nrow(iv)) mean(iv$Fstat) else NA_real_,
      OR=exp(estimate), OR_lower=exp(lower), OR_upper=exp(upper),
      n_conflicting_pQTL_loci_excluded=conflicting, preprocessing_status=reason)]
    # Single-assay output: study-wide FDR/Bonferroni require the full assay panel.
    if (!dir.exists(out) && !dir.create(out, recursive=TRUE)) stop('Cannot create output directory')
    files <- file.path(out, c('MR_results.tsv', 'MR_instruments.tsv.gz'))
    if (any(file.exists(files))) stop('Refusing to overwrite existing MR output')
    if (nrow(iv)) fwrite(iv, files[2], sep='\t', compress='gzip', na='NA') else {
      plain <- file.path(work, 'empty_IVs.tsv'); fwrite(iv, plain, sep='\t', na='NA')
      if (system2('gzip', c('-n', '-c', shQuote(plain)), stdout=files[2]) != 0L) stop('Empty IV gzip export failed')
    }
    fwrite(fits, files[1], sep='\t', na='NA')
    message(reason, '; retained IVs: ', nrow(iv))
    invisible(NULL)
  }

  # GRCh37 gene-body bounds are numeric and strand-independent.
  genes <- fread(config$genes)
  if (!all(c('hgnc_symbol', 'chromosome_name', 'start_position', 'end_position') %in% names(genes))) stop('Invalid GRCh37 gene annotation')
  genes[, chr := suppressWarnings(as.integer(sub('^chr', '', as.character(chromosome_name))))]
  windows <- unique(genes[hgnc_symbol == assay$protein_symbol & chr %in% 1:22,
    .(chr, gene_start=as.numeric(as.character(start_position)), gene_end=as.numeric(as.character(end_position)))])
  if (!nrow(windows)) return(finish(empty_iv(), 'No autosomal target annotation'))
  if (any(!is.finite(windows$gene_start) | !is.finite(windows$gene_end) |
          windows$gene_start < 1 | windows$gene_end < 1 |
          windows$gene_start != floor(windows$gene_start) | windows$gene_end != floor(windows$gene_end))) stop('Invalid numeric gene coordinates')
  windows[, `:=`(cis_start=pmax(1, pmin(gene_start, gene_end)-1000000),
                 cis_end=pmax(gene_start, gene_end)+1000000)]
  in_cis <- function(chr, pos) {
    keep <- rep(FALSE, length(pos))
    for (j in seq_len(nrow(windows))) keep <- keep |
      (chr == windows$chr[j] & pos >= windows$cis_start[j] & pos <= windows$cis_end[j])
    keep & !is.na(keep)
  }

  # deCODE Name/Pos are native GRCh38; Name sequence order is ALT:REF.
  # Frozen, reference-checked maps supply normalized GRCh37 keys and orientation.
  index <- fread(file.path(config$maps, 'Map_chromosome_index.tsv'))
  if (!all(c('cohort', 'source_map_chromosome', 'target_chr37') %in% names(index))) stop('Invalid chromosome mapping index')
  source_chrs <- unique(index[cohort == 'decode' & target_chr37 %in% windows$chr, source_map_chromosome])
  if (!length(source_chrs)) stop('No deCODE source chromosomes mapped to the target interval')
  map_files <- file.path(config$maps, paste0('decode_chr', source_chrs, '.rds'))
  if (!all(file.exists(map_files))) stop('Missing canonical deCODE annotation maps')
  maps <- rbindlist(lapply(map_files, readRDS), fill=TRUE)
  required_maps <- c('row_id', 'normalized_key', 'chr', 'pos', 'ref', 'alt', 'source_ref', 'source_alt',
                     'source_pos38', 'source_chr38',
                     'source_effect', 'source_other', 'source_eaf', 'source_to_alt37', 'alt37_eaf')
  if (!all(required_maps %in% names(maps)) || anyDuplicated(maps$row_id)) stop('Invalid deCODE annotation-map schema or duplicate native IDs')
  maps[, `:=`(chr=as.integer(chr), pos=as.integer(pos))]
  maps <- maps[in_cis(chr, pos)]
  if (!nrow(maps)) stop('No verified source mappings in the gene interval; investigate reference annotation')
  if (anyNA(maps[, ..required_maps]) ||
      any(maps$normalized_key != paste(maps$chr,maps$pos,maps$ref,maps$alt,sep=':')) ||
      any(!maps$source_to_alt37 %in% c(-1,1)) ||
      any(!is.finite(maps$alt37_eaf) | maps$alt37_eaf < 0 | maps$alt37_eaf > 1)) stop('Invalid canonical deCODE allele/orientation map')
  raw <- stream_subset(raw_file, 'Name', maps$row_id, 'raw_decode')
  expected <- c('Chrom', 'Pos', 'Name', 'rsids', 'effectAllele', 'otherAllele', 'Beta', 'Pval', 'minus_log10_pval', 'SE', 'N', 'ImpMAF')
  if (!identical(names(raw), expected)) stop('Unexpected raw deCODE 12-column schema')
  if (!identical(raw_info, file.info(raw_file)[, c('size', 'mtime')])) stop('Raw source changed while reading')
  raw <- unique(raw)
  if (!nrow(raw)) return(finish(empty_iv(), 'No raw deCODE records in the mapped cis interval'))
  if (anyDuplicated(raw$Name)) stop('Conflicting raw records for the same native deCODE Name')
  native <- tstrsplit(raw$Name, ':', fixed=TRUE)
  if (length(native) != 4L || any(paste0('chr', sub('^chr', '', native[[1]])) != raw$Chrom) ||
      any(as.integer(native[[2]]) != raw$Pos)) stop('Raw Name does not match its native GRCh38 coordinates')
  x <- merge(raw, maps, by.x='Name', by.y='row_id', sort=FALSE)
  native <- tstrsplit(x$Name, ':', fixed=TRUE)
  if (any(toupper(native[[3]]) != x$source_alt | toupper(native[[4]]) != x$source_ref) ||
      any(x$Pos != x$source_pos38 | as.integer(sub('^chr','',x$Chrom)) != x$source_chr38)) {
    stop('Native deCODE Name must be GRCh38 chr:position:ALT:REF and match the canonical source annotation')
  }
  same <- toupper(x$effectAllele) == x$source_effect & toupper(x$otherAllele) == x$source_other
  swap <- toupper(x$effectAllele) == x$source_other & toupper(x$otherAllele) == x$source_effect
  if (anyNA(same | swap) || any(!(same | swap))) stop('Raw deCODE alleles disagree with the authoritative annotation')
  sign <- x$source_to_alt37 * ifelse(same, 1, -1)
  if (any(!sign %in% c(-1,1))) stop('Invalid native-to-GRCh37 effect orientation')
  z <- x[, .(normalized_key, chr, pos_bp37=pos, beta_alt37=Beta*sign, se_exposure=SE,
    alt37_eaf, source_minus_log10_pvalue_exposure=minus_log10_pval, pvalue_exposure_numeric_used=Pval,
    source_variant_id=Name, source_pos_bp38=Pos, source_chr38=as.integer(sub('^chr', '', Chrom)),
    N_exposure=N, reference_allele=ref, alternate_allele=alt, ImpMAF)]
  z <- z[in_cis(chr, pos_bp37) & ImpMAF >= .01 & is.finite(beta_alt37) &
           is.finite(se_exposure) & se_exposure > 0 & !is.na(source_minus_log10_pvalue_exposure)]
  if (!nrow(z)) return(finish(empty_iv(), 'No eligible cis variants after source QC'))

  # Read the cleaned meta-GWAS directly; stream only required full-key records.
  g <- stream_subset(config$gwas, 'variant_key', z$normalized_key, 'outcome')
  need_g <- c('variant_key', 'chr', 'pos', 'REF', 'ALT', 'effect_allele', 'other_allele', 'beta', 'se', 'pval', 'eaf', 'rsid')
  if (!all(need_g %in% names(g))) stop('Not the cleaned PanScan1234 + FinnGen R13 GWAS schema')
  if (!nrow(g)) return(finish(empty_iv(), 'No eligible cis variants matched the cleaned outcome'))
  if (anyDuplicated(g$variant_key) || any(g$variant_key != paste(g$chr,g$pos,g$REF,g$ALT,sep=':')) ||
      any(g$effect_allele != g$ALT | g$other_allele != g$REF) ||
      any(!is.finite(g$beta) | !is.finite(g$se) | g$se <= 0)) stop('Invalid, ambiguous, or non-ALT-oriented outcome records')
  outcome <- g[, .(normalized_key=variant_key, chr=as.integer(chr), pos_bp37=as.integer(pos),
    effect_allele, other_allele, beta_outcome=beta, se_outcome=se, freq_outcome=eaf,
    pvalue_outcome=pval, outcome_rsid=rsid)]
  d <- merge(z, outcome, by='normalized_key', sort=FALSE, suffixes=c('', '_outcome'))
  direct <- d$effect_allele == d$alternate_allele & d$other_allele == d$reference_allele
  reverse <- d$effect_allele == d$reference_allele & d$other_allele == d$alternate_allele
  if (any(d$chr != d$chr_outcome | d$pos_bp37 != d$pos_bp37_outcome) ||
      anyNA(direct | reverse) || any(!(direct | reverse))) stop('Outcome/annotation allele-key disagreement')
  d[, `:=`(beta_exposure=ifelse(effect_allele==alternate_allele, beta_alt37, -beta_alt37),
            freq_exposure=ifelse(effect_allele==alternate_allele, alt37_eaf, 1-alt37_eaf), SNPID=normalized_key)]
  dup <- d[duplicated(normalized_key) | duplicated(normalized_key, fromLast=TRUE)]
  audit <- dup[, .(n_associations=uniqueN(paste(signif(beta_exposure,12),signif(se_exposure,12)))), by=normalized_key]
  bad <- audit[n_associations > 1L, normalized_key]
  d <- d[!normalized_key %in% bad]
  setorder(d, source_variant_id)
  d <- d[!duplicated(normalized_key)]
  d[, `:=`(cohort='decode', protein=assay$protein, protein_symbol=assay$protein_symbol,
            Fstat=(beta_exposure/se_exposure)^2, distance_to_gene_bp=Inf)]
  for (j in seq_len(nrow(windows))) {
    distance <- pmax(pmin(windows$gene_start[j],windows$gene_end[j])-d$pos_bp37,
      d$pos_bp37-pmax(windows$gene_start[j],windows$gene_end[j]),0)
    d[, distance_to_gene_bp := pmin(distance_to_gene_bp, ifelse(chr==windows$chr[j], distance, Inf))]
  }
  if (any(!is.finite(d$distance_to_gene_bp) | d$distance_to_gene_bp > 1000000)) stop('Instrument candidate is outside numeric cis bounds')
  candidates <- d[pvalue_exposure_numeric_used <= 5e-8]
  if (!nrow(candidates)) return(finish(empty_iv(), 'No genome-wide significant cis-pQTLs', length(bad)))

  # Full allele keys identify variants. Short PLINK labels are reached ONLY
  # through the allele-verified LD lookup, never by trusting a source rsID.
  if (!dir.exists(out) && !dir.create(out, recursive=TRUE)) stop('Cannot create output directory')
  selected <- character()
  for (ch in unique(candidates$chr)) {
    ref <- file.path(config$ld, paste0('chr', ch))
    required <- paste0(ref, c('.bed','.bim','.fam','.exclude','_clump_lookup.tsv.gz','_CLUMP_LABELS_COMPLETE'))
    if (!all(file.exists(required))) stop('Missing allele-verified, label-matched LD reference for chr', ch)
    lookup <- fread(paste0(ref, '_clump_lookup.tsv.gz'), header=FALSE,
                    col.names=c('chr','pos','normalized_key','clump_id'))
    cc <- merge(candidates[chr==ch], lookup[, .(normalized_key,clump_id)], by='normalized_key', sort=FALSE)
    if (!nrow(cc)) next
    if (anyDuplicated(cc$clump_id) || anyDuplicated(cc$normalized_key)) stop('Ambiguous LD reference lookup')
    setorder(cc, clump_id)
    clump_input <- file.path(work, paste0('clump_chr',ch,'.tsv'))
    fwrite(cc[, .(SNP=clump_id, P=pvalue_exposure_numeric_used)], clump_input, sep='\t', quote=FALSE)
    prefix <- file.path(out, paste0('plink_chr',ch))
    status <- system2(config$plink, c('--bfile',shQuote(ref),'--exclude',shQuote(paste0(ref,'.exclude')),
      '--clump',shQuote(clump_input),'--clump-p1','0.5','--clump-r2','0.1','--clump-kb','2000',
      '--threads','1','--memory','2048','--out',shQuote(prefix)),
      stdout=paste0(prefix,'.stdout'), stderr=paste0(prefix,'.stderr'))
    if (status != 0L) stop('PLINK clumping failed for chr',ch)
    if (file.exists(paste0(prefix,'.clumped'))) {
      ids <- fread(paste0(prefix,'.clumped'))$SNP
      selected <- c(selected, cc[clump_id %in% ids, normalized_key])
    }
  }
  # Keep source-variant ordering for bootstrap reproducibility. F is descriptive;
  # no additional F-statistic filter is imposed by the original analysis.
  iv <- candidates[normalized_key %in% selected]
  finish(iv, if(nrow(iv)) 'Raw-source MR complete' else 'No IVs retained by the LD reference/clumping', length(bad))
}

main()

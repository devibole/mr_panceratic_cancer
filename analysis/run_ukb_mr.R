#!/usr/bin/env Rscript
# Raw UKB-PPP archive -> numeric cis/QC -> GRCh37 allele-key matching ->
# original LD clumping -> MR. Run on allocated compute, not a login node.
# Reference maps/panels are static annotation resources, not previous MR inputs.
# Usage: Rscript analysis/run_ukb_mr.R RAW_UKB_TAR NEW_OUTPUT_DIR

# Set these environment variables, or replace the placeholder paths here.
# VARIANT_MAP_DIR: ukb_chrN.rds, containing source-ID -> normalized-key maps.
# LD_DIR: chrN.{bed,bim,fam,exclude} and chrN_clump_lookup.tsv.gz{,.tbi}.
gwas_file <- Sys.getenv('PAN_MR_GWAS', '/path/to/PanScan1234_FinnGenR13.GRCh37.tsv.gz')
variant_map_dir <- Sys.getenv('PAN_MR_VARIANT_MAP_DIR', '/path/to/variant_maps')
ld_dir <- Sys.getenv('PAN_MR_LD_DIR', '/path/to/allele_verified_EUR_reference')
manifest_file <- Sys.getenv('PAN_MR_MANIFEST', '/path/to/Assay_manifest.tsv')
olink_file <- Sys.getenv('UKB_PROTEIN_MAP', '/path/to/olink_protein_map_3k_v1.tsv')
plink_bin <- Sys.getenv('PLINK_BIN', unname(Sys.which('plink')))

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

main <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (length(args) != 2L) stop('Usage: Rscript analysis/run_ukb_mr.R RAW_UKB_TAR NEW_OUTPUT_DIR')
  packages <- c('data.table','MendelianRandomization','mr.pivw','mr.raps',
                'Rsamtools','GenomicRanges','IRanges')
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if (length(missing)) stop('Missing R packages: ',paste(missing,collapse=', '))
  suppressPackageStartupMessages(library(data.table))
  setDTthreads(1L)
  Sys.setenv(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
  raw_file <- normalizePath(args[1],mustWork=TRUE)
  out <- args[2]
  if ((file.exists(out)&&!dir.exists(out)) ||
      (dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE))))
    stop('Output directory must be new or empty')
  if (!all(file.exists(c(gwas_file,manifest_file,olink_file))) ||
      !dir.exists(variant_map_dir) || !dir.exists(ld_dir))
    stop('Set the reference/input paths in the configuration header')
  if (!nzchar(plink_bin) || !file.exists(plink_bin)) stop('Set PLINK_BIN to the PLINK 1.9 executable')
  protein <- sub('\\.tar$','',basename(raw_file))
  assay_name <- protein
  manifest <- fread(manifest_file)
  if (!all(c('cohort','protein','protein_symbol','bootstrap_seed') %in% names(manifest)))
    stop('Manifest requires cohort, protein, protein_symbol and bootstrap_seed')
  assay <- manifest[cohort=='ukb' & protein==assay_name]
  if (nrow(assay)!=1L || !is.finite(assay$bootstrap_seed) ||
      assay$bootstrap_seed<1 || assay$bootstrap_seed!=floor(assay$bootstrap_seed) ||
      assay$bootstrap_seed>.Machine$integer.max) stop('The raw assay needs exactly one valid saved seed in the manifest')
  if (!dir.create(out,recursive=TRUE,showWarnings=FALSE)&&!dir.exists(out)) stop('Cannot create output directory')
  work <- tempfile('ukb_mr_',tmpdir=tempdir())
  if (!dir.create(work)) stop('Cannot create a private temporary directory')
  on.exit(unlink(work,recursive=TRUE),add=TRUE)
  need <- function(x,fields,label) {
    absent <- setdiff(fields,names(x))
    if (length(absent)) stop(label,' missing columns: ',paste(absent,collapse=', '))
  }
  numeric_bp <- function(x,label) {
    y <- suppressWarnings(as.numeric(as.character(x)))
    if (any(!is.finite(y)|y<1|y!=floor(y))) stop('Invalid numeric coordinates: ',label)
    y
  }
  empty_ivs <- function() data.table(cohort=character(),protein=character(),protein_symbol=character(),
    normalized_key=character(),SNPID=character(),chr=integer(),pos_bp37=integer(),
    effect_allele=character(),other_allele=character(),beta_exposure=numeric(),se_exposure=numeric(),
    beta_outcome=numeric(),se_outcome=numeric(),freq_exposure=numeric(),freq_outcome=numeric(),
    source_minus_log10_pvalue_exposure=numeric(),pvalue_exposure_numeric_used=numeric(),
    source_variant_id=character(),source_pos_bp38=numeric(),Fstat=numeric())
  save_instruments <- function(x) {
    plain <- file.path(work,'MR_instruments.tsv')
    fwrite(x,plain,sep='\t',na='NA')
    if (system2('gzip',c('-n','-c',shQuote(plain)),stdout=file.path(out,'MR_instruments.tsv.gz'))!=0L)
      stop('Could not write the complete instrument gzip')
  }
  finish <- function(ivs,reason) {
    fit <- rbindlist(lapply(c('IVW','Median','Egger','PIVW','RAPS'),function(method)
      fit_original_method(ivs,method,assay$bootstrap_seed)))
    fit[,`:=`(cohort='ukb',protein=protein,protein_symbol=assay$protein_symbol,
      n_IVs=nrow(ivs),mean_F=if(nrow(ivs))mean(ivs$Fstat) else NA_real_,
      OR=exp(estimate),OR_lower=exp(lower),OR_upper=exp(upper),instrument_status=reason)]
    save_instruments(ivs)
    fwrite(fit,file.path(out,'MR_results.tsv'),sep='\t',na='NA')
    message(reason,'; selected IVs: ',nrow(ivs))
  }
  # Always consume the entire gzip stream. A read/decompression error is NOT
  # evidence of no cis-pQTLs. Filtering a large GWAS still checks its full stream.
  read_gzip <- function(path,keys=NULL) {
    filtered <- tempfile('read_',tmpdir=work);err <- paste0(filtered,'.stderr')
    before <- file.info(path)[,c('size','mtime')]
    cmd <- paste('gzip -cd --',shQuote(path))
    if (!is.null(keys)) {
      key_file <- tempfile('keys_',tmpdir=work);writeLines(unique(keys),key_file)
      awk <- 'BEGIN {FS=OFS="\t"} NR==FNR {wanted[$1]=1; next} FNR==1 {for(i=1;i<=NF;i++) if($i=="variant_key") k=i; if(!k) exit 3; print; next} ($k in wanted) {print}'
      cmd <- paste(cmd,'| awk',shQuote(awk),shQuote(key_file),'-')
    }
    script <- paste('set -o pipefail',paste(cmd,'>',shQuote(filtered)),sep='\n')
    status <- system2('/bin/bash',c('-c',shQuote(script)),stderr=err)
    if (status!=0L || !identical(before,file.info(path)[,c('size','mtime')]))
      stop('Incomplete or changed gzip source: ',path,'; ',paste(readLines(err,warn=FALSE),collapse=' | '))
    fread(filtered)
  }

  # 1. Native UKB coordinates and Olink gene annotations are GRCh38.
  ann <- fread(olink_file)
  need(ann,c('UKBPPP_ProteinID','chr','gene_start','gene_end'),'Olink annotation')
  oid <- regmatches(protein,regexpr('OID[0-9]+',protein))
  if (length(oid)!=1L || !nzchar(oid)) stop('Cannot identify the archive Olink assay ID')
  annotation_oid <- vapply(regmatches(ann$UKBPPP_ProteinID,gregexpr('OID[0-9]+',ann$UKBPPP_ProteinID)),
                           function(x)if(length(x)==1L)x else '',character(1))
  windows <- ann[annotation_oid==oid]
  windows[,chr:=suppressWarnings(as.integer(as.character(chr)))];windows <- windows[chr %in% 1:22]
  if (!nrow(windows)) {finish(empty_ivs(),'No autosomal Olink target annotation');return(invisible(NULL))}
  windows[,gene_start:=numeric_bp(gene_start,'Olink gene_start')]
  windows[,gene_end:=numeric_bp(gene_end,'Olink gene_end')]
  windows[,`:=`(cis_start=pmax(1,pmin(gene_start,gene_end)-1000000),
                 cis_end=pmax(gene_start,gene_end)+1000000)]
  native_chrs <- unique(windows$chr)
  members <- utils::untar(raw_file,list=TRUE)
  unsafe <- grepl('^[/\\\\]|^[A-Za-z]:|(^|/)\\.\\.(/|$)|[\r\n]',members)
  if (any(unsafe)) stop('Unsafe path in the raw tar archive')
  selected <- members[grepl(paste0('chr(',paste(native_chrs,collapse='|'),')_'),members)&grepl('\\.gz$',members)]
  if (!length(selected)) stop('Raw archive lacks expected cis-chromosome summary statistics')
  if (utils::untar(raw_file,files=selected,exdir=work)!=0L) stop('Raw tar extraction failed')
  raw <- rbindlist(lapply(file.path(work,selected),read_gzip))
  need(raw,c('ID','N','CHROM','GENPOS','ALLELE0','ALLELE1','BETA','SE','A1FREQ','INFO','LOG10P'),'Raw UKB')
  setkey(raw,ID,N);raw <- raw[!duplicated(ID,fromLast=TRUE)]
  raw[,GENPOS:=numeric_bp(GENPOS,'UKB GENPOS')]
  native_chr <- suppressWarnings(as.integer(as.character(raw$CHROM)))
  if (anyNA(native_chr)) stop('Invalid UKB native chromosomes')
  cis <- rep(FALSE,nrow(raw))
  for (j in seq_len(nrow(windows))) cis <- cis | (native_chr==windows$chr[j]&
    raw$GENPOS>=windows$cis_start[j]&raw$GENPOS<=windows$cis_end[j])
  raw <- raw[cis]
  if (!nrow(raw)) {finish(empty_ivs(),'No raw variants in the numeric cis window');return(invisible(NULL))}

  # 2. Map native source IDs to reference-verified, normalized GRCh37 keys.
  # UKB ID and GENPOS must not be assumed to use the same genome build.
  id_chrs <- suppressWarnings(as.integer(sub(':.*','',raw$ID)))
  if (anyNA(id_chrs)||any(!id_chrs %in% 1:22)) stop('Invalid chromosome in a UKB source ID')
  map_chrs <- unique(c(native_chrs,id_chrs))
  paths <- file.path(variant_map_dir,paste0('ukb_chr',map_chrs,'.rds'))
  if (!all(file.exists(paths))) stop('Missing UKB allele-verified source maps')
  maps <- rbindlist(lapply(paths,readRDS),fill=TRUE)
  need(maps,c('row_id','normalized_key','chr','pos','ref','alt','source_ref','source_alt'),'UKB variant maps')
  if (anyDuplicated(maps$row_id)) stop('Source map has ambiguous raw IDs')
  map_fields <- c('row_id','normalized_key','chr','pos','ref','alt','source_ref','source_alt')
  if (anyNA(maps[, ..map_fields]) ||
      any(maps$normalized_key!=paste(maps$chr,maps$pos,maps$ref,maps$alt,sep=':')))
    stop('Invalid canonical UKB allele-key mapping')
  x <- merge(raw,maps,by.x='ID',by.y='row_id',sort=FALSE)
  same <- toupper(x$ALLELE1)==x$source_alt & toupper(x$ALLELE0)==x$source_ref
  swap <- toupper(x$ALLELE1)==x$source_ref & toupper(x$ALLELE0)==x$source_alt
  if (anyNA(same|swap)||any(!(same|swap))) stop('UKB effect alleles disagree with source-ID alleles')
  sign <- ifelse(same,1,-1)
  z <- x[,.(normalized_key,chr=as.integer(chr),pos_bp37=as.integer(pos),
    beta_alt37=BETA*sign,se_exposure=SE,alt37_eaf=ifelse(sign==1,A1FREQ,1-A1FREQ),
    source_minus_log10_pvalue_exposure=LOG10P,pvalue_exposure_numeric_used=10^(-LOG10P),
    source_variant_id=ID,source_pos_bp38=GENPOS,source_chr38=as.integer(CHROM),
    reference_allele=ref,alternate_allele=alt,INFO)]
  z <- z[is.finite(alt37_eaf)&alt37_eaf>=.01&alt37_eaf<=.99&INFO>.3&
    is.finite(beta_alt37)&is.finite(se_exposure)&se_exposure>0&!is.na(source_minus_log10_pvalue_exposure)]
  if (!nrow(z)) {finish(empty_ivs(),'No allele-verified cis variants after source QC');return(invisible(NULL))}

  # 3. Read the cleaned outcome directly; stream only relevant keys into R.
  g <- read_gzip(gwas_file,unique(z$normalized_key))
  need(g,c('variant_key','chr','pos','REF','ALT','effect_allele','other_allele','beta','se','pval','eaf'),'Cleaned GWAS')
  if (!nrow(g)) {finish(empty_ivs(),'No eligible cis variants matched the cleaned outcome');return(invisible(NULL))}
  if (anyDuplicated(g$variant_key) || any(g$variant_key!=paste(g$chr,g$pos,g$REF,g$ALT,sep=':')) ||
      any(g$effect_allele!=g$ALT|g$other_allele!=g$REF) ||
      any(!is.finite(g$beta)|!is.finite(g$se)|g$se<=0)) stop('Invalid or conflicting cleaned GWAS records')
  outcome <- g[,.(normalized_key=variant_key,chr_outcome=as.integer(chr),pos_outcome=as.integer(pos),
    effect_allele,other_allele,beta_outcome=beta,se_outcome=se,freq_outcome=eaf,pvalue_outcome=pval)]
  d <- merge(z,outcome,by='normalized_key',sort=FALSE)
  if (!nrow(d)) {finish(empty_ivs(),'No eligible cis variants matched the cleaned outcome');return(invisible(NULL))}
  if (any(d$chr!=d$chr_outcome|d$pos_bp37!=d$pos_outcome)) stop('Outcome coordinate mismatch')
  d[,`:=`(beta_exposure=ifelse(effect_allele==alternate_allele,beta_alt37,-beta_alt37),
    freq_exposure=ifelse(effect_allele==alternate_allele,alt37_eaf,1-alt37_eaf),SNPID=normalized_key)]
  duplicate_rows <- d[duplicated(normalized_key)|duplicated(normalized_key,fromLast=TRUE)]
  conflicting <- duplicate_rows[,.(n=uniqueN(paste(signif(beta_exposure,12),signif(se_exposure,12)))),
                               by=normalized_key][n>1L,normalized_key]
  d <- d[!normalized_key %in% conflicting]
  setorder(d,source_variant_id);d <- d[!duplicated(normalized_key)]
  d[,`:=`(cohort='ukb',protein=protein,protein_symbol=assay$protein_symbol,Fstat=(beta_exposure/se_exposure)^2)]
  candidates <- d[pvalue_exposure_numeric_used<=5e-8]

  # 4. Original LD panel/settings and original numeric P-value tie ordering.
  # Long allele keys are bridged to validated short BIM labels, never guessed rsIDs.
  selected_keys <- character()
  for (ch in unique(candidates$chr)) {
    ref <- file.path(ld_dir,paste0('chr',ch))
    required <- paste0(ref,c('.bed','.bim','.fam','.exclude','_clump_lookup.tsv.gz','_clump_lookup.tsv.gz.tbi'))
    if (!all(file.exists(required))) stop('Incomplete allele-verified LD reference for chromosome ',ch)
    cc <- candidates[chr==ch]
    lines <- Rsamtools::scanTabix(paste0(ref,'_clump_lookup.tsv.gz'),param=
      GenomicRanges::GRanges(as.character(ch),IRanges::IRanges(min(cc$pos_bp37),max(cc$pos_bp37))))[[1L]]
    if (!length(lines)) next
    lookup <- fread(text=paste(lines,collapse='\n'),header=FALSE,
                    col.names=c('chr','pos','normalized_key','clump_id'))
    cc <- merge(cc,lookup[,.(normalized_key,clump_id)],by='normalized_key',sort=FALSE)
    if (!nrow(cc)) next
    if (anyDuplicated(cc$clump_id)||anyDuplicated(cc$normalized_key)) stop('Ambiguous LD lookup labels')
    setorder(cc,clump_id)
    clump_input <- file.path(work,paste0('chr',ch,'_clump.tsv'))
    fwrite(cc[,.(SNP=clump_id,P=pvalue_exposure_numeric_used)],clump_input,sep='\t',quote=FALSE)
    prefix <- file.path(work,paste0('chr',ch))
    status <- system2(plink_bin,c('--bfile',shQuote(ref),'--exclude',shQuote(paste0(ref,'.exclude')),
      '--clump',shQuote(clump_input),'--clump-p1','0.5','--clump-r2','0.1','--clump-kb','2000',
      '--threads','1','--memory','2048','--out',shQuote(prefix)),
      stdout=paste0(prefix,'.stdout'),stderr=paste0(prefix,'.stderr'))
    if (status!=0L) stop('PLINK clumping failed for chromosome ',ch)
    if (file.exists(paste0(prefix,'.log'))) file.copy(paste0(prefix,'.log'),
      file.path(out,paste0('Clump_chr',ch,'.log')),overwrite=FALSE)
    if (file.exists(paste0(prefix,'.clumped'))) {
      ids <- fread(paste0(prefix,'.clumped'))$SNP
      selected_keys <- c(selected_keys,cc[clump_id %in% ids,normalized_key])
    }
  }

  # 5. All original MR methods; no new F-statistic filter or assay-level FDR.
  # Cohort-level adjusted P values require collecting the complete positive-IV panel.
  ivs <- candidates[normalized_key %in% selected_keys]
  finish(ivs,if(nrow(ivs))'Raw-source MR complete' else 'No eligible IVs after significance filtering/LD clumping')
}
main()

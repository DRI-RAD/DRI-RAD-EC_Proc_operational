source('run_pipeline.R')
root <- tempfile('uncertainty-schema-'); dir.create(root)
t <- pipeline_time(c('2026-01-01 00:30','2026-01-01 01:00'))
old <- list(data=data.frame(TIMESTAMP=t[1], LE=10),
            units=data.frame(TIMESTAMP=c('TS',NA),LE=c('W m-2','avg')))
for (pattern in pipeline_patterns[c('L3_EC','L4')]) {
  if (pattern == pipeline_patterns[['L4']]) {
    for (name in c('VBR','LE_VBR','LE_VBR_F','LE_VBR_MDS','ET_VBR','ET_VBR_F')) {
      old$data[[name]] <- 99; old$units[[name]] <- c('1','derived')
    }
    old$data$LE_VBR_Ts_MDS <- 7; old$units$LE_VBR_Ts_MDS <- c('W m-2','derived')
  }
  out <- file.path(root, pattern); dir.create(out)
  f <- file.path(out,paste0('ECSM_',pattern,'_old.csv'))
  pipeline_write_table(old$data,old$units,f,'legacy')
  before <- tools::md5sum(f)
  h <- pipeline_history(out,pattern)
  stopifnot(is.na(h$data$LE_rand_err), identical(before,tools::md5sum(f)))
  if (pattern == pipeline_patterns[['L4']]) {
    stopifnot(all(is.na(h$data[pipeline_le_quality_columns()])), h$data$LE_VBR_Ts_MDS==7,
      !any(c('VBR','LE_VBR','LE_VBR_F','LE_VBR_MDS','ET_VBR','ET_VBR_F') %in% names(h$data)),
      identical(h, pipeline_output_schema(h,pattern)))
  }
  fresh <- h; fresh$data$TIMESTAMP <- t[2]; fresh$data$LE_rand_err <- 2
  if (pattern == pipeline_patterns[['L4']]) fresh$data$LE_VBR_Ts_MDS_QC <- 3
  # New stage column order need not match legacy columns after migration.
  order <- c('TIMESTAMP','LE_rand_err',setdiff(names(fresh$data),c('TIMESTAMP','LE_rand_err')))
  fresh$data <- fresh$data[,order]; fresh$units <- fresh$units[,order]
  stage <- file.path(out,'staged'); dir.create(stage)
  pipeline_write_table(fresh$data,fresh$units,file.path(stage,paste0('ECSM_',pattern,'_new.csv')),'new')
  pipeline_level4_figures <- function(...) invisible(NULL)
  w <- list(pending=t[2],read_start=t[2],end=t[2])
  options(ec.pipeline=list(window=w))
  file <- pipeline_publish(stage,out,pattern,w,FALSE)
  saved <- pipeline_read_table(file)
  stopifnot(nrow(saved$data)==2, is.na(saved$data$LE_rand_err[1]), saved$data$LE_rand_err[2]==2,
            identical(names(saved$data),names(saved$units)))
  if (pattern == pipeline_patterns[['L4']]) {
    stopifnot(is.na(saved$data$LE_VBR_Ts_MDS_QC[1]),saved$data$LE_VBR_Ts_MDS_QC[2]==3)
    archive <- list.files(file.path(out,'archive'),pattern='csv$',full.names=TRUE)
    stopifnot(any(vapply(archive,function(p) unname(tools::md5sum(p))==unname(before),logical(1))))
  }
  # Unrelated schema changes are still rejected, not silently padded.
  bad <- saved; bad$data$unrelated <- 1; bad$units$unrelated <- c('1',NA)
  pipeline_write_table(bad$data,bad$units,file.path(out,paste0('ECSM_',pattern,'_bad.csv')),'bad')
  stopifnot(inherits(tryCatch(pipeline_history(out,pattern), error=identity),'error'))
}
options(ec.pipeline=NULL)
unlink(root,recursive=TRUE)
cat('PASS: uncertainty and LE-quality legacy migration, retired VBR removal, archive preservation and strict other schema checks.\n')

# Check launcher settings without reading any network data or running the pipeline.
for(person in c('Larsen','Jasoni')) {
 e <- new.env(parent=baseenv()); e$calls <- list()
 e$source <- function(...) invisible(NULL)
 e$run_pipeline <- function(...) e$calls[[length(e$calls)+1L]] <- list(...)
 for(expr in parse(paste0('run_by_',person,'.R'))) eval(expr,e)
 stopifnot(length(e$calls)==1L, identical(e$calls[[1]]$reprocess,FALSE),
           is.null(e$calls[[1]]$output_dir), length(e$calls[[1]]$sites)==7L,
           startsWith(e$calls[[1]]$base_dir,if(person=='Larsen') 'W:/' else 'Y:/'))
}
# Shared runners resolve the project root from outside the project.
project <- normalizePath('.',winslash='/')
scratch <- tempfile();dir.create(scratch);setwd(scratch)
source(file.path(project,'R','run_level5_ameriflux.R'))
stopifnot(identical(.amf_root,project), file.exists(eval(formals(run_level5_ameriflux)$config_file)))
source(file.path(project,'R','run_pipeline.R'))
stopifnot(identical(.pipeline_root,project),identical(.amf_root,project))
setwd(project);unlink(scratch,recursive=TRUE)
cat('PASS: one incremental run per launcher, correct drives/default outputs/sites, relocated runners from external cwd.\n')

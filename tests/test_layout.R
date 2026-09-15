# Confirm that the operational directory layout resolves independently of cwd.
# Run from the project root: Rscript tests/test_layout.R
project <- normalizePath(".", winslash = "/")
scratch <- tempfile("layout-test-")
dir.create(scratch)
setwd(scratch)
source(file.path(project, "run_pipeline.R"))
stopifnot(identical(.pipeline_root, project))
stopifnot(all(file.exists(file.path(.pipeline_root, pipeline_stages))))
stopifnot(identical(eval(formals(run_pipeline)$sites),
                    c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM")))
config <- eval(formals(run_pipeline)$config_file)
stopifnot(file.exists(config))
dirs <- readr::read_csv(config, show_col_types = FALSE)
options(ec.pipeline = list(dirs = dirs, base_dir = scratch,
                           sites = "ECSM", reference_sites = "ECSM"))
for (stage in c("L1", "L2")) {
  env <- new.env(parent = globalenv())
  # Evaluate startup and metadata reads only; never start a scientific site loop.
  for (expr in parse(file.path(project, pipeline_stages[[stage]]))) {
    if (is.call(expr) && identical(expr[[1]], as.name("for"))) break
    eval(expr, envir = env)
  }
  if (stage == "L1") stopifnot(exists("manual_flag", envir = env), nrow(env$manual_flag) > 0)
  if (stage == "L2") stopifnot(exists("meta_G", envir = env), nrow(env$meta_G) > 0)
}
options(ec.pipeline = NULL)
setwd(project)
unlink(scratch, recursive = TRUE)
cat("PASS: external working directory, default sites, stage paths, and real metadata loading.\n")

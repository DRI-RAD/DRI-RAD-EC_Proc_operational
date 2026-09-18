# Source this file, then call run_level5_ameriflux(). No processing starts on source.
# This exporter reads existing Level 4 snapshots; it does not rerun QAQC or MDS.
.amf_sources <- Filter(Negate(is.null), lapply(sys.frames(), function(x) x$ofile))
.amf_script <- if (length(.amf_sources)) tail(.amf_sources, 1)[[1]] else "run_level5_ameriflux.R"
.amf_root <- dirname(normalizePath(.amf_script, winslash = "/", mustWork = TRUE))
source(file.path(.amf_root, "R", "pipeline_io.R"), encoding = "UTF-8")
source(file.path(.amf_root, "R", "Level_5_AmeriFlux_export.R"), encoding = "UTF-8")

run_level5_ameriflux <- function(
  sites = c("ECDP", "EDVG", "EDVP", "ERVA", "ERVP", "ECSM"),
  base_dir = "Z:/NWI/Task_3_Water_Use_ET_and_Meteoroligical_Monitoring/networks/eddy_stations",
  config_file = file.path(.amf_root, "meta_files", "file_directions_new.csv"),
  start = NULL, end = NULL, output_dir = NULL, dry_run = FALSE,
  unit_overrides = character(), extra_mapping = NULL) {
  # start/end select inclusive LOCAL-STANDARD interval-end observations.
  # NULL bounds mean all saved Level 4 data, without rereading raw/EddyPro files.
  dirs <- readr::read_csv(config_file, show_col_types = FALSE)
  if (!length(sites) || anyNA(sites) || anyDuplicated(sites) || any(!sites %in% dirs$site))
    stop("Select unique sites present in config_file.")
  if (anyDuplicated(dirs$site)) stop("Duplicate site in config_file.")
  if (length(unit_overrides) && (is.null(names(unit_overrides)) || anyDuplicated(names(unit_overrides))))
    stop("unit_overrides must be a named character vector, e.g. c(PA = 'kPa').")
  results <- list()
  for (site in sites) {
    row <- dirs[dirs$site == site, , drop = FALSE]
    site_output <- pipeline_path(base_dir, row$dir_output)
    if (dir.exists(file.path(site_output, ".pipeline-lock"))) stop("Pipeline is writing this site; retry after it finishes: ", site)
    files <- pipeline_files(site_output, "Level_4_post_processed_data")
    if (!length(files)) stop("No Level 4 output for ", site)
    item <- pipeline_history(site_output, "Level_4_post_processed_data")
    prepared <- amf_prepare(item, start = start, end = end,
                            unit_overrides = unit_overrides, extra_mapping = extra_mapping)
    message(site, ": ", nrow(prepared$data), " half-hours; ", ncol(prepared$data) - 2L,
            " variables; ", sum(prepared$audit$status == "omitted"), " unit/mapping omissions.")
    omissions <- prepared$audit[prepared$audit$status == "omitted", c("source", "note"), drop = FALSE]
    if (nrow(omissions)) print(omissions, row.names = FALSE)
    if (dry_run) {
      results[[site]] <- prepared
    } else {
      destination <- if (is.null(output_dir)) file.path(site_output, "AmeriFlux") else file.path(output_dir, site)
      results[[site]] <- amf_write(prepared, site, destination)
      # Record input checksums separately; never modify or rotate Level 4 inputs.
      write.csv(data.frame(file = normalizePath(files, winslash = "/"), md5 = unname(tools::md5sum(files))),
        file.path(dirname(results[[site]]), "metadata", "source_files.csv"), row.names = FALSE)
      message("Saved: ", results[[site]])
    }
  }
  invisible(results)
}

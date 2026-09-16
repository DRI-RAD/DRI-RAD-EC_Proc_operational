# Discover only files belonging to the configured site's LI710 stream.
# Backups are ordered by modification time; the active .dat file wins overlaps.
pipeline_li710_files <- function(row, base_dir) {
  configured <- unlist(lapply(intersect(c("dir_LI710", "dir_LI710_old"), names(row)),
    function(n) pipeline_path(base_dir, row[[n]])), use.names = FALSE)
  configured <- configured[!is.na(configured) & nzchar(configured)]
  folders <- unique(dirname(configured))
  if ("dir_met" %in% names(row) && !is.na(row$dir_met))
    folders <- unique(c(folders, dirname(pipeline_path(base_dir, row$dir_met))))
  candidates <- unique(unlist(lapply(folders, list.files, full.names = TRUE)))
  stems <- unique(tolower(sub("(?i)(li710).*", "\\1", basename(configured), perl = TRUE)))
  stems <- stems[grepl("li710$", stems)]
  names_lower <- tolower(basename(candidates))
  match_stream <- if (length(stems)) Reduce(`|`, lapply(stems, function(s)
    startsWith(names_lower, paste0(s, ".")) | startsWith(names_lower, paste0(s, "_")))) else
    grepl("(^|_)li710[._]", names_lower)
  candidates <- candidates[match_stream & grepl("\\.(dat|backup)$", names_lower)]
  files <- unique(c(configured[file.exists(configured)], candidates))
  files <- files[!dir.exists(files)]
  files <- unique(normalizePath(files, winslash = "/", mustWork = TRUE))
  active <- tolower(basename(files)) %in% paste0(stems, ".dat")
  files[order(active, file.info(files)$mtime, files)]
}

pipeline_li710_read <- function(files, window = NULL) {
  fields <- c("TIMESTAMP", "LE_710", "H_710", "diag", "flow", "tilt", "data_qc")
  empty <- as.data.frame(setNames(lapply(fields, function(n) if (n == "TIMESTAMP")
    as.POSIXct(character(), tz = "UTC") else numeric()), fields))
  data <- empty
  units <- as.data.frame(setNames(lapply(fields, function(n) c(if (n == "TIMESTAMP") "TS" else "unit", "avg")), fields))
  for (file in files) {
    metadata <- readr::read_csv(file, skip = 1, n_max = 2, col_types = readr::cols(.default = "c"))
    required <- setdiff(fields, "data_qc")
    if (!all(required %in% names(metadata))) stop("Missing LI710 fields in ", file, ": ", paste(setdiff(required, names(metadata)), collapse = ", "))
    item <- readr::read_csv(file, skip = 4, col_names = names(metadata),
      col_types = readr::cols(.default = "c"), na = c("", "NA", "NaN", "-9999", "9999999"))
    item$TIMESTAMP <- pipeline_time(item$TIMESTAMP)
    if (!is.null(window)) {
      # Include samples that round down to the requested last half-hour.
      item <- item[item$TIMESTAMP >= window$read_start & item$TIMESTAMP < window$end + 1800, ]
    }
    if (!"data_qc" %in% names(item)) { item$data_qc <- rep("0", nrow(item)); metadata$data_qc <- NA_character_ }
    item <- item[, fields]
    for (name in setdiff(fields, "TIMESTAMP")) {
      converted <- suppressWarnings(as.numeric(item[[name]]))
      if (any(!is.na(item[[name]]) & is.na(converted))) stop("Non-numeric LI710 field ", name, " in ", file)
      item[[name]] <- converted
    }
    data <- dplyr::bind_rows(data, item)
    units <- metadata[, fields]
  }
  # Deduplicate exact sample clocks before averaging to the half-hour grid.
  data <- data[!duplicated(data$TIMESTAMP, fromLast = TRUE), ]
  data <- data[order(data$TIMESTAMP), ]
  list(data = data, units = units, files = files)
}

# Stable optional-input schema keeps EC-only L4 outputs compatible with later
# periods where LI710 exists. Missing measurements and QC remain explicitly NA.
pipeline_empty_li710 <- function(timestamps) {
  fields <- c("LE_710", "H_710", "diag_710", "flow_710", "tilt_710", "data_qc_710")
  # Match the real L3 export's order (all LE flags, then all H flags).
  fields <- c(fields[1:6], paste0("LE_710_", c("QC_diag", "QC_despike", "QC_longrun", "QC_flow")),
              paste0("H_710_", c("QC_diag", "QC_despike", "QC_longrun", "QC_flow")))
  data <- data.frame(TIMESTAMP = timestamps)
  for (n in fields) data[[n]] <- rep(NA_real_, length(timestamps))
  units <- as.data.frame(setNames(lapply(names(data), function(n)
    c(if (n == "TIMESTAMP") "TS" else if (n %in% c("LE_710", "H_710")) "W m-2" else "#", "avg")), names(data)))
  list(data = data, units = units)
}

# Exercise the actual Level 4 export expression without running the pipeline.
# Run from the project root: Rscript tests/test_vbr_mds.R
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
stage_file <- if (length(args)) args[1] else "R/Level_4_gapfilling_ET_calculation_new.R"
expressions <- parse(stage_file)
exports <- list()
find_export <- function(x) {
  if (is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], as.name("MDSout")) &&
    grepl("LE_VBR_Ts_fall", paste(deparse(x), collapse = " "), fixed = TRUE)) {
    exports[[length(exports) + 1L]] <<- x
  }
  if (is.call(x) || is.expression(x)) for (i in seq_along(x)) {
    if (identical(x[[i]], quote(expr = ))) next
    if (!is.symbol(x[[i]])) find_export(x[[i]])
  }
}
find_export(expressions)
stopifnot(length(exports) == 1L)
raw <- data.frame(
  LE_RE_f = 1:5, LE_RE_fall = 11:15, LE_RE_fall_qc = c(1, 2, 1, 3, NA),
  LE_VBR_f = 21:25, LE_VBR_fall = 31:35, LE_VBR_fall_qc = c(1, 2, 3, NA, 1),
  LE_VBR_Ts_f = 41:45, LE_VBR_Ts_fall = 51:55, LE_VBR_Ts_fall_qc = c(2, 1, NA, 1, 3)
)
EProc <- list(sExportResults = function() raw)
EddyData <- data.frame(TIMESTAMP = as.POSIXct("2025-01-01", tz = "UTC") + (1:5) * 1800)
eval(exports[[1]])
stopifnot(
  identical(MDSout$TIMESTAMP, EddyData$TIMESTAMP),
  identical(MDSout$LE_VBR_MDS, c(31L, NA_integer_, NA_integer_, NA_integer_, 35L)),
  identical(MDSout$LE_VBR_Ts_MDS, c(NA_integer_, 52L, NA_integer_, 54L, NA_integer_)),
  identical(MDSout$LE_RE_MDS, c(11L, NA_integer_, 13L, NA_integer_, NA_integer_)),
  identical(MDSout$LE_VBR_F, raw$LE_VBR_f),
  identical(MDSout$LE_VBR_Ts_F, raw$LE_VBR_Ts_f),
  identical(MDSout$LE_RE_F, raw$LE_RE_f)
)
cat("PASS: separate VBR/Ts estimates, independent quality masks, missing QC, and unchanged filled series.\n")

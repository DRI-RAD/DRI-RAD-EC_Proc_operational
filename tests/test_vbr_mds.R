# Exercise actual Level 4 export expressions without running MDS.
# Run from the project root: Rscript tests/test_vbr_mds.R
suppressPackageStartupMessages(library(dplyr))
args <- commandArgs(trailingOnly = TRUE)
stage_file <- if (length(args)) args[1] else "R/Level_4_gapfilling_ET_calculation_new.R"
expressions <- parse(stage_file)
exports <- list()
find_export <- function(x) {
  if (is.call(x) && identical(x[[1]], as.name("<-")) &&
    identical(x[[2]], as.name("MDSout"))) {
    exports[[length(exports) + 1L]] <<- x
  }
  if (is.call(x) || is.expression(x)) for (i in seq_along(x)) {
    if (identical(x[[i]], quote(expr = ))) next
    if (!is.symbol(x[[i]])) find_export(x[[i]])
  }
}
find_export(expressions)
stopifnot(length(exports) == 4L)
raw <- data.frame(
  LE_RE_f = 1:5, LE_RE_fqc = c(0,1,2,3,NA), LE_RE_fall = c(11,12,NA,14,15), LE_RE_fall_qc = c(1, 2, NA, 3, NA),
  LE_VBR_Ts_f = 41:45, LE_VBR_Ts_fqc = c(0,3,2,1,NA), LE_VBR_Ts_fall = 51:55, LE_VBR_Ts_fall_qc = c(2, 1, NA, 1, 3)
)
EProc <- list(sExportResults = function() raw)
EddyData <- data.frame(TIMESTAMP = as.POSIXct("2025-01-01", tz = "UTC") + (1:5) * 1800)
eval(exports[[4]])
stopifnot(
  identical(MDSout$TIMESTAMP, EddyData$TIMESTAMP),
  !any(c('LE_VBR_MDS','LE_VBR_F') %in% names(MDSout)),
  identical(MDSout$LE_VBR_Ts_MDS, raw$LE_VBR_Ts_fall),
  identical(MDSout$LE_RE_MDS, raw$LE_RE_fall),
  identical(MDSout$LE_VBR_Ts_MDS_QC, raw$LE_VBR_Ts_fall_qc),
  identical(MDSout$LE_RE_MDS_QC, raw$LE_RE_fall_qc),
  identical(MDSout$LE_VBR_Ts_F_QC, raw$LE_VBR_Ts_fqc),
  identical(MDSout$LE_RE_F_QC, raw$LE_RE_fqc),
  identical(MDSout$LE_VBR_Ts_F, raw$LE_VBR_Ts_f),
  identical(MDSout$LE_RE_F, raw$LE_RE_f)
)
raw <- data.frame(LE_f=1:5, LE_fqc=c(0,1,2,3,NA), LE_fall=c(11,12,NA,14,15),
                  LE_fall_qc=c(1,2,NA,3,NA), H_f=6:10, NEE_f=11:15)
eval(exports[[2]])
stopifnot(identical(MDSout$LE_MDS,raw$LE_fall), identical(MDSout$LE_MDS_QC,raw$LE_fall_qc),
          identical(MDSout$LE_PI_F,raw$LE_f), identical(MDSout$LE_PI_F_QC,raw$LE_fqc))
mds_li <- raw
eval(exports[[3]])
stopifnot(identical(MDSout$LE_710_MDS,raw$LE_fall), identical(MDSout$LE_710_MDS_QC,raw$LE_fall_qc),
          identical(MDSout$LE_710_PI_F,raw$LE_f), identical(MDSout$LE_710_PI_F_QC,raw$LE_fqc))
mds_li[] <- NA_real_
eval(exports[[3]])
stopifnot(nrow(MDSout)==5, all(is.na(MDSout[-1])))
cat("PASS: EC/LI710/VBR-Ts/RE predictions retain all quality levels, separate QC, missing inputs, unchanged filled series.\n")

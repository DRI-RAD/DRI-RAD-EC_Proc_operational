# Test 100-day windows independently of expensive scientific calculations.
source("run_pipeline.R")
t <- seq(pipeline_time("2026-01-01 00:30"), by = 1800, length.out = 150 * 48)
plan <- pipeline_plan(t, NULL, start = t[140*48+1], end = tail(t, 1))
expanded <- pipeline_mds_window(plan, t)
stopifnot(expanded$read_start == t[50*48+1], expanded$calculation_days == 100,
          identical(plan$pending, expanded$pending), length(expanded$pending) == 10*48)
long <- pipeline_plan(t, NULL)
stopifnot(pipeline_mds_window(long, t)$read_start == t[1])
with_context <- pipeline_plan(t, NULL, start = t[140*48+1], context_days = 110)
stopifnot(pipeline_mds_window(with_context, t)$read_start == t[30*48+1])
error <- tryCatch(pipeline_mds_window(plan, t[-(60*48)]), error = identity)
stopifnot(inherits(error, "error"), grepl("Missing 1 rows", conditionMessage(error)))
error <- tryCatch(pipeline_mds_window(plan, tail(t, 31*48)), error = identity)
stopifnot(inherits(error, "error"), grepl("Process the missing upstream period", conditionMessage(error)))
stopifnot(is.null(pipeline_mds_window(NULL, t)))
# ECDP deployment at 17:30: no artificial rows before real observations.
partial <- t[-(1:34)]
partial_plan <- pipeline_plan(partial, NULL, reprocess=TRUE)
actual <- pipeline_mds_window(partial_plan, partial)
stopifnot(is.null(actual$padding_timestamps),
          identical(actual$pending,partial_plan$pending),actual$read_start==partial[1],
          actual$calculation_days==length(partial)/48)
input <- data.frame(TIMESTAMP=partial, LE=seq_along(partial), processing="EddyPro")
stopifnot(identical(pipeline_slice(input,actual),input))
pipeline_validate_mds(input$TIMESTAMP,100)
# Mid-day end: exactly 100 days, not midnight-rounded extra context.
mid <- pipeline_plan(t,NULL,start=t[140*48+1],end=t[length(t)-13])
mid <- pipeline_mds_window(mid,t)
stopifnot(mid$calculation_days==100,mid$read_start==mid$end-100*86400+1800)
# Interior missing upstream rows still fail; a genuinely short deployment is
# not converted into 100 days by padding nonexistent history.
stopifnot(inherits(tryCatch(pipeline_mds_window(partial_plan,partial[-100]),error=identity),'error'))
short <- tail(t,91*48)
stopifnot(inherits(tryCatch(pipeline_mds_window(pipeline_plan(short,NULL),short,100),error=identity),'error'))
stopifnot(pipeline_mds_window(pipeline_plan(short,NULL),short,90)$calculation_days>=90)
pipeline_validate_mds(tail(t, 100*48), 100)
stopifnot(inherits(tryCatch(pipeline_validate_mds(tail(t, 89*48), 100), error = identity), "error"))
opt <- options(ec.pipeline = list(window = expanded, figure_period = "new"))
df <- data.frame(TIMESTAMP = t, value = seq_along(t))
stopifnot(nrow(pipeline_figure_data(df)) == 480)
options(ec.pipeline = list(window = expanded, figure_period = "calculation"))
stopifnot(nrow(pipeline_figure_data(df)) == 480)
options(opt)
cat("PASS: real 100-day context, mid-day start/end without padding, missing history, regular grid and figure periods.\n")

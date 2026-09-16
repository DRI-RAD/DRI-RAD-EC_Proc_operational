# Regression cases for openeddy's empty/undersized day/night subsets.
source("run_pipeline.R")
n <- 735L
t <- seq(pipeline_time("2026-07-16 17:00"), by = 1800, length.out = n)
sun <- pmax(0, sin(2*pi*(seq_len(n) %% 48)/48))
set.seed(7)
x <- data.frame(timestamp = t, PotRad = sun*500, LE_710 = 50*sun+rnorm(n))
args <- list(var = "LE_710", light = "PotRad", z = 10, var_thr = c(-30, 600), iter = 10)
plain <- do.call(openeddy::despikeLF, c(list(x = x), args))
safe <- do.call(pipeline_li710_despike, c(list(x = x), args))
stopifnot(identical(plain, safe))
missing <- x; missing$LE_710 <- NA_real_
original_error <- tryCatch(do.call(openeddy::despikeLF, c(list(x = missing), args)), error = identity)
stopifnot(inherits(original_error, "error"), grepl("replacement has 1 row", conditionMessage(original_error)))
stopifnot(all(is.na(do.call(pipeline_li710_despike, c(list(x = missing), args)))))
outside <- x; outside$LE_710 <- 1000
stopifnot(all(do.call(pipeline_li710_despike, c(list(x = outside), args)) == 2L))
day <- x; day$LE_710[day$PotRad <= 10] <- NA
day_result <- do.call(pipeline_li710_despike, c(list(x = day), args))
stopifnot(all(is.na(day_result[day$PotRad <= 10])), any(!is.na(day_result[day$PotRad > 10])))
sparse <- x; sparse$LE_710 <- NA_real_; sparse$LE_710[1:5] <- c(1, 2, 3, 4, 1000)
sparse_result <- suppressWarnings(do.call(pipeline_li710_despike, c(list(x = sparse), args)))
stopifnot(all(is.na(sparse_result[1:4])), sparse_result[5] == 2L)
# The installed package and its normal validation must remain untouched.
stopifnot(identical(plain, do.call(openeddy::despikeLF, c(list(x = x), args))))
bad <- x; bad$timestamp[2] <- bad$timestamp[2] + 1
stopifnot(inherits(tryCatch(do.call(pipeline_li710_despike, c(list(x = bad), args)), error = identity), "error"))
cat("PASS: reproduced var_minus error; missing, rejected, daytime-only, sparse and normal inputs; package parity and validation.\n")

source('R/pipeline_io.R')
set.seed(42)
n <- 10*48
t <- pipeline_time('2026-05-01 00:30') + (0:(n-1))*1800
x <- data.frame(TIMESTAMP=t,PotRad=pmax(0,sin(2*pi*(0:(n-1))/48))*500,
                SLE=rnorm(n),SH=rnorm(n),SC=rnorm(n))
x <- x[-(1:34), ] # Actual mid-day deployment; no padding.

before <- x
y <- suppressWarnings(pipeline_storage_qc(x))
for(v in c('SLE','SH','SC')) {
  actual <- x; actual$timestamp <- actual$TIMESTAMP
  flag <- suppressWarnings(pipeline_despike(actual,v,iter=10,light='PotRad',z=10))
  expected <- actual[[v]]; expected[!is.na(flag) & flag>=1] <- NA_real_
  stopifnot(identical(y[[v]],expected))
}
stopifnot(identical(x,before),identical(y$TIMESTAMP,x$TIMESTAMP),identical(y$PotRad,x$PotRad))
# SLE is not a gate for other storage channels.
x$SLE <- NA_real_
z <- suppressWarnings(pipeline_storage_qc(x))
stopifnot(identical(z$SH,y$SH),identical(z$SC,y$SC))
bad <- x;bad$PotRad[100] <- NA_real_
err <- tryCatch(pipeline_storage_qc(bad),error=identity)
stopifnot(inherits(err,'error'),grepl('actual upstream rows',conditionMessage(err)))
x[c('SLE','SH','SC')] <- NA_real_
stopifnot(identical(pipeline_storage_qc(x),x))
invisible(parse('R/Level_4_gapfilling_ET_calculation_new.R'))
cat('PASS: real storage despiking uses actual mid-day observations without padding; all channels, timestamps, missing channels and genuine light gaps handled.\n')

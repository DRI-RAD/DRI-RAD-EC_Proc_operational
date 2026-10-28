# Complete-day aggregation and independent device periods with real plot builds.
suppressPackageStartupMessages({library(ggplot2);library(cowplot)})
source("R/pipeline_level4_figures.R")
t <- as.POSIXct("2026-01-01 00:30:00",tz="UTC")+(0:(6*48-1))*1800
d <- data.frame(TIMESTAMP=t, Rn=100,G=10,LE=60,H=30,
                ET_EC_F=1/48,ET_EC_EBC=0.8/48,
                ET_710_F=NA_real_,ET_710_EBC=NA_real_,P_PI_F=2/48,
                SWC_1_1_1=NA_real_,SWC_1_1_2=NA_real_)
daily <- pipeline_figure_daily(d,c("Rn","G","LE","H"))
stopifnot(nrow(daily)==6,all(daily$complete),all(daily$LE==60))
# Missing row (47), duplicate timestamp despite 48 rows, missing value, and
# off-grid timestamp must each invalidate the corresponding day.
bad <- d[-1,]; bad$TIMESTAMP[49] <- bad$TIMESTAMP[48]
bad$LE[100] <- NA; bad$TIMESTAMP[150] <- bad$TIMESTAMP[150]+1
daily <- pipeline_figure_daily(bad,c("Rn","G","LE","H"))
stopifnot(identical(daily$complete,c(FALSE,FALSE,FALSE,FALSE,TRUE,TRUE)),
          all(is.na(daily$Rn[1:4])))
# Entire missing day stays on the calendar as NA; it must not connect ribbons.
missing_day <- pipeline_figure_daily(d[-(49:96),],c("Rn","G","LE","H"))
stopifnot(nrow(missing_day)==6,!missing_day$complete[2],is.na(missing_day$LE[2]))
keep <- 97:240 # LI710 starts two days later and ends a day earlier than EC.
d$ET_710_F[keep] <- 0.9/48; d$ET_710_EBC[keep] <- 1.1/48
d$ET_EC_F[100] <- NA # middle missing ET interval; cumulative must resume
d$P_PI_F[160] <- NA # same-day comparison excludes this date from all totals
d$SWC_1_1_1[49:nrow(d)] <- 20 # baseline need not exist on first date
ec <- pipeline_water_daily(d,"ET_EC_F","ET_EC_EBC")
li <- pipeline_water_daily(d,"ET_710_F","ET_710_EBC")
stopifnot(nrow(ec)==6,nrow(li)==3,li$DATE[1]==as.Date('2026-01-03'),
          tail(li$DATE,1)==as.Date('2026-01-05'),
          identical(ec$paired,c(TRUE,TRUE,FALSE,FALSE,TRUE,TRUE)),
          is.na(ec$ET_cum[3]),is.na(ec$ET_cum[4]),tail(ec$ET_cum,1)==4,
          tail(ec$P_cum,1)==8,tail(ec$ET_cum_corr,1)==3.2,
          all(ec$cum_min[ec$paired]<=ec$cum_max[ec$paired]),
          !any(is.infinite(as.matrix(ec[vapply(ec,is.numeric,logical(1))]))),
          is.na(ec$SSM[1]),ec$SSM[2]==0)
# No ET is not a zero-flux period; show a placeholder instead of infinite bounds.
empty <- d; empty$ET_710_F <- empty$ET_710_EBC <- NA_real_
stopifnot(nrow(pipeline_water_daily(empty,'ET_710_F','ET_710_EBC'))==0)
artifacts <- Sys.getenv('EC_TEST_ARTIFACT_DIR')
temporary <- !nzchar(artifacts)
if (temporary) artifacts <- tempfile('daily-figures-')
dir.create(artifacts,recursive=TRUE,showWarnings=FALSE)
before <- d
pipeline_water_figure(d,'TEST','EC','ET_EC_F','ET_EC_EBC',file.path(artifacts,'water_EC.jpeg'))
pipeline_water_figure(d,'TEST','LI-710','ET_710_F','ET_710_EBC',file.path(artifacts,'water_LI710.jpeg'))
pipeline_water_figure(empty,'TEST','LI-710','ET_710_F','ET_710_EBC',file.path(artifacts,'water_absent.jpeg'))
pipeline_water_figure(d[1:48,],'TEST','EC','ET_EC_F','ET_EC_EBC',file.path(artifacts,'water_one_day.jpeg'))
pipeline_water_figure(d[1,,drop=FALSE],'TEST','EC','ET_EC_F','ET_EC_EBC',file.path(artifacts,'water_incomplete.jpeg'))
stopifnot(identical(before,d),all(file.info(list.files(artifacts,full.names=TRUE))$size>1000))
if (temporary) unlink(artifacts,recursive=TRUE)
cat('PASS: 48-interval completeness, missing/duplicate days, independent periods, finite ribbons, comparable subtotals, missing SWC and absent/single-day figures.\n')

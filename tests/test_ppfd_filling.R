library(lubridate)
library(dplyr)
source('R/pipeline_io.R')
source('R/Level_5_AmeriFlux_export.R')
expr <- parse('R/Level_2_met_PI_vars_new.R')
# Exercise the very same regression/diurnal functions used by L2.
for(e in expr) if (is.call(e) && identical(e[[1]],as.name('<-')) &&
 as.character(e[[2]])[1] %in% c('fill_long_gap_reg','fill_long_gap_reg_other_sites','fill_mdc_local')) eval(e)
t <- as.POSIXct('2026-06-01',tz='UTC') + seq_len(48*20)*1800
sw <- pmax(0,800*sin(pi*((seq_along(t)%%48)-12)/24))
d <- data.frame(TIMESTAMP=t, SW_IN=sw, PPFD_IN=2*sw)
d$PPFD_IN[300:320] <- NA_real_; d$PPFD_IN_PI_F <- d$PPFD_IN
out <- fill_long_gap_reg(d,'PPFD_IN_PI_F','PPFD_IN','SW_IN',maxgap=4,window=240,min_train=50,force_zero=TRUE)
stopifnot(max(abs(out[300:320]-2*sw[300:320]))<1e-8,identical(out[-(300:320)],d$PPFD_IN[-(300:320)]))
ref <- data.frame(TIMESTAMP=t,PPFD_IN=sw)
out <- fill_long_gap_reg_other_sites(d,'target',list(target=d,other=ref),'PPFD_IN_PI_F','PPFD_IN','PPFD_IN',maxgap=4,window=480,min_train=50,force_zero=TRUE)
stopifnot(max(abs(out[300:320]-2*sw[300:320]))<1e-8)
out <- fill_mdc_local(d,'PPFD_IN_PI_F',timestamp='TIMESTAMP',window_days=5,nmin=3)
stopifnot(all(is.finite(out[300:320])))
# No observations/training cannot invent a filled series.
d$PPFD_IN <- d$PPFD_IN_PI_F <- NA_real_
stopifnot(all(is.na(fill_long_gap_reg(d,'PPFD_IN_PI_F','PPFD_IN','SW_IN',force_zero=TRUE))))
x <- data.frame(TIMESTAMP=t,PPFD_IN=2*sw,PPFD_IN_PI_F=2*sw)
x$PPFD_IN[300:320] <- NA_real_
u <- data.frame(TIMESTAMP=c('TS',''),PPFD_IN=c('umolPhoton m-2 s-1',''),PPFD_IN_PI_F=c('umolPhoton m-2 s-1',''))
for(v in c("LE","H","FC")) {x[[v]]<-0; x[[paste0(v,"_QC_despike")]]<-0; u[[v]]<-c(if(v=="FC") "umolCO2 m-2 s-1" else "W m-2", "")}
p <- amf_prepare(list(data=x,units=u))
stopifnot(identical(p$data$PPFD_IN_F,2*sw),all(is.na(p$data$PPFD_IN[300:320])),tail(names(p$data),1)=='PPFD_IN_F')
old <- pipeline_output_schema(list(data=x['TIMESTAMP'],units=u['TIMESTAMP']),'Level_2_PI_vars_logger_data')
stopifnot(all(is.na(old$data$PPFD_IN_PI_F)),old$units$PPFD_IN_PI_F[1]=='umolPhoton m-2 s-1')
cat('PASS: same-site/other-site/diurnal PPFD filling, observation preservation, no training, photon units, export, legacy schema.\n')



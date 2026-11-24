source('R/pipeline_ppfd_shadow.R')
source('R/Level_5_AmeriFlux_export.R')
day<-rep(1:40,each=48);slot<-rep(1:48,40)
t<-as.POSIXct('2026-06-01',tz='UTC')+seq_along(day)*1800
sw<-pmax(0,800*sin(pi*(slot-12)/24)); sw[slot>36]<-0
d<-data.frame(TIMESTAMP=t,SW_IN=sw,PPFD_IN=2*sw)
five<-day%in%10:14 & slot==24
four<-day%in%20:23 & slot==24
broken<-day%in%c(30,31,33,34,35) & slot==24
d$PPFD_IN[five|four|broken]<-d$PPFD_IN[five|four|broken]*.5
z<-pipeline_ppfd_shadow(d)
stopifnot(abs(z$slope-.5)<1e-10,sum(z$removed)==5,all(z$removed[five]),
 all(!z$removed[four|broken]),all(is.na(z$data$PPFD_IN[five])),
 identical(d$SW_IN,z$data$SW_IN),identical(d$PPFD_IN,z$data$PPFD_IN_pre_shadow),
 all(z$data$PPFD_IN_SW_scaled[five]==d$SW_IN[five]*.5))
# Whole missing day breaks a run, not merely a row of another time slot.
gap<-d[day!=12,];stopifnot(!any(pipeline_ppfd_shadow(gap)$removed))
# Shared cloud reductions and constant calibration offsets are not shadows.
cloud<-d;cloud$SW_IN<-cloud$SW_IN*.5;cloud$PPFD_IN<-cloud$SW_IN*3
stopifnot(!any(pipeline_ppfd_shadow(cloud)$removed))
short<-pipeline_ppfd_shadow(d[1:48,]);stopifnot(all(is.na(short$data$PPFD_IN_QC_shadow)),identical(short$data$PPFD_IN,d$PPFD_IN[1:48]))
empty<-d;empty$PPFD_IN<-NA_real_;stopifnot(!any(pipeline_ppfd_shadow(empty)$removed))
increment<-pipeline_ppfd_shadow(d[day==14,],context=d[day<14,])
stopifnot(nrow(increment$data)==48,sum(increment$removed)==1,
 identical(increment$data$SW_IN,d$SW_IN[day==14]))
# 20% alone without 80 W/m2 absolute deficit must not trigger.
small<-d;small$SW_IN<-150;small$PPFD_IN<-300
small$PPFD_IN[five]<-225;stopifnot(!any(pipeline_ppfd_shadow(small)$removed))
cols<-c('TIMESTAMP_START','TIMESTAMP_END','LE_1_1_1','TA_1_1_10','RH_1_1_1',
 'TA_PI_F','TA_1_1_2','TA_1_1_1','LE_PI_F_1_1_1','LE_1_1_2')
o<-amf_column_order(cols)
stopifnot(identical(o[grep('^TA_',o)],c('TA_1_1_1','TA_1_1_2','TA_1_1_10','TA_PI_F')),
 all(diff(grep('^TA_[0-9]',o))==1),setequal(cols,o))
invisible(parse('R/Level_1_met_QAQC_new.R'))
cat('PASS: five-day/clock persistence, four-day and missing-day negatives, calibration/cloud controls, insufficient data, raw/reference preservation, numeric column order and stage syntax.\n')

# Conservative absolute threshold: a 70 W/m2 deficit survives, 90 is rejected.
for (deficit in c(70,90)) {
 a<-d; a$SW_IN<-300; a$PPFD_IN<-600
 a$PPFD_IN[five]<-2*(300-deficit)
 stopifnot(sum(pipeline_ppfd_shadow(a)$removed)==if(deficit==70) 0 else 5)
}

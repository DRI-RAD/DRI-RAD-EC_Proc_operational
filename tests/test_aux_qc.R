source('R/pipeline_aux_qc.R')
clock <- as.POSIXct('2026-01-01',tz='UTC') + seq_len(40)*1800
x <- data.frame(TIMESTAMP=clock,TA_710=c(rep(20,16),seq_len(24)),
 RH_710=seq(20,59),PA_710=seq(90,93.9,.1),WS_05103_5s_Max=c(seq(1,38),-1,1000))
u<-data.frame(TIMESTAMP=c('timestamp','derived'),TA_710=c('deg C','sensor'),
 RH_710=c('%','sensor'),PA_710=c('kPa','sensor'),WS_05103_5s_Max=c('m/s','sensor'))
q<-pipeline_aux_qc(x,u,'L1')
stopifnot(identical(x$TA_710,q$data$TA_710_raw),all(is.na(q$data$TA_710[1:16])),
 all(q$data$TA_710_QC_aux[1:16]==2),all(is.na(tail(q$data$WS_05103_5s_Max,2))),
 identical(names(q$data),names(q$units)),identical(q$data$PA_710,x$PA_710))
# Fifteen repeats survive; physical limits include their boundaries.
x$TA_710<-c(rep(10,15),-60,60,61,seq_len(22));x$RH_710[1:4]<-c(0,110,-1,111)
x$PA_710[1:4]<-c(50,110,49,111)
q<-pipeline_aux_qc(x,u,'L1')$data
stopifnot(all(q$TA_710[1:15]==10),q$TA_710[16]==-60,q$TA_710[17]==60,is.na(q$TA_710[18]),
 all(is.na(q$RH_710[3:4])),all(is.na(q$PA_710[3:4])))
# A missing timestamp/value breaks a constant run; units are not guessed.
x$TA_710<-c(rep(10,16),seq_len(24));x$TIMESTAMP[9:40]<-x$TIMESTAMP[9:40]+1800
stopifnot(all(!is.na(pipeline_aux_qc(x,u,'L1')$data$TA_710[1:16])))
x$TIMESTAMP<-clock;x$TA_710[8]<-NA_real_
stopifnot(sum(is.na(pipeline_aux_qc(x,u,'L1')$data$TA_710[1:16]))==1)
x$PA_710<-x$PA_710*1000;u$PA_710[1]<-'Pa'
stopifnot(pipeline_aux_qc(x,u,'L1')$data$PA_710[1]==50000)
bad<-u;bad$PA_710[1]<-'%'
stopifnot(inherits(tryCatch(pipeline_aux_qc(x,bad,'L1'),error=identity),'error'))
# Sparse, all missing, nonfinite and zero-IQR inputs are safe.
for(vals in list(rep(NA_real_,40),rep(2,40),c(Inf,rep(NA_real_,39)))) {
 x$WS_05103_5s_Max<-vals;pipeline_aux_qc(x,u,'L1')
}
z<-data.frame(TIMESTAMP=clock,CO2=c(-1,0,1e6,rep(400,37)),H2O=c(-9999,Inf,rep(10,38)),
 CO2_SIGMA=c(-1,rep(2,39)),H2O_SIGMA=c(-1,rep(3,39)),WS_MAX=c(seq(1,38),-1,1000),LE=seq_len(40))
zu<-as.data.frame(lapply(z,function(v)c('1','sensor')))
q<-pipeline_aux_qc(z,zu,'L3')
stopifnot(is.na(q$data$CO2[1]),q$data$CO2[2]==0,q$data$CO2[3]==1e6,
 all(is.na(q$data$H2O[1:2])),is.na(q$data$CO2_SIGMA[1]),
 all(is.na(tail(q$data$WS_MAX,2))),identical(q$data$LE,z$LE),identical(names(q$data),names(q$units)))
for(f in c('R/Level_1_met_QAQC_new.R','R/Level_3_EC_data_QAQC_new.R')) invisible(parse(f))
source('R/pipeline_io.R')
old<-list(data=x,units=u)
migrated<-pipeline_output_schema(old,'Level_1_QAQCed_logger_data')
stopifnot(identical(migrated$data$TA_710,x$TA_710),all(is.na(migrated$data$TA_710_QC_aux)),
 setequal(names(migrated$data),names(pipeline_aux_qc(x,u,'L1')$data)))
cat('PASS: thresholds, run boundaries/gaps, unit conversion, sparse inputs, raw/flag/unit preservation, minimal gas checks, unchanged flux and stage syntax.\n')

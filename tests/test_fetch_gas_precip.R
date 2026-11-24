source('run_level5_ameriflux.R')
source('R/pipeline_aux_qc.R')
stopifnot(identical(pipeline_fetch_export_flag(c(0,2,1,NA_real_)),c(1,0,0,NA_real_)))
raw_flag<-openeddy::fetch_filter(data.frame(FETCH_70=c(10,100,NA),WD=c(0,0,0)),
 'FETCH_70','WD',rep(50,36))
stopifnot(identical(pipeline_fetch_export_flag(raw_flag),c(1,0,NA_real_)))
t<-pipeline_time('2026-01-01 00:30:00')+0:3*1800
x<-data.frame(TIMESTAMP=t,LE=100,H=10,FC=-1,WD=c(10,20,NA,30),FETCH_70=50,FETCH_FILTER=c(1,0,NA,NA),
 LE_PI_F=110,H_PI_F=20,FC_PI_F=-2,LE_PI_CORR=120,H_PI_CORR=30)
for(v in c('LE','H','FC')) {
 x[[paste0(v,'_QC_fetch')]]<-c(0,2,0,NA_real_)
 x[[paste0(v,'_QC_despike')]]<-c(0,0,0,1)
}
u<-as.data.frame(lapply(names(x),function(n)c(if(n=='TIMESTAMP') 'TS' else if(n=='WD') 'Decimal degrees' else if(n=='FETCH_70') 'm' else if(startsWith(n,'FC')) 'umolCO2 m-2 s-1' else 'W m-2','avg')));names(u)<-names(x)
item<-list(data=x,units=u)
for(site in c('EDVG','ERVA')) {
 y<-amf_prepare(item,site=site)$data
 stopifnot(identical(y$FETCH_FILTER,c(1,0,NA_real_,NA_real_)),
  y$LE_1_1_1[2]==100,y$FC[2]==-1,is.na(y$LE_1_1_1[4]),
  identical(y$LE_PI_F_1_1_1,x$LE_PI_F),identical(y$LE_PI_CORR_1_1_1,x$LE_PI_CORR))
}
stopifnot(!'FETCH_FILTER' %in% names(amf_prepare(item,site='ECDP')$data))
# Variable sample tests for k=10 vs k=4 and all-missing/zero-IQR guards.
z<-data.frame(TIMESTAMP=t[1]+seq_len(102)*1800,CO2=c(1:100,500,1e6))
zu<-data.frame(TIMESTAMP=c('TS','avg'),CO2=c('ppm','avg'))
a<-pipeline_aux_qc(z,zu,'L3');b<-pipeline_aux_qc(z,zu,'L3',gas_iqr_k=4)
stopifnot(a$data$CO2[101]==500,is.na(b$data$CO2[101]),is.na(a$data$CO2[102]),
 a$data$CO2_QC_aux[102]==4,identical(a$data$CO2_raw,z$CO2))
z$CO2<-NA_real_;stopifnot(all(is.na(pipeline_aux_qc(z,zu,'L3')$data$CO2)))
z$CO2<-5;stopifnot(all(pipeline_aux_qc(z,zu,'L3')$data$CO2==5))
p<-data.frame(TIMESTAMP=t[1]+0:5*1800,P=c(0,500,50000,-1,NA,Inf),P_Bulk=c(NA,1,501,0,500,2))
pu<-data.frame(TIMESTAMP=c('TS','avg'),P=c('mm','avg'),P_Bulk=c('mm','avg'))
pq<-pipeline_aux_qc(p,pu,'L1')
stopifnot(identical(pq$data$P,c(0,500,NA_real_,NA_real_,NA_real_,NA_real_)),
 is.na(pq$data$P_Bulk[3]),identical(pq$data$P_raw,p$P),pq$data$P_QC_aux[3]==1,
 identical(pipeline_precip_range(c(0,500,501,-1,Inf)),c(0,500,NA_real_,NA_real_,NA_real_)))
source('R/pipeline_io.R')
hist<-pipeline_output_schema(list(data=p,units=pu),'Level_1_QAQCed_logger_data')
stopifnot(setequal(names(hist$data),names(pq$data)))
for(f in c('R/Level_1_met_QAQC_new.R','R/Level_2_met_PI_vars_new.R','R/Level_3_EC_data_QAQC_new.R')) invisible(parse(f))
# Exercise the actual existing same-site regression with screened tipping data.
library(dplyr)
library(zoo)
expr<-parse('R/Level_2_met_PI_vars_new.R')
fun<-Filter(function(e)is.call(e) && identical(e[[1]],as.name('<-')) &&
 identical(e[[2]],as.name('fill_long_gap_reg')),as.list(expr))
stopifnot(length(fun)==1);eval(fun[[1]])
rain<-data.frame(TIMESTAMP=t[1]+0:199*1800,P=rep(1:20,10))
rain$P_PI_F<-2*rain$P;rain$P_PI_F[190:200]<-NA_real_
rain$P[200]<-50000;rain$P<-pipeline_precip_range(rain$P)
filled<-fill_long_gap_reg(rain,y_fill='P_PI_F',y_ref='P_PI_F',x_vars='P',
 maxgap=0,window=4800,min_train=50,force_zero=TRUE)
stopifnot(all(abs(filled[190:199]-2*rain$P[190:199])<1e-8),is.na(filled[200]),
 identical(rain$P_PI_F[1:189],filled[1:189]))
cat('PASS: fetch polarity/site scope/unmasked observed flux and copied final products; gas range; precipitation boundaries, raw provenance and historical schema.\n')

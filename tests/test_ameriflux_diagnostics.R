source('R/run_level5_ameriflux.R')
fail<-function(expr)stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
root<-tempfile('amf-diagnostics-');dir.create(root)
t<-seq(as.POSIXct('2024-01-31 23:30:00',tz='UTC'),by=1800,length.out=6)
d<-data.frame(TIMESTAMP_START=format(t-1800,'%Y%m%d%H%M',tz='UTC'),
 TIMESTAMP_END=format(t,'%Y%m%d%H%M',tz='UTC'),LE=c(1,-9999,3,10000,5,6),
 ALL_NA=-9999,CONSTANT=2,LE_MDS_QC=c(0,1,2,3,-9999,1))
p<-file.path(root,'ECSM_Level_5_ameriflux_ready.csv')
write.csv(d,p,row.names=FALSE,quote=FALSE)
dir.create(file.path(root,'metadata'))
write.csv(data.frame(source=names(d)[-(1:2)],target=names(d)[-(1:2)],status='exported',output_unit='W m-2',note='<script>alert(1)</script>'),file.path(root,'metadata/variable_mapping.csv'),row.names=FALSE)
before<-tools::md5sum(p)
report<-plot_level5_ameriflux(p)
summary<-read.csv(file.path(root,'diagnostics/variable_summary.csv'))
stopifnot(nrow(summary)==4,all(summary$status=='ok'),summary$missing[1]==1,
 summary$max[1]==10000,summary$missing_percent[2]==100,
 summary$min[3]==2,summary$max[3]==2,
 length(list.files(file.path(root,'diagnostics/figures'),pattern='[.]png$'))==4,
 identical(before,tools::md5sum(p)))
html<-paste(readLines(report),collapse='\n')
stopifnot(grepl('&lt;script&gt;',html,fixed=TRUE),!grepl('<script>alert(1)',html,fixed=TRUE),
 grepl('id="search"',html,fixed=TRUE),grepl('loading="lazy"',html,fixed=TRUE))
single<-file.path(root,'one.csv');write.csv(d[1,],single,row.names=FALSE,quote=FALSE)
plot_level5_ameriflux(single,output_dir=file.path(root,'one-diagnostics'))
bad<-d;bad$TIMESTAMP_END[2]<-bad$TIMESTAMP_END[1]
write.csv(bad,file.path(root,'bad.csv'),row.names=FALSE,quote=FALSE)
fail(plot_level5_ameriflux(file.path(root,'bad.csv')))
stopifnot(identical(formals(run_level5_ameriflux)$make_plots,TRUE))
fail(run_level5_ameriflux(make_plots=NA))
source('R/run_pipeline.R')
fail(run_pipeline(ameriflux_plots=NA))
stopifnot(identical(formals(run_pipeline)$ameriflux_plots,TRUE))
# Confirm opt-out, dry run, default integration, and preservation on renderer failure.
site_dir<-file.path(root,'L4');dir.create(site_dir)
x<-data.frame(TIMESTAMP=t,LE=1,H=2,FC=3)
for(v in c('LE','H','FC'))for(s in c('_QC_fetch','_QC_despike'))x[[paste0(v,s)]]<-0
units<-as.data.frame(as.list(setNames(rep('1',ncol(x)),names(x))))
units$LE<-units$H<-'W m-2';units$FC<-'umolCO2 m-2 s-1'
pipeline_write_table(x,rbind(units,units),file.path(site_dir,'ECSM_Level_4_post_processed_data_test.csv'),'test')
config<-file.path(root,'config.csv');write.csv(data.frame(site='ECSM',dir_output=site_dir),config,row.names=FALSE)
args<-list(sites='ECSM',base_dir=root,config_file=config)
preview<-do.call(run_level5_ameriflux,c(args,list(dry_run=TRUE)))
stopifnot(!dir.exists(file.path(site_dir,'AmeriFlux')))
off<-do.call(run_level5_ameriflux,c(args,list(make_plots=FALSE)))
stopifnot(file.exists(off$ECSM),!dir.exists(file.path(dirname(off$ECSM),'diagnostics')))
on<-do.call(run_level5_ameriflux,args)
stopifnot(file.exists(attr(on,'diagnostics')$ECSM))
original_plotter<-plot_level5_ameriflux
plot_level5_ameriflux<-function(...)stop('Simulated graphics failure')
failed<-suppressWarnings(do.call(run_level5_ameriflux,args))
plot_level5_ameriflux<-original_plotter
stopifnot(file.exists(failed$ECSM),is.na(attr(failed,'diagnostics')$ECSM),
 file.exists(file.path(dirname(failed$ECSM),'diagnostic_error.txt')))
unlink(root,recursive=TRUE)
cat('PASS: every submitted variable, all-missing/constant/extreme/flag/single-row inputs, summary reconciliation, safe HTML, invalid clocks and source immutability.\n')

source('run_pipeline.R')
fail<-function(expr) stopifnot(inherits(tryCatch(force(expr),error=identity),'error'))
stopifnot(is.null(pipeline_output_root(NULL)))
for(v in list(NA_character_,'',c('a','b'),TRUE,'C:relative')) fail(pipeline_output_root(v))
root<-tempfile('output-root-');dir.create(root)
original<-file.path(root,'original');dir.create(original)
custom<-file.path(root,'custom results');dir.create(custom)
site<-file.path(custom,'ECSM');dir.create(site)
t<-pipeline_time('2026-01-01 00:30:00')+0:3*1800
d<-data.frame(TIMESTAMP=t,LE=100,H=10,FC=-1)
for(v in c('LE','H','FC')) for(s in c('_QC_despike','_QC_fetch')) d[[paste0(v,s)]]<-0
u<-as.data.frame(lapply(names(d),function(v)c(if(v=='TIMESTAMP') 'TS' else if(v=='FC') 'umolCO2 m-2 s-1' else if(grepl('QC',v)) '1' else 'W m-2','avg')))
names(u)<-names(d)
for(stage in c('L1','L2','L3_EC','L4')) pipeline_write_table(d,u,
 file.path(site,paste0('ECSM_',pipeline_patterns[[stage]],'_test.csv')),'fixture')
old<-d;old$LE<-999
pipeline_write_table(old,u,file.path(original,'ECSM_Level_4_post_processed_data_test.csv'),'fixture')
config<-file.path(root,'config.csv')
write.csv(data.frame(site='ECSM',dir_output=original,dir_eddypro=NA_character_),config,row.names=FALSE)
hash<-tools::md5sum(c(config,list.files(original,full.names=TRUE)))
args<-list(sites='ECSM',stages='L4',base_dir=root,config_file=config,
 start=t[1],end=t[4],output_dir=custom,ameriflux_plots=FALSE)
preview<-do.call(run_pipeline,c(args,list(dry_run=TRUE)))
stopifnot(!dir.exists(file.path(site,'AmeriFlux')))
result<-do.call(run_pipeline,args)
out<-attr(result,'ameriflux')$ECSM
stopifnot(startsWith(normalizePath(out,winslash='/'),normalizePath(custom,winslash='/')),
 all(read.csv(out)$LE_1_1_1==100),!dir.exists(file.path(original,'AmeriFlux')),
 identical(hash,tools::md5sum(names(hash))))
# Standalone L5 and defaults read their respective histories.
stopifnot(all(run_level5_ameriflux('ECSM',root,config,dry_run=TRUE,processed_dir=custom)$ECSM$data$LE_1_1_1==100),
 all(run_level5_ameriflux('ECSM',root,config,dry_run=TRUE)$ECSM$data$LE_1_1_1==999))
caller<-getwd();setwd(root)
stopifnot(identical(pipeline_output_root('custom results'),normalizePath(custom,winslash='/')))
setwd(caller)
missing<-file.path(root,'not-created')
args$output_dir<-missing
fail(do.call(run_pipeline,c(args,list(dry_run=TRUE))))
stopifnot(!dir.exists(missing))
unlink(root,recursive=TRUE)
cat('PASS: custom history and L5 routing, unchanged defaults/config/old outputs, spaces, relative paths, invalid inputs and no dry-run writes.\n')

# Synthetic orchestration plus actual LI710 QAQC; no network data are modified.
source("run_pipeline.R")
root <- tempfile("end-alignment-"); dir.create(root)
real_root <- .pipeline_root
t <- seq(pipeline_time("2026-08-01 00:30"), by = 1800, length.out = 8)
units_for <- function(x) as.data.frame(setNames(lapply(names(x), function(n) c("unit", "avg")), names(x)))
ep_fields <- c("date", "time", "LE", "qc_LE", "H", "qc_H", "co2_flux", "qc_co2_flux", "Tau", "qc_Tau",
 "co2_mole_fraction", "co2_molar_density", "co2_var", "h2o_mole_fraction", "h2o_molar_density", "h2o_var",
 "ts_var", "u_var", "v_var", "w_var", "sonic_temperature", "wind_speed", "max_wind_speed", "wind_dir",
 "u*", "(z-d)/L", "L", "x_90%", "x_70%", "x_50%", "x_30%", "x_10%")
write_ep <- function(folder, times) {
 file <- file.path(folder, "Aug_2026/output/eddypro_site_full_output_2026-09-18T010000_adv.csv")
 dir.create(dirname(file), recursive=TRUE,showWarnings=FALSE)
 x <- as.data.frame(setNames(lapply(ep_fields,function(n)rep("1",length(times))),ep_fields),check.names=FALSE)
 x$date <- format(times,"%Y-%m-%d");x$time<-format(times,"%H:%M")
 x$LE[length(times)]<-"-9999" # A missing flux still has a valid observation clock.
 writeLines("category",file)
 write.table(matrix(names(x),nrow=1),file,append=TRUE,sep=",",row.names=FALSE,col.names=FALSE)
 write.table(matrix(rep("unit",ncol(x)),nrow=1),file,append=TRUE,sep=",",row.names=FALSE,col.names=FALSE)
 readr::write_csv(x,file,append=TRUE,col_names=FALSE)
}
cfg <- data.frame(site=c("ECSM","ECDP"),dir_met=file.path(root,c("one/met.dat","two/met.dat")),
 dir_output=file.path(root,c("one/output","two/output")),dir_eddypro=file.path(root,c("one/ep","two/ep")),
 dir_LI710=file.path(root,c("one/Site_LI710.dat","two/Site_LI710.dat")),dir_LI710_old=NA_character_)
for(i in 1:2) {
 dir.create(dirname(cfg$dir_met[i]),recursive=TRUE)
 logger_times<-if(i==1)t else t[1:5]
 x<-data.frame(TIMESTAMP=logger_times,PotRad=100)
 pipeline_write_table(x,units_for(x),cfg$dir_met[i],"logger")
 li_times<-if(i==1)t[1:2] else t[1:7]
 li<-data.frame(TIMESTAMP=li_times,LE_710=10,H_710=20,diag=0,flow=200,tilt=0)
 pipeline_write_table(li,units_for(li),cfg$dir_LI710[i],"LI710")
 write_ep(cfg$dir_eddypro[i],if(i==1)t[1:4] else t[1:6])
}
config<-file.path(root,"config.csv");readr::write_csv(cfg,config)
fixture<-file.path(root,"scripts");dir.create(file.path(fixture,"R"),recursive=TRUE)
.pipeline_root<-fixture
real_mds<-pipeline_mds_window;pipeline_mds_window<-function(window,available,min_days)window
real_fig<-pipeline_level4_figures;pipeline_level4_figures<-function(data,site,folder)invisible(NULL)
for(stage in names(pipeline_stages)) {
 code<-c('o<-getOption("ec.pipeline"); site<-o$sites',
 'folder<-o$dirs$dir_output[o$dirs$site==site]',
 'x<-data.frame(TIMESTAMP=o$window$pending,PotRad=100)',
 'u<-data.frame(TIMESTAMP=c("TS","avg"),PotRad=c("W m-2","avg"))')
 if(stage=='L3_LI710') code<-c(code,
 'empty<-pipeline_empty_li710(o$window$pending);x<-empty$data;u<-empty$units',
 'j<-match(x$TIMESTAMP,o$li710$data$TIMESTAMP);x$LE_710<-o$li710$data$LE_710[j]')
 code<-c(code,paste0('pipeline_write(x,u,file.path(folder,paste0(site,"_',pipeline_patterns[[stage]],'_fixture.csv")))'))
 writeLines(code,file.path(fixture,pipeline_stages[[stage]]))
}
args<-list(sites=cfg$site,reference_sites=cfg$site,base_dir=root,config_file=config)
preview<-do.call(run_pipeline,c(args,list(dry_run=TRUE)))
for(site in cfg$site) {
 end<-if(site=='ECSM')t[4] else t[6]
 plans<-preview[startsWith(names(preview),paste0(site,":"))]
 stopifnot(all(vapply(plans,function(w)w$end==end,logical(1))),
           length(unique(vapply(plans,pipeline_period_label,character(1))))==1)
}
stopifnot(!any(dir.exists(cfg$dir_output)))
do.call(run_pipeline,args)
for(site in cfg$site) {
 i<-match(site,cfg$site);end<-if(site=='ECSM')t[4] else t[6]
 for(pattern in pipeline_patterns)stopifnot(max(pipeline_history(cfg$dir_output[i],pattern)$data$TIMESTAMP)==end)
 stopifnot(length(list.dirs(file.path(cfg$dir_output[i],"figures"),recursive=FALSE))==1)
}
li<-pipeline_history(cfg$dir_output[1],pipeline_patterns[['L3_LI710']])$data
stopifnot(nrow(li)==4,all(is.na(li$LE_710[3:4])))
stopifnot(all(unlist(do.call(run_pipeline,args))=='skipped'))
# Reprocessing can omit either or both bounds, even with complete saved history.
for (bounds in list(list(), list(start=t[2]), list(end=t[3]), list(start=t[2],end=t[3]))) {
 p<-do.call(run_pipeline,c(args,bounds,list(reprocess=TRUE,dry_run=TRUE)))
 for(site in cfg$site) {
   expected_end<-if(!is.null(bounds$end))bounds$end else if(site=='ECSM')t[4] else t[6]
   for(stage in names(pipeline_stages)) {
     w<-p[[paste(site,stage,sep=':')]]
     # Explicit-end LI710 retains its existing raw-source coverage restriction.
     last<-if(!is.null(bounds$end) && site=='ECSM' && stage=='L3_LI710')min(expected_end,t[2]) else expected_end
     first<-if(is.null(bounds$start))t[1] else bounds$start
     stopifnot(w$start==first,w$end==last,
               identical(as.numeric(w$pending),as.numeric(seq(first,last,by=1800))))
   }
 }
}
# Actual publication replaces saved rows without explicit bounds and retains
# older rows outside the inferred EddyPro cap rather than deleting history.
old<-pipeline_history(cfg$dir_output[1],pipeline_patterns[['L1']])
old$data$PotRad<-999
old$data<-rbind(old$data,data.frame(TIMESTAMP=t[7],PotRad=777))
pipeline_write_table(old$data,old$units,
  file.path(cfg$dir_output[1],paste0('ECSM_',pipeline_patterns[['L1']],'_seed.csv')),'old')
run_pipeline('ECSM',stages='L1',reference_sites='ECSM',base_dir=root,config_file=config,reprocess=TRUE)
updated<-pipeline_history(cfg$dir_output[1],pipeline_patterns[['L1']])$data
stopifnot(all(updated$PotRad[updated$TIMESTAMP<=t[4]]==100),updated$PotRad[updated$TIMESTAMP==t[7]]==777)
# Different historical checkpoints still share one run-level figure period.
uneven<-cfg[1,,drop=FALSE];uneven$dir_output<-file.path(root,'uneven')
dir.create(uneven$dir_output);readr::write_csv(uneven,file.path(root,'uneven.csv'))
for(stage in names(pipeline_stages)) {
 n<-c(L1=4,L2=2,L3_EC=3,L3_LI710=1,L4=3)[[stage]]
 x<-data.frame(TIMESTAMP=t[seq_len(n)],PotRad=100);u<-units_for(x)
 if(stage=='L3_LI710') { item<-pipeline_empty_li710(t[seq_len(n)]);x<-item$data;u<-item$units }
 pipeline_write_table(x,u,file.path(uneven$dir_output,paste0('ECSM_',pipeline_patterns[[stage]],'_old.csv')),'old')
}
p<-run_pipeline('ECSM',reference_sites='ECSM',base_dir=root,config_file=file.path(root,'uneven.csv'),dry_run=TRUE)
stopifnot(p[['ECSM:L1']]=='skipped')
plans<-p[vapply(p,is.list,logical(1))]
stopifnot(length(unique(vapply(plans,pipeline_period_label,character(1))))==1,
          length(unique(vapply(plans,function(w)as.numeric(w$start),numeric(1))))>1)
# Advanced logger history does not let an L2-only automatic run exceed EddyPro.
one<-cfg[1,,drop=FALSE];one$dir_output<-file.path(root,'legacy');readr::write_csv(one,file.path(root,'legacy.csv'))
dir.create(one$dir_output)
x<-data.frame(TIMESTAMP=t,PotRad=100)
pipeline_write_table(x,units_for(x),file.path(one$dir_output,'ECSM_Level_1_QAQCed_logger_data_legacy.csv'),'legacy')
p<-run_pipeline('ECSM',stages='L2',reference_sites='ECSM',base_dir=root,config_file=file.path(root,'legacy.csv'),dry_run=TRUE)
stopifnot(p[[1]]$end==t[4])
# Explicit end retains caller control, even beyond EddyPro for logger-only work.
p<-run_pipeline('ECSM',stages='L1',reference_sites='ECSM',base_dir=root,config_file=config,end=t[8],dry_run=TRUE)
stopifnot(p[[1]]$end==t[8])
# Automatic mode fails before publishing if the EddyPro source is unavailable.
one$dir_eddypro<-file.path(root,'missing_ep');readr::write_csv(one,file.path(root,'missing.csv'))
err<-tryCatch(run_pipeline('ECSM',stages='L2',reference_sites='ECSM',base_dir=root,config_file=file.path(root,'missing.csv'),dry_run=TRUE),error=identity)
stopifnot(inherits(err,'error'))

# Run the REAL LI710 scientific adapter with a short source and then no source.
.pipeline_root<-real_root;pipeline_mds_window<-real_mds;pipeline_level4_figures<-real_fig
actual<-file.path(root,'actual');dir.create(file.path(actual,'figures'),recursive=TRUE)
clock<-seq(t[1],by=1800,length.out=96)
met<-data.frame(TIMESTAMP=clock,PotRad=rep(c(rep(0,24),rep(300,24)),2))
pipeline_write_table(met,units_for(met),file.path(actual,'ECSM_Level_1_QAQCed_logger_data_input.csv'),'met')
dirs<-data.frame(site='ECSM',dir_output=actual)
for(has_data in c(TRUE,FALSE)) {
 input<-pipeline_li710_read(if(has_data)cfg$dir_LI710[1] else character())
 options(ec.pipeline=list(dirs=dirs,base_dir=root,sites='ECSM',li710=input,
   window=list(start=clock[1],end=tail(clock,1),read_start=clock[1],pending=clock),figure_dir=file.path(actual,'figures')))
 old<-pipeline_files(actual,pipeline_patterns[['L3_LI710']]);if(length(old))unlink(old)
 suppressWarnings(sys.source(file.path(real_root,pipeline_stages[['L3_LI710']]),envir=new.env(parent=globalenv())))
 result<-pipeline_read_table(pipeline_files(actual,pipeline_patterns[['L3_LI710']])[1])$data
 stopifnot(nrow(result)==96,max(result$TIMESTAMP)==tail(clock,1))
 rows<-if(has_data)3:96 else 1:96
 stopifnot(all(is.na(result[rows,setdiff(names(result),'TIMESTAMP')])))
}
options(ec.pipeline=NULL)
unlink(root,recursive=TRUE)
cat('PASS: EddyPro caps, per-site ends, padding, unified figure folders, subset runs, explicit end, missing EddyPro, and real LI710 empty/tail QAQC.\n')

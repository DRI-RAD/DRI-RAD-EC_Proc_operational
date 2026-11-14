source('R/pipeline_io.R')
source('R/pipeline_li710.R')
root<-tempfile('li710-met-');dir.create(root)
t<-pipeline_time('2026-01-01 00:00:00')+0:3*1800
write_met<-function(file,idx,ta,pressure=90000,rh=.5,units=c('Pa','deg C','fraction')) {
 d<-data.frame(TIMESTAMP=t[idx],Atm_press=pressure,AirTC=ta,RH_710=rh)
 u<-data.frame(TIMESTAMP=c('TS','avg'),Atm_press=c(units[1],'avg'),AirTC=c(units[2],'avg'),RH_710=c(units[3],'avg'))
 pipeline_write_table(d,u,file.path(root,file),'fixture')
}
write_met('Site_LI710_old.backup',1:2,c(10,11))
write_met('Site_LI710.202601.backup',2:3,c(12,13))
write_met('Site_LI710.dat',3:4,c(14,15)+273.15,90,50,c('kPa','K','%'))
write_met('Another_LI710.dat',1,99)
row<-data.frame(dir_LI710=file.path(root,'Site_LI710.dat'),dir_LI710_old=file.path(root,'Site_LI710_old.backup'))
files<-pipeline_li710_files(row,root)
stopifnot(length(files)==3,basename(tail(files,1))=='Site_LI710.dat')
z<-pipeline_li710_read(files,met=TRUE)
stopifnot(nrow(z$data)==4,all(abs(z$data$AirTC-c(10,12,14,15))<1e-9),
 all(z$data$Atm_press==90),all(z$data$RH_710==50),z$units$Atm_press[1]=='kPa')
w<-pipeline_li710_read(files,list(read_start=t[2],end=t[3]),met=TRUE)
stopifnot(nrow(w$data)==2)
unlink(file.path(root,'Site_LI710.dat'))
row$dir_LI710<-NA_character_
stopifnot(nrow(pipeline_li710_read(pipeline_li710_files(row,root),met=TRUE)$data)==3)
# Missing historical channels are explicit NA; all-NA periods are readable.
d<-data.frame(TIMESTAMP=t[1],LE_710=10)
u<-data.frame(TIMESTAMP=c('TS','avg'),LE_710=c('W m-2','avg'))
f<-file.path(root,'NoMet.dat');pipeline_write_table(d,u,f,'fixture')
z<-suppressWarnings(pipeline_li710_read(f,met=TRUE))
stopifnot(all(is.na(z$data[,-1])),nrow(pipeline_li710_read(character(),met=TRUE)$data)==0)
invisible(parse('R/Level_1_met_QAQC_new.R'))
unlink(root,recursive=TRUE)
cat('PASS: old/backup/current sources, site isolation, active precedence, mixed units, old-only, time windows, missing met and empty input.\n')

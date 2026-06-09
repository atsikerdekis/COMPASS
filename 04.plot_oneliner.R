

############
### INIT ###
############
#file_out  <- paste0(path_plot,gsub(" ","",plot_title),"_",expname1,"-",expname2,"_",sDate,"-",eDate,".png")
figure_box <- FALSE
field_show_box <- FALSE
lonmin <- 0
lonmax <- 360
latmin <- -90
latmax <- +90
projection <- "+proj=robin"
gridlines=10
coastlineWorldFine_lwd=1
seqDate <- format(seq.Date(from=as.Date(sDate, format="%Y%m%d"), to=as.Date(eDate, format="%Y%m%d"), by="day"), "%Y%m%d")


###########
### NCO ### TODO: This can be cleaned and turned into a function or a separate Pre-processing step followed after download...
###########
# 0) ncks -v (only for pl to make them smaller)
# 1) ncrcat
# 2) ncra
# 3) ncwa -lon,lat
### Define files
expfile1_ncrcat <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_ncrcat.nc")
expfile1_ncra   <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_ncra.nc")
expfile1_ncwa   <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_ncwa.nc")
expfile2_ncrcat <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_ncrcat.nc")
expfile2_ncra   <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_ncra.nc")
expfile2_ncwa   <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_ncwa.nc")

if (1==2) {

temp_sfc_ncrcat <- paste0(path_temp,"temp_sfc_ncrcat.nc")
temp_sfc_ncra   <- paste0(path_temp,"temp_sfc_ncra.nc")
temp_sfc_ncwa   <- paste0(path_temp,"temp_sfc_ncwa.nc")
temp_pl_ncrcat  <- paste0(path_temp,"temp_pl_ncrcat.nc")
temp_pl_ncra    <- paste0(path_temp,"temp_pl_ncra.nc")
temp_pl_ncwa    <- paste0(path_temp,"temp_pl_ncwa.nc")
### EXP1
### Surface files
system(paste0(
  path_NCO,"ncrcat ", 
  paste0(path_data,expname1,"/CAMS_",expname1,"_forecast00to21by03_0.7x0.7_sfc_",seqDate,".nc ", collapse=" "), 
  temp_sfc_ncrcat
))
### Presure files
for (d in 1:length(seqDate)) {
  system(paste0(
    path_NCO,"ncks -v BU_AERO,BU_DU,BU_SS,BU_POM,BU_BC,BU_SO4,BU_NI,BU_AM ", 
    paste0(path_data,expname1,"/CAMS_",expname1,"_forecast00to21by03_0.7x0.7_pl_",seqDate[d],".nc ", collapse=" "), 
    path_temp, "temp_pl_",seqDate[d],".nc"
  ))
}
system(paste0(
  path_NCO,"ncrcat ", 
  paste0(path_temp,"temp_pl_",seqDate,".nc ", collapse=" "), 
  temp_pl_ncrcat
))
### Merge files
system(paste0("cp ",temp_sfc_ncrcat," ",expfile1_ncrcat))
system(paste0(path_NCO,"ncks -A ",temp_pl_ncrcat," ",expfile1_ncrcat))
### Spatial (ncwa) and temporal (ncra) mean
system(paste0(path_NCO,"ncwa -O -a longitude,latitude ",expfile1_ncrcat," ",expfile1_ncwa))
system(paste0(path_NCO,"ncra -O ", expfile1_ncrcat," ",expfile1_ncra))
### remove temp  files
system(paste0("rm ",path_temp,"temp*"))
### EXP2
### Surface files
system(paste0(
  path_NCO,"ncrcat ", 
  paste0(path_data,expname2,"/CAMS_",expname2,"_forecast00to21by03_0.7x0.7_sfc_",seqDate,".nc ", collapse=" "), 
  temp_sfc_ncrcat
))
### Presure files
for (d in 1:length(seqDate)) {
  system(paste0(
    path_NCO,"ncks -v BU_AERO,BU_DU,BU_SS,BU_POM,BU_BC,BU_SO4,BU_NI,BU_AM ", 
    paste0(path_data,expname2,"/CAMS_",expname2,"_forecast00to21by03_0.7x0.7_pl_",seqDate[d],".nc ", collapse=" "), 
    path_temp, "temp_pl_",seqDate[d],".nc"
  ))
}
system(paste0(
  path_NCO,"ncrcat ", 
  paste0(path_temp,"temp_pl_",seqDate,".nc ", collapse=" "), 
  temp_pl_ncrcat
))
### Merge files
system(paste0("cp ",temp_sfc_ncrcat," ",expfile2_ncrcat))
system(paste0(path_NCO,"ncks -A ",temp_pl_ncrcat," ",expfile2_ncrcat))
### Spatial (ncwa) and temporal (ncra) mean
system(paste0(path_NCO,"ncwa -O -a longitude,latitude ",expfile2_ncrcat," ",expfile2_ncwa))
system(paste0(path_NCO,"ncra -O ", expfile2_ncrcat," ",expfile2_ncra))
### remove temp  files
system(paste0("rm ",path_temp,"temp*"))

}



for (v in 1:length(variable_title)) {

nvar <- 1
file_out  <- paste0(path_plot,gsub(" ","",plot_title),"_",gsub(" ","",variable_name[v]),"_",expname1,"-",expname2,"_",sDate,"-",eDate,".png")

############
### MAPS ###
############
dpi <- 300
png(file_out, width=(0.2+3*3.9+2*0.5+4)*dpi, height=(0.23+0.15+2*nvar)*dpi)
myl <- layout(mat=matrix(c(1,1,1,1,1,1,1,2:(((2+nvar)*7-6))),2+nvar,7,byrow=T), widths=c(0.2,3.9,3.9,0.5,3.9,0.5,4), heights=c(0.23,0.15,rep(2,nvar)))
#myl <- layout(mat=matrix(c(1,1,1,1,1,1,1,2:(((3+nvar)*6))),2+nvar,7,byrow=T), widths=c(0.2,3.9,3.9,0.4,3.9,0.4,4), heights=c(0.23,0.15,rep(2,nvar)))
layout.show(myl)

par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",plot_title,"   |   Period: ",sDate,"-",eDate), col="grey50", cex=7, family="Century Gothic", srt=0)
abline(h=0, col="grey50", lwd=3)
abline(h=1, col="grey50", lwd=3)
par(mai=c(0,0,0,0));plot.new()
par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0(expname1), col="grey20", cex=5, family="Century Gothic", srt=0)
par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0(expname2), col="grey20", cex=5, family="Century Gothic", srt=0)
par(mai=c(0,0,0,0));plot.new()
par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0(expname1," - ",expname2), col="grey20", cex=5, family="Century Gothic", srt=0)
par(mai=c(0,0,0,0));plot.new()
par(mai=c(0,0,0,0));plot.new()

#for (v in 1:nvar) {
  message(paste0("Plotting: (",v,") ",variable_title[v]))
  ### EXP1 ncra
  file.nc <- nc_open(paste0(expfile1_ncra))
  field_lon  <- ncvar_get(file.nc, "longitude")
  field_lat  <- ncvar_get(file.nc, "latitude")
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { field_var1 <- ncvar_get(file.nc, variable_name[v]) * variable_mul }
  if (variable_name[v] == "ae550to865") { 
    field_temp1 <- ncvar_get(file.nc, "aod550") 
    field_temp2 <- ncvar_get(file.nc, "aod865")
    field_var1  <- get_AngstromExponent(aodW1=field_temp1, aodW2=field_temp2, wave1=550, wave2=865)
    rm(field_temp1, field_temp2)
  }
  if (variable_name[v] == "mec550") {
    field_temp1 <- ncvar_get(file.nc, "aod550")
    field_temp2 <- ncvar_get(file.nc, "BU_AERO") 
    field_var1  <- field_temp1 / field_temp2 # MEC = AOD / MASS
    rm(field_temp1, field_temp2)
  }
  nc_close(file.nc)
  ### EXP1 ncwa
  file.nc   <- nc_open(paste0(expfile1_ncwa))
  tmean_tim <- ncvar_get(file.nc, "time")
  tmean_tim <- as.POSIXct(tmean_tim*3600, origin="1900-01-01 00:00:00", tz="UTC")
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { tmean_var1 <- ncvar_get(file.nc, variable_name[v]) * variable_mul }
  if (variable_name[v] == "ae550to865") { 
    tmean_temp1 <- ncvar_get(file.nc, "aod550") 
    tmean_temp2 <- ncvar_get(file.nc, "aod865")
    tmean_var1  <- get_AngstromExponent(aodW1=tmean_temp1, aodW2=tmean_temp2, wave1=550, wave2=865)
    rm(tmean_temp1, tmean_temp2)
  }
  if (variable_name[v] == "mec550") {
    tmean_temp1 <- ncvar_get(file.nc, "aod550")
    tmean_temp2 <- ncvar_get(file.nc, "BU_AERO") 
    tmean_var1  <- tmean_temp1 / tmean_temp2 # MEC = AOD / MASS
    rm(tmean_temp1, tmean_temp2)
  }
  nc_close(file.nc)
  ### EXP2 ncra
  file.nc <- nc_open(paste0(expfile2_ncra))
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { field_var2 <- ncvar_get(file.nc, variable_name[v]) * variable_mul }
  if (variable_name[v] == "ae550to865") { 
    field_temp1 <- ncvar_get(file.nc, "aod550") 
    field_temp2 <- ncvar_get(file.nc, "aod865") 
    field_var2  <- get_AngstromExponent(aodW1=field_temp1, aodW2=field_temp2, wave1=550, wave2=865)
    rm(field_temp1, field_temp2)
  }
  if (variable_name[v] == "mec550") {
    field_temp1 <- ncvar_get(file.nc, "aod550")
    field_temp2 <- ncvar_get(file.nc, "BU_AERO") 
    field_var2  <- field_temp1 / field_temp2 # MEC = AOD / MASS
    rm(field_temp1, field_temp2)
  }
  nc_close(file.nc)
  ### EXP2 ncwa
  file.nc <- nc_open(paste0(expfile2_ncwa))
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { tmean_var2 <- ncvar_get(file.nc, variable_name[v]) * variable_mul }
  if (variable_name[v] == "ae550to865") { 
    tmean_temp1 <- ncvar_get(file.nc, "aod550") 
    tmean_temp2 <- ncvar_get(file.nc, "aod865")
    tmean_var2  <- get_AngstromExponent(aodW1=tmean_temp1, aodW2=tmean_temp2, wave1=550, wave2=865)
    rm(tmean_temp1, tmean_temp2)
  }
  if (variable_name[v] == "mec550") {
    tmean_temp1 <- ncvar_get(file.nc, "aod550")
    tmean_temp2 <- ncvar_get(file.nc, "BU_AERO") 
    tmean_var2  <- tmean_temp1 / tmean_temp2 # MEC = AOD / MASS
    rm(tmean_temp1, tmean_temp2)
  }
  
  par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, variable_title[v], col="grey20", cex=5, family="Century Gothic", srt=90)
  ### MAP EXP1
  MapNC(filename_topo="",
        figure_box=figure_box,
        field_show_box=field_show_box,
        coastlineWorldFine_lwd=coastlineWorldFine_lwd,
        gridlines=gridlines,
        projection=projection, lonmax=lonmax, lonmin=lonmin, latmax=latmax, latmin=latmin,
        field_value=field_var1, 
        field_lon=field_lon, 
        field_lat=field_lat,
        field_pallete_name=field_pallete_name[v],
        field_breaks=field_breaks[[v]],
        field_units=variable_units[v],
        field_pallete_starting_alpha=100,
        field_show_legend=F)
  ### MAP EXP2
  MapNC(filename_topo="",
        figure_box=figure_box,
        field_show_box=field_show_box,
        coastlineWorldFine_lwd=coastlineWorldFine_lwd,
        gridlines=gridlines,
        projection=projection, lonmax=lonmax, lonmin=lonmin, latmax=latmax, latmin=latmin,
        field_value=field_var2, 
        field_lon=field_lon, 
        field_lat=field_lat,
        field_pallete_name=field_pallete_name[v],
        field_breaks=field_breaks[[v]],
        field_units=variable_units[v],
        field_pallete_starting_alpha=100,
        field_show_legend=T)
  ### MAP EXP1 - EXP2
  MapNC(filename_topo="",
        figure_box=figure_box,
        field_show_box=field_show_box,
        coastlineWorldFine_lwd=coastlineWorldFine_lwd,
        gridlines=gridlines,
        projection=projection, lonmax=lonmax, lonmin=lonmin, latmax=latmax, latmin=latmin,
        field_value=field_var1-field_var2, 
        field_lon=field_lon, 
        field_lat=field_lat,
        field_pallete_name="MNMB",
        field_breaks=field_breaks_diff[[v]],
        field_units=variable_units[v],
        field_pallete_starting_alpha=100,
        field_show_legend=T)
  
  ### TIMESERIES
  par(mai=c(1.2,1.25,0,0.4), family="Century Gothic")
  x     <- 1:length(tmean_tim) # Minus 0.5 to align barplot with line/point plot
  #ymax  <- max(c(tmean_var1,tmean_var2),na.rm=T); ymax <- ymax + ymax*0.1
  #ymin  <- min(c(tmean_var1,tmean_var2),na.rm=T); ymin <- ymin - ymin*0.1
  #ymax  <- max(c(field_breaks[[v]]),na.rm=T); ymax <- ymax + ymax*0.1
  #ymin  <- min(c(field_breaks[[v]]),na.rm=T); ymin <- ymin - ymin*0.1
  #yseq  <- field_breaks[[v]]
  x     <- 1:length(tmean_tim) # Minus 0.5 to align barplot with line/point plot
  ymax  <- max(c(tmean_var1,tmean_var2),na.rm=T); ymax <- ymax + ymax*0.1
  ymin  <- min(c(tmean_var1,tmean_var2),na.rm=T); ymin <- ymin - ymin*0.1
  yseq  <- pretty(diff(seq(ymin, ymax, length.out=11)))[1]
  if (yseq==0) {yseq <- pretty(diff(seq(ymin, ymax, length.out=11)))[2]}
  
  
  plot(x=x,type='n',axes=FALSE,ann=FALSE,ylim=c(ymin,ymax), yaxs="i")
  mtext("Time (3 hourly)", side=1, line=6, cex=2.5)
  IDx_labels <- c(which(format(tmean_tim,"%d %H:%M")=="01 00:00" | format(tmean_tim,"%d %H:%M")=="10 00:00" | format(tmean_tim,"%d %H:%M")=="20 00:00"), length(tmean_tim) )
  axis(side=1, at=x[IDx_labels], labels=paste0(format(tmean_tim[IDx_labels],"%Y-%m-%d"),"\n",format(tmean_tim[IDx_labels],"%H:%M:%S")), cex.axis=2.5, line=2.5, lty=0)
  axis(side=1, at=x[IDx_labels], labels=FALSE, cex.axis=3, lty=1, tck=0.01)
  axis(side=1, at=x[IDx_labels], labels=FALSE, cex.axis=3, lty=1, tck=-0.01)
  #axis(side=2, at=yseq,las=1,cex.axis=2.5)
  axis(side=2, at=format(seq(ymin,ymax,yseq), scientific=TRUE, digits=2),las=1,cex.axis=3)
  box(lwd=2)
  #abline(h=yseq, lwd=0.2, col="grey")
  abline(v=x[IDx_labels], lwd=0.2, col="grey")
  abline(v=x[IDx_labels], lwd=0.2, col="grey")
  lines(x=x,y=tmean_var1,lwd=8);points(x=x,y=tmean_var1,pch=19,cex=2.8)
  lines(x=x,y=tmean_var2,lwd=8);points(x=x,y=tmean_var2,pch=19,cex=2.8)
  lines(x=x,y=tmean_var1,lwd=5,col="blue");points(x=x,y=tmean_var1,pch=19,cex=2.2,col="blue")
  lines(x=x,y=tmean_var2,lwd=5,col="red");points(x=x,y=tmean_var2,pch=19,cex=2.2,col="red")
  legend("top", legend=c(expname1,expname2), lwd=c(7,7), pch=c(19,19), col=c("blue","red"), cex=c(3,3))
  points(x=-0.01630435*length(x),y=mean(tmean_var1,na.rm=T),col="black" , cex=6.8, pch=18)
  points(x=-0.01630435*length(x),y=mean(tmean_var1,na.rm=T),col="blue"   , cex=6.0, pch=18)
  points(x=-0.01630435*length(x),y=mean(tmean_var2,na.rm=T),col="black" , cex=6.8, pch=18)
  points(x=-0.01630435*length(x),y=mean(tmean_var2,na.rm=T),col="red"  , cex=6.0, pch=18)
}
dev.off()

compress(file_in=file_out, file_out=file_out)





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


if (1==1) {

###########
### CDO ### TODO: This can be cleaned and turned into a function or a separate Pre-processing step followed after download...
###########
# 0) ncks -v (only for pl to make them smaller)
# 1) mergetime
# 2) timmean
# 3) fldmean
# 4) fldmean_daymean
### Define files
expfile1_mergetime <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_mergetime.nc")
expfile1_timmean   <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_timmean.nc")
expfile1_fldmean   <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_fldmean.nc")
expfile1_fldmean_daymean   <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_fldmean_daymean.nc")
expfile1_fldmean_dhourmean <- paste0(path_temp,"CAMS_",expname1,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_fldmean_dhourmean.nc")
expfile2_mergetime <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_mergetime.nc")
expfile2_timmean   <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_timmean.nc")
expfile2_fldmean   <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_fldmean.nc")
expfile2_fldmean_daymean   <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_fldmean_daymean.nc")
expfile2_fldmean_dhourmean <- paste0(path_temp,"CAMS_",expname2,"_forecast00to21by03_0.7x0.7_",sDate,"-",eDate,"_fldmean_dhourmean.nc")
temp_sfc_mergetime <- paste0(path_temp,"temp_sfc_mergetime.nc")
temp_sfc_timmean   <- paste0(path_temp,"temp_sfc_timmean.nc")
temp_sfc_fldmean   <- paste0(path_temp,"temp_sfc_fldmean.nc")
temp_pl_mergetime  <- paste0(path_temp,"temp_pl_mergetime.nc")
temp_pl_timmean    <- paste0(path_temp,"temp_pl_timmean.nc")
temp_pl_fldmean    <- paste0(path_temp,"temp_pl_fldmean.nc")

###
### EXP1
###
### Surface files
system(paste0(
  path_CDO,"cdo -b F32 -O mergetime ", 
  paste0(path_data,expname1,"/CAMS_",expname1,"_forecast00to21by03_0.7x0.7_sfc_",seqDate,".nc ", collapse=" "), 
  temp_sfc_mergetime
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
  path_CDO,"cdo -b F32 -O mergetime ", 
  paste0(path_temp,"temp_pl_",seqDate,".nc ", collapse=" "), 
  temp_pl_mergetime
))
### Merge files
system(paste0("cp ",temp_sfc_mergetime," ",expfile1_mergetime))
system(paste0(path_NCO,"ncks -A ",temp_pl_mergetime," ",expfile1_mergetime))
### Spatial (fldmean) and temporal (timmean) mean
system(paste0(path_CDO,"cdo -O fldmean ",expfile1_mergetime," ",expfile1_fldmean))
system(paste0(path_CDO,"cdo -O timmean ", expfile1_mergetime," ",expfile1_timmean))
### Daily mean of global mean (from 3 hourly to daily)
system(paste0(path_CDO,"cdo -O daymean ",expfile1_fldmean," ",expfile1_fldmean_daymean))
### Daily cycle of global mean
system(paste0(path_CDO,"cdo -O dhourmean ",expfile1_fldmean," ",expfile1_fldmean_dhourmean))
### remove temp  files
system(paste0("rm ",path_temp,"temp*"))

###
### EXP2
###
### Surface files
system(paste0(
  path_CDO,"cdo -b F32 -O mergetime ", 
  paste0(path_data,expname2,"/CAMS_",expname2,"_forecast00to21by03_0.7x0.7_sfc_",seqDate,".nc ", collapse=" "), 
  temp_sfc_mergetime
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
  path_CDO,"cdo -b F32 -O mergetime ", 
  paste0(path_temp,"temp_pl_",seqDate,".nc ", collapse=" "), 
  temp_pl_mergetime
))
### Merge files
system(paste0("cp ",temp_sfc_mergetime," ",expfile2_mergetime))
system(paste0(path_NCO,"ncks -A ",temp_pl_mergetime," ",expfile2_mergetime))
### Spatial (fldmean) and temporal (timmean) mean
system(paste0(path_CDO,"cdo -O fldmean ",expfile2_mergetime," ",expfile2_fldmean))
system(paste0(path_CDO,"cdo -O timmean ", expfile2_mergetime," ",expfile2_timmean))
### Daily mean of global mean (from 3 hourly to daily)
system(paste0(path_CDO,"cdo -O daymean ",expfile2_fldmean," ",expfile2_fldmean_daymean))
### Hourly cycle of global mean
system(paste0(path_CDO,"cdo -O dhourmean ",expfile2_fldmean," ",expfile2_fldmean_dhourmean))
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
png(file_out, width=(0.2+3*3.9+0.8+0.8)*dpi, height=(0.23+0.15+2+2.5)*dpi)
#myl <- layout(mat=matrix(c(1,1,1,1,1,1,1,2:(((2+nvar)*7-6))),2+nvar,7,byrow=T), widths=c(0.2,3.9,3.9,0.5,3.9,0.5,4), heights=c(0.23,0.15,rep(2,nvar)))
#myl <- layout(mat=matrix(c(1,1,1,1,1,1,1,2:(((3+nvar)*6))),2+nvar,7,byrow=T), widths=c(0.2,3.9,3.9,0.4,3.9,0.4,4), heights=c(0.23,0.15,rep(2,nvar)))
myl <- layout(mat=matrix(c(1,1,1,1,1,1,2:13,14,14,14,14,15,15),4,6,byrow=T), widths=c(0.2,3.9,3.9,0.8,3.9,0.8), heights=c(0.23,0.15,2,2.5))
layout.show(myl)

par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0("Experiments: ",expname1," VS ",expname2,"   |   Type: ",plot_title,"   |   Period: ",sDate,"-",eDate), col="grey50", cex=6, family="Century Gothic", srt=0)
abline(h=0, col="grey50", lwd=3)
abline(h=1, col="grey50", lwd=3)
par(mai=c(0,0,0,0));plot.new()
par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0(expname1," (",exptype1,")"), col="grey20", cex=4.5, family="Century Gothic", srt=0)
par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0(expname2," (",exptype2,")"), col="grey20", cex=4.5, family="Century Gothic", srt=0)
par(mai=c(0,0,0,0));plot.new()
par(mai=c(0,0,0,0));plot.new();text(0.5,0.5, paste0(expname2," - ",expname1), col="grey20", cex=4.5, family="Century Gothic", srt=0)
par(mai=c(0,0,0,0));plot.new()
#par(mai=c(0,0,0,0));plot.new()

#for (v in 1:nvar) {
  message(paste0("Plotting: (",v,") ",variable_title[v]))
  ### EXP1 timmean
  file.nc <- nc_open(paste0(expfile1_timmean))
  field_lon  <- ncvar_get(file.nc, "longitude")
  field_lat  <- ncvar_get(file.nc, "latitude")
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { field_var1 <- ncvar_get(file.nc, variable_name[v]) * variable_mul[v] }
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
  ### EXP1 fldmean_daymean
  file.nc   <- nc_open(paste0(expfile1_fldmean_daymean))
  tmean_tim <- ncvar_get(file.nc, "time")
  tmean_tim <- as.POSIXct(tmean_tim*3600, origin="1900-01-01 00:00:00", tz="UTC")
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { tmean_var1 <- ncvar_get(file.nc, variable_name[v]) * variable_mul[v] }
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
  ### EXP1 fldmean_dhourmean
  file.nc   <- nc_open(paste0(expfile1_fldmean_dhourmean))
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { dhourmean_var1 <- ncvar_get(file.nc, variable_name[v]) * variable_mul[v] }
  if (variable_name[v] == "ae550to865") {
    tmean_temp1 <- ncvar_get(file.nc, "aod550")
    tmean_temp2 <- ncvar_get(file.nc, "aod865")
    dhourmean_var1  <- get_AngstromExponent(aodW1=tmean_temp1, aodW2=tmean_temp2, wave1=550, wave2=865)
    rm(tmean_temp1, tmean_temp2)
  }
  if (variable_name[v] == "mec550") {
    tmean_temp1 <- ncvar_get(file.nc, "aod550")
    tmean_temp2 <- ncvar_get(file.nc, "BU_AERO")
    dhourmean_var1 <- tmean_temp1 / tmean_temp2 # MEC = AOD / MASS
    rm(tmean_temp1, tmean_temp2)
  }
  nc_close(file.nc)

  ### EXP2 timmean
  file.nc <- nc_open(paste0(expfile2_timmean))
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { field_var2 <- ncvar_get(file.nc, variable_name[v]) * variable_mul[v] }
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
  ### EXP2 fldmean_daymean
  file.nc <- nc_open(paste0(expfile2_fldmean_daymean))
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { tmean_var2 <- ncvar_get(file.nc, variable_name[v]) * variable_mul[v] }
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
  ### EXP1 fldmean_dhourmean
  file.nc   <- nc_open(paste0(expfile2_fldmean_dhourmean))
  if (variable_name[v] != "ae550to865" && variable_name[v] != "mec550") { dhourmean_var2 <- ncvar_get(file.nc, variable_name[v]) * variable_mul[v] }
  if (variable_name[v] == "ae550to865") {
    tmean_temp1 <- ncvar_get(file.nc, "aod550")
    tmean_temp2 <- ncvar_get(file.nc, "aod865")
    dhourmean_var2  <- get_AngstromExponent(aodW1=tmean_temp1, aodW2=tmean_temp2, wave1=550, wave2=865)
    rm(tmean_temp1, tmean_temp2)
  }
  if (variable_name[v] == "mec550") { 
    tmean_temp1 <- ncvar_get(file.nc, "aod550")
    tmean_temp2 <- ncvar_get(file.nc, "BU_AERO")
    dhourmean_var2 <- tmean_temp1 / tmean_temp2 # MEC = AOD / MASS
    rm(tmean_temp1, tmean_temp2)
  }
  nc_close(file.nc)
  
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
        field_show_legend=T, 
        field_legend_mai_right=1.8)
  ### MAP EXP1 - EXP2
  MapNC(filename_topo="",
        figure_box=figure_box,
        field_show_box=field_show_box,
        coastlineWorldFine_lwd=coastlineWorldFine_lwd,
        gridlines=gridlines,
        projection=projection, lonmax=lonmax, lonmin=lonmin, latmax=latmax, latmin=latmin,
        field_value=field_var2-field_var1, 
        field_lon=field_lon, 
        field_lat=field_lat,
        field_pallete_name="MNMB",
        field_breaks=field_breaks_diff[[v]],
        field_units=variable_units[v],
        field_pallete_starting_alpha=100,
        field_show_legend=T,
        field_legend_mai_right=1.8)
  
  ### TIMESERIES
  par(mai=c(2.0,2.0,0,0.4), family="Century Gothic")
  #x     <- 1:length(tmean_tim) # Minus 0.5 to align barplot with line/point plot
  #ymax  <- max(c(tmean_var1,tmean_var2),na.rm=T); ymax <- ymax + ymax*0.1
  #ymin  <- min(c(tmean_var1,tmean_var2),na.rm=T); ymin <- ymin - ymin*0.1
  #ymax  <- max(c(field_breaks[[v]]),na.rm=T); ymax <- ymax + ymax*0.1
  #ymin  <- min(c(field_breaks[[v]]),na.rm=T); ymin <- ymin - ymin*0.1
  #yseq  <- field_breaks[[v]]
  x     <- 1:length(tmean_tim) # Minus 0.5 to align barplot with line/point plot
  #ymax  <- quantile(c(tmean_var1,tmean_var2),0.99,na.rm=T); ymax <- ymax + ymax*0.1
  #ymin  <- quantile(c(tmean_var1,tmean_var2),0.01,na.rm=T); ymin <- ymin - ymin*0.1
  #yseq  <- pretty(diff(seq(ymin, ymax, length.out=11)))[1]
  #if (yseq==0) {yseq <- pretty(diff(seq(ymin, ymax, length.out=11)))[2]}
  #axis_ticks <- function(x, n = 5) {
  #  rng <- range(x, na.rm = TRUE, finite = TRUE)
  #  ticks <- pretty(rng, n = n)
  #  list( breaks = ticks, labels = scales::label_number()(ticks) )
  #}
  
axis_ticks <- function(x, n = 5, small_thresh = 0.01, large_thresh = 100) {
  rng <- range(x, na.rm = TRUE, finite = TRUE)
  ticks <- pretty(rng, n = n)
  
  max_abs <- max(abs(ticks))
  min_abs <- if (all(ticks == 0)) 0 else min(abs(ticks[ticks != 0]))
  
  if (max_abs >= large_thresh || (min_abs > 0 && min_abs < small_thresh)) {
    # Scientific notation
    labels <- sprintf("%.2e", ticks)
  } else {
    # Adaptive decimals: at least 2, but more if needed to distinguish
    ndigits <- max(2, -floor(log10(min(diff(ticks)))))
    labels <- formatC(ticks, format = "f", digits = ndigits)
  }
  
  list(breaks = ticks, labels = labels)
}



  yseq <- axis_ticks(c(tmean_var1,tmean_var2), n=10)
  ymin <- min(yseq$breaks)
  ymax <- max(yseq$breaks)
  pad  <- diff(range(yseq$breaks))*0.05

  plot(x=x,type='n',axes=FALSE,ann=FALSE,ylim=c(ymin-pad,ymax+pad), yaxs="i")
  mtext("Time (Day)", side=1, line=12, cex=3.5)
#  IDx_labels <- c(which(format(tmean_tim,"%d %H:%M")=="01 00:00" |
#                        format(tmean_tim,"%d %H:%M")=="05 00:00" | 
#                        format(tmean_tim,"%d %H:%M")=="10 00:00" | 
#                        format(tmean_tim,"%d %H:%M")=="15 00:00" |
#                        format(tmean_tim,"%d %H:%M")=="20 00:00" |
#                        format(tmean_tim,"%d %H:%M")=="25 00:00"), length(tmean_tim) )
  IDx_labels <- c(which(format(tmean_tim,"%d")=="01" |
                        format(tmean_tim,"%d")=="05" |
                        format(tmean_tim,"%d")=="10" |
                        format(tmean_tim,"%d")=="15" |
                        format(tmean_tim,"%d")=="20" |
                        format(tmean_tim,"%d")=="25"), length(tmean_tim) )
#  axis(side=1, at=x[IDx_labels], labels=paste0(format(tmean_tim[IDx_labels],"%Y-%m-%d"),"\n",format(tmean_tim[IDx_labels],"%H:%M:%S")), cex.axis=4.0, line=5.5, lty=0)
  axis(side=1, at=x[IDx_labels], labels=paste0(format(tmean_tim[IDx_labels],"%Y-%m-%d")), cex.axis=4.0, line=4, lty=0)
  axis(side=1, at=x[IDx_labels], labels=FALSE, cex.axis=4, lty=1, tck=0.01)
  axis(side=1, at=x[IDx_labels], labels=FALSE, cex.axis=4, lty=1, tck=-0.01)
  axis(side=2, at=yseq$breaks, labels=yseq$labels, las=1,cex.axis=3)
  box(lwd=2)
  abline(h=yseq$breaks, lwd=1.0, col="grey")
  abline(v=x[IDx_labels], lwd=1.0, col="grey")
  lines(x=x,y=tmean_var1,lwd=8);points(x=x,y=tmean_var1,pch=19,cex=2.8)
  lines(x=x,y=tmean_var2,lwd=8);points(x=x,y=tmean_var2,pch=19,cex=2.8)
  lines(x=x,y=tmean_var1,lwd=5,col="red");points(x=x,y=tmean_var1,pch=19,cex=2.2,col="red")
  lines(x=x,y=tmean_var2,lwd=5,col="blue");points(x=x,y=tmean_var2,pch=19,cex=2.2,col="blue")
  legend("top", legend=c(paste0(expname1,"(",exptype1,")"),paste0(expname2,"(",exptype2,")")), lwd=c(7,7), pch=c(19,19), col=c("red","blue"), cex=c(4,4))
  points(x=-0.01630435*length(x),y=mean(tmean_var1,na.rm=T),col="black" , cex=6.8, pch=18)
  points(x=-0.01630435*length(x),y=mean(tmean_var1,na.rm=T),col="red"   , cex=6.0, pch=18)
  points(x=-0.01630435*length(x),y=mean(tmean_var2,na.rm=T),col="black" , cex=6.8, pch=18)
  points(x=-0.01630435*length(x),y=mean(tmean_var2,na.rm=T),col="blue"  , cex=6.0, pch=18)

  ### DAILY CYCLE
  yseq <- axis_ticks(c(dhourmean_var1,dhourmean_var2), n = 10)
  ymin <- min(yseq$breaks)
  ymax <- max(yseq$breaks)
  pad  <- diff(range(yseq$breaks))*0.05
  print(yseq)
  par(mai=c(2.0,2.0,0,0.4), family="Century Gothic")
  plot(x=1:8, type='n',axes=FALSE,ann=FALSE,ylim=c(ymin-pad,ymax+pad), yaxs="i")
  mtext("Time (3 hourly UTC)", side=1, line=12, cex=3.5)
  axis(side=1, at=1:8, labels=c("00","03","06","09","12","15","18","21"), cex.axis=4.0, line=4, lty=0)
  axis(side=1, at=1:8, labels=FALSE, cex.axis=4, lty=1, tck=0.01)
  axis(side=1, at=1:8, labels=FALSE, cex.axis=4, lty=1, tck=-0.01)
  axis(side=2, at=yseq$breaks, labels=yseq$labels,las=1,cex.axis=3)
  box(lwd=2)
  abline(h=yseq$breaks, lwd=1.0, col="grey")
  abline(v=1:8,, lwd=1.0, col="grey")
  # normalized daily cycle lines and points (inflates the cycle to the range of the yy axis)
  ynormalized1 = Normalization(x=dhourmean_var1, ymin=ymin, ymax=ymax)
  ynormalized2 = Normalization(x=dhourmean_var2, ymin=ymin, ymax=ymax)
  lines(x=1:8,,y=ynormalized1,lwd=7, lty=2, col=rgb(1,0,0,0.5));points(x=1:8,,y=ynormalized1,pch=19,cex=3.2,col=rgb(1,0,0,0.5))
  lines(x=1:8,,y=ynormalized2,lwd=7, lty=2, col=rgb(0,0,1,0.5));points(x=1:8,,y=ynormalized2,pch=19,cex=3.2,col=rgb(0,0,1,0.5))
  # daily cycle lines and points
  lines(x=1:8,,y=dhourmean_var1,lwd=10);points(x=1:8,,y=dhourmean_var1,pch=19,cex=3.8)
  lines(x=1:8,,y=dhourmean_var2,lwd=10);points(x=1:8,,y=dhourmean_var2,pch=19,cex=3.8)
  lines(x=1:8,,y=dhourmean_var1,lwd=7,col="red");points(x=1:8,,y=dhourmean_var1,pch=19,cex=3.2,col="red")
  lines(x=1:8,,y=dhourmean_var2,lwd=7,col="blue");points(x=1:8,,y=dhourmean_var2,pch=19,cex=3.2,col="blue")
}
dev.off()

compress(file_in=file_out, file_out=file_out)



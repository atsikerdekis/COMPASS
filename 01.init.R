############
### USER ###
############
username <- system("whoami", intern = TRUE)

############
### PATH ###
############
path_function <- paste0("/perm/",username,"/cams2_35/COMPASS/code/function/")
path_plot     <- paste0("/perm/",username,"/cams2_35/COMPASS/plot/",expname1,"_",expname2,"_",sDate,"-",eDate,"/")
path_log      <- paste0("/perm/",username,"/cams2_35/COMPASS/log/")
path_data     <- paste0("/scratch/",username,"/cams2_35/data/")
path_temp     <- paste0("/scratch/",username,"/cams2_35/temp/")
path_NCO      <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")
path_CDO      <- paste0("/usr/local/apps/cdo/2.5.1/bin/")
path_PYTHON   <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")
path_R        <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")

###################
### CREATE PATH ###
###################
suppressWarnings(dir.create(paste0(path_plot), recursive=T, showWarnings=T))
suppressWarnings(dir.create(paste0(path_log) , recursive=T, showWarnings=T))
suppressWarnings(dir.create(paste0(path_temp), recursive=T, showWarnings=T))
suppressWarnings(dir.create(paste0(path_data,"/",expname1), recursive=T, showWarnings=T))
suppressWarnings(dir.create(paste0(path_data,"/",expname2), recursive=T, showWarnings=T))

#################
### LIBRARIES ###
#################
suppressWarnings( library("ncdf4"   , quietly=T) )
suppressWarnings( library("showtext", quietly=T) )
font_add(family = "Century Gothic", regular = "/home/nktt/fonts/CenturyGothic/centurygothic.ttf")
showtext_auto()


################
### FUNCTION ###
################
source(paste0(path_function,"MapNC.R"))
source(paste0(path_function,"get_AngstromExponent.R"))
source(paste0(path_function,"compress.R"))
source(paste0(path_function,"SubmitJob.R"))
source(paste0(path_function,"Normalization.R"))
Sys.setenv(PROJ_LIB = "/etc/ecmwf/nfs/dh1_perm_b/nktt/miniforge3/envs/COMPASS/share/proj") # So sf can stop complaining...

############
### INIT ###
############
#nvar <- length(variable_name) 

################
### METADATA ### TODO: Could add a big list of supported dataset instead of a config file
################       To make this more flexible... it will get quite complicated... 
### TODO: Could add a big list of supported dataset instead of a config file
###       To make this more flexible... it will get quite complicated... 
###       Take MEC e.g. it will require sfc and pl stream
###       Lots of changes will have to be done (switches) in the download and pre-processing step...
###       Maybe in the future could be more generalizaed. 
###       For now use fixed config files for predefined figures variables
# Supported parameters from surface (sfc) stream
#param_sfc <- list(
#  VARNAME=c(
#    "AOD550",
#    "AAOD550",
#    "SSA550",
#    "AE550865",
#    "MEC550"
#  ), 
#  ID=c(
#    "207.210",
#    "104.215"
#  ),
#  VARTITLE=c(
#    "AOD 550nm",
#    "AE 550-865nm"
#  )
#)

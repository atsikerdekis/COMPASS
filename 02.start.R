###
### Created : Thanos Tsikerdekis (KNMI) | May 2025 |
### Contact : thanos.tsikerdekis@knmi.nl
### Purpose : Start a comparison
### Example : Rscript 01.start.R
### Envirom : ...
###

###################
### START TIMER ###
###################
stime <- Sys.time()
message(paste0("---> Starting... "))

##################
### INITIALIZE ### Pretty dum for now add some checks...
##################
args      <- commandArgs(trailingOnly=TRUE)
expname1  <- args[1]
exptype1  <- args[2]
expname2  <- args[3]
exptype2  <- args[4]
sDate     <- args[5]
eDate     <- args[6]
path_code <- paste0(dirname(normalizePath(sub("--file=", "", args[grep("--file=", args)]))), "/")
source(paste0(path_code,"config.R"))
source(paste0(path_code,"01.init.R"))

########################
### CHECK & DOWNLOAD ###
########################
message(paste0("---> Downloading ",expname1," (",exptype1,")..."))
system(paste0(path_PYTHON,"python ",path_code,"03.download_sfc.py ",expname1," ",sDate," ",eDate," ",path_data))
if (exptype1 == "HAM") { system(paste0(path_PYTHON,"python ",path_code,"03.download_pl_M7.py ",expname1," ",sDate," ",eDate," ",path_data)) }
if (exptype1 == "AER") { system(paste0(path_PYTHON,"python ",path_code,"03.download_pl.py ",expname1," ",sDate," ",eDate," ",path_data)) }
message(paste0("---> Downloading ",expname2," (",exptype2,")..."))
system(paste0(path_PYTHON,"python ",path_code,"03.download_sfc.py ",expname2," ",sDate," ",eDate," ",path_data))
if (exptype2 == "HAM") { system(paste0(path_PYTHON,"python ",path_code,"03.download_pl_M7.py ",expname2," ",sDate," ",eDate," ",path_data)) }
if (exptype2 == "AER") { system(paste0(path_PYTHON,"python ",path_code,"03.download_pl.py ",expname2," ",sDate," ",eDate," ",path_data)) }

############
### PLOT ###
############
message(paste0("--> Ploting... "))
source(paste0(path_code,"04.plot.R"))
  
#################
### END TIMER ###
#################
etime <- Sys.time()
message(paste0("---> Completed in ",round(difftime(etime,stime, units="mins"),1)," minutes."))



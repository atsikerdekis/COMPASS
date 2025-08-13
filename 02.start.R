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
seqDate   <- format(seq.Date(from=as.Date(sDate, format="%Y%m%d"), to=as.Date(eDate, format="%Y%m%d"), by="day"), "%Y%m%d")
path_code <- paste0(dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))), "/")
source(paste0(path_code,"config.R"))
source(paste0(path_code,"01.init.R"))

if (1==2) {
########################
### CHECK & DOWNLOAD ###
########################
for (d in 1:length(seqDate)) {
  for (e in 1:2) {
    if (e == 1) {expname=expname1;exptype=exptype1}
    if (e == 2) {expname=expname2;exptype=exptype2}
    message(paste0("---> Downloading ",expname," (",exptype,")..."))
    # Download surface (sfc)
    SubmitJob(
      JOB_name     = paste0("download_sfc_",expname,"_",seqDate[d]),
      JOB_out      = paste0(path_log,"download_sfc_",expname,"_",seqDate[d],".out"),
      JOB_err      = paste0(path_log,"download_sfc_",expname,"_",seqDate[d],".err"),
      PATH_program = paste0(path_PYTHON,"python"),
      PATH_script  = paste0(path_code,"03.download_sfc.py"),
      SCRIPT_flag  = paste0(expname," ",seqDate[d]," ",path_data)
    )
    # Download pressure level (HAM)
    if (exptype == "HAM") {
        SubmitJob(
          JOB_name     = paste0("download_pl_M7_",expname,"_",seqDate[d]),
          JOB_out      = paste0(path_log,"download_pl_M7_",expname,"_",seqDate[d],".out"),
          JOB_err      = paste0(path_log,"download_pl_M7_",expname,"_",seqDate[d],".err"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_pl_M7.py"),
          SCRIPT_flag  = paste0(expname," ",seqDate[d]," ",path_data)
      )
    }
    # Download pressure level (AER)
    if (exptype == "AER") {
        SubmitJob( 
          JOB_name     = paste0("download_pl_",expname,"_",seqDate[d]),
          JOB_out      = paste0(path_log,"download_pl_",expname,"_",seqDate[d],".out"),
          JOB_err      = paste0(path_log,"download_pl_",expname,"_",seqDate[d],".err"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_pl.py"),
          SCRIPT_flag  = paste0(expname," ",seqDate[d]," ",path_data)
      )   
    } 
  } # END LOOP experiments (e)
} # END LOOP dates (d)
}

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



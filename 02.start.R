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
path_code <- paste0(dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))), "/")
source(paste0(path_code,"config.R"))
source(paste0(path_code,"01.init.R"))
### for downloading...
seqDate    <- format(seq.Date(from=as.Date(sDate, format="%Y%m%d"), to=as.Date(eDate, format="%Y%m%d"), by="day"), "%Y%m%d")
step       <- round(length(seqDate)/NumberOfDownloadJobs,0) ; if (step < 1) {step <- 1}
sDate_temp <- as.Date(sDate, format="%Y%m%d")
eDate_temp <- as.Date(eDate, format="%Y%m%d")
vdate1     <- seq.Date(sDate_temp, eDate_temp, paste0(step," day"))
vdate1     <- gsub(pattern="-", replacement="", x=vdate1)
sDate_temp <- as.Date(sDate_temp)-1+step
eDate_temp <- as.Date(eDate_temp)-1+step
vdate2     <- seq.Date(sDate_temp, eDate_temp, paste0(step," day"))
vdate2     <- gsub(pattern="-", replacement="", x=vdate2)
#stop()

########################
### CHECK & DOWNLOAD ###
########################
if (runtype == "download") {
  message("---> Download mode...")
  for (d in 1:length(vdate1)) {
    for (e in 1:2) {
      if (e == 1) {expname=expname1;exptype=exptype1;expclass=expclass1}
      if (e == 2) {expname=expname2;exptype=exptype2;expclass=expclass2}
      message(paste0("---> Downloading ",expname," (",exptype,")..."))

#if (1==2) {
      # Download surface (sfc)
      SubmitJob(
        JOB_name     = paste0("download_sfc_",expname,"_",vdate1[d],"_",vdate2[d]),
        JOB_out      = paste0(path_log,"download_sfc_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
        JOB_err      = paste0(path_log,"download_sfc_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
        PATH_program = paste0(path_PYTHON,"python"),
        PATH_script  = paste0(path_code,"03.download_sfc.py"),
        SCRIPT_flag  = paste0(expname," ",expclass," ",vdate1[d]," ",vdate2[d]," ",path_data)
      )
      # Download pressure level (HAM)
      if (exptype == "HAM") {
        SubmitJob(
          JOB_name     = paste0("download_pl_M7_",expname,"_",vdate1[d],"_",vdate2[d]),
          JOB_out      = paste0(path_log,"download_pl_M7_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          JOB_err      = paste0(path_log,"download_pl_M7_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_pl_M7.py"),
          SCRIPT_flag  = paste0(expname," ",expclass," ",vdate1[d]," ",vdate2[d]," ",path_data)
        )
      }
      # Download pressure level (HAM)
      if (exptype == "HAM_NI_CS") {
        SubmitJob(
          JOB_name     = paste0("download_pl_M7_",expname,"_",vdate1[d],"_",vdate2[d]),
          JOB_out      = paste0(path_log,"download_pl_M7_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          JOB_err      = paste0(path_log,"download_pl_M7_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_pl_M7_Add_NI_CS.py"),
          SCRIPT_flag  = paste0(expname," ",expclass," ",vdate1[d]," ",vdate2[d]," ",path_data)
        )
      }
      # Download pressure level (AER)
      if (exptype == "AER") {
        SubmitJob(
          JOB_name     = paste0("download_pl_",expname,"_",vdate1[d],"_",vdate2[d]),
          JOB_out      = paste0(path_log,"download_pl_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          JOB_err      = paste0(path_log,"download_pl_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_pl.py"),
          SCRIPT_flag  = paste0(expname," ",expclass," ",vdate1[d]," ",vdate2[d]," ",path_data)
        )
      }
#}

    } # END LOOP experiments (e)
  } # END LOOP dates (d)
} # END IF runtype == "download"


############
### PLOT ###
############
if (runtype == "plot") {
  message("---> Plot mode...")
  source(paste0(path_code,"04.plot.R"))
} # END IF runtype == "plot"

#################
### END TIMER ###
#################
etime <- Sys.time()
message(paste0("---> Completed in ",round(difftime(etime,stime, units="mins"),1)," minutes."))



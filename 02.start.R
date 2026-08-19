###
### Created : Thanos Tsikerdekis (KNMI) | May 2025
### Contact : thanos.tsikerdekis@knmi.nl
### Purpose : Start a comparison
###

###################
### START TIMER ###
###################
stime <- Sys.time()
message("---> Starting...")

##################
### INITIALIZE ###
##################
path_code <- paste0(dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))), "/")
source(paste0(path_code,"config.R"))
source(paste0(path_code,"01.init.R"))

######################
### DOWNLOAD DATES ###
######################
seqDate <- format(seq.Date(from=as.Date(sDate,format="%Y%m%d"),to=as.Date(eDate,format="%Y%m%d"),by="day"),"%Y%m%d")

step <- round(length(seqDate)/NumberOfDownloadJobs,0)
if (step < 1) step <- 1

sDate_temp <- as.Date(sDate,format="%Y%m%d")
eDate_temp <- as.Date(eDate,format="%Y%m%d")

vdate1 <- gsub("-","",seq.Date(sDate_temp,eDate_temp,paste0(step," day")))

sDate_temp <- sDate_temp-1+step
eDate_temp <- eDate_temp-1+step

vdate2 <- gsub("-","",seq.Date(sDate_temp,eDate_temp,paste0(step," day")))

########################
### CHECK & DOWNLOAD ###
########################
if (runtype == "download") {

  message("---> Download mode...")

  for (d in seq_along(vdate1)) {

    for (e in 1:2) {

      if (e == 1) {
        expname  <- expname1
        exptype  <- exptype1
        expclass <- expclass1
        expvars  <- variables_exp1
      } else {
        expname  <- expname2
        exptype  <- exptype2
        expclass <- expclass2
        expvars  <- variables_exp2
      }

      message("---> Downloading ",expname," (",exptype,")...")

      #########################
      ### SELECT PARAMETERS ###
      #########################

      ### Pressure-level diagnostics
      params_pl <- unique(expvars$grib[expvars$grib_column == "grib"])
      params_pl <- params_pl[!is.na(params_pl) & params_pl != ""]
      params_pl <- paste(params_pl,collapse="/")

      ### Surface diagnostics
      surface_columns <- c("gribddp","gribsdm","gribwdl","gribwdc")
      params_sfc <- unique(expvars$grib[expvars$grib_column %in% surface_columns])
      params_sfc <- params_sfc[!is.na(params_sfc) & params_sfc != ""]
      params_sfc <- paste(params_sfc,collapse="/")

      message("---> PL parameters: ",ifelse(params_pl == "","none",params_pl))
      message("---> SFC parameters: ",ifelse(params_sfc == "","none",params_sfc))

      ########################
      ### SURFACE DOWNLOAD ###
      ########################
      if (params_sfc != "") {

        SubmitJob(
          JOB_name     = paste0("download_sfc_",expname,"_",vdate1[d],"_",vdate2[d]),
          JOB_out      = paste0(path_log,"download_sfc_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          JOB_err      = paste0(path_log,"download_sfc_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_sfc.py"),
          SCRIPT_flag  = paste(expname,expclass,vdate1[d],vdate2[d],path_data,params_sfc)
        )

      }

      ###############################
      ### PRESSURE LEVEL DOWNLOAD ###
      ###############################
      if (params_pl != "") {

        SubmitJob(
          JOB_name     = paste0("download_pl_",expname,"_",vdate1[d],"_",vdate2[d]),
          JOB_out      = paste0(path_log,"download_pl_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          JOB_err      = paste0(path_log,"download_pl_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_pl.py"),
          SCRIPT_flag  = paste(expname,expclass,vdate1[d],vdate2[d],path_data,params_pl)
        )

      }

    }

  }

}

############
### PLOT ###
############
if (runtype == "plot") {
  message("---> Plot mode...")
  source(paste0(path_code,"05.plot.R"))
}

#################
### END TIMER ###
#################
etime <- Sys.time()
message("---> Completed in ",round(difftime(etime,stime,units="mins"),1)," minutes.")

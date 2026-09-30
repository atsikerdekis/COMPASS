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

#############################
### ARCHIVE OLD DOWNLOADS ###
#############################
next_archive_dir <- function(exp_dir) {

  existing_dirs <- list.dirs(exp_dir,full.names=FALSE,recursive=FALSE)
  archive_dirs <- existing_dirs[grepl("^archive[0-9]{3}$",existing_dirs)]

  if (length(archive_dirs) == 0) {
    next_id <- 1
  } else {
    ids <- suppressWarnings(as.integer(sub("^archive","",archive_dirs)))
    ids <- ids[is.finite(ids)]
    next_id <- if (length(ids) == 0) 1 else max(ids)+1
  }

  file.path(exp_dir,paste0("archive",sprintf("%03d",next_id)))
}

archive_existing_download_files <- function(
  expname,dates,
  archive_sfc=FALSE,
  archive_precip=FALSE,
  archive_optics=FALSE,
  archive_mec=FALSE,
  archive_columns=FALSE,
  archive_pl=FALSE,
  archive_ml=FALSE
) {

  exp_dir <- file.path(path_data,expname)
  dir.create(exp_dir,recursive=TRUE,showWarnings=FALSE)

  targets <- character(0)

  for (date in dates) {

    if (archive_sfc)
      targets <- c(targets,file.path(
        exp_dir,paste0("CAMS_",expname,"_forecast00to21by03_0.7x0.7_sfc_",date,".nc")
      ))

    if (archive_precip)
      targets <- c(targets,file.path(
        exp_dir,paste0("CAMS_",expname,"_forecast03to24by03_0.7x0.7_sfc_precip_",date,".nc")
      ))

    if (archive_optics)
      targets <- c(targets,file.path(
        exp_dir,paste0("CAMS_",expname,"_forecast00to21by03_0.7x0.7_sfc_optics_",date,".nc")
      ))

    if (archive_mec)
      targets <- c(targets,file.path(
        exp_dir,paste0("CAMS_",expname,"_forecast00to21by03_0.7x0.7_sfc_mec_",date,".nc")
      ))

    if (archive_columns)
      targets <- c(targets,file.path(
        exp_dir,paste0("CAMS_",expname,"_forecast00to21by03_0.7x0.7_sfc_columns_",date,".nc")
      ))

    if (archive_pl)
      targets <- c(targets,file.path(
        exp_dir,paste0("CAMS_",expname,"_forecast00to21by03_0.7x0.7_pl_",date,".nc")
      ))

    if (archive_ml) {
      targets <- c(
        targets,
        file.path(exp_dir,paste0("CAMS_",expname,"_forecast00to21by03_0.7x0.7_ml_",date,".nc")),
        file.path(exp_dir,paste0("CAMS_",expname,"_forecast00to21by03_0.7x0.7_lnsp_",date,".nc"))
      )
    }
  }

  existing <- unique(targets[file.exists(targets)])

  if (length(existing) == 0) {
    message("---> No existing download files need archiving for ",expname,".")
    return(invisible(NULL))
  }

  archive_dir <- next_archive_dir(exp_dir)
  dir.create(archive_dir,recursive=TRUE,showWarnings=FALSE)

  message(
    "---> Existing files detected for ",expname,": ",length(existing),
    " file(s). Moving them to ",archive_dir," before fresh download."
  )

  for (src in existing) {
    dst <- file.path(archive_dir,basename(src))
    ok <- file.rename(src,dst)
    if (!ok)
      stop("Could not archive existing file: ",src," -> ",dst)
    message("     archived: ",basename(src))
  }

  message("---> Archive complete: ",archive_dir)
  invisible(archive_dir)
}

########################
### CHECK & DOWNLOAD ###
########################
if (runtype == "download") {

  message("---> Download mode...")

  if (exists("massdiag_compare") && massdiag_compare) {
    for (e in 1:2) {
      if (e == 1) expname <- expname1 else expname <- expname2

      SubmitJob(
        JOB_name     = paste0("download_massdiag_",expname),
        JOB_out      = paste0(path_log,"download_massdiag_",expname,".out"),
        JOB_err      = paste0(path_log,"download_massdiag_",expname,".out"),
        PATH_program = paste0(path_R,"Rscript"),
        PATH_script  = paste0(path_code,"03.download_massdiag.R"),
        SCRIPT_flag  = paste(expname,sDate,eDate,path_data)
      )
    }
  }

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

      ### Pressure-level aerosol diagnostics
      params_pl <- unique(expvars$grib[expvars$grib_column == "grib"])
      params_pl <- params_pl[!is.na(params_pl) & params_pl != ""]
      params_pl <- paste(params_pl,collapse="/")

      ### Surface aerosol diagnostics
      surface_columns <- c("gribddp","gribsdm","gribwdl","gribwdc","gribmss","gribod","gribngt")
      params_sfc <- unique(expvars$grib[expvars$grib_column %in% surface_columns])
      params_sfc <- params_sfc[!is.na(params_sfc) & params_sfc != ""]
      params_sfc <- paste(params_sfc,collapse="/")

      ### Model levels required for derived RH diagnostics
      levels_ml <- ""
      if (length(rh_variables) > 0) {
        levels_ml <- unique(sapply(rh_variables,function(x) get_rh_definition(x)$model_level))
        levels_ml <- paste(sort(levels_ml),collapse="/")
      }

      ### Surface precipitation diagnostics
      params_precip <- ""
      if (length(precip_variables) > 0) {
        params_precip <- unique(sapply(precip_variables,function(x) get_precip_definition(x)$grib))
        params_precip <- paste(params_precip,collapse="/")
      }


### Optical-property diagnostics
params_optics <- ""
if (length(optical_variables) > 0) {
  params_optics <- paste(get_optical_required_gribs(optical_variables),collapse="/")
}


### Surface total-column gas diagnostics
params_columns <- ""
if (length(column_variables) > 0) {
  params_columns <- unique(sapply(column_variables,function(x) get_column_definition(x)$grib))
  params_columns <- paste(params_columns,collapse="/")
}

### MEC needs total aerosol column burden. Use all available mss GRIBs.
params_mec <- ""
if ("mec550" %in% strip_time_aggregation(optical_variables)) {
  params_mec <- unique(expvars$grib[expvars$grib_column == "gribmss"])
  params_mec <- params_mec[!is.na(params_mec) & params_mec != ""]
  ### If no MSS variable was explicitly requested, resolve the full aerosol burden.
  if (length(params_mec) == 0) {
    mec_names <- paste0("mss_",c("ss","du","pom","bc","so4","ni","am","soa"))
    mec_table <- data.frame()
    for (nm in mec_names) {
      tmp <- tryCatch(resolve_variables(exptype,nm),error=function(e) data.frame())
      if (nrow(tmp) > 0) mec_table <- rbind(mec_table,tmp)
    }
    if (nrow(mec_table) > 0) params_mec <- unique(mec_table$grib[mec_table$grib_column == "gribmss"])
  }
  params_mec <- paste(params_mec,collapse="/")
}

      message("---> PL parameters: ",ifelse(params_pl == "","none",params_pl))
      message("---> SFC parameters: ",ifelse(params_sfc == "","none",params_sfc))
      message("---> RH model levels: ",ifelse(levels_ml == "","none",levels_ml))
      message("---> Precipitation parameters: ",ifelse(params_precip == "","none",params_precip))
      message("---> Optical parameters: ",ifelse(params_optics == "","none",params_optics))
      message("---> MEC burden parameters: ",ifelse(params_mec == "","none",params_mec))
      message("---> Total-column parameters: ",ifelse(params_columns == "","none",params_columns))

      ### Archive existing target files once per experiment for this run.
      ### Doing this before any child jobs are submitted avoids race conditions
      ### between parallel download jobs choosing archiveNNN directories.
      if (d == 1) {
        archive_existing_download_files(
          expname=expname,
          dates=seqDate,
          archive_sfc=params_sfc != "",
          archive_precip=params_precip != "",
          archive_optics=params_optics != "",
          archive_mec=params_mec != "",
          archive_columns=params_columns != "",
          archive_pl=params_pl != "",
          archive_ml=levels_ml != ""
        )
      }

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

      ################################
      ### PRECIPITATION DOWNLOAD ###
      ################################
      if (params_precip != "") {

        SubmitJob(
          JOB_name     = paste0("download_precip_",expname,"_",vdate1[d],"_",vdate2[d]),
          JOB_out      = paste0(path_log,"download_precip_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          JOB_err      = paste0(path_log,"download_precip_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_sfc.py"),
          SCRIPT_flag  = paste(expname,expclass,vdate1[d],vdate2[d],path_data,params_precip,"3/6/9/12/15/18/21/24","sfc_precip")
        )

      }


############################
### OPTICAL DOWNLOAD ###
############################
if (params_optics != "") {
  SubmitJob(
    JOB_name     = paste0("download_optics_",expname,"_",vdate1[d],"_",vdate2[d]),
    JOB_out      = paste0(path_log,"download_optics_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
    JOB_err      = paste0(path_log,"download_optics_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
    PATH_program = paste0(path_PYTHON,"python"),
    PATH_script  = paste0(path_code,"03.download_sfc.py"),
    SCRIPT_flag  = paste(expname,expclass,vdate1[d],vdate2[d],path_data,params_optics,"0/3/6/9/12/15/18/21","sfc_optics")
  )
}

############################
### MEC BURDEN DOWNLOAD ###
############################
if (params_mec != "") {
  SubmitJob(
    JOB_name     = paste0("download_mec_",expname,"_",vdate1[d],"_",vdate2[d]),
    JOB_out      = paste0(path_log,"download_mec_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
    JOB_err      = paste0(path_log,"download_mec_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
    PATH_program = paste0(path_PYTHON,"python"),
    PATH_script  = paste0(path_code,"03.download_sfc.py"),
    SCRIPT_flag  = paste(expname,expclass,vdate1[d],vdate2[d],path_data,params_mec,"0/3/6/9/12/15/18/21","sfc_mec")
  )
}


################################
### TOTAL-COLUMN GAS DOWNLOAD ###
################################
if (params_columns != "") {
  SubmitJob(
    JOB_name     = paste0("download_columns_",expname,"_",vdate1[d],"_",vdate2[d]),
    JOB_out      = paste0(path_log,"download_columns_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
    JOB_err      = paste0(path_log,"download_columns_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
    PATH_program = paste0(path_PYTHON,"python"),
    PATH_script  = paste0(path_code,"03.download_sfc.py"),
    SCRIPT_flag  = paste(expname,expclass,vdate1[d],vdate2[d],path_data,params_columns,"0/3/6/9/12/15/18/21","sfc_columns")
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

      ############################
      ### RH MODEL-LEVEL DOWNLOAD ###
      ############################
      if (levels_ml != "") {

        SubmitJob(
          JOB_name     = paste0("download_ml_rh_",expname,"_",vdate1[d],"_",vdate2[d]),
          JOB_out      = paste0(path_log,"download_ml_rh_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          JOB_err      = paste0(path_log,"download_ml_rh_",expname,"_",vdate1[d],"_",vdate2[d],".out"),
          PATH_program = paste0(path_PYTHON,"python"),
          PATH_script  = paste0(path_code,"03.download_ml.py"),
          SCRIPT_flag  = paste(expname,expclass,vdate1[d],vdate2[d],path_data,levels_ml)
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

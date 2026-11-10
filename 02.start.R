###
### Created : Thanos Tsikerdekis (KNMI) | May 2025
### Purpose : Start a comparison; batch MARS downloads preserving daily NetCDF files
###

stime <- Sys.time()
message("---> Starting...")

path_code <- paste0(dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))), "/")
source(paste0(path_code,"config.R"))
source(paste0(path_code,"01.init.R"))

seqDate <- format(seq.Date(as.Date(sDate,"%Y%m%d"),as.Date(eDate,"%Y%m%d"),by="day"),"%Y%m%d")
if (length(seqDate) < 1) stop("Invalid download date range: ",sDate,"-",eDate)

### Maximum number of consecutive calendar days in one queued MARS job.
### DownloadBatchDays can be added to config.R; 30 is the default.
if (!exists("DownloadBatchDays")) DownloadBatchDays <- 30L
DownloadBatchDays <- as.integer(DownloadBatchDays)
if (length(DownloadBatchDays) != 1 || is.na(DownloadBatchDays) || DownloadBatchDays < 1 || DownloadBatchDays > 30)
  stop("DownloadBatchDays must be an integer between 1 and 30")

batch_ids <- split(seq_along(seqDate),ceiling(seq_along(seqDate)/DownloadBatchDays))

### Check existing output paths before scheduling; the worker also verifies
### NetCDF timestep counts and skips completed outputs within partial batches.
batch_missing <- function(expname,prefix,dates,params) {
  files <- file.path(path_data,expname,paste0(prefix,"_",dates,".nc"))
  expected <- vapply(strsplit(params,"/",fixed=TRUE)[[1]],function(code) {
    bits <- strsplit(code,".",fixed=TRUE)[[1]]
    paste0("p",sprintf("%03d",as.integer(bits[2])),sprintf("%03d",as.integer(bits[1])))
  },character(1))
  for (file in files) {
    if (!file.exists(file) || is.na(file.info(file)$size) || file.info(file)$size <= 0)
      return(TRUE)
    hdr <- tryCatch(system2("ncdump",c("-h",shQuote(file)),stdout=TRUE,stderr=FALSE),
                    error=function(e) character(0))
    if (length(hdr) == 0 || !all(vapply(expected,function(x)
      any(grepl(paste0(x,"\\("),hdr)),logical(1)))) return(TRUE)
    ntime <- tryCatch(system2("/usr/local/apps/cdo/2.5.1/bin/cdo",
                                   c("-s","ntime",shQuote(file)),
                                   stdout=TRUE,stderr=FALSE),
                      error=function(e) character(0))
    if (length(ntime) != 1 || !isTRUE(suppressWarnings(as.integer(ntime[1])) == 8L))
      return(TRUE)
  }
  FALSE
}

submit_batch <- function(job_tag,script,expname,expclass,date1,date2,
                         params,steps="",out_tag="",levels="") {
  if (is.null(params) || !nzchar(params)) return(invisible(NULL))
  base <- paste0("CAMS_",expname,"_forecast",
                 if (job_tag == "precip") "03to24" else "00to21",
                 "by03_0.7x0.7_")
  suffix <- switch(job_tag,
                   pl="pl",sfc="sfc",precip="sfc_precip",
                   optics="sfc_optics",mec="sfc_mec",columns="sfc_columns",
                   ml_rh="ml_rh",ml_wat="ml_wat",job_tag)
  ### RH additionally creates a dedicated lnsp file. Its existence alone is
  ### not enough to suppress a needed RH job, and vice versa.
  dates <- format(seq.Date(as.Date(date1,"%Y%m%d"),as.Date(date2,"%Y%m%d"),by="day"),"%Y%m%d")
  missing <- batch_missing(expname,paste0(base,suffix),dates,params)
  if (job_tag == "ml_rh") missing <- missing ||
    batch_missing(expname,paste0(base,"lnsp"),dates,"152.128")
  if (!missing) {
    message("---> All daily files exist for ",expname," / ",job_tag," / ",date1,"-",date2,"; skipping job.")
    return(invisible(NULL))
  }

  name <- paste0("download_",job_tag,"_",expname,"_",date1,"_",date2)
  flags <- if (job_tag %in% c("ml_rh","ml_wat")) {
    paste(expname,expclass,date1,date2,path_data,levels,params,out_tag)
  } else if (job_tag == "pl") {
    paste(expname,expclass,date1,date2,path_data,params)
  } else if (job_tag == "sfc") {
    paste(expname,expclass,date1,date2,path_data,params)
  } else {
    paste(expname,expclass,date1,date2,path_data,params,steps,out_tag)
  }
  message("---> Submitting ",name," (up to ",length(dates)," days, missing outputs only)")
  SubmitJob(
    JOB_name=name,
    JOB_out=paste0(path_log,name,".out"),
    JOB_err=paste0(path_log,name,".out"),
    PATH_program=paste0(path_PYTHON,"python"),
    PATH_script=paste0(path_code,script),
    SCRIPT_flag=flags
  )
}

if (runtype == "download") {
  message("---> Download mode: GRIB batches of up to ",DownloadBatchDays," days.")
  message("---> Existing valid daily NetCDF files will not be archived or overwritten.")

  if (exists("massdiag_compare") && isTRUE(massdiag_compare)) {
    for (expname in unique(c(expname1,expname2))) {
      SubmitJob(
        JOB_name=paste0("download_massdiag_",expname),
        JOB_out=paste0(path_log,"download_massdiag_",expname,".out"),
        JOB_err=paste0(path_log,"download_massdiag_",expname,".out"),
        PATH_program=paste0(path_R,"Rscript"),
        PATH_script=paste0(path_code,"03.download_massdiag.R"),
        SCRIPT_flag=paste(expname,sDate,eDate,path_data)
      )
    }
  }

  for (e in 1:2) {
    if (e == 1) {
      expname <- expname1; exptype <- exptype1
      expclass <- expclass1; expvars <- variables_exp1
    } else {
      expname <- expname2; exptype <- exptype2
      expclass <- expclass2; expvars <- variables_exp2
    }

    ### Pressure-level aerosol diagnostics.
    params_pl <- unique(expvars$grib[expvars$grib_column == "grib"])
    params_pl <- params_pl[!is.na(params_pl) & nzchar(params_pl)]
    params_pl <- paste(params_pl,collapse="/")

    ### Surface aerosol diagnostics, including additional species OD/MSS fields.
    surface_columns <- c("gribddp","gribsdm","gribwdl","gribwdc","gribmss","gribod","gribngt")
    params_sfc <- unique(expvars$grib[expvars$grib_column %in% surface_columns])
    params_sfc <- params_sfc[!is.na(params_sfc) & nzchar(params_sfc)]
    panel_bases <- tolower(strip_time_aggregation(panel_species_variables))

    if (any(strip_time_aggregation(optical_variables) %in% c("od_sum_species","aod550_species")) ||
        any(panel_bases %in% c("aod_per_species","aodratio_per_species",
                                  "aod_per_species_residual_wat",
                                  "aodratio_per_species_residual_wat","mec_per_species"))) {
      species_vars <- paste0("od_",get_od_species_suffixes(include_soa=TRUE))
      tab <- resolve_variables(exptype,species_vars)
      extra <- unique(tab$grib[tab$grib_column == "gribod"])
      params_sfc <- unique(c(params_sfc,extra[!is.na(extra) & nzchar(extra)]))
    }

    if (any(panel_bases %in% c("mass_per_species","massratio_per_species","mec_per_species"))) {
      for (nm in paste0("mss_",c("ss","du","pom","bc","so4","ni","am","soa"))) {
        tab <- tryCatch(resolve_variables(exptype,nm),error=function(err) data.frame())
        if (nrow(tab)) params_sfc <- c(params_sfc,tab$grib[tab$grib_column == "gribmss"])
      }
    }
    params_sfc <- unique(params_sfc)
    params_sfc <- paste(params_sfc[!is.na(params_sfc) & nzchar(params_sfc)],collapse="/")

    levels_ml <- ""
    if (length(rh_variables)) {
      levels_ml <- paste(sort(unique(sapply(rh_variables,function(x)
        get_rh_definition(x)$model_level))),collapse="/")
    }
    params_ml_wat <- ""
    levels_ml_wat <- ""
    if (length(wat_variables) || (exptype == "HAM" &&
       any(panel_bases %in% c("aod_per_species","aodratio_per_species")))) {
      levels_ml_wat <- "2/3/4"
      params_ml_wat <- "22.210"
    }

    params_precip <- ""
    if (length(precip_variables))
      params_precip <- paste(unique(sapply(precip_variables,function(x)
        get_precip_definition(x)$grib)),collapse="/")

    params_optics <- ""
    if (length(optical_variables) || length(panel_species_variables)) {
      optical_need <- optical_variables
      if (any(panel_bases %in% c("aod_per_species","aodratio_per_species",
                                  "aod_per_species_residual_wat",
                                  "aodratio_per_species_residual_wat","mec_per_species")))
        optical_need <- unique(c(optical_need,"aod550"))
      if (length(panel_species_variables))
        optical_need <- unique(c(optical_need,get_panel_validation_definition()$model_variable))
      params_optics <- paste(get_optical_required_gribs(optical_need),collapse="/")
    }

    params_columns <- ""
    if (length(column_variables))
      params_columns <- paste(unique(sapply(column_variables,function(x)
        get_column_definition(x)$grib)),collapse="/")

    params_mec <- ""
    if ("mec550" %in% strip_time_aggregation(optical_variables) ||
        "mec_per_species" %in% panel_bases) {
      mec_gribs <- unique(expvars$grib[expvars$grib_column == "gribmss"])
      mec_gribs <- mec_gribs[!is.na(mec_gribs) & nzchar(mec_gribs)]
      if (!length(mec_gribs)) {
        for (nm in paste0("mss_",c("ss","du","pom","bc","so4","ni","am","soa"))) {
          tab <- tryCatch(resolve_variables(exptype,nm),error=function(err) data.frame())
          if (nrow(tab)) mec_gribs <- c(mec_gribs,tab$grib[tab$grib_column == "gribmss"])
        }
      }
      params_mec <- paste(unique(mec_gribs[!is.na(mec_gribs) & nzchar(mec_gribs)]),collapse="/")
    }

    message("---> ",expname," (",exptype,") PL: ",ifelse(nzchar(params_pl),params_pl,"none"))
    message("---> ",expname," SFC: ",ifelse(nzchar(params_sfc),params_sfc,"none"))
    message("---> ",expname," OPTICS: ",ifelse(nzchar(params_optics),params_optics,"none"))

    for (ids in batch_ids) {
      date1 <- seqDate[min(ids)]
      date2 <- seqDate[max(ids)]

      submit_batch("sfc","03.download_sfc.py",expname,expclass,date1,date2,params_sfc)
      submit_batch("precip","03.download_sfc.py",expname,expclass,date1,date2,
                   params_precip,"3/6/9/12/15/18/21/24","sfc_precip")
      submit_batch("optics","03.download_sfc.py",expname,expclass,date1,date2,
                   params_optics,"0/3/6/9/12/15/18/21","sfc_optics")
      submit_batch("mec","03.download_sfc.py",expname,expclass,date1,date2,
                   params_mec,"0/3/6/9/12/15/18/21","sfc_mec")
      submit_batch("columns","03.download_sfc.py",expname,expclass,date1,date2,
                   params_columns,"0/3/6/9/12/15/18/21","sfc_columns")
      submit_batch("pl","03.download_pl.py",expname,expclass,date1,date2,params_pl)
      if (nzchar(levels_ml)) submit_batch("ml_rh","03.download_ml.py",expname,expclass,date1,date2,
                   "130.128/133.128",out_tag="rh",levels=levels_ml)
      ### Only request RH model levels when configured.
      ### Avoid a zero-level RH request when no RH variable is selected.
      submit_batch("ml_wat","03.download_ml.py",expname,expclass,date1,date2,
                   params_ml_wat,out_tag="wat",levels=levels_ml_wat)
    }
  }
}

if (runtype == "plot") {
  message("---> Plot mode...")
  source(paste0(path_code,"05.plot.R"))
}

etime <- Sys.time()
message("---> Completed in ",round(difftime(etime,stime,units="mins"),1)," minutes.")

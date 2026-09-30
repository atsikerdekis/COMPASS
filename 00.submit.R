###
### Created : Thanos Tsikerdekis (KNMI) | Jun 2024 |
### Contact : thanos.tsikerdekis@knmi.nl
### Purpose : Submit COMPASS work
###

############
### INIT ###
############
path_code <- paste0(
  dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))),
  "/"
)

source(paste0(path_code,"config.R"))

username <- system("whoami",intern=TRUE)
path_log <- paste0("/perm/",username,"/cams2_35/COMPASS/log/")
path_R   <- paste0("/etc/ecmwf/nfs/dh1_perm_b/",username,"/miniforge3/envs/COMPASS/bin/")

dir.create(path_log,recursive=TRUE,showWarnings=FALSE)

path_script <- paste0(path_code,"02.start.R")

########################################
### DOWNLOAD: SUBMIT CHILD JOBS HERE ###
########################################
### 02.start.R is lightweight in download mode: it prepares/archive-checks
### files and submits the actual MARS download jobs to SLURM. Running it here
### means the child job IDs and stdout/stderr paths are visible immediately
### in the terminal instead of being hidden inside a parent SLURM .out file.
if (runtype == "download") {

  message("--> Download mode: submitting download jobs directly...")

  status <- system2(
    paste0(path_R,"Rscript"),
    args=shQuote(path_script)
  )

  if (!identical(status,0L))
    stop("02.start.R failed while submitting download jobs. Exit status: ",status)

  message("--> All requested download jobs have been submitted.")
  quit(save="no",status=0)
}

####################################
### OTHER MODES: SUBMIT ONE JOB ###
####################################
script_name <- paste0(
  expname1," ",exptype1," ",
  expname2," ",exptype2," ",
  sDate," ",eDate," ",runtype
)

script_flag <- ""
file.slurm <- paste0(path_log,"Job_",gsub(" ","_",script_name),".slurm")
file.out   <- paste0(path_log,"Job_",gsub(" ","_",script_name),".out")

slurm <- readLines(paste0(path_code,"function/TemplateJob.slurm"))
slurm <- gsub(pattern="JobName",      replacement=paste0("Job_",gsub(" ","_",script_name)),slurm)
slurm <- gsub(pattern="SCRIPT_flag",  replacement=script_flag,slurm)
slurm <- gsub(pattern="Job.out",      replacement=file.out,slurm)
slurm <- gsub(pattern="Job.err",      replacement=file.out,slurm)
slurm <- gsub(pattern="PATH_program", replacement=paste0(path_R,"Rscript"),slurm)
slurm <- gsub(pattern="PATH_script",  replacement=path_script,slurm)

writeLines(slurm,file.slurm)

message("--> Submit: ",file.slurm)
message("--> Output: ",file.out)

submit_output <- system2(
  "sbatch",
  args=shQuote(file.slurm),
  stdout=TRUE,
  stderr=TRUE
)

status <- attr(submit_output,"status")
if (!is.null(status) && status != 0)
  stop("sbatch failed:\n",paste(submit_output,collapse="\n"))

message(paste(submit_output,collapse="\n"))

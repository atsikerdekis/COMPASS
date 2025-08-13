###
### Created : Thanos Tsikerdekis (KNMI) | Jun 2024 |
### Contact : thanos.tsikerdekis@knmi.nl
### Purpose : Submits a Job per day given 2 or 3 arguments
### Example : Rscript 00.SubmitJob.R MOD 20240101 20240102
### Envirom : ...
###

############
### INIT ###
############
args      <- commandArgs(trailingOnly=TRUE)
expname1  <- args[1]
exptype1  <- args[2]
expname2  <- args[3]
exptype2  <- args[4]
sDate     <- args[5]
eDate     <- args[6]
path_code <- paste0(dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))), "/")
source(paste0(path_code,"01.init.R"))

#############
### SLURM ### Create slurm file and submit job
#############
path_script  <- paste0(path_code,"02.start.R")
script_flag  <- paste0(expname1," ",exptype1," ",expname2," ",exptype2," ",sDate," ",eDate)
file.slurm   <- paste0(path_log,"Job_",gsub(" ","_",script_flag),".slurm")
file.out     <- paste0(path_log,"Job_",gsub(" ","_",script_flag),".out")
slurm <- readLines(paste0(path_code,"function/TemplateJob.slurm"))
slurm <- gsub(pattern="JobName",      replacement=paste0("Job_",gsub(" ","_",script_flag)), slurm)
slurm <- gsub(pattern="SCRIPT_flag",  replacement=script_flag, slurm)
slurm <- gsub(pattern="Job.out",      replacement=file.out, slurm)
slurm <- gsub(pattern="Job.err",      replacement=file.out, slurm)
slurm <- gsub(pattern="PATH_program", replacement=paste0(path_R,"Rscript"), slurm)
slurm <- gsub(pattern="PATH_script",  replacement=path_script, slurm)
writeLines(slurm, file.slurm)
message(paste0("--> Submit: sbatch ",file.slurm))
system(paste0("sbatch ",file.slurm))


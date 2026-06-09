###
### Created : Thanos Tsikerdekis (KNMI) | Jun 2024 |
### Contact : thanos.tsikerdekis@knmi.nl
### Purpose : ...
### Example : Rscript 00.submit.R b2t0 HAM b2t0 HAM 20230601 20230602 download
### Envirom : ...
###

############
### INIT ###
############
path_code <- paste0(dirname(normalizePath(sub("--file=", "", commandArgs(FALSE)[grep("--file=", commandArgs(FALSE))]))), "/")
source(paste0(path_code,"config.R"))
source(paste0(path_code,"01.init.R"))

#############
### SLURM ### Create slurm file and submit job
#############
path_script  <- paste0(path_code,"02.start.R")
script_name  <- paste0(expname1," ",exptype1," ",expname2," ",exptype2," ",sDate," ",eDate," ",runtype)
script_flag  <- ""
file.slurm   <- paste0(path_log,"Job_",gsub(" ","_",script_name),".slurm")
file.out     <- paste0(path_log,"Job_",gsub(" ","_",script_name),".out")
slurm <- readLines(paste0(path_code,"function/TemplateJob.slurm"))
slurm <- gsub(pattern="JobName",      replacement=paste0("Job_",gsub(" ","_",script_name)), slurm)
slurm <- gsub(pattern="SCRIPT_flag",  replacement=script_flag, slurm)
slurm <- gsub(pattern="Job.out",      replacement=file.out, slurm)
slurm <- gsub(pattern="Job.err",      replacement=file.out, slurm)
slurm <- gsub(pattern="PATH_program", replacement=paste0(path_R,"Rscript"), slurm)
slurm <- gsub(pattern="PATH_script",  replacement=path_script, slurm)
writeLines(slurm, file.slurm)
message(paste0("--> Submit: ",file.slurm))
message(paste0("--> Output: ",file.out))
system(paste0("sbatch ",file.slurm))


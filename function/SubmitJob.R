### SUBMIT JOB
SubmitJob <- function(JOB_name, JOB_out, JOB_err, PATH_program, PATH_script, SCRIPT_flag) {
  ### Create a slurm file
  writeLines(con=paste0(path_log,JOB_name,".slurm"),
    text=paste0("#! /bin/sh\n\n#SBATCH --job-name=",JOB_name,"\n#SBATCH --output=",JOB_out,"\n#SBATCH --error=",JOB_err,"\n#SBATCH --time=47:00:00\n#SBATCH --mem=8G\n\n",PATH_program," ",PATH_script," ",SCRIPT_flag,"\n"
  ))
  ### Submit the Job
  system(paste0("sbatch ",path_log,JOB_name,".slurm"))
}


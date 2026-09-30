### SUBMIT JOB
SubmitJob <- function(JOB_name,JOB_out,JOB_err,PATH_program,PATH_script,SCRIPT_flag) {

  slurm_file <- paste0(path_log,JOB_name,".slurm")

  ### Create SLURM file
  writeLines(
    con=slurm_file,
    text=paste0(
      "#! /bin/sh\n\n",
      "#SBATCH --job-name=",JOB_name,"\n",
      "#SBATCH --output=",JOB_out,"\n",
      "#SBATCH --error=",JOB_err,"\n",
      "#SBATCH --time=47:00:00\n",
      "#SBATCH --mem=8G\n\n",
      PATH_program," ",PATH_script," ",SCRIPT_flag,"\n"
    )
  )

  ### Submit and capture SLURM response
  submit_output <- system2(
    "sbatch",
    args=shQuote(slurm_file),
    stdout=TRUE,
    stderr=TRUE
  )

  status <- attr(submit_output,"status")
  if (!is.null(status) && status != 0) {
    stop(
      "Failed to submit ",JOB_name,":\n",
      paste(submit_output,collapse="\n")
    )
  }

  response <- if (length(submit_output) > 0) submit_output[1] else ""
  job_id <- sub("^Submitted batch job[[:space:]]+","",response)

  message(
    "---> Submitted ",JOB_name,
    " | Job ID: ",job_id,
    "\n     SLURM:  ",slurm_file,
    "\n     stdout: ",JOB_out,
    "\n     stderr: ",JOB_err
  )

  invisible(job_id)
}

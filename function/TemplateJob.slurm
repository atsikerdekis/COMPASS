#! /bin/sh

#SBATCH --job-name=JobName
#SBATCH --output=Job.out
#SBATCH --error=Job.err
#SBATCH --time=24:00:00
#SBATCH --mem=8G

PATH_program PATH_script SCRIPT_flag 


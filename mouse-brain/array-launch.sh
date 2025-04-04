#!/bin/bash
#SBATCH -p short # partition (queue)
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code)
#SBATCH -n 1 # number of cores
#SBATCH --array=1-1000 # array job indices
#SBATCH --mem=8192 # total memory
#SBATCH --job-name="runLR" # job name
#SBATCH -o /home/aeoz001/a-s-omics/data-outputs/out-err/slurm-%x-%A-%4a.out # STDOUT %a job array index number
#SBATCH -e /home/aeoz001/a-s-omics/data-outputs/out-err/slurm-%x-%A-%4a.err # STDERR %a job array index number
#SBATCH --mail-user=aeoz001@bucknell.edu # address to email
#SBATCH --mail-type=NONE # mail events (NONE, BEGIN, END, FAIL, ALL)
Rscript --vanilla /home/aeoz001/a-s-omics/code-spatial-smoothing/array-run-job.R 2024-03-14 0
## launch from a-s-omics with
## sbatch ./code-spatial-smoothing/array-launch.sh

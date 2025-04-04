#!/bin/bash
#SBATCH -p short # partition (queue)
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code)
#SBATCH -n 1 # number of cores
#SBATCH --array=1-1000 # array job indices
#SBATCH --mem=8192 # total memory
#SBATCH --job-name="runLR" # job name
#SBATCH -o /home/aeoz001/a-s-omics/data-outputs/mouse-embryo/out-err/slurm-%x-%A-%4a.out # STDOUT %a job array index number
#SBATCH -e /home/aeoz001/a-s-omics/data-outputs/mouse-embryo/out-err/slurm-%x-%A-%4a.err # STDERR %a job array index number
#SBATCH --mail-user=aeoz001@bucknell.edu # address to email
#SBATCH --mail-type=NONE # mail events (NONE, BEGIN, END, FAIL, ALL)
## arg 1: timestamp for output folder
## arg 2: index of embryo timepoint to process (1-8)
## arg 3: array maxes out at 1000, this will argument is added to array index to process more than 1000 jobs. eg if set to 1000, index 0 will be interpretted as 1001 job, etc
Rscript --vanilla /home/aeoz001/a-s-omics/code-spatial-smoothing/mouse-embryo/array-run-job.R 2024-07-16 7 0
## launch from a-s-omics with
## sbatch ./code-spatial-smoothing/array-launch.sh

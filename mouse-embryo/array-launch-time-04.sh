#!/bin/bash
#SBATCH -p short # partition (queue)
#SBATCH -N 1 # (leave at 1 unless using multi-node specific code)
#SBATCH -n 4 # number of cores
#SBATCH --array=1-1619 # array job indices # 1619 for embryo LR cut at 5% thresh
#SBATCH --mem=32768 # total memory
#SBATCH --job-name="runLR-4" # job name
#SBATCH -o /home/aeoz001/a-s-omics/data-outputs/mouse-embryo/out-err/slurm-%x-%A-%4a.out # STDOUT %a job array index number
#SBATCH -e /home/aeoz001/a-s-omics/data-outputs/mouse-embryo/out-err/slurm-%x-%A-%4a.err # STDERR %a job array index number
#SBATCH --mail-user=aeoz001@bucknell.edu # address to email
#SBATCH --mail-type=NONE # mail events (NONE, BEGIN, END, FAIL, ALL)
## arg 1: timestamp for output folder
## arg 2: index of embryo timepoint to process (1-8)
## arg 3: array maxes out at 3000, this will argument is added to array index to process more than 3000 jobs. eg if set to 3000, index 0 will be interpretted as 1001 job, etc
Rscript --vanilla /home/aeoz001/a-s-omics/code-spatial-smoothing/mouse-embryo/array-run-job.R 2024-08-27 4 0
## launch from a-s-omics with
## sbatch ./code-spatial-smoothing/array-launch.sh

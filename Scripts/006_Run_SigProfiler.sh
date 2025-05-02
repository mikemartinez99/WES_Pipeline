#!/bin/bash
#SBATCH --job-name=Subtract
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=Subtract_%j.out

#----- START
echo "Starting job: $SLURM_JOB_NAME (Job ID: $SLURM_JOB_ID)"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo -e "#-----------------------------------------------------------#\n"

#----- Source and activate conda
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate  /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/SigProfilerAssignment

#----- Move into SigProfiler directory
cd /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/curiel/20250206_Mouse_Exomes/code/snakemake/SigProfilerAssignment

#----- Run SigProfilerAssignment
for i in `cat mafs.list`; do python sig_profiler_assignment.py -i "$i" -o "$i"-outputs --write-results-per-sample --genome-build mm10;done

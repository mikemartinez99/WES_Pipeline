#!/bin/bash
#SBATCH --job-name=WES_Pipeline
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=WES_Pipeline_%j.out

source /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/etc/profile.d/conda.sh
conda activate /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/snakemake

#snakemake -s Snakefile \
#    --use-conda \
#    --conda-frontend conda \
#    --conda-prefix /dartfs/rc/nosnapshots/G/GMBSR_refs/envs/DAC-RNAseq-pipeline \
#	--rerun-incomplete \
#    --keep-going \
 #   --profile cluster_profile 


snakemake --rulegraph | dot -Tpng > rulegraph.png


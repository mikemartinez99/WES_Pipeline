#! /usr/bin/bash

#!/bin/bash
#SBATCH --job-name=GATK
#SBATCH --nodes=1
#SBATCH --partition=standard
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=GATK_%j.out

#----- Set paths
GATK_PATH="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/sullivan/tools/gatk/gatk-4.5.0.0/gatk"
JAVA_PATH="/usr/bin/java"
TARGETS="Twist_Mouse_Exome_Target_Rev1_7APR20.bed"
REF="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/mouse/UCSC_mm10/mm10.fa"

#----- Source conda
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/gatk4

#----- Create sequence dictionary
gatk CreateSequenceDictionary -R "$REF" -O mm10.dict

#----- Convert to interval list format
gatk BedToIntervalList \
    -I "$TARGETS" \
    -O mm10.interval_list \
    -SD mm10.dict

#----- Split bed file into regions, add regions to list
split -l 2155 "$TARGETS"
ls x* > bed.list

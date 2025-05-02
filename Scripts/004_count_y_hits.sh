#!/bin/bash
#SBATCH --job-name=Y_hits
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=Y_hits_%j.out

#----- START
echo "Starting job: $SLURM_JOB_NAME (Job ID: $SLURM_JOB_ID)"
echo "Running on node: $(hostname)"
echo "Start time: $(date)"
echo -e "#-----------------------------------------------------------#\n"

#----- Set directories
VCF_DIR="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/curiel/20250206_Mouse_Exomes/code/snakemake/freebayes"
OUTPUT_FILE="chrY_hits_summary.txt"

#----- Source and activate conda
source /optnfs/common/miniconda3/etc/profile.d/conda.sh
conda activate  /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/bcftools

#----- Script
# Print the header of the output file
echo -e "Sample\tHits" > "$OUTPUT_FILE"

# Iterate through each VCF file in the specified directory
for VCF_FILE in "$VCF_DIR"/*.merged.raw.vcf; do
    # Get the base name of the sample (file name without extension)
    SAMPLE_NAME=$(basename "$VCF_FILE" .merged.raw.vcf)
    
    # Count the number of hits on the Y chromosome in the VCF file
    # Assuming Y chromosome is named either 'Y' or 'chrY' in the VCF
    Y_HITS=$(bcftools view -H "$VCF_FILE" | grep -P "chrY" | wc -l)
    
    # Output results to file in the format "Sample Name: Y Hits"
    echo -e "$SAMPLE_NAME\t$Y_HITS" >> "$OUTPUT_FILE"
done

echo "Summary written to $OUTPUT_FILE"
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
conda activate  /dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/owen/sharedconda/miniconda/envs/bcftools

#----- Set directories
VCF_DIR="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/curiel/20250206_Mouse_Exomes/code/snakemake/filtered"
OUTPUT_DIR="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/curiel/20250206_Mouse_Exomes/code/snakemake/somatic_mutations"
TMP_DIR="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/curiel/20250206_Mouse_Exomes/code/snakemake/somatic_tmp"

#----- Make output directory
mkdir -p "$OUTPUT_DIR"
mkdir -p "$TMP_DIR"

#----- Make sure all VCFs have an index
for VCF_FILE in "$VCF_DIR"/*.vcf.gz; do
    bcftools index "$VCF_FILE"
done
echo "Indexing worked"

#----- Define VCF lists for Group 1 and Group 2
GROUP1_VCFS=("NCH1.filt.ann.vcf.gz" "NCH2.filt.ann.vcf.gz" "NCH3.filt.ann.vcf.gz" "NCN1.filt.ann.vcf.gz" "NCN3.filt.ann.vcf.gz" "NCN4.filt.ann.vcf.gz")  # Modify this list with actual VCF filenames
GROUP2_VCFS=("NFH1.filt.ann.vcf.gz" "NFH2.filt.ann.vcf.gz" "NFH3.filt.ann.vcf.gz" "NFN1.filt.ann.vcf.gz" "NFN2.filt.ann.vcf.gz" "NFN3.filt.ann.vcf.gz")  # Modify this list with actual VCF filenames

#----- Subtraction for Group 1 (using reference TNQ61R_CTRL)
REFERENCE1="$VCF_DIR/TNQ61R_CTRL.filt.ann.vcf.gz"
for VCF_FILE in "${GROUP1_VCFS[@]}"; do
    SAMPLE_NAME=$(basename "$VCF_FILE" .vcf.gz)
    echo "$SAMPLE_NAME"
    
    # Perform isec operation
    bcftools isec "$REFERENCE1" "$VCF_DIR/$VCF_FILE" -p "$TMP_DIR"
    
    # Move and rename files
    mv "$TMP_DIR/0001.vcf" "$OUTPUT_DIR/${SAMPLE_NAME}.somatic.filt.ann.vcf"
    
    # Remove intermediate files
    rm "$TMP_DIR/*.vcf"
    rm "$TMP_DIR*.vcf.gz"
    rm "$TMP_DIR/*.vcf.gz.tbi"
    rm "$TMP_DIR/README.txt"
    rm "$TMP_DIR/sites.txt"
done

#----- Subtraction for Group 2 (using reference for Group 2)
REFERENCE2="$VCF_DIR/TNQ61R_CD274_flfl.filt.ann.vcf.gz"
for VCF_FILE in "${GROUP2_VCFS[@]}"; do
    SAMPLE_NAME=$(basename "$VCF_FILE" .vcf.gz)
    echo "$SAMPLE_NAME"
    
    # Perform isec operation
    bcftools isec "$REFERENCE2" "$VCF_DIR/$VCF_FILE" -p "$TMP_DIR"
    
    # Move and rename files
    mv "$TMP_DIR/0001.vcf" "$OUTPUT_DIR/${SAMPLE_NAME}.somatic.filt.ann.vcf"
    
    # Remove intermediate files
    rm "$TMP_DIR/*.vcf"
    rm "$TMP_DIR*.vcf.gz"
    rm "$TMP_DIR/*.vcf.gz.tbi"
    rm "$TMP_DIR/README.txt"
    rm "$TMP_DIR/sites.txt"
done

#----- Double subtraction
bcftools isec "$VCF_DIR/NCH1.filt.ann.vcf.gz" "$VCF_DIR/NCH1P.filt.ann.vcf.gz" -p "$TMP_DIR"
mv "$TMP_DIR/0001.vcf" "$OUTPUT_DIR/NCH1P.somatic.filt.ann.vcf"
rm "$TMP_DIR/*.vcf"
rm "$TMP_DIR*.vcf.gz"
rm "$TMP_DIR/*.vcf.gz.tbi"
rm "$TMP_DIR/README.txt"
rm "$TMP_DIR/sites.txt"

#----- Triple subtraction
bgzip "$OUTPUT_DIR/NCH1P.somatic.filt.ann.vcf"
bcftools index "$OUTPUT_DIR/NCH1P.somatic.filt.ann.vcf.gz"
bcftools isec "$OUTPUT_DIR/NCH1P.somatic.filt.ann.vcf.gz" "$VCF_DIR/NCH1Lu.filt.ann.vcf.gz" -p "$TMP_DIR"
mv "$TMP_DIR/0001.vcf" "$OUTPUT_DIR/NCH1Lu.somatic.filt.ann.vcf"

#----- Unzip the zipped file
bgzip -d $OUTPUT_DIR/NCH1P.somatic.filt.ann.vcf.gz

#----- Clean up
rm "$OUTPUT_DIR/*.csi"
rm "$TMP_DIR/*

echo "Job complete. Summary written to $OUTPUT_DIR"
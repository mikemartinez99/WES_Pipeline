#! /usr/bin/bash
#SBATCH --job-name=WES_vcf2maf
#SBATCH --nodes=1
#SBATCH --partition=preempt1
#SBATCH --account=dac
#SBATCH --time=60:00:00
#SBATCH --mail-user=f007qps@dartmouth.edu
#SBATCH --mail-type=FAIL
#SBATCH --output=WES_vcf2maf_%j.out

# Define paths
vcf2maf="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/misc/shared-software/tools/vcf2maf/mskcc-vcf2maf-f6d0c40/vcf2maf.pl"
ref="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/genomic_references/mouse/UCSC_mm10/mm10.fa"
input_dir="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/curiel/20250206_Mouse_Exomes/code/snakemake/somatic_mutations"  # Change this to your actual VCF directory
output_dir="/dartfs-hpc/rc/lab/G/GMBSR_bioinfo/Labs/curiel/20250206_Mouse_Exomes/code/snakemake/somatic_MAFs"  # Change this to your actual MAF output directory

# Create the output directory if it doesn't exist
mkdir -p "$output_dir"

# Load module for tabix (required)
module load htslib

# Loop through all VCF files in the input directory
for vcf in "$input_dir"/*.vcf; do
    # Extract filename without extension
    sample_name=$(basename "$vcf" .vcf)
    
    # Define output MAF file path
    maf="$output_dir/${sample_name}.maf"
    
    # Run vcf2maf
    perl "$vcf2maf" --input-vcf "$vcf" \
        --output-maf "$maf" \
        --inhibit-vep \
        --ref-fasta "$ref" \
        --tumor-id "$sample_name"
    echo "Processed: $vcf -> $maf"
done


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Whole exome sequencing - somatic variant calling pipeline
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
#
# NOTES:
#   - the picard.yaml is a place holder conda environment for the freebayes and snpeff rules. Might want to run without conda
#   - double check that all freebayes variants calls inherently pass its quality filter if it has one.
#   - Additional filters in freebayes commands
#   - Check mitochondrial mapping methods, might be worth it to remove mito, X, and Y from results
#   - remove overlaps with repeat regions, which can also cause artifacts. you can get the repeatmasker file for mm10 from the uscs table browser. good to check on a test sample how many
#----- Import modules
import pandas as pd

#----- Set config file
configfile: "config.yaml"

#----- read in sample data
samples_df = pd.read_table(config["sample_tsv"]).set_index("sample_id", drop=False)
sample_list = list(samples_df['sample_id'])

#----- read in the bed list
bed_list = [line.strip() for line in open("bed.list")]

#----- Execute pipeline
rule all:
    input:
        expand("trimming/{sample}.R1.trim.fastq.gz", sample = sample_list),
        expand("trimming/{sample}.R2.trim.fastq.gz", sample = sample_list),
        expand("alignment/{sample}.srt.bam", sample = sample_list),
        expand("alignment/{sample}.srt.bam.bai", sample = sample_list),
        expand("markdup/{sample}.mkdup.bam", sample = sample_list),
        expand("markdup/{sample}.mkdup.log.txt", sample = sample_list),
        expand("markdup/{sample}.mkdup.bai", sample = sample_list),
        expand("metrics/{sample}.HSMetrics.txt", sample = sample_list),
        expand("freebayes/{sample}/{sample}.{region}.raw.vcf", sample = sample_list, region=bed_list),
        expand("freebayes/{sample}.merged.raw.vcf", sample = sample_list),
        expand("snpeff/{sample}.ann.vcf", sample = sample_list),
        "all_freebayes_commands.txt",
        expand("filtered/{sample}.filt.ann.vcf", sample = sample_list),
        expand("snpeff/{sample}.exclude_scaffolds.txt", sample = sample_list),
        expand("MAFs/{sample}.filt.ann.maf", sample = sample_list)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# RULES
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#----- Trimming
rule trim_reads:
    output:
        trimmed_fastq1 = "trimming/{sample}.R1.trim.fastq.gz",
        trimmed_fastq2 = "trimming/{sample}.R2.trim.fastq.gz",
        report = "trimming/{sample}.cutadapt.report"
    conda:
        "env_config/cutadapt.yaml",
    params:
        sample = lambda wildcards:  wildcards.sample,
        cutadapt = config["cutadapt_path"],
        fastq_file_1 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_1"],
        fastq_file_2 = lambda wildcards: samples_df.loc[wildcards.sample, "fastq_2"],
        adapterF = config["adapterF"],
        adapterR = config["adapterR"],
        nextseq_flag = "--nextseq-trim=20",
        min_length = "-m 18",
        max_n = "--max-n 0.8"
    resources: cpus="10", maxtime="2:00:00", mem_mb="60gb",
    shell: """
        cutadapt \
            -a {params.adapterF} -A {params.adapterR} \
            -o {output.trimmed_fastq1} \
            -p {output.trimmed_fastq2} \
            {params.fastq_file_1} {params.fastq_file_2} \
            {params.nextseq_flag} \
            {params.min_length} \
            {params.max_n} \
            --trim-n > {output.report}
    """

#----- Alignment
rule alignment:
    input: 
        trimmed_fastq_1 = "trimming/{sample}.R1.trim.fastq.gz",
        trimmed_fastq_2 = "trimming/{sample}.R2.trim.fastq.gz"
    output:
        sorted_bam = "alignment/{sample}.srt.bam",
        bai = "alignment/{sample}.srt.bam.bai"
    conda:
        "env_config/alignment.yaml",
    params:
        sample = lambda wildcards:  wildcards.sample,
        samtools = config["samtools_path"],
        reference = config["reference"],
        bwa = config["bwa_path"],
        rg = r"-R '@RG\tID:{sample}\tSM:{sample}'",
    resources: threads="6", maxtime="8:00:00", mem_mb="40gb",
    shell: """
        #----- Map with bwa mem and generated bam file
        {params.bwa} mem \
        {params.reference} \
        {input.trimmed_fastq_1} \
        {input.trimmed_fastq_2} \
        -t 8 \
        {params.rg} | {params.samtools} view -@ 2 -b | {params.samtools} sort -T \
        /scratch/samtools_{params.sample} \
        -@ 4 -m 512M 1> {output.sorted_bam} 

        #----- Sort sam file, convert to bam
        {params.samtools} index -@ 4 {output.sorted_bam}
    """

#----- Duplicate marking
rule mark_duplicates:
    input:
        bam_srt = "alignment/{sample}.srt.bam"
    output:
        dedup_bam = "markdup/{sample}.mkdup.bam",
        dedup_log = "markdup/{sample}.mkdup.log.txt",
        dedup_bai = "markdup/{sample}.mkdup.bai",
    conda:
        "env_config/picard.yaml"
    params:
        reference = config["reference"],
        picard_path = config["picard_path"],
        bedtools = config["bedtools_path"],
        targets = config["targets"]
    resources: cpus=2, maxtime="8:00:00", mem_mb=128000,
    shell: """
        #----- Mark duplicates
        picard -Xmx32g MarkDuplicates \
            I={input.bam_srt} \
            O={output.dedup_bam} \
            M={output.dedup_log} \
            OPTICAL_DUPLICATE_PIXEL_DISTANCE=100 \
            CREATE_INDEX=true \
            MAX_RECORDS_IN_RAM=2000000 \
            ASSUME_SORTED=true \
            MAX_FILE_HANDLES=768

            #----- Target coverages
            bedtools coverage \
            -a {params.targets} \
            -b {output.dedup_bam}
    """

#----- Metrics
rule collect_metrics:
    input: 
        mkdup_bam = "markdup/{sample}.mkdup.bam"
    output: 
        metrics = "metrics/{sample}.HSMetrics.txt"
    conda:
        "env_config/picard.yaml"
    params:
        reference = config["reference"],
        picard_path = config["picard_path"],
        targets = config["targets"],
        intervals = config["intervals"] # Generated with bed_to_interval_list.sh
    resources: threads="4", maxtime="8:00:00", mem_mb="12gb",
    shell: """
        #----- Collate metrics
        #picard CollectAlignmentSummaryMetrics \
        #    R={params.reference}.fa \
         #   I={input.mkdup_bam} \
         #   O={output.metrics}

        picard CollectHsMetrics \
            I={input.mkdup_bam} \
            O={output.metrics} \
            R={params.reference}.fa \
            BAIT_INTERVALS={params.intervals} \
            TARGET_INTERVALS={params.intervals}

    """

#----- Prepare FreeBayes commands
rule prepare_freebayes:
    input: 
    	bed_list = "bed.list",
        dedup_bams = expand("markdup/{sample}.mkdup.bam", sample=sample_list),
        dedup_bais = expand("markdup/{sample}.mkdup.bai", sample=sample_list)
    output: "all_freebayes_commands.txt"
    conda:
	    "freebayes"
    params:
        sample_list=" ".join(sample_list),  
        reference_fasta = config["reference_fa"],
    resources: cpus=2, maxtime="8:00:00", mem_mb="12gb",
    shell: """
        > {output}
        for sample in {params.sample_list}; do
            mkdir -p freebayes/$sample
            for i in `cat {input.bed_list}`; do 
                echo "freebayes -K -= --min-coverage 20 --min-alternate-fraction 1 -m 15 -f {params.reference_fasta} -t $i markdup/$sample.mkdup.bam > freebayes/$sample/$sample.$i.raw.vcf"
            done
        done >> {output}
    """

#----- Run FreeBayes commands in parallel
rule run_freebayes_parallel:
    input:
        commands = "all_freebayes_commands.txt"
    output:
        vcf_output = expand("freebayes/{sample}/{sample}.{region}.raw.vcf", sample = sample_list, region=bed_list)
    conda:
        "freebayes"
    resources: cpus=2, maxtime="8:00:00", mem_mb=128000,
    params:
        parallel_exec = config["parallel_path"],
        num_jobs = 16
    shell: """
        {params.parallel_exec} -j {params.num_jobs} < {input.commands}

    """

#----- Annotate mutations
rule snpeff:
    input:
        rawVCFs = expand("freebayes/{sample}/{sample}.{region}.raw.vcf", sample = sample_list, region=bed_list)
    output:
        mergedVCF = "freebayes/{sample}.merged.raw.vcf",
        annoVCF = "snpeff/{sample}.ann.vcf"
    conda:
        "bcftools"
    resources: cpus=2, maxtime="8:00:00", mem_mb=128000,
    params:
        sample = lambda wildcards: wildcards.sample,
        snpeff_jar = config["snpeff_path"],
        snpeff_genome = config["snpeff_genome"],
        java = config["java_path"],
    resources: threads="2", maxtime="4:00:00", memory="4gb"
    shell: """
        for sample_dir in freebayes; do
            # Merge per-region VCFs before annotation
            bcftools concat $(ls freebayes/{params.sample}/*.raw.vcf) -o {output.mergedVCF} -O v; done
        
        # Annotate the merged VCF using SnpEff
        {params.java} -jar {params.snpeff_jar} {params.snpeff_genome} {output.mergedVCF} > {output.annoVCF}
    """

#----- Filter mitochondrial and sex chromosomes, and unplaced scaffolds
rule filter_vcf:
    input:
        ann_vcf = "snpeff/{sample}.ann.vcf"
    output:
        filt_vcf = "filtered/{sample}.filt.ann.vcf",
        unplaced_scaffolds = "snpeff/{sample}.exclude_scaffolds.txt"
    conda:
        "bcftools"
    resources:
        cpus=2, maxtime="8:00:00", mem_mb=64000,
    params:
        sample = lambda wildcards: wildcards.sample,
        targets = config["targets"]
    shell: """
        # Extract the unplaced scaffolds (chromosome, start, end) from the BED file
        cut -f1,2,3 {params.targets} | grep -E "random|GL|chrUn|X|Y|M" | sort -u > {output.unplaced_scaffolds}
        
        # Remove variants called on unplaced scaffolds
        bcftools view -T ^{output.unplaced_scaffolds} {input.ann_vcf} -o {output.filt_vcf}
    """

#----- Convert VCFs to MAFs
rule generate_mafs:
    input:
        vcf_file = "filtered/{sample}.filt.ann.vcf"
    output:
        maf_files = "MAFs/{sample}.filt.ann.maf"
    conda:
        "bcftools"
    resources:
        cpus=2, maxtime="8:00:00", mem_mb=64000,
    params:
        sample = lambda wildcards: wildcards.sample,
        vcf2maf_path = config["vcf2maf_path"],
        reference_fasta = config["reference_fa"]
    shell: """
        # Run vcf2maf
        perl {params.vcf2maf_path} \
        --input-vcf {input.vcf_file} \
        --inhibit-vep \
        --ref-fasta {params.reference_fasta} \
        --tumor-id {params.sample}

    """
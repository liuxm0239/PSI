#########################################
# CNV calling using QDNAseq for lowpass WGS samples
# by Xingmin Liu
#########################################
from datetime import date
import json
from os.path import join, basename, dirname

configfile: 'config.yaml'

# Handle sample file names from samples.json
FILES = json.load(open(config['SAMPLES_JSON']))
SAMPLES = sorted(FILES.keys())

#SAMPLES={"MS22081258P1C"}
REFERENCE="/HD101TB/conda/Genome/hg38.fa"
REFBED="/HD101TB/conda/Genome/hg38.target.bed"
BWAINDEX="/HD101TB/conda/Genome/bwaIndex/hg38"
GATK="/SlurmDatabase/Bladder_cancer/DNA/rawdata/tools/gatk-4.3.0.0/gatk"
DBSNP="/HD101TB/conda/Genome/dbsnp_138.hg38.vcf.gz"
OMNI="/HD101TB/conda/Genome/1000G_omni2.5.hg38.vcf.gz"
PHASE1="/HD101TB/conda/Genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
MILLS="/HD101TB/conda/Genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
PON="/HD101TB/conda/Genome/cnvkit/PoN/Reference.cnn"
REFFLT="/HD101TB/conda/Genome/refFlat.txt"
QDNASEQ="/HD101TB/conda/Genome/QDNAseq/qdnaseq.R"
QDNASEG="/HD101TB/conda/Genome/QDNAseq/qdnaseq_segment_calls.R"
TIME=date.today().isoformat()
#TIME="2023-07-10"

rule all:
    input:
        expand("{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam", sample=SAMPLES, time=TIME ),
        expand("{time}/s03_variant/QDNAseq/{sample}/{sample}/100kb-bins/{sample}-copyNumbersCalled-100kb-bins.tab", sample=SAMPLES, time=TIME)

rule FASTP:
    input:
        R1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        R2 = lambda wildcards: FILES[wildcards.sample]['R2']
        #R1="raw.fastq/{sample}_1.fq.gz",
        #R2="raw.fastq/{sample}_2.fq.gz"
    output:
        "{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        "{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    threads: 2
    resources: mem_mb=16000
    conda:
        "Genomics"
    log:
        "{time}/s01_fastqc/{sample}.fastp_fastqc.log"
    shell: """
        #######################
        #trimmomatic PE -threads {threads} -trimlog {log} {input.R1} {input.R2} {wildcards.time}/s01_fastqc/{wildcards.sample}.trimmed_R1.fastq.gz {wildcards.time}/s01_fastqc/{wildcards.sample}.trimmed.unpaired_R1.fastq.gz {wildcards.time}/s01_fastqc/{wildcards.sample}.trimmed_R2.fastq.gz {wildcards.time}/s01_fastqc/{wildcards.sample}.trimmed.unpaired_R2.fastq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 CROP:200 HEADCROP:15 MINLEN:36 
	#echo "Working on raw data fastqc"
	#fastqc -t {threads} {input.R1} {input.R2} -o {wildcards.time}/s01_fastqc/ > {log} 2>&1
	
        echo "Working on fastp triming"
        fastp --in1 {input.R1} --in2 {input.R2} --out1 {output[0]} --out2 {output[1]} -5 -c -p -q 20 --json {wildcards.time}/s01_fastqc/{wildcards.sample}.fastp.Clean.json --html {TIME}/s01_fastqc/{wildcards.sample}.fastp.Clean.html --thread  {threads} --trim_poly_g --detect_adapter_for_pe >> {log} 2>&1 

        echo "Working on clean data fastqc"
        fastqc -t  {threads} {output} -o {wildcards.time}/s01_fastqc/ >> {log} 2>&1
	#######################
    """

rule BWA:
    input: 
        R1="{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        R2="{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    output:
        raw_bam=temp("{time}/s02_alignment/s020_bwa/{sample}.bam"),
        bam=temp("{time}/s02_alignment/s020_bwa/{sample}.sort.bam"),
        bai=temp("{time}/s02_alignment/s020_bwa/{sample}.sort.bam.bai")
    threads: 20
    resources: mem_mb=45000
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_alignment.log"
    params:
        platform ="Illumina"
    shell: """
        #######################	
        echo "Working on BWA alignment"
	bwa mem -M -Y -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}' {BWAINDEX} {input.R1} {input.R2} | samtools view -Sbh - -o {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}
        echo "Working on samtools sorting"
        samtools sort -m 3G -@ {threads} -T {wildcards.time}/s02_alignment/s020_bwa/tmp -o {output.bam} {TIME}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}

        echo "Working on bam flagstat"
	samtools index -@ {threads} {output.bam} 2>>{log} 
        samtools flagstat -@ {threads} {output.bam} > {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.sort.bam.stat 2>>{log}
	#######################
	"""

rule MarkDup:
    input:
        bam="{time}/s02_alignment/s020_bwa/{sample}.sort.bam",
        bai="{time}/s02_alignment/s020_bwa/{sample}.sort.bam.bai"
    output:
        bam=temp("{time}/s02_alignment/s021_Dedup/{sample}.dd.bam"),
        index=temp("{time}/s02_alignment/s021_Dedup/{sample}.dd.bam.bai"),
        metrics="{time}/s02_alignment/s021_Dedup/{sample}.markdup.txt"
    threads: 8
    resources: mem_mb=20000
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_s021_Dedup.log"
    shell: """
	#######################
        echo "Working on MarkDuplicates"
        {GATK} --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4  -Djava.io.tmpdir={wildcards.time}/s02_alignment/s021_Dedup/" MarkDuplicates -AS true -M {output.metrics} -I {input.bam} -O {output.bam} --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT  --CREATE_INDEX true >>{log} 2>&1
        mv {wildcards.time}/s02_alignment/s021_Dedup/{wildcards.sample}.dd.bai {output.index}
	#######################
	"""

rule GATK_BQSR:
    input:
        bam="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam",
        bai="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam.bai"
    output:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam",
        index="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam.bai"
    threads: 20
    resources: mem_mb=60000
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_s022_Brecal.log"
    shell: """
	#######################
        echo "Working on SplitIntervals for parallel BQSR"
        {GATK} SplitIntervals -R {REFERENCE} -L {REFBED} --scatter-count 20 -O {wildcards.time}/s02_alignment/s022_Brecal/interval-files

        echo "Working on GATK BQSR"
        for sam in {wildcards.sample}
        do 
            for i in `seq -f '%04g' 0 19`
            do
                outfile={wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_data_${{i}}.table
                {GATK} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
                BaseRecalibrator -L {wildcards.time}/s02_alignment/s022_Brecal/interval-files/${{i}}-scattered.interval_list \
                -R {REFERENCE} \
                --known-sites {DBSNP} --known-sites {OMNI} --known-sites {PHASE1} --known-sites {MILLS} \
                -I {input.bam} -O ${{outfile}} &
            done
            wait
            echo "Working on ApplyBQSR"
            for i in `seq -f '%04g' 0 19`
            do 
                bqfile={wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_data_${{i}}.table
                outputf={wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_${{i}}.bam
                {GATK} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
                ApplyBQSR -L {wildcards.time}/s02_alignment/s022_Brecal/interval-files/${{i}}-scattered.interval_list \
                -R {REFERENCE} -bqsr ${{bqfile}} \
                --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
                -I {input.bam} -O ${{outputf}} &
            done
            wait
        done
        
        echo "Working on merge BQSR bam"
        samtools merge -f {output.bam} {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal*.bam &&\
        samtools index -@ {threads} {output.bam}
        samtools flagstat {output.bam} > {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.BQSR.bam.stat 2>>{log}
        rm -f {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_*.ba*
	#######################
	"""

############
#   QDNAseq
############
rule QDNAseq:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam",
        bai="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam.bai"
    output:
        "{time}/s03_variant/QDNAseq/{sample}/{sample}/100kb-bins/{sample}-copyNumbersCalled-100kb-bins.tab"
    threads: 2
    resources: mem_mb=4000
    conda:
        "Genomics"
    log:
        "{time}/s03_variant/{sample}_QDNAseq.log"
    shell: """
	#######################
        #!/bin/bash
        echo "Working on QDNAseq"
        #mkdir -p {wildcards.time}/s03_variant/QDNAseq/{wildcards.sample}
        ln -s $(realpath {input.bam}) {wildcards.time}/s03_variant/qDNAseq/{wildcards.sample}
        ln -s $(realpath {input.bai}) {wildcards.time}/s03_variant/qDNAseq/{wildcards.sample}

        /HD101TB/conda/miniconda3/envs/Genomics/bin/Rscript {QDNASEQ} {wildcards.time}/s03_variant/qDNAseq/{wildcards.sample} {wildcards.time}/s03_variant/qDNAseq/{wildcards.sample} 100 {wildcards.sample}
        /HD101TB/conda/miniconda3/envs/Genomics/bin/Rscript {QDNASEG} {wildcards.time}/s03_variant/qDNAseq/{wildcards.sample}/{wildcards.sample}-binSize-100/{wildcards.sample}-copyNumbers-100kb-bins.rds {wildcards.time}/s03_variant/QDNAseq/{wildcards.sample} {wildcards.sample}

	#######################
	"""
############

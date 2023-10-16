#########################################
# CNV calling using cnvkit for WGS samples
# by Xingmin Liu
#########################################
from datetime import date

SAMPLES={"Ajing", "SHP",  "wangfeiyunB",  "wanghuajun"}
REFERENCE="/root/conda/Genome/hg38.fa"
REFBED="/root/conda/Genome/hg38.target.bed"
BWAINDEX="/root/conda/Genome/bwaIndex/hg38"
GATK="/SlurmDatabase/Bladder_cancer/DNA/rawdata/tools/gatk-4.3.0.0/gatk"
DBSNP="/root/conda/Genome/dbsnp_138.hg38.vcf.gz"
OMNI="/root/conda/Genome/1000G_omni2.5.hg38.vcf.gz"
PHASE1="/root/conda/Genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
MILLS="/root/conda/Genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
PON="/root/conda/Genome/cnvkit/PoN/Reference.cnn"
REFFLT="/root/conda/Genome/refFlat.txt"
TIME=date.today().isoformat()
#TIME="2023-06-30"

rule all:
    input:
        expand("{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam", sample=SAMPLES, time=TIME ),
        expand("{time}/s03_somatic_variant/cnvkit/{sample}.marked.BQSR.call.cns", sample=SAMPLES, time=TIME)

rule FASTP:
    input:
        R1="raw.fastq/{sample}_raw_1.fq.gz",
        R2="raw.fastq/{sample}_raw_2.fq.gz"
    output:
        "{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        "{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    resources: mem_mb=6000, cpus=2
    conda:
        "Genomics"
    log:
        "{time}/s01_fastqc/{sample}.fastp_fastqc.log"
    shell: """
        #######################
	#echo "Working on raw data fastqc"
	#fastqc -t {resources.cpus} {input.R1} {input.R2} -o {wildcards.time}/s01_fastqc/ > {log} 2>&1
	
        echo "Working on fastp triming"
        fastp --in1 {input.R1} --in2 {input.R2} --out1 {output[0]} --out2 {output[1]} -5 -c -p -q 20 --json {wildcards.time}/s01_fastqc/{wildcards.sample}.fastp.Clean.json --html {TIME}/s01_fastqc/{wildcards.sample}.fastp.Clean.html --thread  {resources.cpus} --trim_poly_g --detect_adapter_for_pe >> {log} 2>&1 

        echo "Working on clean data fastqc"
        fastqc -t  {resources.cpus} {output} -o {wildcards.time}/s01_fastqc/ >> {log} 2>&1
	#######################
    """

rule BWA:
    input: 
        R1="{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        R2="{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    output:
        raw_bam="{time}/s02_alignment/s020_bwa/{sample}.bam",
        bam="{time}/s02_alignment/s020_bwa/{sample}.sort.bam",
        bai="{time}/s02_alignment/s020_bwa/{sample}.sort.bam.bai"
    resources: mem_mb=45000, cpus=20
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_alignment.log"
    params:
        platform ="Illumina"
    shell: """
        #######################	
        echo "Working on BWA alignment"
	bwa mem -M -Y -t {resources.cpus} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}' {BWAINDEX} {input.R1} {input.R2} | samtools view -Sbh - -o {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}
        echo "Working on samtools sorting"
        samtools sort -m 3G -@ {resources.cpus} -T {wildcards.time}/s02_alignment/s020_bwa/tmp -o {output.bam} {TIME}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}

        echo "Working on bam flagstat"
	samtools index -@ {resources.cpus} {output.bam} 2>>{log} 
        samtools flagstat -@ {resources.cpus} {output.bam} > {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.sort.bam.stat 2>>{log}
	#######################
	"""

rule MarkDup:
    input:
        bam="{time}/s02_alignment/s020_bwa/{sample}.sort.bam",
        bai="{time}/s02_alignment/s020_bwa/{sample}.sort.bam.bai"
    output:
        bam="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam",
        index="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam.bai",
        metrics="{time}/s02_alignment/s021_Dedup/{sample}.markdup.txt"
    resources: mem_mb=20000, cpus=8
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_s021_Dedup.log"
    shell: """
	#######################
        echo "Working on MarkDuplicates"
        {GATK} --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4  -Djava.io.tmpdir={wildcards.time}/s02_alignment/s021_Dedup/" MarkDuplicates -AS true -M {output.metrics} -I {input.bam} -O {output.bam} --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT  --CREATE_INDEX true >>{log} 2>&1
        mv {wildcards.time}/s02_alignment/s021_Dedup/{wildcards.sample}.dd.bai {output.index}
        samtools flagstat -@ {resources.cpus} {output.bam} > {wildcards.time}/s02_alignment/s021_Dedup/{wildcards.sample}.dedup.bam.stat
	#######################
	"""

rule GATK_BQSR:
    input:
        bam="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam",
        bai="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam.bai"
    output:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam",
        index="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam.bai"
    resources: mem_mb=60000, cpus=20
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
        samtools index -@ {resources.cpus} {output.bam}
        samtools flagstat -@ {resources.cpus} {output.bam} > {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.BQSR.bam.stat 2>>{log}
        rm -f {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_*.ba*
	#######################
	"""

############
#       cnvkit
############
rule cnvkitPooled:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam"
    output:
        "{time}/s03_somatic_variant/cnvkit/{sample}.marked.BQSR.call.cns"
    resources: mem_mb=4000, cpus=2
    conda:
        "Genomics"
    log:
        "{time}/s03_somatic_variant/{sample}_cnvkit.log"
    shell: """
	#######################
        #!/bin/bash
        echo "Working on CNVkit"
        #source activate Genomics
	cnvkit.py batch {input.bam} -r {PON} -p {resources.cpus} --method wgs --drop-low-coverage --scatter --diagram -d {wildcards.time}/s03_somatic_variant/cnvkit/ >{log} 2>&1
	#cnvkit.py scatter -s {wildcards.time}/s03_somatic_variant/cnvkit/{wildcards.sample}.cn{{s,r}} -o s03_somatic_variant/cnvkit/{wildcards.sample}.pdf
	#######################
	"""
############

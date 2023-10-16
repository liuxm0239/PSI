#########################################
# CNV calling using cnvkit for WGS samples
# by Xingmin Liu
#########################################
from datetime import date
import json
from os.path import join, basename, dirname

configfile: 'config.yaml'

# Handle sample file names from samples.json
FILES = json.load(open(config['SAMPLES_JSON']))
SAMPLES = sorted(FILES.keys())

PLATFORM="MGISEQ"
REFERENCE="/SlurmDatabase/Clinical/ClinSV/reference/mm10/mm10.fasta"
REFBED="/SlurmDatabase/Clinical/ClinSV/reference/mm10/mm10.bed"
#BWAINDEX="/HD101TB/conda/Genome/bwaIndex/hg38"
GATK="/SlurmDatabase/Bladder_cancer/DNA/rawdata/tools/gatk-4.3.0.0/gatk"
#HAPMAP="/HD101TB/conda/Genome/hapmap_3.3.hg38.vcf.gz"
DBSNP="/SlurmDatabase/Clinical/ClinSV/reference/mm10/mgp.v5.snp.pass.sorted.vcf.gz"
DBINDEL="/SlurmDatabase/Clinical/ClinSV/reference/mm10/mgp.v5.indels.pass.sorted.vcf.gz"
#OMNI="/HD101TB/conda/Genome/1000G_omni2.5.hg38.vcf.gz"
#PHASE1="/HD101TB/conda/Genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
#MILLS="/HD101TB/conda/Genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
PON="/HD101TB/conda/Genome/cnvkit/PoN/mouse_reference.cnn"
TIME=date.today().isoformat()
#TIME="2023-07-19"

rule all:
    input:
        expand("{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam", sample=SAMPLES, time=TIME ),
        expand("{time}/s03_variant/cnvkit/{sample}.marked.BQSR.call.cns", sample=SAMPLES, time=TIME),
        expand("{time}/s03_variant/GATK/{sample}_HC_calls.vcf.gz" , sample=SAMPLES, time=TIME),
        expand("{time}/s03_variant/GATK/{sample}_HC.CNNscore.filtered.vcf", sample=SAMPLES, time=TIME)

rule FASTP:
    input:
        R1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        R2 = lambda wildcards: FILES[wildcards.sample]['R2']
        #R1="raw.fastq/{sample}_raw_1.fq.gz",
        #R2="raw.fastq/{sample}_raw_2.fq.gz"
    output:
        "{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        "{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    threads: 2
    resources: mem_mb=6000
    #conda:        "Genomics"
    log:
        "{time}/s01_fastqc/{sample}.fastp_fastqc.log"
    shell: """
        #######################
	#echo "Working on raw data fastqc"
	#fastqc -t {threads} {input.R1} {input.R2} -o {wildcards.time}/s01_fastqc/ > {log} 2>&1
	
        echo "Working on fastp triming"
        fastp --in1 {input.R1} --in2 {input.R2} \
            --out1 {output[0]} --out2 {output[1]} \
            -5 -c -p -q 20 \
            --json {wildcards.time}/s01_fastqc/{wildcards.sample}.fastp.Clean.json \
            --html {TIME}/s01_fastqc/{wildcards.sample}.fastp.Clean.html \
            --thread  {threads} \
            --trim_poly_g --detect_adapter_for_pe >> {log} 2>&1 

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
    resources: mem_mb=16000
    conda:        "base"
    log:
        "{time}/s02_alignment/{sample}_alignment.log"
    params:
        platform ="Illumina"
    shell: """
        #######################	
        echo "Working on BWA alignment"
	/usr/local/bioinfo/bin/bwa mem -M -Y -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{PLATFORM}' {REFERENCE} {input.R1} {input.R2} | /usr/local/bioinfo/bin/samtools view -Sbh - -o {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}
        echo "Working on samtools sorting"
        /usr/local/bioinfo/bin/samtools sort -m 3G -@ {threads} -T {wildcards.time}/s02_alignment/s020_bwa/tmp_{wildcards.sample} -o {output.bam} {TIME}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}

        echo "Working on bam flagstat"
	/usr/local/bioinfo/bin/samtools index -@ {threads} {output.bam} 2>>{log} 
        /usr/local/bioinfo/bin/samtools flagstat -@ {threads} {output.bam} > {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.sort.bam.stat 2>>{log}
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
    conda:        "base"
    log:
        "{time}/s02_alignment/{sample}_s021_Dedup.log"
    shell: """
	#######################
        echo "Working on MarkDuplicates"
        {GATK} --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4  -Djava.io.tmpdir={wildcards.time}/s02_alignment/s021_Dedup/tmp_{wildcards.sample}" MarkDuplicates -AS true -M {output.metrics} -I {input.bam} -O {output.bam} --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT  --CREATE_INDEX true >>{log} 2>&1
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
    conda:        "base"
    log:
        "{time}/s02_alignment/{sample}_s022_Brecal.log"
    shell: """
	#######################
        echo "Working on SplitIntervals for parallel BQSR"
        {GATK} SplitIntervals -R {REFERENCE} -L {REFBED} --scatter-count 20 -O {wildcards.time}/s02_alignment/s022_Brecal/interval-files 2>&1 > {log}

        echo "Working on GATK BQSR"
        for sam in {wildcards.sample}
        do 
            for i in `seq -f '%04g' 0 19`
            do
                outfile={wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_data_${{i}}.table
                {GATK} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
                BaseRecalibrator -L {wildcards.time}/s02_alignment/s022_Brecal/interval-files/${{i}}-scattered.interval_list \
                -R {REFERENCE} \
                --known-sites {DBSNP} --known-sites {DBINDEL} \
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
        /usr/local/bioinfo/bin/samtools merge -f {output.bam} {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal*.bam &&\
        /usr/local/bioinfo/bin/samtools index -@ {threads} {output.bam}
        /usr/local/bioinfo/bin/samtools flagstat {output.bam} > {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.BQSR.bam.stat 2>>{log}
        rm -f {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_*.ba*
        rm -f {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}*.table
	#######################
	"""

############  
# HaplotypeCaller gVCF
############ 
rule gVCF:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam"
    output:
        "{time}/s03_variant/GATK/{sample}.g.vcf.gz"
    threads: 20
    resources:
        mem_mb=30000 
    #conda:         "Genomics"
    log:
        "{time}/s03_variant/{sample}_GATK.log"
    shell:"""
        #######################
        #!/bin/bash
        echo "Working on HaplotypeCaller gVCF"
        {GATK} SplitIntervals -R {REFERENCE} -L {REFBED} --scatter-count 10\
            -O {wildcards.time}/s03_variant/interval-files

        for sam in {wildcards.sample}
        do 
            for i in `seq -f '%04g' 0 9`
            do 
                outfile={wildcards.time}/s03_variant/GATK/{wildcards.sample}_${{i}}.g.vcf.gz
                {GATK} --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" HaplotypeCaller \
                    -R {REFERENCE} -I {input.bam} -O ${{outfile}} \
                    -L {wildcards.time}/s03_variant/interval-files/${{i}}-scattered.interval_list \
                    -ERC GVCF &
            done
            wait
        done

        find {wildcards.time}/s03_variant/GATK/ -type f -name "{wildcards.sample}_*.g.vcf.gz" > {wildcards.time}/s03_variant/GATK/{wildcards.sample}_input.list
        {GATK} --java-options '-Xmx60g' CombineGVCFs \
                -R {REFERENCE} -V {wildcards.time}/s03_variant/GATK/{wildcards.sample}_input.list \
                --create-output-variant-index -O {output[0]}
        #######################
        """
############

rule GenotypeGVCFs:
    input:
        "{time}/s03_variant/GATK/{sample}.g.vcf.gz"
    output:
        "{time}/s03_variant/GATK/{sample}_HC_calls.vcf.gz"
    #conda:        "Genomics"
    log:
        '{time}/s03_variant/GATK/{sample}_GenotypeGVCFs.log'
    threads: 4
    resources: mem_mb=8000
    shell:"""
    #!/bin/bash
    echo "Working on HaplotypeCaller GenotypeGVCFs"
    {GATK} --java-options '-Xms2G -Xmx8G -XX:ParallelGCThreads=4' GenotypeGVCFs \
        -R {REFERENCE}  -V {input[0]}  -O {output[0]}
    """

rule CNNScoreVariants:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam",
        vcf="{time}/s03_variant/GATK/{sample}_HC_calls.vcf.gz"
    output:
        anno="{time}/s03_variant/GATK/{sample}_HC.CNNscore.vcf.gz",
        filt="{time}/s03_variant/GATK/{sample}_HC.CNNscore.filtered.vcf.gz"
    conda:
        "gatk"
    log:
        '{time}/s03_variant/GATK/{sample}_CNNScoreVariants.log'
    threads: 15
    resources: mem_mb=50000
    shell:"""
    #!/bin/bash
    echo "Working on CNNScoreVariants annotation"
    for sam in {wildcards.sample}
    do
        for i in `seq -f '%04g' 0 11`
        do
            outfile={wildcards.time}/s03_variant/GATK/{wildcards.sample}_${{i}}.CNNscore.vcf.gz
            {GATK} --java-options '-Xms2G -Xmx4G -XX:ParallelGCThreads=2'  CNNScoreVariants \
                    -I {input.bam} -V {input.vcf} -R {REFERENCE} \
                    -L {wildcards.time}/s03_variant/interval-files/${{i}}-scattered.interval_list \
                    -O ${{outfile}} \
                    -tensor-type read_tensor &
        done
        wait
    done
    find {wildcards.time}/s03_variant/GATK/ -type f -name "{wildcards.sample}_*.CNNscore.vcf.gz" >\
            {wildcards.time}/s03_variant/GATK/{wildcards.sample}_CNNscore.list
    {GATK} --java-options '-Xmx60g' MergeVcfs \
            -I {wildcards.time}/s03_variant/GATK/{wildcards.sample}_CNNscore.list \
            -O {output.anno}

    echo "Working on FilterVariantTranches"
    {GATK} --java-options '-Xms2G -Xmx4G -XX:ParallelGCThreads=2' FilterVariantTranches \
        -V {output.anno} -O {output.filt} \
        -resource {DBSNP}  -resource {DBINDEL}
    """

############
#       cnvkit
############
rule cnvkitPooled:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam"
    output:
        "{time}/s03_variant/cnvkit/{sample}.marked.BQSR.call.cns"
    threads: 2
    resources: mem_mb=4000
    conda:
        "Genomics"
    log:
        "{time}/s03_variant/{sample}_cnvkit.log"
    shell: """
	#######################
        #!/bin/bash
        echo "Working on CNVkit"

	cnvkit.py batch {input.bam} -r {PON} -p {threads} --method wgs \
                --drop-low-coverage --scatter --diagram -d {wildcards.time}/s03_variant/cnvkit/ >{log} 2>&1
	#cnvkit.py scatter -s {wildcards.time}/s03_variant/cnvkit/{wildcards.sample}.cn{{s,r}} -o s03_variant/cnvkit/{wildcards.sample}.pdf
	#######################
	"""
############

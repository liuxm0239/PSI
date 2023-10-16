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

REFERENCE="/SlurmDatabase/Clinical/ClinSV/reference/hg38/hg38.fasta"
REFBED="/SlurmDatabase/Clinical/ClinSV/reference/hg38/hg38.bed"
GATK="/opt/gatk/gatk-4.3.0.0/gatk"
plotBAF="/HD101TB/bioinfo/liuxingmin/Clinical/snk/bin/plotBAF.R"
HAPMAP="/HD101TB/conda/Genome/hapmap_3.3.hg38.vcf.gz"
DBSNP="/HD101TB/conda/Genome/dbsnp_138.hg38.vcf.gz"
OMNI="/HD101TB/conda/Genome/1000G_omni2.5.hg38.vcf.gz"
PHASE1="/HD101TB/conda/Genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
MILLS="/HD101TB/conda/Genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
PON="/HD101TB/conda/Genome/cnvkit/PoN/Reference.cnn"
#TIME=date.today().isoformat()
TIME=config['TIME']

rule all:
    input:
        expand("{time}/s01_fastqc/{sample}_clean_R2_fastqc.html", sample=SAMPLES, time=TIME),
        expand("{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam", sample=SAMPLES, time=TIME ),
        #expand("{time}/s03_variant/cnvkit/{sample}.marked.BQSR.call.cns", sample=SAMPLES, time=TIME),
        expand("{time}/s03_variant/GATK/{sample}_merged_VQSR_SNP_INDEL.vcf.gz", sample=SAMPLES, time=TIME)
        #expand("{time}/s03_variant/GATK/{sample}_HC.CNNscore.filtered.vcf.gz", sample=SAMPLES, time=TIME),
        #expand("{time}/s03_variant/clinSV/{sample}/SVs/joined/SV-CNV.PASS.vcf" , sample=SAMPLES, time=TIME),
        #expand("{time}/s03_variant/clinSV/{sample}/SVs/joined/CNV.RARE_PASS_GENE.txt", sample=SAMPLES, time=TIME)

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
    resources: mem_mb=6000
    #conda:        "Genomics"
    log:
        "{time}/s01_fastqc/{sample}.fastp.log"
    shell: """
        #######################
	#echo "Working on raw data fastqc"
	#fastqc -t {threads} {input.R1} {input.R2} -o {wildcards.time}/s01_fastqc/ 2>&1 > {log} 
	
        echo "Working on fastp triming"
        # !!! note --fix_mgi_id cause memory leak, do not use
        fastp --in1 {input.R1} --in2 {input.R2} \
                --adapter_sequence CTGTCTCTTATACACATCT \
                --adapter_sequence_r2 AGATGTGTATAAGAGACAG \
                --out1 {output[0]} --out2 {output[1]} \
                --detect_adapter_for_pe --trim_poly_g --trim_poly_x \
                --length_required 50 \
                --json {wildcards.time}/s01_fastqc/{wildcards.sample}.fastp.json \
                --html {TIME}/s01_fastqc/{wildcards.sample}.fastp.html \
                --thread {threads} -p 2>&1  >> {log}

        #echo "Working on clean data fastqc"
        #fastqc -t  {threads} {output} -o {wildcards.time}/s01_fastqc/ 2>&1 >> {log}
	#######################
    """

rule fastQC:
    input:
        R1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        R2 = lambda wildcards: FILES[wildcards.sample]['R2'],
        R1_clean =  "{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        R2_clean =  "{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
        #R1="raw.fastq/{sample}_1.fq.gz",
        #R2="raw.fastq/{sample}_2.fq.gz"
    output:
        "{time}/s01_fastqc/{sample}_clean_R1_fastqc.html",
        "{time}/s01_fastqc/{sample}_clean_R2_fastqc.html"
    threads: 4
    resources: mem_mb=6000
    #conda:        "Genomics"
    log:
        "{time}/s01_fastqc/{sample}.fastqc.log"
    shell: """
        #######################
        echo "Working on clean data fastqc"
        fastqc -t {threads} {input} -o {wildcards.time}/s01_fastqc/ 2>&1 >>{log}
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
    resources: mem_mb=20000
    #conda:        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_alignment.log"
    params:
        platform = config['PLATFORM']
    shell: """
        #######################	
        echo "Working on BWA alignment"
	bwa mem -M -Y -t {threads} \
                -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}' \
                {REFERENCE} {input.R1} {input.R2} | samtools view -Sbh -F 256 - \
                -o {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}
        echo "Working on samtools sorting"
        samtools sort -m 1G -@ {threads} \
                -T {wildcards.time}/s02_alignment/s020_bwa/tmp_{wildcards.sample} \
                -o {output.bam} \
                {TIME}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>>{log}

        echo "Working on bam flagstat"
	samtools index -@ {threads} {output.bam} 2>>{log} 
        samtools flagstat -@ {threads} {output.bam} \
                > {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.sort.bam.stat 2>>{log}
        {GATK} --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
                ValidateSamFile \
                -I {output.bam} -MODE SUMMARY \
                -O {wildcards.time}/s02_alignment/s020_bwa/{wildcards.sample}.sort.bam.Val.txt
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
    #conda:        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_s021_Dedup.log"
    shell: """
	#######################
        echo "Working on MarkDuplicates"
        {GATK} --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4 -Djava.io.tmpdir={wildcards.time}/s02_alignment/s021_Dedup/tmp_{wildcards.sample}" \
                MarkDuplicates -AS true -M {output.metrics} -I {input.bam} -O {output.bam} \
                --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT  \
                --CREATE_INDEX true 2>&1 >{log}

        mv {wildcards.time}/s02_alignment/s021_Dedup/{wildcards.sample}.dd.bai {output.index}
        {GATK} --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
                ValidateSamFile -I {output.bam} -MODE SUMMARY \
                -O {wildcards.time}/s02_alignment/s021_Dedup/{wildcards.sample}.dd.bam.Val.txt
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
    #conda:        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_s022_Brecal.log"
    shell: """
	#######################
        echo "Working on SplitIntervals for parallel BQSR"
        if [ ! -d "{wildcards.time}/s02_alignment/s022_Brecal/interval-files" ] 
        then
            {GATK} SplitIntervals -R {REFERENCE} -L {REFBED} --scatter-count 20 \
                    -O {wildcards.time}/s02_alignment/s022_Brecal/interval-files 2>&1 > {log}
        fi

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
                -I {input.bam} -O ${{outfile}} 2>>{log}  & 
            done
            wait

            find {wildcards.time}/s02_alignment/s022_Brecal/ \
                    -type f -name "{wildcards.sample}*_*.table" \
                    > {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.recal.list 2>>{log}

            {GATK} GatherBQSRReports -I {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.recal.list \
                    -O {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.recal.table 2>>{log} 
        done
        
        #samtools view -b -f 4 {input.bam} > {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.unmapped.bam

        bqfile={wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.recal.table
        {GATK} --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
            ApplyBQSR -R {REFERENCE} -bqsr ${{bqfile}} \
            --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
            -I {input.bam} -O {output.bam} 2>>{log} 
        mv {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.marked.BQSR.bai {output.index}
        #samtools index -@ {threads} {output.bam} 2>>{log}

        # bam validation
        {GATK} --java-options "-Xmx20G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
                ValidateSamFile -I {output.bam} -MODE SUMMARY \
                -O {wildcards.sample}.marked.BQSR.bam.Val.txt

        echo "Working on CollectInsertSizeMetrics"
        {GATK} CollectInsertSizeMetrics -I {output.bam} \
                -O {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_insert_size_metrics.txt \
                -H {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_insert_size_metrics.pdf -M 0.01 2>>{log}

        samtools flagstat {output.bam} > {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.BQSR.bam.stat 2>>{log}
        rm -f {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_dedup_recal_*.ba*
        rm -f {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}*.table

        echo "Working on Mosdepth"
        mkdir -p  {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_mosdepth
        mosdepth -t {threads} -b 5000 -Q 20 \
                {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}_mosdepth/{wildcards.sample}_depth {output.bam} 2>>{log}

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
        mem_mb=36000 
    #conda:         "Genomics"
    log:
        "{time}/s03_variant/GATK/{sample}_gVCF.log"
    shell:"""
        #######################
        #!/bin/bash
        echo "Working on HaplotypeCaller gVCF"
        if [ ! -d "{wildcards.time}/s03_variant/interval-files" ]
        then 
            {GATK} SplitIntervals -R {REFERENCE} -L {REFBED} --scatter-count 20\
                    -O {wildcards.time}/s03_variant/interval-files 2>>{log}
        fi

        for sam in {wildcards.sample}
        do 
            for i in `seq -f '%04g' 0 19`
            # '0 11' according to {GATK} SplitIntervals  --scatter-count 20
            do 
                outfile={wildcards.time}/s03_variant/GATK/{wildcards.sample}_${{i}}.g.vcf.gz
                {GATK} --java-options "-Xms4G -Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" HaplotypeCaller \
                    -R {REFERENCE} -I {input.bam} -O ${{outfile}} \
                    -L {wildcards.time}/s03_variant/interval-files/${{i}}-scattered.interval_list \
                    --native-pair-hmm-threads 6 \
                    -ERC GVCF &
            done
            wait
        done

        find {wildcards.time}/s03_variant/GATK/ -type f -name "{wildcards.sample}_*.g.vcf.gz" \
                > {wildcards.time}/s03_variant/GATK/{wildcards.sample}_input.list
        find {wildcards.time}/s03_variant/GATK/ -type f -name "{wildcards.sample}_*.g.vcf.gz.tbi" \
                > {wildcards.time}/s03_variant/GATK/{wildcards.sample}_input.tbi.list

        {GATK} --java-options '-Xmx10g -XX:+UseParallelGC -XX:ParallelGCThreads=2' CombineGVCFs \
                -R {REFERENCE} -V {wildcards.time}/s03_variant/GATK/{wildcards.sample}_input.list \
                --create-output-variant-index -O {output[0]} 2>>{log}

        xargs rm < {wildcards.time}/s03_variant/GATK/{wildcards.sample}_input.list
        xargs rm < {wildcards.time}/s03_variant/GATK/{wildcards.sample}_input.tbi.list
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
    threads: 8
    resources: mem_mb=8000
    shell:"""
    #!/bin/bash
    echo "Working on HaplotypeCaller GenotypeGVCFs"
    {GATK} --java-options '-Xms2G -Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4' GenotypeGVCFs \
        -R {REFERENCE}  -V {input[0]}  -O {output[0]}  2>>{log}
    """

rule VariantRecalibrator:
    input:
        "{time}/s03_variant/GATK/{sample}_HC_calls.vcf.gz"
    output:
        snprecal=temp("{time}/s03_variant/GATK/{sample}_SNP.recal"),
        snptranches=temp("{time}/s03_variant/GATK/{sample}_SNP.tranches"),              
        snpplot=temp("{time}/s03_variant/GATK/{sample}_SNP.plots.R"),
        indelrecal=temp("{time}/s03_variant/GATK/{sample}_INDEL.recal"),
        indeltranches=temp("{time}/s03_variant/GATK/{sample}_INDEL.tranches"),
        indelplot=temp("{time}/s03_variant/GATK/{sample}_INDEL.plots.R"),
        snpvqsr=temp("{time}/s03_variant/GATK/{sample}_merged_VQSR_SNP.vcf.gz"),
        indelvqsr="{time}/s03_variant/GATK/{sample}_merged_VQSR_SNP_INDEL.vcf.gz",
        vqsrdown=temp("{time}/s03_variant/GATK/{sample}_merged_VQSR_SNP_INDEL.downsample.vcf.gz")
    #conda:        "Genomics"
    log:
        '{time}/s03_variant/GATK/{sample}_VQSR.log'
    threads: 8
    resources: mem_mb=16000
    shell:"""
    #!/bin/bash   
    echo "Working on VQSR"

    {GATK} --java-options "-Xms4G -Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2" VariantRecalibrator \
        -tranche 100.0 -tranche 99.95 -tranche 99.9 \
        -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
        -tranche 95.0 -tranche 94.0 \
        -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
        -R {REFERENCE} -V {input[0]} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {HAPMAP} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 {OMNI} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {PHASE1} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP -O {output.snprecal} --tranches-file {output.snptranches} \
        --rscript-file {output.snpplot} 2>&1 >>{log}

    {GATK} --java-options "-Xms4G -Xmx4G -XX:+UseParallelGC  -XX:ParallelGCThreads=2" VariantRecalibrator \
        -tranche 100.0 -tranche 99.95 -tranche 99.9 \
        -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
        -tranche 95.0 -tranche 94.0 \
        -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
        -R {REFERENCE} -V {input[0]} \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 {MILLS} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {DBSNP} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
        -mode INDEL -O {output.indelrecal} --tranches-file {output.indeltranches} \
        --rscript-file {output.indelplot} 2>&1 >>{log}

    echo "Working on ApplyVQSR"
    {GATK} --java-options "-Xms2G -Xmx2G -XX:+UseParallelGC -XX:ParallelGCThreads=2" ApplyVQSR \
        -V {input[0]} \
        -mode SNP \
        --recal-file {output.snprecal} \
        --tranches-file {output.indeltranches} \
        --truth-sensitivity-filter-level 99.9 \
        --create-output-variant-index true \
        -O {output.snpvqsr} 2>&1 >>{log}

    {GATK} --java-options "-Xms2G -Xmx2G -XX:+UseParallelGC -XX:ParallelGCThreads=2" ApplyVQSR \
        -V {output.snpvqsr} \
        -mode INDEL \
        --recal-file {output.indelrecal} \
        --tranches-file {output.indeltranches} \
        --truth-sensitivity-filter-level 99.9 \
        --create-output-variant-index true \
        -O {output.indelvqsr} 2>&1 >>{log}
    
    echo "Working on plotBAF"
    sleep 120

    (zcat {output.indelvqsr} | head -n 10000 | zgrep ^# ; zgrep -v ^# {output.indelvqsr} | \
            shuf -n 150000 | LC_ALL=C sort -k1,1V -k2,2n) | bgzip > {output.vqsrdown}

    tabix -p vcf {output.vqsrdown}
    Rscript {plotBAF} {output.vqsrdown} {wildcards.sample}
    """

rule CNNScoreVariants:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam",
        vcf="{time}/s03_variant/GATK/{sample}_HC_calls.vcf.gz"
    output:
        anno="{time}/s03_variant/GATK/{sample}_HC.CNNscore.vcf.gz",
        filt="{time}/s03_variant/GATK/{sample}_HC.CNNscore.filtered.vcf.gz"
    conda: "gatk"
    log:
        '{time}/s03_variant/GATK/{sample}_CNNScoreVariants.log'
    threads: 15
    resources: mem_mb=48000
    shell:"""
    #!/bin/bash

    echo "Working on CNNScoreVariants annotation"
        if [ ! -d "{wildcards.time}/s03_variant/GATK/interval-files" ]
        then 
            {GATK} SplitIntervals -R {REFERENCE} -L {REFBED} --scatter-count 12\
                    -O {wildcards.time}/s03_variant/GATK/interval-files 2>>{log}
        fi
    for sam in {wildcards.sample}
    do
        for i in `seq -f '%04g' 0 11`
        # '0 11' according to {GATK} SplitIntervals  --scatter-count 12
        do
            outfile={wildcards.time}/s03_variant/GATK/{wildcards.sample}_${{i}}.CNNscore.vcf.gz
            {GATK} --java-options '-Xms2G -Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4'  CNNScoreVariants \
                    -I {input.bam} -V {input.vcf} -R {REFERENCE} \
                    -L {wildcards.time}/s03_variant/GATK/interval-files/${{i}}-scattered.interval_list \
                    -O ${{outfile}} \
                    -tensor-type read_tensor &
        done
        wait
    done

    find {wildcards.time}/s03_variant/GATK/ -type f -name "{wildcards.sample}_*.CNNscore.vcf.gz" >\
            {wildcards.time}/s03_variant/GATK/{wildcards.sample}_CNNscore.list
    find {wildcards.time}/s03_variant/GATK/ -type f -name "{wildcards.sample}_*.CNNscore.vcf.gz.tbi" >\
            {wildcards.time}/s03_variant/GATK/{wildcards.sample}_CNNscore.tbi.list
    {GATK} --java-options '-Xmx10g -XX:+UseParallelGC  -XX:ParallelGCThreads=2' MergeVcfs \
            -I {wildcards.time}/s03_variant/GATK/{wildcards.sample}_CNNscore.list \
            -O {output.anno}

    echo "Working on FilterVariantTranches"
    {GATK} --java-options '-Xms2G -Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=2' FilterVariantTranches \
        -V {output.anno} -O {output.filt} \
        -resource {HAPMAP} -resource {OMNI} -resource {PHASE1} \
        -resource {DBSNP} -resource {MILLS} \
        --info-key CNN_1D --snp-tranche 99.95 --indel-tranche 99.4 2>>{log}

    # temp files clenup
    xargs rm < {wildcards.time}/s03_variant/GATK/{wildcards.sample}_CNNscore.list
    xargs rm < {wildcards.time}/s03_variant/GATK/{wildcards.sample}_CNNscore.tbi.list
    """

############
# clinSV
# Minoche, A.E., Lundie, B., Peters, G.B. et al. ClinSV: clinical grade structural and copy number variant detection from whole genome sequencing data. Genome Med 13, 32 (2021). https://doi.org/10.1186/s13073-021-00841-x
############
rule clinSV:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam"
    output:
        vcf="{time}/s03_variant/clinSV/{sample}/SVs/joined/SV-CNV.PASS.vcf",
        txtpass="{time}/s03_variant/clinSV/{sample}/SVs/joined/SV-CNV.RARE_PASS_GENE.txt",
        txt="{time}/s03_variant/clinSV/{sample}/SVs/joined/SV-CNV.txt"
    threads: 8
    resources: mem_mb=36000
    log:
        "{time}/s03_variant/clinSV/{sample}_clinSV.log"
    singularity:
        "/HD101TB/bioinfo/liuxingmin/Clinical/WGS_cnv/clinsv_1.0.sif"
    #container:
    #    image("/HD101TB/bioinfo/liuxingmin/Clinical/WGS_cnv/clinsv_1.0.sif", singularity_args="--fakeroot -B $PWD:/root")
    shell: """
        #######################
        clinsv -i "/root/{wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.marked.BQSR.bam" \
                -ref /root/clinsv/refdata-b38 -p /root/{wildcards.time}/s03_variant/clinSV/{wildcards.sample}
        #######################
    """

############
# ClassfyCNV
# Gurbich, T.A., Ilinsky, V.V. ClassifyCNV: a tool for clinical annotation of copy-number variants. Sci Rep 10, 20375 (2020). https://doi.org/10.1038/s41598-020-76425-3
############
rule ClassfyCNV:
    input:
        "{time}/s03_variant/clinSV/{sample}/SVs/joined/SV-CNV.txt",
        "{time}/s03_variant/clinSV/{sample}/SVs/joined/SV-CNV.RARE_PASS_GENE.txt"
    output:
        cnv="{time}/s03_variant/clinSV/{sample}/SVs/joined/CNV.txt",
        cnvpass="{time}/s03_variant/clinSV/{sample}/SVs/joined/CNV.RARE_PASS_GENE.txt"
    threads: 
        2
    resources: mem_mb=4000
    log:
        "{time}/s03_variant/clinSV/{sample}_ClassfyCNV.log"
    shell:"""
        #######################
        grep -v LOW {input[0]} | cut -f18,19 | sed '1d;s/:/\t/;s/-/\t/' | awk '$4 =="DUP" || $4 =="DEL" ' > {output.cnv}
        /HD101TB/bioinfo/software/ClassifyCNV/ClassifyCNV.py \
            --infile {output.cnv} --GenomeBuild hg38 \
            --outdir {wildcards.time}/s03_variant/clinSV/{wildcards.sample}/SVs/joined/ClassifyCNV_results/

        grep -v LOW {input[1]} | cut -f18,19 | sed '1d;s/:/\t/;s/-/\t/' | awk '$4 =="DUP" || $4 =="DEL" ' > {output.cnvpass}
        /HD101TB/bioinfo/software/ClassifyCNV/ClassifyCNV.py \
            --infile {output.cnvpass} --GenomeBuild hg38 \
            --outdir {wildcards.time}/s03_variant/clinSV/{wildcards.sample}/SVs/joined/ClassifyCNV_results_RARE_PASS_GENE
        #######################
        """

#
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
        "{time}/s03_variant/cnvkit/{sample}_cnvkit.log"
    shell: """
	#######################
        #!/bin/bash
        echo "Working on CNVkit"
        #source activate Genomics
	cnvkit.py batch {input.bam} -r {PON} -p {threads} --method wgs \
                --drop-low-coverage --scatter --diagram -d {wildcards.time}/s03_variant/cnvkit/ 2>&1 >{log}
	#cnvkit.py scatter -s {wildcards.time}/s03_variant/cnvkit/{wildcards.sample}.cn{{s,r}} \
        #        -o s03_variant/cnvkit/{wildcards.sample}.pdf
	#######################
	"""
############

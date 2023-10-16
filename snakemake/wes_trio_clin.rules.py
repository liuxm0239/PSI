#########################################
# SNP calling using GATK joint-calling  for Trio samples
# by Xingmin Liu
# https://hpc.nih.gov/training/gatk_tutorial/
# https://www.ibm.com/downloads/cas/ZJQD0QAL 
#########################################
from datetime import date

SAMPLES={"MS22081258F1", "MS22081258M1", "MS22081258P1"}
FAMILY={"MS22081258"}
TARGET="/SlurmDatabase/Clinical/ClinSV/reference/WES_target_region/SeqCap_EZ_MedExomePlusMito_hg38_capture_targets.bed"
REFERENCE="/SlurmDatabase/Clinical/ClinSV/reference/hg38/hg38.fasta"
REFBED="/root/conda/Genome/hg38.target.bed"
BWAINDEX="/SlurmDatabase/Clinical/ClinSV/reference/hg38/hg38.fasta"
GATK="/SlurmDatabase/Bladder_cancer/DNA/rawdata/tools/gatk-4.3.0.0/gatk"
HAPMAP="/root/conda/Genome/hapmap_3.3.hg38.vcf.gz"
DBSNP="/root/conda/Genome/dbsnp_146.hg38.vcf.gz"
OMNI="/root/conda/Genome/1000G_omni2.5.hg38.vcf.gz"
PHASE1="/root/conda/Genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
MILLS="/root/conda/Genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
TIME=date.today().isoformat()
#TIME="2023-07-04"

rule all:
    input:
        expand("{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam", sample=SAMPLES, time=TIME ),
        expand("{time}/s03_variant/GATK/{family}_HC_calls_combined.g.vcf.gz", family=FAMILY, time=TIME),
        expand("{time}/s03_variant/GATK/{family}_merged_VQSR_SNP_INDEL.vcf.gz", family=FAMILY, time=TIME)

rule FASTP:
    input:
        R1="raw.fastq/{sample}_1.fq.gz",
        R2="raw.fastq/{sample}_2.fq.gz"
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
        raw_bam="{time}/s02_alignment/s020_bwa/{sample}.bam",
        bam="{time}/s02_alignment/s020_bwa/{sample}.sort.bam",
        bai="{time}/s02_alignment/s020_bwa/{sample}.sort.bam.bai"
    threads: 20
    resources: mem_mb=25000
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
        samtools sort -m 3G -@ {threads} -T {wildcards.time}/s02_alignment/s020_bwa/tmp_{wildcards.sample} -o {output.bam} {TIME}/s02_alignment/s020_bwa/{wildcards.sample}.bam 2>{log}

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
        bam="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam",
        index="{time}/s02_alignment/s021_Dedup/{sample}.dd.bam.bai",
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
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_s022_Brecal.log"
    shell: """
	#######################
        echo "Working on SplitIntervals for parallel BQSR"
        {GATK} SplitIntervals -R {REFERENCE} -L {TARGET} --scatter-count 20 -O {wildcards.time}/s02_alignment/s022_Brecal/interval-files

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
            find {wildcards.time}/s02_alignment/s022_Brecal/ -type f -name "{wildcards.sample}*.table" > {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.recal.list
            {GATK} GatherBQSRReports -I {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.recal.list -O {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}.recal.table && \
                    rm -f {wildcards.time}/s02_alignment/s022_Brecal/{wildcards.sample}*.table

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
# HaplotypeCaller gVCF                                                                                                                                                
############                                                                                                                                                          
rule gVCF:
    input:
        bam="{time}/s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam"
    output:
        "{time}/s03_variant/GATK/{sample}.g.vcf.gz"
    threads: 2
    resources:
        mem_mb=20000 
    conda: 
        "Genomics"
    log:
        "{time}/s03_variant/{sample}_GATK.log"
    shell:"""
        #######################
        #!/bin/bash
        echo "Working on HaplotypeCaller gVCF"
        gatk --java-options "-Xms20G -Xmx20G -XX:ParallelGCThreads=2" HaplotypeCaller\
          -R {REFERENCE} -I {input.bam} -L {TARGET} \
          -O {output[0]} -ERC GVCF
        #######################
        """
############
# HaplotypeCaller combineGVCFs
############
rule combineGVCFs:
    input:
        gvcfs = lambda wildcards:expand('{time}/s03_variant/GATK/{sample}.g.vcf.gz', sample=SAMPLES, time=TIME)
    output:
        combined = '{time}/s03_variant/GATK/{family}_HC_calls_combined.g.vcf.gz'
    conda:
        "Genomics"
    log:
        '{time}/s03_variant/GATK/{family}_HC_calls_combined.g.vcf.log'
    threads: 8
    shell:"""
    #!/bin/bash
    gatk --java-options '-Xmx60g' CombineGVCFs \
        -R {REFERENCE} \
        -V {input.gvcfs[0]} \
        -V {input.gvcfs[1]} \
        -V {input.gvcfs[2]} \
        -L {TARGET} \
        --create-output-variant-index \
        -O {output.combined}
    """


rule GenotypeGVCFs:
    input:
        "{time}/s03_variant/GATK/{family}_HC_calls_combined.g.vcf.gz"
    output:
        "{time}/s03_variant/GATK/{family}_HC_calls_combined.vcf.gz"
    conda:
        "Genomics"
    log:
        '{time}/s03_variant/GATK/{family}_GenotypeGVCFs.log'
    threads: 8
    shell:"""
    #!/bin/bash
    gatk --java-options '-Xms2G -Xmx2G -XX:ParallelGCThreads=2' GenotypeGVCFs \
        -R {REFERENCE} \
        -V {input[0]} \
        -O {output[0]}
    """

rule VariantRecalibrator:
    input:
        "{time}/s03_variant/GATK/{family}_HC_calls_combined.vcf.gz"
    output:
        snprecal="{time}/s03_variant/GATK/{family}_SNP.recal",
        snptranches="{time}/s03_variant/GATK/{family}_SNP.tranches",
        snpplot="{time}/s03_variant/GATK/{family}_SNP.plots.R",
        indelrecal="{time}/s03_variant/GATK/{family}_INDEL.recal",
        indeltranches="{time}/s03_variant/GATK/{family}_INDEL.tranches",
        indelplot="{time}/s03_variant/GATK/{family}_INDEL.plots.R",
        snpvqsr="{time}/s03_variant/GATK/{family}_merged_VQSR_SNP.vcf.gz",
        indelvqsr="{time}/s03_variant/GATK/{family}_merged_VQSR_SNP_INDEL.vcf.gz"
    conda:
        "Genomics"
    log:
        '{time}/s03_variant/GATK/{family}_VQSR.log'
    threads: 8
    shell:"""
    #!/bin/bash
    gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
        -tranche 100.0 -tranche 99.95 -tranche 99.9 \
        -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
        -tranche 95.0 -tranche 94.0 \
        -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
        -R {REFERENCE} \
        -V {input[0]} \
        --resource:hapmap,known=false,training=true,truth=true,prior=15.0 {HAPMAP} \
        --resource:omni,known=false,training=true,truth=false,prior=12.0 {OMNI} \
        --resource:1000G,known=false,training=true,truth=false,prior=10.0 {PHASE1} \
        -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
        -mode SNP -O {output.snprecal} --tranches-file {output.snptranches} \
        --rscript-file {output.snpplot}

    gatk --java-options "-Xms4G -Xmx4G -XX:ParallelGCThreads=2" VariantRecalibrator \
        -tranche 100.0 -tranche 99.95 -tranche 99.9 \
        -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 \
        -tranche 95.0 -tranche 94.0 \
        -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
        -R {REFERENCE} \
        -V {input[0]} \
        --resource:mills,known=false,training=true,truth=true,prior=12.0 {MILLS} \
        --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 {DBSNP} \
        -an QD -an MQRankSum -an ReadPosRankSum -an FS -an SOR -an DP \
        -mode INDEL -O {output.indelrecal} --tranches-file {output.indeltranches} \
        --rscript-file {output.indelplot}

    gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
        -V {input[0]} \
        -mode SNP \
        --recal-file {output.snprecal} \
        --tranches-file {output.indeltranches} \
        --truth-sensitivity-filter-level 99.9 \
        --create-output-variant-index true \
        -O {output.snpvqsr}

    gatk --java-options "-Xms2G -Xmx2G -XX:ParallelGCThreads=2" ApplyVQSR \
        -V {output.snpvqsr} \
        -mode INDEL \
        --recal-file {output.indelrecal} \
        --tranches-file {output.indeltranches} \
        --truth-sensitivity-filter-level 99.9 \
        --create-output-variant-index true \
        -O {output.indelvqsr}


    """


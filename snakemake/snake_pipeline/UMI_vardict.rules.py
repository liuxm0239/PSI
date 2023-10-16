#########################################
# SNV calling using vardict for UMI samples
# by 
#########################################
from datetime import date

SAMPLES={"U1240_P1", "U1240_BC" }
SPARA=2
BWAINDEX="/SlurmDatabase/reference/xukejin/SNP/NEW/01.align/index/hg38"
REFERENCE="/SlurmDatabase/reference/xukejin/SNP/NEW/01.align/index/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
REFBED="/SlurmDatabase/reference/xukejin/SNP/NEW/01.align/Nonacus_Pan_Cancer_524_TMB_Panel_GRCh38_1x_v1.0_target_merged.bed"
TSOMATIC="/root/conda/miniconda3/envs/Genomics/share/vardict-java-1.8.3-0/bin/testsomatic.R"
VARPAIR="/root/conda/miniconda3/envs/Genomics/share/vardict-java-1.8.3-0/bin/var2vcf_paired.pl"
REFSD="/SlurmDatabase/reference/xukejin/SNP/NEW/01.align/index/Homo_sapiens.GRCh38.dna.primary_assembly.dict"
#TIME=date.today().isoformat()
TIME="2023-07-02"

rule all:
    input:
        expand("{time}/s02_alignment/s021_gencore/{sample}.gencore.{spara}.sorted.bam", sample=SAMPLES, time=TIME, spara=SPARA ),
        expand("{time}/s03_somatic_variant/vardict/{sample}.gencore.{spara}.vcf", sample=SAMPLES, time=TIME, spara=SPARA)

rule FASTP:
    input:
        R1="fastq_data/{sample}_R1.fastq.gz",
        R2="fastq_data/{sample}_R2.fastq.gz"
    output:
        "{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        "{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    threads: 2    
    resources: mem_mb=6000
    conda:
        "Genomics"
    log:
        "{time}/s01_fastqc/{sample}.fastp_fastqc.log"
    shell: """
        #######################
        #echo "Working on raw data fastqc"
        #fastqc -t {threads} {input.R1} {input.R2} -o {wildcards.time}/s01_fastqc/ > {log} 2>&1

        echo "Working on fastp triming"
        #fastp --in1 {input.R1} --in2 {input.R2} --out1 {output[0]} --out2 {output[1]} -5 -c -p -q 20 --json {wildcards.time}/s01_fastqc/{wildcards.sample}.fastp.Clean.json --html {wildcards.time}/s01_fastqc/{wildcards.sample}.fastp.Clean.html --thread  {threads} --trim_poly_g --detect_adapter_for_pe >> {log} 2>&1 
        R1=$(realpath {input.R1})
        R2=$(realpath {input.R2})
        ln -s  ${{R1}} {output[0]}
        ln -s  ${{R2}} {output[1]}

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
        bai="{time}/s02_alignment/s020_bwa/{sample}.sort.bam.bai",
        bcount="{time}/s02_alignment/s020_bwa/{sample}.sort.bam.txt"
    threads: 20    
    resources: mem_mb=16000
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_alignment.log"
    params:
        platform ="Illumina"
    shell: """
        #######################
        echo "Working on BWA alignment"
        bwa mem -M -t {threads} -R '@RG\\tID:{wildcards.sample}\\tSM:{wildcards.sample}\\tLB:{wildcards.sample}\\tPL:{params.platform}' {BWAINDEX} {input.R1} {input.R2} | samtools view - -S -b -F 0x900 -o {output.raw_bam} 2>&1 >{log}

        echo "Working on samtools sorting"
        samtools sort -m 4G -@ {threads} -T {wildcards.time}/s02_alignment/s020_bwa/tmp -o {output.bam} {output.raw_bam} 2>>{log}
        samtools index -@ {threads} {output.bam}

        echo "Working on bam counting"
        samtools view -@ {threads} -c {output.bam} > {output.bcount} 2>>{log}
        #######################
        """

rule gencore:
    input:
        bam="{time}/s02_alignment/s020_bwa/{sample}.sort.bam",
        bai="{time}/s02_alignment/s020_bwa/{sample}.sort.bam.bai"
    output:
        bam="{time}/s02_alignment/s021_gencore/{sample}.gencore.{spara}.bam",
        stat="{time}/s02_alignment/s021_gencore/{sample}.gencore.{spara}.bam.txt",
        sortedbam="{time}/s02_alignment/s021_gencore/{sample}.gencore.{spara}.sorted.bam",
        index="{time}/s02_alignment/s021_gencore/{sample}.gencore.{spara}.sorted.bam.bai"
    threads: 8      
    resources: mem_mb=20000
    conda:
        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_s021_gencore.{spara}.log"
    shell: """
        #######################
        echo "Working on gencore"
        /SlurmDatabase/reference/xukejin/SNP/NEW/tools/gencore -i {input.bam} -o  {output.bam} -r {REFERENCE} -s {wildcards.spara} --score_threshold {threads} --high_qual 30 
        echo "Working on umi_dedup_bam_static"

        samtools sort -@ {threads} -m 4G -T {wildcards.time}/s02_alignment/s021_gencore/tmp -o {output.sortedbam} {output.bam}
        samtools index -@ {threads} {output.sortedbam}
        samtools view -@ {threads} -c {output.sortedbam} > {output.stat}
        #######################
        """


############
#    vardict
############
rule vardict:
    input:
        bam="{time}/s02_alignment/s021_gencore/{sample}.gencore.{spara}.sorted.bam"
    output:
        tmp="{time}/s03_somatic_variant/vardict/{sample}.tmp.gencore.{spara}.vcf",
        vcf="{time}/s03_somatic_variant/vardict/{sample}.gencore.{spara}.vcf"
    threads: 4    
    resources: mem_mb=28000
    conda:
        "Genomics"
    log:
        "{time}/s03_somatic_variant/{sample}.gencore.{spara}.log"
    shell: """
        #######################
        #!/bin/bash
        echo "Working on vardict"
        vardict-java -G {REFERENCE} -N {wildcards.sample} -f 0.01 -th {threads} -b {input.bam} -O 15 -q 30 -Q 40 -r 4 -z -c 1 -S 2 -E 3 -g 4 {REFBED} | \
                {TSOMATIC} |\
                {VARPAIR} -N "{wildcards.sample}|U1183_BC" -f 0.01 |\
                awk '{if ($1 ~ /^#/) print; else print }' > {output.tmp}
        picard -Xmx128g SortVcf -I {output.tmp} -O {output.vcf} -SD {REFSD}
        #######################
        """
############

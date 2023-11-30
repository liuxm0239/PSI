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
TIME=config['TIME']

#SAMPLES={ "TB1", "TB2"}
REFERENCE="/HD101TB/bioinfo/reference/hg38/"
REFBED="/HD101TB/bioinfo/reference/hg38/hg38.bed"
#BWAINDEX="/HD101TB/conda/Genome/bwaIndex/hg38"
GATK="/opt/gatk/gatk-4.3.0.0/gatk"
HAPMAP="/HD101TB/conda/Genome/hapmap_3.3.hg38.vcf.gz"
DBSNP="/HD101TB/conda/Genome/dbsnp_138.hg38.vcf.gz"
OMNI="/HD101TB/conda/Genome/1000G_omni2.5.hg38.vcf.gz"
PHASE1="/HD101TB/conda/Genome/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
MILLS="/HD101TB/conda/Genome/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
PON="/HD101TB/conda/Genome/cnvkit/PoN/Reference.cnn"
#TIME="2023-08-20"

def double_threads(wildcards, threads):
    return threads * 2

rule all:
    input:
        expand("{time}/s02_alignment/s021_dedup/{sample}.bismark.deduplicated.bam", sample=SAMPLES, time=TIME),
        expand("{time}/s03_methylation/{sample}.bismark_report.html", sample=SAMPLES, time=TIME),
        expand("{time}/s03_methylation/cgmaptools/{sample}.mstat.data", sample=SAMPLES, time=TIME)

rule trim_galore:
    input:
        R1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        R2 = lambda wildcards: FILES[wildcards.sample]['R2']
        #R1="raw.fastq/{sample}_R1.fastq.gz",
        #R2="raw.fastq/{sample}_R2.fastq.gz"
    output:
        "{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        "{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    threads: 2
    resources: mem_mb=2000
    #conda:        "Genomics"
    log:
        "{time}/s01_fastqc/{sample}.trim_galore.log"
    shell: """
        #######################
        if [ -f {wildcards.time}/s01_fastqc/{wildcards.sample}_R1.fq.gz ]; then
            ls {wildcards.time}/s01_fastqc/{wildcards.sample}* >\
                    {wildcards.time}/s01_fastqc/{wildcards.sample}_old
            xargs rm < {wildcards.time}/s01_fastqc/{wildcards.sample}_old
        fi

        echo "Working on trim_galore triming"
        ln -s $(realpath {input.R1}) {wildcards.time}/s01_fastqc/{wildcards.sample}_R1.fq.gz
        ln -s $(realpath {input.R2}) {wildcards.time}/s01_fastqc/{wildcards.sample}_R2.fq.gz
        trim_galore --rrbs --paired -j 4 -o {wildcards.time}/s01_fastqc/\
                {wildcards.time}/s01_fastqc/{wildcards.sample}_R1.fq.gz\
                {wildcards.time}/s01_fastqc/{wildcards.sample}_R2.fq.gz

        mv {wildcards.time}/s01_fastqc/{wildcards.sample}_R1_val_1.fq.gz {output[0]}
        mv {wildcards.time}/s01_fastqc/{wildcards.sample}_R2_val_2.fq.gz {output[1]}
	#######################
    """

rule fastQC:
    input:
        R1 = lambda wildcards: FILES[wildcards.sample]['R1'],
        R2 = lambda wildcards: FILES[wildcards.sample]['R2'],
        R1_clean =  "{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        R2_clean =  "{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
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

rule Bismark:
    input: 
        R1="{time}/s01_fastqc/{sample}_clean_R1.fastq.gz",
        R2="{time}/s01_fastqc/{sample}_clean_R2.fastq.gz"
    output:
        bam=temp("{time}/s02_alignment/s020_bismark/{sample}.bismark.bam")
    threads: 6
    resources: mem_mb=50000
    #conda:        "Genomics"
    log:
        "{time}/s02_alignment/{sample}_alignment.log"
    params:
        platform = config['PLATFORM']
    shell: """
        #######################	
        echo "Working on Bismark alignment"
        #'--score_min L,0,-0.2' for bowtie2 end 2 end alignment
        bismark --parallel {threads} --bowtie2\
            --rg_tag --rg_id {wildcards.sample} --rg_sample {wildcards.sample}\
            --score_min L,0,-0.2 {REFERENCE} \
            -1 {input.R1} -2 {input.R2} \
            --temp_dir {wildcards.time}/s02_alignment/s020_bismark \
            --output_dir {wildcards.time}/s02_alignment/s020_bismark \
            2>logs/s020_bismark.log &&\
        reads_1=$(echo {input.R1} | cut -d/ -f3 | sed 's/.fq.gz//;s/.fastq.gz//;s/.fq//;s/.fastq//')
        mv {wildcards.time}/s02_alignment/s020_bismark/${{reads_1}}_bismark_bt2_pe.bam \
            {wildcards.time}/s02_alignment/s020_bismark/{wildcards.sample}.bismark.bam &&\
        mv {wildcards.time}/s02_alignment/s020_bismark/${{reads_1}}_bismark_bt2_PE_report.txt \
        {wildcards.time}/s02_alignment/s020_bismark/{wildcards.sample}.bismark_report.txt
	#######################
	"""

rule methylation_extractor:
    input:
        bam="{time}/s02_alignment/s020_bismark/{sample}.bismark.bam"
    output:
        bam="{time}/s03_methylation/{sample}.bismark_report.html",
        cxreport="{time}/s03_methylation/{sample}.bismark.CX_report.txt.gz"
    threads: 6
    resources: mem_mb=16000
    #conda:        "Genomics"
    log:
        "{time}/s03_methylation/{sample}_s03_methylation.log"
    shell: """
	#######################
        echo "Working on bismark_methylation_extractor"
        bismark_methylation_extractor -p --no_overlap --parallel {threads} --buffer_size 8G \
            --ignore_r2 3 --cytosine_report --CX_context --comprehensive \
            --gzip -bedGraph --genome_folder {REFERENCE} -o {wildcards.time}/s03_methylation \
            {input.bam}  2>{log}

        bismark2report --alignment_report {wildcards.time}/s02_alignment/s020_bismark/{wildcards.sample}.bismark_report.txt \
            --splitting_report  {wildcards.time}/s03_methylation/{wildcards.sample}.bismark_splitting_report.txt \
            --mbias_report {wildcards.time}/s03_methylation/{wildcards.sample}.bismark.M-bias.txt \
            -o {wildcards.time}/s03_methylation/{wildcards.sample}.bismark_report.html 2>>{log}
	#######################
	"""

rule cgmaptools:
    input:
        cxreport="{time}/s03_methylation/{sample}.bismark.CX_report.txt.gz"
    output:
        cgmap="{time}/s03_methylation/cgmaptools/{sample}.CGmap.gz",
        mbin="{time}/s03_methylation/cgmaptools/{sample}.mbin.data",
        mecs="{time}/s03_methylation/cgmaptools/{sample}.mec_stat.data",
        mstat="{time}/s03_methylation/cgmaptools/{sample}.mstat.data"
    threads: 2
    resources: mem_mb=6000
    conda: "py27"
    log:
        "{time}/s03_methylation/{sample}_s031_cgmaptools.log"
    shell:"""
    #######################
    cgmaptools convert bismark2cgmap -i {input.cxreport} -o {output.cgmap} 2>{log}

    cgmaptools mbin -i {output.cgmap} -B 100000 -c 5 -f png \
        -p {wildcards.time}/s03_methylation/cgmaptools/{wildcards.sample} \
        -t {wildcards.sample} > \
        {output.mbin} 2>>{log}

    cgmaptools mec stat -i {output.cgmap} \
        -p {wildcards.time}/s03_methylation/cgmaptools/{wildcards.sample}  -f png >\
        {output.mecs} 2>>{log} 

    cgmaptools mstat -i {output.cgmap} -c 5 -f png \
        -p {wildcards.time}/s03_methylation/cgmaptools/{wildcards.sample} \
        -t {wildcards.sample} > \
        {output.mstat} 2>>{log} 
    #######################
    """



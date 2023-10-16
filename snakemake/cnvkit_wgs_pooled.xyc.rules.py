#########################################
# CNV calling using cnvkit for WGS samples
# by Xingmin Liu
#########################################

SAMPLES={"xiaoyecui"}
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

rule all:
    input:
        expand("s03_somatic_variant/cnvkit/{sample}.marked.BQSR.call.cns", sample=SAMPLES)


############
#       cnvkit
############
rule cnvkitPooled:
    input:
        bam="s02_alignment/s022_Brecal/{sample}.marked.BQSR.bam"
    output:
        "s03_somatic_variant/cnvkit/{sample}.marked.BQSR.call.cns"
    resources: mem_mb=4000, cpus=2
    conda:
        "Genomics"
    log:
        "s03_somatic_variant/{sample}_cnvkit.log"
    shell: """
	#######################
        #!/bin/bash
        echo "Working on CNVkit"
        #source activate Genomics
	cnvkit.py batch {input.bam} -r {PON} -p {resources.cpus} --method wgs --drop-low-coverage --scatter --diagram -d s03_somatic_variant/cnvkit/ >{log} 2>&1
	#cnvkit.py scatter -s s03_somatic_variant/cnvkit/{wildcards.sample}.cn{{s,r}} -o s03_somatic_variant/cnvkit/{wildcards.sample}.pdf
	#######################
	"""
############

nohup snakemake -p --use-singularity --singularity-args " --fakeroot --bind $PWD:/root --bind /HD101TB/bioinfo/liuxingmin/Clinical/WGS_cnv/clinsv/:/root/clinsv/ " --rerun-incomplete --profile slurm --jobs 28 -s /HD101TB/bioinfo/liuxingmin/Clinical/snk/wgs_clin_v1.0.0.rules.py &
snakemake --unlock ./ -s /HD101TB/bioinfo/liuxingmin/Clinical/snk/wgs_matepair_v1.0.0.rules.py
snakemake -s /HD101TB/bioinfo/liuxingmin/Clinical/snk/wgs_clin_v1.0.0.rules.py --rulegraph | dot -Tpng > wgs_clin.rules.png

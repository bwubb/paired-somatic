import os

### INIT ###

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

### FUNCTIONS ###

def paired_bams(wildcards):
    ref=config['reference']['key']
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

### SNAKEMAKE ###

wildcard_constraints:
    work_dir=f"data/work/{config['resources']['targets_key']}"

rule collect_sequenza:
    input:
        expand("data/work/{lib}/{tumor}/sequenza/{tumor}_segments.txt",lib=f"{config['resources']['targets_key']}",tumor=PAIRS.keys())


#This can be parallelized; divided into chr.seqz.gz
#probably combined after binning. each file has same header
#zcat sample_chr1.seqz.gz sample_chr1.seqz.gz | \
#    gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip > \
#    sample.seqz.gz
#tabix -f -s 1 -b 2 -e 2 -S 1 sample.seqz.gz
#I have done this already in another script, but it has not been added.

#awk 'FNR==1 && NR!=1 { while (/^HRD/) getline; } 1 {print}' *hrd.txt



rule Sequenza_bam2seqz:
    input:
        unpack(paired_bams)
    params:
        ref=config['reference']['fasta'],
        gc=config['reference']['gc_wiggle']
        #"$HOME/resources/Genomes/Human/GRCh37/custom-GRCh37.gc50Base.txt.gz"#different file for tcga Homo_sapien_assembly19
    output:
        "{work_dir}/{tumor}/sequenza/seqz.gz"
    shell:
        #"sequenza-utils.py pileup2seqz -gc {params.gc} -n {input.normal} -t {input.tumor} | gzip > {output}"
        "sequenza-utils bam2seqz -F {params.ref} -gc {params.gc} -n {input.normal} -t {input.tumor} | gzip > {output}"
        #CHECK IF THE FILE IS EMPTY
#Make workflow to create gc file

rule Sequenza_bin:
    input:
        "{work_dir}/{tumor}/sequenza/seqz.gz"
    output:
        "{work_dir}/{tumor}/sequenza/seqz.small.txt"
    params:
        bin=50
        #50 for exome 200 for genome
    shell:
        "sequenza-utils seqz_binning -w {params.bin} -s {input} -o {output}"

#need to remove Y MT M
#At somepoint sequenza started having trouble with gzipped files.
rule Sequenza_extract:
    input:
        "{work_dir}/{tumor}/sequenza/seqz.small.txt"
    output:
        "{work_dir}/{tumor}/sequenza/{tumor}_confints_CP.txt",
        "{work_dir}/{tumor}/sequenza/{tumor}_segments.txt"
    params:
        outdir="{work_dir}/{tumor}/sequenza"
    threads:
        4
    shell:
        "Rscript sequenza-snakemake.R --id {wildcards.tumor} --input {input} --outdir {params.outdir} --threads {threads}"

rule Sequenza_hrd:
    input:
        segments="{work_dir}/{tumor}/sequenza/{tumor}_segments.txt",
        confints="{work_dir}/{tumor}/sequenza/{tumor}_confints_CP.txt",
    output:
        "{tumor}_hrd.txt"
    script:
        "Rscript $HOME/software/BRADtools/beta/hrd_wrapper.R {tumor} {input.segments} {input.confints}"

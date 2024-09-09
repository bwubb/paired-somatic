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

CHR=[f'chr{c}' for c in range(1,23)]+['chrX']
print(CHR)

### FUNCTIONS ###

def paired_bams(wildcards):
    ref=config['reference']['key']
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

### SNAKEMAKE ###

wildcard_constraints:
    work_dir=f"data/work/{config['resources']['targets_key']}",
    chr="chr[1-2][0-9]|chr[1-9]|chrX"

rule sequenza_all:
    input:
        expand("data/work/{targets}/{tumor}/sequenza/{tumor}.sequenza.annotsv.gene_split.report.csv",targets=f"{config['resources']['targets_key']}",tumor=PAIRS.keys()),
        expand("data/work/{targets}/{tumor}/sequenza/{tumor}.sequenza.hrd.txt",targets=f"{config['resources']['targets_key']}",tumor=PAIRS.keys()),
        "sequenza_purity.table",
        "sequenza_ploidy.table"

#Make workflow to create gc file
#
rule sequenza_bam2seqz_byChr:
    input:
        unpack(paired_bams)
    output:
        temp("data/work/{targets}/{tumor}/sequenza/seqz.{chr}.gz")
    params:
        ref=config['reference']['fasta'],
        gc=config['reference']['gc_wiggle'],
        chr=lambda wildcards: wildcards.chr.lstrip('chr')
    shell:
        "sequenza-utils bam2seqz -F {params.ref} -gc {params.gc} -n {input.normal} -t {input.tumor} -C {params.chr} | gzip > {output}"
#CHECK IF THE FILE IS EMPTY

rule sequenza_bin_byChr:
    input:
        "data/work/{targets}/{tumor}/sequenza/seqz.{chr}.gz"
    output:
        temp("data/work/{targets}/{tumor}/sequenza/seqz.{chr}.small.gz")
    params:
        #50 for exome 200 for genome
        bin=50
    shell:
        """
        sequenza-utils seqz_binning -w {params.bin} -s {input} -o - | gzip > {output}
        """

rule sequenza_combine_seqz:
    input:
        ["data/work/{{targets}}/{{tumor}}/sequenza/seqz.{chr}.small.gz".format(chr=chr) for chr in CHR]
    output:
        "data/work/{targets}/{tumor}/sequenza/seqz.small.all.gz"
    shell:
        """
        zcat {input} | gawk '{{if (NR!=1 && $1 != \"chromosome\") {{print $0}}}}' | bgzip > {output}
        tabix -f -s 1 -b 2 -e 2 -S 1 {output}
        """

rule sequenza_extract:
    input:
        "data/work/{targets}/{tumor}/sequenza/seqz.small.all.gz"
    output:
        "data/work/{targets}/{tumor}/sequenza/{tumor}_confints_CP.txt",
        "data/work/{targets}/{tumor}/sequenza/{tumor}_segments.txt"
    params:
        outdir="data/work/{targets}/{tumor}/sequenza"
    threads:
        4
    shell:
        """
        Rscript sequenza-snakemake.R --id {wildcards.tumor} --input {input} --outdir {params.outdir} --threads {threads}
        """

rule sequenza_2bed:
    input:
        "data/work/{targets}/{tumor}/sequenza/{tumor}_segments.txt"
    output:
        "data/work/{targets}/{tumor}/sequenza/{tumor}_segments.bed"
    shell:
        """
        python cnv_to_bed.py -c sequenza {input}
        """

rule sequenza_AnnotSV:
    input:
        "data/work/{targets}/{tumor}/sequenza/{tumor}_segments.bed"
    output:
        "data/work/{targets}/{tumor}/sequenza/{tumor}.sequenza.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule sequenza_AnnotSV_parser:
    input:
        "data/work/{targets}/{tumor}/sequenza/{tumor}.sequenza.annotsv.gene_split.tsv"
    output:
        "data/work/{targets}/{tumor}/sequenza/{tumor}.sequenza.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

rule sequenza_HRDex:
    input:
        "data/work/{targets}/{tumor}/sequenza/{tumor}_segments.bed"
    output:
        "data/work/{targets}/{tumor}/sequenza/{tumor}.sequenza.hrd.txt"
    params:
        tumor="{tumor}",
        build="grch38"
    shell:
        """
        Rscript runHRDex.R -i {input} -o {output} --build {params.build} --tumor {params.tumor}
        """

rule sequenza_purity_table:
    input:
        expand("data/work/{targets}/{tumor}/sequenza/{tumor}_confints_CP.txt",targets=config['resources']['targets_key'],tumor=PAIRS.keys())
    output:
        "sequenza_purity.table"
    run:
        with open(output[0],'w') as outfile:
            for f in input:
                sampleid=os.path.basename(f).split('_')[0]
                with open(f,'r') as infile:
                    line=infile.readline()#header
                    line=infile.readline()#line1
                    line=infile.readline()#line2
                    cellularity=line.rstrip().split('\t')[0]
                outfile.write(f'{sampleid}\t{cellularity}\n')

rule sequenza_ploidy_table:
    input:
        expand("data/work/{targets}/{tumor}/sequenza/{tumor}_confints_CP.txt",targets=config['resources']['targets_key'],tumor=PAIRS.keys())
    output:
        "sequenza_ploidy.table"
    run:
        with open(output[0],'w') as outfile:
            for f in input:
                sampleid=os.path.basename(f).split('_')[0]
                with open(f,'r') as infile:
                    line=infile.readline()#header
                    line=infile.readline()#line1
                    line=infile.readline()#line2
                    ploidy=line.rstrip().split('\t')[1]
                outfile.write(f'{sampleid}\t{ploidy}\n')

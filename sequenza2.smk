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

rule collect_sequenza:
    input:
        expand("data/work/{lib}/{tumor}/sequenza/{tumor}_segments.txt",lib=f"{config['resources']['targets_key']}",tumor=PAIRS.keys())

rule collect_hrd:
    input:
        expand("data/work/{lib}/{tumor}/sequenza/{tumor}_hrd.txt",lib=f"{config['resources']['targets_key']}",tumor=PAIRS.keys())

rule collect_purity:
    input:
        "purity.table"

rule collect_ploidy:
    input:
        "ploidy.table"

rule collect_annotsv:
    input:
        expand("data/work/{lib}/{tumor}/sequenza/annotsv_gene_split.report.csv",lib=f"{config['resources']['targets_key']}",tumor=PAIRS.keys())
#This can be parallelized; divided into chr.seqz.gz
#probably combined after binning. each file has same header
#zcat sample_chr1.seqz.gz sample_chr1.seqz.gz | \
#    gawk '{if (NR!=1 && $1 != "chromosome") {print $0}}' | bgzip > \
#    sample.seqz.gz
#tabix -f -s 1 -b 2 -e 2 -S 1 sample.seqz.gz

#awk 'FNR==1 && NR!=1 { while (/^HRD/) getline; } 1 {print}' *hrd.txt

#Make workflow to create gc file

#
rule Sequenza_bam2seqz_byChr:
    input:
        unpack(paired_bams)
    output:
        temp("{work_dir}/{tumor}/sequenza/seqz.{chr}.gz")
    params:
        ref=config['reference']['fasta'],
        gc=config['reference']['gc_wiggle'],
        chr=lambda wildcards: wildcards.chr.lstrip('chr')
    shell:
        "sequenza-utils bam2seqz -F {params.ref} -gc {params.gc} -n {input.normal} -t {input.tumor} -C {params.chr} | gzip > {output}"
#CHECK IF THE FILE IS EMPTY

rule Sequenza_bin_byChr:
    input:
        "{work_dir}/{tumor}/sequenza/seqz.{chr}.gz"
    output:
        temp("{work_dir}/{tumor}/sequenza/seqz.{chr}.small.gz")
    params:
        #50 for exome 200 for genome
        bin=50
    shell:
        """
        sequenza-utils seqz_binning -w {params.bin} -s {input} -o - | gzip > {output}
        """

rule Sequenza_combine_seqz:
    input:
        ["{{work_dir}}/{{tumor}}/sequenza/seqz.{chr}.small.gz".format(chr=chr) for chr in CHR]
    output:
        "{work_dir}/{tumor}/sequenza/seqz.small.all.gz"
    shell:
        """
        zcat {input} | gawk '{{if (NR!=1 && $1 != \"chromosome\") {{print $0}}}}' | bgzip > {output}
        tabix -f -s 1 -b 2 -e 2 -S 1 {output}
        """

rule Sequenza_extract:
    input:
        "{work_dir}/{tumor}/sequenza/seqz.small.all.gz"
    output:
        "{work_dir}/{tumor}/sequenza/{tumor}_confints_CP.txt",
        "{work_dir}/{tumor}/sequenza/{tumor}_segments.txt"
    params:
        outdir="{work_dir}/{tumor}/sequenza"
    threads:
        4
    shell:
        """
        Rscript sequenza-snakemake.R --id {wildcards.tumor} --input {input} --outdir {params.outdir} --threads {threads}
        """

rule Sequenza_2bed:
    input:
        "{work_dir}/{tumor}/sequenza/{tumor}_segments.txt"
    output:
        "{work_dir}/{tumor}/sequenza/{tumor}_segments.bed"
    shell:
        """
        python sequenza2bed.py {input}
        """
        #change it to write to stdout if output arg isnt given?

rule Sequenza_AnnotSV:
    input:
        "{work_dir}/{tumor}/sequenza/{tumor}_segments.bed"
    output:
        "{work_dir}/{tumor}/sequenza/annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule Sequenza_AnnotSV_parser:
    input:
        "{work_dir}/{tumor}/sequenza/annotsv.gene_split.tsv"
    output:
        "{work_dir}/{tumor}/sequenza/annotsv_gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

rule Sequenza_hrd:
    input:
        segments="{work_dir}/{tumor}/sequenza/{tumor}_segments.txt"
    output:
        "{work_dir}/{tumor}/sequenza/{tumor}_hrd.txt"
    params:
        build=config['reference']['key'].lower()
    shell:
        """
        Rscript $HOME/software/HRDex/R/run_HRDex.R -i {input} -o {output} --tumor {wildcards.tumor} --build {params.build}
        """

rule purity_table:
    input:
        expand("data/work/{lib}/{tumor}/sequenza/{tumor}_confints_CP.txt",lib=config['resources']['targets_key'],tumor=PAIRS.keys())
    output:
        "purity.table"
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

rule ploidy_table:
    input:
        expand("data/work/{lib}/{tumor}/sequenza/{tumor}_confints_CP.txt",lib=config['resources']['targets_key'],tumor=PAIRS.keys())
    output:
        "ploidy.table"
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

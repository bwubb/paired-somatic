
##Author: Brad Wubbenhorst

import os

### INIT ###

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('project',{}).get('pair_table','pair.tble'),'r') as p:
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

localrules: lancet_sample_name

rule run_lancet:
    input: expand("data/final/{tumor}/{tumor}.lancet.somatic.final.bcf",tumor=PAIRS.keys())

rule processed_lancet:
    input: expand("data/work/{tumor}/lancet/somatic.norm.clean.vcf.gz",tumor=PAIRS.keys())

rule unprocessed_lancet:
    input: expand("data/work/{tumor}/lancet/somatic.vcf.gz",tumor=PAIRS.keys())

rule lancet_main:
    input:
        unpack(paired_bams)
    output:
        "data/work/{tumor}/lancet/somatic.vcf.gz"
    params:
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed']
    threads:
        4
    shell:
        """
        lancet --tumor {input.tumor} --normal {input.normal} --ref {params.ref} --bed {params.bed} --num-threads {threads} | bcftools view -W=tbi -Oz -o {output}
        """

rule lancet_somatic_normalized:
    input:
        "data/work/{tumor}/lancet/somatic.vcf.gz"
    output:
        norm="data/work/{tumor}/lancet/somatic.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -W=tbi -Oz -o {output.norm}
        """

rule lancet_sample_name:
    output:
        "data/work/{tumor}/lancet/sample.name"
    params:
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        echo -e "TUMOR\\t{wildcards.tumor}\\nNORMAL\\t{params.normal}" > {output}
        echo -e "tumor\\t{wildcards.tumor}\\nnormal\\t{params.normal}" >> {output}
        """

rule lancet_somatic_clean:
    input:
        name="data/work/{tumor}/lancet/sample.name",
        vcf="data/work/{tumor}/lancet/somatic.norm.vcf.gz"
    output:
        clean="data/work/{tumor}/lancet/somatic.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("data/work/{tumor}/lancet/temp.h.vcf.gz")
    shell:
        """
        bcftools reheader -f {params.fai} -s {input.name} -o {params.vcf} {input.vcf}
        bcftools index {params.vcf}

        bcftools view -s {wildcards.tumor},{params.normal} -e 'ALT="*"' -R {params.regions} {params.vcf} | \
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output.clean}
        """

rule lancet_somatic_final:
    input:
        "data/work/{tumor}/lancet/somatic.vcf.gz",
        "data/work/{tumor}/lancet/somatic.norm.clean.vcf.gz"
    output:
        "data/final/{tumor}/{tumor}.lancet.somatic.bcf",
        "data/final/{tumor}/{tumor}.lancet.somatic.final.bcf"
    shell:
        """
        bcftools view -W=csi -Ob -o {output[0]} {input[0]}
        bcftools view -W=csi -Ob -o {output[1]} {input[1]}
        """

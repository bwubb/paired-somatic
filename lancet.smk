
#lancet.snake
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
localrules: collect_lancet

wildcard_constraints:
    work_dir=f"data/work/{config['resources']['targets_key']}"

rule processed_lancet:
    input:
        expand("data/work/{lib}/{tumor}/lancet/somatic.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule unprocessed_lancet:
    input:
        expand("data/work/{lib}/{tumor}/lancet/somatic.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule run_lancet:
    input:
        unpack(paired_bams)
    output:
        "{work_dir}/{tumor}/lancet/somatic.vcf.gz"
    params:
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed']
    threads:
        4
    shell:
        """
        lancet --tumor {input.tumor} --normal {input.normal} --ref {params.ref} --bed {params.bed} --num-threads {threads} | bcftools view -O z -o {output}
        tabix -fp vcf {output}
        """
        ##Reheader?

rule lancet_somatic_normalized:
    input:
        "{work_dir}/{tumor}/lancet/somatic.vcf.gz"
    output:
        norm="{work_dir}/{tumor}/lancet/somatic.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output.norm}
        tabix -fp vcf {output.norm}
        """

rule lancet_somatic_clean:
    input:
        "{work_dir}/{tumor}/lancet/somatic.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/lancet/sample.name",
        clean="{work_dir}/{tumor}/lancet/somatic.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/lancet/temp.h.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

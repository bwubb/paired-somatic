import os
import csv

include:"./sequenza2.smk"

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
    work_dir=f"data/work/{config['resources']['targets_key']}",
    results_dir=f"data/final/{config['project']['name']}"


rule filter_applied_germline_varscan2:
    input:
        expand("data/work/{lib}/{tumor}/varscan2/loh.fpfilter.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{lib}/{tumor}/varscan2/germline.fpfilter.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule bamreadcount_loh_regions:
    input:
        "{work_dir}/{tumor}/varscan2/loh.snp.vcf",
        "{work_dir}/{tumor}/varscan2/loh.indel.vcf"
    output:
        "{work_dir}/{tumor}/varscan2/loh.snp.regions",
        "{work_dir}/{tumor}/varscan2/loh.indel.regions"
    script:
        "bam-readcount_regions.py"

rule bamreadcount_loh_recounts:
    input:
        bam=lambda wildcards: BAMS[PAIRS[wildcards.tumor]],
        snp="{work_dir}/{tumor}/varscan2/loh.snp.regions",
        indel="{work_dir}/{tumor}/varscan2/loh.indel.regions"
    output:
        snp="{work_dir}/{tumor}/varscan2/loh.snp.readcounts",
        indel="{work_dir}/{tumor}/varscan2/loh.indel.readcounts"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        $HOME/software/bam-readcount/bin/bam-readcount -b 15 -q 1 -w 1 -f {params.ref} -l {input.snp} {input.bam} > {output.snp}
        $HOME/software/bam-readcount/bin/bam-readcount -i -b 15 -q 1 -w 1 -f {params.ref} -l {input.indel} {input.bam} > {output.indel}
        """

rule VarScan2_loh_fpfilter:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/loh.snp.vcf",
        indel_vcf="{work_dir}/{tumor}/varscan2/loh.indel.vcf",
        snp_readcount="{work_dir}/{tumor}/varscan2/loh.snp.readcounts",
        indel_readcount="{work_dir}/{tumor}/varscan2/loh.indel.readcounts"
    output:
        snps="{work_dir}/{tumor}/varscan2/loh.snps.fpfilter.vcf.gz",
        indels="{work_dir}/{tumor}/varscan2/loh.indels.fpfilter.vcf.gz"
    params:
        snp_out="{work_dir}/{tumor}/varscan2/loh.snps.fpfilter.vcf",
        indel_out="{work_dir}/{tumor}/varscan2/loh.indels.fpfilter.vcf"
    shell:
        """
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.snp_vcf} {input.snp_readcount} --output-file {params.snp_out}
        bgzip {params.snp_out}
        tabix -f -p vcf {output.snps}
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.indel_vcf} {input.indel_readcount} --output-file {params.indel_out}
        bgzip {params.indel_out}
        tabix -f -p vcf {output.indels}
        """

rule VarScan2_loh_merge:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/loh.snps.fpfilter.vcf.gz",
        indel_vcf="{work_dir}/{tumor}/varscan2/loh.indels.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.vcf.gz"
    shell:
        """
        bcftools concat -a {input} | bcftools sort -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_loh_normalized:
    input:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.norm.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_loh_clean:
    input:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/varscan2/sample.loh.name",
        clean="{work_dir}/{tumor}/varscan2/loh.fpfilter.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/varscan2/temp.loh.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

rule bamreadcount_germline_regions:
    input:
        "{work_dir}/{tumor}/varscan2/germline.snp.vcf",
        "{work_dir}/{tumor}/varscan2/germline.indel.vcf"
    output:
        "{work_dir}/{tumor}/varscan2/germline.snp.regions",
        "{work_dir}/{tumor}/varscan2/germline.indel.regions"
    script:
        "bam-readcount_regions.py"

rule bamreadcount_germline_recounts:
    input:
        bam=lambda wildcards: BAMS[PAIRS[wildcards.tumor]],
        snp="{work_dir}/{tumor}/varscan2/germline.snp.regions",
        indel="{work_dir}/{tumor}/varscan2/germline.indel.regions"
    output:
        snp="{work_dir}/{tumor}/varscan2/germline.snp.readcounts",
        indel="{work_dir}/{tumor}/varscan2/germline.indel.readcounts"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        $HOME/software/bam-readcount/bin/bam-readcount -b 15 -q 1 -w 1 -f {params.ref} -l {input.snp} {input.bam} > {output.snp}
        $HOME/software/bam-readcount/bin/bam-readcount -i -b 15 -q 1 -w 1 -f {params.ref} -l {input.indel} {input.bam} > {output.indel}
        """

rule VarScan2_germline_fpfilter:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/germline.snp.vcf",
        indel_vcf="{work_dir}/{tumor}/varscan2/germline.indel.vcf",
        snp_readcount="{work_dir}/{tumor}/varscan2/germline.snp.readcounts",
        indel_readcount="{work_dir}/{tumor}/varscan2/germline.indel.readcounts"
    output:
        snps="{work_dir}/{tumor}/varscan2/germline.snps.fpfilter.vcf.gz",
        indels="{work_dir}/{tumor}/varscan2/germline.indels.fpfilter.vcf.gz"
    params:
        snp_out="{work_dir}/{tumor}/varscan2/germline.snps.fpfilter.vcf",
        indel_out="{work_dir}/{tumor}/varscan2/germline.indels.fpfilter.vcf"
    shell:
        """
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.snp_vcf} {input.snp_readcount} --output-file {params.snp_out}
        bgzip {params.snp_out}
        tabix -f -p vcf {output.snps}
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.indel_vcf} {input.indel_readcount} --output-file {params.indel_out}
        bgzip {params.indel_out}
        tabix -f -p vcf {output.indels}
        """

rule VarScan2_germline_merge:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/germline.snps.fpfilter.vcf.gz",
        indel_vcf="{work_dir}/{tumor}/varscan2/germline.indels.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.vcf.gz"
    shell:
        """
        bcftools concat -a {input} | bcftools sort -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_germline_normalized:
    input:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.norm.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_germline_clean:
    input:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/varscan2/sample.g.name",
        clean="{work_dir}/{tumor}/varscan2/germline.fpfilter.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/varscan2/temp.g.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

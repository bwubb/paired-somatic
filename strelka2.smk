##Author: Brad Wubbenhorst

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

rule run_strelka2:
    input: expand("data/final/{tumor}/{tumor}.strelka2.somatic.final.bcf",tumor=PAIRS.keys())
    #ADD MANTA

rule manta_write_workflow:
    input:
        unpack(paired_bams)
    output:
        temp("data/work/{tumor}/manta/runWorkflow.py")
    params:
        runDir="data/work/{tumor}/manta",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz']
    shell:
        "$HOME/software/manta/bin/configManta.py --normalBam {input.normal} --tumorBam {input.tumor} --referenceFasta {params.reference} --callRegions {params.bedgz} --exome --runDir {params.runDir}"

rule manta_main:
    input:
        "data/work/{tumor}/manta/runWorkflow.py"
    output:
        "data/work/{tumor}/manta/results/variants/candidateSmallIndels.vcf.gz",
        "data/work/{tumor}/manta/results/variants/candidateSV.vcf.gz",
        "data/work/{tumor}/manta/results/variants/diploidSV.vcf.gz",
        "data/work/{tumor}/manta/results/variants/somaticSV.vcf.gz"
    threads:
        4
    shell:
        "{input} -m local -j {threads}"

rule strelka2_write_workflow:
    input:
        unpack(paired_bams),
        indels="data/work/{tumor}/manta/results/variants/candidateSmallIndels.vcf.gz"
    output:
        temp("data/work/{tumor}/strelka2/runWorkflow.py")
    params:
        runDir="data/work/{tumor}/strelka2",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz']
    shell:
        "$HOME/software/strelka/bin/configureStrelkaSomaticWorkflow.py --normalBam {input.normal} --tumorBam {input.tumor} --indelCandidates {input.indels} --referenceFasta {params.reference} --callRegions {params.bedgz} --exome --runDir {params.runDir}"

rule strelka2_main:
    input:
        "data/work/{tumor}/strelka2/runWorkflow.py"
    output:
        snvs="data/work/{tumor}/strelka2/results/variants/somatic.snvs.vcf.gz",#not snps
        indels="data/work/{tumor}/strelka2/results/variants/somatic.indels.vcf.gz"
    params:
        workdir="data/work/{tumor}/strelka2",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz']
    threads:
        4
    shell:
        "{input} -m local -j {threads}"
        #can add -m, --max-mem <float>[kMG]    maximum memory to use [768M] for sort

rule strelka2_concat:
    input:
        snvs="data/work/{tumor}/strelka2/results/variants/somatic.snvs.vcf.gz",
        indels="data/work/{tumor}/strelka2/results/variants/somatic.indels.vcf.gz"
    output:
        "data/work/{tumor}/strelka2/somatic.vcf.gz"
    shell:
        """
        bcftools concat -a {input.snvs} {input.indels} | bcftools sort -W=tbi -Oz -o {output}
        """

rule strelka2_somatic_normalized:
    input:
        "data/work/{tumor}/strelka2/somatic.vcf.gz"
    output:
        norm="data/work/{tumor}/strelka2/somatic.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -W=tbi -Oz -o {output.norm}
        """

rule strelka2_sample_name:
    output:
        "data/work/{tumor}/strelka2/sample.name"
    params:
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        echo -e "TUMOR\\t{wildcards.tumor}\\nNORMAL\\t{params.normal}" > {output}
        echo -e "tumor\\t{wildcards.tumor}\\nnormal\\t{params.normal}" >> {output}
        """

rule strelka2_somatic_clean:
    input:
        name="data/work/{tumor}/strelka2/sample.name",
        vcf="data/work/{tumor}/strelka2/somatic.norm.vcf.gz"
    output:
        clean="data/work/{tumor}/strelka2/somatic.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("data/work/{tumor}/strelka2/temp.h.vcf.gz")
    shell:
        """
        bcftools reheader -f {params.fai} -s {input.name} -o {params.vcf} {input.vcf}
        bcftools index {params.vcf}

        bcftools view -s {wildcards.tumor},{params.normal} -e 'ALT="*"' -R {params.regions} {params.vcf} | \
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output.clean}
        """

rule strelka2_somatic_final:
    input:
        "data/work/{tumor}/strelka2/somatic.vcf.gz",
        "data/work/{tumor}/strelka2/somatic.norm.clean.vcf.gz"
    output:
        "data/final/{tumor}/{tumor}.strelka2.somatic.bcf",
        "data/final/{tumor}/{tumor}.strelka2.somatic.final.bcf"
    shell:
        """
        bcftools view -W=csi -Ob -o {output[0]} {input[0]}
        bcftools view -W=csi -Ob -o {output[1]} {input[1]}
        """

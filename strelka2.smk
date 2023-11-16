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

rule collect_strelka2:
    input:
        expand("data/work/{lib}/{tumor}/strelka2/somatic.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule write_Manta:
    input:
        unpack(paired_bams)
    output:
        temp("{work_dir}/{tumor}/manta/runWorkflow.py")
    params:
        runDir="{work_dir}/{tumor}/manta",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz']
    shell:
        "$HOME/software/manta/bin/configManta.py --normalBam {input.normal} --tumorBam {input.tumor} --referenceFasta {params.reference} --callRegions {params.bedgz} --exome --runDir {params.runDir}"

rule run_Manta:
    input:
        "{work_dir}/{tumor}/manta/runWorkflow.py"
    output:
        "{work_dir}/{tumor}/manta/results/variants/candidateSmallIndels.vcf.gz",
        "{work_dir}/{tumor}/manta/results/variants/candidateSV.vcf.gz",
        "{work_dir}/{tumor}/manta/results/variants/diploidSV.vcf.gz",
        "{work_dir}/{tumor}/manta/results/variants/somaticSV.vcf.gz"
    threads:
        4
    shell:
        "{input} -m local -j {threads}"

rule write_Strelka2:
    input:
        unpack(paired_bams),
        indels="{work_dir}/{tumor}/manta/results/variants/candidateSmallIndels.vcf.gz"
    output:
        temp("{work_dir}/{tumor}/strelka2/runWorkflow.py")
    params:
        runDir="{work_dir}/{tumor}/strelka2",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz']
    shell:
        "$HOME/software/strelka/bin/configureStrelkaSomaticWorkflow.py --normalBam {input.normal} --tumorBam {input.tumor} --indelCandidates {input.indels} --referenceFasta {params.reference} --callRegions {params.bedgz} --exome --runDir {params.runDir}"

rule run_Strelka2:#You sholuld put all concats in a differnt rule. Strelka wont work if the runWorkflow.py script already exists.
    input:
        "{work_dir}/{tumor}/strelka2/runWorkflow.py"
    output:
        snvs="{work_dir}/{tumor}/strelka2/results/variants/somatic.snvs.vcf.gz",#not snps
        indels="{work_dir}/{tumor}/strelka2/results/variants/somatic.indels.vcf.gz"
    params:
        workdir="{work_dir}/{tumor}/strelka2",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz']
    threads:
        4
    shell:
        "{input} -m local -j {threads}"
        #can add -m, --max-mem <float>[kMG]    maximum memory to use [768M] for sort

rule Strelka2_concat:
    input:
        snvs="{work_dir}/{tumor}/strelka2/results/variants/somatic.snvs.vcf.gz",
        indels="{work_dir}/{tumor}/strelka2/results/variants/somatic.indels.vcf.gz"
    output:
        "{work_dir}/{tumor}/strelka2/results/variants/somatic.vcf.gz"
    shell:
        """
        bcftools concat -a {input.snvs} {input.indels} | bcftools sort -O z -o {output}
        tabix -f -p vcf {output}
        """

rule Strelka2_somatic_normalized:
    input:
        "{work_dir}/{tumor}/strelka2/results/variants/somatic.vcf.gz"
    output:
        norm="{work_dir}/{tumor}/strelka2/somatic.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output.norm}
        tabix -f -p vcf {output.norm}
        """

rule Strelka2_somatic_clean:
    input:
        "{work_dir}/{tumor}/strelka2/somatic.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/strelka2/sample.name",
        clean="{work_dir}/{tumor}/strelka2/somatic.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/strelka2/temp.h.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

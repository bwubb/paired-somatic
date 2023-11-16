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

rule filter_applied_mutect2:
    input: expand("data/work/{lib}/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule unprocessed_mutect2:
    input: expand("data/work/{lib}/{tumor}/mutect2/somatic.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())
#rule three steps to create pon



rule run_Mutect2:
    input:
        unpack(paired_bams)
    output:#--bamOutput {output.bam}
        raw="{work_dir}/{tumor}/mutect2/somatic.vcf.gz",
        snps="{work_dir}/{tumor}/mutect2/somatic.snps.vcf.gz",
        indels="{work_dir}/{tumor}/mutect2/somatic.indels.vcf.gz",
        stats="{work_dir}/{tumor}/mutect2/somatic.vcf.gz.stats",
        f1r2="{work_dir}/{tumor}/mutect2/f1r2.tar.gz"
    params:
        ref=config['reference']['fasta'],
        intervals=config['resources']['targets_intervals'],
        tumor=lambda wildcards: wildcards.tumor,
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        memory='16g'
    shell:
        """
        gatk --java-options '-Xmx{params.memory}' Mutect2 -R {params.ref} -I {input.tumor} -I {input.normal} -tumor {params.tumor} -normal {params.normal} -L {params.intervals} -O {output.raw} --f1r2-tar-gz {output.f1r2} --genotype-germline-sites true --genotype-pon-sites tru
        gatk --java-options '-Xmx{params.memory}' SelectVariants -R {params.ref} -V {output.raw} -O {output.snps} -L {params.intervals} -select-type SNP
        gatk --java-options '-Xmx{params.memory}' SelectVariants -R {params.ref} -V {output.raw} -O {output.indels} -L {params.intervals} -select-type INDEL
        """
#If GRCh38 consider adding --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter, see documentation for details
#        -germline-resource af-only-gnomad.vcf \
#        -pon panel_of_normals.vcf   \

rule Mutect2_LearnReadOrientationModel:
    input:
        "{work_dir}/{tumor}/mutect2/f1r2.tar.gz"
    output:
        "{work_dir}/{tumor}/mutect2/read-orientation-model.tar.gz"
    shell:
        "gatk LearnReadOrientationModel -I {input} -O {output}"

rule Mutect2_GetPileupSummaries:
    input:
        unpack(paired_bams)
    output:
        pileup="{work_dir}/{tumor}/mutect2/getpileupsummaries.table"
    params:
        #allele=f"$HOME/resources/Vcf_files/gnomad.exomes.r2.1.1.sites.{config['reference']['key']}.{config['resources']['targets_key']}.common_biallelic_snps.simplified.vcf.gz",
        allele=config['resources']['common_snps'],
        intervals=config['resources']['targets_intervals']
    shell:
        "gatk GetPileupSummaries -I {input.tumor} -V {params.allele} -L {params.intervals} -O {output.pileup}"

rule Mutect2_CalculateContamination:
    input:
        pileup="{work_dir}/{tumor}/mutect2/getpileupsummaries.table"
    output:
        contamination="{work_dir}/{tumor}/mutect2/calculatecontamination.table",
        segments="{work_dir}/{tumor}/mutect2/segments.table"
    shell:
        "gatk CalculateContamination -I {input.pileup} -tumor-segmentation {output.segments} -O {output.contamination}"

rule Mutect2_FilterMutectCalls:
    input:
        vcf="{work_dir}/{tumor}/mutect2/somatic.vcf.gz",
        stats="{work_dir}/{tumor}/mutect2/somatic.vcf.gz.stats",
        contamination="{work_dir}/{tumor}/mutect2/calculatecontamination.table",
        segments="{work_dir}/{tumor}/mutect2/segments.table",
        model="{work_dir}/{tumor}/mutect2/read-orientation-model.tar.gz"
    output:
        "{work_dir}/{tumor}/mutect2/somatic.filtered.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        "gatk FilterMutectCalls -R {params.ref} -V {input.vcf} --tumor-segmentation {input.segments} --contamination-table {input.contamination} --ob-priors {input.model} -O {output}"

rule Mutect2_somatic_normalized:
    input:
        "{work_dir}/{tumor}/mutect2/somatic.filtered.vcf.gz"
    output:
        norm="{work_dir}/{tumor}/mutect2/somatic.filtered.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output.norm}
        tabix -f -p vcf {output.norm}
        """

rule Mutect2_somatic_clean:
    input:
        "{work_dir}/{tumor}/mutect2/somatic.filtered.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/mutect2/sample.name",
        clean="{work_dir}/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/mutect2/temp.h.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

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

rule run_mutect2:
    input: expand("data/final/{tumor}/{tumor}.mutect2.somatic.final.bcf",tumor=PAIRS.keys())

rule filter_applied_mutect2:
    input: expand("data/work/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz",tumor=PAIRS.keys())

rule unprocessed_mutect2:
    input: expand("data/work/{tumor}/mutect2/somatic.vcf.gz",tumor=PAIRS.keys())

rule vep_mutect2:
    input: expand("data/work/{tumor}/mutect2/somatic.filtered.norm.clean.vep.vcf",tumor=PAIRS.keys())
#Attempting the pon_db is a waste of time.
#Its takes FOREVER, a lot of disk space, and then doesn't work

rule mutect2_main:
    input:
        unpack(paired_bams)
    output:
        raw="data/work/{tumor}/mutect2/somatic.vcf.gz",
        snps="data/work/{tumor}/mutect2/somatic.snps.vcf.gz",
        indels="data/work/{tumor}/mutect2/somatic.indels.vcf.gz",
        stats="data/work/{tumor}/mutect2/somatic.vcf.gz.stats",
        f1r2="data/work/{tumor}/mutect2/f1r2.tar.gz"
    params:
        ref=config['reference']['fasta'],
        intervals=config['resources']['targets_intervals'],
        tumor=lambda wildcards: wildcards.tumor,
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        memory='16g'
    shell:
        """
        gatk --java-options '-Xmx{params.memory}' Mutect2 -R {params.ref} -I {input.tumor} -I {input.normal} -tumor {params.tumor} -normal {params.normal} -L {params.intervals} -O {output.raw} --f1r2-tar-gz {output.f1r2}
        gatk --java-options '-Xmx{params.memory}' SelectVariants -R {params.ref} -V {output.raw} -O {output.snps} -L {params.intervals} -select-type SNP
        gatk --java-options '-Xmx{params.memory}' SelectVariants -R {params.ref} -V {output.raw} -O {output.indels} -L {params.intervals} -select-type INDEL
        """
#--genotype-germline-sites true
#--genotype-pon-sites true

rule mutect2_LearnReadOrientationModel:
    input:
        "data/work/{tumor}/mutect2/f1r2.tar.gz"
    output:
        "data/work/{tumor}/mutect2/read-orientation-model.tar.gz"
    shell:
        """
        gatk LearnReadOrientationModel -I {input} -O {output}
        """

rule mutect2_GetPileupSummaries:
    input:
        unpack(paired_bams)
    output:
        pileup="data/work/{tumor}/mutect2/getpileupsummaries.table"
    params:
        allele=config['resources']['common_snps'],
        intervals=config['resources']['targets_intervals']
    shell:
        """
        gatk GetPileupSummaries -I {input.tumor} -V {params.allele} -L {params.intervals} -O {output.pileup}
        """

rule mutect2_CalculateContamination:
    input:
        pileup="data/work/{tumor}/mutect2/getpileupsummaries.table"
    output:
        contamination="data/work/{tumor}/mutect2/calculatecontamination.table",
        segments="data/work/{tumor}/mutect2/segments.table"
    shell:
        """
        gatk CalculateContamination -I {input.pileup} -tumor-segmentation {output.segments} -O {output.contamination}
        """

rule mutect2_FilterMutectCalls:
    input:
        vcf="data/work/{tumor}/mutect2/somatic.vcf.gz",
        stats="data/work/{tumor}/mutect2/somatic.vcf.gz.stats",
        contamination="data/work/{tumor}/mutect2/calculatecontamination.table",
        segments="data/work/{tumor}/mutect2/segments.table",
        model="data/work/{tumor}/mutect2/read-orientation-model.tar.gz"
    output:
        "data/work/{tumor}/mutect2/somatic.filtered.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        gatk FilterMutectCalls -R {params.ref} -V {input.vcf} --tumor-segmentation {input.segments} --contamination-table {input.contamination} --ob-priors {input.model} -O {output}
        """

rule mutect2_somatic_normalized:
    input:
        "data/work/{tumor}/mutect2/somatic.filtered.vcf.gz"
    output:
        norm="data/work/{tumor}/mutect2/somatic.filtered.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -W=tbi -Oz -o {output.norm}
        """

rule mutect2_sample_name:
    output:
        "data/work/{tumor}/mutect2/sample.name"
    params:
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        echo -e "TUMOR\\t{wildcards.tumor}\\nNORMAL\\t{params.normal}" > {output}
        echo -e "tumor\\t{wildcards.tumor}\\nnormal\\t{params.normal}" >> {output}
        """

rule mutect2_somatic_clean:
    input:
        name="data/work/{tumor}/mutect2/sample.name",
        vcf="data/work/{tumor}/mutect2/somatic.filtered.norm.vcf.gz"
    output:
        clean="data/work/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("data/work/{tumor}/mutect2/temp.h.vcf.gz")
    shell:
        """
        bcftools reheader -f {params.fai} -s {input.name} -o {params.vcf} {input.vcf}
        bcftools index {params.vcf}

        bcftools view -s {wildcards.tumor},{params.normal} -e 'ALT="*"' -R {params.regions} {params.vcf} | \
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output.clean}
        """

rule mutect2_somatic_final:
    input:
        "data/work/{tumor}/mutect2/somatic.vcf.gz",
        "data/work/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    output:
        "data/final/{tumor}/{tumor}.mutect2.somatic.bcf",
        "data/final/{tumor}/{tumor}.mutect2.somatic.final.bcf"
    shell:
        """
        bcftools view -W=csi -Ob -o {output[0]} {input[0]}
        bcftools view -W=csi -Ob -o {output[1]} {input[1]}
        """

rule mutect2_somatic_vep:
    input:
        "data/work/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    output:
        "data/work/{tumor}/mutect2/somatic.filtered.norm.clean.vep.vcf"
    shell:
        """
        singularity run -H $PWD:/home \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i {input} \
        -o {output} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf --everything --canonical \
        --assembly GRCh38 \
        --species homo_sapiens \
        --fasta /opt/vep/resources/Genomes/Human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin REVEL,/opt/vep/.vep/revel/revel_grch38.tsv.gz \
        --plugin SpliceAI,snv=/opt/vep/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz,indel=/opt/vep/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz \
        --plugin gnomADc,/opt/vep/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz \
        --plugin UTRAnnotator,/opt/vep/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt \
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar_20250106.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
        """
        
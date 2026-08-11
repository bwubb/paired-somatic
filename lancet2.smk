##Author: Brad Wubbenhorst

### INIT ###

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

### FUNCTIONS ###

def paired_bams(wildcards):
    ref=config['reference']['key']
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

### SNAKEMAKE ###

localrules: lancet2_sample_name

rule run_lancet2:
    input: expand("data/final/{tumor}/{tumor}.lancet2.somatic.final.bcf",tumor=PAIRS.keys())

rule processed_lancet2:
    input: expand("data/work/{tumor}/lancet2/somatic.norm.clean.vcf.gz",tumor=PAIRS.keys())

rule unprocessed_lancet2:
    input: expand("data/work/{tumor}/lancet2/somatic.vcf.gz",tumor=PAIRS.keys())

# IT shipped Lancet2 as a SIF (no local build).
# Bind /project + /home/bwubb (projects often symlink into /project). Never --pwd.
# /scratch exists on compute nodes (submitted jobs), not on cybertron2 head node.
LANCET2_SIF=config.get('lancet2',{}).get('container','/appl/containers/lancet_v2.8.7.sif')

rule lancet2_main:
    input:
        unpack(paired_bams)
    output:
        "data/work/{tumor}/lancet2/somatic.vcf.gz"
    params:
        ref=config['reference']['fasta'],
        region=config['resources']['targets_bed'],
        sif=LANCET2_SIF
    threads: 4
    resources:
        mem_mb=32768
    shell:
        """
        WORKDIR="$(pwd -P)"
        OUT_VCF="$WORKDIR/{output}"
        mkdir -p "$(dirname "$OUT_VCF")"

        singularity exec \
          --bind /project:/project \
          --bind /home/bwubb:/home/bwubb \
          --bind /scratch:/scratch \
          --bind "$WORKDIR:$WORKDIR" \
          {params.sif} \
          /usr/bin/Lancet2 pipeline \
            --tumor {input.tumor} \
            --normal {input.normal} \
            --reference {params.ref} \
            --region {params.region} \
            --num-threads {threads} \
            --out-vcfgz "$OUT_VCF"
        """

rule lancet2_somatic_normalized:
    input:
        "data/work/{tumor}/lancet2/somatic.vcf.gz"
    output:
        norm="data/work/{tumor}/lancet2/somatic.norm.vcf.gz"
    resources:
        mem_mb=6144
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -W=tbi -Oz -o {output.norm}
        """

rule lancet2_sample_name:
    output:
        "data/work/{tumor}/lancet2/sample.name"
    params:
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        echo -e "CASE\\t{wildcards.tumor}\\nCTRL\\t{params.normal}" > {output}
        """

rule lancet2_somatic_clean:
    input:
        name="data/work/{tumor}/lancet2/sample.name",
        vcf="data/work/{tumor}/lancet2/somatic.norm.vcf.gz"
    output:
        clean="data/work/{tumor}/lancet2/somatic.norm.clean.vcf.gz"
    resources:
        mem_mb=6144
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("data/work/{tumor}/lancet2/temp.h.vcf.gz")
    shell:
        """
        bcftools reheader -f {params.fai} -s {input.name} -o {params.vcf} {input.vcf}
        bcftools index {params.vcf}

        bcftools view -s {wildcards.tumor},{params.normal} -e 'ALT="*"' -R {params.regions} {params.vcf} | \
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output.clean}
        """

rule lancet2_somatic_final:
    input:
        "data/work/{tumor}/lancet2/somatic.vcf.gz",
        "data/work/{tumor}/lancet2/somatic.norm.clean.vcf.gz"
    output:
        "data/final/{tumor}/{tumor}.lancet2.somatic.bcf",
        "data/final/{tumor}/{tumor}.lancet2.somatic.final.bcf"
    resources:
        mem_mb=6144
    shell:
        """
        bcftools view -W=csi -Ob -o {output[0]} {input[0]}
        bcftools view -W=csi -Ob -o {output[1]} {input[1]}
        """

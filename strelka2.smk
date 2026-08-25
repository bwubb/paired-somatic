##Author: Brad Wubbenhorst

### INIT ###

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

# Illumina Manta + Strelka2 are archived and locked to Python 2. RHEL9 has no
# python2 module, so host $HOME/software/.../config*.py fails with
# "python2: No such file or directory". Run configure AND the generated
# runWorkflow.py via Singularity BioContainers images that bundle python2.
# Pull examples:
#   singularity pull manta_1.6.0.sif \
#     docker://quay.io/biocontainers/manta:1.6.0--py27h9948957_6
#   singularity pull strelka_2.9.10.sif \
#     docker://quay.io/biocontainers/strelka:2.9.10--0
MANTA_SIF=config.get('strelka2',{}).get('manta_container','/home/bwubb/containers/manta_1.6.0.sif')
STRELKA_SIF=config.get('strelka2',{}).get('strelka_container','/home/bwubb/containers/strelka_2.9.10.sif')

# Same bind pattern as lancet2 / cnvkit_singularity: physical cwd under /project,
# never --pwd (symlink /home/.../projects breaks chdir).
def singularity_exec(sif):
    return (
        'WORKDIR="$(pwd -P)"; cd "$WORKDIR"; '
        'singularity exec '
        '--bind /project:/project '
        '--bind /home/bwubb:/home/bwubb '
        '--bind /scratch:/scratch '
        '--bind "$WORKDIR:$WORKDIR" '
        f'{sif}'
    )

### FUNCTIONS ###

def paired_bams(wildcards):
    ref=config['reference']['key']
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

### SNAKEMAKE ###

localrules: strelka2_sample_name

rule run_strelka2:
    input: expand("data/final/{tumor}/{tumor}.strelka2.somatic.final.bcf",tumor=PAIRS.keys())

rule manta_write_workflow:
    input:
        unpack(paired_bams)
    output:
        temp("data/work/{tumor}/manta/runWorkflow.py")
    resources:
        mem_mb=6144
    params:
        runDir="data/work/{tumor}/manta",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz'],
        exec=singularity_exec(MANTA_SIF)
    shell:
        """
        {params.exec} configManta.py \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --referenceFasta {params.reference} \
            --callRegions {params.bedgz} \
            --exome \
            --runDir {params.runDir}
        """

rule manta_main:
    input:
        "data/work/{tumor}/manta/runWorkflow.py"
    output:
        "data/work/{tumor}/manta/results/variants/candidateSmallIndels.vcf.gz",
        "data/work/{tumor}/manta/results/variants/candidateSV.vcf.gz",
        "data/work/{tumor}/manta/results/variants/diploidSV.vcf.gz",
        "data/work/{tumor}/manta/results/variants/somaticSV.vcf.gz"
    threads: 4
    resources:
        mem_mb=18432
    params:
        exec=singularity_exec(MANTA_SIF)
    shell:
        """
        {params.exec} {input} -m local -j {threads}
        """

rule strelka2_write_workflow:
    input:
        unpack(paired_bams),
        indels="data/work/{tumor}/manta/results/variants/candidateSmallIndels.vcf.gz"
    output:
        temp("data/work/{tumor}/strelka2/runWorkflow.py")
    resources:
        mem_mb=6144
    params:
        runDir="data/work/{tumor}/strelka2",
        reference=config['reference']['fasta'],
        bedgz=config['resources']['targets_bedgz'],
        exec=singularity_exec(STRELKA_SIF)
    shell:
        """
        {params.exec} configureStrelkaSomaticWorkflow.py \
            --normalBam {input.normal} \
            --tumorBam {input.tumor} \
            --indelCandidates {input.indels} \
            --referenceFasta {params.reference} \
            --callRegions {params.bedgz} \
            --exome \
            --runDir {params.runDir}
        """

rule strelka2_main:
    input:
        "data/work/{tumor}/strelka2/runWorkflow.py"
    output:
        snvs="data/work/{tumor}/strelka2/results/variants/somatic.snvs.vcf.gz",
        indels="data/work/{tumor}/strelka2/results/variants/somatic.indels.vcf.gz"
    threads: 4
    resources:
        mem_mb=18432
    params:
        exec=singularity_exec(STRELKA_SIF)
    shell:
        """
        {params.exec} {input} -m local -j {threads}
        """

rule strelka2_concat:
    input:
        snvs="data/work/{tumor}/strelka2/results/variants/somatic.snvs.vcf.gz",
        indels="data/work/{tumor}/strelka2/results/variants/somatic.indels.vcf.gz"
    output:
        "data/work/{tumor}/strelka2/somatic.vcf.gz"
    resources:
        mem_mb=6144
    shell:
        """
        bcftools concat -a {input.snvs} {input.indels} | bcftools sort -W=tbi -Oz -o {output}
        """

rule strelka2_somatic_normalized:
    input:
        "data/work/{tumor}/strelka2/somatic.vcf.gz"
    output:
        norm="data/work/{tumor}/strelka2/somatic.norm.vcf.gz"
    resources:
        mem_mb=6144
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
    resources:
        mem_mb=6144
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
        bcftools annotate --set-id '%CHROM_%POS_%REF_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output.clean}
        """

rule strelka2_somatic_final:
    input:
        "data/work/{tumor}/strelka2/somatic.vcf.gz",
        "data/work/{tumor}/strelka2/somatic.norm.clean.vcf.gz"
    output:
        "data/final/{tumor}/{tumor}.strelka2.somatic.bcf",
        "data/final/{tumor}/{tumor}.strelka2.somatic.final.bcf"
    resources:
        mem_mb=6144
    shell:
        """
        bcftools view -W=csi -Ob -o {output[0]} {input[0]}
        bcftools view -W=csi -Ob -o {output[1]} {input[1]}
        """

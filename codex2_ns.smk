import errno
import os

with open(config['project']['sample_list'],'r') as i:
    SAMPLES=i.read().splitlines()

with open(config['project']['bam_list'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

#ALL Codex R code needs a total rewrite.
#Get GC and mappability once.
#Multiplex get coverage for each sample
#merge data.
#Run analysis by Chr



rule all:
    input:
        f"project_data/work/{config['resources'][targets_key]}/codex2/coverageQC.csv"

rule CODEX2_CoverageQC:
    input:
        [BAMS[sample] for sample in SAMPLES]
    output:
        Y_qc="project_data/work/{lib}/codex2/coverageQC.csv",
        gc_qc="project_data/work/{lib}/codex2/gc_qc.csv",
        N="project_data/work/{lib}/codex2/library_size_factor.csv"
    params:
        project="project_data/work/{lib}/codex2",
        bed=config['resources']['targets_bed']
    script:
        "codex2_paired_getCoverage.snakemake.R"

rule CODEX2_ns:
    input:
        Y_qc="project_data/work/{lib}/codex2/coverageQC.csv",
        gc_qc="project_data/work/{lib}/codex2/gc_qc.csv",
        N="project_data/work/{lib}/codex2/library_size_factor.csv"
    output:
        "project_data/work/{lib}/codex2/chr{chr}.codex2.filtered_segments.txt"
    params:
        project="project_data/work/{lib}/codex2",
        chr=lambda wildcards: wildcards.chr,
        normals=config['project']['normal_list']
    script:
        "codex2_paired_runChr.snakemake.R"

#rule CODEX2_merge:
#    input:
#    output:
#    shell:
    
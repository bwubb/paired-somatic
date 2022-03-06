import os
#include ./lancet.smk
#include ./mutect2.smk
#include ./strelka2.smk
#include ./vardictjava.smk
#include ./varscan2.smk

##INIT
with open(config.get('project',{}).get('sample_list','sample.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

TUMORS=PAIRS.keys()

##PYTHON
def map_vcf(wildcards):
    V={'lancet':f'data/final/{wildcards.lib}/{wildcards.tumor}/lancet/{wildcards.tumor}.somatic.norm.clean.vcf.gz',
    'mutect2':f'data/final/{wildcards.lib}/{wildcards.tumor}/mutect2/{wildcards.tumor}.somatic.filtered.norm.clean.vcf.gz',
    'strelka2':f'data/final/{wildcards.lib}/{wildcards.tumor}/strelka2/{wildcards.tumor}.somatic.norm.clean.vcf.gz',
    'vardict':f'data/final/{wildcards.lib}/{wildcards.tumor}/vardict/{wildcards.tumor}.somatic.twice_filtered.norm.clean.vcf.gz',
    'varscan2':f'data/final/{wildcards.lib}/{wildcards.tumor}/varscan2/{wildcards.tumor}.somatic.fpfilter.norm.clean.vcf.gz'}
    return V[wildcards.caller]

def map_preprocess(wildcards):
    return {'bam':BAMS[wildcards.sample],'bcf':f'data/work/{wildcards.lib}/{wildcards.tumor}/varlociraptor/candidate.bcf','aln':f'data/work/{wildcards.lib}/{wildcards.tumor}/varlociraptor/{wildcards.sample}.alignment-properties.json'}

def map_varlociraptor_paired(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':f'data/work/{wildcards.lib}/{tumor}/varlociraptor/{tumor}.observations.bcf','normal':f'data/work/{wildcards.lib}/{tumor}/varlociraptor/{normal}.observations.bcf'}

def genome_size(wildcards):
    G={'S04380110':'5.0e7','S07604715':'6.6e7','S31285117':'4.9e7','xgen-exome-research-panel-targets-grch37':'3.9e7'}
    return G[config['resources']['targets_key']]

##TARGET RULES
rule local_fdr:
    input:
        expand("data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.bcf",lib=config['resources']['targets_key'],tumor=TUMORS)

rule collect_vep:
    input:
        expand("data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.vep.vcf.gz",lib=config['resources']['targets_key'],tumor=TUMORS)

##SNAKEMAKE
rule candidate_tsv:
    input:
        map_vcf
    output:
        "data/work/{lib}/{tumor}/{caller}/cadidate.tsv"
    shell:
        """
        bcftools query -i 'FILTER=\"PASS\"' -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\t.\n' {input} > {output}
        """
    #In the future a slower python version will be needed to make a standardized format INFO field.
    #Drop GTs still though

rule candidate_bcf:
    input:
        "data/work/{lib}/{tumor}/lancet/cadidate.tsv",
        "data/work/{lib}/{tumor}/mutect2/cadidate.tsv",
        "data/work/{lib}/{tumor}/strelka2/cadidate.tsv",
        "data/work/{lib}/{tumor}/vardict/cadidate.tsv",
        "data/work/{lib}/{tumor}/varscan2/cadidate.tsv"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/candidate.bcf"
    params:
        tsv=temp("data/work/{lib}/{tumor}/varlociraptor/candidate.tsv"),
        header="/home/bwubb/resources/Vcf_files/simplified-header-w_contig.vcf"
    shell:
        """
        cat {input} | sort | uniq > {params.tsv}
        cat {params.header} {params.tsv} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

#if {{ conda env list | grep 'varlociraptor'; }} >/dev/null 2>&1; then source activate varlociraptor; fi
rule estimate_properties:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{sample}.alignment-properties.json"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        varlociraptor estimate alignment-properties {params.ref} --bam {input.bam} > {output}
        """

rule preprocess_sample:
    input:
        unpack(map_preprocess)
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{sample}.observations.bcf"
    params:
        ref='/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta'
    shell:
        """
        varlociraptor preprocess variants {params.ref} --alignment-properties {input.aln} --bam {input.bam} --candidates {input.bcf} > {output}
        """

rule varlociraptor_paired:
    input:
        unpack(map_varlociraptor_paired)
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.bcf"
    params:
        purity='0.75'
    shell:
        """
        varlociraptor call variants tumor-normal --purity {params.purity} --tumor {input.tumor} --normal {input.normal} > {output}
        bcftools index {output}
        """

rule varlociraptor_local_fdr:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.bcf"
    params:
        fdr="0.05"
    shell:
        """
        varlociraptor filter-calls control-fdr --local {input} --events SOMATIC_TUMOR --fdr {params.fdr} > {output}
        bcftools index {output}
        """

rule varlociraptor_fdr_snv:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.snv.bcf"
    params:
        fdr='0.05'
    shell:
        """
        varlociraptor filter-calls control-fdr {input} --events SOMATIC_TUMOR --fdr {params.fdr} --var SNV > {output}
        bcftools index {output}
        """

rule varlociraptor_fdr_mnv:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.mnv.bcf"
    params:
        fdr='0.03'
    shell:
        """
        varlociraptor filter-calls control-fdr {input} --events SOMATIC_TUMOR --fdr {params.fdr} --var MNV > {output}
        bcftools index {output}
        """

rule varlociraptor_fdr_ins:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.ins.bcf"
    params:
        fdr='0.03'
    shell:
        """
        varlociraptor filter-calls control-fdr {input} --events SOMATIC_TUMOR --fdr {params.fdr} --var INS > {output}
        bcftools index {output}
        """

rule varlociraptor_fdr_del:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.del.bcf"
    params:
        fdr='0.03'
    shell:
        """
        varlociraptor filter-calls control-fdr {input} --events SOMATIC_TUMOR --fdr {params.fdr} --var DEL > {output}
        bcftools index {output}
        """

rule merge_fdr_bcf:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.snv.bcf",
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.mnv.bcf",
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.ins.bcf",
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.fdr-filtered.del.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.global-fdr.bcf"
    shell:
        """
        bcftools concat -a {input} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

rule varlociraptor_vep:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.bcf"
    output:
        vcf="data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.vep.vcf.gz",
        bcf="data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.vep.bcf"
    params:
        in_vcf='data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.vcf.gz',
        out_vcf='data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.local-fdr.vep.vcf',
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        ref_fa="data/work/{lib}/{tumor}/varlociraptor/reference.fa",
        mut_fa="data/work/{lib}/{tumor}/varlociraptor/mutated.fa"
    shell:
        """
        bcftools view -O v -o {params.in_vcf} {input} && tabix -fp vcf {params.in_vcf}

        if {{ conda env list | grep 'vep'; }} >/dev/null 2>&1; then source activate vep; fi

        vep -i {params.in_vcf} -o {params.out_vcf} \
        --force_overwrite \
        --offline \
        --cache \
        --format vcf \
        --vcf \
        --everything \
        --canonical \
        --assembly {params.assembly} \
        --species homo_sapiens \
        --fasta {params.fa} \
        --plugin NMD \
        --plugin ProteinSeqs,{params.ref_fa},{params.mut_fa} \
        --plugin Downstream \
        --plugin REVEL,{params.revel} \
        --plugin SpliceAI,snv={params.splice_snv},indel={params.splice_indel} \
        --plugin gnomADc,{params.gnomAD} \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

        bgzip {params.out_vcf} && tabix -fp vcf {output.vcf}

        bcftools view -O b -o {output.bcf} {output.vcf}
        bcftools index -f {output.bcf}
        """
        #--plugin LoF,loftee_path:{params.loftee}
        #--plugin UTRannotator \
        #move ref mut fa to output after testing is complete.
        #possibly bgzip/index output fasta file

rule varlociraptor_estimate_mutational_burden:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.posterior-odds.strong.vep.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.posterior-odds.strong.vep.tmb.json",
        "data/work/{lib}/{tumor}/varlociraptor/{tumor}.somatic.posterior-odds.strong.vep.tmb.svg"
    params:
        size=genome_size
        #--vaf-cutoff <cutoff>
    shell:
        """
        varlociraptor estimate mutational-burden --plot-mode curve --coding-genome-size {params.size} --sample tumor --events SOMATIC_TUMOR < {input} > {output[0]}
        vl2svg {output[0]} {output[1]}
        """
rule bcftools_isec:
    input:
        lancet="data/work/{lib}/{tumor}/lancet/somatic.norm.clean.vcf.gz",
        mutect2="data/work/{lib}/{tumor}/mutect2/somatic.filtered.norm.clean.std.vcf.gz",
        strelka2="data/work/{lib}/{tumor}/strelka2/somatic.norm.clean.std.vcf.gz",
        vardict="data/work/{lib}/{tumor}/vardict/somatic.twice_filtered.norm.clean.std.vcf.gz",
        varscan2="data/work/{lib}/{tumor}/varscan2/somatic.fpfilter.norm.clean.std.vcf.gz"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/sites.txt"
    shell:
        """
        bcftools isec {input.lancet} {input.mutect2} {input.strelka2} {input.vardict} {input.varscan2} > {output}
        """

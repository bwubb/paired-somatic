import os

#include:"./mutect2.snake"
#include:"./strelka2.snake"
#include:"./vardictjava.snake"
#include:"./varscan2.snake"

def annovar_samples(wildcards):
    return [wildcards.tumor,PAIRS[wildcards.tumor]]

### INIT ###

with open(config['project']['pair_table'],'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

TUMORS=PAIRS.keys()
#clean up rule
#perform du first to see all files and size removed.

### SNAKEMAKE

rule all_ensemble:
    input:
        expand("data/work/{lib}/{tumor}/annovar/{variant}.ensemble.hg19_multianno.report.tsv",variant=['somatic','germline-loh'],lib=config['resources']['targets_key'],tumor=TUMORS)

rule somatic_ensemble:
    input:
        expand("data/work/{lib}/{tumor}/annovar/somatic.ensemble.hg19_multianno.report.tsv",lib=config['resources']['targets_key'],tumor=TUMORS)

rule germline_ensemble:
    input:
        expand("data/work/{lib}/{tumor}/annovar/germline-loh.ensemble.hg19_multianno.report.tsv",lib=config['resources']['targets_key'],tumor=TUMORS)

rule Bcftools_somatic_isec:
    input:
        mutect2="data/work/{lib}/{tumor}/mutect2/somatic.filtered.norm.clean.std.vcf.gz",
        strelka2="data/work/{lib}/{tumor}/strelka2/somatic.norm.clean.std.vcf.gz",
        vardict="data/work/{lib}/{tumor}/vardict/somatic.twice_filtered.norm.clean.std.vcf.gz",
        varscan2="data/work/{lib}/{tumor}/varscan2/somatic.fpfilter.norm.clean.std.vcf.gz"
    params:
        outdir="data/work/{lib}/{tumor}/bcftools/somatic"
    output:
        "data/work/{lib}/{tumor}/bcftools/somatic/sites.txt"
    shell:
        "bcftools isec -n+2 -p {params.outdir} {input.mutect2} {input.strelka2} {input.vardict} {input.varscan2}"

rule concordant_somatic_calls:
    input:
        sites="data/work/{lib}/{tumor}/bcftools/somatic/sites.txt",
        mutect2="data/work/{lib}/{tumor}/mutect2/somatic.filtered.norm.clean.std.vcf.gz",
        strelka2="data/work/{lib}/{tumor}/strelka2/somatic.norm.clean.std.vcf.gz",
        vardict="data/work/{lib}/{tumor}/vardict/somatic.twice_filtered.norm.clean.std.vcf.gz",
        varscan2="data/work/{lib}/{tumor}/varscan2/somatic.fpfilter.norm.clean.std.vcf.gz"
    output:
        "data/work/{lib}/{tumor}/bcftools/somatic/somatic.ensemble.vcf.gz"
    params:
        tumor=lambda wildcards: wildcards.tumor,
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        lib=config['resources']['targets_key'],
        outfile="data/work/{lib}/{tumor}/bcftools/somatic/somatic.ensemble.vcf"
    shell:
        """
        python tumor-normal.ensemble-vcf.py -b {input.sites} -o {params.outfile} --mutect2_vcf {input.mutect2} --strelka2_vcf {input.strelka2} --vardict_vcf {input.vardict} --varscan2_vcf {input.varscan2} -T {params.tumor} -N {params.normal} -L {params.lib}
        """

rule ANNOVAR_somatic_ensemble_vcf:
    input:
        'data/work/{lib}/{tumor}/bcftools/somatic/somatic.ensemble.vcf.gz'
    output:
        'data/work/{lib}/{tumor}/annovar/somatic.ensemble.hg19_multianno.vcf.gz'
    params:
        humandb='$HOME/resources/annovar/humandb/',
        output_p='data/work/{lib}/{tumor}/annovar/somatic.ensemble'
    shell:
        """
        tabix -fp vcf {input}
        table_annovar.pl {input} {params.humandb} --buildver hg19 --vcfinput --outfile {params.output_p} --protocol refGene,cytoband,genomicSuperDups,dbscsnv11,avsnp150,dbnsfp35a,mcap,revel,popfreq_max_20150413,exac03,exac03nontcga,gnomad211_exome,intervar_20180118,icgc21,cosmic84_coding,cosmic84_noncoding,clinvar_20190305 --operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -remove
        bgzip {params.output_p}.hg19_multianno.vcf
        tabix -fp vcf {params.output_p}.hg19_multianno.vcf.gz
        """

rule run_somatic_annovartools:
    input:
        'data/work/{lib}/{tumor}/annovar/somatic.ensemble.hg19_multianno.vcf.gz'
    output:
        'data/work/{lib}/{tumor}/annovar/somatic.ensemble.hg19_multianno.report.tsv'
    params:
        header='/home/bwubb/resources/annovar/ensemble.somatic-annotation-header.20200827.txt',
        mode='paired',
        samples=annovar_samples
    shell:
        'python annovartools.v03.py -I {input} -O {output} --header {params.header} -m {params.mode} {params.samples}'


rule Bcftools_germline_loh_concat:
    input:
        vardict_germline="data/work/{lib}/{tumor}/vardict/germline.twice_filtered.norm.clean.std.vcf.gz",
        vardict_loh="data/work/{lib}/{tumor}/vardict/loh.twice_filtered.norm.clean.std.vcf.gz",
        varscan2_germline="data/work/{lib}/{tumor}/varscan2/germline.fpfilter.norm.clean.std.vcf.gz",
        varscan2_loh="data/work/{lib}/{tumor}/varscan2/loh.fpfilter.norm.clean.std.vcf.gz"
    output:
        vardict_out="data/work/{lib}/{tumor}/vardict/germline-loh.twice_filtered.norm.clean.std.vcf.gz",
        varscan2_out="data/work/{lib}/{tumor}/varscan2/germline-loh.fpfilter.norm.clean.std.vcf.gz"
    params:
    shell:
        """
        bcftools concat -a {input.vardict_germline} {input.vardict_loh} | bcftools sort -O z -o {output.vardict_out}
        tabix -f -p vcf {output.vardict_out}
        bcftools concat -a {input.varscan2_germline} {input.varscan2_loh} | bcftools sort -O z -o {output.varscan2_out}
        tabix -f -p vcf {output.varscan2_out}
        """

rule Bcftools_germline_loh_isec:
    input:
        vardict="data/work/{lib}/{tumor}/vardict/germline-loh.twice_filtered.norm.clean.std.vcf.gz",
        varscan2="data/work/{lib}/{tumor}/varscan2/germline-loh.fpfilter.norm.clean.std.vcf.gz"
    output:
        "data/work/{lib}/{tumor}/bcftools/germline-loh/sites.txt"
    params:
        outdir="data/work/{lib}/{tumor}/bcftools/germline-loh"
    shell:
        """
        tabix -f -p vcf {input.vardict}
        tabix -f -p vcf {input.varscan2}
        bcftools isec -n+1 -p {params.outdir} {input.vardict} {input.varscan2}
        """
        #Additionally created vcfs are private calls.
        #Can be renamed and used.

rule concordant_germline_loh_calls:
    input:
        sites="data/work/{lib}/{tumor}/bcftools/germline-loh/sites.txt",
        vardict="data/work/{lib}/{tumor}/vardict/germline-loh.twice_filtered.norm.clean.std.vcf.gz",
        varscan2="data/work/{lib}/{tumor}/varscan2/germline-loh.fpfilter.norm.clean.std.vcf.gz"
    output:
        'data/work/{lib}/{tumor}/bcftools/germline-loh/germline-loh.ensemble.vcf.gz'
    params:
        tumor=lambda wildcards: wildcards.tumor,
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        lib=config['resources']['targets_key'],
        outfile="data/work/{lib}/{tumor}/bcftools/germline-loh/germline-loh.ensemble.vcf"
    shell:
        """
        tabix -f -p vcf {input.vardict}
        tabix -f -p vcf {input.varscan2}
        python tumor-normal.ensemble-vcf.py -b {input.sites} -o {params.outfile} --vardict_vcf {input.vardict} --varscan2_vcf {input.varscan2} -T {params.tumor} -N {params.normal} -L {params.lib}
        """

rule ANNOVAR_germline_loh_ensemble_vcf:
    input:
        'data/work/{lib}/{tumor}/bcftools/germline-loh/germline-loh.ensemble.vcf.gz'
    output:
        'data/work/{lib}/{tumor}/annovar/germline-loh.ensemble.hg19_multianno.vcf.gz'
    params:
        humandb='$HOME/resources/annovar/humandb/',
        output_p='data/work/{lib}/{tumor}/annovar/germline-loh.ensemble'
    shell:
        """
        tabix -fp vcf {input}
        table_annovar.pl {input} {params.humandb} --buildver hg19 --vcfinput --outfile {params.output_p} --protocol refGene,cytoband,genomicSuperDups,dbscsnv11,avsnp150,dbnsfp35a,mcap,revel,popfreq_max_20150413,exac03,exac03nontcga,gnomad211_exome,intervar_20180118,icgc21,cosmic84_coding,cosmic84_noncoding,clinvar_20210501 --operation g,r,r,f,f,f,f,f,f,f,f,f,f,f,f,f,f -remove
        bgzip {params.output_p}.hg19_multianno.vcf
        tabix -fp vcf {params.output_p}.hg19_multianno.vcf.gz
        """
        #modify: params.build [hg19,hg38]

rule run_germline_loh_annovartools:
    input:
        'data/work/{lib}/{tumor}/annovar/germline-loh.ensemble.hg19_multianno.vcf.gz'
    output:
        'data/work/{lib}/{tumor}/annovar/germline-loh.ensemble.hg19_multianno.report.tsv'
    params:
        header='/home/bwubb/resources/annovar/ensemble.somatic-annotation-header.20200827.txt',
        mode='paired',
        samples=annovar_samples
    shell:
        'python annovartools.v03.py -I {input} -O {output} --header {params.header} -m {params.mode} {params.samples}'

#rule clean_up:
#    input:
#        pass

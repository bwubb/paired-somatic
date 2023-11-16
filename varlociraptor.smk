#https://github.com/bwubb

import os
import yaml
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

with open(config.get('analysis',{}).get('purity_table','purity.table'),'r') as u:
    PURITY=dict(line.split('\t') for line in u.read().splitlines())

TUMORS=PAIRS.keys()

with open(f"{config['project']['name']}.{config['resources']['targets_key']}.sample_name",'w') as name_file:
    name_file.write("TUMOR tumor\n")
    for i in TUMORS:
        name_file.write(f"{i} tumor\n")

##PYTHON
def map_vcf(wildcards):
    V={'lancet':f'data/work/{wildcards.lib}/{wildcards.tumor}/lancet/somatic.norm.clean.vcf.gz',
    'mutect2':f'data/work/{wildcards.lib}/{wildcards.tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz',
    'strelka2':f'data/work/{wildcards.lib}/{wildcards.tumor}/strelka2/somatic.norm.clean.vcf.gz',
    'vardict':f'data/work/{wildcards.lib}/{wildcards.tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz data/work/{wildcards.lib}/{wildcards.tumor}/vardict/loh.twice_filtered.norm.clean.vcf.gz data/work/{wildcards.lib}/{wildcards.tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz',
    'varscan2':f'data/work/{wildcards.lib}/{wildcards.tumor}/varscan2/somatic.fpfilter.norm.clean.vcf.gz data/work/{wildcards.lib}/{wildcards.tumor}/varscan2/loh.fpfilter.norm.clean.vcf.gz data/work/{wildcards.lib}/{wildcards.tumor}/varscan2/germline.fpfilter.norm.clean.vcf.gz'}
    return V[wildcards.caller]

def map_preprocess(wildcards):
    return {'bam':BAMS[wildcards.sample],'bcf':f'data/work/{wildcards.lib}/{wildcards.tumor}/varlociraptor/pass_candidates.bcf','aln':f'data/work/{wildcards.lib}/{wildcards.tumor}/varlociraptor/{wildcards.sample}.alignment-properties.json'}

def map_varlociraptor_scenario(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    filename=os.path.basename(config['analysis']['vlr'])
    #relapse1..n
    return {'tumor':f'data/work/{wildcards.lib}/{tumor}/varlociraptor/{tumor}.observations.bcf','normal':f'data/work/{wildcards.lib}/{tumor}/varlociraptor/{normal}.observations.bcf','scenario':f'data/work/{wildcards.lib}/{tumor}/varlociraptor/{filename}'}
#change {sample}'s to tumor/normal'

def genome_size(wildcards):
    G={'S04380110':'5.0e7','S07604715':'6.6e7','S31285117':'4.9e7','xgen-exome-research-panel-targets-grch37':'3.9e7','':'4.9e7','XGEN-EXOME-HYB-V2':'3.4e7'}
    return G[config['resources']['targets_key']]

##TARGET RULES
rule local_fdr:
    input:
        expand("data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.bcf",lib=config['resources']['targets_key'],tumor=TUMORS,scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0])
        #need better assignment of scenario wildcard
        #AmbiguousRuleException without constraints

rule collect_vep:
    input:
        expand("data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.report.csv",lib=config['resources']['targets_key'],tumor=TUMORS,scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0])

rule tmb:
    input:
        expand("data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.tmb.json",lib=config['resources']['targets_key'],tumor=TUMORS,scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0])

rule sites:
    input:
        expand("data/work/{lib}/{tumor}/varlociraptor/sites.txt",lib=config['resources']['targets_key'],tumor=TUMORS)

##SNAKEMAKE
wildcard_constraints:
    scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0]


rule candidate_tsv:
    input:
        map_vcf
    output:
        "data/work/{lib}/{tumor}/{caller}/candidates.tsv"
    shell:
        """
        bcftools concat -a {} | bcftools view -s tumor | bcftools query -f '%CHROM\t%POS\t.\t%REF\t%ALT\t.\t%FILTER\t.\n' | perl -lne '@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/\REJECT/; print join ("\\t", @row)' > {output}
        """
#bcftools concat {input} | bcftools

#add manta?
#
rule candidate_bcf:
    input:
        "data/work/{lib}/{tumor}/lancet/candidates.tsv",
        "data/work/{lib}/{tumor}/mutect2/candidates.tsv",
        "data/work/{lib}/{tumor}/strelka2/candidates.tsv",
        "data/work/{lib}/{tumor}/vardict/candidates.tsv",
        "data/work/{lib}/{tumor}/varscan2/candidates.tsv"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/pass_candidates.bcf"
    params:
        tsv=temp("data/work/{lib}/{tumor}/varlociraptor/candidates.tsv"),
        header="$HOME/resources/Vcf_files/headers/simplified-header-w_contig.grch38.vcf"#needs genome versioning
    shell:
        """
        cat {input} | grep PASS | sort -u > {params.tsv}
        cat {params.header} {params.tsv} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

rule varlociraptor_estimate_properties:
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

#bam
#bcf=pass_candidates.bcf
#aln=alignment-properties.json
rule varlociraptor_preprocess_sample:
    input:
        unpack(map_preprocess)
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{sample}.observations.bcf"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        varlociraptor preprocess variants {params.ref} --alignment-properties {input.aln} --bam {input.bam} --candidates {input.bcf} > {output}
        """

rule write_somatic_scenario:
    input:
        yaml=config['analysis']['vlr']
    output:
        yaml="data/work/{lib}/{tumor}/varlociraptor/{scenario}.yml"
    params:
        contamination=lambda wildcards: f"{1-float(PURITY.get(wildcards.tumor,'1.0')):.2f}"
    run:
        with open(input.yaml,'r') as i:
            scenario=yaml.load(i,Loader=yaml.SafeLoader)
            scenario['samples']['tumor']['contamination']['fraction']=float(params.contamination)
        with open(output.yaml,'w') as o:
            yaml.dump(scenario,o,default_flow_style=False)

rule varlociraptor_scenario:
    input:
        unpack(map_varlociraptor_scenario)
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.bcf"
    shell:
        """
        varlociraptor call variants generic --scenario {input.scenario} --obs tumor={input.tumor} normal={input.normal} > {output}
        """

rule varlociraptor_local_fdr:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.bcf"
    params:
        fdr=0.05
    shell:
        """
        varlociraptor filter-calls control-fdr --local {input} --events SOMATIC_TUMOR --fdr {params.fdr} > {output}
        bcftools index {output}
        """

rule varlociraptor_alt_depth_filter:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.filter.bcf"
    params:
        cutoff=8.0
    shell:
        """
        bcftools filter -i '(FORMAT/DP[1] * FORMAT/AF[1:0]) > {params.cutoff}' -O b -o {output} {input}
        bcftools index {output}
        """

rule bcftools_called_sites:
    input:
        "data/work/{lib}/{tumor}/{caller}/candidates.tsv"
    output:
        "data/work/{lib}/{tumor}/{caller}/called_sites.bcf"
    params:
        file=temp("data/work/{lib}/{tumor}/{caller}/tmp.tsv"),
        header="$HOME/resources/Vcf_files/headers/{caller}.contig_header.grch38.vcf"
        #how to do code base path?
    shell:
        """
        cat {input} | perl -lne '@row=split /\\t/; $row[7] = uc("{wildcards.caller}=").$row[6]; print join ("\\t", @row)' > {params.file}
        cat {params.header} {params.file} | bcftools view -O b -o {output}
        bcftools index {output}
        """

rule bcftools_merge_sites:
    input:
        lancet="data/work/{lib}/{tumor}/lancet/called_sites.bcf",
        mutect2="data/work/{lib}/{tumor}/mutect2/called_sites.bcf",
        strelka2="data/work/{lib}/{tumor}/strelka2/called_sites.bcf",
        vardict="data/work/{lib}/{tumor}/vardict/called_sites.bcf",
        varscan2="data/work/{lib}/{tumor}/varscan2/called_sites.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/paired_somatic.bcf"
    shell:
        """
        bcftools merge -m none {input} | bcftools sort -O b -o {output}
        bcftools index {output}
        """

rule bcftools_annotate:
    input:
        bcf="data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.bcf",
        sites="data/work/{lib}/{tumor}/varlociraptor/paired_somatic.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.bcf"
    shell:
        """
        bcftools annotate -a {input.sites} -c LANCET,MUTECT2,STRELKA2,VARDICT,VARSCAN2 -O b -o {output} {input.bcf}
        """

rule varlociraptor_vep:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.vcf.gz"
    params:
        in_vcf=temp('data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vcf'),
        out_vcf='data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.vcf',
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='$HOME/.vep/Plugins/loftee',#check
        utr=config['resources']['utrannotator'],
        ref_fa="data/work/{lib}/{tumor}/varlociraptor/{scenario}.reference.fa",
        mut_fa="data/work/{lib}/{tumor}/varlociraptor/{scenario}.mutated.fa"
    shell:
        """
        bcftools view -O v -o {params.in_vcf} {input}

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
        --vcf_info_field ANN \
        --plugin NMD \
        --plugin ProteinSeqs,{params.ref_fa},{params.mut_fa} \
        --plugin Downstream \
        --plugin REVEL,{params.revel} \
        --plugin SpliceAI,snv={params.splice_snv},indel={params.splice_indel} \
        --plugin gnomADc,{params.gnomAD} \
        --plugin UTRannotator,{params.utr} \
        --plugin LoF,loftee_path:$HOME/software/loftee,human_ancestor_fa:$HOME/.vep/Plugins/loftee/human_ancestor.fa.gz \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

        bgzip {params.out_vcf}
        tabix -fp vcf {output}
        """
        #loftee cut again until after lab meeting
        #--plugin LoF,loftee_path:$HOME/software/loftee,human_ancestor_fa:$HOME/.vep/Plugins/loftee/human_ancestor.fa.gz
        #move ref mut fa to output after testing is complete.
        #possibly bgzip/index output fasta file

#I need to go back to old snakemake parsing for proper github pathing
rule varlociraptor_vep_report:
    input:
        vcf="data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.vcf.gz"
    output:
        csv="data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.report.csv"
    params:
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        python vep_vcf_parser.py -i {input.vcf} -o {output.csv} --mode tumor_normal,{wildcards.tumor},{params.normal} everything
        """

rule varlociraptor_estimate_mutational_burden:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.vcf.gz"
    output:
        bcf="data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.bcf",
        json="data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.tmb.json"
        #"data/work/{lib}/{tumor}/varlociraptor/{scenario}.local-fdr.paired_somatic.vep.tmb.svg"
    params:
        size=genome_size
    shell:
        """
        bcftools view -O b -o {output.bcf} {input}
        varlociraptor estimate mutational-burden --plot-mode curve --coding-genome-size {params.size} --sample tumor --events SOMATIC_TUMOR < {output.bcf} > {output.json}
        """

rule bcftools_isec:
    input:
        lancet="data/work/{lib}/{tumor}/lancet/somatic.norm.clean.vcf.gz",
        mutect2="data/work/{lib}/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz",
        strelka2="data/work/{lib}/{tumor}/strelka2/somatic.norm.clean.vcf.gz",
        vardict="data/work/{lib}/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz",
        varscan2="data/work/{lib}/{tumor}/varscan2/somatic.fpfilter.norm.clean.vcf.gz"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/somatic_sites.txt"
    shell:
        """
        bcftools isec -n+1 {input.lancet} {input.mutect2} {input.strelka2} {input.vardict} {input.varscan2} > {output}
        """

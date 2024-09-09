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

##PYTHON
def map_vcf(wildcards):
    V={'lancet':f'data/work/{wildcards.lib}/{wildcards.tumor}/lancet/somatic.norm.clean.vcf.gz',
    'mutect2':f'data/work/{wildcards.lib}/{wildcards.tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz',
    'strelka2':f'data/work/{wildcards.lib}/{wildcards.tumor}/strelka2/somatic.norm.clean.vcf.gz',
    'vardict':[f'data/work/{wildcards.lib}/{wildcards.tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz',f'data/work/{wildcards.lib}/{wildcards.tumor}/vardict/loh.twice_filtered.norm.clean.vcf.gz',f'data/work/{wildcards.lib}/{wildcards.tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz'],
    'varscan2':[f'data/work/{wildcards.lib}/{wildcards.tumor}/varscan2/somatic.fpfilter.norm.clean.vcf.gz',f'data/work/{wildcards.lib}/{wildcards.tumor}/varscan2/loh.fpfilter.norm.clean.vcf.gz',f'data/work/{wildcards.lib}/{wildcards.tumor}/varscan2/germline.fpfilter.norm.clean.vcf.gz']}
    #'vardict':f'data/work/{wildcards.lib}/{wildcards.tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz',
    #'varscan2':f'data/work/{wildcards.lib}/{wildcards.tumor}/varscan2/somatic.fpfilter.norm.clean.vcf.gz'}
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
rule vlr_paired_somatic:
    input:
        expand("data/work/{lib}/{tumor}/varlociraptor/{scenario}.paired_somatic.vep.vcf",lib=config['resources']['targets_key'],tumor=TUMORS,scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0])

##SNAKEMAKE
wildcard_constraints:
    scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0]

#header file paths need to be altered when repository is working
rule lancet_vlr_input:
    input:
        vcf="data/work/{targets}/{tumor}/lancet/somatic.norm.clean.vcf.gz"
    output:
        tsv="data/work/{targets}/{tumor}/varlociraptor/lancet.candidates.tsv",
        bcf="data/work/{targets}/{tumor}/varlociraptor/lancet.candidates.bcf"
    params:
        tsv=temp("data/work/{targets}/{tumor}/varlociraptor/lancet1.tsv"),
        header="/home/bwubb/resources/Vcf_files/headers/simplified-header-w_contig.grch38.lancet.vcf"
    shell:
        """
        bcftools view -G -Ov {input.vcf} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=LANCET|SOMATIC;LANCET=%FILTER\\n' > {params.tsv}

        cat {params.header} {params.tsv} | bcftools sort -Ob -W=csi -o {output.bcf}
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\t.\\n' {output.bcf} > {output.tsv}
        """

rule mutect2_vlr_input:
    input:
        vcf="data/work/{targets}/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    output:
        tsv="data/work/{targets}/{tumor}/varlociraptor/mutect2.candidates.tsv",
        bcf="data/work/{targets}/{tumor}/varlociraptor/mutect2.candidates.bcf"
    params:
        tsv=temp("data/work/{targets}/{tumor}/varlociraptor/mutect2_1.tsv"),
        header="/home/bwubb/resources/Vcf_files/headers/simplified-header-w_contig.grch38.mutect2.vcf"
    shell:
        """
        bcftools view -G -Ov {input.vcf} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=MUTECT2|SOMATIC;MUTECT2=%FILTER\\n' > {params.tsv}

        cat {params.header} {params.tsv} | bcftools sort -Ob -W=csi -o {output.bcf}
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\t.\\n' {output.bcf} > {output.tsv}
        """

rule strelka2_vlr_input:
    input:
        vcf="data/work/{targets}/{tumor}/strelka2/somatic.norm.clean.vcf.gz"
    output:
        tsv="data/work/{targets}/{tumor}/varlociraptor/strelka2.candidates.tsv",
        bcf="data/work/{targets}/{tumor}/varlociraptor/strelka2.candidates.bcf"
    params:
        tsv=temp("data/work/{targets}/{tumor}/varlociraptor/strelka2_1.tsv"),
        header="/home/bwubb/resources/Vcf_files/headers/simplified-header-w_contig.grch38.strelka2.vcf"
    shell:
        """
        bcftools view -G -Ov {input.vcf} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=STRELKA2|SOMATIC;STRELKA2=%FILTER\\n' > {params.tsv}

        cat {params.header} {params.tsv} | bcftools sort -Ob -W=csi -o {output.bcf}
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\t.\\n' {output.bcf} > {output.tsv}
        """

rule vardict_vlr_input:
    input:
        vcf1="data/work/{targets}/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz",
        vcf2="data/work/{targets}/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz"
    output:
        tsv="data/work/{targets}/{tumor}/varlociraptor/vardict.candidates.tsv",
        bcf="data/work/{targets}/{tumor}/varlociraptor/vardict.candidates.bcf"
    params:
        tsv1=temp("data/work/{targets}/{tumor}/varlociraptor/vardict1.tsv"),
        tsv2=temp("data/work/{targets}/{tumor}/varlociraptor/vardict2.tsv"),
        header="/home/bwubb/resources/Vcf_files/headers/simplified-header-w_contig.grch38.vardict.vcf"
    shell:
        """
        bcftools view -G -Ov {input.vcf1} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=VARDICT|SOMATIC;VARDICT=%FILTER\\n' > {params.tsv1}

        bcftools view -G -Ov {input.vcf2} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=VARDICT|GERMLINE;VARDICT=%FILTER\\n' > {params.tsv2}

        cat {params.header} {params.tsv1} {params.tsv2} | bcftools sort -Ob -W=csi -o {output.bcf}
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\t.\\n' {output.bcf} > {output.tsv}
        """

rule varscan2_vlr_input:
    input:
        vcf1="data/work/{targets}/{tumor}/varscan2/somatic.fpfilter.fix.norm.clean.vcf.gz",
        vcf2="data/work/{targets}/{tumor}/varscan2/germline.fpfilter.fix.norm.clean.vcf.gz"
    output:
        tsv="data/work/{targets}/{tumor}/varlociraptor/varscan2.candidates.tsv",
        bcf="data/work/{targets}/{tumor}/varlociraptor/varscan2.candidates.bcf"
    params:
        tsv1=temp("data/work/{targets}/{tumor}/varlociraptor/varscan2_1.tsv"),
        tsv2=temp("data/work/{targets}/{tumor}/varlociraptor/varscan2_2.tsv"),
        header="/home/bwubb/resources/Vcf_files/headers/simplified-header-w_contig.grch38.varscan2.vcf"
    shell:
        """
        bcftools view -G -Ov {input.vcf1} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=VARSCAN2|SOMATIC;VARSCAN2=%FILTER\\n' > {params.tsv1}

        bcftools view -G -Ov {input.vcf2} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=VARSCAN2|GERMLINE;VARSCAN2=%FILTER\\n' > {params.tsv2}

        cat {params.header} {params.tsv1} {params.tsv2} | bcftools sort -Ob -W=csi -o {output.bcf}
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\t.\\n' {output.bcf} > {output.tsv}
        """

rule candidate_bcf:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/lancet.candidates.tsv",
        "data/work/{lib}/{tumor}/varlociraptor/mutect2.candidates.tsv",
        "data/work/{lib}/{tumor}/varlociraptor/strelka2.candidates.tsv",
        "data/work/{lib}/{tumor}/varlociraptor/vardict.candidates.tsv",
        "data/work/{lib}/{tumor}/varlociraptor/varscan2.candidates.tsv"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/pass_candidates.bcf"
    params:
        tsv=temp("data/work/{lib}/{tumor}/varlociraptor/pass_candidates.tsv"),
        header="$HOME/resources/Vcf_files/headers/simplified-header-w_contig.grch38.vcf"#needs genome versioning
    shell:
        """
        cat {input} | grep PASS | sort -u > {params.tsv}
        cat {params.header} {params.tsv} | bcftools sort -Ob -W=csi -o {output}
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

rule write_scenario:
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

rule varlociraptor_call_scenario:
    input:
        unpack(map_varlociraptor_scenario)
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.bcf"
    shell:
        """
        varlociraptor call variants generic --scenario {input.scenario} --obs tumor={input.tumor} normal={input.normal} > {output}
        """

rule bcftools_merge_sites:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/lancet.candidates.bcf",
        "data/work/{lib}/{tumor}/varlociraptor/mutect2.candidates.bcf",
        "data/work/{lib}/{tumor}/varlociraptor/strelka2.candidates.bcf",
        "data/work/{lib}/{tumor}/varlociraptor/vardict.candidates.bcf",
        "data/work/{lib}/{tumor}/varlociraptor/varscan2.candidates.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/paired_somatic.bcf"
    shell:
        """
        bcftools merge -m none --info-rules CATEGORY:join {input} | bcftools sort -Ob -W=csi -o {output}
        """

rule bcftools_annotate_sites:
    input:
        sites="data/work/{lib}/{tumor}/varlociraptor/paired_somatic.bcf",
        bcf="data/work/{lib}/{tumor}/varlociraptor/{scenario}.bcf"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.paired_somatic.vcf.gz"
    shell:
        """
        bcftools index -f {input.bcf}
        bcftools index -f {input.sites}
        bcftools annotate -a {input.sites} -c LANCET,MUTECT2,STRELKA2,VARDICT,VARSCAN2,CATEGORY -Oz -W=tbi -o {output} {input.bcf}
        """

#VEP appears to run a lot slower if the input is compressed.
rule vep_annotation:
    input:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.paired_somatic.vcf.gz"
    output:
        "data/work/{lib}/{tumor}/varlociraptor/{scenario}.paired_somatic.vep.vcf"
    params:
        vcf=temp("data/work/{lib}/{tumor}/varlociraptor/{scenario}.paired_somatic.vcf")
    shell:
        """
        bcftools view -Ov -o {params.vcf} {input}

        singularity run -H $PWD:/home \
        --bind /home/bwubb/resources:/opt/vep/resources \
        --bind /home/bwubb/.vep:/opt/vep/.vep \
        /appl/containers/vep112.sif vep \
        --dir /opt/vep/.vep \
        -i {params.vcf} \
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
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
        """

##Author: Brad Wubbenhorst
# https://github.com/bwubb
#
# Varlociraptor is a Rust binary. Host cargo builds need libclang (not on RHEL9
# modules here), so we run the BioContainers SIF instead:
#   singularity pull $HOME/containers/varlociraptor.sif \
#     docker://quay.io/biocontainers/varlociraptor:8.9.5--h24073b4_0
# Override with varlociraptor.container in config if needed.

import os
import yaml

#include ./lancet2.smk
#include ./mutect2.smk
#include ./strelka2.smk
#include ./vardictjava.smk
#include ./varscan2.smk

### INIT ###

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

with open(config.get('analysis',{}).get('purity_table','purity.table'),'r') as u:
    PURITY=dict(line.split('\t') for line in u.read().splitlines())

TUMORS=PAIRS.keys()

VLR_SIF=config.get('varlociraptor',{}).get('container','/home/bwubb/containers/varlociraptor.sif')

# Same bind pattern as lancet2 / strelka2 / cnvkit_singularity.
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

def map_preprocess(wildcards):
    return {'bam':BAMS[wildcards.sample],
            'bcf':f'data/work/{wildcards.tumor}/varlociraptor/pass_candidates.bcf',
            'aln':f'data/work/{wildcards.tumor}/varlociraptor/{wildcards.sample}.alignment-properties.json'}

def map_varlociraptor_scenario(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    filename=os.path.basename(config['analysis']['vlr'])
    return {'tumor':f'data/work/{tumor}/varlociraptor/{tumor}.observations.bcf',
            'normal':f'data/work/{tumor}/varlociraptor/{normal}.observations.bcf',
            'scenario':f'data/work/{tumor}/varlociraptor/{filename}'}

def genome_size(wildcards):
    G={'S04380110':'5.0e7','S07604715':'6.6e7','S31285117':'4.9e7','xgen-exome-research-panel-targets-grch37':'3.9e7','':'4.9e7','XGEN-EXOME-HYB-V2':'3.4e7'}
    return G[config['resources']['targets_key']]

### TARGET RULES ###

rule vlr_paired_somatic:
    input:
        expand("data/work/{tumor}/varlociraptor/{scenario}.paired_somatic.vep.vcf",
               tumor=TUMORS,
               scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0])

### SNAKEMAKE ###

wildcard_constraints:
    scenario=os.path.splitext(os.path.basename(config['analysis']['vlr']))[0]

#header file paths need to be altered when repository is working
rule lancet2_vlr_input:
    input:
        vcf="data/work/{tumor}/lancet2/somatic.norm.clean.vcf.gz"
    output:
        tsv="data/work/{tumor}/varlociraptor/lancet2.candidates.tsv",
        bcf="data/work/{tumor}/varlociraptor/lancet2.candidates.bcf"
    resources:
        mem_mb=6144
    params:
        tsv=temp("data/work/{tumor}/varlociraptor/lancet2_1.tsv"),
        header="$HOME/resources/Vcf_files/headers/lancet2.contig_header.grch38.vcf"
    shell:
        """
        bcftools view -G -Ov {input.vcf} |
        perl -lne 'if (/^#/) {{print}} else {{@row=split /\\t/; $row[6] =~ s/^(?!PASS).*/FAIL/; print join ("\\t", @row)}}' |
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\tCATEGORY=LANCET2|SOMATIC;LANCET2=%FILTER\\n' > {params.tsv}

        cat {params.header} {params.tsv} | bcftools sort -Ob -W=csi -o {output.bcf}
        bcftools query -f'%CHROM\\t%POS\\t.\\t%REF\\t%ALT\\t.\\t%FILTER\\t.\\n' {output.bcf} > {output.tsv}
        """

rule mutect2_vlr_input:
    input:
        vcf="data/work/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    output:
        tsv="data/work/{tumor}/varlociraptor/mutect2.candidates.tsv",
        bcf="data/work/{tumor}/varlociraptor/mutect2.candidates.bcf"
    resources:
        mem_mb=6144
    params:
        tsv=temp("data/work/{tumor}/varlociraptor/mutect2_1.tsv"),
        header="$HOME/resources/Vcf_files/headers/mutect2.contig_header.grch38.vcf"
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
        vcf="data/work/{tumor}/strelka2/somatic.norm.clean.vcf.gz"
    output:
        tsv="data/work/{tumor}/varlociraptor/strelka2.candidates.tsv",
        bcf="data/work/{tumor}/varlociraptor/strelka2.candidates.bcf"
    resources:
        mem_mb=6144
    params:
        tsv=temp("data/work/{tumor}/varlociraptor/strelka2_1.tsv"),
        header="$HOME/resources/Vcf_files/headers/strelka2.contig_header.grch38.vcf"
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
        vcf1="data/work/{tumor}/vardict/somatic.filter2.norm.clean.vcf.gz",
        vcf2="data/work/{tumor}/vardict/germline.filter2.norm.clean.vcf.gz"
    output:
        tsv="data/work/{tumor}/varlociraptor/vardict.candidates.tsv",
        bcf="data/work/{tumor}/varlociraptor/vardict.candidates.bcf"
    resources:
        mem_mb=6144
    params:
        tsv1=temp("data/work/{tumor}/varlociraptor/vardict1.tsv"),
        tsv2=temp("data/work/{tumor}/varlociraptor/vardict2.tsv"),
        header="$HOME/resources/Vcf_files/headers/vardict.contig_header.grch38.vcf"
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
        vcf1="data/work/{tumor}/varscan2/somatic.fpfilter.processed.norm.clean.vcf.gz",
        vcf2="data/work/{tumor}/varscan2/germline.fpfilter.processed.norm.clean.vcf.gz"
    output:
        tsv="data/work/{tumor}/varlociraptor/varscan2.candidates.tsv",
        bcf="data/work/{tumor}/varlociraptor/varscan2.candidates.bcf"
    resources:
        mem_mb=6144
    params:
        tsv1=temp("data/work/{tumor}/varlociraptor/varscan2_1.tsv"),
        tsv2=temp("data/work/{tumor}/varlociraptor/varscan2_2.tsv"),
        header="$HOME/resources/Vcf_files/headers/varscan2.contig_header.grch38.vcf"
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
        "data/work/{tumor}/varlociraptor/lancet2.candidates.tsv",
        "data/work/{tumor}/varlociraptor/mutect2.candidates.tsv",
        "data/work/{tumor}/varlociraptor/strelka2.candidates.tsv",
        "data/work/{tumor}/varlociraptor/vardict.candidates.tsv",
        "data/work/{tumor}/varlociraptor/varscan2.candidates.tsv"
    output:
        "data/work/{tumor}/varlociraptor/pass_candidates.bcf"
    resources:
        mem_mb=6144
    params:
        tsv=temp("data/work/{tumor}/varlociraptor/pass_candidates.tsv"),
        header="$HOME/resources/Vcf_files/headers/paired_somatic.contig_header.grch38.vcf"
    shell:
        """
        cat {input} | grep PASS | sort -u > {params.tsv}
        cat {params.header} {params.tsv} | bcftools sort -Ob -W=csi -o {output}
        """

rule varlociraptor_estimate_properties:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "data/work/{tumor}/varlociraptor/{sample}.alignment-properties.json"
    resources:
        mem_mb=18432
    params:
        ref=config['reference']['fasta'],
        exec=singularity_exec(VLR_SIF)
    shell:
        """
        {params.exec} varlociraptor estimate alignment-properties {params.ref} --bams {input.bam} > {output}
        """

rule varlociraptor_preprocess_sample:
    input:
        unpack(map_preprocess)
    output:
        "data/work/{tumor}/varlociraptor/{sample}.observations.bcf"
    resources:
        mem_mb=18432
    params:
        ref=config['reference']['fasta'],
        exec=singularity_exec(VLR_SIF)
    shell:
        """
        {params.exec} varlociraptor preprocess variants {params.ref} --alignment-properties {input.aln} --bam {input.bam} --candidates {input.bcf} > {output}
        """

rule write_scenario:
    input:
        yaml=config['analysis']['vlr']
    output:
        yaml="data/work/{tumor}/varlociraptor/{scenario}.yml"
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
        "data/work/{tumor}/varlociraptor/{scenario}.bcf"
    resources:
        mem_mb=18432
    params:
        exec=singularity_exec(VLR_SIF)
    shell:
        """
        {params.exec} varlociraptor call variants generic --scenario {input.scenario} --obs tumor={input.tumor} normal={input.normal} > {output}
        """

rule bcftools_merge_sites:
    input:
        "data/work/{tumor}/varlociraptor/lancet2.candidates.bcf",
        "data/work/{tumor}/varlociraptor/mutect2.candidates.bcf",
        "data/work/{tumor}/varlociraptor/strelka2.candidates.bcf",
        "data/work/{tumor}/varlociraptor/vardict.candidates.bcf",
        "data/work/{tumor}/varlociraptor/varscan2.candidates.bcf"
    output:
        "data/work/{tumor}/varlociraptor/paired_somatic.bcf"
    resources:
        mem_mb=6144
    shell:
        """
        bcftools merge -m none --info-rules CATEGORY:join {input} | bcftools sort -Ob -W=csi -o {output}
        """

rule bcftools_annotate_sites:
    input:
        sites="data/work/{tumor}/varlociraptor/paired_somatic.bcf",
        bcf="data/work/{tumor}/varlociraptor/{scenario}.bcf"
    output:
        "data/work/{tumor}/varlociraptor/{scenario}.paired_somatic.vcf.gz"
    resources:
        mem_mb=6144
    shell:
        """
        bcftools index -f {input.bcf}
        bcftools index -f {input.sites}
        bcftools annotate -a {input.sites} -c LANCET2,MUTECT2,STRELKA2,VARDICT,VARSCAN2,CATEGORY -Oz -W=tbi -o {output} {input.bcf}
        """

#VEP appears to run a lot slower if the input is compressed.
rule vep_annotation:
    input:
        "data/work/{tumor}/varlociraptor/{scenario}.paired_somatic.vcf.gz"
    output:
        "data/work/{tumor}/varlociraptor/{scenario}.paired_somatic.vep.vcf"
    resources:
        mem_mb=32768
    params:
        vcf=temp("data/work/{tumor}/varlociraptor/{scenario}.paired_somatic.vcf")
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
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
        """

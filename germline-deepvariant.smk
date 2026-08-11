##Author: Brad Wubbenhorst

### INIT ###

with open(config.get('project',{}).get('germline_list','germline.list'),'r') as s:
    SAMPLES=s.read().splitlines()

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

### FUNCTIONS ###

def bam_input(wildcards):
    return BAMS[wildcards.sample]

### SNAKEMAKE ###

rule germline_deepvariant:
    input:
        expand("data/work/{sample}/deepvariant/germline.norm.vep.vcf",sample=SAMPLES)

rule run_deepvariant:
    input:
        bam=bam_input
    output:
        vcf="data/work/{sample}/deepvariant/germline.vcf.gz",
        g_vcf="data/work/{sample}/deepvariant/{sample}.g.vcf.gz"
    threads: 4
    resources:
        mem_mb=18432
    params:
        bed=config['resources']['targets_bed'],
        ref=config['reference']['fasta']
    shell:
        """
        WORKDIR="$(pwd)"
        OUT_VCF="$WORKDIR/{output.vcf}"
        OUT_GVCF="$WORKDIR/{output.g_vcf}"
        SCRATCH="/scratch/bwubb/deepvariant/{wildcards.sample}"

        mkdir -p "$(dirname "$OUT_VCF")" "$SCRATCH/tmp" "$SCRATCH/intermediate"
        export TMPDIR="$SCRATCH/tmp"

        # /home/bwubb/projects/... often symlinks into /project/...; bind both so
        # absolute paths and symlink targets resolve inside the container.
        # Do not use --pwd here: Singularity realpath() of a /home -> /project
        # symlink fails hard if it chdirs before/without the /project bind.
        singularity run \
          --bind /usr/lib/locale:/usr/lib/locale \
          --bind /project:/project \
          --bind /home/bwubb:/home/bwubb \
          --bind /scratch:/scratch \
          --bind "$WORKDIR:$WORKDIR" \
          /appl/containers/deepvariant_1.4.0.sif \
          run_deepvariant \
          --model_type=WES \
          --ref={params.ref} \
          --reads={input.bam} \
          --regions={params.bed} \
          --output_vcf="$OUT_VCF" \
          --output_gvcf="$OUT_GVCF" \
          --intermediate_results_dir="$SCRATCH/intermediate" \
          --num_shards={threads} \
          --sample_name={wildcards.sample}
        """

rule deepvariant_norm:
    input:
        "data/work/{sample}/deepvariant/germline.vcf.gz"
    output:
        "data/work/{sample}/deepvariant/germline.norm.vcf.gz"
    resources:
        mem_mb=6144
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -W=tbi -Oz -o {output}
        """

rule deepvariant_vep:
    input:
        "data/work/{sample}/deepvariant/germline.norm.vcf.gz"
    output:
        "data/work/{sample}/deepvariant/germline.norm.vep.vcf"
    resources:
        mem_mb=32768
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
        --custom /opt/vep/.vep/clinvar/vcf_GRCh38/clinvar.autogvp.vcf.gz,ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN,AutoGVP \
        --plugin AlphaMissense,file=/opt/vep/.vep/alphamissense/AlphaMissense_GRCh38.tsv.gz \
        --plugin MaveDB,file=/opt/vep/.vep/mavedb/MaveDB_variants.tsv.gz
        """

rule parse_deepvariant_vep:
    input:
        "data/work/{sample}/deepvariant/germline.norm.vep.vcf"
    output:
        "data/work/{sample}/deepvariant/germline.norm.vep.report.csv"
    shell:
        "python vep_vcf_parser2.py -i {input} -o {output} -m single,{wildcards.sample}"

rule report_germline_results:
    input:
        expand("data/work/{sample}/deepvariant/germline.norm.vep.report.csv",sample=SAMPLES)
    output:
        "germline.report.html"
    shell:
        ""
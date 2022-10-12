


with open(config.get('project',{}).get('germline_list','germline.list'),'r') as s:
    SAMPLES=s.read().splitlines()

with open(config.get('project',{}).get('bam_table','bam.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

rule all:
    input:
        expand("data/work/{library}/{sample}/vardict/vardict.vep.vcf.gz",sample=SAMPLES,library=config['resources']['targets_key'])

rule run_germline_vardict:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "data/work/{library}/{sample}/vardict/vardict.vcf.gz"
    params:
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed'],
        path="$HOME/software/VarDictJava/VarDict",
        AF_THR=0.01
    threads:
        4
    shell:
        "$HOME/software/VarDictJava/build/install/VarDict/bin/VarDict -th {threads} -G {params.ref} -f {params.AF_THR} -N {wildcards.sample} -b {input.bam} -z -c 1 -S 2 -E 3 -g 4 {params.bed} | {params.path}/teststrandbias.R | {params.path}/var2vcf_valid.pl -N {wildcards.sample} -f {params.AF_THR} | bgzip -c > {output}"

rule run_germline_vep:
    input:
        "data/work/{library}/{sample}/vardict/vardict.vcf.gz"
    output:
        "data/work/{library}/{sample}/vardict/vardict.vep.vcf.gz"
    params:
        out_vcf="data/work/{library}/{sample}/vardict/vardict.vep.vcf",
        assembly=config['reference']['key'],
        fa=config['reference']['fasta'],
        splice_snv=config['resources']['splice_snv'],
        splice_indel=config['resources']['splice_indel'],
        gnomAD=config['resources']['gnomAD'],
        clinvar=config['resources']['clinvar'],
        revel=config['resources']['revel'],
        loftee='$HOME/.vep/Plugins/loftee',#check
        utr=config['resources']['utr'],
        ref_fa="data/work/{library}/{sample}/vardict/reference.vardict.fa",
        mut_fa="data/work/{library}/{sample}/vardict/mutated.vardict.fa"
    shell:
        """
        vep -i {input} -o {params.out_vcf} \
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
        --plugin UTRannotator,{params.utr} \
        --custom {params.clinvar},ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

        bgzip {params.out_vcf} && tabix -fp vcf {output}
        """
        #--plugin LoF,loftee_path:{params.loftee} \

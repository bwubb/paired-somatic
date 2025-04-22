import os

### INIT ###

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as s:
    SAMPLES=s.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

### FUNCTIONS ###

def paired_bams(wildcards):
    ref=config['reference']['key']
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

### SNAKEMAKE ###

localrules: vardict_sample_name

#should be bcf format...
rule run_vardict:
    input:
        expand("data/final/{tumor}/{tumor}.vardict.somatic.vcf.gz", tumor=PAIRS.keys()),
        expand("data/final/{tumor}/{tumor}.vardict.germline.vcf.gz", tumor=PAIRS.keys())

#switch to 1filter/2filter
rule filter_applied_vardict_somatic:
    input:
        expand("data/work/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz",tumor=PAIRS.keys())

rule filter_applied_vardict:
    input:
        expand("data/work/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz",tumor=PAIRS.keys()),
        expand("data/work/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz",tumor=PAIRS.keys())

rule vardict_sample_name:
    output:
        "data/work/{tumor}/vardict/sample.name"
    params:
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        echo -e "TUMOR\\t{wildcards.tumor}\\nNORMAL\\t{params.normal}" > {output}
        echo -e "tumor\\t{wildcards.tumor}\\nnormal\\t{params.normal}" >> {output}
        """

#First a more lenient -P val, not sure what
rule vardict_paired_main:
    input:
        unpack(paired_bams)
    output:
        "data/work/{tumor}/vardict/variants.vcf.gz"
    params:
        init=temp("data/work/{tumor}/vardict/output.vcf.gz"),
        bam=lambda wildcards: BAMS[wildcards.tumor],
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed'],
        path="$HOME/software/VarDictJava/VarDict",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        AF_THR=0.01
    threads:
        4
    shell:
        """
        VarDict -th {threads} -G {params.ref} -f {params.AF_THR} -N {wildcards.tumor} -b '{input.tumor}|{input.normal}' -c 1 -S 2 -E 3 -g 4 {params.bed} | {params.path}/testsomatic.R | {params.path}/var2vcf_paired.pl -N '{wildcards.tumor}|{params.normal}' -f {params.AF_THR} | bgzip -c > {params.init}
        tabix -p vcf {params.init}
        gatk UpdateVCFSequenceDictionary -V {params.init} --source-dictionary {params.bam} --output {output}
        """

rule vardict_filter:
    input:
        "data/work/{tumor}/vardict/variants.vcf.gz"
    output:
        once="data/work/{tumor}/vardict/variants.once_filtered.vcf.gz",
        twice="data/work/{tumor}/vardict/variants.twice_filtered.vcf.gz"
    shell:
        """
        bcftools filter --threads {threads} -e '((FORMAT/AF[0] * FORMAT/DP[0] < 6) && ((FORMAT/MQ[0] < 55.0 && FORMAT/NM[0] > 1.0) || (FORMAT/MQ[0] < 60.0 && FORMAT/NM[0] > 2.0) || (FORMAT/DP[0] < 10) || (QUAL < 45)))' -s filter_1 -m + -W=tbi -O z {input} > {output.once}
        bcftools filter --threads {threads} -e 'FORMAT/AF[0] < 0.2 && FORMAT/QUAL[0] < 55 && INFO/SSF[0] > 0.06' -s filter_2 -m + -W=tbi -O z {output.once} > {output.twice}
        """

#Maybe filter is only good for somatic?
#What germline/loh are filtered out?
rule vardict_split:
    input:
        raw="data/work/{tumor}/vardict/variants.vcf.gz",
        filtered="data/work/{tumor}/vardict/variants.twice_filtered.vcf.gz"
    output:
        somatic1="data/work/{tumor}/vardict/somatic.vcf.gz",
        germline1="data/work/{tumor}/vardict/germline.vcf.gz",
        somatic2="data/work/{tumor}/vardict/somatic.twice_filtered.vcf.gz",
        germline2="data/work/{tumor}/vardict/germline.twice_filtered.vcf.gz"
    shell:
        """
        bcftools view -i 'INFO/STATUS==\"StrongSomatic\" || INFO/STATUS==\"LikelySomatic\"' -W=tbi -O z -o {output.somatic1} {input.raw}
        bcftools view -i 'INFO/STATUS==\"Germline\" || INFO/STATUS==\"StrongLOH\" || INFO/STATUS==\"LikelyLOH\"' -W=tbi -O z -o {output.germline1} {input.raw}

        bcftools view -i 'INFO/STATUS==\"StrongSomatic\" || INFO/STATUS==\"LikelySomatic\"' -W=tbi -O z -o {output.somatic2} {input.filtered}
        bcftools view -i 'INFO/STATUS==\"Germline\" || INFO/STATUS==\"StrongLOH\" || INFO/STATUS==\"LikelyLOH\"' -W=tbi -O z -o {output.germline2} {input.filtered}
        """

rule vardict_somatic_normalized:
    input:
        "data/work/{tumor}/vardict/somatic.twice_filtered.vcf.gz"
    output:
        norm="data/work/{tumor}/vardict/somatic.twice_filtered.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -W=tbi -O z -o {output.norm}
        """

rule vardict_somatic_clean:
    input:
        name="data/work/{tumor}/vardict/sample.name",
        vcf="data/work/{tumor}/vardict/somatic.twice_filtered.norm.vcf.gz"
    output:
        clean="data/work/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("data/work/{tumor}/vardict/temp.h.vcf.gz")
    shell:
        """
        bcftools reheader -f {params.fai} -s {input.name} -o {params.vcf} {input.vcf}
        bcftools index {params.vcf}

        bcftools view -s {wildcards.tumor},{params.normal} -e 'ALT="*"' -R {params.regions} {params.vcf} | \
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output.clean}
        """

rule vardict_somatic_final:
    input:
        "data/work/{tumor}/vardict/somatic.vcf.gz",
        "data/work/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz"
    output:
        "data/final/{tumor}/{tumor}.vardict.somatic.vcf.gz",
        "data/final/{tumor}/{tumor}.vardict.somatic.final.vcf.gz"
    shell:
        """
        bcftools view -W=tbi -Oz -o {output[0]} {input[0]}
        bcftools view -W=tbi -Oz -o {output[1]} {input[1]}
        """

rule vardict_germline_normalized:
    input:
        "data/work/{tumor}/vardict/germline.twice_filtered.vcf.gz"
    output:
        norm="data/work/{tumor}/vardict/germline.twice_filtered.norm.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -W=tbi -Oz -o {output.norm}
        """

rule vardict_germline_clean:
    input:
        name="data/work/{tumor}/vardict/sample.name",
        vcf="data/work/{tumor}/vardict/germline.twice_filtered.norm.vcf.gz"
    output:
        clean="data/work/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("data/work/{tumor}/vardict/temp.g.vcf.gz")
    shell:
        """
        bcftools reheader -f {params.fai} -s {input.name} -o {params.vcf} {input.vcf}
        bcftools index {params.vcf}

        bcftools view -s {wildcards.tumor},{params.normal} -e 'ALT="*"' -R {params.regions} {params.vcf} | \
        bcftools annotate --set-id '%CHROM\_%POS\_%REF\_%ALT' | \
        bcftools sort -W=tbi -Oz -o {output.clean}
        """

rule vardict_germline_final:
    input:
        "data/work/{tumor}/vardict/germline.vcf.gz",
        "data/work/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz"
    output:
        "data/final/{tumor}/{tumor}.vardict.germline.vcf.gz",
        "data/final/{tumor}/{tumor}.vardict.germline.final.vcf.gz"
    shell:
        """
        bcftools view -W=tbi -Oz -o {output[0]} {input[0]}
        bcftools view -W=tbi -Oz -o {output[1]} {input[1]}
        """

rule vardict_unpaired_main:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "data/work/{sample}/vardict/unpaired_variants.vcf.gz"
    params:
        init="data/work/{sample}/vardict/output2.vcf.gz",
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed'],
        path="$HOME/software/VarDictJava/VarDict",
        AF_THR=0.01
    threads:
        4
    shell:
        """
        VarDict -G {params.ref} -f {params.AF_THR} -N {wildcards.sample} -b {input.bam} -c 1 -S 2 -E 3 -g 4 {params.bed} | {params.path}/teststrandbias.R | {params.path}/var2vcf_valid.pl -N {wildcards.sample} -f {params.AF_THR} | bgzip -c > {params.init}
        tabix -fp vcf {params.init}
        gatk UpdateVCFSequenceDictionary -V {params.init} --source-dictionary {input} --output {output}
        """

rule vardict_unpaired_final:
    input:
        "data/work/{sample}/vardict/unpaired_variants.vcf.gz"
    output:
        "data/final/{sample}/{sample}.vardict.unpaired.vcf.gz"
    shell:
        "bcftools view -W=tbi -Oz -o {output} {input}"
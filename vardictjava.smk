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

wildcard_constraints:
    work_dir=f"data/work/{config['resources']['targets_key']}"

#switch to 1filter/2filter
rule filter_applied_vardict_somatic:
    input:
        expand("data/work/{lib}/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule filter_applied_vardict:
    input:
        expand("data/work/{lib}/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{lib}/{tumor}/vardict/loh.twice_filtered.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{lib}/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule unprocessed_vardict_unpaired:
    input:
        expand("data/work/{lib}/{sample}/vardict/unpaired_variants.vcf.gz",lib=config['resources']['targets_key'],sample=SAMPLES)

rule run_VarDictJava:#First a more lenient -P val, not sure what
    input:
        unpack(paired_bams)
    output:
        "{work_dir}/{tumor}/vardict/variants.vcf.gz"
    params:
        init="{work_dir}/{tumor}/vardict/output.vcf.gz",
        bam=lambda wildcards: BAMS[wildcards.tumor],
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed'],
        path="$HOME/software/VarDictJava/VarDict",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        AF_THR=0.01
    threads:
        4
    shell:#split to snps and indels
        """
        VarDict -th {threads} -G {params.ref} -f {params.AF_THR} -N {wildcards.tumor} -b '{input.tumor}|{input.normal}' -c 1 -S 2 -E 3 -g 4 {params.bed} | {params.path}/testsomatic.R | {params.path}/var2vcf_paired.pl -N '{wildcards.tumor}|{params.normal}' -f {params.AF_THR} | bgzip -c > {params.init}
        tabix -p vcf {params.init}
        gatk UpdateVCFSequenceDictionary -V {params.init} --source-dictionary {params.bam} --output {output}
        """

rule run_unpaired_VarDictJava:#First a more lenient -P val, not sure what
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        "{work_dir}/{sample}/vardict/unpaired_variants.vcf.gz"
    params:
        init="{work_dir}/{sample}/vardict/output2.vcf.gz",
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed'],
        path="/home/bwubb/software/VarDictJava/VarDict",
        AF_THR=0.01
    shell:
        """
        VarDict -G {params.ref} -f {params.AF_THR} -N {wildcards.sample} -b {input.bam} -c 1 -S 2 -E 3 -g 4 {params.bed} | {params.path}/teststrandbias.R | {params.path}/var2vcf_valid.pl -N {wildcards.sample} -f {params.AF_THR} | bgzip -c > {params.init}
        tabix -fp vcf {params.init}
        gatk UpdateVCFSequenceDictionary -V {params.init} --source-dictionary {input} --output {output}
        """

#the part where it makes somatic is not needed any long
rule VarDict_filter:
    input:
        "{work_dir}/{tumor}/vardict/variants.vcf.gz"
    output:
        somatic="{work_dir}/{tumor}/vardict/somatic.vcf.gz",
        once="{work_dir}/{tumor}/vardict/variants.once_filtered.vcf.gz",
        twice="{work_dir}/{tumor}/vardict/variants.twice_filtered.vcf.gz"
    shell:
        """
        bcftools view -i 'INFO/STATUS==\"StrongSomatic\" || INFO/STATUS==\"LikelySomatic\"' -O z -o {output.somatic} {input}
        tabix -p vcf {output.somatic}
        bcftools filter --threads {threads} -e '((FORMAT/AF[0] * FORMAT/DP[0] < 6) && ((FORMAT/MQ[0] < 55.0 && FORMAT/NM[0] > 1.0) || (FORMAT/MQ[0] < 60.0 && FORMAT/NM[0] > 2.0) || (FORMAT/DP[0] < 10) || (QUAL < 45)))' -s filter_1 -m + -O z {input} > {output.once}
        tabix -p vcf {output.once}
        bcftools filter --threads {threads} -e 'FORMAT/AF[0] < 0.2 && FORMAT/QUAL[0] < 55 && INFO/SSF[0] > 0.06' -s filter_2 -m + -O z {output.once} > {output.twice}
        tabix -p vcf {output.twice}
        """

#Maybe filter is only good for somatic?
#What germline/loh are filtered out?
rule VarDict_split:
    input:
        "{work_dir}/{tumor}/vardict/variants.twice_filtered.vcf.gz"
    output:
        somatic="{work_dir}/{tumor}/vardict/somatic.twice_filtered.vcf.gz",
        loh="{work_dir}/{tumor}/vardict/loh.twice_filtered.vcf.gz",
        germline="{work_dir}/{tumor}/vardict/germline.twice_filtered.vcf.gz"
    shell:
        """
        bcftools view -i 'INFO/STATUS==\"StrongSomatic\" || INFO/STATUS==\"LikelySomatic\"' -O z -o {output.somatic} {input} && tabix -p vcf {output.somatic}
        bcftools view -i 'INFO/STATUS==\"StrongLOH\" || INFO/STATUS==\"LikelyLOH\"' -O z -o {output.loh} {input} && tabix -p vcf {output.loh}
        bcftools view -i 'INFO/STATUS==\"Germline\"' -O z -o {output.germline} {input} && tabix -p vcf {output.germline}
        """

rule VarDict_somatic_normalized:
    input:
        "{work_dir}/{tumor}/vardict/somatic.twice_filtered.vcf.gz"
    output:
        norm="{work_dir}/{tumor}/vardict/somatic.twice_filtered.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output.norm}
        tabix -f -p vcf {output.norm}
        """

rule VarDict_somatic_clean:
    input:
        "{work_dir}/{tumor}/vardict/somatic.twice_filtered.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/vardict/sample.name",
        clean="{work_dir}/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/vardict/temp.h.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

rule VarDict_loh_normalized:
    input:
        "{work_dir}/{tumor}/vardict/loh.twice_filtered.vcf.gz"
    output:
        norm="{work_dir}/{tumor}/vardict/loh.twice_filtered.norm.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output.norm}
        tabix -f -p vcf {output.norm}
        """

rule VarDict_loh_clean:
    input:
        "{work_dir}/{tumor}/vardict/loh.twice_filtered.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/vardict/sample.l.name",
        clean="{work_dir}/{tumor}/vardict/loh.twice_filtered.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/vardict/temp.l.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

rule VarDict_germline_normalized:
    input:
        "{work_dir}/{tumor}/vardict/germline.twice_filtered.vcf.gz"
    output:
        norm="{work_dir}/{tumor}/vardict/germline.twice_filtered.norm.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output.norm}
        tabix -f -p vcf {output.norm}
        """

rule VarDict_germline_clean:
    input:
        "{work_dir}/{tumor}/vardict/germline.twice_filtered.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/vardict/sample.g.name",
        clean="{work_dir}/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/vardict/temp.g.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

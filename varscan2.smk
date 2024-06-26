import os
import csv

#include:"./sequenza2.smk"

### INIT ###

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
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
    work_dir=f"data/work/{config['resources']['targets_key']}",
    results_dir=f"data/final/{config['project']['name']}"

rule filter_applied_varscan2:
    input:
        expand("data/work/{lib}/{tumor}/varscan2/somatic.fpfilter.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{lib}/{tumor}/varscan2/loh.fpfilter.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{lib}/{tumor}/varscan2/germline.fpfilter.norm.clean.vcf.gz",lib=config['resources']['targets_key'],tumor=PAIRS.keys())

rule pair_mpileup:
    input:
        unpack(paired_bams)
    output:
        temp("{work_dir}/{tumor}/varscan2/normal_tumor.mpileup")
    params:
        ref=config['reference']['fasta'],
        bed=config['resources']['targets_bed']
    shell:
        "samtools mpileup -ABR -f {params.ref} -l {params.bed} -o {output} -Q 20 {input.normal} {input.tumor}"

rule run_VarScan2:#OR return to mpileup.gz and zcat into VarScan2
    input:
        pileup="{work_dir}/{tumor}/varscan2/normal_tumor.mpileup",
        #CP="{work_dir}/{tumor}/sequenza/{tumor}_confints_CP.txt"
    output:
        indel="{work_dir}/{tumor}/varscan2/variants.indel.vcf",
        snp="{work_dir}/{tumor}/varscan2/variants.snp.vcf"
        #Note the none plural
        #I can use --output-snp and indel to name files
    params:
        prefix="{work_dir}/{tumor}/varscan2/variants",
        args='--min-coverage 3 --min-var-freq 0.08 --p-value 0.10 --somatic-p-value 0.05 --strand-filter 0',
        memory='16g'
        #old settings: args="--min-coverage 5 --p-value 0.98 --strand-filter 1"
    shell:
        """
        java -Xmx{params.memory} -jar $HOME/software/varscan/VarScan.v2.4.4.jar somatic {input.pileup} {params.prefix} --mpileup 1 --output-vcf 1 {params.args}
        """
    #run:
    #    with open(input["CP"],'r') as file:
    #        lines=file.read().splitlines()
    #        purity=lines[3].split('\t')[0]#3 for max
    #    shell(f"java -Xmx{params.memory} -jar $HOME/software/varscan/VarScan.v2.4.4.jar somatic {input.pileup} {params.prefix} --mpileup 1 --output-vcf 1 --tumor-purity {purity} {params.args}")

rule VarScan2_processSomatic:
    input:
        snp="{work_dir}/{tumor}/varscan2/variants.snp.vcf",
        indel="{work_dir}/{tumor}/varscan2/variants.indel.vcf"
    output:
        "{work_dir}/{tumor}/varscan2/somatic.snp.vcf",
        "{work_dir}/{tumor}/varscan2/somatic.indel.vcf",
        "{work_dir}/{tumor}/varscan2/loh.snp.vcf",
        "{work_dir}/{tumor}/varscan2/loh.indel.vcf",
        "{work_dir}/{tumor}/varscan2/germline.snp.vcf",
        "{work_dir}/{tumor}/varscan2/germline.indel.vcf"
    params:
        bam=lambda wildcards: BAMS[wildcards.tumor],
        snp="{work_dir}/{tumor}/varscan2/variants.snp.2.vcf",
        indel="{work_dir}/{tumor}/varscan2/variants.indel.2.vcf"
        #--min-tumor-freq - Minimum variant allele frequency in tumor [0.10]
        #--max-normal-freq - Maximum variant allele frequency in normal [0.05]
        #--p-value - P-value for high-confidence calling [0.07]
    shell:
        """
        gatk UpdateVCFSequenceDictionary -V {input.snp} --source-dictionary {params.bam} --output {params.snp}
        gatk UpdateVCFSequenceDictionary -V {input.indel} --source-dictionary {params.bam} --output {params.indel}
        java -jar ~/software/varscan/VarScan.v2.4.4.jar processSomatic {params.snp}
        java -jar ~/software/varscan/VarScan.v2.4.4.jar processSomatic {params.indel}
        mv {wildcards.work_dir}/{wildcards.tumor}/varscan2/variants.snp.2.Somatic.vcf {wildcards.work_dir}/{wildcards.tumor}/varscan2/somatic.snp.vcf
        mv {wildcards.work_dir}/{wildcards.tumor}/varscan2/variants.indel.2.Somatic.vcf {wildcards.work_dir}/{wildcards.tumor}/varscan2/somatic.indel.vcf
        mv {wildcards.work_dir}/{wildcards.tumor}/varscan2/variants.snp.2.LOH.vcf {wildcards.work_dir}/{wildcards.tumor}/varscan2/loh.snp.vcf
        mv {wildcards.work_dir}/{wildcards.tumor}/varscan2/variants.indel.2.LOH.vcf {wildcards.work_dir}/{wildcards.tumor}/varscan2/loh.indel.vcf
        mv {wildcards.work_dir}/{wildcards.tumor}/varscan2/variants.snp.2.Germline.vcf {wildcards.work_dir}/{wildcards.tumor}/varscan2/germline.snp.vcf
        mv {wildcards.work_dir}/{wildcards.tumor}/varscan2/variants.indel.2.Germline.vcf {wildcards.work_dir}/{wildcards.tumor}/varscan2/germline.indel.vcf
        """
        #RENAME

rule bamreadcount_somatic_regions:
    input:
        "{work_dir}/{tumor}/varscan2/somatic.snp.vcf",
        "{work_dir}/{tumor}/varscan2/somatic.indel.vcf"
    output:
        "{work_dir}/{tumor}/varscan2/somatic.snp.regions",
        "{work_dir}/{tumor}/varscan2/somatic.indel.regions"
    script:
        "bam-readcount_regions.py"

rule bamreadcount_somatic_readcounts:
    input:
        bam=lambda wildcards: BAMS[wildcards.tumor],
        snp="{work_dir}/{tumor}/varscan2/somatic.snp.regions",
        indel="{work_dir}/{tumor}/varscan2/somatic.indel.regions"
    output:
        snp="{work_dir}/{tumor}/varscan2/somatic.snp.readcounts",
        indel="{work_dir}/{tumor}/varscan2/somatic.indel.readcounts"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        $HOME/software/bam-readcount/bin/bam-readcount -b 15 -q 1 -w 1 -f {params.ref} -l {input.snp} {input.bam} > {output.snp}
        $HOME/software/bam-readcount/bin/bam-readcount -i -b 15 -q 1 -w 1 -f {params.ref} -l {input.indel} {input.bam} > {output.indel}
        """

#the precision of variant and mutation calling by removing artifacts associated with short-read alignment.
#-For somatic mutations, generate bam-readcounts with the Tumor BAM. For LOH and Germline, generate readcounts with the Normal BAM
#-For de novo mutations (trio calling), generate readcounts with the child BAM.
#The filter requires the bam-readcount utility: https://github.com/genome/bam-readcount

rule VarScan2_somatic_fpfilter:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/somatic.snp.vcf",
        indel_vcf="{work_dir}/{tumor}/varscan2/somatic.indel.vcf",
        snp_readcount="{work_dir}/{tumor}/varscan2/somatic.snp.readcounts",
        indel_readcount="{work_dir}/{tumor}/varscan2/somatic.indel.readcounts"
    output:
        snp_out="{work_dir}/{tumor}/varscan2/somatic.fpfilter.snp.out",
        indel_out="{work_dir}/{tumor}/varscan2/somatic.fpfilter.indel.out"
    shell:
        """
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.snp_vcf} {input.snp_readcount} --output-file {output.snp_out} --keep-failures
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.indel_vcf} {input.indel_readcount} --output-file {output.indel_out} --keep-failures
        """
        ##
        #For VarScan fpfilter (--dream-3-settings):
        #--min-var-count = 3
        #--min-var-count-lc = 1
        #--min-strandedness = 0
        #--min-var-basequal = 30
        #--min-ref-readpos = 0.20
        #--min-ref-dist3 = 0.20
        #--min-var-readpos = 0.15
        #--min-var-dist3 = 0.15
        #--max-rl-diff = 0.05
        #--max-mapqual-diff = 10
        #--min-ref-mapqual = 20
        #--min-var-mapqual = 30
        #--max-var-mmqs = 100
        #--max-ref-mmqs = 50

rule Varscan2_somatic_fix_output:
    input:
        "{work_dir}/{tumor}/varscan2/somatic.fpfilter.snp.out",
        "{work_dir}/{tumor}/varscan2/somatic.fpfilter.indel.out"
    output:
        "{work_dir}/{tumor}/varscan2/somatic.fpfilter.snps.vcf",
        "{work_dir}/{tumor}/varscan2/somatic.fpfilter.indels.vcf"
    script:
        'fix_varscan2_out.py'
        #revisit

rule VarScan2_somatic_merge:
    input:
        snp="{work_dir}/{tumor}/varscan2/somatic.fpfilter.snps.vcf",
        indel="{work_dir}/{tumor}/varscan2/somatic.fpfilter.indels.vcf"
    output:
        "{work_dir}/{tumor}/varscan2/somatic.fpfilter.vcf.gz"
    params:
        bgzsnp="{work_dir}/{tumor}/varscan2/somatic.snps.fpfilter.vcf.gz",
        bgzindel="{work_dir}/{tumor}/varscan2/somatic.indels.fpfilter.vcf.gz"
    shell:
        """
        bgzip -c {input.snp} > {params.bgzsnp}
        tabix -f -p vcf {params.bgzsnp}
        bgzip -c {input.indel} > {params.bgzindel}
        tabix -f -p vcf {params.bgzindel}
        bcftools concat -a {params.bgzsnp} {params.bgzindel}| bcftools sort -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_somatic_normalized:
    input:
        "{work_dir}/{tumor}/varscan2/somatic.fpfilter.vcf.gz"
    output:
        norm="{work_dir}/{tumor}/varscan2/somatic.fpfilter.norm.vcf.gz"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output.norm}
        tabix -f -p vcf {output.norm}
        """

rule VarScan2_somatic_clean:
    input:
        "{work_dir}/{tumor}/varscan2/somatic.fpfilter.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/varscan2/sample.name",
        clean="{work_dir}/{tumor}/varscan2/somatic.fpfilter.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/varscan2/temp.h.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

rule bamreadcount_loh_regions:
    input:
        "{work_dir}/{tumor}/varscan2/loh.snp.vcf",
        "{work_dir}/{tumor}/varscan2/loh.indel.vcf"
    output:
        "{work_dir}/{tumor}/varscan2/loh.snp.regions",
        "{work_dir}/{tumor}/varscan2/loh.indel.regions"
    script:
        "bam-readcount_regions.py"

rule bamreadcount_loh_recounts:
    input:
        bam=lambda wildcards: BAMS[PAIRS[wildcards.tumor]],
        snp="{work_dir}/{tumor}/varscan2/loh.snp.regions",
        indel="{work_dir}/{tumor}/varscan2/loh.indel.regions"
    output:
        snp="{work_dir}/{tumor}/varscan2/loh.snp.readcounts",
        indel="{work_dir}/{tumor}/varscan2/loh.indel.readcounts"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        $HOME/software/bam-readcount/bin/bam-readcount -b 15 -q 1 -w 1 -f {params.ref} -l {input.snp} {input.bam} > {output.snp}
        $HOME/software/bam-readcount/bin/bam-readcount -i -b 15 -q 1 -w 1 -f {params.ref} -l {input.indel} {input.bam} > {output.indel}
        """

rule VarScan2_loh_fpfilter:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/loh.snp.vcf",
        indel_vcf="{work_dir}/{tumor}/varscan2/loh.indel.vcf",
        snp_readcount="{work_dir}/{tumor}/varscan2/loh.snp.readcounts",
        indel_readcount="{work_dir}/{tumor}/varscan2/loh.indel.readcounts"
    output:
        snps="{work_dir}/{tumor}/varscan2/loh.snps.fpfilter.vcf.gz",
        indels="{work_dir}/{tumor}/varscan2/loh.indels.fpfilter.vcf.gz"
    params:
        snp_out="{work_dir}/{tumor}/varscan2/loh.snps.fpfilter.vcf",
        indel_out="{work_dir}/{tumor}/varscan2/loh.indels.fpfilter.vcf"
    shell:
        """
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.snp_vcf} {input.snp_readcount} --output-file {params.snp_out}
        bgzip {params.snp_out}
        tabix -f -p vcf {output.snps}
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.indel_vcf} {input.indel_readcount} --output-file {params.indel_out}
        bgzip {params.indel_out}
        tabix -f -p vcf {output.indels}
        """

rule VarScan2_loh_merge:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/loh.snps.fpfilter.vcf.gz",
        indel_vcf="{work_dir}/{tumor}/varscan2/loh.indels.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.vcf.gz"
    shell:
        """
        bcftools concat -a {input} | bcftools sort -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_loh_normalized:
    input:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.norm.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_loh_clean:
    input:
        "{work_dir}/{tumor}/varscan2/loh.fpfilter.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/varscan2/sample.loh.name",
        clean="{work_dir}/{tumor}/varscan2/loh.fpfilter.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/varscan2/temp.loh.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

rule bamreadcount_germline_regions:
    input:
        "{work_dir}/{tumor}/varscan2/germline.snp.vcf",
        "{work_dir}/{tumor}/varscan2/germline.indel.vcf"
    output:
        "{work_dir}/{tumor}/varscan2/germline.snp.regions",
        "{work_dir}/{tumor}/varscan2/germline.indel.regions"
    script:
        "bam-readcount_regions.py"

rule bamreadcount_germline_recounts:
    input:
        bam=lambda wildcards: BAMS[PAIRS[wildcards.tumor]],
        snp="{work_dir}/{tumor}/varscan2/germline.snp.regions",
        indel="{work_dir}/{tumor}/varscan2/germline.indel.regions"
    output:
        snp="{work_dir}/{tumor}/varscan2/germline.snp.readcounts",
        indel="{work_dir}/{tumor}/varscan2/germline.indel.readcounts"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        $HOME/software/bam-readcount/bin/bam-readcount -b 15 -q 1 -w 1 -f {params.ref} -l {input.snp} {input.bam} > {output.snp}
        $HOME/software/bam-readcount/bin/bam-readcount -i -b 15 -q 1 -w 1 -f {params.ref} -l {input.indel} {input.bam} > {output.indel}
        """

rule VarScan2_germline_fpfilter:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/germline.snp.vcf",
        indel_vcf="{work_dir}/{tumor}/varscan2/germline.indel.vcf",
        snp_readcount="{work_dir}/{tumor}/varscan2/germline.snp.readcounts",
        indel_readcount="{work_dir}/{tumor}/varscan2/germline.indel.readcounts"
    output:
        snps="{work_dir}/{tumor}/varscan2/germline.snps.fpfilter.vcf.gz",
        indels="{work_dir}/{tumor}/varscan2/germline.indels.fpfilter.vcf.gz"
    params:
        snp_out="{work_dir}/{tumor}/varscan2/germline.snps.fpfilter.vcf",
        indel_out="{work_dir}/{tumor}/varscan2/germline.indels.fpfilter.vcf"
    shell:
        """
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.snp_vcf} {input.snp_readcount} --output-file {params.snp_out}
        bgzip {params.snp_out}
        tabix -f -p vcf {output.snps}
        java -jar $HOME/software/varscan/VarScan.v2.4.4.jar fpfilter {input.indel_vcf} {input.indel_readcount} --output-file {params.indel_out}
        bgzip {params.indel_out}
        tabix -f -p vcf {output.indels}
        """

rule VarScan2_germline_merge:
    input:
        snp_vcf="{work_dir}/{tumor}/varscan2/germline.snps.fpfilter.vcf.gz",
        indel_vcf="{work_dir}/{tumor}/varscan2/germline.indels.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.vcf.gz"
    shell:
        """
        bcftools concat -a {input} | bcftools sort -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_germline_normalized:
    input:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.vcf.gz"
    output:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.norm.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        ref=config['reference']['fasta']
    shell:
        """
        bcftools norm -m-both {input} | bcftools norm -f {params.ref} -O z -o {output}
        tabix -f -p vcf {output}
        """

rule VarScan2_germline_clean:
    input:
        "{work_dir}/{tumor}/varscan2/germline.fpfilter.norm.vcf.gz"
    output:
        name="{work_dir}/{tumor}/varscan2/sample.g.name",
        clean="{work_dir}/{tumor}/varscan2/germline.fpfilter.norm.clean.vcf.gz"
    params:
        regions=config['resources']['targets_bedgz'],
        fai=f"{config['reference']['fasta']}.fai",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        vcf=temp("{work_dir}/{tumor}/varscan2/temp.g.vcf.gz")
    shell:
        """
        echo -e "TUMOR\\ttumor\\nNORMAL\\tnormal" > {output.name}
        echo -e "{wildcards.tumor}\\ttumor\\n{params.normal}\\tnormal" >> {output.name}

        bcftools reheader -f {params.fai} -s {output.name} -o {params.vcf} {input}
        bcftools index {params.vcf}

        bcftools view -e 'ALT~\"*\"' -R {params.regions} {params.vcf} | bcftools sort -O z -o {output.clean}
        tabix -f -p vcf {output.clean}
        """

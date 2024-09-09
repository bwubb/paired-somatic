import errno
import os

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno==errno.EEXIST and os.path.isdir(path):
            pass

#Pair Table.
with open(config['project']['pair_table'],'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

def normal_cnn(wildcards):
    #Provide the *.targetcoverage.cnn and *.antitargetcoverage.cnn files created by the coverage command:
    cnn=[]
    for n in list(PAIRS.values()):
        if f"data/work/{config['resources']['targets_key']}/{n}/cnvkit/targetcoverage.cnn" not in cnn:
            cnn.append(f"data/work/{config['resources']['targets_key']}/{n}/cnvkit/{n}.targetcoverage.cnn")
            cnn.append(f"data/work/{config['resources']['targets_key']}/{n}/cnvkit/{n}.antitargetcoverage.cnn")
    return cnn

rule cnvkit_all:
    input:
        expand('data/work/{targets}/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.report.csv',targets=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{targets}/{tumor}/cnvkit/{tumor}.cnvkit.hrd.txt",targets=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand('data/work/{targets}/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.report.csv',targets=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{targets}/{tumor}/purecn/{tumor}.purecn.hrd.txt",targets=config['resources']['targets_key'],tumor=PAIRS.keys())


rule cnvkit_access:
    output:
        "/home/bwubb/resources/Bed_files/cnvkit.access-excludes.GRCh38.bed"
    params:
        ref=config['reference']['fasta'],
        exclude="/home/bwubb/resources/Bed_files/SV_blacklist.10xGenomics.GRCh38.bed"
    shell:
        """
        cnvkit.py access {params.ref} -x {params.exclude} -o {output}
        """
        #This had all the KT*** stuff and is causing errors in PureCN...

#rule cnvkit_target:
#    input:
#        bed=config['resources']['targets_bed'],
#        ref="/home/bwubb/resources/refGene/refFlat.grch38.txt"
#    output:
#        bed="{targets}.cnv-targets.bed"


#I think this is supose to use the output from access as -g
#rule cnvkit_antitarget:
#    input:
#        "/home/bwubb/resources/Bed_files/cnvkit.access-excludes.GRCh38.bed"
#    output:
#        "/home/bwubb/resources/Bed_files/{targets}.GRCh38.antitargets.bed"
#    params:
#        bed=config['resources']['targets_bed']
#    shell:
#        """
#        cnvkit.py antitarget {params.bed} -g {input} -o {output}
#        """

#From documentation: "If multiple BAMs are given, use the BAM with median file size.""
rule cnvkit_autobin:
    input:
        access="/home/bwubb/resources/Bed_files/cnvkit.access-excludes.GRCh38.bed",
        bams=BAMS.values()
    output:
        targets="data/cnvkit/{project}.{targets}.cnvkit-targets.bed",
        antitargets="data/cnvkit/{project}.{targets}.cnvkit-antitargets.bed"
    params:
        bed=config['resources']['targets_bed'],
        ref_flat="/home/bwubb/resources/refGene/refFlat.grch38.txt"
    shell:
        """
        cnvkit.py autobin {input.bams} -t {params.bed} -g {input.access} --annotate {params.ref_flat} --short-names --target-output-bed {output.targets} --antitarget-output-bed {output.antitargets}
        """
#This one is confusing, it appears to be a different way to generate the antitarget.bed, and maybe a targets.bed?
#But taking bam files into account, but its not clear if normals AND tumors should be used or what.
#It might replace some stuff, but there is no -o and I can parallel coverage.


rule cnvkit_target_coverage:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample],
        bed=lambda wildcards: f"data/cnvkit/{config['project']['name']}.{wildcards.targets}.cnvkit-targets.bed"
    output:
        "data/work/{targets}/{sample}/cnvkit/{sample}.targetcoverage.cnn"
    shell:
        """
        cnvkit.py coverage {input.bam} {input.bed} -o {output}
        """

rule cnvkit_antitarget_coverage:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample],
        bed=lambda wildcards: f"data/cnvkit/{config['project']['name']}.{wildcards.targets}.cnvkit-antitargets.bed"
    output:
        'data/work/{targets}/{sample}/cnvkit/{sample}.antitargetcoverage.cnn'
    shell:
        """
        cnvkit.py coverage {input.bam} {input.bed} -o {output}
        """

#To analyze a cohort sequenced on a single platform, we recommend combining all normal samples into a pooled reference,
#even if matched tumor-normal pairs were sequenced – our benchmarking showed that a pooled reference performed slightly
#better than constructing a separate reference for each matched tumor-normal pair. Furthermore, even matched normals
#from a cohort sequenced together can exhibit distinctly different copy number biases
#(see Plagnol et al. 2012 and Backenroth et al. 2014);
#reusing a pooled reference across the cohort provides some consistency to help diagnose such issues.


rule cnvkit_reference:
    input:
        normal_cnn
    output:
        "data/cnvkit/reference.cnn"
    params:
        ref=config['reference']['fasta']
    shell:
        """
        cnvkit.py reference {input} --fasta {params.ref} -o {output}
        """
#file names needs to be sample.{,anti}targetcoverage.cnn
rule cnvkit_fix:
    input:
        targets="data/work/{targets}/{tumor}/cnvkit/{tumor}.targetcoverage.cnn",
        antitargets="data/work/{targets}/{tumor}/cnvkit/{tumor}.antitargetcoverage.cnn",
        reference="data/cnvkit/reference.cnn"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.cnr"
    shell:
        """
        cnvkit.py fix {input.targets} {input.antitargets} {input.reference} -o {output}
        """

rule cnvkit_input_vcf:
    input:
        "data/work/{targets}/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz"
    output:
        "data/work/{targets}/{tumor}/cnvkit/vardict.snps.clean.vcf.gz"
    shell:
        """
        bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y --type snps -f PASS -Oz -W=tbi -o {output} {input}
        """

#can use threads with -p after you have a proper cluser config
rule cnvkit_segment:
    input:
        cnr="data/work/{targets}/{tumor}/cnvkit/{tumor}.cnr",
        vcf="data/work/{targets}/{tumor}/cnvkit/vardict.snps.clean.vcf.gz"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.cns"
    params:
        tumor="{tumor}",
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        cnvkit.py segment {input.cnr} -v {input.vcf} -i {params.tumor} -n {params.normal} -o {output}
        """
#Native VarScan2 calls lack FORMAT:AF

#export to seg
rule cnvkit_export_seg:
    input:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.cns"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.seg"
    shell:
        """
        cnvkit.py export seg {input} -o {output}
        """

rule purecn_input_vcf:
    input:
        "data/work/{targets}/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz",
        "data/work/{targets}/{tumor}/vardict/germline.twice_filtered.norm.clean.vcf.gz"
    output:
        "data/work/{targets}/{tumor}/purecn/vardict.snps.clean.vcf.gz"
    params:
        vcf="data/work/{targets}/{tumor}/purecn/vardict.snps.vcf.gz"
    shell:
        """
        bcftools concat -a {input} | bcftools view --type snps -f PASS | bcftools sort -W=tbi -Oz -o {params.vcf}
        bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y -W=tbi -Oz -o {output} {params.vcf}
        """

rule purecn_run:
    input:
        snps="data/work/{targets}/{tumor}/purecn/vardict.snps.clean.vcf.gz",
        seg="data/work/{targets}/{tumor}/cnvkit/{tumor}.seg",
        cnr="data/work/{targets}/{tumor}/cnvkit/{tumor}.cnr"
    output:
        "data/work/{targets}/{tumor}/purecn/{tumor}.csv",
        "data/work/{targets}/{tumor}/purecn/{tumor}_loh.csv"
    params:
        seg="data/work/{targets}/{tumor}/purecn/input.seg",
        cnr="data/work/{targets}/{tumor}/purecn/input.cnr",
        outdir="data/work/{targets}/{tumor}/purecn"
    shell:
        """
        awk '$2 !~ /^[GK]/' {input.seg} > {params.seg}
        awk '$1 !~ /^[GK]/' {input.cnr} > {params.cnr}

        Rscript /home/bwubb/software/PureCN/inst/extdata/PureCN.R \
        --out {params.outdir} \
        --sampleid {wildcards.tumor} \
        --tumor {params.cnr} \
        --seg-file {params.seg} \
        --vcf {input.snps} \
        --genome grch38 \
        --fun-segmentation Hclust \
        --force
        """

rule purecn_to_bed:
    input:
        "data/work/{targets}/{tumor}/purecn/{tumor}_loh.csv"
    output:
        "data/work/{targets}/{tumor}/purecn/{tumor}_loh.bed"
    shell:
        """
        python cnv_to_bed.py -c purecn {input}
        """

#I would like to add --sex
rule purecn_HRDex:
    input:
        "data/work/{targets}/{tumor}/purecn/{tumor}_loh.bed"
    output:
        "data/work/{targets}/{tumor}/purecn/{tumor}.purecn.hrd.txt"
    params:
        tumor="{tumor}",
        build="grch38"
    shell:
        """
        Rscript runHRDex.R -i {input} -o {output} --build {params.build} --tumor {params.tumor}
        """

rule purecn_annotsv:
    input:
        "data/work/{targets}/{tumor}/purecn/{tumor}_loh.bed"
    output:
        "data/work/{targets}/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule purecn_annotsv_parser:
    input:
        "data/work/{targets}/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.tsv"
    output:
        "data/work/{targets}/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

#--ploidy requires int
#NEED sex
#remove -y? change -x to female
rule cnvkit_call:
    input:
        csv="data/work/{targets}/{tumor}/purecn/{tumor}.csv",
        vcf="data/work/{targets}/{tumor}/cnvkit/vardict.snps.clean.vcf.gz",
        cns="data/work/{targets}/{tumor}/cnvkit/{tumor}.cns"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.call.cns"
    params:
        tumor="{tumor}",
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    shell:
        """
        purity=`grep {params.tumor} {input.csv} | cut -d, -f2`

        cnvkit.py call {input.cns} -y -x male -m clonal --purity $purity -v {input.vcf} -i {params.tumor} -n {params.normal} -o {output}
        """

rule cnvkit_to_bed:
    input:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.call.cns"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.call.bed"
    shell:
        """
        python cnv_to_bed.py -c cnvkit {input}
        """

rule cnvkit_HRDex:
    input:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.call.bed"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.cnvkit.hrd.txt"
    params:
        tumor="{tumor}",
        build="grch38"
    shell:
        """
        Rscript runHRDex.R -i {input} -o {output} --build {params.build} --tumor {params.tumor}
        """

rule cnvkit_annotsv:
    input:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.call.bed"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule cnvkit_annotsv_parser:
    input:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.tsv"
    output:
        "data/work/{targets}/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

rule purecn_purity_table:
    input:
        expand("data/work/{targets}/{tumor}/purecn/{tumor}.csv",targets=config['resources']['targets_key'],tumor=PAIRS.keys())
    output:
        "purecn_purity.table","purecn_purity.table"
    run:
        with open(output[0],'w') as outfile0, open(output[1],'w') as outfile1:
            for f in input:
                with open(f,'r') as infile:
                    line=infile.readline()#header
                    line=infile.readline()
                    sampleid,purity,ploidy=line.replace('"','').rstrip().split(',')[0:3]
                outfile0.write(f'{sampleid}\t{purity}\n')
                outfile1.write(f'{sampleid}\t{ploidy}\n')

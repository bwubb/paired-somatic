import csv

#ASCAT proto-pipeline

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()

#with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
#ASCAT has dropped several samples.
#This probably wont be an issue in the future and should be changed back.
#We included pairs with germline <30X
with open('ascat_pair.table','r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config.get('project',{}).get('bam_table','bams.table'),'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

#I need FULL PATHS to bams
def bam_input(wildcards):
    bams=[]
    for t,n in PAIRS.items():
        bams+=[BAMS[t],BAMS[n]]
    return list(set(bams))

def paired_bams(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

rule ascat_all:
    input:
        expand("data/work/{targets}/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.report.csv",targets=config['resources']['targets_key'],tumor=PAIRS.keys()),
        expand("data/work/{targets}/{tumor}/ascat/{tumor}.ascat.hrd.txt",targets=config['resources']['targets_key'],tumor=PAIRS.keys())

rule ASCAT_write_worksheet:
    input:
        bam_input
    output:
        "data/work/ASCAT/worksheet.tsv"
    run:
        with open(output[0],'w') as f:
            writer=csv.writer(f,delimiter="\t")
            writer.writerow(['Patient_ID','Normal_ID','Normal_file','Gender'])
            for t,n in PAIRS.items():
                writer.writerow([t,n,BAMS[n],"XX"])


#output not directly used, picked X out of all any chr number.
rule ASCAT_prepareTargetedSeq:
    input:
        worksheet="data/work/ASCAT/worksheet.tsv"
    output:
        "data/work/ASCAT/alleleData/Cleaned/alleleData_chrX.txt",
        "data/work/ASCAT/alleleData/Cleaned/loci_chrX.txt"
    params:
        outdir="data/work/ASCAT",
        alleles="/home/bwubb/projects/021-POSH_FFPE/G1000_allelesAll_hg38/G1000_alleles_hg38_chr",
        bed=config['resources']['targets_bed']
    threads:
        8
    script:
        "ASCAT_prepareTargetedSeq.R"


rule ASCAT_prepareHTS:
    input:
        unpack(paired_bams),
        allele="data/work/ASCAT/alleleData/Cleaned/alleleData_chrX.txt",
        loci="data/work/ASCAT/alleleData/Cleaned/loci_chrX.txt"
    output:
        tumorLogR="data/work/{targets}/{tumor}/ascat/Tumor_LogR.txt",
        tumorBAF="data/work/{targets}/{tumor}/ascat/Tumor_BAF.txt",
        germlineLogR="data/work/{targets}/{tumor}/ascat/Germline_LogR.txt",
        germlineBAF="data/work/{targets}/{tumor}/ascat/Germline_BAF.txt"
    params:
        bed=config["resources"]["targets_bed"],
        alleles="data/work/ASCAT/alleleData/Cleaned/alleleData_chr",
        loci="data/work/ASCAT/alleleData/Cleaned/loci_chr",
        tumor='{tumor}',
        normal=lambda wildcards: PAIRS[wildcards.tumor]
    threads:
        1
    script:
        "ASCAT_prepareHTS.R"

rule ASCAT_runAscat:
    input:
        tumorLogR="data/work/{targets}/{tumor}/ascat/Tumor_LogR.txt",
        tumorBAF="data/work/{targets}/{tumor}/ascat/Tumor_BAF.txt",
        germlineLogR="data/work/{targets}/{tumor}/ascat/Germline_LogR.txt",
        germlineBAF="data/work/{targets}/{tumor}/ascat/Germline_BAF.txt"
    output:
        segments="data/work/{targets}/{tumor}/ascat/{tumor}.segments.txt",
        rdata="data/work/{targets}/{tumor}/ascat/ASCAT_objects.Rdata"
    params:
        outdir="data/work/{targets}/{tumor}/ascat",
        gc_file="GC_G1000_hg38.txt"
    script:
        "ASCAT_runAscat.R"

rule ASCAT_to_bed:
    input:
        "data/work/{targets}/{tumor}/ascat/{tumor}.segments.txt"
    output:
        "data/work/{targets}/{tumor}/ascat/{tumor}.segments.bed"
    shell:
        "python cnv_to_bed.py -c ascat {input}"

rule ASCAT_annotsv:
    input:
        "data/work/{targets}/{tumor}/ascat/{tumor}.segments.bed"
    output:
        "data/work/{targets}/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule ASCAT_annotsv_parser:
    input:
        "data/work/{targets}/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.tsv"
    output:
        "data/work/{targets}/{tumor}/ascat/{tumor}.ascat.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

rule ascat_HRDex:
    input:
        "data/work/{targets}/{tumor}/ascat/{tumor}.segments.bed"
    output:
        "data/work/{targets}/{tumor}/ascat/{tumor}.ascat.hrd.txt"
    params:
        tumor="{tumor}",
        build="grch38"
    shell:
        """
        Rscript runHRDex.R -i {input} -o {output} --build {params.build} --tumor {params.tumor}
        """

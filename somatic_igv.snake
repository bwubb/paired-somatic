import csv
import os
import errno
from collections import defaultdict

#Read tumor/normals names/bams
#Read mutation file
#Make intervals
#Make bed +500 each side
#IGV thing to make bat, -nosnap, Can I make my own with better settings???
#Run xvfb with modified settings
#???
#Profit

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno==errno.EEXIST and os.path.isdir(path):
            pass

def paired_bams(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

def get_position(wildcards):
    position=MUTS[wildcards.mut]
    chr,start=position.split(':')
    return f'{chr}:{int(start)-100}-{int(start)+100}'

#TumorID', 'MayoBrca2Br34
MUTS={}
SAMPLE_MUT=defaultdict(list)
with open('SecondaryMuts_in_question.csv','r') as file:
    #dict['GENE_pA123B']=1:23345
    reader=csv.DictReader(file,delimiter=',')
    for row in reader:
        if row['AAChange.refGene']=='.':
            X=row['NTChange.refGene'].replace('.','')
            #If not AAChange.refGene use NTChange.refGene
            #Need to check for special characters
        else:
            X=row['AAChange.refGene'].replace('.','')
        key=f"{row['Gene.refGene']}_{X}"
        
        MUTS[key]=f"{row['Chr']}:{row['Start']}"
        SAMPLE_MUT[row['TumorID']].append(key)
        #mkdir_p(f"analysis/work/{row['TumorID']}/somatic_variants/igv_png")
        with open(f"analysis/work/{row['TumorID']}/somatic_variants/igv_png/{key}.bed",'w') as bed_file:
            bed_row=[row['Chr']]
            bed_row.append(str(int(row['Start'])-100))
            bed_row.append(str(int(row['Start'])+100))
            bed_row.append(key)
            out="\t".join(bed_row)
            bed_file.write(f'{out}\n')

TARGET_FILES=[]
for tumor,muts in SAMPLE_MUT.items():
    for mut in muts:
        TARGET_FILES.append(f'analysis/work/{tumor}/somatic_variants/igv_png/{tumor}_{mut}.png')

BAMS={}
with open('bams.table','r') as file:
    for line in file:
        sample,bam=line.rstrip().split('\t')
        BAMS[sample]=bam

PAIRS={}
with open('pairs.tsv','r') as file:
    for line in file:
        tumor,normal=line.rstrip().split('\t')
        PAIRS[tumor]=normal

rule all:
    input:
        TARGET_FILES

rule Mutect2_bamout:
    input:
        unpack(paired_bams)
    output:
        bam="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.bam",
        vcf="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2_activeregion.vcf.gz"
    params:
        ref='~/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta',
        tumor=lambda wildcards: wildcards.tumor,
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        intervals=lambda wildcards: MUTS[wildcards.mut]
    log:
        "analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.log"
    shell:
        """
        gatk Mutect2 -R {params.ref} -I {input.tumor} -I {input.normal} \
        -tumor {params.tumor} -normal {params.normal} -L {params.intervals} \
        -ip 1000 -bamout {output.bam} -O {output.vcf}
        
        samtools index {output.bam}
        """

#rule annotate_activeregion_vcf:
#    input:
#        "analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2_activeregion.vcf.gz"
#    output:
#        "analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.vcf.gz"
#    shell:
#        """
#        bcftools view
#        """

rule write_bat:
    input:
        bam="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.bam"
    output:
        bat="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.bat"
    params:
        dir="analysis/work/{tumor}/somatic_variants/igv_png/",
        name="{tumor}_{mut}",
        position=get_position
    run:
        with open(output['bat'],'w') as file:
            file.write(f"new\n")
            file.write(f"genome b37\n")
            file.write(f"snapshotDirectory {params['dir']}\n")
            file.write(f"load {input['bam']}\n")
            file.write(f"maxPanelHeight 2000\n")
            file.write(f"goto {params['position']}\n")
            file.write(f"group strand\n")
            file.write(f"snapshot {params['name']}\n")
            file.write(f"exit\n")

rule run_bat:
    input:
        bat="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.bat"
    output:
        bat="analysis/work/{tumor}/somatic_variants/igv_png/{tumor}_{mut}.png"
    shell:
        "xvfb-run --auto-servernum --server-num=1 --server-args='-screen 0, 1900x1200x24' java -Xmx10240m -jar /usr/local/software/IGV_2.4.4/igv.jar -b {input}"
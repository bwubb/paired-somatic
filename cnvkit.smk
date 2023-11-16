import errno
import os

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno==errno.EEXIST and os.path.isdir(path):
            pass

def tumors_and_normals(wildcards):
    tumor_bams=[]
    for t in TUMORS:
        tumor_bams.append(BAMS[t])
    norm_bams=[]
    for n in NORMALS:
        norm_bams.append(BAMS[n])
    #print(' '.join(tumor_bams))
    #print(' '.join(norm_bams))
    return {"tumors":TUMORS,"normals":NORMALS}

with open(config['project']['bam_list'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

with open('tumors.list','r') as file:
    TUMORS=file.read().splitlines()

with open('normals.list','r') as file:
    NORMALS=file.read().splitlines()

rule all:
    input:
        'reference.cnn'

rule cnvkit_batch:
    input:
        unpack(tumors_and_normals)
    output:
        reference.cnn
    params:
        targets="/home/bwubb/resources/Bed_files/SureSelect-Exon_v5.S04380110.cnvkit.targets.bed",
        antitargets="/home/bwubb/resources/Bed_files/SureSelect-Exon_v5.S04380110.cnvkit.antitargets.bed",
        ref="/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta",
        refFlat="/home/bwubb/resources/refFlat.GRCh37.txt",
        access="/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.access-excludes.bed"
    threads:
        8
    shell:
        # From baits and tumor/normal BAMs
        '''
        cnvkit.py batch {input.tumors} --normal {input.normals} \
        --targets {params.cnv_targets} --annotate {params.refFlat} \
        --fasta {params.ref} --antitargets {params.cnv_antitargets} \
        --access {params.access} -p {threads} \
        --output-reference reference.cnn --output-dir results/ \
        --diagram --scatter
        '''

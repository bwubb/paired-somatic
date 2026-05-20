import csv
import os

# Layout: analysis/work/{sample}/somatic_variants/igv_png/{mut}.(bwa|mutect2).bam[.bai]
# Set IGV jar in rule run_bat if your install path differs.

# (tumor_or_sample, mut_name) -> MUTS entry dict with name, range, samples
MUT_LOOKUP={}

def get_position(wildcards):
    entry=MUT_LOOKUP[f'{wildcards.tumor}_{wildcards.mut}']
    return entry['range']

def igv_snapshot_bam(wildcards):
    """Mutect2 excerpt when tumor/normal pair exists; else samtools bwa excerpt."""
    entry=MUT_LOOKUP[f'{wildcards.tumor}_{wildcards.mut}']
    base=f"analysis/work/{wildcards.tumor}/somatic_variants/igv_png/{wildcards.mut}"
    if len(entry['samples'])>=2:
        return f"{base}.mutect2.bam"
    return f"{base}.bwa.bam"

def igv_snapshot_bai(wildcards):
    entry=MUT_LOOKUP[f'{wildcards.tumor}_{wildcards.mut}']
    base=f"analysis/work/{wildcards.tumor}/somatic_variants/igv_png/{wildcards.mut}"
    if len(entry['samples'])>=2:
        return f"{base}.mutect2.bam.bai"
    return f"{base}.bwa.bam.bai"

def paired_bams(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

MUTS={}
RANGES={}
with open(config['input'],'r') as file:
    reader=csv.DictReader(file,delimiter='\t')
    for n,row in enumerate(reader):
        k=f'n{n}'
        if 'Tumor.ID' in row.keys():
            samples=[row['Tumor.ID'],row['Normal.ID']]
            for sid in samples:
                os.makedirs(f'analysis/work/{sid}/somatic_variants/igv_png',exist_ok=True)
        else:
            samples=[row['Sample.ID']]
            os.makedirs(f'analysis/work/{row["Sample.ID"]}/somatic_variants/igv_png',exist_ok=True)
        chr=row['Chr']
        start=row['Start']
        gene=row['Gene']
        if row.get('HGVSp',f'variant{n}')=='.':
            hgvs='splice'
        else:
            hgvs=row.get('HGVSp',f'variant{n}').replace('p.','').replace('?','')
        mut_name=f'chr{chr}-{start}.{gene}.{hgvs}'
        MUTS[k]={'name':mut_name,'range':f'{chr}:{int(start)-100}-{int(start)+100}','samples':samples}
        for sample in samples:
            bwa_path=f'analysis/work/{sample}/somatic_variants/igv_png/{mut_name}.bwa.bam'
            RANGES[bwa_path]=MUTS[k]['range']
        for sample in samples:
            MUT_LOOKUP[f'{sample}_{mut_name}']=MUTS[k]

TARGET_FILES=[]
SCREENSHOT_FILES=[]
for k,v in MUTS.items():
    mut_name=v['name']
    samples=v['samples']
    for sample in samples:
        TARGET_FILES.append(f'analysis/work/{sample}/somatic_variants/igv_png/{mut_name}.bwa.bam')
    # Mutect2 bamout is one file per variant (tumor sample path when paired)
    if len(samples)>=2:
        tumor_sample=samples[0]
        TARGET_FILES.append(f'analysis/work/{tumor_sample}/somatic_variants/igv_png/{mut_name}.mutect2.bam')
        SCREENSHOT_FILES.append(f'analysis/work/{tumor_sample}/somatic_variants/igv_png/{mut_name}.bat')
        SCREENSHOT_FILES.append(f'analysis/work/{tumor_sample}/somatic_variants/igv_png/{tumor_sample}_{mut_name}.png')
    else:
        s=samples[0]
        # no tumor/normal pair: no Mutect2; IGV loads samtools bwa excerpt only
        SCREENSHOT_FILES.append(f'analysis/work/{s}/somatic_variants/igv_png/{mut_name}.bat')
        SCREENSHOT_FILES.append(f'analysis/work/{s}/somatic_variants/igv_png/{s}_{mut_name}.png')

BAMS={}
with open('bam.table','r') as file:
    for line in file:
        sample,bam=line.rstrip().split('\t')
        BAMS[sample]=bam

PAIRS={}
with open('pair.table','r') as file:
    for line in file:
        tumor,normal=line.rstrip().split('\t')
        PAIRS[tumor]=normal

wildcard_constraints:
    mut='[^/]+'

rule all:
    input:
        TARGET_FILES+SCREENSHOT_FILES

rule samtools_bamout:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample]
    output:
        bam="analysis/work/{sample}/somatic_variants/igv_png/{mut}.bwa.bam",
        bai="analysis/work/{sample}/somatic_variants/igv_png/{mut}.bwa.bam.bai",
    params:
        range=lambda wildcards: RANGES[f"analysis/work/{wildcards.sample}/somatic_variants/igv_png/{wildcards.mut}.bwa.bam"]
    shell:
        """
        samtools view -O BAM -o {output.bam} {input.bam} {params.range}
        samtools index {output.bam}
        """

rule mutect2_bamout:
    input:
        unpack(paired_bams)
    output:
        bam="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.bam",
        bai="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.bam.bai",
        vcf="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.vcf"
    params:
        ref=config['reference']['fasta'],
        tumor=lambda wildcards: wildcards.tumor,
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        range=lambda wildcards: RANGES[f"analysis/work/{wildcards.tumor}/somatic_variants/igv_png/{wildcards.mut}.bwa.bam"]
    log:
        "analysis/work/{tumor}/somatic_variants/igv_png/{mut}.mutect2.log"
    shell:
        """
        gatk Mutect2 -R {params.ref} \
        -I {input.tumor} -I {input.normal} \
        -tumor {params.tumor} -normal {params.normal} \
        -L {params.range} \
        -ip 1000 -bamout {output.bam} -O {output.vcf}

        samtools index {output.bam}
        """

rule write_bat:
    input:
        bam=igv_snapshot_bam,
        bai=igv_snapshot_bai,
    output:
        bat="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.bat"
    params:
        snapshot_dir=lambda wildcards: f"analysis/work/{wildcards.tumor}/somatic_variants/igv_png",
        snapshot_name=lambda wildcards: f"{wildcards.tumor}_{wildcards.mut}",
        position=get_position,
    run:
        with open(output['bat'],'w') as file:
            file.write(f"new\n")
            file.write(f"genome {config['reference']['key']}\n")
            file.write(f"snapshotDirectory {params.snapshot_dir}/\n")
            file.write(f"load {input['bam']}\n")
            file.write(f"maxPanelHeight 2000\n")
            file.write(f"goto {params.position}\n")
            file.write(f"group strand\n")
            file.write(f"snapshot {params.snapshot_name}\n")
            file.write(f"exit\n")

rule run_bat:
    input:
        bat="analysis/work/{tumor}/somatic_variants/igv_png/{mut}.bat",
        bam=igv_snapshot_bam,
        bai=igv_snapshot_bai,
    output:
        png="analysis/work/{tumor}/somatic_variants/igv_png/{tumor}_{mut}.png"
    shell:
        "xvfb-run --auto-servernum --server-num=1 --server-args='-screen 0, 1900x1200x24' java -Xmx4g -jar /usr/local/software/IGV_2.4.4/igv.jar -b {input.bat}"

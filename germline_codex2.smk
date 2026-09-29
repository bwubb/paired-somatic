import os

#Germline-only CODEX2. Outdir hardcoded: data/work/codex2.
#Samples: project.germline_list; BAM paths: project.bam_table.
#BED: resources.targets_bed; genome: reference.key (GRCh38 / hg38).
#Exclude list: out of norm_index, still in coverage matrix.
#
#Optional:
#  codex2.mapq (20)
#  codex2.mappability none|getmapp (default none)
#  codex2.Kmax (10)

with open(config['project']['germline_list'],'r') as s:
    SAMPLES=[ln.strip() for ln in s if ln.strip() and not ln.startswith('#')]

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines() if line.strip())

CODEX=config.get('codex2',{})
OUTDIR='data/work/codex2'
BED=config['resources']['targets_bed']
BAM_TABLE=config['project']['bam_table']
GENOME_KEY=config['reference']['key']
if GENOME_KEY not in ('GRCh38','hg38'):
    raise ValueError(f"reference.key must be GRCh38 or hg38 for CODEX2 (got {GENOME_KEY})")
MAPQ=CODEX.get('mapq',20)
#Accept old key "mapp: ones" as alias for none
MAPP=CODEX.get('mappability',CODEX.get('mapp','none'))
if MAPP=='ones':
    MAPP='none'
if MAPP not in ('none','getmapp'):
    raise ValueError(f"codex2.mappability must be none or getmapp (got {MAPP})")
KMAX=CODEX.get('Kmax',10)

EXCLUDE=set()
_excl=CODEX.get('exclude_normals',config.get('cnvkit',{}).get('exclude_normals','cnvkit_exclude_normals.list'))
if os.path.exists(_excl):
    with open(_excl,'r') as e:
        EXCLUDE={line.strip() for line in e if line.strip() and not line.startswith('#')}

_missing=[n for n in SAMPLES if n not in BAMS]
if _missing:
    raise ValueError(f"germline_list samples missing from bam_table: {_missing}")

CONTROLS=sorted(n for n in SAMPLES if n not in EXCLUDE)
CASES=sorted(n for n in SAMPLES if n in EXCLUDE)
if not CONTROLS:
    raise ValueError("No control samples left for CODEX2 norm_index after exclude list")

def codex2_chr_list(wildcards):
    if CODEX.get('chrs'):
        return [str(c) for c in CODEX['chrs']]
    ck=checkpoints.codex2_prep.get(**wildcards)
    with open(ck.output.chrs,'r') as fh:
        return [ln.strip() for ln in fh if ln.strip()]

rule codex2_all:
    input:
        os.path.join(OUTDIR,'codex2.segments.filtered.txt'),
        expand(os.path.join(OUTDIR,'{sample}.codex2.annotsv.gene_split.report.csv'),sample=SAMPLES),
        os.path.join(OUTDIR,'samples.list'),
        os.path.join(OUTDIR,'controls.list'),
        os.path.join(OUTDIR,'cases.list')

rule codex2_segments_all:
    input:
        os.path.join(OUTDIR,'codex2.segments.filtered.txt')

rule codex2_sample_lists:
    input:
        bam_table=BAM_TABLE,
        germline=config['project']['germline_list']
    output:
        samples=os.path.join(OUTDIR,'samples.list'),
        controls=os.path.join(OUTDIR,'controls.list'),
        cases=os.path.join(OUTDIR,'cases.list')
    run:
        os.makedirs(OUTDIR,exist_ok=True)
        with open(output.samples,'w') as o:
            for n in SAMPLES:
                o.write(n+'\n')
        with open(output.controls,'w') as o:
            for n in CONTROLS:
                o.write(n+'\n')
        with open(output.cases,'w') as o:
            for n in CASES:
                o.write(n+'\n')

checkpoint codex2_prep:
    input:
        bam_table=BAM_TABLE,
        samples=os.path.join(OUTDIR,'samples.list'),
        bed=BED
    output:
        qc=os.path.join(OUTDIR,'coverageQC.csv'),
        refqc=os.path.join(OUTDIR,'ref_qc.rds'),
        N=os.path.join(OUTDIR,'library_size_factor.csv'),
        chrs=os.path.join(OUTDIR,'chrs.txt')
    params:
        outdir=OUTDIR,
        genome_key=GENOME_KEY,
        mapq=MAPQ,
        mappability=MAPP
    threads: 1
    shell:
        """
        Rscript codex2_prep.R \
          --bam-table {input.bam_table} \
          --samples {input.samples} \
          --bed {input.bed} \
          --outdir {params.outdir} \
          --genome {params.genome_key} \
          --mapq {params.mapq} \
          --mappability {params.mappability}
        """

rule codex2_chr:
    input:
        qc=os.path.join(OUTDIR,'coverageQC.csv'),
        refqc=os.path.join(OUTDIR,'ref_qc.rds'),
        N=os.path.join(OUTDIR,'library_size_factor.csv'),
        controls=os.path.join(OUTDIR,'controls.list')
    output:
        raw=os.path.join(OUTDIR,'chr{chr}.codex2.segments.txt'),
        filt=os.path.join(OUTDIR,'chr{chr}.codex2.segments.filtered.txt')
    params:
        outdir=OUTDIR,
        kmax=KMAX
    threads: 1
    shell:
        """
        Rscript codex2_run_chr.R \
          --outdir {params.outdir} \
          --chr {wildcards.chr} \
          --controls {input.controls} \
          --Kmax {params.kmax}
        """

def codex2_merge_input(wildcards):
    return expand(os.path.join(OUTDIR,'chr{chr}.codex2.segments.filtered.txt'),chr=codex2_chr_list(wildcards))

rule codex2_merge:
    input:
        codex2_merge_input
    output:
        os.path.join(OUTDIR,'codex2.segments.filtered.txt')
    params:
        outdir=OUTDIR
    shell:
        """
        Rscript codex2_merge.R --outdir {params.outdir} --output {output}
        """

rule codex2_to_bed:
    input:
        segs=os.path.join(OUTDIR,'codex2.segments.filtered.txt'),
        samples=os.path.join(OUTDIR,'samples.list')
    output:
        expand(os.path.join(OUTDIR,'{sample}.codex2.bed'),sample=SAMPLES)
    shell:
        """
        python cnv_to_bed.py -c codex2 --samples {input.samples} {input.segs}
        """

rule codex2_annotsv:
    input:
        os.path.join(OUTDIR,'{sample}.codex2.bed')
    output:
        os.path.join(OUTDIR,'{sample}.codex2.annotsv.gene_split.tsv')
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule codex2_annotsv_parser:
    input:
        os.path.join(OUTDIR,'{sample}.codex2.annotsv.gene_split.tsv')
    output:
        os.path.join(OUTDIR,'{sample}.codex2.annotsv.gene_split.report.csv')
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.sample}
        """

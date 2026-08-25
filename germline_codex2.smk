import os

#Germline-only CODEX2 (GRCh38, no "chr" prefix). No tumors.
#Negative-control normalization: known exon-del carriers excluded from norm_index
#(same spirit as cnvkit_exclude_normals), but still present in the coverage matrix for calling.
#
#Config (examples):
#  project.pair_table / project.bam_table  (reuse)
#  codex2.bed                              required Covered BED (GRCh38, no chr)
#  codex2.exclude_normals                  default cnvkit_exclude_normals.list
#  codex2.outdir                           default data/work/codex2
#  codex2.mapq                             default 20
#  codex2.mapp                             ones|getmapp (default ones)
#  codex2.Kmax                             default 10
#  codex2.chrs                             optional hard list (else from prep chrs.txt)

with open(config['project']['pair_table'],'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines() if line.strip())

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines() if line.strip())

CODEX=config.get('codex2',{})
OUTDIR=CODEX.get('outdir','data/work/codex2')
BED=CODEX.get('bed',config.get('resources',{}).get('targets_bed'))
if not BED:
    raise ValueError("Set codex2.bed (or resources.targets_bed) to a GRCh38 Covered BED without chr prefix")
MAPQ=CODEX.get('mapq',20)
MAPP=CODEX.get('mapp','ones')
KMAX=CODEX.get('Kmax',10)

EXCLUDE=set()
_excl=CODEX.get('exclude_normals',config.get('cnvkit',{}).get('exclude_normals','cnvkit_exclude_normals.list'))
if os.path.exists(_excl):
    with open(_excl,'r') as e:
        EXCLUDE={line.strip() for line in e if line.strip() and not line.startswith('#')}

NORMALS=sorted(set(PAIRS.values()))
CONTROLS=sorted(n for n in NORMALS if n not in EXCLUDE)
CASES=sorted(n for n in NORMALS if n in EXCLUDE)

_missing=[n for n in NORMALS if n not in BAMS]
if _missing:
    raise ValueError(f"normals missing from bam_table: {_missing}")
if not CONTROLS:
    raise ValueError("No control normals left for CODEX2 norm_index after exclude list")

def codex2_chr_list(wildcards):
    if CODEX.get('chrs'):
        return [str(c) for c in CODEX['chrs']]
    ck=checkpoints.codex2_prep.get(**wildcards)
    with open(ck.output.chrs,'r') as fh:
        return [ln.strip() for ln in fh if ln.strip()]

rule codex2_all:
    input:
        os.path.join(OUTDIR,'codex2.segments.filtered.txt'),
        os.path.join(OUTDIR,'controls.list'),
        os.path.join(OUTDIR,'cases.list')

rule codex2_sample_lists:
    output:
        bams=os.path.join(OUTDIR,'bams.list'),
        controls=os.path.join(OUTDIR,'controls.list'),
        cases=os.path.join(OUTDIR,'cases.list')
    run:
        os.makedirs(OUTDIR,exist_ok=True)
        with open(output.bams,'w') as o:
            for n in NORMALS:
                o.write(BAMS[n]+'\n')
        with open(output.controls,'w') as o:
            for n in CONTROLS:
                o.write(n+'\n')
        with open(output.cases,'w') as o:
            for n in CASES:
                o.write(n+'\n')

checkpoint codex2_prep:
    input:
        bams=os.path.join(OUTDIR,'bams.list'),
        bed=BED
    output:
        qc=os.path.join(OUTDIR,'coverageQC.csv'),
        refqc=os.path.join(OUTDIR,'ref_qc.rds'),
        N=os.path.join(OUTDIR,'library_size_factor.csv'),
        chrs=os.path.join(OUTDIR,'chrs.txt')
    params:
        outdir=OUTDIR,
        mapq=MAPQ,
        mapp=MAPP
    threads: 1
    shell:
        """
        Rscript codex2_prep.R {input.bams} {input.bed} {params.outdir} {params.mapq} {params.mapp}
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
        Rscript codex2_run_chr.R {params.outdir} {wildcards.chr} {input.controls} {params.kmax}
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
        Rscript codex2_merge.R {params.outdir} {output}
        """

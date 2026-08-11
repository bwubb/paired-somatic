import errno
import os

def mkdir_p(path):
    try:
        os.makedirs(path)
    except OSError as exc: # Python >2.5
        if exc.errno==errno.EEXIST and os.path.isdir(path):
            pass

# Host cnvkit.py dies on RHEL9 OpenSSL FIPS (pysam). Run via Singularity BioContainers image.
# Override with cnvkit.container in config if needed.
CNVKIT_SIF=config.get('cnvkit',{}).get('container','/home/bwubb/containers/cnvkit_0.9.10.sif')
# Shell preamble: pwd -P so we sit on /project/... not a broken /home/.../projects symlink.
# No --pwd (same DeepVariant gotcha). Relative paths work after cd "$WORKDIR".
CNVKIT_EXEC=(
    'WORKDIR="$(pwd -P)"; cd "$WORKDIR"; '
    'singularity exec '
    '--bind /project:/project '
    '--bind /home/bwubb:/home/bwubb '
    '--bind /scratch:/scratch '
    '--bind "$WORKDIR:$WORKDIR" '
    f'{CNVKIT_SIF} cnvkit.py'
)

#Pair Table.
with open(config['project']['pair_table'],'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

#Normals for pooled reference; drop carriers with known germline exon dels (one ID per line).
EXCLUDE_NORMALS=set()
if os.path.exists(config.get('cnvkit',{}).get('exclude_normals','cnvkit_exclude_normals.list')):
    with open(config.get('cnvkit',{}).get('exclude_normals','cnvkit_exclude_normals.list'),'r') as e:
        EXCLUDE_NORMALS={line.strip() for line in e if line.strip() and not line.startswith('#')}

#sample<TAB>female|male|F|M|XX|XY ; missing => female (breast cohort default)
SEX={}
if os.path.exists(config.get('project',{}).get('sex_table','sex.table')):
    with open(config.get('project',{}).get('sex_table','sex.table'),'r') as s:
        for line in s:
            if line.strip() and not line.startswith('#'):
                sample,sex=line.rstrip().split('\t')[:2]
                SEX[sample]=sex

def sample_sex(sample):
    # Default female: all-female breast cancer project; set sex_table only for exceptions.
    sex=SEX.get(sample,config.get('cnvkit',{}).get('default_sex','female')).lower()
    if sex in ('male','m','xy','1'):
        return {'cnvkit':'male','purecn':'M'}
    return {'cnvkit':'female','purecn':'F'}

def reference_normals():
    return sorted({n for n in PAIRS.values() if n not in EXCLUDE_NORMALS})

def normal_cnn(wildcards):
    #Provide the *.targetcoverage.cnn and *.antitargetcoverage.cnn files created by the coverage command:
    cnn=[]
    for n in reference_normals():
        tgt=f"data/work/{n}/cnvkit/{n}.targetcoverage.cnn"
        anti=f"data/work/{n}/cnvkit/{n}.antitargetcoverage.cnn"
        if tgt not in cnn:
            cnn.append(tgt)
            cnn.append(anti)
    return cnn

def normal_bams(wildcards):
    return [BAMS[n] for n in reference_normals()]

#Excluded normals present in this pair_table → germline exon del/dup vs pooled reference.
#Requires a matched tumor's paired VarDict (vardictjava.smk), not germline-vardict.smk paths.
NORMAL_TO_TUMOR={}
for tumor,normal in PAIRS.items():
    NORMAL_TO_TUMOR.setdefault(normal,tumor)

_orphan_excl=sorted(EXCLUDE_NORMALS - set(PAIRS.values()))
if _orphan_excl:
    print(f"WARNING: cnvkit exclude_normals not in pair_table (skipped for germline CN): {_orphan_excl}")

GERMLINE_CN=sorted(EXCLUDE_NORMALS & set(PAIRS.values()))

def germline_pair_vcf_inputs(wildcards):
    g=wildcards.germline
    if g not in NORMAL_TO_TUMOR:
        raise ValueError(
            f"{g} is in cnvkit exclude_normals / GERMLINE_CN but is not a normal in pair_table"
        )
    t=NORMAL_TO_TUMOR[g]
    return {
        'vcf': f"data/work/{t}/vardict/germline.vcf.gz",
        'name': f"data/work/{t}/vardict/sample.name",
    }

rule cnvkit_all:
    input:
        expand("data/work/{tumor}/cnvkit/{tumor}.call.cns",tumor=PAIRS.keys()),
        expand("data/work/{tumor}/purecn/{tumor}.csv",tumor=PAIRS.keys()),
        "purecn_purity.table",
        "purecn_ploidy.table",
        expand("data/work/{germline}/cnvkit/{germline}.germline.call.cns",germline=GERMLINE_CN)

#Optional downstream; needs AnnotSV / HRDex installed and maintained separately.
rule cnvkit_report_all:
    input:
        expand('data/work/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.report.csv',tumor=PAIRS.keys()),
        expand("data/work/{tumor}/cnvkit/{tumor}.cnvkit.hrd.txt",tumor=PAIRS.keys()),
        expand('data/work/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.report.csv',tumor=PAIRS.keys()),
        expand("data/work/{tumor}/purecn/{tumor}.purecn.hrd.txt",tumor=PAIRS.keys())

rule cnvkit_germline_all:
    input:
        expand("data/work/{germline}/cnvkit/{germline}.germline.call.cns",germline=GERMLINE_CN)

rule cnvkit_germline_report_all:
    input:
        expand('data/work/{germline}/cnvkit/{germline}.germline.annotsv.gene_split.report.csv',germline=GERMLINE_CN)

rule cnvkit_access:
    output:
        "/home/bwubb/resources/Bed_files/cnvkit.access-excludes.GRCh38.bed"
    params:
        ref=config['reference']['fasta'],
        exclude="/home/bwubb/resources/Bed_files/SV_blacklist.10xGenomics.GRCh38.bed",
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} access {params.ref} -x {params.exclude} -o {output}
        """
        #This had all the KT*** stuff and is causing errors in PureCN...

#From documentation: "If multiple BAMs are given, use the BAM with median file size."
#batch uses normals only for autobin; use Covered BED via cnvkit.baits_bed when set.
rule cnvkit_autobin:
    input:
        access="/home/bwubb/resources/Bed_files/cnvkit.access-excludes.GRCh38.bed",
        bams=normal_bams
    output:
        targets=f"data/cnvkit/{config['project']['name']}.{config['resources']['targets_key']}.cnvkit-targets.bed",
        antitargets=f"data/cnvkit/{config['project']['name']}.{config['resources']['targets_key']}.cnvkit-antitargets.bed"
    params:
        bed=config.get('cnvkit',{}).get('baits_bed',config['resources']['targets_bed']),
        ref_flat="/home/bwubb/resources/refGene/refFlat.grch38.txt",
        cnvkit=CNVKIT_EXEC
    shell:
        """
        mkdir -p data/cnvkit
        {params.cnvkit} autobin {input.bams} -t {params.bed} -g {input.access} --annotate {params.ref_flat} --short-names --target-output-bed {output.targets} --antitarget-output-bed {output.antitargets}
        """

rule cnvkit_target_coverage:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample],
        bed=lambda wildcards: f"data/cnvkit/{config['project']['name']}.{config['resources']['targets_key']}.cnvkit-targets.bed"
    output:
        "data/work/{sample}/cnvkit/{sample}.targetcoverage.cnn"
    params:
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} coverage {input.bam} {input.bed} -o {output}
        """

rule cnvkit_antitarget_coverage:
    input:
        bam=lambda wildcards: BAMS[wildcards.sample],
        bed=f"data/cnvkit/{config['project']['name']}.{config['resources']['targets_key']}.cnvkit-antitargets.bed"
    output:
        'data/work/{sample}/cnvkit/{sample}.antitargetcoverage.cnn'
    params:
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} coverage {input.bam} {input.bed} -o {output}
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
        ref=config['reference']['fasta'],
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} reference {input} --fasta {params.ref} -o {output}
        """
#file names needs to be sample.{,anti}targetcoverage.cnn
rule cnvkit_fix:
    input:
        targets="data/work/{tumor}/cnvkit/{tumor}.targetcoverage.cnn",
        antitargets="data/work/{tumor}/cnvkit/{tumor}.antitargetcoverage.cnn",
        reference="data/cnvkit/reference.cnn"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.cnr"
    params:
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} fix {input.targets} {input.antitargets} {input.reference} -o {output}
        """

# Shared germline heterozygous SNPs for CNVkit segment/call and PureCN.
# Source: unfiltered VarDict germline (not filter2). Tumor=FORMAT[0], normal=FORMAT[1].
# Female cohort: autosomes + X; no Y.
rule germline_het_snps:
    input:
        germline="data/work/{tumor}/vardict/germline.vcf.gz"
    output:
        "data/work/{tumor}/snps/germline.hets.vcf.gz"
    resources:
        mem_mb=6144
    params:
        ref=config['reference']['fasta'],
        regions=config['resources']['targets_bedgz'],
        common=config['resources']['common_snps'],
        chroms="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
        tmp="data/work/{tumor}/snps/germline.hets.tmp.vcf.gz"
    shell:
        """
        mkdir -p data/work/{wildcards.tumor}/snps

        bcftools view -v snps -r {params.chroms} -R {params.regions} {input.germline} | \
        bcftools norm -m-both -f {params.ref} | \
        bcftools view -i 'FORMAT/DP[0]>=15 && FORMAT/DP[1]>=15 && FORMAT/AF[1]>=0.3 && FORMAT/AF[1]<=0.7' \
            -W=tbi -Oz -o {params.tmp}

        bcftools annotate -a {params.common} -m DB {params.tmp} | \
        bcftools sort -W=tbi -Oz -o {output}

        rm -f {params.tmp} {params.tmp}.tbi
        """

#can use threads with -p after you have a proper cluser config
rule cnvkit_segment:
    input:
        cnr="data/work/{tumor}/cnvkit/{tumor}.cnr",
        vcf="data/work/{tumor}/snps/germline.hets.vcf.gz"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.cns"
    params:
        tumor="{tumor}",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} segment {input.cnr} -v {input.vcf} -i {params.tumor} -n {params.normal} -o {output}
        """
#Native VarScan2 calls lack FORMAT:AF

#export to seg
rule cnvkit_export_seg:
    input:
        "data/work/{tumor}/cnvkit/{tumor}.cns"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.seg"
    params:
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} export seg {input} -o {output}
        """

rule purecn_run:
    input:
        snps="data/work/{tumor}/snps/germline.hets.vcf.gz",
        seg="data/work/{tumor}/cnvkit/{tumor}.seg",
        cnr="data/work/{tumor}/cnvkit/{tumor}.cnr"
    output:
        "data/work/{tumor}/purecn/{tumor}.csv",
        "data/work/{tumor}/purecn/{tumor}_loh.csv"
    params:
        seg="data/work/{tumor}/purecn/input.seg",
        cnr="data/work/{tumor}/purecn/input.cnr",
        outdir="data/work/{tumor}/purecn",
        sex=lambda wildcards: sample_sex(wildcards.tumor)['purecn']
    shell:
        """
        # Drop decoy/alt contigs (GL*, KI*, etc.); PureCN --genome grch38 wants primary chroms.
        awk 'NR==1 || $2 ~ /^([1-9]|1[0-9]|2[0-2]|[XY])$/' {input.seg} > {params.seg}
        awk 'NR==1 || $1 ~ /^([1-9]|1[0-9]|2[0-2]|[XY])$/' {input.cnr} > {params.cnr}

        Rscript /home/bwubb/software/PureCN/inst/extdata/PureCN.R \
        --out {params.outdir} \
        --sampleid {wildcards.tumor} \
        --tumor {params.cnr} \
        --seg-file {params.seg} \
        --vcf {input.snps} \
        --genome grch38 \
        --fun-segmentation Hclust \
        --sex {params.sex} \
        --force
        """

rule purecn_to_bed:
    input:
        "data/work/{tumor}/purecn/{tumor}_loh.csv"
    output:
        "data/work/{tumor}/purecn/{tumor}_loh.bed"
    shell:
        """
        python cnv_to_bed.py -c purecn {input}
        """

rule purecn_HRDex:
    input:
        "data/work/{tumor}/purecn/{tumor}_loh.bed"
    output:
        "data/work/{tumor}/purecn/{tumor}.purecn.hrd.txt"
    params:
        tumor="{tumor}",
        build="grch38"
    shell:
        """
        Rscript runHRDex.R -i {input} -o {output} --build {params.build} --tumor {params.tumor}
        """

rule purecn_annotsv:
    input:
        "data/work/{tumor}/purecn/{tumor}_loh.bed"
    output:
        "data/work/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule purecn_annotsv_parser:
    input:
        "data/work/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.tsv"
    output:
        "data/work/{tumor}/purecn/{tumor}.purecn.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

#--ploidy requires int
rule cnvkit_call:
    input:
        csv="data/work/{tumor}/purecn/{tumor}.csv",
        vcf="data/work/{tumor}/snps/germline.hets.vcf.gz",
        cns="data/work/{tumor}/cnvkit/{tumor}.cns"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.call.cns"
    params:
        tumor="{tumor}",
        normal=lambda wildcards: PAIRS[wildcards.tumor],
        sex=lambda wildcards: sample_sex(wildcards.tumor)['cnvkit'],
        cnvkit=CNVKIT_EXEC
    shell:
        """
        purity=`grep {params.tumor} {input.csv} | cut -d, -f2`
        ploidy=`grep {params.tumor} {input.csv} | cut -d, -f3 | cut -d. -f1`

        {params.cnvkit} call {input.cns} -x {params.sex} -m clonal --purity $purity --ploidy $ploidy -v {input.vcf} -i {params.tumor} -n {params.normal} -o {output}
        """

rule cnvkit_to_bed:
    input:
        "data/work/{tumor}/cnvkit/{tumor}.call.cns"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.call.bed"
    shell:
        """
        python cnv_to_bed.py -c cnvkit {input}
        """

rule cnvkit_HRDex:
    input:
        "data/work/{tumor}/cnvkit/{tumor}.call.bed"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.cnvkit.hrd.txt"
    params:
        tumor="{tumor}",
        build="grch38"
    shell:
        """
        Rscript runHRDex.R -i {input} -o {output} --build {params.build} --tumor {params.tumor}
        """

rule cnvkit_annotsv:
    input:
        "data/work/{tumor}/cnvkit/{tumor}.call.bed"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule cnvkit_annotsv_parser:
    input:
        "data/work/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.tsv"
    output:
        "data/work/{tumor}/cnvkit/{tumor}.cnvkit.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

rule purecn_purity_table:
    input:
        expand("data/work/{tumor}/purecn/{tumor}.csv",tumor=PAIRS.keys())
    output:
        "purecn_purity.table","purecn_ploidy.table"
    run:
        with open(output[0],'w') as outfile0, open(output[1],'w') as outfile1:
            for f in input:
                with open(f,'r') as infile:
                    line=infile.readline()#header
                    line=infile.readline()
                    sampleid,purity,ploidy=line.replace('"','').rstrip().split(',')[0:3]
                outfile0.write(f'{sampleid}\t{purity}\n')
                outfile1.write(f'{sampleid}\t{ploidy}\n')

#Germline exon del/dup: excluded normals vs clean pooled reference (not in reference.cnn).
#Purity/ploidy fixed at 1/2. BAF from matched tumor's unfiltered paired VarDict germline.
rule cnvkit_germline_fix:
    input:
        targets="data/work/{germline}/cnvkit/{germline}.targetcoverage.cnn",
        antitargets="data/work/{germline}/cnvkit/{germline}.antitargetcoverage.cnn",
        reference="data/cnvkit/reference.cnn"
    output:
        "data/work/{germline}/cnvkit/{germline}.germline.cnr"
    params:
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} fix {input.targets} {input.antitargets} {input.reference} -o {output}
        """

rule cnvkit_germline_input_vcf:
    input:
        unpack(germline_pair_vcf_inputs)
    output:
        "data/work/{germline}/cnvkit/vardict.snps.clean.vcf.gz"
    resources:
        mem_mb=6144
    params:
        ref=config['reference']['fasta'],
        regions=config['resources']['targets_bedgz'],
        chroms="1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X",
        tmp="data/work/{germline}/cnvkit/temp.reheader.vcf.gz"
    shell:
        """
        bcftools reheader -s {input.name} -o {params.tmp} {input.vcf}
        bcftools index -f {params.tmp}

        bcftools view -s {wildcards.germline} -v snps -r {params.chroms} -R {params.regions} {params.tmp} | \
        bcftools norm -m-both -f {params.ref} | \
        bcftools view -i 'FORMAT/DP>=15 && FORMAT/AF>=0.3 && FORMAT/AF<=0.7' | \
        bcftools sort -W=tbi -Oz -o {output}

        rm -f {params.tmp} {params.tmp}.tbi
        """

rule cnvkit_germline_segment:
    input:
        cnr="data/work/{germline}/cnvkit/{germline}.germline.cnr",
        vcf="data/work/{germline}/cnvkit/vardict.snps.clean.vcf.gz"
    output:
        "data/work/{germline}/cnvkit/{germline}.germline.cns"
    params:
        germline="{germline}",
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} segment {input.cnr} -v {input.vcf} -i {params.germline} -n {params.germline} -o {output}
        """

rule cnvkit_germline_call:
    input:
        vcf="data/work/{germline}/cnvkit/vardict.snps.clean.vcf.gz",
        cns="data/work/{germline}/cnvkit/{germline}.germline.cns"
    output:
        "data/work/{germline}/cnvkit/{germline}.germline.call.cns"
    params:
        germline="{germline}",
        sex=lambda wildcards: sample_sex(wildcards.germline)['cnvkit'],
        cnvkit=CNVKIT_EXEC
    shell:
        """
        {params.cnvkit} call {input.cns} -x {params.sex} -m clonal --purity 1 --ploidy 2 -v {input.vcf} -i {params.germline} -n {params.germline} -o {output}
        """

rule cnvkit_germline_to_bed:
    input:
        "data/work/{germline}/cnvkit/{germline}.germline.call.cns"
    output:
        "data/work/{germline}/cnvkit/{germline}.germline.call.bed"
    shell:
        """
        python cnv_to_bed.py -c cnvkit {input}
        """

rule cnvkit_germline_annotsv:
    input:
        "data/work/{germline}/cnvkit/{germline}.germline.call.bed"
    output:
        "data/work/{germline}/cnvkit/{germline}.germline.annotsv.gene_split.tsv"
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule cnvkit_germline_annotsv_parser:
    input:
        "data/work/{germline}/cnvkit/{germline}.germline.annotsv.gene_split.tsv"
    output:
        "data/work/{germline}/cnvkit/{germline}.germline.annotsv.gene_split.report.csv"
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.germline}
        """

### INIT ###

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())

### FUNCTIONS ###

def paired_bams(wildcards):
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

# Pileup shard file names always use chr1..chrX; SNP VCF contents match BAM contigs.
CHR=[f'chr{c}' for c in range(1,23)]+['chrX']

wildcard_constraints:
    work_dir=f"data/work/{config['resources']['targets_key']}",
    chr="chr[1-2][0-9]|chr[1-9]|chrX"

### SNAKEMAKE ###

rule collect_facets:
    input:
        expand("data/work/{lib}/{tumor}/facets/annotsv_gene_split.report.csv",lib=f"{config['resources']['targets_key']}",tumor=PAIRS.keys())

rule facets_chr_pileup:
    input:
        unpack(paired_bams)
    output:
        temp("{work_dir}/{tumor}/facets/pileup.{chr}.csv.gz")
    resources:
        mem_mb=6144
    params:
        snp=lambda wildcards: f"/home/bwubb/resources/Vcf_files/dbsnp156.GRCh38.snps.{wildcards.chr}.20240405.vcf.gz"
    shell:
        """
        snp-pileup -g -q15 -Q20 -P100 -r25,0 {params.snp} {output} {input.normal} {input.tumor}
        """

rule facets_merge_pileup:
    input:
        [f"{{work_dir}}/{{tumor}}/facets/pileup.{chr}.csv.gz" for chr in CHR]
    output:
        "{work_dir}/{tumor}/facets/pileup.csv.gz"
    resources:
        mem_mb=6144
    shell:
        """
        zcat {input} | awk 'NR == 1 || !/^Chromosome/' | bgzip -c > {output}
        """

# Lower cval -> higher sensitivity for small changes (procSample).
rule run_facets:
    input:
        "{work_dir}/{tumor}/facets/pileup.csv.gz"
    output:
        "{work_dir}/{tumor}/facets/segmentation_cncf.csv",
        "{work_dir}/{tumor}/facets/purity_ploidy.csv",
        "{work_dir}/{tumor}/facets/copynumber_profile.pdf",
        "{work_dir}/{tumor}/facets/fit_diagnostic.pdf",
        "{work_dir}/{tumor}/facets/{tumor}_segments.txt",
        "{work_dir}/{tumor}/facets/flags.txt"
    resources:
        mem_mb=18432
    params:
        cval=150,
        ndepth=25,
        gbuild="hg38"
    shell:
        """
        Rscript facets-snakemake.R --id {wildcards.tumor} --input {input} --cval {params.cval} --ndepth {params.ndepth} --gbuild {params.gbuild}
        """

rule facets_2bed:
    input:
        "{work_dir}/{tumor}/facets/segmentation_cncf.csv"
    output:
        "{work_dir}/{tumor}/facets/segmentation_cncf.bed"
    resources:
        mem_mb=6144
    shell:
        """
        python cnv_to_bed.py -c facets {input}
        """

rule facets_AnnotSV:
    input:
        "{work_dir}/{tumor}/facets/segmentation_cncf.bed"
    output:
        "{work_dir}/{tumor}/facets/annotsv.gene_split.tsv"
    resources:
        mem_mb=8192
    params:
        build=config['reference']['key']
    shell:
        """
        AnnotSV -SVinputFile {input} -annotationMode split -genomeBuild {params.build} -tx ENSEMBL -outputFile {output}
        """

rule facets_AnnotSV_parser:
    input:
        "{work_dir}/{tumor}/facets/annotsv.gene_split.tsv"
    output:
        "{work_dir}/{tumor}/facets/annotsv_gene_split.report.csv"
    resources:
        mem_mb=6144
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

#"Reference and variant allele read counts were extraced from the bams files for germline polymorphic
#sites cataloagues in the dvSNP and 1000genome databases."

with open(config.get('project',{}).get('sample_list','samples.list'),'r') as i:
    SAMPLES=i.read().splitlines()
    for sample in SAMPLES:
        os.makedirs(f'logs/cluster/{sample}',exist_ok=True)

with open(config.get('project',{}).get('pair_table','pair.table'),'r') as p:
    PAIRS=dict(line.split('\t') for line in p.read().splitlines())

with open(config['project']['bam_table'],'r') as b:
    BAMS=dict(line.split('\t') for line in b.read().splitlines())


def paired_bams(wildcards):
    ref=config['reference']['key']
    tumor=wildcards.tumor
    normal=PAIRS[wildcards.tumor]
    return {'tumor':BAMS[wildcards.tumor],'normal':BAMS[normal]}

#gnomad is not enough
#need dbsnp

#def gnomad_snp(wildcards):
#    #return f"/home/bwubb/resources/Vcf_files/gnomad.exomes.r2.0.2.sites.{config['resources']['targets_key']}.common_biallelic_snps.simplified.vcf.gz"
#    return f"/home/bwubb/resources/Vcf_files/dbsnp151.snps.20180423.vcf.gz"

wildcard_constraints:
    work_dir=f"data/work/{config['resources']['targets_key']}",
    chr="chr[1-2][0-9]|chr[1-9]|chrX"

rule collect_facets:
    input:
        expand("data/work/{lib}/{tumor}/facets/segmentation_cncf.csv",lib=f"{config['resources']['targets_key']}",tumor=PAIRS.keys())

#input quality can be adjusted
rule facets_pileup:
    input:
        unpack(paired_bams)
    output:
        "{work_dir}/{tumor}/facets/pileup.csv.gz"
    params:
        #snp=gnomad_snp
        #snp=config['resources']['dbsnp']
        snp="$HOME/resources/Vcf_files/dbsnp156.GCF_000001405.40.GRCh38.contig.sort.vcf.gz"
    shell:
        """
        snp-pileup -g -q15 -Q20 -P100 -r25,0 {params.snp} {output} {input.normal} {input.tumor}
        """

#Lower cval lead to higher sensitivity for small changes.
rule run_facets:
    input:
        "{work_dir}/{tumor}/facets/pileup.csv.gz"
    output:
        "{work_dir}/{tumor}/facets/segmentation_cncf.csv"
    params:
        cval=150
    shell:
        """
        Rscript facets-snakemake.R --id {wildcards.tumor} --input {input} --cval {params.cval}
        """

rule facets_2bed:
    input:
        "{work_dir}/{tumor}/facets/segmentation_cncf.csv"
    output:
        "{work_dir}/{tumor}/facets/segmentation_cncf.bed"
    shell:
        """
        python facets2bed.py {input}
        """

rule facets_AnnotSV:
    input:
        "{work_dir}/{tumor}/facets/segmentation_cncf.bed"
    output:
        "{work_dir}/{tumor}/facets/annotsv.gene_split.tsv"
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
    shell:
        """
        python annotsv_parser.py -i {input} -o {output} --tumor {wildcards.tumor}
        """

rule facets_hrd:
    input:
        segments="{work_dir}/{tumor}/facets/{tumor}_segments.txt"
    output:
        "{work_dir}/{tumor}/facets/{tumor}_hrd.txt"
    params:
        build=config['reference']['key'].lower()
    shell:
        """
        Rscript $HOME/software/HRDex/R/run_HRDex.R -i {input} -o {output} --tumor {wildcards.tumor} --build {params.build}
        """

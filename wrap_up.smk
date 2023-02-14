


#wrap up all the files I need.
#data/final/tumor/tumor.caller.unprocessed.vcf.gz
#...
#data/final/tumor/tumor.caller.filter_applied.vcf.gz
#clean/write vcf header fix for each. Contigs, signatures, filters
#
#vlr annotated data
#purity/ploidy table
#sequenza figures
#tmb figures
#run some sort of vcf stats

rule wrap_up:
    input:
        expand("data/final/{tumor}/{tumor}.lancet.normalized.vcf.gz",tumor=TUMORS),
        expand("data/final/{tumor}/{tumor}.mutect2.filter_applied.vcf.gz",tumor=TUMORS),
        expand("data/final/{tumor}/{tumor}.strelka2.normalized.vcf.gz",tumor=TUMORS),
        expand("data/final/{tumor}/{tumor}.vardict.filter_applied.vcf.gz",tumor=TUMORS),
        expand("data/final/{tumor}/{tumor}.varscan2.filter_applied.vcf.gz",tumor=TUMORS),
        expand("data/final/{tumor}/{tumor}.sequenza.segments.txt",tumor=TUMORS),
        expand()

rule copy_lancet:
    input:
        u="data/work/{config['resources']['targets_key']}/{tumor}/lancet/somatic.vcf.gz",
        f="data/work/{config['resources']['targets_key']}/{tumor}/lancet/somatic.filtered.norm.clean.vcf.gz"
    output:
        u="data/final/{tumor}/{tumor}.lancet.unprocessed.vcf.gz",
        f="data/final/{tumor}/{tumor}.lancet.normalized.vcf.gz"
    shell:
        """
        rsync -vr {input.u} {output.u}
        rsync -vr {input.f} {output.f}
        """

rule copy_mutect2:
    input:
        u="data/work/{config['resources']['targets_key']}/{tumor}/mutect2/somatic.vcf.gz",
        f="data/work/{config['resources']['targets_key']}/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    output:
        u="data/final/{tumor}/{tumor}.mutect2.unprocessed.vcf.gz",
        f="data/final/{tumor}/{tumor}.mutect2.filter_applied.vcf.gz"
    shell:
        """
        rsync -vr {input.u} {output.u}
        rsync -vr {input.f} {output.f}
        """

#paths need to be fixed for strelka2
rule copy_strelka2:
    input:
        u="data/work/{config['resources']['targets_key']}/{tumor}/strelka2/result/variants/somatic.vcf.gz",
        f="data/work/{config['resources']['targets_key']}/{tumor}/strelka2/results/variants/somatic.filtered.norm.clean.vcf.gz"
    output:
        u="data/final/{tumor}/{tumor}.strelka2.unprocessed.vcf.gz",
        f="data/final/{tumor}/{tumor}.strelka2.normalized.vcf.gz"
    shell:
        """
        rsync -vr {input.u} {output.u}
        rsync -vr {input.f} {output.f}
        """

rule copy_vardict:
    input:
        u="data/work/{config['resources']['targets_key']}/{tumor}/vardict/somatic.vcf.gz",
        f="data/work/{config['resources']['targets_key']}/{tumor}/vardict/somatic.twice_filtered.norm.clean.vcf.gz"
    output:
        u="data/final/{tumor}/{tumor}.vardict.unprocessed.vcf.gz",
        f="data/final/{tumor}/{tumor}.vardict.filter_applied.vcf.gz"
    shell:
        """
        rsync -vr {input.u} {output.u}
        rsync -vr {input.f} {output.f}
        """

rule copy_varscan2:
    input:
        u="data/work/{config['resources']['targets_key']}/{tumor}/mutect2/somatic.vcf.gz",
        f="data/work/{config['resources']['targets_key']}/{tumor}/mutect2/somatic.filtered.norm.clean.vcf.gz"
    output:
        u="data/final/{tumor}/{tumor}.varscan2.unprocessed.vcf.gz",
        f="data/final/{tumor}/{tumor}.varscan2.filter_applied.vcf.gz"
    shell:
        """
        rsync -vr {input.u} {output.u}
        rsync -vr {input.f} {output.f}
        """

#lacks flexibility if its not an ffpe scenario
rule copy_varlociraptor:
    input:
        u="data/work/{config['resources']['targets_key']}/{tumor}/varlociraptor/ffpe_scenario.bcf",
        f="data/work/{config['resources']['targets_key']}/{tumor}/varlociraptor/ffpe_scenario.local-fdr.vcf"
    output:
        u="data/final/{tumor}/{tumor}.varlociraptor.unprocessed.vcf.gz",
        f="data/final/{tumor}/{tumor}.varlociraptor.filter_applied.vcf.gz",
    shell:
        """
        bcftools view -Oz -o {output.u} {input.u}
        bgzip -c {input.f} > {output.f}
        """

rule copy_vep:
    input:
    output:
    shell:
        """
        """

rule copy_sequenza:
    input:
        "data/work/{config['resources']['targets_key']}/{tumor}/sequenza/{tumor}_segments.txt",
        "data/work/{config['resources']['targets_key']}/{tumor}/sequenza/{tumor}_confints_CP.txt",
        "data/work/{config['resources']['targets_key']}/{tumor}/sequenza/{tumor}_model_fit.pdf"
    output:
        "data/final/{tumor}/{tumor}.sequenza.segments.txt",
        "data/final/{tumor}/{tumor}.sequenza.confints_CP.txt",
        "data/final/{tumor}/{tumor}.sequenza.model_fit.pdf"
    shell:
        """
        rsync -vr {input[0]} {output[0]}
        rsync -vr {input[1]} {output[1]}
        rsync -vr {input[2]} {output[2]}
        """

rule copy_facets:
    input: pass

rule copy_cnvkit:
    input: pass

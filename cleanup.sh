#

##copy and clean
##I like my sugar with coffee and cream.

#find data/final/$lib/$name/ -iname *.vcf.gz and exec tabix -f -p vcf

#RENAME all somatic.* to callername
#for F in data/final/$lib/$name/manta/somatic.*; do ${F} ; done

name=$1
lib=$2

#manta
mkdir -p data/final/$lib/$name/manta
rsync -v data/work/$lib/$name/manta/results/variants/candidateSmallIndels.vcf.gz data/final/$lib/$name/manta/$name.candidateSmallIndels.vcf.gz
rsync -v data/work/$lib/$name/manta/results/variants/candidateSV.vcf.gz data/final/$lib/$name/manta/$name.candidateSV.vcf.gz
rsync -v data/work/$lib/$name/manta/results/variants/diploidSV.vcf.gz data/final/$lib/$name/manta/$name.diploidSV.vcf.gz
rsync -v data/work/$lib/$name/manta/results/variants/somaticSV.vcf.gz data/final/$lib/$name/manta/$name.somaticSV.vcf.gz
rsync -vr data/work/$lib/$name/manta/results/stats/* data/final/$lib/$name/manta/

#strelka2
mkdir -p data/final/$lib/$name/strelka2
rsync -vr data/work/$lib/$name/strelka2/results/variants/somatic.raw.vcf.gz data/final/$lib/$name/strelka2/$name.somatic.raw.vcf.gz
rsync -vr data/work/$lib/$name/strelka2/results/variants/somatic.raw.norm.clean.std.vcf.gz data/final/$lib/$name/strelka2/$name.somatic.raw.norm.clean.std.vcf.gz
#strelka2/stats, not much there.

#mutect2
mkdir -p data/final/$lib/$name/mutect2
rsync -v data/work/$lib/$name/mutect2/somatic.raw.vcf.gz data/final/$lib/$name/mutect2/$name.somatic.raw.vcf.gz
rsync -v data/work/$lib/$name/mutect2/somatic.filtered.vcf.gz data/final/$lib/$name/mutect2/$name.somatic.filtered.vcf.gz
rsync -v data/work/$lib/$name/mutect2/somatic.filtered.norm.clean.std.vcf.gz data/final/$lib/$name/mutect2/$name.somatic.filtered.norm.clean.std.vcf.gz
rsync -vr data/work/$lib/$name/mutect2/*_metrics.txt data/final/$lib/$name/mutect2/

#vardict
mkdir -p data/final/$lib/$name/vardict
#rsync -vr data/work/$lib/$name/vardict/variants.raw.vcf.gz data/final/$lib/$name/vardict/variants.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/raw.vcf.gz data/final/$lib/$name/vardict/$name.variants.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/somatic.raw.vcf.gz data/final/$lib/$name/vardict/$name.somatic.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/loh.raw.vcf.gz data/final/$lib/$name/vardict/$name.loh.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/germline.raw.vcf.gz data/final/$lib/$name/vardict/$name.germline.raw.vcf.gz

rsync -vr data/work/$lib/$name/vardict/variants.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.variants.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/somatic.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.somatic.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/germline.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.germline.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/loh.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.loh.twice_filtered.vcf.gz


rsync -vr data/work/$lib/$name/vardict/somatic.twice_filtered.norm.clean.std.vcf.gz data/final/$lib/$name/vardict/$name.somatic.twice_filtered.norm.clean.std.vcf.gz
rsync -vr data/work/$lib/$name/vardict/germline.twice_filtered.norm.clean.std.vcf.gz data/final/$lib/$name/vardict/$name.germline.twice_filtered.norm.clean.std.vcf.gz
rsync -vr data/work/$lib/$name/vardict/loh.twice_filtered.norm.clean.std.vcf.gz data/final/$lib/$name/vardict/$name.loh.twice_filtered.norm.clean.std.vcf.gz

#varscan2
mkdir -p data/final/$lib/$name/varscan2

rsync -vr data/work/$lib/$name/varscan2/raw.Germline.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.snp.Germline.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.indel.Germline.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/raw.LOH.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.snp.LOH.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.indel.LOH.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.indel.vcf.gz


#rsync -vr data/work/$lib/$name/varscan2/raw.snp.Somatic.vcf.gz data/work/$lib/$name/varscan2/$name.somatic.raw.snp.vcf.gz
#rsync -vr data/work/$lib/$name/varscan2/raw.indel.Somatic.vcf.gz data/work/$lib/$name/varscan2/$name.somatic.raw.indel.vcf.gz
if [ ! -e data/work/$lib/$name/varscan2/somatic.raw.vcf.gz ]; then bcftools concat -a data/work/$lib/$name/varscan2/raw.snp.Somatic.vcf.gz data/work/$lib/$name/varscan2/raw.indel.Somatic.vcf.gz | bcftools sort -O z -o data/work/$lib/$name/varscan2/somatic.raw.vcf.gz; fi

rsync -vr data/work/$lib/$name/varscan2/somatic.raw.vcf.gz data/final/$lib/$name/varscan2/$name.somatic.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.Germline.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.LOH.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/somatic.raw.snp.vcf.gz data/final/$lib/$name/varscan2/$name.somatic.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/somatic.raw.indel.vcf.gz data/final/$lib/$name/varscan2/$name.somatic.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/germline.raw.snp.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/germline.raw.indel.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/loh.raw.snp.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/loh.raw.indel.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/somatic.fpfilter.vcf.gz data/final/$lib/$name/varscan2/$name.somatic.fpfilter.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/somatic.fpfilter.norm.clean.std.vcf.gz data/final/$lib/$name/varscan2/$name.somatic.fpfilter.norm.clean.std.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/germline.fpfilter.vcf.gz data/final/$lib/$name/varscan2/$name.germline.fpfilter.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/germline.fpfilter.norm.clean.std.vcf.gz data/final/$lib/$name/varscan2/$name.germline.fpfilter.norm.clean.std.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/loh.fpfilter.vcf.gz data/final/$lib/$name/varscan2/$name.loh.fpfilter.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/loh.fpfilter.norm.clean.std.vcf.gz data/final/$lib/$name/varscan2/$name.loh.fpfilter.norm.clean.std.vcf.gz

#sequenza
mkdir -p data/final/$lib/$name/sequenza
rsync -vr data/work/$lib/$name/sequenza/*.pdf data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*_confints_CP.txt data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*_mutations.txt data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*_segments.txt data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*.RData data/final/$lib/$name/sequenza/

#bcftools
mkdir -p data/final/$lib/$name/bcftools
bgzip -c data/work/$lib/$name/bcftools/somtic/0000.vcf > data/final/$lib/$name/bcftools/$name.somatic.0000.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/somatic/0001.vcf > data/final/$lib/$name/bcftools/$name.somatic.0001.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/somatic/0002.vcf > data/final/$lib/$name/bcftools/$name.somatic.0002.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/somatic/0003.vcf > data/final/$lib/$name/bcftools/$name.somatic.0003.vcf.gz
rsync -vr data/work/$lib/$name/bcftools/somatic/sites.txt data/final/$lib/$name/bcftools/$name.somatic.sites.txt



#What files are missing?

#missing vardict/variants.twice_filtered.vcf.gz
#missing loh and germline

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
rsync -vr data/work/$lib/$name/manta/results/variants/*.vcf.gz data/final/$lib/$name/manta/
#manta/stats

#strelka2
mkdir -p data/final/$lib/$name/strelka2
rsync -vr data/work/$lib/$name/strelka2/results/variants/*.vcf.gz data/final/$lib/$name/strelka2/
rsync -vr data/work/$lib/$name/strelka2/somatic.raw.vcf.gz data/final/$lib/$name/strelka2/somatic.raw.vcf.gz #temp
#strelka2/stats, not much there.


#mutect2
mkdir -p data/final/$lib/$name/mutect2
rsync -vr data/work/$lib/$name/mutect2/somatic.raw.vcf.gz data/final/$lib/$name/mutect2/somatic.raw.vcf.gz
rsync -vr data/work/$lib/$name/mutect2/somatic.twice_filtered.vcf.gz data/final/$lib/$name/mutect2/somatic.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/mutect2/*_metrics.txt data/final/$lib/$name/mutect2/

#vardict
mkdir -p data/final/$lib/$name/vardict
#rsync -vr data/work/$lib/$name/vardict/variants.raw.vcf.gz data/final/$lib/$name/vardict/variants.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/raw.vcf.gz data/final/$lib/$name/vardict/variants.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/somatic.raw.vcf.gz data/final/$lib/$name/vardict/somatic.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/twice_filtered.vcf.gz data/final/$lib/$name/vardict/variants.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/somatic.twice_filtered.vcf.gz data/final/$lib/$name/vardict/somatic.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/germline.twice_filtered.vcf.gz data/final/$lib/$name/vardict/germline.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/loh.twice_filtered.vcf.gz data/final/$lib/$name/vardict/loh.twice_filtered.vcf.gz

#varscan2
mkdir -p data/final/$lib/$name/varscan2

cp data/work/$lib/$name/varscan2/raw.Germline.vcf.gz data/work/$lib/$name/varscan2/germline.raw.vcf.gz
cp data/work/$lib/$name/varscan2/raw.snp.Germline.vcf.gz data/work/$lib/$name/varscan2/germline.raw.snp.vcf.gz
cp data/work/$lib/$name/varscan2/raw.indel.Germline.vcf.gz data/work/$lib/$name/varscan2/germline.raw.indel.vcf.gz

cp data/work/$lib/$name/varscan2/raw.LOH.vcf.gz data/work/$lib/$name/varscan2/loh.raw.vcf.gz
cp data/work/$lib/$name/varscan2/raw.snp.LOH.vcf.gz data/work/$lib/$name/varscan2/loh.raw.snp.vcf.gz
cp data/work/$lib/$name/varscan2/raw.indel.LOH.vcf.gz data/work/$lib/$name/varscan2/loh.raw.indel.vcf.gz

if [ ! -e data/work/$lib/$name/varscan2/somatic.raw.snp.vcf.gz ]; then cp data/work/$lib/$name/varscan2/raw.Somatic.vcf.gz data/work/$lib/$name/varscan2/somatic.raw.vcf.gz; fi
cp data/work/$lib/$name/varscan2/raw.snp.Somatic.vcf.gz data/work/$lib/$name/varscan2/somatic.raw.snp.vcf.gz
cp data/work/$lib/$name/varscan2/raw.indel.Somatic.vcf.gz data/work/$lib/$name/varscan2/somatic.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/somatic.raw.vcf.gz data/final/$lib/$name/varscan2/somatic.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.Germline.vcf.gz data/final/$lib/$name/varscan2/germline.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.LOH.vcf.gz data/final/$lib/$name/varscan2/loh.raw.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/somatic.raw.snp.vcf.gz data/final/$lib/$name/varscan2/somatic.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/somatic.raw.indel.vcf.gz data/final/$lib/$name/varscan2/somatic.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/germline.raw.snp.vcf.gz data/final/$lib/$name/varscan2/germline.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/germline.raw.indel.vcf.gz data/final/$lib/$name/varscan2/germline.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/loh.raw.snp.vcf.gz data/final/$lib/$name/varscan2/loh.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/loh.raw.indel.vcf.gz data/final/$lib/$name/varscan2/loh.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/somatic.fpfilter.vcf.gz data/final/$lib/$name/varscan2/somatic.fpfilter.vcf.gz

#NEED
#loh.fpfilter.vcf.gz
#germline.fpfilter.vcf.gz

#sequenza
mkdir -p data/final/$lib/$name/sequenza
rsync -vr data/work/$lib/$name/sequenza/*.pdf data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*_confints_CP.txt data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*_mutations.txt data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*_segments.txt data/final/$lib/$name/sequenza/
rsync -vr data/work/$lib/$name/sequenza/*.RData data/final/$lib/$name/sequenza/

#bcftools
mkdir -p data/final/$lib/$name/bcftools
bgzip -c data/work/$lib/$name/bcftools/0000.vcf > data/final/$lib/$name/bcftools/somatic.0000.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/0001.vcf > data/final/$lib/$name/bcftools/somatic.0001.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/0002.vcf > data/final/$lib/$name/bcftools/somatic.0002.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/0003.vcf > data/final/$lib/$name/bcftools/somatic.0003.vcf.gz
rsync -vr data/work/$lib/$name/bcftools/sites.txt data/final/$lib/$name/bcftools/somatic.sites.txt



#What files are missing?


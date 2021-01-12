
name=$1
lib=$2

#bcftools
mkdir -p data/final/$lib/$name/bcftools
bgzip -c data/work/$lib/$name/bcftools/somatic/0000.vcf > data/final/$lib/$name/bcftools/$name.somatic.0000.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/somatic/0001.vcf > data/final/$lib/$name/bcftools/$name.somatic.0001.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/somatic/0002.vcf > data/final/$lib/$name/bcftools/$name.somatic.0002.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/somatic/0003.vcf > data/final/$lib/$name/bcftools/$name.somatic.0003.vcf.gz
rsync -vr data/work/$lib/$name/bcftools/somatic/sites.txt data/final/$lib/$name/bcftools/$name.somatic.sites.txt
rsync -vr data/work/$lib/$name/bcftools/somatic/somatic.ensemble.vcf.gz data/final/$lib/$name/bcftools/$name.somatic.ensemble.vcf.gz



bgzip -c data/work/$lib/$name/bcftools/germline-loh/0000.vcf > data/final/$lib/$name/bcftools/$name.germline-loh.0000.vcf.gz
bgzip -c data/work/$lib/$name/bcftools/germline-loh/0001.vcf > data/final/$lib/$name/bcftools/$name.germline-loh.0001.vcf.gz
rsync -vr data/work/$lib/$name/bcftools/germline-loh/sites.txt data/final/$lib/$name/bcftools/$name.germline-loh.sites.txt
rsync -vr data/work/$lib/$name/bcftools/germline-loh/germline-loh.ensemble.vcf.gz data/final/$lib/$name/bcftools/$name.germline-loh.ensemble.vcf.gz





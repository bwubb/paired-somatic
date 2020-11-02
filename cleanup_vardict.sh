
name=$1
lib=$2

#vardict
mkdir -p data/final/$lib/$name/vardict
#rsync -vr data/work/$lib/$name/vardict/variants.raw.vcf.gz data/final/$lib/$name/vardict/variants.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/raw.vcf.gz data/final/$lib/$name/vardict/$name.variants.raw.vcf.gz
rsync -vr data/work/$lib/$name/vardict/somatic.raw.vcf.gz data/final/$lib/$name/vardict/$name.somatic.raw.vcf.gz

rsync -vr data/work/$lib/$name/vardict/variants.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.variants.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/somatic.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.somatic.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/germline.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.germline.twice_filtered.vcf.gz
rsync -vr data/work/$lib/$name/vardict/loh.twice_filtered.vcf.gz data/final/$lib/$name/vardict/$name.loh.twice_filtered.vcf.gz

rsync -vr data/work/$lib/$name/vardict/somatic.twice_filtered.norm.clean.std.vcf.gz data/final/$lib/$name/vardict/$name.somatic.twice_filtered.norm.clean.std.vcf.gz
rsync -vr data/work/$lib/$name/vardict/germline.twice_filtered.norm.clean.std.vcf.gz data/final/$lib/$name/vardict/$name.germline.twice_filtered.norm.clean.std.vcf.gz
rsync -vr data/work/$lib/$name/vardict/loh.twice_filtered.norm.clean.std.vcf.gz data/final/$lib/$name/vardict/$name.loh.twice_filtered.norm.clean.std.vcf.gz

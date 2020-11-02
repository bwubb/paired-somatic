
name=$1
lib=$2

#mutect2
mkdir -p data/final/$lib/$name/mutect2
rsync -v data/work/$lib/$name/mutect2/somatic.raw.vcf.gz data/final/$lib/$name/mutect2/$name.somatic.raw.vcf.gz
rsync -v data/work/$lib/$name/mutect2/somatic.filtered.vcf.gz data/final/$lib/$name/mutect2/$name.somatic.filtered.vcf.gz
rsync -v data/work/$lib/$name/mutect2/somatic.filtered.norm.clean.std.vcf.gz data/final/$lib/$name/mutect2/$name.somatic.filtered.norm.clean.std.vcf.gz
rsync -vr data/work/$lib/$name/mutect2/*_metrics.txt data/final/$lib/$name/mutect2/
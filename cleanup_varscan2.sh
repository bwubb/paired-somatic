

#varscan2
mkdir -p data/final/$lib/$name/varscan2

rsync -vr data/work/$lib/$name/varscan2/raw.Germline.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.snp.Germline.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.indel.Germline.vcf.gz data/final/$lib/$name/varscan2/$name.germline.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/raw.LOH.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.snp.LOH.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.indel.LOH.vcf.gz data/final/$lib/$name/varscan2/$name.loh.raw.indel.vcf.gz

rsync -vr data/work/$lib/$name/varscan2/raw.snp.Somatic.vcf.gz data/work/$lib/$name/varscan2/$name.somatic.raw.snp.vcf.gz
rsync -vr data/work/$lib/$name/varscan2/raw.indel.Somatic.vcf.gz data/work/$lib/$name/varscan2/$name.somatic.raw.indel.vcf.gz
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

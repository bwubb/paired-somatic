#mkdir -p data/work/S04380110/"$1"/bcftools

ref="/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta"

#bcftools isec -n +1 -p data/work/S04380110/"$1"/bcftools data/final/BasserExome_S04380110/"$1".S04380110.mutect2.somatic.twice_filtered.pass.vcf.gz data/final/BasserExome_S04380110/"$1".S04380110.strelka2.somatic.pass.vcf.gz data/final/BasserExome_S04380110/"$1".S04380110.vardict.somatic.twice_filtered.pass.vcf.gz data/final/BasserExome_S04380110/"$1".S04380110.varscan2.somatic.fpfilter.pass.vcf.gz
#bcftools isec -n +2 -p data/work/S04380110/"$1"/bcftools \

#bcftools norm -m-both data/work/S04380110/"$1"/mutect2/somatic.twice_filtered.vcf.gz | bcftools norm -f "$ref" -O z -o data/work/S04380110/"$1"/mutect2/somatic.twice_filtered.norm.vcf.gz
#tabix -p vcf data/work/S04380110/"$1"/mutect2/somatic.twice_filtered.norm.vcf.gz

#bcftools norm -m-both data/work/S04380110/"$1"/strelka2/somatic.raw.vcf.gz | bcftools norm -f "$ref" -O z -o data/work/S04380110/"$1"/strelka2/somatic.raw.norm.vcf.gz
#tabix -p vcf data/work/S04380110/"$1"/strelka2/somatic.raw.norm.vcf.gz

#bcftools norm -m-both data/work/S04380110/"$1"/vardict/somatic.twice_filtered.vcf.gz | bcftools norm -f "$ref" -O z -o data/work/S04380110/"$1"/vardict/somatic.twice_filtered.norm.vcf.gz
#tabix -p vcf data/work/S04380110/"$1"/vardict/somatic.twice_filtered.norm.vcf.gz

#bcftools norm -m-both data/work/S04380110/"$1"/varscan2/somatic.fpfilter.vcf.gz | bcftools norm -f "$ref" -O z -o data/work/S04380110/"$1"/varscan2/somatic.fpfilter.norm.vcf.gz
#tabix -p vcf data/work/S04380110/"$1"/varscan2/somatic.fpfilter.norm.vcf.gz

tabix -p vcf data/work/S04380110/"$1"/mutect2/somatic.twice_filtered.norm.std.vcf.gz
tabix -p vcf data/work/S04380110/"$1"/strelka2/somatic.raw.norm.std.vcf.gz
tabix -p vcf data/work/S04380110/"$1"/vardict/somatic.twice_filtered.norm.std.vcf.gz
tabix -p vcf data/work/S04380110/"$1"/varscan2/somatic.fpfilter.norm.vcf.gz

bcftools isec -n +2 -p data/work/S04380110/"$1"/bcftools \
data/work/S04380110/"$1"/mutect2/somatic.twice_filtered.norm.std.vcf.gz \
data/work/S04380110/"$1"/strelka2/somatic.raw.norm.std.vcf.gz \
data/work/S04380110/"$1"/vardict/somatic.twice_filtered.norm.std.vcf.gz \
data/work/S04380110/"$1"/varscan2/somatic.fpfilter.norm.vcf.gz

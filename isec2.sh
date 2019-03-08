


echo data/work/"$1"/S07604715/mutect/somatic.twice_filtered.vcf.gz
tabix -f -p vcf data/work/"$1"/S07604715/mutect/somatic.twice_filtered.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/mutect/somatic.twice_filtered.vcf.gz | bcftools norm -m-both | bcftools view -e 'ALT~"*"' -O z -o data/work/"$1"/S07604715/mutect/somatic.twice_filtered.norm.vcf.gz && tabix -f -p vcf data/work/"$1"/S07604715/mutect/somatic.twice_filtered.norm.vcf.gz

echo data/work/"$1"/S07604715/strelka/somatic.raw.vcf.gz
tabix -f -p vcf data/work/"$1"/S07604715/strelka/somatic.raw.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/strelka/somatic.raw.vcf.gz | bcftools norm -m-both | bcftools view -e 'ALT~"*"' -O z -o data/work/"$1"/S07604715/strelka/somatic.raw.norm.vcf.gz && tabix -f -p vcf data/work/"$1"/S07604715/strelka/somatic.raw.norm.vcf.gz

echo data/work/"$1"/S07604715/vardict/somatic.twice_filtered.vcf.gz
tabix -f -p vcf data/work/"$1"/S07604715/vardict/somatic.twice_filtered.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/vardict/somatic.twice_filtered.vcf.gz | bcftools norm -m-both | bcftools view -e 'ALT~"*"' -O z -o data/work/"$1"/S07604715/vardict/somatic.twice_filtered.norm.vcf.gz  && tabix -f -p vcf data/work/"$1"/S07604715/vardict/somatic.twice_filtered.norm.vcf.gz

echo data/work/"$1"/S07604715/varscan/somatic.fpfilter.vcf.gz
tabix -f -p vcf data/work/"$1"/S07604715/varscan/somatic.fpfilter.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/varscan/somatic.fpfilter.vcf.gz | bcftools norm -m-both | bcftools view -e 'ALT~"*"' -O z -o data/work/"$1"/S07604715/varscan/somatic.fpfilter.norm.vcf.gz && tabix -f -p vcf data/work/"$1"/S07604715/varscan/somatic.fpfilter.norm.vcf.gz


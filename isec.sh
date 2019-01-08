


#for each vcf file
tabix -f -p vcf data/work/"$1"/S07604715/mutect/somatic.twice_filtered.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/mutect/somatic.twice_filtered.vcf.gz | bcftools norm -m-both -O z -o data/work/"$1"/S07604715/mutect/somatic.twice_filtered.norm.vcf.gz && tabix -f -p vcf data/work/"$1"/S07604715/mutect/somatic.twice_filtered.norm.vcf.gz

tabix -f -p vcf data/work/"$1"/S07604715/strelka/somatic.raw.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/strelka/somatic.raw.vcf.gz | bcftools norm -m-both -O z -o data/work/"$1"/S07604715/strelka/somatic.raw.norm.vcf.gz && tabix -f -p vcf data/work/"$1"/S07604715/strelka/somatic.raw.norm.vcf.gz

tabix -f -p vcf data/work/"$1"/S07604715/vardict/somatic.twice_filtered.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/vardict/somatic.twice_filtered.vcf.gz | bcftools norm -m-both -O z -o data/work/"$1"/S07604715/vardict/somatic.twice_filtered.norm.vcf.gz  && tabix -f -p vcf data/work/"$1"/S07604715/vardict/somatic.twice_filtered.norm.vcf.gz

tabix -f -p vcf data/work/"$1"/S07604715/varscan/somatic.fpfilter.vcf.gz
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/"$1"/S07604715/varscan/somatic.fpfilter.vcf.gz | bcftools norm -m-both -O z -o data/work/"$1"/S07604715/varscan/somatic.fpfilter.norm.vcf.gz && tabix -f -p vcf data/work/"$1"/S07604715/varscan/somatic.fpfilter.norm.vcf.gz


bcftools isec -n +2 -f PASS -p data/work/"$1"/S07604715/concordant_calls/ data/work/"$1"/S07604715/mutect/somatic.twice_filtered.norm.vcf.gz data/work/"$1"/S07604715/strelka/somatic.raw.norm.vcf.gz data/work/"$1"/S07604715/vardict/somatic.twice_filtered.norm.vcf.gz data/work/"$1"/S07604715/varscan/somatic.fpfilter.norm.vcf.gz

bcftools isec -n =1 -f PASS -p data/work/"$1"/S07604715/private_calls/ data/work/"$1"/S07604715/mutect/somatic.twice_filtered.norm.vcf.gz data/work/"$1"/S07604715/strelka/somatic.raw.norm.vcf.gz data/work/"$1"/S07604715/vardict/somatic.twice_filtered.norm.vcf.gz data/work/"$1"/S07604715/varscan/somatic.fpfilter.norm.vcf.gz

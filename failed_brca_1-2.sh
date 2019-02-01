


mkdir -p data/work/$1/S07604715/failed_calls

echo -e "$2\n$1\n" > data/work/$1/S07604715/names.txt

bcftools view -e 'FILTER=="PASS"' -r 13:32889611-32973805,17:41196312-41277500 -O z -o data/work/$1/S07604715/failed_calls/mutect.brca_1-2.failed.vcf.gz data/work/$1/S07604715/mutect/somatic.twice_filtered.vcf.gz
tabix -p vcf data/work/$1/S07604715/failed_calls/mutect.brca_1-2.failed.vcf.gz

#rename reorder
bcftools reheader -s data/work/$1/S07604715/names.txt -o data/work/$1/S07604715/strelka/rename.vcf.gz data/work/$1/S07604715/strelka/somatic.raw.vcf.gz
tabix -p vcf data/work/$1/S07604715/strelka/rename.vcf.gz

bcftools view -O z -o data/work/$1/S07604715/strelka/rename.reorder.vcf.gz -s "$1","$2" data/work/$1/S07604715/strelka/rename.vcf.gz
tabix -p vcf data/work/$1/S07604715/strelka/rename.reorder.vcf.gz

python strelka_fix.py $1
bgzip data/work/$1/S07604715/strelka/rename.reorder.fix.vcf
tabix -p vcf data/work/$1/S07604715/strelka/rename.reorder.fix.vcf.gz

bcftools view -e 'FILTER=="PASS"' -r 13:32889611-32973805,17:41196312-41277500 -O z -o data/work/$1/S07604715/failed_calls/strelka.brca_1-2.failed.vcf.gz data/work/$1/S07604715/strelka/rename.reorder.fix.vcf.gz
tabix -p vcf data/work/$1/S07604715/failed_calls/strelka.brca_1-2.failed.vcf.gz

bcftools view -e 'FILTER=="PASS"' -r 13:32889611-32973805,17:41196312-41277500 -O z -o data/work/$1/S07604715/failed_calls/vardict.brca_1-2.failed.vcf.gz data/work/$1/S07604715/vardict/somatic.twice_filtered.vcf.gz
tabix -p vcf data/work/$1/S07604715/failed_calls/vardict.brca_1-2.failed.vcf.gz


#rename reorder
bcftools reheader -s data/work/$1/S07604715/names.txt -o data/work/$1/S07604715/varscan/rename.vcf.gz data/work/$1/S07604715/varscan/somatic.fpfilter.vcf.gz
tabix -p vcf data/work/$1/S07604715/varscan/rename.vcf.gz

bcftools view -O z -o data/work/$1/S07604715/varscan/rename.reorder.vcf.gz -s "$1","$2" data/work/$1/S07604715/varscan/rename.vcf.gz
tabix -p vcf data/work/$1/S07604715/varscan/rename.reorder.vcf.gz

bcftools view -e 'FILTER=="PASS"' -r 13:32889611-32973805,17:41196312-41277500 -O z -o data/work/$1/S07604715/failed_calls/varscan.brca_1-2.failed.vcf.gz data/work/$1/S07604715/varscan/rename.reorder.vcf.gz
tabix -p vcf data/work/$1/S07604715/failed_calls/varscan.brca_1-2.failed.vcf.gz


####
bcftools concat --allow-overlaps -O z -o data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.vcf.gz data/work/$1/S07604715/failed_calls/mutect.brca_1-2.failed.vcf.gz data/work/$1/S07604715/failed_calls/vardict.brca_1-2.failed.vcf.gz data/work/$1/S07604715/failed_calls/varscan.brca_1-2.failed.vcf.gz
tabix -p vcf data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.vcf.gz

bcftools sort -m 10G -O z -o data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.vcf.gz data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.vcf.gz
tabix -p vcf data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.vcf.gz

bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.vcf.gz | bcftools norm -m-both -O z -o data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.norm.vcf.gz
tabix -f -p vcf data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.norm.vcf.gz

table_annovar.pl --buildver hg19 --vcfinput --outfile data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.norm --protocol refGene,cytoband,gwasCatalog,genomicSuperDups,dbscsnv11,dbnsfp33a,popfreq_max_20150413,exac03,exac03nontcga,gnomad_exome,avsnp150,cosmic84_coding,cosmic84_noncoding,clinvar_20180603 --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f -remove data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.norm.vcf.gz $HOME/resources/annovar/humandb/

python annotated_somatic_vcf_v4.py --no_filter -I data/work/$1/S07604715/failed_calls/concat.brca_1-2.failed.sort.norm.hg19_multianno.vcf --tumor_id $1 --normal_id $2
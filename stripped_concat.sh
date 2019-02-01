

#strelka2
bcftools reheader -s data/work/$1/S07604715/names.txt -o data/work/$1/S07604715/concordant_calls/0001.rename.vcf data/work/$1/S07604715/concordant_calls/0001.vcf
bcftools view -O v -o data/work/$1/S07604715/concordant_calls/0001.rename.reorder.vcf -s "$1","$2" data/work/$1/S07604715/concordant_calls/0001.rename.vcf
python strelka_fix2.py $1


#varscan
bcftools reheader -s data/work/$1/S07604715/names.txt -o data/work/$1/S07604715/concordant_calls/0003.rename.vcf data/work/$1/S07604715/concordant_calls/0003.vcf
bcftools view -O v -o data/work/$1/S07604715/concordant_calls/0003.rename.reorder.vcf -s "$1","$2" data/work/$1/S07604715/concordant_calls/0003.rename.vcf

bgzip -c data/work/$1/S07604715/concordant_calls/0000.vcf > data/work/$1/S07604715/concordant_calls/0000.vcf.gz
tabix -p vcf data/work/$1/S07604715/concordant_calls/0000.vcf.gz

bgzip -c data/work/$1/S07604715/concordant_calls/0001.rename.reorder.fix.vcf > data/work/$1/S07604715/concordant_calls/0001.rename.reorder.fix.vcf.gz
tabix -p vcf data/work/$1/S07604715/concordant_calls/0001.rename.reorder.fix.vcf.gz

bgzip -c data/work/$1/S07604715/concordant_calls/0002.vcf > data/work/$1/S07604715/concordant_calls/0002.vcf.gz
tabix -p vcf data/work/$1/S07604715/concordant_calls/0002.vcf.gz

bgzip -c data/work/$1/S07604715/concordant_calls/0003.rename.reorder.vcf > data/work/$1/S07604715/concordant_calls/0003.rename.reorder.vcf.gz
tabix -p vcf data/work/$1/S07604715/concordant_calls/0003.rename.reorder.vcf.gz

bcftools concat -a -O z -o data/work/$1/S07604715/concordant_calls/stripped_variants_concat.vcf.gz data/work/$1/S07604715/concordant_calls/0000.vcf.gz data/work/$1/S07604715/concordant_calls/0001.rename.reorder.fix.vcf.gz data/work/$1/S07604715/concordant_calls/0002.vcf.gz data/work/$1/S07604715/concordant_calls/0003.rename.reorder.vcf.gz

bcftools sort -O z -m 5G -o data/work/$1/S07604715/concordant_calls/stripped_variants_concat.sort.vcf.gz data/work/$1/S07604715/concordant_calls/stripped_variants_concat.vcf.gz

bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta data/work/$1/S07604715/concordant_calls/stripped_variants_concat.sort.vcf.gz | bcftools norm -m-both -O z -o data/work/$1/S07604715/concordant_calls/stripped_variants_concat.sort.norm.vcf.gz
tabix -f -p vcf data/work/$1/S07604715/concordant_calls/stripped_variants_concat.sort.norm.vcf.gz

table_annovar.pl --buildver hg19 --vcfinput --outfile data/work/$1/S07604715/concordant_calls/stripped_variants_concat.sort.norm --protocol refGene,cytoband,gwasCatalog,genomicSuperDups,dbscsnv11,dbnsfp33a,popfreq_max_20150413,exac03,exac03nontcga,gnomad_exome,avsnp150,cosmic84_coding,cosmic84_noncoding,clinvar_20180603 --operation g,r,r,r,f,f,f,f,f,f,f,f,f,f -remove data/work/$1/S07604715/concordant_calls/stripped_variants_concat.sort.norm.vcf.gz $HOME/resources/annovar/humandb/

python annotated_somatic_vcf_v4.py --no_filter -I data/work/$1/S07604715/concordant_calls/stripped_variants_concat.sort.norm.hg19_multianno.vcf --tumor_id $1 --normal_id $2

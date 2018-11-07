

#for each vcf file
bcftools norm -f "$HOME"/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta {caller}.vcf.gz bcftools norm -m-both -O z -o {normalized_caller}.vcf.gz

bcftools isec -n +2 -f PASS -p data/work/{tumor}/{targets}/concordant_calls/ [inputs] > {output.sites.txt}

bcftools isec -n -1 -f PASS -p data/work/{tumor}/{targets}/private_calls/ [inputs] > output.sites.txt


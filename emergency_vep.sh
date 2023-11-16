
input=$1
DIR=`dirname $1`
base1=`basename $1 .bcf`

vcf1="${DIR}/${base1}.vcf"
vcf2="${DIR}/${base1}.vep.vcf"
fa1="${DIR}/ffpe_scenario.reference.fa"
fa2="${DIR}/ffpe_scenario.mutated.fa"


bcftools view -O v -o $vcf1 $input

vep -i $vcf1 -o $vcf2 \
--force_overwrite \
--offline \
--cache \
--format vcf \
--vcf \
--everything \
--canonical \
--assembly GRCh38 \
--species homo_sapiens \
--fasta $HOME/resources/Genomes/Human/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
--vcf_info_field ANN \
--plugin NMD \
--plugin ProteinSeqs,"$fa1","$fa2" \
--plugin Downstream \
--plugin REVEL,"$HOME/.vep/revel/revel_grch38.tsv.gz" \
--plugin SpliceAI,snv="$HOME/.vep/spliceai/spliceai_scores.raw.snv.hg38.vcf.gz",indel="$HOME/.vep/spliceai/spliceai_scores.raw.indel.hg38.vcf.gz" \
--plugin gnomADc,"$HOME/.vep/gnomAD/gnomad.v3.1.1.hg38.genomes.gz" \
--plugin UTRannotator,"$HOME/.vep/Plugins/UTRannotator/uORF_5UTR_GRCh38_PUBLIC.txt" \
--plugin LoF,loftee_path:"$HOME/.vep/Plugins/loftee" \
--custom "$HOME/.vep/clinvar/vcf_GRCh38/clinvar.vcf.gz",ClinVar,vcf,exact,0,CLNSIG,CLNREVSTAT,CLNDN

bgzip $vcf2
tabix -fp vcf ""$vcf2".gz

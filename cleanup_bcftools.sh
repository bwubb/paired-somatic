
name=$1
lib=$2

#bcftools
mkdir -p data/final/$lib/$name/bcftools
if [ ! -e data/final/$lib/$name/bcftools/$name.somatic.sites.txt ]
then
    #bgzip -c data/work/$lib/$name/bcftools/somtic/0000.vcf > data/final/$lib/$name/bcftools/$name.somatic.0000.vcf.gz
    #bgzip -c data/work/$lib/$name/bcftools/somatic/0001.vcf > data/final/$lib/$name/bcftools/$name.somatic.0001.vcf.gz
    #bgzip -c data/work/$lib/$name/bcftools/somatic/0002.vcf > data/final/$lib/$name/bcftools/$name.somatic.0002.vcf.gz
    #bgzip -c data/work/$lib/$name/bcftools/somatic/0003.vcf > data/final/$lib/$name/bcftools/$name.somatic.0003.vcf.gz
    #rsync -vr data/work/$lib/$name/bcftools/somatic/sites.txt data/final/$lib/$name/bcftools/$name.somatic.sites.txt
    echo $lib $name
fi


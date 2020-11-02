

##clean up file check
name=$1
lib=$2


#mutect2
if [ -e data/final/$lib/$name/mutect2/$name.somatic.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/mutect2/$name.somatic.raw.vcf.gz
    #echo "data/final/$lib/$name/mutect2/$name.somatic.raw.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/mutect2/$name.somatic.raw.vcf.gz"
fi


#manta
if [ -e data/final/$lib/$name/manta/$name.candidateSmallIndels.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/manta/$name.candidateSmallIndels.vcf.gz
    #echo "data/final/$lib/$name/manta/$name.candidateSmallIndels.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/manta/$name.candidateSmallIndels.vcf.gz"
fi

if [ -e data/final/$lib/$name/manta/$name.candidateSV.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/manta/$name.candidateSV.vcf.gz
    #echo "data/final/$lib/$name/manta/$name.candidateSV.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/manta/$name.candidateSV.vcf.gz"
fi

if [ -e data/final/$lib/$name/manta/$name.diploidSV.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/manta/$name.diploidSV.vcf.gz
    #echo "data/final/$lib/$name/manta/$name.diploidSV.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/manta/$name.diploidSV.vcf.gz"
fi

if [ -e data/final/$lib/$name/manta/$name.somaticSV.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/manta/$name.somaticSV.vcf.gz
    #echo "data/final/$lib/$name/manta/$name.somaticSV.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/manta/$name.somaticSV.vcf.gz"
fi


#strelka2
if [ -e data/final/$lib/$name/strelka2/$name.somatic.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/strelka2/$name.somatic.raw.vcf.gz
    #echo "data/final/$lib/$name/strelka2/$name.somatic.raw.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/strelka2/$name.somatic.raw.vcf.gz"
fi


#vardict

if [ -e data/final/$lib/$name/vardict/$name.variants.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/vardict/$name.variants.raw.vcf.gz
else
    echo Missing "data/final/$lib/$name/vardict/$name.variants.raw.vcf.gz"
fi

if [ -e data/final/$lib/$name/vardict/$name.somatic.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/vardict/$name.somatic.raw.vcf.gz
else
    echo Missing "data/final/$lib/$name/vardict/$name.somatic.raw.vcf.gz"
fi

if [ -e data/final/$lib/$name/vardict/$name.loh.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/vardict/$name.loh.raw.vcf.gz
else
    echo Missing "data/final/$lib/$name/vardict/$name.loh.raw.vcf.gz"
fi

if [ -e data/final/$lib/$name/vardict/$name.germline.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/vardict/$name.germline.raw.vcf.gz
else
    echo Missing "data/final/$lib/$name/vardict/$name.germline.raw.vcf.gz"
fi

#varscan2
if [ -e data/final/$lib/$name/varscan2/$name.germline.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.germline.raw.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.germline.raw.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.germline.raw.snp.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.germline.raw.snp.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.germline.raw.snp.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.germline.raw.indel.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.germline.raw.indel.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.germline.raw.indel.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.loh.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.loh.raw.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.loh.raw.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.loh.raw.snp.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.loh.raw.snp.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.loh.raw.snp.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.loh.raw.indel.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.loh.raw.indel.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.loh.raw.indel.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.somatic.raw.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.somatic.raw.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.somatic.raw.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.somatic.raw.snp.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.somatic.raw.snp.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.somatic.raw.snp.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.somatic.raw.indel.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.somatic.raw.indel.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.somatic.raw.indel.vcf.gz"
fi

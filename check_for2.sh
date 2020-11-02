

##clean up file check
name=$1
lib=$2


#mutect2
if [ -e data/final/$lib/$name/mutect2/$name.somatic.filtered.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/mutect2/$name.somatic.filtered.norm.clean.std.vcf.gz
    #echo "data/final/$lib/$name/mutect2/$name.somatic.raw.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/mutect2/$name.somatic.filtered.norm.clean.std.vcf.gz" 
fi


#strelka2
if [ -e data/final/$lib/$name/strelka2/$name.somatic.raw.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/strelka2/$name.somatic.raw.norm.clean.std.vcf.gz
    #echo "data/final/$lib/$name/strelka2/$name.somatic.raw.vcf.gz" check.
else
    echo Missing "data/final/$lib/$name/strelka2/$name.somatic.raw.norm.clean.std.vcf.gz"
fi


#vardict

if [ -e data/final/$lib/$name/vardict/$name.somatic.twice_filtered.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/vardict/$name.somatic.twice_filtered.norm.clean.std.vcf.gz
else
    echo Missing "data/final/$lib/$name/vardict/$name.somatic.twice_filtered.norm.clean.std.vcf.gz"
fi

if [ -e data/final/$lib/$name/vardict/$name.germline.twice_filtered.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/vardict/$name.germline.twice_filtered.norm.clean.std.vcf.gz
else
    echo Missing "data/final/$lib/$name/vardict/$name.germline.twice_filtered.norm.clean.std.vcf.gz"
fi

if [ -e data/final/$lib/$name/vardict/$name.loh.twice_filtered.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/vardict/$name.loh.twice_filtered.norm.clean.std.vcf.gz
else
    echo Missing "data/final/$lib/$name/vardict/$name.loh.twice_filtered.norm.clean.std.vcf.gz"
fi

#varscan2

if [ -e data/final/$lib/$name/varscan2/$name.somatic.fpfilter.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.somatic.fpfilter.norm.clean.std.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.somatic.fpfilter.norm.clean.std.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.germline.fpfilter.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.germline.fpfilter.norm.clean.std.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.germline.fpfilter.norm.clean.std.vcf.gz"
fi

if [ -e data/final/$lib/$name/varscan2/$name.loh.fpfilter.norm.clean.std.vcf.gz ]
then
    tabix -f -p vcf data/final/$lib/$name/varscan2/$name.loh.fpfilter.norm.clean.std.vcf.gz
else
    echo Missing "data/final/$lib/$name/varscan2/$name.loh.fpfilter.norm.clean.std.vcf.gz"
fi

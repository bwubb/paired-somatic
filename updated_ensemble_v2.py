

import vcfpy
import operator
import argparse

#header.add_filter_line(vcfpy.OrderedDict([('ID', 'DP10'), ('Description', 'total DP < 10')]))

#Read sites.txt
#vcf reader get enries, test this out
#merge records
#write record

#Mutect2 is not easy to fine germline calls that pass filtering.

def set_status(info):
    if info.get('SOMATIC',False):
        return 'somatic'
    elif info.get('STATUS',False):
        return info['STATUS'].lower()#Its probably a list
    elif info.get('SS',False):
        return {'0':'Reference','1':'Germline','2':'Somatic','3':'LOH','5':'Unknown'}[info['SS']].lower()
    else:
        return 'somatic'#set as default for now, mutect doesnt have anything.

def set_build_info(record,name):
    return vcfpy.OrderedDict([('STATUS',[f"{':'.join([set_status(record.INFO),name])}"]),('CALLER',[name])])

def set_build_call(call,tumor,normal):
    #Fuck you Varscan
    if type(call.data.get('AD',False))==int:
        #print('VarScan call')
        RD=call.data['RD']
        AD=call.data['AD']
        call.data['AD']=[RD,AD]
    #Strelka2 is worse and needs to be examined/compared to other callers where in common
    #data=vcfpy.OrderedDict([('GT',call.data.get['GT']),('DP',sum(call.data['AD'])),('AD',call.data['AD'])])
    data=vcfpy.OrderedDict([('GT',call.data.get('GT','./.')),('DP',sum(call.data.get('AD',[-1]))),('AD',call.data.get('AD',[-1,-1]))])#strelka ignore
    if call.sample=='TUMOR':
        return vcfpy.Call(tumor,data)
    elif call.sample=='NORMAL':
        return vcfpy.Call(normal,data)
    else:
        return vcfpy.Call(call.sample,data)

#Strelka info extraction.
#Somatic SNVs:

#refCounts = Value of FORMAT column $REF + “U” (e.g. if REF="A" then use the value in FOMRAT/AU)
#altCounts = Value of FORMAT column $ALT + “U” (e.g. if ALT="T" then use the value in FOMRAT/TU)
#tier1RefCounts = First comma-delimited value from $refCounts
#tier1AltCounts = First comma-delimited value from $altCounts
#Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)

#Somatic indels:

#tier1RefCounts = First comma-delimited value from FORMAT/TAR
#tier1AltCounts = First comma-delimited value from FORMAT/TIR
#Somatic allele freqeuncy is $tier1AltCounts / ($tier1AltCounts + $tier1RefCounts)



def main(argv=None):
    tumor=argv.tumor
    normal=argv.normal
    lib=argv.lib
    with open('data/work/{tumor}/{lib}/concordant_calls/sites.txt'.format(tumor=tumor,lib=lib),'r') as file:
        #head=[next(file).rstrip().split('\t') for x in range(1000)]
        head=[x.split('\t') for x in file.read().splitlines()]#I should change this
    
    #Make your own header.vcf
    #then add header_line for tumor_sample and normal_sample
    #Add header line for ##reference=file:/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta\n
    with vcfpy.Reader.from_path('/home/bwubb/resources/Vcf_files/ensemble-header.GRCh37.vcf') as reader:
        pass
    
#    tumor='TB5050-T1'
#    normal='TB5050'
    
    writer=vcfpy.Writer.from_path('data/work/{tumor}/{lib}/concordant_calls/somatic.n2.pass.vcf.gz'.format(tumor=tumor,lib=lib),vcfpy.Header(reader.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    #writer header=vcfpy.Header(header_vcf.header.lines,vcfpy.SampleInfos([tumor.name,normal.name]))
    
    mutect_reader=vcfpy.Reader.from_path('data/work/{tumor}/{lib}/mutect/somatic.twice_filtered.norm.vcf.gz'.format(tumor=tumor,lib=lib))#PASS=Somatic
    strelka_reader=vcfpy.Reader.from_path('data/work/{tumor}/{lib}/strelka/somatic.raw.norm.vcf.gz'.format(tumor=tumor,lib=lib))#SOMATIC in INFO, no LOH
    varscan_reader=vcfpy.Reader.from_path('data/work/{tumor}/{lib}/varscan/somatic.fpfilter.norm.vcf.gz'.format(tumor=tumor,lib=lib))#SOMATIC in INFO, but also SS in INFO. Look for SPV and SSC
    vardict_reader=vcfpy.Reader.from_path('data/work/{tumor}/{lib}/vardict/somatic.twice_filtered.norm.vcf.gz'.format(tumor=tumor,lib=lib))#Status=text
    
    #Can I read this from a yaml or json?
    #change this to be properly ordered
    reader_dict={0:{'reader':mutect_reader,'name':'mutect2'},1:{'reader':strelka_reader,'name':'strelka2'},2:{'reader':vardict_reader,'name':'vardictjava'},3:{'reader':varscan_reader,'name':'varscan2'}}
    #1101
    serialize=operator.methodcaller('serialize')
    for line in head:#line is already a list
        N=list(line[-1])
        build=None
        for i,j in enumerate(N):
            if j=='1':#I didnt check the PASS ever
                try:
                    for raw in reader_dict[i]['reader'].fetch('{0}:{1}-{1}'.format(*line[:2])): #chr,pos
                        '''
                        Fetch has less intuitive behavior.
                        It will jump to the position and limit the iterable to the end.
                        If no end is supplied it will carry to the end of the vcf.
                        '''
                        if any([raw.REF!=line[2],','.join(map(serialize,raw.ALT))!=line[3]]):#REF ALTcheck
                            continue
                        if build:
                            if any(C.data.get('GT','./.')=='./.' for C in build.calls):#strelka correction
                                #print('Correcting',build.calls,build.CHROM,build.POS)
                                #I looked more at getting the GT DP AND AD, the latter two being very annoying.
                                #two tiers of analysis and you have to parse different key:vals depending on which is present.
                                FILTER=['PASS']
                                FORMAT=['GT','DP','AD']
                                CALLS=[set_build_call(C,tumor,normal) for C in raw.calls]
                                QUAL='.'
                                INFO=vcfpy.OrderedDict([('STATUS',[f"{':'.join(['strelka','somatic'])}",f"{':'.join([reader_dict[i]['name'],set_status(raw.INFO)])}"]),('CALLER',['strelka2',reader_dict[i]['name']])])
                                build=vcfpy.Record(raw.CHROM,raw.POS,raw.ID,raw.REF,raw.ALT,QUAL,FILTER,INFO,FORMAT,CALLS)
                                #build.calls=[set_build_call(C,tumor,normal) for C in raw.calls]
                                #print('New calls',build.calls)
                            else:
                                build.INFO['CALLER']+=[reader_dict[i]['name']]
                                build.INFO['STATUS']+=[f"{':'.join([reader_dict[i]['name'],set_status(raw.INFO)])}"]
                        else:
                            FILTER=['PASS']
                            INFO=set_build_info(raw,reader_dict[i]['name'])
                            FORMAT=['GT','DP','AD']
                            CALLS=[set_build_call(C,tumor,normal) for C in raw.calls]
                            QUAL='.'
                            build=vcfpy.Record(raw.CHROM,raw.POS,raw.ID,raw.REF,raw.ALT,QUAL,FILTER,INFO,FORMAT,CALLS)
                except ValueError:
                    print(f"ValueError: IN {reader_dict[i]['reader'].path} AT '{line[0]}:{line[1]}-{line[1]}'")
                    print("Skipping variant")
                    break
        #try except into break?
        #Do I have a reason to break?
        #ValueError: could not create iterator for region '12:10602624-10602624'
        else:#Always wanted to use else with for loop
            #if any(C.data.get('GT','./.')=='./.' for C in build.calls):
            #    print('Writing bad call',build.CHROM,build.POS,build.calls)
            writer.write_record(build)
            #if str(build.POS)=='21186881':
            #    print(build.calls)
        
    writer.close()

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('-t','--tumor',help='Tumor name')
    p.add_argument('-n','--normal',help='Normal name')
    p.add_argument('--lib',default='SureSelect-Exon_v6+COSMIC',help='Lib name in /work dir')
    argv=p.parse_args()
    print('Arguments Initialized...')
    #LOG the run conditions
    for k,v in vars(argv).items():
        print(k,':',v)
    main(argv)
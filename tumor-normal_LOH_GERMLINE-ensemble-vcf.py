import vcfpy
import operator
import argparse
import csv
from collections import defaultdict

#header.add_filter_line(vcfpy.OrderedDict([('ID', 'DP10'), ('Description', 'total DP < 10')]))

#Read sites.txt
#vcf reader get enries, test this out
#merge records
#write record

def set_build_calls(data,sample):
    temp_GT=[x['GT'] for x in data]
    _DP=max([x['DP'] for x in data])
    _AD=max([x['AD'] for x in data])
    _AAF=max([x['AAF'] for x in data])
    if len(set(temp_GT))>2:
        print('Discordant genotypes',temp_GT)
        _GT='./.'
    else:
        _GT=temp_GT[0]
    return vcfpy.Call(sample,vcfpy.OrderedDict([('GT',_GT),('DP',_DP),('AD',_AD),('AAF',_AAF)]))

def get_Args():
    p=argparse.ArgumentParser()
    p.add_argument('-T','--tumor',help='Tumor name')
    p.add_argument('-N','--normal',help='Normal name')
    p.add_argument('-L','--lib',default='S04380110',help='Lib name in /work dir')
    argv=p.parse_args()
    print('Arguments Initialized...')
    #LOG the run conditions
    for k,v in vars(argv).items():
        print(k,':',v)
    return argv

def main(argv=None):
    try:
        snakemake
        #samples file in params
    except NameError:
        args=get_Args()
        tumor=args.tumor
        normal=args.normal
        lib=args.lib
    else:
        #Better snakemake
        tumor=snakemake.params['tumor']
        normal=snakemake.params['normal']
        lib=snakemake.params['lib']
    
    with vcfpy.Reader.from_path('/home/bwubb/resources/Vcf_files/loh_germline.ensemble-header.GRCh37.20200102.vcf') as _VCFH:
        pass
    #Need custom lines
    #Make your own header.vcf
    #files
    #script version
    #then add header_line for tumor_sample and normal_sample
    #Add header line for ##reference=file:/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta\n
    
    #high_conf_writer=vcfpy.Writer.from_path(f'data/work/{lib}/{tumor}/bcftools/somatic.high_confidence.vcf.gz',vcfpy.Header(_VCFH.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    #low_conf_writer=vcfpy.Writer.from_path(f'data/work/{lib}/{tumor}/bcftools/somatic.low_confidence.vcf.gz',vcfpy.Header(_VCFH.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    #failed_writer=vcfpy.Writer.from_path(f'data/work/{lib}/{tumor}/bcftools/somatic.failed_confidence.vcf.gz',vcfpy.Header(_VCFH.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    outfile_writer=vcfpy.Writer.from_path(f'data/work/{lib}/{tumor}/bcftools/loh_germline/loh_germline_variants.vcf.gz',vcfpy.Header(_VCFH.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    ####READERS####
    vardict_reader=vcfpy.Reader.from_path(f'data/work/{lib}/{tumor}/vardict/loh_germline.twice_filtered.norm.clean.std.vcf.gz')
    varscan_reader=vcfpy.Reader.from_path(f'data/work/{lib}/{tumor}/varscan2/loh_germline.fpfilter.norm.clean.std.vcf.gz')
    #Can I read this from a yaml or json?
    reader_dict={0:{'reader':vardict_reader,'name':'vardict'},1:{'reader':varscan_reader,'name':'varscan2'}}
    #1101
    serialize=operator.methodcaller('serialize')
    with open(f'data/work/{lib}/{tumor}/bcftools/loh_germline/sites.txt','r') as file:
        sites_reader=csv.DictReader(file,delimiter='\t',fieldnames=['CHROM','POS','REF','ALT','BIN'])
        for row in sites_reader:
            bin=list(row['BIN'])
            FORMAT=['GT','DP','AD','AAF']
            build=defaultdict(list)
            for i,b in enumerate(bin):
                if b=='1':
                    for record in reader_dict[i]['reader'].fetch(f"{row['CHROM']}:{row['POS']}-{row['POS']}"):
                        #Fetch has less intuitive behavior.
                        #It will jump to the position and limit the iterable to the end.
                        #If no end is supplied it will carry to the end of the vcf.
                        if any([record.REF!=row['REF'],record.ALT[0].value!=row['ALT']]):#REF ALTcheck
                            #Good place for custom exceptions
                            print(f"REF ALT do not match: {row['CHROM']} {row['POS']} {record.REF}:{row['REF']};{record.ALT[0].value}:{row['ALT']}")
                            continue
                        #INFO
                        build[record.INFO['CALLER'].upper()]=record.INFO['OFS']#Is this a list?
                        build['OFS']+=record.INFO['OFS']
                        #Check if calls exist and if any FORMAT if tumor DP is 0.
                        build['SS']+=[record.INFO['SS']]
                        if build[tumor]==[]:
                            for call in record.calls:
                                data=vcfpy.OrderedDict([(k,call.data[k]) for k in FORMAT])
                                build[call.sample]+=[data]
            else:#Always wanted to use else with for loop, Will likely switch to try, except though
                ### INFO ###
                INFO=vcfpy.OrderedDict([(caller,build.get(caller,'.')) for caller in ['VARDICTJAVA','VARSCAN2']])
                INFO['STATUS']='|'.join(build.get('SS'))
                #### CALLS ####
                CALLS=[set_build_calls(build[tumor],tumor),set_build_calls(build[normal],normal)]
                ##### QUAL FILTER#####
                QUAL='.'
                FILTER='.'
                out_record=vcfpy.Record(record.CHROM,record.POS,record.ID,record.REF,record.ALT,QUAL,FILTER,INFO,FORMAT,CALLS)
                
                outfile_writer.write_record(out_record)
                #if build['OFS'].count('PASS')>=2:
                #    high_conf_writer.write_record(out_record)
                #elif build['OFS'].count('PASS')==1:
                #    low_conf_writer.write_record(out_record)
                #else:
                #    failed_writer.write_record(out_record)
    outfile_writer.close()
    #high_conf_writer.close()
    #low_conf_writer.close()
    #failed_writer.close()

if __name__=='__main__':
    main()
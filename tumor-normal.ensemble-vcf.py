import vcfpy
import operator
import argparse
import csv
import os
from collections import defaultdict

#vcfpy HeaderLine class
##fileformat=VCFv4.3
##fileDate=20090805
##source=myImputationProgramV3.1
##reference=file:///seq/references/1000GenomesPilot-NCBI36.fasta

#Read sites.txt
#vcf reader get enries, test this out
#merge records
#write record

class ConcordantCall(object):
    def __init__(self,fields,tumor,normal):
        '''
        from bcftools sites.txt
        fieldnames=['CHROM','POS','REF','ALT','BIN']
        '''
        self.CHROM=fields['CHROM']
        self.POS=fields['POS']
        self.REF=fields['REF']
        self.ALT=self.AltRecord(fields['REF'],fields['ALT'])
        self.BIN=fields['BIN']
        self.MUTECT2='.'
        self.STRELKA2='.'
        self.VARDICT='.'
        self.VARSCAN2='.'
        self.NumPASS=0
        self.ID='.'
        self.QUAL='.'
        self.FILTER=['.']
        self.INFO=vcfpy.OrderedDict([])
        self.FORMAT=['GT','DP','AD','AAF']
        self.GT_CONFLICT='.'
        self.GT_check={}
        self.tumor=vcfpy.Call(tumor,vcfpy.OrderedDict([('GT','./.'),('DP',0),('AD',[0,0]),('AAF',0.0)]))
        self.normal=vcfpy.Call(normal,vcfpy.OrderedDict([('GT','./.'),('DP',0),('AD',[0,0]),('AAF',0.0)]))
        #May need one for Tumor and Normal, record conflicts in both?
    
    @staticmethod
    def AltRecord(REF,ALT):
        '''
        retun list of vcfpy.Substitution
        depending on what REF and ALT look like
        Work on DUP MNV later
        '''
        if all([len(REF)==1,len(ALT)==1]):
            return [vcfpy.Substitution('SNV',ALT)]
        elif len(REF)>len(ALT):
            return [vcfpy.Substitution('DEL',ALT)]
        elif len(REF)<len(ALT):
            return [vcfpy.Substitution('INS',ALT)]
        elif len(REF)==len(ALT):
            return [vcfpy.Substitution('MNV',ALT)]
            #This may be incorrect; not sure how MNV is written in vcf format
        else:
            raise
    
    def update(self,caller,ofs):
        if getattr(self,caller,'.')=='.' and ofs==['PASS']:
            self.NumPASS+=1
        setattr(self,caller,ofs)
    
    def add_gt(self,caller,calls):
        self.GT_check[caller]=calls[0].data['GT']
        #If DP is highest and AD[1]>0
        ##Not sure if this is a good check
        if calls[0].data['DP']>self.tumor.data['DP'] and calls[0].data['AD'][1]>0:
            for k in self.FORMAT:
                self.tumor.data[k]=calls[0].data[k]
                self.normal.data[k]=calls[1].data[k]
        elif self.tumor.data['GT']=='0/0' and calls[0].data['GT']!='0/0':
            #Overwrite HOM_REF is possible
            for k in self.FORMAT:
                self.tumor.data[k]=calls[0].data[k]
                self.normal.data[k]=calls[1].data[k]
        else:
            pass
    
    def set_GT_CONFLICT(self):
        if len(set(self.GT_check.values()))>1:
            d=defaultdict(list)
            for k,v in self.GT_check.items():
                d[v].append(k)
            #A,B:C
            self.GT_CONFLICT=':'.join([','.join(v) for v in d.values()])
        else:
            self.GT_CONFLICT='NO_CONFLICT'
    
    def set_filter(self,n):
        if self.NumPASS>=n:
            self.FILTER=['PASS']
        else:
            self.FILTER=['REJECT']
    
    def set_INFO(self):
        INFO=vcfpy.OrderedDict([(x,getattr(self,x,'.')) for x in ['MUTECT2','STRELKA2','VARDICT','VARSCAN2']])
        INFO['NumPASS']=self.NumPASS
        INFO['GT_CONFLICT']=[str(self.GT_CONFLICT)]
        self.INFO=INFO
        
    def outrecord(self):
        self.set_GT_CONFLICT()
        self.set_INFO()
        return vcfpy.Record(self.CHROM,int(self.POS),self.ID,self.REF,self.ALT,self.QUAL,self.FILTER,self.INFO,self.FORMAT,[self.tumor,self.normal])

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

def parse_arguments():
    #later on add header options,
    ##library, other things
    ##META
    ##SAMPLE
    p=argparse.ArgumentParser()
    p.add_argument('-b','--sites',help='bcftools isec sites.txt')
    p.add_argument('-T','--tumor',help='Tumor name')
    p.add_argument('-N','--normal',help='Normal name')
    p.add_argument('-o','--out_fp',help='Outfile path and name')
    p.add_argument('--mutect2_vcf',help='path to mutect2.norm.clean.std.vcf.gz')
    p.add_argument('--strelka2_vcf',help='path to strelka2.norm.clean.std.vcf.gz')
    p.add_argument('--vardict_vcf',help='path to vardict.norm.clean.std.vcf.gz')
    p.add_argument('--varscan2_vcf',help='path to varscan2.norm.clean.std.vcf.gz')
    p.add_argument('-L','--lib',default='S04380110',help='Lib name in work dir')
    #p.add_Argument('--yaml',help='yaml info')
    return vars(p.parse_args())

def main(argv=None):
    for k,v in argv.items():
        print(k,':',v)
    with vcfpy.Reader.from_path('/home/bwubb/resources/Vcf_files/ensemble-header.GRCh37.20200716_v3.vcf') as VCFH:
        pass
    #files
    #script version
    #then add header_line for tumor_sample and normal_sample
    #Add header line for ##reference=file:human_g1k_v37.fasta\n
    #Others mentioned above
    tumor=argv['tumor']
    normal=argv['normal']
    OUTFILE=argv['out_fp']
    passed_writer=vcfpy.Writer.from_path(OUTFILE,vcfpy.Header(VCFH.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    #failed_writer=vcfpy.Writer.from_path(f'FAILFILE',vcfpy.Header(VCFH.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    ####READERS####
    '''
    Here we look for the files given.
    If the file exists, a reader object is made, eg. mutect2_reader
    At the end we define a dictionary of reader objects.
    ##yaml option??##
    '''
    AVAILABLE_METHODS=['mutect2','strelka2','vardict','varscan2']
    CALLERS=[]
    for CALLER in AVAILABLE_METHODS:
        if eval(f'argv["{CALLER}_vcf"]')!=None:
            assert(os.path.isfile(eval(f'argv["{CALLER}_vcf"]')))
            #This assert can be replaced with FileNotFound or something. Or not.
            exec(f'{CALLER}_reader=vcfpy.Reader.from_path(argv["{CALLER}_vcf"])')
            CALLERS.append(CALLER)
    reader_dict={}
    for i,n in enumerate(CALLERS):
        reader_dict[i]={'reader':eval(f'{n}_reader'),'name':n}
    #reader_dict={0:{'reader':mutect_reader,'name':'mutect2'},1:{'reader':strelka_reader,'name':'strelka2'},2:{'reader':vardict_reader,'name':'vardict'},3:{'reader':varscan_reader,'name':'varscan2'}}
    #1101
    with open(argv['sites'],'r') as file:
        sites_reader=csv.DictReader(file,delimiter='\t',fieldnames=['CHROM','POS','REF','ALT','BIN'])
        for row in sites_reader:
            bin=list(row['BIN'])
            Build=ConcordantCall(row,tumor,normal)
            for i,b in enumerate(bin):
                if b=='1':
                    #Can I yield these?
                    for record in reader_dict[i]['reader'].fetch(f"{row['CHROM']}:{row['POS']}-{row['POS']}"):
                        if any([record.REF!=row['REF'],record.ALT[0].value!=row['ALT']]):#REF ALTcheck
                            continue
                        #Good place for custom exceptions
                        #Verbose output
                        Build.update(record.INFO['CALLER'].upper(),record.INFO['OFS'])
                        Build.add_gt(record.INFO['CALLER'].upper(),record.calls)
                        #Check if calls exist and if any FORMAT if tumor DP is 0.
            else:#Always wanted to use else with for loop
                Build.set_filter(2)#Hard coded to 2 for now
                passed_writer.write_record(Build.outrecord())
#                if Build.FILTER=='PASS':
#                    passed_writer.write_record(Build.outrecord())
#                else:
#                    failed_writer.write_record(Build.outrecord())
    passed_writer.close()
    #failed_writer.close()
    #I guess I should use subprocess, or at least os.popen?
    os.system(f'bgzip -f {OUTFILE}')
    os.system(f'tabix -fp vcf {OUTFILE}.gz')

if __name__=='__main__':
    try:
        snakemake
    except NameError:
        main(parse_arguments())
    else:
        main(parse_snakemake())
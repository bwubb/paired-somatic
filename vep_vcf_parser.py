
#VEP parser
import vcfpy
import os
import csv
import argparse
from collections import defaultdict

#class
class VEPannotation(object):
    def __init__(self,vcfpy_record,tumor_normal=False,tumor_id=None):
        self.fields=defaultdict(str)
        self.fields['Chr']=f"{vcfpy_record.CHROM}"
        self.fields['Start']=f"{vcfpy_record.POS}"
        self.fields['REF']=f"{vcfpy_record.REF}"
        self.fields['ALT']=f"{vcfpy_record.ALT[0].value}"
        self.fields['FILTER']=f"{';'.join(vcfpy_record.FILTER)}".replace("MONOALLELIC",'.')
        if tumor_normal:
            self.calls=self.__tumor_normal__(vcfpy_record.calls,tumor_id)
        else:
            self.calls=self.__get_calls__(vcfpy_record.calls)
        self.call_count=len(self.calls)

    def in_region(self,intervals):
        #integer start,end
        if any([start<=int(self.fields['Start'])<=end for start,end in intervals]):
            return True
        return False

    @staticmethod
    def __get_calls__(record_calls):
        #PMBB data has some errors with Ref/Alt genotypes split across two records
        #'AD': [None, 26],
        def AD_check(AD):
            for i,v in enumerate(AD):
                if v in [None,'None','.']:
                    AD[i]=0
            return AD

        def ZYG_check(GT,AD):
            if GT in ['.','./.']:
                return '.'
            elif '.' in GT:
                f=AD[1]/sum(AD)
                if f>=0.985:
                    return 'HOM_ALT'
                elif f<=0.015:
                    return 'HOM_REF'
                else:
                    return 'HET_ALT'
            else:
                #If it is not in these options, I want it reported as is.
                return {'0/0':'HOM_REF','0/1':'HET_ALT','1/0':'HET_ALT','1/1':'HOM_ALT'}.get(GT,GT)

        _calls=[]
        index_error=False
        zerodivision_error=False
        #One day I want to write failed output to a separate vcf.
        for c in record_calls:
            #Custom errors for metrics
            #This is all messed up in germline VLR.
            #There is an AF but no AD
            x=defaultdict(str)
            x['Sample.ID']=c.sample
            x['Sample.Depth']=f"{c.data.get('DP',0)}"
            x['Sample.AltDepth']='.'
            x['Sample.AltFrac']=f"{c.data.get('AF',[0.0])[0]}"
            x['Sample.Zyg']={'0.0':'HOM_REF','0.5':'HET_ALT','1.0':'HOM_ALT'}.get(x['Sample.AltFrac'],'.')
            _calls.append(x)
        return _calls

    @staticmethod
    def __tumor_normal__(record_calls,tumor_id):
        x=defaultdict(str)
        for c in record_calls:
            #Custom errors for metrics
            if c.sample.lower()=='tumor' or c.sample==tumor_id:
                x['Tumor.ID']=c.sample
                x['Tumor.Depth']=f"{c.data.get('DP','.')}"
                x['Tumor.Zyg']=f"{c.data.get('GT','./.').replace('/',';')}"
                x['Tumor.AltDepth']=f"{c.data.get('AD','.')}"
                x['Tumor.AltFrac']=f"{c.data.get('AF',['.'])[0]:.3f}"
                #VLR does not like to give AD/AF. This is an unchecked estimation.
                if x['Tumor.AltDepth']=='.' and x['Tumor.AltFrac']!='.':
                        x['Tumor.AltDepth']=f"{int(c.data.get('DP','0'))*float(c.data.get('AF',['0.0'])[0]):.0f}"
            else:
                x['Normal.ID']=c.sample
                x['Normal.Depth']=f"{c.data.get('DP','.')}"
                x['Normal.Zyg']=f"{c.data.get('GT','./.').replace('/',';')}"
                x['Normal.AltDepth']=f"{c.data.get('AD','.')}"
                x['Normal.AltFrac']=f"{c.data.get('AF',['.'])[0]:.3f}"
                if x['Normal.AltDepth']=='.' and x['Normal.AltFrac']!='.':
                        x['Normal.AltDepth']=f"{int(c.data.get('DP','0'))*float(c.data.get('AF',['0.0'])[0]):.0f}"
        return [x]

    @staticmethod
    def variant_type():
        ##Ways to get this
        #1 INFO Field information; preferred but non-standard for many vcf generation
        #2 commandline argument; if vcf files are split by type
        #3 filename; similar to above but less robust.
        #4 Conditional logic?
        return 'Imprecise_variant'

    def info(self,CSQ):
        #INFOS
        self.fields['Gene']=CSQ['SYMBOL']
        self.fields['Gene.Accession']=CSQ['Gene']
        #self.fields['Genomeic.Region']=
        self.fields['Variant.Class']=CSQ['VARIANT_CLASS']
        self.fields['Variant.Type']=self.variant_type()
        self.fields['Variant.Consequence']=CSQ['Consequence']
        self.fields['HGVSc']=CSQ['HGVSc'].split(':')[-1]
        self.fields['HGVSp']=CSQ['HGVSp'].split(':')[-1]
        self.fields['Feature.Type']=CSQ['Feature_type']
        self.fields['Feature.Accession']=CSQ['Feature']
        self.fields['Bio.type']=CSQ['BIOTYPE']
        self.fields['Existing.variation']=CSQ['Existing_variation'] #trim to rs?
        self.fields['EXON']=CSQ['EXON'].replace('/','|')
        self.fields['INTRON']=CSQ['INTRON'].replace('/','|')
        self.fields['STRAND']=CSQ['STRAND']
        self.fields['cDNA.position']=CSQ['cDNA_position']
        self.fields['CDS.position']=CSQ['CDS_position']
        self.fields['Protein.position']=CSQ['Protein_position']
        self.fields['Amino.acids']=CSQ['Amino_acids']
        self.fields['Codons']=CSQ['Codons']

    def splice_ai(self,CSQ):
        #SpliceAI
        self.fields['SpliceAI.DS_AG']=CSQ['SpliceAI_pred_DS_AG']
        self.fields['SpliceAI.DS_AL']=CSQ['SpliceAI_pred_DS_AL']
        self.fields['SpliceAI.DS_DG']=CSQ['SpliceAI_pred_DS_DG']
        self.fields['SpliceAI.DS_DL']=CSQ['SpliceAI_pred_DS_DL']

    def snv_prediction(self,CSQ):
        #Prediction
        self.fields['SIFT']=CSQ['SIFT']
        self.fields['PolyPhen']=CSQ['PolyPhen']
        self.fields['REVEL']=CSQ['REVEL']

    def gnomAD(self,CSQ):
        #gnomAD
        self.fields['gnomAD.AF']=CSQ['gnomAD_AF']
        self.fields['gnomAD.AFR']=CSQ['gnomAD_AFR_AF']
        self.fields['gnomAD.AMR']=CSQ['gnomAD_AMR_AF']
        self.fields['gnomAD.ASJ']=CSQ['gnomAD_ASJ_AF']
        self.fields['gnomAD.EAS']=CSQ['gnomAD_EAS_AF']
        self.fields['gnomAD.FIN']=CSQ['gnomAD_FIN_AF']
        self.fields['gnomAD.NFE']=CSQ['gnomAD_NFE_AF']
        self.fields['gnomAD.OTH']=CSQ['gnomAD_OTH_AF']
        self.fields['gnomAD.SAS']=CSQ['gnomAD_SAS_AF']
        self.fields['gnomAD.MAX_AF']=CSQ['MAX_AF']
        self.fields['gnomAD.MAX_POPS']=CSQ['MAX_AF_POPS'] #.replace('_','.')

    def clinvar(self,CSQ):
        #ClinVar
        self.fields['ClinVar']=CSQ['ClinVar']#ClinVar
        self.fields['ClinVar.SIG']=CSQ['ClinVar_CLNSIG']#ClinVar.SIG
        self.fields['ClinVar.REVSTAT']=CSQ['ClinVar_CLNREVSTAT']#ClinVar.REVSTAT
        self.fields['ClinVar.DN']=CSQ['ClinVar_CLNDN']#ClinVar.CLNDN

    def loftee(self):
        pass

    def utrannotator(self):
        pass

    def fill_values(self,header):
        for h in header:
            if self.fields[h]=='':
                self.fields[h]='.'

    def print(self,header=['Chr','Start','REF','ALT','FILTER']):
        return ','.join([self.fields[h] for h in header])

    def report(self,writer):
        #need other report if 'none' is used.
        for o in self.calls:
            writer.writerow({**self.fields,**o})
                #reverse dict order had incorrect fields.

def report_header(annotations,tumor_normal=False):
    if tumor_normal:
        header=['Tumor.ID','Normal.ID','Chr','Start','REF','ALT','FILTER']
    else:
        header=['Sample.ID','Chr','Start','REF','ALT','FILTER']
    header+=['Gene','Gene.Accession','Variant.Class','Variant.Type','Variant.Consequence']
    header+=['HGVSc','HGVSp','Feature.Type','Feature.Accession','Bio.type','Existing.variation']
    header+=['EXON','INTRON','STRAND','cDNA.position','CDS.position','Protein.position','Amino.acids','Codons']

    def splice_ai(tumor_normal=False):
        return ['SpliceAI.DS_AG','SpliceAI.DS_AL','SpliceAI.DS_DG','SpliceAI.DS_DL']

    def snv_predition(tumor_normal=False):
        #excluded SIFT and Polyphen
        return ['REVEL']

    def gnomAD(tumor_normal=False):
        return ['gnomAD.AF','gnomAD.AFR','gnomAD.AMR','gnomAD.ASJ','gnomAD.EAS','gnomAD.FIN','gnomAD.NFE','gnomAD.OTH','gnomAD.SAS','gnomAD.MAX_AF','gnomAD.MAX_POPS']

    def clinvar(tumor_normal=False):
        return ['ClinVar','ClinVar.SIG','ClinVar.REVSTAT','ClinVar.DN']

    def loftree(tumor_normal=False):
        pass

    def genotype(tumor_normal=False):
        if tumor_normal:
            return ['Tumor.Zyg','Tumor.Depth','Tumor.AltDepth','Tumor.AltFrac','Normal.Zyg','Normal.Depth','Normal.AltDepth','Normal.AltFrac']
        return ['Sample.Zyg','Sample.Depth','Sample.AltDepth','Sample.AltFrac']

    if 'none' in annotations:
        return header
    elif 'everything' in annotations:
        for F in [splice_ai,snv_predition,gnomAD,clinvar,genotype]:
            header+=F(tumor_normal)
    else:
        for x in ['splice_ai','snv_predition','gnomAD','clinvar','genotype']:
            if x in annotations:
                header+=eval(f"{x}({tumor_normal})")
    return header

def annotation_check(annotations):
    choices=['everything','snv_prediction','gnomAD','splice_ai','clinvar','genotype','none']
    for x in annotations:
        if x not in choices:
            raise NameError

def gene_list_check(fp):
    #receive error for NoneType
    try:
        if os.path.isfile(fp):
            with open(fp,'r') as file:
                gene_list=file.read().splitlines()
                #print(gene_list)
                return True,gene_list
    except TypeError:
        pass
    return False,[]

def bed_region_check(fp):
    #receive error for NoneType
    try:
        if os.path.isfile(fp):
            bed_regions=defaultdict(list)
            with open(fp,'r') as file:
                reader=csv.reader(file,delimiter='\t')
                for row in reader:
                    bed_regions[row[0]].append([int(row[1]),int(row[2])])
                return True,bed_regions
    except TypeError:
        pass
    return False,{}

def get_args(argv):
    #make single tumor_normal cohort as subparser?
    #scenrios may need some sort of sample id mapping file.
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input_vcf',help='Input: vcf file from ensembl-vep.')
    p.add_argument('-o','--output_csv',help='Output: csv file trimmed for specific design.')
    p.add_argument('-g','--gene_list',help='Optional list of genes to include only.')
    p.add_argument('-r','--bed_region',help='Bed file format to subset regions.')
    p.add_argument('-t','--variant_type',help='Optional value for column field.')
    p.add_argument('-m','--mode',default='cohort',help='Run mode determines how calls are reported. Single should be "single,{sample.id}. Tumor/Normal should be "tumor_normal,{tumor.id},{normal.id}"')
    p.add_argument('annotations',nargs=argparse.REMAINDER,default='everything',choices=['everything','snv_prediction','gnomAD','splice_ai','clinvar','genotype','none'],help='annotation blocks to include')
    return p.parse_args(argv)

def main(argv=None):
    args=get_args(argv)
    gene_filter,gene_list=gene_list_check(args.gene_list)
    region_filter,bed_regions=bed_region_check(args.bed_region)
    annotation_check(args.annotations)
    #need better mode check
    if args.mode.startswith('tumor_normal'):
        #def
        z=args.mode.split(',')
        tumor_normal=True
        tumor=z[1]
        normal=z[2]
        #index error default tumor normal
    else:
        tumor_normal=False
        tumor=None

    if args.mode.startswith('single'):
        #def
        z=args.mode.split(',')
        single=True
        sample=z[1]
    else:
        single=False
        sample=None
    #this is becoming a lot of if/else stuff.
    #Need to clean up this main()
    #default everything is not working properly
    header=report_header(args.annotations,tumor_normal)
    #test data
    #VcfReader=vcfpy.Reader.from_path('data/vcf/FLCN/PMBB-Release-2020-2.0_genetic_exome_FLCN_NF.norm.vep.vcf.gz')
    VcfReader=vcfpy.Reader.from_path(f'{args.input_vcf}')
    csq_keys=VcfReader.header.get_info_field_info('CSQ').serialize().split(' ')[-1].rstrip('">').split('|')
    with open(args.output_csv,'w') as outfile:
        writer=csv.DictWriter(outfile,fieldnames=header,delimiter=',',restval='.',extrasaction='ignore',quoting=csv.QUOTE_NONNUMERIC,dialect='excel')
        writer.writeheader()
        for record in VcfReader:
            if len(record.ALT)>1:
                print(f"Warning! : record.ALT length is {len(record.ALT)}. Not currently supported")
            vep_data=VEPannotation(record,tumor_normal,tumor)
            if 'CSQ' not in record.INFO:
                # This is from those fucking * calls in haplotype caller.
                #print(f'Missing CSQ: {record}')
                continue
            for csq_i in record.INFO['CSQ']:
                csq_dict=dict(zip(csq_keys,csq_i.split('|')))
                if csq_dict['SYMBOL']!='' and csq_dict['CANONICAL']=='YES':
                    #checks
                    vep_data.info(csq_dict)
                    vep_data.snv_prediction(csq_dict)
                    vep_data.splice_ai(csq_dict)
                    vep_data.gnomAD(csq_dict)
                    vep_data.clinvar(csq_dict)
                    #more checks
                    if tumor_normal:
                        vep_data.calls[0]['Tumor.ID']=tumor
                        vep_data.calls[0]['Normal.ID']=normal

                    if single:
                        vep_data.calls[0]['Sample.ID']=sample

                    vep_data.fill_values(header)
                    #print(vep_data.print(header))
                    if gene_filter and vep_data.fields['Gene'] in gene_list:
                        vep_data.report(writer)
                    elif region_filter and vep_data.in_region(bed_regions[record.CHROM]):
                        vep_data.report(writer)
                    elif not gene_filter and not region_filter:
                        vep_data.report(writer)
                    else:
                        continue

if __name__=='__main__':
    main()


#Need to recover
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  PMBB1495195945216
#17      17216394        17_17216394_T_TG        T       TG      641     .       AF=6.8e-05;AQ=641;AC=0;AN=1     GT:DP:AD:SBPV:GQ:PL:FT:RNC      ./0:50:24,0:0.3414:99:.:.:-.
#17      17216394        17_17216394_TG_T        TG      T       435     MONOALLELIC     AF=1.1e-05;AQ=435;AC=1;AN=1     GT:DP:AD:SBPV:GQ:PL:FT:RNC      ./1:50:.,26:.:435:.:.:1.

#Add check to make sure custom filter columns are there


#Existing_variation > Existing.variation [rs780588085]
#IMPACT ? [Look up what this value is from. WE determine impact]

#Variant call method coutn
#VLR filter
#LOFTEE ?
#Custom errors for get_calls to report rate of messed up AD fields, etc.

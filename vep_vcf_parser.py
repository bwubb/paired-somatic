
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
        #I could put all the callers stuff here.
        #There creating a few empty attributes on default is no big deal.

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
        #make a first attempt here, if it exists. Then later make a double check of AD/AF/DP
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
            #This really needs to be sorted, ZYG keeps throwing problems
            #Generally need to exclude HOM_REF
            #x['Sample.AltDepth']='.'
            #x['Sample.AltFrac']=f"{c.data.get('AF',[0.0])[0]}"
            #x['Sample.Zyg']={'0.0':'HOM_REF','0.5':'HET_ALT','1.0':'HOM_ALT'}.get(x['Sample.AltFrac'],'.')
            #_calls.append(x)
            #This keeps HOM_REF
            try:
                x['Sample.AltDepth']=f"{AD_check(c.data.get('AD',[0,0]))[1]}"
            except IndexError as e:
                #report error
                index_error=True
                x['Sample.AltDepth']='.'
            try:
                x['Sample.AltFrac']=f"{AD_check(c.data.get('AD',[0,0]))[1]/c.data.get('DP',0):.3f}"#round decimals
                x['Sample.Zyg']=f"{ZYG_check(c.data.get('GT','./.'),AD_check(c.data.get('AD',[0,0])))}"
            except (ZeroDivisionError,IndexError,TypeError) as e:
                #report error
                zerodivision_error=True
                x['Sample.AltFrac']='.'
                x['Sample.Zyg']='.'
            if x['Sample.Zyg'] not in ['.','HOM_REF']:
                #print(x)
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
                if type(c.data.get('AF',['.']))==list:#some callers have AF as a float rather than list of floats
                    x['Tumor.AltFrac']=f"{c.data.get('AF',['.'])[0]:.3f}"
                else:
                    x['Tumor.AltFrac']=f"{c.data['AF']:.3f}"
                #VLR does not like to give AD/AF. This is an unchecked estimation.
                if x['Tumor.AltDepth']=='.' and x['Tumor.AltFrac']!='.':
                        x['Tumor.AltDepth']=f"{int(c.data.get('DP','0'))*float(c.data.get('AF',['0.0'])[0]):.0f}"
                if x['Tumor.Zyg']=='.;.':
                    if 0.0<float(x['Tumor.AltFrac'])<1.0:
                        x['Tumor.Zyg']='HET'
                    elif float(x['Tumor.AltFrac'])==1.0:
                        x['Tumor.Zyg']='HOM_ALT'
                    else:
                        x['Tumor.Zyg']='HOM_REF'
            else:
                x['Normal.ID']=c.sample
                x['Normal.Depth']=f"{c.data.get('DP','.')}"
                x['Normal.Zyg']=f"{c.data.get('GT','./.').replace('/',';')}"
                x['Normal.AltDepth']=f"{c.data.get('AD','.')}"
                if type(c.data.get('AF',['.']))==list:
                    x['Normal.AltFrac']=f"{c.data.get('AF',['.'])[0]:.3f}"
                else:
                    x['Normal.AltFrac']=f"{c.data['AF']:.3f}"
                if x['Normal.AltDepth']=='.' and x['Normal.AltFrac']!='.':
                        x['Normal.AltDepth']=f"{int(c.data.get('DP','0'))*float(c.data.get('AF',['0.0'])[0]):.0f}"
                if x['Normal.Zyg']=='.;.':
                    if 0.0<float(x['Normal.AltFrac'])<1.0:
                        x['Normal.Zyg']='HET'
                    elif float(x['Normal.AltFrac'])==1.0:
                        x['Normal.Zyg']='HOM_ALT'
                    else:
                        x['Normal.Zyg']='HOM_REF'
        return [x]

    def info(self,CSQ):
        #INFOS
        self.fields['Gene']=CSQ['SYMBOL']
        self.fields['Gene.Accession']=CSQ['Gene']
        #self.fields['Genomeic.Region']=
        self.fields['Variant.Class']=CSQ['VARIANT_CLASS']
        self.fields['Variant.LoF_level']='.'
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
        self.fields['gnomAD.AF']=CSQ.get('gnomADg_AF',CSQ.get('gnomADe_AF','.'))
        self.fields['gnomAD.AFR']=CSQ.get('gnomADg_AFR_AF',CSQ.get('gnomADe_AFR_AF','.'))
        self.fields['gnomAD.AMR']=CSQ.get('gnomADg_AMR_AF',CSQ.get('gnomADe_AMR_AF','.'))
        self.fields['gnomAD.ASJ']=CSQ.get('gnomADg_ASJ_AF',CSQ.get('gnomADe_ASJ_AF','.'))
        self.fields['gnomAD.EAS']=CSQ.get('gnomADg_EAS_AF',CSQ.get('gnomADe_EAS_AF','.'))
        self.fields['gnomAD.FIN']=CSQ.get('gnomADg_FIN_AF',CSQ.get('gnomADe_FIN_AF','.'))
        self.fields['gnomAD.NFE']=CSQ.get('gnomADg_NFE_AF',CSQ.get('gnomADe_NFE_AF','.'))
        self.fields['gnomAD.OTH']=CSQ.get('gnomADg_OTH_AF',CSQ.get('gnomADe_OTH_AF','.'))
        self.fields['gnomAD.SAS']=CSQ.get('gnomADg_SAS_AF',CSQ.get('gnomADe_SAS_AF','.'))
        self.fields['gnomAD.MAX_AF']=CSQ.get('MAX_AF','.')
        self.fields['gnomAD.MAX_POPS']=CSQ.get('MAX_AF_POPS','.')

    def clinvar(self,CSQ):
        #ClinVar
        self.fields['ClinVar']=CSQ['ClinVar']#ClinVar
        self.fields['ClinVar.SIG']=CSQ['ClinVar_CLNSIG']#ClinVar.SIG
        self.fields['ClinVar.REVSTAT']=CSQ['ClinVar_CLNREVSTAT']#ClinVar.REVSTAT
        self.fields['ClinVar.DN']=CSQ['ClinVar_CLNDN']#ClinVar.CLNDN

    def loftee(self,CSQ):
        #LOFTEE
        self.fields['LOFTEE.lof']=CSQ.get('LoF','.')
        self.fields['LOFTEE.filter']=CSQ.get('LoF_filter','.')
        self.fields['LOFTEE.flags']=CSQ.get('LoF_flags','.')
        self.fields['LOFTEE.info']=CSQ.get('LoF_info','.')

    def alphamissense(self,CSQ):
        self.fields['AM.class']=CSQ.get('am_class','.')#AlphaMissense pathogenicity prediction
        self.fields['AM.pathogenicity']=CSQ.get('am_pathogenicity','.')#AlphaMissense pathogenicity score

    def mavedb(self,CSQ):
        self.fields['MaveDB.nt']=CSQ.get('MaveDB_nt','.')
        self.fields['MaveDB.pro']=CSQ.get('MaveDB_pro','.')
        self.fields['MaveDB.score']=CSQ.get('MaveDB_score','.')
        self.fields['MaveDB.urn']=CSQ.get('MaveDB_urn','.')
        ##MaveDB_nt=MaveDB HGVS (nucleotide); column from MaveDB_variants.tsv.gz
        ##MaveDB_pro=MaveDB HGVS (protein); column from MaveDB_variants.tsv.gz
        ##MaveDB_score=MaveDB score - see MaveDB for interpretation of scores; column from MaveDB_variants.tsv.gz
        ##MaveDB_urn=MaveDB database identifier; column from MaveDB_variants.tsv.gz

    def utrannotator(self):
        pass

    def lof_level(self):
        def check_level_one():
            def condition_one(cv,am):
                #ClinVar.SIG
                #AM.Class
                return ('pathogenic' in cv.lower() and 'conflicting' not in cv.lower()) or 'pathogenic' in am.lower()
            def condition_two(vc):
                #Variant.Consequence
                return 'frameshift' in vc.lower() or 'stop_gained' in vc.lower()
            def condition_three(ai_vals):
                return any(v>0.5 for v in ai_vals)
            def format_spliceai(vals):
                #SpliceAI.DS_AG,SpliceAI.DS_AL,SpliceAI.DS_DG,SpliceAI.DS_DL
                return [0.0 if x in ['.',''] else float(x) for x in vals]
            def condition_four(r):
                return r not in ['.',''] and float(r)>0.5
            return (condition_one(self.fields['ClinVar.SIG'],self.fields['AM.Class']) or condition_two(self.fields['Variant.Consequence']) or condition_three(format_spliceai([self.fields[v] for v in ['SpliceAI.DS_AG','SpliceAI.DS_AL','SpliceAI.DS_DG','SpliceAI.DS_DL']])) or condition_four(self.fields['REVEL']))

        def check_level_two():
            #ClinVar.SIG
            #AM.Class
            def condition_one(cv,am):
                return 'benign' not in cv.lower() and 'benign' not in am.lower()
            #Variant.Consequence
            def condition_two(vc):
                return any(x in vc.lower() for x in ['protein_altering','missense','inframe','start_lost'])
            #gnomadAD_AF_MAX
            def condition_three(ad):
                return ad in ['.',''] or float(ad)<0.01
            return (condition_one(self.fields['ClinVar.SIG'],self.fields['AM.Class']) and condition_two(self.fields['Variant.Consequence']) and condition_three(self.fields['gnomAD.MAX_AF']))

        if check_level_one():
            self.fields['Variant.LoF_level']='1'
        elif check_level_two():
            self.fields['Variant.LoF_level']='2'
        else:
            self.fields['Variant.LoF_level']='3'

    def fill_values(self,header):
        for h in header:
            if self.fields[h]=='':
                self.fields[h]='.'

    def print(self,header=['Chr','Start','REF','ALT','FILTER']):
        return ','.join([self.fields[h] for h in header])

    def report(self,writer):
        if len(self.calls)>0:
            for o in self.calls:
                writer.writerow({**self.fields,**o})
                #reverse dict order had incorrect fields.
        else:
            writer.writerow({**self.fields})

#END CLASS

def report_header(annotations,tumor_normal=False):
    if tumor_normal:
        header=['Tumor.ID','Normal.ID','Chr','Start','REF','ALT','FILTER']
    else:
        header=['Sample.ID','Chr','Start','REF','ALT','FILTER']
    header+=['Gene','Gene.Accession','Variant.LoF_level','Variant.Category','Variant.Class','Variant.Consequence']
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

    def loftee(tumor_normal=False):
        return ['LOFTEE.lof','LOFTEE.filter','LOFTEE.flags','LOFTEE.info']

    def alphamissense(tumor_normal=False):
        return ['AM.class','AM.pathogenicity']

    def mavedb(tumor_normal=False):
        return ['MaveDB.nt','MaveDB.pro','MaveDB.score','MaveDB.urn']

    def genotype(tumor_normal=False):
        if tumor_normal:
            return ['Tumor.Zyg','Tumor.Depth','Tumor.AltDepth','Tumor.AltFrac','Normal.Zyg','Normal.Depth','Normal.AltDepth','Normal.AltFrac']
        return ['Sample.Zyg','Sample.Depth','Sample.AltDepth','Sample.AltFrac']

    def vlr_values(tumor_normal=False):
        return ['PROB_SOMATIC_TUMOR','PROB_GERMLINE','PROB_SOMATIC_NORMAL','PROB_FFPE_ARTIFACT','PROB_ARTIFACT','PROB_ABSENT']

    #Cut OR rework
    if 'none' in annotations:
        return header
    elif 'everything' in annotations:
        for F in [splice_ai,snv_predition,gnomAD,loftee,clinvar,alphamissense,mavedb,genotype]:#,gnomAD, removed temp from before clinvar, the symbols are changed. loftee, too before genotype
            header+=F(tumor_normal)
    else:
        for x in ['splice_ai','snv_predition','gnomAD','loftee','clinvar','alphamissense','mavedb','genotype']:#here too'gnomAD''loftee',
            if x in annotations:
                header+=eval(f"{x}({tumor_normal})")
    if tumor_normal:
        header+=['LANCET','MUTECT2','STRELKA2','VARDICT','VARSCAN2']
    return header

def annotation_check(annotations):
    choices=['everything','snv_prediction','gnomAD','splice_ai','clinvar','alphamissense','mavedb','genotype','none']
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

def phred_to_probability(phred_score):
    #phred_score is a list
    if phred_score==["NA"] or phred_score==[]:
        return "NA"
    elif phred_score[0]=='inf' or phred_score[0]==float('inf'):
        return "0.0"
    return f"{10**(-float(phred_score[0])/10):.3e}"

def get_args(argv):
    #make single tumor_normal cohort as subparser?
    #scenrios may need some sort of sample id mapping file.
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input_vcf',help='Input: vcf file from ensembl-vep.')
    p.add_argument('-o','--output_csv',help='Output: csv file trimmed for specific design.')
    p.add_argument('-g','--gene_list',help='Optional list of genes to include only.')
    p.add_argument('-r','--bed_region',help='Bed file format to subset regions.')
    p.add_argument('-V','--include_vlr',action='store_true',help='Include transformed VLR probability values.')
    p.add_argument('-m','--mode',default='cohort',help='Run mode determines how calls are reported. Single should be "single,{sample.id}. Tumor/Normal should be "tumor_normal,{tumor.id},{normal.id}" OR use "no_sample"')
    p.add_argument('annotations',nargs=argparse.REMAINDER,default='everything',choices=['everything','snv_prediction','gnomAD','splice_ai','clinvar','alphamissense','mavedb','genotype','none'],help='annotation blocks to include')
    return p.parse_args(argv)
    #remove the whole everything choice. No one else is annotating, and we always do everything and subset columns after.

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
    if args.mode=='no_sample':
        for x in ['Sample.ID','Sample.Zyg','Sample.Depth','Sample.AltDepth','Sample.AltFrac']:
            header.remove(x)

    #I did it like this simply to control the column order
    if args.include_vlr:
        for x in ['PROB_SOMATIC_TUMOR','PROB_GERMLINE','PROB_SOMATIC_NORMAL','PROB_FFPE_ARTIFACT','PROB_ARTIFACT','PROB_ABSENT']:
            header.append(x)
    #test data
    #VcfReader=vcfpy.Reader.from_path('data/vcf/FLCN/PMBB-Release-2020-2.0_genetic_exome_FLCN_NF.norm.vep.vcf.gz')
    VcfReader=vcfpy.Reader.from_path(f'{args.input_vcf}')
    #VLR wants ANN tag instead of CSQ.
    #Changed CSQ to ANN
    csq_keys=VcfReader.header.get_info_field_info('ANN').serialize().split(' ')[-1].rstrip('">').split('|')
    with open(args.output_csv,'w') as outfile:
        writer=csv.DictWriter(outfile,fieldnames=header,delimiter=',',restval='.',extrasaction='ignore',quoting=csv.QUOTE_NONNUMERIC,dialect='excel')
        writer.writeheader()
        for record in VcfReader:
            if len(record.ALT)>1:
                print(f"Warning! : record.ALT length is {len(record.ALT)}. Not currently supported")
            vep_data=VEPannotation(record,tumor_normal,tumor)

            if 'ANN' not in record.INFO:
                continue
            if 'RefCall' in record.FILTER:
                continue

            if tumor_normal:
                tools=['LANCET','MUTECT2','STRELKA2','VARDICT','VARSCAN2']
                categories=[y for x in tools for y in record.INFO.get("CATEGORY",['NA']) if x in y and record.INFO.get(x,['NA'])==['PASS']]#Im a monster
                #print(categories)
                for x in tools:
                    vep_data.fields[x]=record.INFO.get(x,['NA'])[0]
                vep_data.fields["Variant.Category"]=";".join(categories) if categories else 'NA'


            if args.include_vlr:
                for x in VcfReader.header.info_ids():
                    if x.startswith("PROB_"):
                        vep_data.fields[x]=phred_to_probability(record.INFO.get(x,['NA']))
                        #Maybe this goes into class

            for csq_i in record.INFO['ANN']:
                csq_dict=dict(zip(csq_keys,csq_i.split('|')))
                if csq_dict['SYMBOL']!='' and csq_dict['CANONICAL']=='YES':
                    #checks
                    vep_data.info(csq_dict)
                    vep_data.snv_prediction(csq_dict)
                    vep_data.splice_ai(csq_dict)
                    vep_data.gnomAD(csq_dict)
                    vep_data.loftee(csq_dict)
                    vep_data.clinvar(csq_dict)
                    vep_data.alphamissense(csq_dict)
                    vep_data.mavedb(csq_dict)
                    vep_data.lof_level()

                    #more checks
                    if tumor_normal:
                        vep_data.calls[0]['Tumor.ID']=tumor
                        vep_data.calls[0]['Normal.ID']=normal

                    if single:
                        try:
                            vep_data.calls[0]['Sample.ID']=sample
                        except IndexError as e:
                            print(vep_data.print(),vep_data.calls)
                            print(record.INFO['ANN'])
                            print(record.calls)
                            continue

                    vep_data.fill_values(header)
                    if gene_filter and vep_data.fields['Gene'] in gene_list:
                        vep_data.report(writer)
                    elif region_filter and vep_data.in_region(bed_regions[record.CHROM]):
                        vep_data.report(writer)
                    elif not gene_filter and not region_filter:
                        vep_data.report(writer)
                    else:
                        continue
    print(f"{outfile.name} written.")

if __name__=='__main__':
    main()

#try:
#        snakemake
#    except NameError:
#        main(parse_arguments())
#    else:
#        main(parse_snakemake())

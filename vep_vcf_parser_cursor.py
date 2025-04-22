import vcfpy
import os
import csv
import argparse
from collections import defaultdict
import logging
import re

class SampleAnnot:
    @staticmethod
    def AD_check(AD):
        for i,v in enumerate(AD):
            if v in [None,'None','.']:
                AD[i]=0
        return AD

    @staticmethod
    def ZYG_check(GT,AD):
        if GT in ['.','./.']:
            #No GT info, try to make call from AD
            if not AD or sum(AD)==0:
                return '.'
            f=AD[1]/sum(AD)
            if f>=0.85:
                return 'HOM_ALT'
            elif f<=0.15:
                return 'HOM_REF'
            else:
                return 'HET_ALT'
        elif '.' in GT:
            #Partial GT info, use AD to refine
            #May need to add a haploid output.
            if not AD or sum(AD)==0:
                return GT
            f=AD[1]/sum(AD)
            if f>=0.85:
                return 'HOM_ALT'
            elif f<=0.15:
                return 'HOM_REF'
            else:
                return 'HET_ALT'
        else:
            #Full GT info
            return {'0/0':'HOM_REF','0/1':'HET_ALT','1/0':'HET_ALT','1/1':'HOM_ALT'}.get(GT,GT)

    @staticmethod
    def __get_calls__(record_calls):
        _calls=[]
        errors=defaultdict(int)

        def safe_process_call(c):
            x=defaultdict(str)
            x['Sample.ID']=c.sample
            x['Sample.Depth']=f"{c.data.get('DP',0)}"

            try:
                ad=SampleAnnot.AD_check(c.data.get('AD',[0,0]))
                x['Sample.AltDepth']=f"{ad[1]}"
                dp=c.data.get('DP',0)
                x['Sample.AltFrac']=f"{ad[1]/dp:.3f}" if dp else '.'
                x['Sample.Zyg']=f"{SampleAnnot.ZYG_check(c.data.get('GT','./.'),ad)}"

                return x if x['Sample.Zyg'] not in ['.','HOM_REF'] else None

            except (IndexError,ZeroDivisionError,TypeError) as e:
                errors[type(e).__name__]+=1
                return None

        _calls=[x for c in record_calls if (x:=safe_process_call(c))]

        if errors:
            logging.debug(f"Encountered errors processing calls: {dict(errors)}")

        return _calls

class TumorNormalAnnot:
    @staticmethod
    def __tumor_normal__(record_calls,tumor_id):
        x=defaultdict(str)
        for c in record_calls:
            #Set prefix based on whether this is tumor or normal
            prefix='Tumor.' if c.sample.lower()=='tumor' or c.sample==tumor_id else 'Normal.'

            #Basic sample information
            x[f'{prefix}ID']=c.sample
            x[f'{prefix}Depth']=f"{c.data.get('DP','.')}"
            x[f'{prefix}Zyg']=f"{c.data.get('GT','./.').replace('/',';')}"

            #First try to use AD if available
            if 'AD' in c.data and c.data.get('AD') not in [None,'.',['.'],[],['None']]:
                ad=c.data.get('AD')
                if isinstance(ad,list) and len(ad)>1:
                    x[f'{prefix}AltDepth']=f"{ad[1]}"
                else:
                    x[f'{prefix}AltDepth']=f"{ad}"

            #Then try SAOBS/SROBS from VLR
            elif 'SAOBS' in c.data and 'SROBS' in c.data:
                try:
                    #Extract counts from SAOBS (alt allele) and SROBS (ref allele)
                    saobs=c.data.get('SAOBS')
                    srobs=c.data.get('SROBS')
                    if isinstance(saobs,list) and saobs: saobs=saobs[0]
                    if isinstance(srobs,list) and srobs: srobs=srobs[0]

                    #Sum all numeric prefixes in SAOBS and SROBS
                    alt_count=0
                    ref_count=0
                    for match in re.finditer(r'(\d+)[A-Za-z]',saobs):
                        alt_count+=int(match.group(1))
                    for match in re.finditer(r'(\d+)[A-Za-z]',srobs):
                        ref_count+=int(match.group(1))

                    x[f'{prefix}AltDepth']=f"{alt_count}"

                    #Calculate AltFrac from SAOBS/SROBS counts
                    total_depth=alt_count+ref_count
                    if total_depth>0:
                        x[f'{prefix}AltFrac']=f"{alt_count/total_depth:.3f}"
                except (AttributeError,IndexError,ValueError,TypeError):
                    x[f'{prefix}AltDepth']='.'
            else:
                x[f'{prefix}AltDepth']='.'

            #Handle allele fraction if not already set from SAOBS/SROBS
            if f'{prefix}AltFrac' not in x or x[f'{prefix}AltFrac']=='':
                try:
                    if type(c.data.get('AF',['.']))==list:
                        af_value=c.data.get('AF',['.'])[0]
                        x[f'{prefix}AltFrac']=f"{float(af_value):.3f}" if af_value!='.' else '.'
                    else:
                        af_value=c.data.get('AF','.')
                        x[f'{prefix}AltFrac']=f"{float(af_value):.3f}" if af_value!='.' else '.'
                except (ValueError,TypeError,IndexError):
                    x[f'{prefix}AltFrac']='.'

            #Calculate AltFrac from AltDepth/Depth if not already set
            if (x[f'{prefix}AltFrac']=='.' or x[f'{prefix}AltFrac']=='0.000') and x[f'{prefix}AltDepth']!='.' and x[f'{prefix}Depth']!='.':
                try:
                    alt_depth=float(x[f'{prefix}AltDepth'])
                    depth=float(x[f'{prefix}Depth'])
                    if depth>0 and alt_depth>0:
                        x[f'{prefix}AltFrac']=f"{alt_depth/depth:.3f}"
                except (ValueError,TypeError):
                    pass

            #Calculate AltDepth from AF*DP if needed
            if x[f'{prefix}AltDepth']=='.' and x[f'{prefix}AltFrac']!='.' and x[f'{prefix}Depth']!='.':
                try:
                    x[f'{prefix}AltDepth']=f"{int(float(x[f'{prefix}Depth'])*float(x[f'{prefix}AltFrac']))}"
                except (ValueError,TypeError):
                    pass

            #Handle zygosity
            #Seems weird and overly complicated...
            if x[f'{prefix}Zyg']=='.;.':
                try:
                    if x[f'{prefix}AltFrac']!='.':
                        alt_frac=float(x[f'{prefix}AltFrac'])
                        if 0.0<alt_frac<1.0:
                            x[f'{prefix}Zyg']='HET'
                        elif alt_frac==1.0:
                            x[f'{prefix}Zyg']='HOM_ALT'
                        else:
                            x[f'{prefix}Zyg']='HOM_REF'
                    else:
                        x[f'{prefix}Zyg']='.'
                except (ValueError,TypeError):
                    pass
            else:
                x[f'{prefix}Zyg']={'0;0':'HOM_REF','0;1':'HET_ALT','1;0':'HET_ALT','1;1':'HOM_ALT'}.get(x[f'{prefix}Zyg'],x[f'{prefix}Zyg'])

        return [x]

class GnomadAnnot:
    def gnomAD(self,CSQ):
        populations=[
            ('AF',''),
            ('AFR','_AFR_AF'),
            ('AMR','_AMR_AF'),
            ('ASJ','_ASJ_AF'),
            ('EAS','_EAS_AF'),
            ('FIN','_FIN_AF'),
            ('NFE','_NFE_AF'),
            ('OTH','_OTH_AF'),
            ('SAS','_SAS_AF')
        ]
        for pop,suffix in populations:
            field_name=f'gnomAD.{pop}'
            g_value=CSQ.get(f'gnomADg{suffix}','.')
            e_value=CSQ.get(f'gnomADe{suffix}','.')
            self.fields[field_name]=g_value if g_value!='.' else e_value
        self.fields['gnomAD.MAX_AF']=CSQ.get('MAX_AF','.')
        self.fields['gnomAD.MAX_POPS']=CSQ.get('MAX_AF_POPS','.')

class ClinvarAnnot:
    def clinvar(self,CSQ):
        fields=[
            ('ClinVar','ClinVar'),
            ('ClinVar.SIG','ClinVar_CLNSIG'),
            ('ClinVar.REVSTAT','ClinVar_CLNREVSTAT'),
            ('ClinVar.DN','ClinVar_CLNDN'),
            ('AutoGVP','ClinVar_AutoGVP')
        ]
        for field,csq_key in fields:
            self.fields[field]=CSQ.get(csq_key,'.')

class AutoGVPAnnot:
    def autogvp(self,CSQ):
        fields=[('AutoGVP','ClinVar_AutoGVP')]
        for field,csq_key in fields:
            self.fields[field]=CSQ.get('AutoGVP','.')

class SpliceAIAnnot:
    def splice_ai(self,CSQ):
        metrics=['AG','AL','DG','DL']
        for m in metrics:
            self.fields[f'SpliceAI.DS_{m}']=CSQ[f'SpliceAI_pred_DS_{m}']

class PredictionAnnot:
    def snv_prediction(self,CSQ):
        predictors=[
            ('SIFT','SIFT'),
            ('PolyPhen','PolyPhen'),
            ('REVEL','REVEL')
        ]
        for field,csq_key in predictors:
            self.fields[field]=CSQ[csq_key]

class AlphaMissenseAnnot:
    def alphamissense(self,CSQ):
        fields=[
            ('AM.class','am_class'),
            ('AM.pathogenicity','am_pathogenicity')
        ]
        for field,csq_key in fields:
            self.fields[field]=CSQ.get(csq_key,'.')

class MaveDBAnnot:
    def mavedb(self,CSQ):
        fields=[
            ('MaveDB.nt','MaveDB_nt'),
            ('MaveDB.pro','MaveDB_pro'),
            ('MaveDB.score','MaveDB_score'),
            ('MaveDB.urn','MaveDB_urn')
        ]
        for field,csq_key in fields:
            self.fields[field]=CSQ.get(csq_key,'.')

class BasicInfoAnnot:
    def info(self,CSQ):
        fields=[
            ('Gene','SYMBOL'),
            ('Gene.Accession','Gene'),
            ('Variant.Class','VARIANT_CLASS'),
            ('Variant.Consequence','Consequence'),
            ('HGVSc','HGVSc'),
            ('HGVSp','HGVSp'),
            ('Feature.Type','Feature_type'),
            ('Feature.Accession','Feature'),
            ('Bio.type','BIOTYPE'),
            ('Existing.variation','Existing_variation'),
            ('EXON','EXON'),
            ('INTRON','INTRON'),
            ('STRAND','STRAND'),
            ('cDNA.position','cDNA_position'),
            ('CDS.position','CDS_position'),
            ('Protein.position','Protein_position'),
            ('Amino.acids','Amino_acids'),
            ('Codons','Codons')
        ]
        for field,csq_key in fields:
            value=CSQ[csq_key]
            if field in ['HGVSc','HGVSp']:
                value=value.split(':')[-1]
            elif field in ['EXON','INTRON']:
                value=value.replace('/','|')
            self.fields[field]=value
        self.fields['Variant.LoF_level']='.'

class LofLevelAnnot:
    def lof_level(self):
        def check_level_one():
            def condition_zero(bt):
                return 'protein_coding' in bt.lower()
            def condition_one(autogvp):#AutoGVP
                return ('pathogenic' in autogvp.lower())
            def condition_two(vc,en):
                if '|' not in en:
                    return False
                i,j=en.split('|')
                return ('frameshift' in vc or 'stop_gained' in vc) and i!=j
            def condition_three(af):
                if af in ['.','']:
                    return True
                return float(af)<=0.01
            def condition_four(vc,hgvs):
                p=re.compile(r'^c\.\d+([-+][12])([ACGT>]+|del|ins|dup)?$')
                if 'splice' in vc:
                    return bool(p.fullmatch(hgvs))
                return False

            return (condition_zero(self.fields['Bio.type']) and
                   (condition_one(self.fields['AutoGVP']) or
                    (condition_two(self.fields['Variant.Consequence'],self.fields['EXON']) and
                     condition_three(self.fields['gnomAD.MAX_AF'])) or
                    (condition_four(self.fields['Variant.Consequence'],self.fields['HGVSc']) and
                     condition_three(self.fields['gnomAD.MAX_AF']))))

        def check_level_two():
            def condition_zero(bt):
                return 'protein_coding' in bt.lower()
            def condition_one(autogvp):
                #May need to include AM not benign.
                return (autogvp.lower() in ['','.','uncertain_significance'])
            def condition_two(vc):
                return any(x in vc for x in ['protein_altering','missense','inframe','start_lost','frameshift','splice'])
            def condition_three(revel):
                if revel in ['.','']:
                    return True
                return float(revel)>0.5
            def condition_four(vc,ai):
                if 'splice' in vc:
                    return any(float(i)>0.5 for i in ai if i not in ['','.'])
                return True
            def condition_five(vc,hgvsp):
                if 'inframe' not in vc:
                    return True
                # Match Gln/Glu deletions/duplications
                p=re.compile(r'p\.[GQ]l[nu][0-9]+(?:_[GQ]l[nu][0-9]+)?(?:del|dup)')
                return not bool(p.match(hgvsp))

            return (condition_zero(self.fields['Bio.type']) and
                   condition_one(self.fields['AutoGVP']) and
                   condition_two(self.fields['Variant.Consequence']) and
                   condition_three(self.fields['REVEL']) and
                   condition_four(self.fields['Variant.Consequence'],[self.fields['SpliceAI.DS_AG'],self.fields['SpliceAI.DS_AL'],self.fields['SpliceAI.DS_DG'],self.fields['SpliceAI.DS_DL']]) and
                   condition_five(self.fields['Variant.Consequence'],self.fields['HGVSp']))

        def check_level_four():
            return ('benign' in self.fields['AutoGVP'].lower() or
                   'benign' in self.fields['ClinVar.SIG'].lower() or
                   'protein_coding' not in self.fields['Bio.type'].lower())

        if check_level_four():
            self.fields['Variant.LoF_level']='4'
        elif check_level_one():
            if any(x in self.fields['Variant.Consequence'] for x in ["stream","UTR","intron"]):
                self.fields['Variant.LoF_level']='2'
                #Special case for stream, UTR, intron that are "Pathogenic"
            else:
                self.fields['Variant.LoF_level']='1'
        elif check_level_two():
            self.fields['Variant.LoF_level']='2'
        else:
            self.fields['Variant.LoF_level']='3'
        #Need to know onco genes vs tumor suppressors from OncoKB, but this is a bad thing to have to run.
        #I can make a simple lookup table


class VEPannotation(BasicInfoAnnot,GnomadAnnot,ClinvarAnnot,SpliceAIAnnot,
                   PredictionAnnot,AlphaMissenseAnnot,MaveDBAnnot,
                   SampleAnnot,TumorNormalAnnot,LofLevelAnnot):
    """
    Main class for parsing VEP-annotated VCF files.

    Supports multiple modes:
    - Cohort analysis
    - Tumor/Normal paired analysis
    - Single sample analysis
    - No sample mode

    Attributes:
        fields (defaultdict): Stores variant annotations
        calls (list): Sample genotype information
        call_count (int): Number of calls processed
    """
    def __init__(self,vcfpy_record,tumor_normal=False,tumor_id=None):
        self.fields=defaultdict(str)
        self.fields['Chr']=f"{vcfpy_record.CHROM}"
        self.fields['Start']=f"{vcfpy_record.POS}"
        self.fields['REF']=f"{vcfpy_record.REF}"
        self.fields['ALT']=f"{vcfpy_record.ALT[0].value}"
        self.fields['FILTER']=f"{';'.join(vcfpy_record.FILTER)}".replace("MONOALLELIC",'.')
        self.fields['ID']=f"{vcfpy_record.ID[0]}" if len(vcfpy_record.ID)>0 else "."
        if tumor_normal:
            self.calls=self.__tumor_normal__(vcfpy_record.calls,tumor_id)
        else:
            self.calls=self.__get_calls__(vcfpy_record.calls)
        self.call_count=len(self.calls)

    def in_region(self,intervals):
        if any([start<=int(self.fields['Start'])<=end for start,end in intervals]):
            return True
        return False

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
        else:
            writer.writerow({**self.fields})

#END CLASS

def report_header(tumor_normal=False,no_caller=False):
    # Base headers
    if tumor_normal:
        header=['Tumor.ID','Normal.ID','Chr','Start','REF','ALT','FILTER','ID']
    else:
        header=['Sample.ID','Chr','Start','REF','ALT','FILTER','ID']

    # Core annotation headers
    header+=['Gene','Gene.Accession','Variant.LoF_level','Variant.Category','Variant.Class','Variant.Consequence',
            'HGVSc','HGVSp','Feature.Type','Feature.Accession','Bio.type','Existing.variation',
            'EXON','INTRON','STRAND','cDNA.position','CDS.position','Protein.position','Amino.acids','Codons']

    # All additional annotations
    header+=['SpliceAI.DS_AG','SpliceAI.DS_AL','SpliceAI.DS_DG','SpliceAI.DS_DL',  # SpliceAI
            'REVEL',  # SNV prediction
            'gnomAD.AF','gnomAD.AFR','gnomAD.AMR','gnomAD.ASJ','gnomAD.EAS','gnomAD.FIN','gnomAD.NFE','gnomAD.OTH','gnomAD.SAS','gnomAD.MAX_AF','gnomAD.MAX_POPS',  # gnomAD
            'ClinVar','ClinVar.SIG','ClinVar.REVSTAT','ClinVar.DN',  # ClinVar
            'AutoGVP',  # AutoGVP
            'AM.class','AM.pathogenicity',  # AlphaMissense
            'MaveDB.nt','MaveDB.pro','MaveDB.score','MaveDB.urn']  # MaveDB
            #'LOFTEE.lof','LOFTEE.filter','LOFTEE.flags','LOFTEE.info',  # LOFTEE doesnt work!

    # Genotype information
    if tumor_normal:
        header+=['Tumor.Zyg','Tumor.Depth','Tumor.AltDepth','Tumor.AltFrac',
                'Normal.Zyg','Normal.Depth','Normal.AltDepth','Normal.AltFrac']
    else:
        header+=['Sample.Zyg','Sample.Depth','Sample.AltDepth','Sample.AltFrac']

    # Tumor/Normal specific headers - only include if no_caller is False
    if tumor_normal and not no_caller:
        header+=['LANCET','MUTECT2','STRELKA2','VARDICT','VARSCAN2']

    return header


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
    return f"{10**(-float(phred_score[0])/10):.5f}"

def get_args(argv):
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input_vcf',help='Input: vcf file from ensembl-vep.')
    p.add_argument('-o','--output_csv',help='Output: csv file trimmed for specific design.')
    p.add_argument('-g','--gene_list',help='Optional list of genes to include only.')
    p.add_argument('-r','--bed_region',help='Bed file format to subset regions.')
    p.add_argument('-N','--no-caller',action='store_true',help='Do not include caller information.')
    p.add_argument('-R','--include-ref',action='store_true',help='Include RefCall in output.')
    p.add_argument('-V','--include_vlr',action='store_true',help='Include transformed VLR probability values.')
    p.add_argument('-m','--mode',default='cohort',help='Run mode determines how calls are reported. Single should be "single,{sample.id}. Tumor/Normal should be "tumor_normal,{tumor.id},{normal.id}" OR use "no_sample"')
    return p.parse_args(argv)

def main(argv=None):
    args=get_args(argv)
    gene_filter,gene_list=gene_list_check(args.gene_list)
    region_filter,bed_regions=bed_region_check(args.bed_region)
    tumor=None

    # Mode handling
    if args.mode.startswith('tumor_normal'):
        try:
            mode,tumor,normal=args.mode.split(',')
            tumor_normal=True
            single=False
            cohort=False
        except ValueError:
            raise ValueError("Tumor/Normal mode requires format: tumor_normal,{tumor.id},{normal.id}")
    elif args.mode.startswith('single'):
        try:
            mode,sample=args.mode.split(',')
            single=True
            tumor_normal=False
            cohort=False
        except ValueError:
            raise ValueError("Single mode requires format: single,{sample.id}")
    elif args.mode=='no_sample':
        single=False
        tumor_normal=False
        cohort=False
    elif args.mode=='cohort':  #Explicit cohort mode handling
        cohort=True
        single=False
        tumor_normal=False
        tumor=None
        sample=None
    else:
        raise ValueError(f"Unsupported mode: {args.mode}")

    header=report_header(tumor_normal,args.no_caller)

    if args.mode=='no_sample':
        for x in ['Sample.ID','Sample.Zyg','Sample.Depth','Sample.AltDepth','Sample.AltFrac']:
            header.remove(x)

    #I did it like this simply to control the column order
    if args.include_vlr:
        for x in ['PROB_SOMATIC_TUMOR','PROB_GERMLINE','PROB_SOMATIC_NORMAL','PROB_FFPE_ARTIFACT','PROB_ARTIFACT','PROB_ABSENT']:
            header.append(x)
    
    if args.include_ref:
        ref_call=True
    else:
        ref_call=False

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
            if 'RefCall' in record.FILTER and not ref_call:
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
                    vep_data.clinvar(csq_dict)
                    #vep_data.autogvp(csq_dict)
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

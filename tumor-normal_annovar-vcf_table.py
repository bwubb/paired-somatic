

import argparse
import vcfpy
import csv
import os
import sys
import datetime
from collections import defaultdict

class TranscriptError(Exception):
    def __init__(self,*args,**kwargs):
        Exception.__init__(self,*args,**kwargs)
#Sometimes I get non-exiting complaint about not having ##FILTER <ID=PASS...>

#GeneDetail.refGene "," sep, then ":" sep

def __zyg__(call):
    if not call.called:
        return '.'  # not called
    elif all(a == 0 for a in call.gt_alleles):
        return 'HOM_REF'
    elif len(set(call.gt_alleles)) == 1:
        return 'HOM_ALT'
    elif len(set(call.gt_alleles)) == 2:
        return 'HET_ALT'
    elif len(set(call.gt_alleles)) > 2:
        return 'COMPLEX'
    else:
        return 'HET'

def __AltFrac__(call):
    AF=[]
    if call.data.get('AD',False):
        #print(call)
        try:
            DP=float(call.data.get('DP'))
        except TypeError:#DP=. is NoneType
            DP=float(sum(call.data.get('AD')))
        for AD in call.data.get('AD')[1:]:#sometimes AD is a list...
            try:
                if AD>0:
                    AF.append('{0:.3f}'.format(float(AD)/DP))
            except (ZeroDivisionError,TypeError,ValueError) as e:
                print('{0}:'.format(e),'DP =',str(DP),'AD =',str(AD))
    if len(AF)==0:
        AF=['0.0']
    return ','.join(AF)

def __exonic__(annotation,gene_key,tx_id):#no multi gene in this
    def exonic_update(d,v):
        #print(v)
        try:#Nonframeshift sub have no aa, <INSERT> custom fix
            d['AAChange.refGene'].append(v[4])
        except IndexError:
            d['AAChange.refGene'].append('.')
        try:
            d['NTChange.refGene'].append(v[3])
            d['Exon.refGene'].append(v[2].replace('exon','Exon '))
        except IndexError:
            if v[2]=='wholegene':
                d['NTChange.refGene'].append('.')
                d['Exon.refGene'].append(v[2])
        d['Transcript.refGene'].append(v[1])
        return d
    pref_details=defaultdict(list)
    details=defaultdict(list)
    details['Gene.refGene']=[gene_key]
    for A in annotation.get('AAChange.refGene'):
        details['GeneDetail.refGene'].append(A.replace('%s:'%gene_key,'.'))
        v=A.split(':')
        if tx_id in v:
            pref_details.update(exonic_update(pref_details,v))
            pref_details['Gene.refGene']=[gene_key]
        else:
            details.update(exonic_update(details,v))
    if pref_details['Gene.refGene']:
        return dict((k,','.join(v)) for k,v in pref_details.items())
    else:
        return dict((k,','.join(v)) for k,v in details.items())

def __splicing__(annotation,gene_key,tx_id):
    def splicing_update(d,v):
        d['NTChange.refGene'].append(v[2])
        d['Exon.refGene'].append(v[1].replace('exon','Exon '))
        d['Transcript.refGene'].append(v[0])
        return d
    pref_details=defaultdict(list)
    details=defaultdict(list)
    details['Gene.refGene']=[gene_key]
    for A in annotation.get('GeneDetail.refGene'):
        details['GeneDetail.refGene'].append(A.replace('%s:'%gene_key,'.'))
        v=A.split(':')
        if tx_id in v:
            pref_details.update(splicing_update(pref_details,v))
            pref_details['Gene.refGene']=[gene_key]
        else:
            details.update(splicing_update(details,v))
    if pref_details['Gene.refGene']:
        return dict((k,','.join(v)) for k,v in pref_details.items())
    else:
        return dict((k,','.join(v)) for k,v in details.items())

def __exonicsplicing__(annotation,gene_key,tx_id):
    pref_details=defaultdict(list)
    try:
        exonic=__exonic__(annotation,gene_key,tx_id)#not a defaultdict, by design
    except IndexError:
        exonic={}
    splicing=__splicing__(annotation,gene_key,tx_id)#I mean kinda, it allows me to do try except below
    for d in [exonic,splicing]:
        if d.get('Transcript.refGene','.')==tx_id:
            pref_details=d.copy()
    details=exonic.copy()
    for k,v in splicing.items():
        try:
            details[k]=','.join([exonic[k],splicing[k]])
        except KeyError:
            details[k]=splicing[k]
    if pref_details['Gene.refGene']:
        return pref_details
    else:
        return details

def __UTR__(annotation,gene_key,tx_id):
    def UTR_update(d,v):
        d['NTChange.refGene'].append(v[1])
        d['Transcript.refGene'].append(v[0])
        return d
    pref_details=defaultdict(list)
    details=defaultdict(list)
    details['Gene.refGene']=[gene_key]
    for A in annotation.get('GeneDetail.refGene'):
        details['GeneDetail.refGene'].append(A.replace('%s:'%gene_key,'.'))
        v=A.split(':')
        if tx_id in v:
            pref_details.update(UTR_update(pref_details,v))
            pref_details['Gene.refGene']=[gene_key]
        else:
            details.update(UTR_update(details,v))
    if pref_details['Gene.refGene']:
        return dict((k,','.join(v)) for k,v in pref_details.items())
    else:
        return dict((k,','.join(v)) for k,v in details.items())


#GeneDetail.refGene=NM_001008536:c.*93A>C

def __intergenic__():
    pass

def serialize(ALT):
    return ','.join([alt.value for alt in ALT])

def __preferredTX__(tx_file):
    pTX={}
    with open(tx_file,'r') as file:
        for line in file:
            try:
                k,v=line.rstrip().split('\t')
                pTX[k]=v
            except IndexError:
                pass
    return pTX

def trimALT(ref,call):
    alts=[]
    print(ref,call)
    for a in call.gt_bases:
        if a!=ref and a not in alts:
            alts.append(a)
    return ','.join(alts)

def report_sample_variant(writer,fields):
    #filter checks
    if fields.get('SampleID',False):
        writer.writerow(fields)

def BasicCallFilter(call):
    #if call.data.get('GT','./.') in ['0/0']:#removed ./. until strelka is fixed
    #    return True
    if not call.data.get('DP',False):
        print('No Depth',call)
        return True
    elif set(call.data.get('AD',[0]))==set([0]):
        print('No Allele Depth',call)
        return True
    else:
        return False

def somatic_out_fielder(fields,record,header,argv):#need to change for tumor normal
    iT=next(i for i,x in enumerate(record.calls) if x.sample==argv.tumor_id or x.sample=='TUMOR')
    iN=next(i for i,x in enumerate(record.calls) if x.sample==argv.normal_id or x.sample=='NORMAL')
    
    fields['TumorID']=argv.tumor_id
    fields['NormalID']=argv.normal_id
    fields['Chr']=record.CHROM#Could set earlier
    fields['Start']=str(record.POS)#Could set earlier
    fields['Ref']=record.REF
    fields['Alt']=serialize(record.ALT)
    fields['FILTER']=','.join(record.FILTER)
    
    fields['Tumor_Genotype']='[{0}]'.format(record.calls[iT].data.get('GT').replace('/',','))
    
    #Some tumor records from GGA are DP=.
    #I need to look for those, try to use AD otherwise set some defaults.
    if record.calls[iT].data.get('DP',None)==None:
        fields['Tumor_Total_Depth']=str(sum(record.calls[iT].data.get('AD',0)))
    else:
        try:
            fields['Tumor_Total_Depth']=str(record.calls[iT].data.get('DP',0))#if I use .get() w/o default it does not raise KeyError; returns None
        except KeyError:
            fields['Tumor_Total_Depth']=str(sum(record.calls[iT].data.get('AD',0)))#Need defaults and checks for other FORMATs
    
    try:
        fields['Tumor_ALT_AlleleDepth']=','.join(str(d) for d in record.calls[iT].data.get('AD')[1:])#Could also be AO
    except TypeError:#quick fix
        fields['Tumor_ALT_AlleleDepth']='0.0'
    if fields['Tumor_ALT_AlleleDepth']=='.':
        fields['Tumor_ALT_AlleleDepth']='0.0'
    if record.calls[iT].data.get('AAF',False):
        fields['Tumor_ALT_AlleleFrac']=record.calls[iT].data['AAF']
    else:
        fields['Tumor_ALT_AlleleFrac']=__AltFrac__(record.calls[iT])
    fields['Tumor_Zyg']=__zyg__(record.calls[iT])
    
    fields['Normal_Genotype']='[{0}]'.format(record.calls[iN].data.get('GT').replace('/',','))
    
    if record.calls[iN].data.get('DP',None)==None:
        fields['Normal_Total_Depth']=str(sum(record.calls[iN].data.get('AD',0)))
    else:
        try:#pretty sure this is not needed
            fields['Normal_Total_Depth']=str(record.calls[iN].data.get('DP',0))
        except KeyError:
            fields['Normal_Total_Depth']=str(sum(record.calls[iN].data.get('AD',0)))
    try:
        fields['Normal_ALT_AlleleDepth']=','.join(str(d) for d in record.calls[iN].data.get('AD')[1:])#Could also be AO
    except TypeError:
        fields['Normal_ALT_AlleleDepth']='0.0'
    if record.calls[iN].data.get('AAF',False):
        fields['Normal_ALT_AlleleFrac']=record.calls[iN].data['AAF']
    else:
        fields['Normal_ALT_AlleleFrac']=__AltFrac__(record.calls[iN])
    fields['Normal_Zyg']=__zyg__(record.calls[iN])
    
    for k,v in record.INFO.items():
        if k in header and k not in fields.keys():
            try:#nci60 behaves weird, its either None or float, but not in a list as normal
                fields[k]=','.join(v).replace('\\x3d','=').replace('\\x3b',';')
            except TypeError:
                if v==None:
                    fields[k]='.'
                else:
                    fields[k]=str(v).replace('\\x3d','=').replace('\\x3b',';')
    if 'Somatic_LOH_pval' in header:
        fields['Somatic_LOH_pval']=str(record.INFO.get('SPV','.'))
    if 'Germline_pval' in header:
        fields['Germline_pval']=str(record.INFO.get('GPV','.'))
    #Quick fix
    for k,v in fields.items():
        if v=='':
            fields[k]='.'
    return fields

def main(argv=None):
    COUNT=defaultdict(int)
    with open(argv.annovar_header,'r') as file:
        summary_header=file.read().splitlines()
        #print(summary_header)
    reader=vcfpy.Reader.from_path(os.path.abspath(argv.input_fp))
    pTX=__preferredTX__(argv.preferred_transcripts)
    
    #if argv.output_fp:
    #   outfile=open('{0}.somesortofname.tsv'.format(argv.output_fp),'w')
    #else:
    #outfile=sys.stdout
    
    output=str(argv.input_fp).replace('.vcf','.report.tsv')#Fucking switch to bgzip tbi already
    outfile=open(output,'w')
    
    writer=csv.DictWriter(outfile,delimiter='\t',restval='.',extrasaction='ignore',fieldnames=summary_header)
    writer.writeheader()
    #Need to calculate my own END
    #if END use that
    #
    #Dup
    #else end=start+len(ref)-1
    #($newref, $newalt) = ($ref, '-')
    
    #Runner
    for record in reader:
        try:
            if '\\x3b' in record.INFO.get('Gene.refGene',['.'])[0]:
                record.INFO['Gene.refGene']=list(set(record.INFO.get('Gene.refGene')[0].split('\\x3b')))
        except IndexError:
            pass
        try:
            if '\\x3b' in record.INFO.get('GeneDetail.refGene',['.'])[0]:#I dont fail to get this, instead is an empty list which gives a error for index 0.
                record.INFO['GeneDetail.refGene']=list(set(record.INFO.get('GeneDetail.refGene')[0].split('\\x3b')))
        except IndexError:
            pass
        
        #record.INFO['genomicSuperDups']=[x.replace('\\x3b',';').replace('\\x3d','=') for x in record.INFO.get('genomicSuperDups',[])]
        fields={}#one init
        
        #Some Filtering#
        if record.ALT=='*':
            print('Star Alt')
            continue
        elif any('UNKNOWN' in record.INFO.get(x,'.') for x in ['Gene.refGene','Func.refGene','AAChange.refGene']):
            print('UNKNOWN',[record.INFO.get(x,'.') for x in ['Gene.refGene','Func.refGene','AAChange.refGene']])
            continue
        
        gene_key=record.INFO.get('Gene.refGene')[0]
        tx_id=pTX.get(gene_key,'NONE!')
        if record.INFO.get('Func.refGene','.')==['exonic']:
            try:
                fields=__exonic__(record.INFO,gene_key,tx_id)
                fields['GenomicRegion.refGene']='exonic'
                writer.writerow(somatic_out_fielder(fields,record,summary_header,argv))
                COUNT['exonic']+=1
            except (KeyError,IndexError):
                print(record.CHROM,record.POS,record.INFO.get('Func.refGene'),record.INFO.get('Gene.refGene'),record.INFO.get('AAChange.refGene'),'is not a supported record')
        elif record.INFO.get('Func.refGene','.')==['exonic\\x3bsplicing']:
            #print('running exonic;splicing!!!!!!!!')
            #previous method, splits Gene_refGene by ; with 0 being used in exonic and 1 used in splicing
            #print(record.INFO.get('AAChange.refGene'),record.INFO.get('GeneDetail.refGene'))
            fields=__exonicsplicing__(record.INFO,gene_key,tx_id)
            fields['GenomicRegion.refGene']='exonic;splicing'
            writer.writerow(somatic_out_fielder(fields,record,summary_header,argv))
            COUNT['exonic;splicing']+=1
        elif record.INFO.get('Func.refGene','.')==['splicing']:
            fields=__splicing__(record.INFO,gene_key,tx_id)
            fields['GenomicRegion.refGene']='splicing'
            writer.writerow(somatic_out_fielder(fields,record,summary_header,argv))
            COUNT['splicing']+=1
        elif record.INFO.get('Func.refGene') in [['UTR3'],['UTR5']]:
            #print(record.INFO['GeneDetail.refGene'])
            fields=__UTR__(record.INFO,gene_key,tx_id)
            fields['GenomicRegion.refGene']=record.INFO.get('Func.refGene')[0]
            writer.writerow(somatic_out_fielder(fields,record,summary_header,argv))
            COUNT['UTR']+=1
        elif record.INFO.get('Func.refGene','.')==['intronic']:
            #would also work for ncRNA_exonic,ncRNA_intronic
            fields={}
            fields['Gene.refGene']=record.INFO.get('Gene.refGene')[0]
            fields['Transcript.refGene']=pTX.get(record.INFO.get('Gene.refGene')[0],'.')#intronic is no longer annovar'd well. No transcript given.
            if any(record.INFO.get(i)!=[] for i in ['dbscSNV_ADA_SCORE','dbscSNV_RF_SCORE']):#['value']
                fields['GenomicRegion.refGene']='intronic;splicing'
            else:
                fields['GenomicRegion.refGene']='intronic'
            writer.writerow(somatic_out_fielder(fields,record,summary_header,argv))
            COUNT[fields['GenomicRegion.refGene']]+=1
            #print(count,line)
        else:
            try:
                fields={}
                fields['Gene.refGene']=record.INFO.get('Gene.refGene')[0]
                fields['Transcript.refGene']=pTX.get(record.INFO.get('Gene.refGene')[0],'.')
                fields['GenomicRegion.refGene']=record.INFO.get('Func.refGene')[0].replace('\\x3b',';')
                writer.writerow(somatic_out_fielder(fields,record,summary_header,argv))
                COUNT[fields['GenomicRegion.refGene']]+=1
            except (KeyError,IndexError):
                print(record.CHROM,record.POS,record.REF,record.ALT,record.INFO.get('Func.refGene'),record.INFO.get('Gene.refGene'),'is not a supported record')
    outfile.close()
    for k,v in COUNT.items():
        print(k,v)

if __name__=='__main__':
    csv.field_size_limit(sys.maxsize)
    p=argparse.ArgumentParser()
    p.add_argument('-I','--input_fp',required=True,help='Annovar vcf.gz')
    p.add_argument('--tumor_id',required=True,help='Tumor id as it appears in input')
    p.add_argument('--normal_id',required=True,help='Normal id as it appears in input')
    p.add_argument('--preferred_transcripts',default='/home/bwubb/resources/preferred_transcripts.20170505.txt',help='File with preferred refGene transcript IDs')
    p.add_argument('--annovar_header',default='/home/bwubb/resources/annovar/somatic-annotation-header.20181218.txt',help='File with annovar columns to keep')
    argv=p.parse_args()
    print('Arguments Initialized...')
    #LOG the run conditions
    for k,v in vars(argv).items():
        print(k,':',v)
    main(argv)

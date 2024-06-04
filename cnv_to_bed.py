import argparse
import os
import csv
from collections import defaultdict

def sequenza_row(row,ucsc=False):
    try:
        start=int(row['start.pos'])-1
        end=int(row['end.pos'])
        type=[]
        if row['A']=='0' or row['B']=='0':
            type.append('loh')
        elif row['A']!=row['B']:
            type.append('imbalance')
        else:
            type.append('nonloh')
        if int(row['CNt'])>4:
            type.append('amp')
        elif 2<int(row['CNt'])<=4:
            type.append('gain')
        elif 1<=int(row['CNt'])<2:
            type.append('loss')
        elif int(row['CNt'])<1:
            type.append('del')
        elif int(row['CNt'])==2:
            type.append('neutral')
        else:
            type.append('unknown')
        name=f"{end-start+1}bp;{';'.join(type)};A{row['A']};B{row['B']}"
        score=f"{int(row['CNt'])*100}"
        strand="+"
        bed_row={'chrom':row['chromosome'],'chromStart':f"{start}",'chromEnd':f"{end}",'name':name,'score':score,'strand':strand}
    except ValueError:
        return None
    return bed_row

def cnvkit_row(row,ucsc=False):
    try:
        start=int(row['start'])-1
        end=int(row['end'])
        type=[]
        #check
        if int(row['cn1'])<0:
            print(f"Warning: cn1={row['cn1']}")
            row['cn1']='0'
        if int(row['cn2'])<0:
            print(f"Warning: cn2={row['cn2']}")
            row['cn2']='0'
        if int(row['cn'])<0:
            print(f"Warning: cn={row['cn']}")
            row['cn']='0'

        if row['cn1']=='0' or row['cn2']=='0':
            type.append('loh')
        elif row['cn1']!=row['cn2']:
            type.append('imbalance')
        else:
            type.append('nonloh')
        if int(row['cn'])>4:
            type.append('amp')
        elif 2<int(row['cn'])<=4:
            type.append('gain')
        elif 1<=int(row['cn'])<2:
            type.append('loss')
        elif int(row['cn'])<1:
            type.append('del')
        elif int(row['cn'])==2:
            type.append('neutral')
        else:
            type.append('unknown')
        name=f"{end-start+1}bp;{';'.join(type)};A{row['cn1']};B{row['cn2']}"
        score=f"{int(row['cn'])*100}"
        strand="+"
        bed_row={'chrom':row['chromosome'],'chromStart':f"{start}",'chromEnd':f"{end}",'name':name,'score':score,'strand':strand}
    except ValueError:
        return None
    return bed_row

#Sampleid,chr,start,end,arm,C,M
def purecn_row(row,ucsc=False):
    try:
        start=int(row['start'])-1
        end=int(row['end'])
        type=[]
        CN=int(row['C'])
        B=int(row['M'])
        A=CN-B
        if any([A==0,B==0]):
            type.append('loh')
        elif A!=B:
            type.append('imbalance')
        else:
            type.append('nonloh')
        if int(CN)>4:
            type.append('amp')
        elif 2<int(CN)<=4:
            type.append('gain')
        elif 1<=int(CN)<2:
            type.append('loss')
        elif int(CN)<1:
            type.append('del')
        elif int(CN)==2:
            type.append('neutral')
        else:
            type.append('unknown')
        name=f"{end-start+1}bp;{';'.join(type)};A{A};B{B}"
        score=f"{int(CN)*100}"
        strand="+"
        bed_row={'chrom':row['chr'],'chromStart':f"{start}",'chromEnd':f"{end}",'name':name,'score':score,'strand':strand}
    except ValueError:
        return None
    return bed_row

def ascat_row(row,ucsc=False):
    try:
        start=int(row['start'])-1
        end=int(row['end'])
        type=[]
        CN=f"{int(row['nMajor'])+int(row['nMinor'])}"
        if any([row['nMajor']=='0',row['nMinor']=='0']):
            type.append('loh')
        elif row['nMajor']!=row['nMinor']:
            type.append('imbalance')
        else:
            type.append('nonloh')
        if int(CN)>4:
            type.append('amp')
        elif 2<int(CN)<=4:
            type.append('gain')
        elif 1<=int(CN)<2:
            type.append('loss')
        elif int(CN)<1:
            type.append('del')
        elif int(CN)==2:
            type.append('neutral')
        else:
            type.append('unknown')
        name=f"{end-start+1}bp;{';'.join(type)};A{row['nMajor']};B{row['nMinor']}"
        score=f"{int(CN)*100}"
        strand="+"
        bed_row={'chrom':row['chromosome'],'chromStart':f"{start}",'chromEnd':f"{end}",'name':name,'score':score,'strand':strand}
    except ValueError:
        return None
    return bed_row

#There is a simpler segments.txt file?
def facets_row(row,ucsc=False):
    try:
        start=int(row['start'])-1
        end=int(row['end'])
        try:
            length=end-start+1
        except ValueError:
            length='.'
        type=[]
        if any([row['lcn.em']=='0',row['lcn.em']==row['tcn.em']]):
            type+=['loh']
        elif all([row['lcn.em'] not in ['NA','0'],row['lcn.em']!=row['tcn.em'],row['tcn.em']!='2']):
            type+=['imbalance']
        else:
            type+=['nonloh']
        if int(row['tcn.em'])>4:
            type+=['amp']
        elif int(row['tcn.em'])>2 and int(row['tcn.em'])<=4:
            type+=['gain']
        elif int(row['tcn.em'])<2 and int(row['tcn.em'])>=1:
            type+=['loss']
        elif int(row['tcn.em'])<1:
            type+=['del']
        elif int(row['tcn.em'])==2:
            type+=['neutral']
        else:
            type+=['unknown']

        if row["lcn.em"]=="NA":
            A="NA"
            B="NA"
        else:
            A=f"{int(row['tcn.em'])-int(row['lcn.em'])}"
            B=row["lcn.em"]
        name=f"{length}bp;{';'.join(type)};A{A};B{B}"
        score=f"{int(row['tcn.em'])*100}"
        strand="+"
        if row['chrom']=='23':
            chrom='X'
        else:
            chrom=row['chrom']
        bed_row={'chrom':chrom,'chromStart':f"{start}",'chromEnd':f"{end}",'name':f"{name}",'score':f"{score}",'strand':strand}
        #rint(bed_row)
    except ValueError:
        return None
    return bed_row

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('-c','--caller',choices=['sequenza','cnvkit','ascat','purecn','facets'],default='sequenza',help='Which ASCNV caller is this data from?')
    p.add_argument('--ucsc',action='store_true',default=False,help='Format outfiles for UCSC browser')
    p.add_argument('input_fp',nargs=argparse.REMAINDER,help='One or more input files')
    argv=p.parse_args()
    return vars(argv)

def main(argv=None):
    bed_fields=['chrom','chromStart','chromEnd','name','score','strand']
    argv=get_args() if argv is None else argv
    for i,file in enumerate(argv['input_fp']):
        base,ext=os.path.splitext(file)
        outfile=f"{base}.bed"
        with open(file,'r') as infile, open(outfile,'w',newline='') as bed_file:
            if ext=='.csv':
                reader=csv.DictReader(infile,delimiter=',')
            else:
                reader=csv.DictReader(infile,delimiter='\t')
            writer=csv.DictWriter(bed_file,delimiter='\t',fieldnames=bed_fields)
            #writer.writeheader()
            if argv['caller'] == 'cnvkit':
                run_row=cnvkit_row
            elif argv['caller'] == 'sequenza':
                run_row=sequenza_row
            elif argv['caller'] == 'ascat':
                run_row=ascat_row
            elif argv['caller']=='purecn':
                run_row=purecn_row
            elif argv['caller']=='facets':
                run_row=facets_row
            else:
                raise ValueError(f"Unsupported caller: {argv['caller']}")
            for row in reader:
                out_row=run_row(row,argv['ucsc'])
                if out_row is not None and out_row['chrom'] not in ['24','Y','chrY']:
                    writer.writerow(out_row)

if __name__=='__main__':
    main()

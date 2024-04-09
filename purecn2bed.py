import argparse
import csv
from collections import defaultdict

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('--ucsc',action='store_true',default=False,help='Format outfiles for UCSC browser')
    p.add_argument('input_fp',nargs=argparse.REMAINDER,help='One or more input files')
    argv=p.parse_args()
    #converts argparse.Namespace() to dict
    return vars(argv)

def format_row(row,ucsc=False):
    try:
        start=int(row['start'])-1
        end=int(row['end'])
        try:
            length=end-start+1
        except ValueError:
            length='.'
        type=[]
        if any([row['type']=='LOH']):
            type+=['loh']
        elif all([row['M'] not in ['NA','0'],int(row['M'])!=(int(row['C'])-int(row['M'])),row['C']!='2']):
            type+=['imbalance']
        else:
            type+=['nonloh']
        if int(row['C'])>4:
            type+=['amp']
        elif int(row['C'])>2 and int(row['C'])<=4:
            type+=['gain']
        elif int(row['C'])<2 and int(row['C'])>=1:
            type+=['loss']
        elif int(row['C'])<1:
            type+=['del']
        elif int(row['C'])==2:
            type+=['neutral']
        else:
            type+=['unknown']

        if row["M"]=="NA":
            A="NA"
            B="NA"
        else:
            A=f"{int(row['C'])-int(row['M'])}"
            B=row["M"]
        name=f"{length}bp;{';'.join(type)};A{A};B{B}"
        score=f"{int(row['C'])*100}"
        strand="+"
        if row['chr']=='23':
            chrom='X'
        else:
            chrom=row['chr']
        bed_row={'chrom':chrom,'chromStart':f"{start}",'chromEnd':f"{end}",'name':f"{name}",'score':f"{score}",'strand':strand}
        #rint(bed_row)
    except ValueError:
        return None
    return bed_row

def main(argv=None):
    bed_fields=['chrom','chromStart','chromEnd','name','score','strand']
    #d={'chromosome':'Chr','start.pos':'Start','end.pos':'End','A':'A','B':'B'}
    for i,file in enumerate(argv['input_fp']):
        outfile=file.replace('csv','bed')
        #print(outfile)
        with open(file,'r') as infile,open(outfile,'w') as bed_file:
            reader=csv.DictReader(infile,delimiter=',')
            writer=csv.DictWriter(bed_file,delimiter='\t',fieldnames=bed_fields,lineterminator='\n')
            for row in reader:
                out_row=format_row(row)
                if out_row!=None:
                    writer.writerow(out_row)

if __name__=='__main__':
    main(get_args())

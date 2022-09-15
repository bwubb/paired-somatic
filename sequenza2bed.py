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

def get_snake():
    ucsc=snakemake.params.get('ucsc','False') in ['True','true','1','yes','Yes']
    lib=snakemake.params['targets_key']
    argv={'lib':lib,'ucsc':ucsc,'input_fp':list(snakemake.input)}
    return argv

def format_row(row,ucsc=False):
    #print(row)
    try:
        start=int(row['start.pos'])-1
        end=int(row['end.pos'])
        try:
            length=end-start+1
        except ValueError:
            length='.'
        type=[]
        if any([row['A']=='0',row['B']=='0']):
            type+=['loh']
        elif row['A']!=row['B']:
            type+=['imbalance']
        else:
            type+=['nonloh']
        if int(row['CNt'])>4:
            type+=['amp']
        elif int(row['CNt'])>2 and int(row['CNt'])<=4:
            type+=['gain']
        elif int(row['CNt'])<2 and int(row['CNt'])>=1:
            type+=['loss']
        elif int(row['CNt'])<1:
            type+=['del']
        elif int(row['CNt'])==2:
            type+=['neutral']
        else:
            type+=['unknown']
        name=f"{length}bp;{';'.join(type)};A{row['A']};B{row['B']}"
        score=f"{int(row['CNt'])*100}"
        strand="+"
        bed_row={'chrom':row['chromosome'],'chromStart':f"{start}",'chromEnd':f"{end}",'name':f"{name}",'score':f"{score}",'strand':strand}
        #print(bed_row)
    except ValueError:
        return None
    return bed_row

def main(argv=None):
    bed_fields=['chrom','chromStart','chromEnd','name','score','strand']
    #d={'chromosome':'Chr','start.pos':'Start','end.pos':'End','A':'A','B':'B'}
    for i,file in enumerate(argv['input_fp']):
        outfile=file.replace('txt','bed')
        print(outfile)
        with open(file,'r') as infile,open(outfile,'w') as bed_file:
            reader=csv.DictReader(infile,delimiter='\t')
            writer=csv.DictWriter(bed_file,delimiter='\t',fieldnames=bed_fields,lineterminator='\n')
            for row in reader:
                out_row=format_row(row)
                if out_row!=None:
                    writer.writerow(out_row)

###Additional Feature, LPP filtering

if __name__=='__main__':
    try:
        snakemake
    except NameError:
        main(get_args())
    else:
        main(get_snake())

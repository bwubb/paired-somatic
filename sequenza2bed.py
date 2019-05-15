import os
import sys
import csv
from collections import defaultdict


try:
    snakemake
except NameError:
    input=sys.argv[1:]
    output=[i.replace('txt','bed') for i in input if i.endswith('txt')]
else:
    input=snakemake.input
    output=[i.replace('txt','bed') for i in input if i.endswith('txt')]

bed_fields=['chrom','chromStart','chromEnd','name','score','strand']
#d={'chromosome':'Chr','start.pos':'Start','end.pos':'End','A':'A','B':'B'}
for i,file in enumerate(input):
    with open(file,'r') as segments_file,open(output[i],'w') as bed_file:
        print(f'Reading {segments_file.name}...')
        reader=csv.DictReader(segments_file,delimiter='\t')
        writer=csv.DictWriter(bed_file,delimiter='\t',fieldnames=bed_fields,lineterminator='\n')
        for row in reader:
            start=int(row['start.pos'])-1
            end=int(row['end.pos'])
            try:
                length=end-start+1
            except ValueError:
                length='.'
            type=[]
            if any([row['A']=='0',row['B']=='0']):
                type+=['loh']
            if int(row['CNt'])>2:
                type+=['gain']
            elif int(row['CNt'])<2:
                type+=['loss']
            elif int(row['CNt'])==2:
                type+=['neutral']
            else:
                type+=['unknown']
            name=f"{length}bp;{';'.join(type)}_A{row['A']};B{row['B']}"
            score=f"{int(row['CNt'])*100}"
            strand="+"
            bed_row={'chrom':row['chromosome'],'chromStart':f"{start}",'chromEnd':f"{end}",'name':f"{name}",'score':f"{score}",'strand':strand}
            writer.writerow(bed_row)

###Additional Feature, LPP filtering
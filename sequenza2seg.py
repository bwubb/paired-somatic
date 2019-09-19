import os
import sys
import csv
import math

#need samples

sample=''

try:
    snakemake
except NameError:
    input=sys.argv[1:]
    output=[i.replace('txt','seg') for i in input if i.endswith('txt')]
else:
    input=snakemake.input
    output=[i.replace('txt','seg') for i in input if i.endswith('txt')]

sample='5093-Brca1Ov13'

seg_fields=['ID','chrom','loc.start','loc.end','num.mark','seg.mean']
#d={'chromosome':'Chr','start.pos':'Start','end.pos':'End','A':'A','B':'B'}
for i,file in enumerate(input):
    with open(file,'r') as segments_file,open(output[i],'w') as seg_file:
        print(f'Reading {segments_file.name}...')
        reader=csv.DictReader(segments_file,delimiter='\t')
        writer=csv.DictWriter(seg_file,delimiter='\t',fieldnames=seg_fields,lineterminator='\n')
        for row in reader:
            outrow={'ID':sample,'chrom':row['chromosome'],'loc.start':row['start.pos'],'loc.end':row['end.pos'],'num.mark':row['N.BAF'],'seg.mean':f"{math.log(int(row['CNt'])/2,2)}"}
            writer.writerow(outrow)
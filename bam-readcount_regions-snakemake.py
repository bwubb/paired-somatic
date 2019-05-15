
import csv
import sys
from itertools import dropwhile

def is_comment(s):
    return s.startswith('#')

try:
    snakemake
except NameError:
    input=sys.argv[1:]
    output=[i.replace('vcf','regions') for i in input if i.endswith('vcf')]
else:
    input=snakemake.input
    output=snakemake.output

for i,file in enumerate(input):
    with open(input[i],'r') as psuedo_vcf,open(output[i],'w') as outfile:
        writer=csv.writer(outfile,delimiter='\t')
        for line in dropwhile(is_comment,psuedo_vcf):
            if any([line.startswith('M'),line.startswith('GL')]):
                break#hard stop
            xline=line.split('\t')[:5]
            #snp and mnp
            if xline[3]!='-' and xline[4]!='-' and len(xline[3])==len(xline[4]):
                writer.writerow([xline[0],xline[1],str(int(xline[1])+len(xline[3])-1)])
            elif xline[4]=='-':#deletion
                writer.writerow([xline[0],xline[1],str(int(xline[1])+len(xline[3])-1)])
            else:#insertion
                writer.writerow([xline[0],xline[1],str(int(xline[1])+1)])

import csv
from itertools import dropwhile

def is_comment(s):
    return s.startswith('#')

with open(snakemake.input[0],'r') as psuedo_vcf,open(snakemake.output[0],'w') as outfile:
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
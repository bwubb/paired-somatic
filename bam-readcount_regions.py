import csv
import sys
from itertools import dropwhile

#USAGE: python bam-readcount_regions.py A.vcf [B.vcf ...]

def is_comment(s):
    return s.startswith('#')

def parse_arguments():
    argv={}
    argv['input']=sys.argv[1:]
    argv['output']=[i.replace('vcf','regions') for i in sys.argv[1:] if i.endswith('vcf')]
    return argv

def parse_snakemake():
    argv={}
    argv['input']=snakemake.input
    argv['output']=snakemake.output
    return argv

def main(argv=None):
    input=argv['input']
    output=argv['output']
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

if __name__=='__main__':
    try:
        snakemake
    except NameError:
        main(parse_arguments())
    else:
        main(parse_snakemake())
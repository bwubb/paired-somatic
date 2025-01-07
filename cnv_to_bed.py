"""
CNV to BED Format Converter

Converts Copy Number Variation (CNV) output from various callers into BED format.
Handles output from: Sequenza, CNVkit, ASCAT, PureCN, and FACETS.

BED Format Output:
    - chrom: Chromosome number
    - chromStart: 0-based start position
    - chromEnd: End position
    - name: Format is "{length}bp;{type};A{cn_a};B{cn_b}"
        - type includes: loh/nonloh/imbalance + amp/gain/loss/del/neutral
    - score: Total copy number * 100
    - strand: Always "+"

Usage:
    python cnv_to_bed.py -c [caller] [input_files]
    
Required Arguments:
    input_files: One or more input CNV files
    
Optional Arguments:
    -c/--caller: CNV caller (sequenza,cnvkit,ascat,purecn,facets)
    --ucsc: Format output for UCSC browser
"""

import argparse
import os
import csv
from collections import defaultdict

class CNVProcessor:
    def __init__(self,ucsc=False):
        self.ucsc=ucsc
    
    def get_copy_number_type(self,cn_total,cn_a,cn_b):
        type=[]
        if cn_a==0 or cn_b==0:
            type.append('loh')
        elif cn_a!=cn_b:
            type.append('imbalance')
        else:
            type.append('nonloh')
        if cn_total>4:
            type.append('amp')
        elif 2<cn_total<=4:
            type.append('gain')
        elif 1<=cn_total<2:
            type.append('loss')
        elif cn_total<1:
            type.append('del')
        elif cn_total==2:
            type.append('neutral')
        else:
            type.append('unknown')
        return type
    
    def format_bed_row(self,chrom,start,end,cn_total,cn_a,cn_b):
        try:
            start=int(start)-1
            end=int(end)
            type=self.get_copy_number_type(cn_total,cn_a,cn_b)
            name=f"{end-start+1}bp;{';'.join(type)};A{cn_a};B{cn_b}"
            score=f"{int(cn_total)*100}"
            return {
                'chrom':chrom,
                'chromStart':f"{start}",
                'chromEnd':f"{end}",
                'name':name,
                'score':score,
                'strand':"+"
            }
        except ValueError:
            return None

class SequenzaProcessor(CNVProcessor):
    def process_row(self,row):
        try:
            cn_total=int(row['CNt'])
            cn_a=int(row['A'])
            cn_b=int(row['B'])
            return self.format_bed_row(
                row['chromosome'],
                row['start.pos'],
                row['end.pos'],
                cn_total,
                cn_a,
                cn_b
            )
        except ValueError:
            return None

class CNVkitProcessor(CNVProcessor):
    def process_row(self,row):
        try:
            if int(row['cn1'])<0:
                print(f"Warning: cn1={row['cn1']}")
                row['cn1']='0'
            if int(row['cn2'])<0:
                print(f"Warning: cn2={row['cn2']}")
                row['cn2']='0'
            if int(row['cn'])<0:
                print(f"Warning: cn={row['cn']}")
                row['cn']='0'
            cn_total=int(row['cn'])
            cn_a=int(row['cn1'])
            cn_b=int(row['cn2'])
            return self.format_bed_row(
                row['chromosome'],
                row['start'],
                row['end'],
                cn_total,
                cn_a,
                cn_b
            )
        except ValueError:
            return None

class PureCNProcessor(CNVProcessor):
    def process_row(self,row):
        try:
            cn_total=int(row['C'])
            cn_b=int(row['M'])
            cn_a=cn_total-cn_b
            return self.format_bed_row(
                row['chr'],
                row['start'],
                row['end'],
                cn_total,
                cn_a,
                cn_b
            )
        except ValueError:
            return None

class ASCATProcessor(CNVProcessor):
    def process_row(self,row):
        try:
            cn_a=int(row['nMajor'])
            cn_b=int(row['nMinor'])
            cn_total=cn_a+cn_b
            return self.format_bed_row(
                row['chromosome'],
                row['start'],
                row['end'],
                cn_total,
                cn_a,
                cn_b
            )
        except ValueError:
            return None

class FacetsProcessor(CNVProcessor):
    def process_row(self,row):
        try:
            if row['lcn.em']=="NA":
                cn_a="NA"
                cn_b="NA"
            else:
                cn_total=int(row['tcn.em'])
                cn_b=int(row['lcn.em'])
                cn_a=cn_total-cn_b
            chrom='X' if row['chrom']=='23' else row['chrom']
            return self.format_bed_row(
                chrom,
                row['start'],
                row['end'],
                int(row['tcn.em']),
                cn_a,
                cn_b
            )
        except ValueError:
            return None

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
    
    processors={
        'sequenza':SequenzaProcessor,
        'cnvkit':CNVkitProcessor,
        'ascat':ASCATProcessor,
        'purecn':PureCNProcessor,
        'facets':FacetsProcessor
    }
    
    for i,file in enumerate(argv['input_fp']):
        base,ext=os.path.splitext(file)
        outfile=f"{base}.bed"
        with open(file,'r') as infile,open(outfile,'w',newline='') as bed_file:
            if ext=='.csv':
                reader=csv.DictReader(infile,delimiter=',')
            else:
                reader=csv.DictReader(infile,delimiter='\t')
            writer=csv.DictWriter(bed_file,delimiter='\t',fieldnames=bed_fields)
            
            processor=processors.get(argv['caller'])(ucsc=argv['ucsc'])
            if not processor:
                raise ValueError(f"Unsupported caller: {argv['caller']}")
            
            for row in reader:
                out_row=processor.process_row(row)
                if out_row is not None and out_row['chrom'] not in ['24','Y','chrY']:
                    writer.writerow(out_row)

if __name__=='__main__':
    main()

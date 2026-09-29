"""
CNV to BED Format Converter

Converts CNV caller output into BED for AnnotSV.
Callers: sequenza, cnvkit, ascat, purecn, facets, codex2

BED name (user#1): "{length}bp;{zyg};{type};A{cn_a};B{cn_b}"
  zyg: loh / imbalance / nonloh / unknown (unknown when no allele CN)
  type: amp / gain / loss / del / neutral
  A/B: allele copy numbers, or NA for total-CN-only callers (CODEX2)
BED score (user#2): total CN * 100

Usage:
  python cnv_to_bed.py -c cnvkit sample.call.cns
  python cnv_to_bed.py -c codex2 data/work/codex2/codex2.segments.filtered.txt
    (writes one BED per sample_name next to the input)
"""

import argparse
import os
import csv
from collections import defaultdict

class CNVProcessor:
    def __init__(self,ucsc=False):
        self.ucsc=ucsc

    def get_zygosity(self,cn_a,cn_b):
        if cn_a in ('NA',None) or cn_b in ('NA',None):
            return 'unknown'
        cn_a=int(cn_a)
        cn_b=int(cn_b)
        if cn_a==0 or cn_b==0:
            return 'loh'
        if cn_a!=cn_b:
            return 'imbalance'
        return 'nonloh'

    def get_cn_class(self,cn_total):
        cn_total=float(cn_total)
        if cn_total>4:
            return 'amp'
        if cn_total>2:
            return 'gain'
        if cn_total==2:
            return 'neutral'
        if cn_total>=1:
            return 'loss'
        if cn_total<1:
            return 'del'
        return 'unknown'

    def format_bed_row(self,chrom,start,end,cn_total,cn_a,cn_b):
        try:
            start=int(start)-1
            end=int(end)
            cn_total=float(cn_total)
            zyg=self.get_zygosity(cn_a,cn_b)
            cn_class=self.get_cn_class(cn_total)
            name=f"{end-start+1}bp;{zyg};{cn_class};A{cn_a};B{cn_b}"
            score=f"{int(cn_total)*100}"
            return {
                'chrom':chrom,
                'chromStart':f"{start}",
                'chromEnd':f"{end}",
                'name':name,
                'score':score,
                'strand':"+"
            }
        except (ValueError,TypeError):
            return None

class SequenzaProcessor(CNVProcessor):
    def process_row(self,row):
        try:
            return self.format_bed_row(
                row['chromosome'],
                row['start.pos'],
                row['end.pos'],
                int(row['CNt']),
                int(row['A']),
                int(row['B'])
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
            return self.format_bed_row(
                row['chromosome'],
                row['start'],
                row['end'],
                int(row['cn']),
                int(row['cn1']),
                int(row['cn2'])
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
            return self.format_bed_row(
                row['chr'],
                row['startpos'],
                row['endpos'],
                cn_a+cn_b,
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
                cn_b=int(row['lcn.em'])
                cn_a=int(row['tcn.em'])-cn_b
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

class Codex2Processor(CNVProcessor):
    #Total-CN only (integer mode). Alleles packed as NA for annotsv_parser.
    def process_row(self,row):
        try:
            return self.format_bed_row(
                row['chr'],
                row['st_bp'],
                row['ed_bp'],
                int(float(row['copy_no'])),
                'NA',
                'NA'
            )
        except (ValueError,KeyError,TypeError):
            return None

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('-c','--caller',choices=['sequenza','cnvkit','ascat','purecn','facets','codex2'],default='sequenza',help='CNV caller')
    p.add_argument('--ucsc',action='store_true',default=False,help='Format outfiles for UCSC browser')
    p.add_argument('--samples',default=None,help='Optional sample list (codex2): write empty BEDs for samples with no calls')
    p.add_argument('input_fp',nargs=argparse.REMAINDER,help='One or more input files')
    argv=p.parse_args()
    return vars(argv)

def write_bed(path,rows,bed_fields):
    with open(path,'w',newline='') as bed_file:
        writer=csv.DictWriter(bed_file,delimiter='\t',fieldnames=bed_fields)
        for out_row in rows:
            writer.writerow(out_row)

def main(argv=None):
    bed_fields=['chrom','chromStart','chromEnd','name','score','strand']
    argv=get_args() if argv is None else argv
    processors={
        'sequenza':SequenzaProcessor,
        'cnvkit':CNVkitProcessor,
        'ascat':ASCATProcessor,
        'purecn':PureCNProcessor,
        'facets':FacetsProcessor,
        'codex2':Codex2Processor
    }
    processor=processors[argv['caller']](ucsc=argv['ucsc'])

    for file in argv['input_fp']:
        base,ext=os.path.splitext(file)
        delim=',' if ext=='.csv' else '\t'
        with open(file,'r') as infile:
            reader=csv.DictReader(infile,delimiter=delim)
            if argv['caller']=='codex2':
                #Merged CODEX2 table -> one BED per sample next to input.
                by_sample=defaultdict(list)
                for row in reader:
                    out_row=processor.process_row(row)
                    if out_row is None or out_row['chrom'] in ['24','Y','chrY']:
                        continue
                    sample=row.get('sample_name') or row.get('sample')
                    if not sample:
                        raise ValueError(f"codex2 row missing sample_name: {row}")
                    by_sample[sample].append(out_row)
                outdir=os.path.dirname(file) or '.'
                want=set(by_sample)
                if argv.get('samples'):
                    with open(argv['samples'],'r') as sf:
                        want|={ln.strip() for ln in sf if ln.strip() and not ln.startswith('#')}
                for sample in sorted(want):
                    rows=by_sample.get(sample,[])
                    outfile=os.path.join(outdir,f"{sample}.codex2.bed")
                    write_bed(outfile,rows,bed_fields)
                    print(f"Wrote {outfile} ({len(rows)} intervals)")
            else:
                outfile=f"{base}.bed"
                rows=[]
                for row in reader:
                    out_row=processor.process_row(row)
                    if out_row is not None and out_row['chrom'] not in ['24','Y','chrY']:
                        rows.append(out_row)
                write_bed(outfile,rows,bed_fields)

if __name__=='__main__':
    main()



#import vcf
import csv
import argparse
import os
from itertools import dropwhile

"""
To refGene annoate pindel calls
Long format events do not translate to non-"standard vcf entries
ANNOVAR wont like it
"""
#pindel.vcf > pindel.avinput
#pindel.hg19_multianno.txt > pindel.hg19_multianno.report.tsv


#def _format_filter(self, flt):
#    if flt == []:
#        return 'PASS'
#    return self._stringify(flt, none='.', delim=';')

#def unused():
#    row=[record.CHROM,record.POS,record.end]
#    if isinstance(record.ALT[0],vcf.model._Substitution):
#        row.append(record.REF)
#        row.append(str(record.ALT[0]))
#    elif isinstance(record.ALT[0],vcf.model._SV):
#        if record.ALT[0].startswith('<'):
#            row+=['0','0']
#        else:
#            row.append(record.REF)
#            row.append(str(record.ALT[0]))
#    else:
#        raise
    #Now would add all VCF info, but Im only including INFO,FORMAT,SAMPLES
    #pyvcf would sort info by key_order from vcf write that way.
#    row.append(';'.join(_stringify_pair(x,y) for x,y in record.INFO.items()))
#    row.append(record.FORMAT)

#def _format_info(self, info):
#    if not info:
#        return '.'
#    def order_key(field):
#        # Order by header definition first, alphabetically second.
#        return self.info_order[field], field
#    return ';'.join(self._stringify_pair(f, info[f]) for f in
#                    sorted(info, key=order_key))

#def _format_sample(self, fmt, sample):
#    if hasattr(sample.data, 'GT'):
#        gt = sample.data.GT
#    else:
#        gt = './.' if 'GT' in fmt else ''
#    result = [gt] if gt else []
#    # Following the VCF spec, GT is always the first item whenever it is present.
#    for field in sample.data._fields:
#        value = getattr(sample.data,field)
#        if field == 'GT':
#            continue
#        if field == 'FT':
#            result.append(self._format_filter(value))
#        else:
#            result.append(self._stringify(value))
#    return ':'.join(result)

#Maybe I can put some of these inside each other


#def _stringify_pair(x,y,none='.',delim=','):
#    #From pyvcf
#    def _map(func,iterable,none='.'):
#        '''``map``, but make None values none.'''
#        return [func(x) if x is not None else none for x in iterable]
#    def _stringify(x,none='.',delim=','):
#        if type(x)==type([]):
#            return delim.join(_map(str,x,none))
#        return str(x) if x is not None else none
#    if isinstance(y,bool):
#        return str(x) if y else ""
#    return "%s=%s" % (str(x),_stringify(y,none=none,delim=delim))

def is_comment(s):
    return s.startswith('#')

def format_record(line):
    outrow=[]
    record=line.rstrip().split('\t')
    outrow+=record[:2]
    info=split_info(record[7])
    outrow.append(info['END'])
    if record[4].startswith('<') and record[4].endswith('>'):
        outrow+=['0','0']
    else:
        outrow+=record[3:5]
    outrow+=record[5:]
    return outrow

def split_info(info):
    _info={}
    items=info.split(';')
    for item in items:
        i=item.split('=')
        if len(i)==2:
            _info[i[0]]=i[1]
        elif len(i)==0:
            _info[i[0]]='True'
        else:
            print('INFO :',i)
    return _info


def get_snake():
    argv={}
    argv['input']=snakemake.input
    argv['output']=snakemake.output
    return argv

def get_args():
    argv={}
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input',help="VCF file made from 'pindel2vcf'.")#need gzip support
    args=p.parse_args()
    argv['input']=args.input
    argv['output']=os.path.splitext(args.input)[0]+'.avinput'
    return argv

def main(argv=None):
    with open(argv['output'],'w') as OUTFILE,open(argv['input'],'r') as vcfFILE:
        writer=csv.writer(OUTFILE,delimiter=' ')
        for record in dropwhile(is_comment,vcfFILE):
            writer.writerow(format_record(record))

if __name__=='__main__':
    try:
        snakemake
    except NameError:
        main(get_args())
    else:
        main(get_snake())
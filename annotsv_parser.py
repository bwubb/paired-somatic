import sys
import argparse
import csv
import os

csv.field_size_limit(sys.maxsize)

def name_info(v):
    y=v.split(';')
    # Check and convert A and B values safely
    A=y[3][1:] if y[3][1:].isdigit() or y[3][1:]=="NA" else "0"
    B=y[4][1:] if y[4][1:].isdigit() or y[4][1:]=="NA" else "0"
    A=int(A) if A!="NA" else A
    B=int(B) if B!="NA" else B
    x={'Segment.Length':y[0],'Segment.Zyg':y[1],'Segment.Type':y[2],'Segment.CNa':A,'Segment.CNb':B}
    x['Segment.CN']=x['Segment.CNa']+x['Segment.CNb'] if isinstance(x['Segment.CNa'], int) and isinstance(x['Segment.CNb'], int) else "NA"
    return x

def get_header():
    return 'TumorID,Gene,SV.Chrom,SV.Start,SV.End,Segment.Length,Segment.Zyg,Segment.Type,Segment.CNa,Segment.CNb,Segment.CN,Overlapped_CDS.PCT,Frameshift,Location1,Location2'.split(',')

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('-i','--input_fp',help='Input file')
    p.add_argument('-o','--output_fp',help='Output file')
    p.add_argument('--tumor',help='')
    argv=p.parse_args()
    return vars(argv)

def main(argv=None):
    header=get_header()
    argv=get_args()
    with open(argv['input_fp'],'r') as in_file, open(argv['output_fp'],'w',newline='') as out_file:
        reader=csv.DictReader(in_file,delimiter='\t')
        writer=csv.DictWriter(out_file,delimiter=',',fieldnames=header)
        writer.writeheader()
        for row in reader:
            outrow={'SV.Chrom':row['SV_chrom'],'SV.Start':row['SV_start'],'SV.End':row['SV_end'],'Gene':row['Gene_name']}
            outrow.update({'Overlapped_CDS.PCT':row['Overlapped_CDS_percent'],'Frameshift':row['Frameshift'].upper(),'Location1':row['Location'],'Location2':row['Location2']})
            outrow['TumorID']=argv['tumor']
            outrow.update(name_info(row['user#1']))
            if outrow['Segment.CN']=='NA':
                outrow['Segment.CN']=f"{int(row['user#2'])/100}"
            writer.writerow(outrow)

if __name__=='__main__':
    main()

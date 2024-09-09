import vcfpy
import argparse
import warnings

def get_args():
    p=argparse.ArgumentParser()
    p.add_argument('input',nargs='+',help='One or more Varscan2 fpfilter output files')
    argv=p.parse_args()
    return vars(argv)

def parse_snakemake():
    argv['input']=snakemake.input
    return argv

def preprocess_vcf(file):
    lines=[]
    with open(file,'r') as f:
        for line in f:
            if not line.startswith('#'):
                parts=line.split('\t')
                parts[6]=parts[6].replace(',', ';')
                lines.append('\t'.join(parts))
            else:
                lines.append(line)
    temp_file=file.replace('.out','.preprocessed.vcf')
    with open(temp_file,'w') as f:
        f.writelines(lines)
    return temp_file

def process_vcf(file):
    filters=["RefReadPos","RefDist3","RefMapQual","RefMMQS","RefAvgRL",
             "SomaticP","VarCount","VarFreq","VarReadPos","VarDist3",
             "VarMMQS","VarMapQual","RefBaseQual","VarBaseQual",
             "VarAvgRL","Strand","MaxBAQdiff","MMQSdiff","MinMMQSdiff",
             "MapQualDiff","ReadLenDiff","NoReadCounts"]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        reader=vcfpy.Reader.from_path(file)

        new_header=reader.header.copy()
        new_header.add_format_line(vcfpy.OrderedDict([('ID','AF'),('Number','A'),('Type','Float'),('Description','Allele Frequency')]))

        # Add missing filter lines
        for filter_id in filters:
            if not reader.header.has_header_line("FILTER",filter_id):
                new_header.add_filter_line(vcfpy.OrderedDict([('ID',filter_id),('Description','No Description provided by VarScan2 authors')]))

        output_file=file.replace('.out','.processed.vcf')
        writer=vcfpy.Writer.from_path(output_file,new_header)

        for record in reader:
            for call in record.calls:
                if 'AD' in call.data and 'DP' in call.data:
                    ad=call.data['AD']
                    dp=call.data['DP']
                    if dp>0:
                        af=[round(alt/dp,2) for alt in ad] if ad else [0.0]
                        call.data['AF']=af
            if len(record.FILTER)==1 and ',' in record.FILTER[0]:
                record.FILTER=record.FILTER[0].split(',')
            writer.write_record(record)

        writer.close()

def main(argv=None):
    argv=get_args()if argv is None else argv
    for file in argv['input']:
        #preprocessed_file=preprocess_vcf(file)
        process_vcf(file)

if __name__=="__main__":
    main()

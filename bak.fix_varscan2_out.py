
import csv

MISSING_FILTERS=["RefReadPos","RefDist3","RefMapQual","RefMMQS","RefAvgRL","SomaticP","VarCount","VarFreq","VarReadPos","VarDist3","VarMMQS","VarMapQual","RefBaseQual","VarBaseQual","VarAvgRL","Strand","MaxBAQdiff","MMQSdiff","MinMMQSdiff","MapQualDiff","ReadLenDiff","NoReadCounts"]

for i,file in enumerate(snakemake.input):
    with open(file,'r') as IN,open(snakemake.output[i],'w') as OUT:
        reader=csv.reader(IN,delimiter='\t')
        writer=csv.writer(OUT,delimiter='\t',quoting=csv.QUOTE_NONE,quotechar='',)
        for row in reader:
            if row[0].startswith('#CHROM'):
                for NAME in MISSING_FILTERS:
                    writer.writerow([f'##FILTER=<ID={NAME},Description="No Description provided by VarScan2 authors">'])
                writer.writerow(row)
            elif row[0].startswith('##'):
                writer.writerow(row)
            else:
                new_row=row
                new_row[6]=row[6].replace(',',';')
                writer.writerow(new_row)

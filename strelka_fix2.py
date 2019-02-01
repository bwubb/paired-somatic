import sys
import vcfpy

tumor=str(sys.argv[1])

myvcf=vcfpy.Reader.from_path(f"data/work/{tumor}/S07604715/concordant_calls/0001.rename.reorder.vcf")
outheader=myvcf.header
outheader.add_format_line(vcfpy.OrderedDict([('ID', 'GT'),('Number','1'),('Type','String'),('Description', 'Genotype')]))
writer=vcfpy.Writer.from_path(f"data/work/{tumor}/S07604715/concordant_calls/0001.rename.reorder.fix.vcf",outheader)
for record in myvcf:
    record.FORMAT.insert(0,'GT')
    for call in record.calls:
        call.data['GT']='0/1'
    writer.write_record(record)
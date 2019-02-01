
import sys
import vcfpy

def _tumor_normal_allele_freqs(ref,alt,calls):#unclear if Im doing this after establishing genotypes. For now assume one one variant
    """
    Retreive AD, AF fields for Strelka vcf files
    Makes best guess attempt
    """
    def _snv_freq(call):
        try:
            used_depth=call.data['DP']-call.data['FDP']#-call.data['SDP']
            _ad=call.data[f'{alt}U'][0]
            _af=float(_ad)/float(used_dp)
        except (ZeroDivisionError,TypeError) as e:
            print(f'{e}: {ref} {alt} {call}')
            _ad=0
            _af=0.00
        return _ad,float(f"{_af:.2f}")
    #DUP,DEL,INS,SNV,MNV
    def _indel_freq(call):
        try:
            used_Depth=call.data['DP']
            _ad=call.data['TIR']
            _af=float(_ad)/float(used_dp)
        except (ZeroDivisionError,TypeError) as e:
            print(f'{e}: {ref} {alt} {call}')
            _ad=0
            _af=0.00
        return _ad,float(f"{_af:.2f}")
    
    for call in calls:
        if alt[0].type=='SNV':
            _ad,_af=_snv_freq(call)
        elif alt[0].type in ['DEL','DUP','INS']:
            _ad,_af=_indel_freq(call)
        else:
            print(f'{alt[0].type} type UNKNOWN')
        call.data['AD']=_ad
        call.data['AF']=_af
    return calls

def _tumor_normal_genotypes(ref, alt, info):#, fname, coords):#record
    """Retrieve standard 0/0, 0/1, 1/1 style genotypes from INFO field.
    Normal -- NT field (ref, het, hom, conflict)
    Tumor -- SGT field
      - for SNPs specified as GG->TT for the normal and tumor diploid alleles. These
        can also represent more complex alleles in which case we set at heterozygotes
        pending longer term inclusion of genotypes in Strelka2 directly
        (https://github.com/Illumina/strelka/issues/16)
      - For indels, uses the ref, het, hom convention
    """
    known_names=set(["het","hom","ref","conflict"])
    def name_to_gt(val):
        if val.lower()=="het":
            return "0/1"
        elif val.lower()=="hom":
            return "1/1"
        elif val.lower() in set(["ref","confict"]):#Check examples after running bcftools norm. #Expand all variants and make separate 0/1 calls
            return "0/0"
        else:
            return "0/1"
    def alleles_to_gt(val):
        gt_indices = {gt.upper(): i for i, gt in enumerate([ref]+alt)} #combo list of ref and alt, enumerated for GT integers
        tumor_gts = [gt_indices[x.upper()] for x in val if x in gt_indices]
        if tumor_gts and val not in known_names:
            if max(tumor_gts)==0:
                tumor_gt="0/0"
            elif 0 in tumor_gts:
                tumor_gt=f"0/{min([x for x in tumor_gts if x > 0])}" 
            else:
                tumor_gt=f"{min(tumor_gts)}/{max(tumor_gts)}"
        else:
            tumor_gt=name_to_gt(val)
        return tumor_gt
    nt_val=info['NT']
    normal_gt=name_to_gt(nt_val)
    sgt_val=info.get('SGT',False)
    if not sgt_val:
        tumor_gt = "0/0"
    else:
        sgt_val=sgt_val.split("->")[-1]
        tumor_gt=alleles_to_gt(sgt_val)
return tumor_gt,normal_gt



def main(argv=None):
    myvcf=vcfpy.Reader.from_path(f"data/work/{tumor}/S07604715/strelka/rename.reorder.vcf.gz")
    outheader=myvcf.header
    outheader.add_format_line(vcfpy.OrderedDict([('ID', 'GT'),('Number','1'),('Type','String'),('Description', 'Genotype')]))
    writer=vcfpy.Writer.from_path(f"data/work/{tumor}/S07604715/strelka/rename.reorder.fix.vcf",outheader)#,vcfpy.Header(reader.header.lines,vcfpy.SamplesInfos([tumor,normal])))
    for record in myvcf:
        record.FORMAT.insert(0,'GT')
        for call in record.calls:
            call.data['GT']='0/1'
        writer.write_record(record)

if __name__=='__main__':
    p=argparse.ArgumentParser()
    p.add_argument('-I','--input',help='Input file')
    p.add_argument('-c','--caler',help='Caller vcf to reformat')
    argv=p.parse_args()
    main(args)
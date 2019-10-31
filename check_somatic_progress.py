
import os
import yaml
import argparse
import datetime

p=argparse.ArgumentParser()
p.add_argument('-c', '--configfile',help='snakemake config yaml')
args=p.parse_args()
#Check Somatic progess

missing_bam=[]
missing_vcf=[]
missing_pair=[]

with open(args.configfile,'r') as file:
    config=yaml.load(file)

with open(config['project']['sample_list']) as file:
    SAMPLES=file.read().splitlines()

for sample in SAMPLES:
    if not os.path.exists(f'bam_input/final/{sample}/{config["reference"]["key"]}/{sample}.ready.bam'):
        missing_bam.append(sample)

TUMORS=[t for t in SAMPLES if 'germline' not in t]
print('TUMORS')
print(TUMORS)
print()

for tumor in TUMORS:
    if not any([os.path.exists(f'data/final/{tumor}.{caller}.somatic.raw.vcf.gz') for caller in ['mutect2','strelka2','vardict','varscarn2']]):
        if tumor not in missing_vcf:
            missing_vcf.append(tumor)

with open(f'{config["resources"]["targets_key"]}.{datetime.date.today().strftime("%Y%m%d")}_pairs.list','w') as ofile:
    for tumor in TUMORS:
        root=tumor.split('-')[0]
        if f'{root}-germline' in SAMPLES:
            pass
            #ofile.write(f'{tumor}\t{root}-germline\n')
        elif os.path.exists(f'bam_input/final/{root}-germline/GRCh37/{root}-germline.ready.bam'):
            pass
        else:
            missing_pair.append(tumor)

print('missing_bam')
print(missing_bam)
print()

print('missing_vcf')
print(missing_vcf)
print()

print('missing_pair')
print(missing_pair)
print()
with open(f'{config["resources"]["targets_key"]}.{datetime.date.today().strftime("%Y%m%d")}_missing-pairs.list','w') as ofile:
    for s in missing_pair:
        ofile.write(f'{s}\n')
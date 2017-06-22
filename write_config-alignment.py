#!/usr/bin/env python

#20170217

import os
import glob
import argparse
import datetime
import errno

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno==errno.EEXIST and os.path.isdir(path):
			pass

def _get_ID_PU(fq):
	x=os.path.basename(fq).rstrip('.fastq.gz').split('_')
	return '.'.join([x[1],x[2],x[4]]),x[4]

def get_fastqs(sample,dir):
	path=os.path.abspath(dir)
	fqs=glob.glob('{0}/{1}_*.gz'.format(path,sample))#You see the problem with this.
	if len(fqs)<2:
		return False,None
	else:
		return sorted(fqs)

def ref_by_build(build):
	REF={'hg19':'','GRCh37':'/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta','hg38':'','GRCh38':'','GRCh38_GDC':'/home/bwubb/GRCh38_GDC/GRCh38.d1.vd1.fa'}
	return REF[build]

def write_script(sample,fqs,LB,ref,no_merge):#Project?
	with open('bam_input/config/{0}/aln.sh'.format(sample),'wb') as script:
		script.write('#!/bin/bash\n\n')
		script.write('ref="{0}"\n'.format(ref))#standardize?
		#script.write('PICARD=\n\n')
		script.write('module load picard-1.96\n\n')
		#script.write('SCRIPT=$(readlink -f "$0")\n')
		#script.write('cd $(dirname "$SCRIPT")\n\n')
		for i,j in enumerate(range(1,len(fqs),2)):
			ID,PU=_get_ID_PU(fqs[j-1])
			script.write('bwa mem -M -t 16 $ref {0} {1} | samtools view -bS - > {2}.bam\n'.format(fqs[j-1],fqs[j],ID))
			script.write('samtools sort -o {0}.COsorted.bam {0}.bam \n'.format(ID))
			script.write('if [ -e {0}.COsorted.bam ]; then rm {0}.bam; fi\n\n'.format(ID))
			#Add Readgroup
			#script.write('java -jar $PICARD AddOrReplaceReadGroups \\\n')#figure out picard tools commandline?
			script.write('AddOrReplaceReadGroups \\\n')
			script.write(' I={0}.COsorted.bam \\\n'.format(ID))
			script.write(' O={0}.COsorted.wRG.bam \\\n'.format(ID))
			script.write(' RGID={0} \\\n'.format(ID))
			script.write(' RGLB={0} \\\n'.format(LB))
			script.write(' RGPL=illumina \\\n')
			script.write(' RGPU={0} \\\n'.format(PU))
			script.write(' RGSM={0} \n\n'.format(sample))
		bam='{0}.{1}.input.bam'.format(sample,len(fqs)/2)
		if len(fqs)>2:
			script.write('samtools merge -f {0} *COsorted.wRG.bam\n\n'.format(bam))
		else:
			script.write('mv {0}.COsorted.wRG.bam {1}\n\n'.format(ID,bam))
		
		
		
		#Add verbosity option, or just always have it in output
		#script.write('samtools index {0}\n\n'.format(bam))
		script.write('if [ -e {0} ]; then rm *.COsorted*; fi\n\n'.format(bam))
		#script.write('samtools view -h {0} | samblaster -M | samtools view -Sb - > {1}\n'.format(bam,)#more bam names
		

def get_args():
	'''Parse sys.argv'''
	parser=argparse.ArgumentParser()
	parser.add_argument('-i','--infile', help='')
	parser.add_argument('--fastq-dir',default='FASTQ',help='Directory where fastq files may be found.') #Display default
	parser.add_argument('--lib',default='unknown',help='Name associated with capture library. Will become LB tag')
	parser.add_argument('--build',default='GRCh37',help='')
	parser.add_argument('--no-merge',action='store_true',default=False,help='Do no merge')
	#parser.add_argument('-w','--template',required=True,help='Name of bcbio template')
	#parser.add_argument('-r','--resource_file',help='resource file to include')
	return parser.parse_args()

def main(argv=None):
	args=get_args()#need to check values
	for a,b in vars(args).items():
		print '{}: {}'.format(a,b)
	with open(args.infile,'rb') as file:
		samples=file.read().splitlines()
	
	mkdir_p('bam_input')#baminput/config/{0} and final and work and work/{0}?
	mkdir_p('bam_input/final')
	mkdir_p('bam_input/config')
	mkdir_p('bam_input/work')
	for sample in samples:
		mkdir_p('bam_input/work/{0}'.format(sample))
		mkdir_p('bam_input/config/{0}'.format(sample))
		fqs=get_fastqs(sample,args.fastq_dir)
		
		ref=ref_by_build(args.build)
		if ref:
			write_script(sample,fqs,args.lib,ref,args.no_merge)
			print sample
#		if check:
#			write_yaml(sample,fqs)


if __name__=='__main__':
	main()

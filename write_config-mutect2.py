import errno
import os
import glob
import argparse
import csv


def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass

def write_script(tumor,normal,x,y):
	tbam='/media/Nandi/BAM_FILES/MSKCC/'+x
	nbam='/media/Nandi/BAM_FILES/MSKCC/'+y
	#for dir in ['BRCA','MORE_BRCA','MIXED']:
	#	if os.path.isfile('{0}/{1}.bam'.format(dir,tumor)):
	#		tbam=os.path.abspath('{0}/{1}.bam'.format(dir,tumor))
	#	if os.path.isfile('{0}/{1}.bam'.format(dir,normal)):
	#		nbam=os.path.abspath('{0}/{1}.bam'.format(dir,normal))
	with open('mutect/{0}/mutect.sh'.format(tumor), 'wb') as script:
		script.write('#!/bin/bash\n\n')
		script.write('SCRIPT=$(readlink -f "$0")\n')
		script.write('cd $(dirname "$SCRIPT")\n\n')
		script.write('ref="/media/Shakti/Genomes/Human/b37_genomes/human_g1k_v37.fasta"\n')
		script.write('dbsnp="/media/Shakti/Genomes/Human/b37_genomes/dbsnp_137.b37.vcf"\n')
		script.write('cosmic="/media/Shakti/Genomes/Human/b37_genomes/CosmicCodingMuts_v69_b37.vcf"\n\n')
		script.write('java -Xmx100g -jar /software/$GATK/GenomeAnalysisTK.jar \\\n')
		script.write('-T MuTect \\\n')
		script.write('-R $ref \\\n')
		script.write('--dbsnp $dbsnp \\\n')
		script.write('--cosmic $cosmic \\\n')
		#script.write('--intervals {} \\\n'.format(intervals_file))
		script.write('--input_file:normal {} \\\n'.format(nbam))
		script.write('--input_file:tumor {} \\\n'.format(tbam))
		script.write('--out {}.call_stats.txt \\\n'.format(tumor))
		script.write('--coverage_file {}.coverage.wig.txt \\\n'.format(tumor))
		script.write('--vcf {}.mutect.vcf\n'.format(tumor))
	print 'mutect/{}/mutect.sh was written.'.format(tumor)

def make_dirs(tumor,normal,):
	mkdir_p('mutect/{}'.format(tumor))

def main(argv=None):
	p = argparse.ArgumentParser()
	p.add_argument('-i', '--infile', help='Tab-delimited file; Tumor Normal pair IDs')
	args = p.parse_args()
	for a, b in vars(args).items():
		print '{}: {}'.format(a,b)
	
	with open(args.infile, 'rb') as file:
		pairs = [line.split('\t') for line in file.read().splitlines()]
		mkdir_p('mutect')
		for pair in pairs:
			mkdir_p('mutect/{0}'.format(pair[0]))
			write_script(*pair[:4 ])


if __name__ =='__main__':
	main()

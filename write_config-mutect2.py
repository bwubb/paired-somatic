import errno
import os
import glob
import argparse
import csv

class BAM():
	def __init__(self,name,bam):
		self.id=name
		if bam:
			self.bam=bam
		elif os.path.isfile('bam_input/final/{0}.ready.bam'.format(self.id)):
			self.bam=os.path.abspath('bam_input/final/{0}.ready.bam'.format(self.id))
		else:
			raise OSError('No file found: {0}\nTry including full path to bam in input file'.format(os.path.abspath('bam_input/final/{0}.ready.bam'.format(self.id))))

def mkdir_p(path):
	try:
		os.makedirs(path)
	except OSError as exc: # Python >2.5
		if exc.errno == errno.EEXIST and os.path.isdir(path):
			pass

def write_script(t,n):
	with open('data/final/{0}/mutect2.sh'.format(t.id),'wb') as script:
		script.write('#!/bin/bash\n\n')
		script.write('module load java-sdk-1.8.0\n\n')
		script.write('ref="/home/bwubb/resources/Genomes/Human/GRCh37/human_g1k_v37.fasta"\n')
		script.write('dbsnp="/home/bwubb/b37_genomes/dbsnp_137.b37.vcf"\n')
		script.write('cosmic="/home/bwubb/b37_genomes/CosmicCodingMuts_v69_b37.vcf"\n')
		script.write('GATK="/home/bwubb/software/GenomeAnalysisTK-3.7"\n\n')
		script.write('module load java-sdk-1.8.0\n\n')
		script.write('java -Xmx30g -jar /software/$GATK/GenomeAnalysisTK.jar \\\n')
		script.write('-T MuTect2 \\\n')
		script.write('-R $ref \\\n')
		script.write('--dbsnp $dbsnp \\\n')
		script.write('--cosmic $cosmic \\\n')
		#script.write('--intervals {} \\\n'.format(intervals_file))
		script.write('--input_file:normal {} \\\n'.format(n.bam))
		script.write('--input_file:tumor {} \\\n'.format(t.bam))
		script.write('--out {}.call_stats.txt \\\n'.format(t.id))
		script.write('--coverage_file {}.coverage.wig.txt \\\n'.format(t.id))
		script.write('--vcf {0}\n'.format(os.path.abspath('data/work/{0}/{0}.mutect.vcf\n'.format(t.id))))
	print 'data/config/{0}/mutect.sh was written.'.format(t.id)

def make_dirs(tumor,normal):
	mkdir_p('data/work/{}'.format(tumor))

def main(argv=None):
	p = argparse.ArgumentParser()
	p.add_argument('-i','--infile',help='Tab-delimited file; Tumor Normal pair IDs')
	p.add_argument('-L','--intervals',help='')
	#p.add_argument('--gatk',help='Specify a specific GATK version (if it exists in $HOME/software/)')
	#p.add_argument('--build',help='')#Add json reading for files of versions
	args = p.parse_args()
	for a, b in vars(args).items():
		print '{}: {}'.format(a,b)
	
	with open(args.infile, 'rb') as file:
		reader=csv.DictReader(file,delimiter='\t',fieldnames=['t_id','n_id','t_bam','n_bam'],restval=False)
		for row in reader:
			if any(v=k for k,v in row.items()):
				continue
			t=BAM(row['t_id'],row['t_bam'])
			n=BAM(row['n_id'],row['n_bam'])
			write_script(t,n)
			


if __name__ =='__main__':
	main()

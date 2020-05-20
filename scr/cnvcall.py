#!/usr/bin/env python
# coding=utf-8
import os
import click
from pathlib import Path
from configparser import ConfigParser

if os.path.exists('config.ini'):
	cfg = ConfigParser()
	cfg.read('config.ini')
	cnvcaller_path = cfg.get('configuration','cnvcaller_path')
	blasr_path = cfg.get('configuration','blasr_path')
	sawriter_path = cfg.get('configuration','sawriter_path')
	java8_path = cfg.get('configuration','java8_path')
	jar_picard = cfg.get('configuration','jar_picard')
	samtools_path = cfg.get('configuration','samtools_path')
else:
	cnvcaller_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/cnvcaller'
	blasr_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/blasr'
	sawriter_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/sawriter'
	java8_path ='/ifs4/BC_PUB/biosoft/pipeline/newblc/03.Soft_ALL/jdk1.8.0_131/bin/java'
	jar_picard = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/picard.jar'
	samtools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/samtools'

@click.command("cnvcall", help="cnv calling piplines")
@click.option(
    "-r", "--reference", type=click.Path(exists=True), help="the file of reference genome"
    )
@click.option(
    "-b", "--listfile", type=click.Path(exists=True), help="bam list"
    )
@click.option(
    "-o", "--outfile", default="cnv.sh", type=str, help="Output file", show_default=True)
def cnvcall(reference,listfile,outfile):
	refername = Path(reference).stem
	listfile = os.path.abspath(listfile)
	try:
		if not os.path.exists("06.Cnv_call"):
			os.mkdir('06.Cnv_call')
			os.chdir('06.Cnv_call')
		else: 
			print("the Folder of realign already exists")
			os._exit()
		with open(listfile) as IN,open('index_'+outfile,'w') as OUT1, open('work_'+outfile,'w') as OUT2:
			createDB = 'perl {}/bin/CNVReferenceDB.pl {} -w 800'.format(cnvcaller_path,reference)
			kmer = 'python {}/bin/0.1.Kmer_Generate.py {} 800 kmer.fa'.format(cnvcaller_path,reference)
			sawriter = '{} {}'.format(sawriter_path,reference)
			align = '{} kmer.fa {} --sa {}.sa --out kmer.aln -m 5 --noSplitSubreads --minMatch 15 --maxMatch 20 '\
					'--advanceHalf --advanceExactMatches 10 --fastMaxInterval --fastSDP --aggressiveIntervalCut --bestn 10'.format(blasr_path,reference,reference)
			link = 'python {}/bin/0.2.Kmer_Link.py kmer.aln 800 {}.link'.format(cnvcaller_path,refername)
			scr1= createDB+"\n"+kmer+"\n"+sawriter+"\n"+align+"\n"+link
			OUT1.write(scr1)

			tmp_scr = []
			count = 0
			for bam in IN:
				bamname = Path(Path(bam).stem).stem
				addrefer = '{} view -bT {} {} > {}.bam'.format(samtools_path,reference,bam.rstrip('\n'),bamname)
				AddOrReplaceReadGroups = '{} -Xmx5g  -jar {} AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT '\
										'I={}.bam O={}.head.bam ID={}  LB={}  PL=illumina SM={} PU={}'.format(java8_path,jar_picard,bamname,bamname,bamname,bamname,bamname,bamname)
				Individual = '{}/Individual.Process.sh -b {}.head.bam -h {} -d {}.link -s none'.format(cnvcaller_path,bamname,bamname,reference)
				tmp1_scr = addrefer+"\n"+AddOrReplaceReadGroups+"\n"+Individual
				tmp_scr.append(tmp1_scr)
				count += 1 #Count the number of files
			scr2 = '\n'.join(tmp_scr)
			if count <= 30: #Determine Pearson's coefficient
				pearson = 0.5
			elif 30< count <= 50:
				pearson = 0.4
			elif 50< count <= 100:
				pearson = 0.3
			elif 100< count <= 200:
				pearson = 0.2
			elif 200< count <= 500:
				pearson = 0.15
			else:
				pearson = 0.1
			RD_normalized_dict = str(Path.cwd())+'/RD_normalized'
			result_list = 'python /ldfssz1/MS_OP/USER/xutengfei1/script_py/abpath.py {} >list.txt'.format(RD_normalized_dict)
			exclude_list = "touch exclude_list.txt"
			discovery = 'bash {}/CNV.Discovery.sh -l list.txt -e exclude_list.txt  -f 0.1 -h 3 -r {} -p primaryCNVR -m mergeCNVR'.format(cnvcaller_path,pearson)
			genotype = 'python {}/Genotype.py --cnvfile mergeCNVR --outprefix {} --nproc 12'.format(cnvcaller_path,refername)
			scr3 = result_list+"\n"+exclude_list+"\n"+discovery+"\n"+genotype
			OUT2.write(scr2+"\n")
			OUT2.write(scr3)
	except TypeError :
			print("File not found,please input existing file or use option --help")
	except FileNotFoundError :
			print("File not found,please input existing file or use option --help")
	else:
		print("successfully! The script has been saved in 06.Cnv_call.")
#!/usr/bin/env python
# coding=utf-8

import os
import click
from pathlib import Path
import toml

if os.path.exists('./input.toml'):
	Configuration = toml.load("./input.toml")
	samtools_path = Configuration['configuration']['samtools_path']
	sambamba_path = Configuration['configuration']['sambamba_path']
	java7_path = Configuration['configuration']['java7_path']
	java8_path = Configuration['configuration']['java8_path']
	jar_gatk = Configuration['configuration']['jar_gatk']
	jar_picard = Configuration['configuration']['jar_picard']
else:
	samtools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/samtools'
	java7_path = '/ifs4/BC_PUB/biosoft/pipeline/Package/jre1.7.0_55/bin/java'
	java8_path ='/ifs4/BC_PUB/biosoft/pipeline/newblc/03.Soft_ALL/jdk1.8.0_131/bin/java'
	sambamba_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/sambamba-0.6.9-linux-static'
	jar_gatk = '/ifs4/BC_PUB/biosoft/pipeline/Package/GATK-3.3.0/GenomeAnalysisTK.jar'
	jar_picard = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/picard.jar'
@click.command("realign", help="revision of bam file")
@click.option(
    "-d", "--chrbam_dic", type=click.Path(exists=True), help="the dictionary of chrbam"
    )
@click.option(
    "-rd", "--refer_dic", type=click.Path(exists=True), help="the dictionary of split refer genome"
    )
@click.option(
    "-o", "--outfile", default="work_realign.sh", type=str, help="Output file", show_default=True)
def realign(chrbam_dic,refer_dic,outfile):
	try:
		if not os.path.exists("04.Realign"):
			os.mkdir('04.Realign')
			os.chdir('04.Realign')
		else: 
			print("the Folder of realign already exists")
			os._exit()
		current_dict = os.path.abspath('.')
		with open(outfile,'w') as OUT1,open('bam_list.sh','w') as OUT2:
			all_scr = []
			bamlist = []
			all_name = []
			#chrbam_dic = chrbam_dic.rstrip('/') if chrbam_dic.endswith('/') else chrbam_dic
			chrbam_dic = os.path.abspath(chrbam_dic)
			if chrbam_dic.endswith('/'):
				chrbam_dic = chrbam_dic.rstrip('/')
			else:
				pass
			chrbam_dic_path = Path(chrbam_dic)
			#for root,dirs,files in os.walk(chrbam_dic):
			tmp_iter = sorted([dirs for dirs in chrbam_dic_path.iterdir()]) #sort file list
			for dirs in tmp_iter:
				os.mkdir(Path(dirs).name)
				for fil in Path(dirs).iterdir():
					tmp = Path(fil).stem
					name = Path(tmp).stem
					all_name.append(name)
					fname2=os.path.join(dirs,fil)
					upper_directory = os.path.abspath(os.path.dirname(fname2))
					dic2 = Path(upper_directory).name
					dic = Path(upper_directory).name+"/"+Path(upper_directory).name+"."+name
					add_header = 'gunzip -c {} | {} view -bT {}/{}.fa - > ./{}.bam'.format(fname2,samtools_path,refer_dic,name,dic)
					sortbam = '{} sort -m 3G -t 8 --tmpdir ./ -o ./{}.sort.bam  ./{}.bam'.format(sambamba_path,dic,dic)
					add_RG = '{} -Xmx5g  -jar {} AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT   I=./{}.sort.bam O=./{}.head.bam '\
						     'ID={}  LB={}  PL=illumina SM={} PU={}'.format(java8_path,jar_picard,dic,dic,dic2,dic2,dic2,dic2)
					samtools_index = '{} index ./{}.head.bam'.format(samtools_path,dic)
					RealignerTargetCreator = '{} -Xmx5g  -jar {} -T RealignerTargetCreator -R {}/{}.fa -I ./{}.head.bam '\
						                     '-o ./{}.realn.intervals -allowPotentiallyMisencodedQuals'.format(java7_path,jar_gatk,refer_dic,name,dic,dic)
					IndelRealigner = '{} -Xmx5g  -jar {} -T IndelRealigner -R {}/{}.fa -I ./{}.head.bam  '\
					                 '-targetIntervals ./{}.realn.intervals -o ./{}.realign.bam  -allowPotentiallyMisencodedQuals'.format(java7_path,jar_gatk,refer_dic,name,dic,dic,dic)
					src = add_header+"\n"+sortbam+"\n"+add_RG+"\n"+samtools_index+"\n"+RealignerTargetCreator+"\n"+IndelRealigner
					all_scr.append(src)
			for  i in set(all_name):
				tmp_bamlist = 'ls {}/*/*{}.realign.bam >{}.realign.bam.list'.format(current_dict,i,i)
				bamlist.append(tmp_bamlist)
					
			final_scr = '\n'.join(all_scr)
			final_bamlist = '\n'.join(bamlist)
			OUT1.write(final_scr+"\n")
			OUT2.write(final_bamlist+"\n")

	except TypeError:
			print("File not found,please input existing file or use option --help")
	except FileNotFoundError :
			print("File not found,please input existing file or use option --help")
	else:
		print("successfully! The script has been saved in 04.Realign/{}.".format(outfile))

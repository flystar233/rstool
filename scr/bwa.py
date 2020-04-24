#!/usr/bin/env python
# coding=utf-8
import os
import click
from pathlib import Path

@click.command("bwa", help="The sequence will be aligned to reference genome by bwa")
@click.option(
    "-r", "--reference", type=click.Path(exists=True), help="the soft chain file of reference genome"
    )	
@click.option(
    "-l", "--listfile", type=click.Path(exists=True), help="the list file of clean data"
    )
@click.option(
    "-o", "--outfile", default="work_bwa.sh", type=str, help="Output file ", show_default=True)
def bwa(listfile,reference,outfile):
	#file format
	#clean_data	clean_data_1.fq.gz	clean_data_2.fq.gz
	#clean_data	clean_data2_1.fq.gz	clean_data2_2.fq.gz
	listfile = os.path.abspath(listfile)
	reference = os.path.abspath(reference)
	try:
		with open(listfile,'r') as IN1:
			all_scr1 = []
			all_scr2 = []
			ref_stem = Path(reference).stem
			dirname = str(Path(reference).parent)
			split_directory = dirname +'/{}_cut'.format(ref_stem)
			for _,_,files in os.walk(split_directory):
				for ref_chr in files:
					ref_chr_stem = Path(ref_chr).stem
					chr_ = Path(ref_chr_stem).stem
					gzip = 'gzip -d {}/{}'.format(split_directory,ref_chr)
					dic_split_chr = '{} -jar {} CreateSequenceDictionary R={}/{} O={}/{}.dict'.format(java8_path,jar_picard,split_directory,ref_chr_stem,split_directory,chr_)
					samtools_chr_index = '{} faidx {}/{}'.format(samtools_path,split_directory,ref_chr_stem)
					scr = gzip+"\n"+dic_split_chr+"\n"+samtools_chr_index
					all_scr1.append(scr)

			if not os.path.exists("03.Bwa"):
				os.mkdir('03.Bwa')
				os.chdir('03.Bwa')
				os.mkdir('chrbam')
			else: 
				print("the Folder of bwa already exists")
				os._exit()
			
			#file = IN1.readlines()
			for i in IN1:
				name,fq1,fq2 = i.split()
				os.mkdir(name)
				os.mkdir('chrbam/'+name)
				#out1 = Path(fq1).name
				#out2 = Path(fq2).name
				bwa = '{} mem -t 12 -M {} {} {} | {} view -Sb - -o {}/{}.bam'.format(bwa_path,reference,fq1,fq2,samtools_path,name,name)
				sortbam = '{} sort -m 3G -t 8 --tmpdir ./ -o {}/{}.sort.bam  {}/{}.bam'.format(sambamba_path,name,name,name,name)
				markduplicates = '{} -Xmx5g -Djava.io.tmpdir=tmp -jar {} MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=8000 INPUT={}/{}.sort.bam '\
					             'OUTPUT={}/{}.dedup.bam METRICS_FILE={}/{}.dedup.metrics VALIDATION_STRINGENCY=SILENT'.format(java8_path,jar_picard,name,name,name,name,name,name)
				samtools_index = '{} index {}/{}.dedup.bam'.format(samtools_path,name,name)
				bam_stat = '{} bamqc -bam {}/{}.dedup.bam -outfile {}.pdf -outformat PDF'.format(qualimap_path,name,name,name)
				uniq_map = '{} view -F 4 -q 30 {}/{}.dedup.bam |awk \'$7=="="\'|{} '\
						   'view -SbT {}  - > {}/{}.uniq.bam'.format(samtools_path,name,name,samtools_path,reference,name,name)
				chrbam = '{} Xamtools split -Ref {} -InFile {}/{}.uniq.bam -Bam -OutDir ./chrbam/{}'.format(iTools_path,reference,name,name,name)

				scr2 = bwa+"\n"+sortbam+"\n"+markduplicates+"\n"+samtools_index+"\n"+bam_stat+"\n"+uniq_map+"\n"+chrbam
				all_scr2.append(scr2)		

			with open ('index.sh','w') as OUT1,open(outfile,'w') as OUT2:
				final_scr1 = "\n".join(all_scr1)
				OUT1.write(final_scr1+"\n")
				final_scr2 = "\n".join(all_scr2)
				OUT2.write(final_scr2+"\n")

	except TypeError:
			print("File not found,please input existing file or use option --help")
	except FileNotFoundError :
			print("File not found,please input existing file or use option --help")
	else:
		print("successfully! The script has been saved in 03.Bwa/{}.".format(outfile))
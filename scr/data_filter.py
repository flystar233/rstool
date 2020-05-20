#!/usr/bin/env python
# coding=utf-8
import os
import click
from pathlib import Path
from configparser import ConfigParser

if os.path.exists('config.ini'):
	cfg = ConfigParser()
	cfg.read('config.ini')
	cut_adapter_path = cfg.get('configuration','cut_adapter_path')
	sickle_path = cfg.get('configuration','sickle_path')
else:
	cut_adapter_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/cutadapt'
	sickle_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/sickle-1.33/sickle'

@click.command("data_filter", help="raw data filter by cutadapt and sickle")
@click.option(
    "-l", "--listfile", type=click.Path(exists=True), help="the list file of raw data"
    )
@click.option(
    "-a", "--adapter", type=str,default=False,help=" the sequence of adapter. If don't cut adapter, please set -a F or False", show_default=True
    )	
@click.option(
    "-o", "--outfile", default="work_filter.sh", type=str, help="Output file ", show_default=True)
def data_filter(listfile,adapter,outfile):
	#file format
	#rawdata	raw_data_1.fq.gz	raw_data_2.fq.gz
	#raw_data2	raw_data2_1.fq.gz	raw_data2_2.fq.gz
	listfile = os.path.abspath(listfile) # Get absolute path
	try:
		with open(listfile,'r') as IN1:
			#file = IN1.readlines()

			if not os.path.exists("01.Data_filter"):
				os.mkdir('01.Data_filter')
				os.chdir('01.Data_filter')
			else: 
				print("the Folder of data_filter already exists")

				os._exit()
			with open(outfile,'w') as OUT1:
				for i in IN1: #meet big file
					name,fq1,fq2 = i.split()
					out1 = Path(fq1).name
					out2 = Path(fq2).name
					try:
						if  adapter:
							adapter = adapter
							r_adapter = adapter.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1] # reverse complementary sequence
							cut_adapter = '{} -a {} -A {} -e 0.1 -O 5 -m 50 -q 20,20 ' \
							'--discard-trimmed --pair-filter=any  -o trim_{} -p trim_{} {} {} >& {}_log.txt'.format(cut_adapter_path,adapter,r_adapter,out1,out2,fq1,fq2,name)
							cut_low_reads = '{} pe -f trim_{} -r trim_{} -t sanger '\
							'-o clean_{} -p clean_{} -s singles_{}.gz -l 50 -q 20 -g '.format(sickle_path,out1,out2,out1,out2,name)
							OUT1.write(cut_adapter+'\n')
							OUT1.write(cut_low_reads+'\n')
						else:
							cut_low_reads = '{} pe -f {} -r {} -t sanger '\
							'-o clean_{} -p clean_{} -s singles_{}.gz -l 50 -q 20 -g '.format(sickle_path,fq1,fq2,out1,out2,name)
							OUT1.write(cut_low_reads+'\n')
					except ValueError:
						print("Please check the list file for blank lines(delete it)")
						os._exit()

			with open('make_list.sh','w') as OUT2:
				pwd = os.getcwd()
				all_name  = 'cut -f1 {} |sort >all_name.txt'.format(listfile)
				allcleanfile = 'python /ldfssz1/MS_OP/USER/xutengfei1/script_py/abpath.py {} | grep clean_ |sort >allcleanfile.txt'.format(pwd)
				FQ1 = "sed -n '1~2p' allcleanfile.txt >fq1.txt"
				FQ2 = "sed -n '0~2p' allcleanfile.txt >fq2.txt"
				Clean_list = 'paste all_name.txt fq1.txt fq2.txt > clean_list.txt'
				all_scr = all_name+"\n"+allcleanfile+"\n"+FQ1+"\n"+FQ2+"\n"+Clean_list
				OUT2.write(all_scr+"\n")

	except TypeError:
			print("File not found,please input existing file or use option --help")
	except FileNotFoundError :
			print("File not found,please input existing file or use option --help")
	else:
		print("successfully! The script has been saved in 01.Data_filter/{}.".format(outfile))
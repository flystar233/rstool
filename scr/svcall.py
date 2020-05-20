#!/usr/bin/env python
# coding=utf-8
import os
import click
from pathlib import Path
from configparser import ConfigParser

#if os.path.exists('config.ini'):
if Path('config.ini').exists():
	cfg = ConfigParser()
	cfg.read('config.ini')
	samtools_path = cfg.get('configuration','samtools_path')
	breakdancer_path = cfg.get('configuration','breakdancer_path')
	java8_path = cfg.get('configuration','java8_path')
	jar_picard = cfg.get('configuration','jar_picard')
	smoove_path = cfg.get('configuration','smoove_path')
	sambamba_path = cfg.get('configuration','sambamba_path')
	manta_path = cfg.get('configuration','manta_path')
else:
	sambamba_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/sambamba-0.6.9-linux-static'
	samtools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/samtools'
	breakdancer_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin'
	java8_path ='/ifs4/BC_PUB/biosoft/pipeline/newblc/03.Soft_ALL/jdk1.8.0_131/bin/java'
	jar_picard = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/picard.jar'
	smoove_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/smoove'
	manta_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/manta-1.6.0.centos6_x86_64/bin/configManta.py'
@click.command("svcall", help="sv calling piplines (breakdancer+lumpy+manta)")
@click.option(
    "-r", "--reference", type=click.Path(exists=True), help="the file of reference genome"
    )
@click.option(
    "-b", "--listfile", type=click.Path(exists=True), help="bam list"
    )
@click.option(
    "-o", "--outfile", default="SV.sh", type=str, help="Output file", show_default=True)
def svcall(listfile,outfile,reference):
	refername = os.path.abspath(reference)
	listfile = os.path.abspath(listfile)
	try:
		if not Path('07.SV_call').exists():
			os.mkdir('07.SV_call')
			os.chdir('07.SV_call')
			os.mkdir('SV_RESULT')
		else: 
			print("the Folder of realign already exists")
			os._exit()
		with open(listfile) as IN,open('breakdancer_'+outfile,'w') as OUT1,open('smoove_'+outfile,'w') as OUT2,open('manta_'+outfile,'w') as OUT3:
			env1 = 'export PATH=/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT:$PATH'
			env2 = 'export PATH=/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/envs/python27/bin:$PATH'
			tmp_scr1 = []
			count = 0
			file1 = IN.readlines()
			for bam in file1:
				bamname = Path(Path(bam).stem).stem
				AddOrReplaceReadGroups = '{} -Xmx5g  -jar {} AddOrReplaceReadGroups VALIDATION_STRINGENCY=LENIENT '\
										'I={} O={}.head.bam ID={}  LB={}  PL=illumina SM={} PU={}'.format(java8_path,jar_picard,bam.strip('\n'),bamname,bamname,bamname,bamname,bamname)
				sort_bam = '{} sort  -m 3G -t 8 --tmpdir ./ -o ./{}.sort.bam  ./{}.head.bam'.format(sambamba_path,bamname,bamname)
				index = '{} index {}.sort.bam'.format(samtools_path,bamname)
				sample_bai = bamname+'.sort.bam.bai'
				rm_head = 'rm {}.head.bam'.format(bamname)
				tmp = AddOrReplaceReadGroups+"\n"+sort_bam+"\n"+index+"\n"+rm_head
				tmp_scr1.append(tmp)
				count += 1
			scr1 = '\n'.join(tmp_scr1)
			bam2cfg = 'perl {}/bam2cfg.pl -v 50 -gh *.sort.bam >./SV_RESULT/all_sample.breakdancer.cfg'.format(breakdancer_path,bam)
			breakdancer = '{}/breakdancer-max -q 35 ./SV_RESULT/all_sample.breakdancer.cfg >./SV_RESULT/all.breakdancer.ctx'.format(breakdancer_path)
			filt_sv = 'perl /ldfssz1/MS_OP/USER/xutengfei1/bin/filt_sv.pl -m 100 -x 1000000 -s 30 -d 5 -i ./SV_RESULT/all.breakdancer.ctx -o ./SV_RESULT/all.filter.breakdancer.ctx'
			scr2 = bam2cfg+"\n"+breakdancer+"\n"+filt_sv

			if count > 40:
				smoove_tmp = []
				merge_tmp = []
				genotype_tmp = []
				for bam in file1:
					bamname = Path(Path(bam).stem).stem
					smoove = '{} call --outdir ./SV_RESULT/results_smoove/  --name {} --fasta {} -p 1 --genotype {}.sort.bam'.format(smoove_path,bamname,refername,bamname)
					merge = '{} merge --name merged -f {} --outdir ./SV_RESULT/ ./SV_RESULT/results_smoove/*.genotyped.vcf.gz'.format(smoove_path,refername)
					genotype = '{} genotype -d -x -p 1 --name {}-joint --outdir ./SV_RESULT/results_genotped/ --fasta {} --vcf ./SV_RESULT/merged.sites.vcf.gz {}.sort.bam'.format(smoove_path,bamname,reference,bamname)
					smoove_tmp.append(smoove)
					merge_tmp.append(merge)
					genotype_tmp.append(genotype)
				merge_tmp = list(set(merge_tmp))[0]
				paste = '{} paste --name paste ./SV_RESULT/results_genotyped/*.vcf.gz'.format(smoove_path)
				scr_smoove = '\n'.join(smoove_tmp)
				scr_genotype = '\n'.join(genotype_tmp)
				scr3 = scr_smoove+"\n"+merge_tmp+"\n"+scr_genotype+"\n"+paste
			else:
				smoove = '{} call -x --name ./SV_RESULT/result_smoove  --fasta {} -p 4 --genotype *.sort.bam'.format(smoove_path,refername)
				scr3 = smoove

			manta_file = ' '.join(['--bam '+ bam.strip("\n") for bam in file1])
			manta_config = 'python {} {} --referenceFasta {} --runDir ./SV_RESULT/'.format(manta_path,manta_file,refername)
			manta_flow = 'echo "python runWorkflow.py" >./SV_RESULT/manta_flow.sh'
			manta_qsub = '#qsub -cwd -l vf=10G,num_proc=1 -q af.q ./SV_RESULT/manta_flow.sh'
			scr4 = manta_config+"\n"+manta_flow+"\n"+manta_qsub

			OUT1.write(scr1+"\n")
			OUT1.write(scr2+"\n")
			OUT2.write(env1+"\n"+env2+"\n")
			OUT2.write(scr3+"\n")
			OUT3.write(env2+"\n")
			OUT3.write(scr4+"\n")
	except TypeError :
			print("File not found,please input existing file or use option --help")
	except FileNotFoundError :
			print("File not found,please input existing file or use option --help")
	else:
		print("successfully! The script has been saved in 07.SV_call.")

#!/usr/bin/env python
# coding=utf-8
import os
import click
import re
from pathlib import Path
from configparser import ConfigParser

if os.path.exists('config.ini'):
	cfg = ConfigParser()
	cfg.read('config.ini')
	java7_path = cfg.get('configuration','java7_path')
	vcffilter_path = cfg.get('configuration','vcffilter_path')
	vcftools_path = cfg.get('configuration','vcftools_path')
	freebayes_path = cfg.get('configuration','freebayes_path')
	bcftools_path = cfg.get('configuration','bcftools_path')
	beagle_path = cfg.get('configuration','beagle_path')
	jar_gatk = cfg.get('configuration','jar_gatk')
	iTools_path = cfg.get('configuration','iTools_path')


else:
	java7_path = '/ifs4/BC_PUB/biosoft/pipeline/Package/jre1.7.0_55/bin/java'
	vcffilter_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/freebayes/vcflib/bin/vcffilter'
	vcftools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/vcftools_0.1.13/bin/vcftools'
	freebayes_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/freebayes/bin/freebayes'
	bcftools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/bcftools'
	beagle_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/beagle.28Sep18.793.jar'
	jar_gatk = '/ifs4/BC_PUB/biosoft/pipeline/Package/GATK-3.3.0/GenomeAnalysisTK.jar'
	iTools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/iTools_Code/iTools'

@click.command("snpcall", help="GATK + freebayes + bcftools snp calling piplines")
@click.option(
    "-r", "--reference", type=click.Path(exists=True), help="the folder of reference genome"
    )
@click.option(
    "-l", "--listfile", type=click.Path(exists=True), help="all bam list folder"
    )
@click.option(
    "-i","--imputation", default="false",type=click.Choice(['false','f','F','true','t','T'])
    )
@click.option(
    "-d","--dict_path", default="05.Snp_call",type=str,help="the home directory of snpcall", show_default=True)
@click.option(
    "-of", "--outfile_f", default="workf.sh", type=str, help="Output file", show_default=True)
@click.option(
    "-ob", "--outfile_b", default="workb.sh", type=str, help="Output file", show_default=True)
@click.option(
    "-og", "--outfile_g", default="workg.sh", type=str, help="Output file", show_default=True)
def snpcall(reference,listfile,outfile_f,outfile_b,outfile_g,imputation,dict_path):
	reference = os.path.abspath(reference)
	listfile = os.path.abspath(listfile)
	path1 = Path(reference)
	path2 = Path(listfile)
	falist = [child for child in path1.iterdir()]
	falist = [i for i in falist if Path(i).suffix=='.fa']
	falist = sorted(falist)
	filelist2 = [child for child in path2.iterdir()]
	filelist2 = [i for i in filelist2 if Path(i).suffix=='.list']
	filelist2 = sorted(filelist2)

	all1 = zip(falist,filelist2)
	all2 = zip(falist,filelist2)
	all3 = zip(falist,filelist2)
	try:
		if not os.path.exists(dict_path):
			os.mkdir(dict_path)
			os.chdir(dict_path)
		else:
			print("the Folder of Snp_call already exists")
			os._exit()
		if not os.path.exists("Freebayes_call"):
			os.mkdir('Freebayes_call')
			os.chdir('Freebayes_call')
		else:
			print("the Folder of freebayes_call already exists")
			os._exit()
		
		with open(outfile_f,'w') as OUT1:
			# freebayes piplines
			for i in all1:
				chrom = re.findall(r'chr[0-9]+',i[0].name,re.IGNORECASE)[0]
				snpcalling_f = '{} '\
				'--use-best-n-alleles 4 -p 2 -f {} -L {} >freebayes.{}.vcf'.format(freebayes_path,i[0],i[1],chrom)
				out1file_f = 'freebayes.{}.vcf'.format(chrom)
				vcffilter_f = '{} '\
				'-f "TYPE = snp & QUAL >20  & DP > 10 & AC > 0 & SAF > 0 & SAR > 0 & RPL >0 & RPR >0" {} >filter.{}.vcf'.format(vcffilter_path,out1file_f,chrom)
				out2file_f = 'filter.{}.vcf'.format(chrom)
				vcftools_f = '{}  '\
				'--vcf {} --min-alleles 2.0 --max-alleles 2.0 --max-missing 0.95 --non-ref-af 0.05 --max-non-ref-af 0.95 --recode -c >finnal.{}.vcf'.format(vcftools_path,out2file_f,chrom) # maf 0.05
				out3file_f = 'finnal.{}.vcf'.format(chrom)
				change_spot_f = 'perl -pe "s/\\t.:/\\t.\/.:/g" {} >finnal2.{}.vcf'.format(out3file_f,chrom)
				out4file_f = 'finnal2.{}.vcf'.format(chrom)

				final_scrf = snpcalling_f+"\n"+vcffilter_f+"\n"+vcftools_f+"\n"+change_spot_f
				OUT1.write(final_scrf+"\n")

				if imputation in ['true','t','T']:
					impute_f = 'java -jar {} gt={} ne=1000 out=imputation.{}.vcf'.format(beagle_path,out4file_f,chrom)
					OUT1.write(impute_f+"\n")
				else:
					pass

		os.chdir(os.path.abspath(os.path.dirname(os.getcwd())))
		if not os.path.exists("Bcftools_call"):
			os.mkdir('Bcftools_call')
			os.chdir('Bcftools_call')
		else:
			print("the Folder of bcftools_call already exists")
			os._exit()

		with open(outfile_b,'w') as OUT2:
			# bcftools piplines
			for i in all2:
				chrom = re.findall(r'chr[0-9]+',i[0].name,re.IGNORECASE)[0]
				snpcalling_b = '{} mpileup -Ou --skip-indels -f {} -b {} '\
				'|{} call -P 1.1e-5 -Ob -mv -o bcftools.{}.bcf'.format(bcftools_path,i[0],i[1],bcftools_path,chrom)
				out1file_b = 'bcftools.{}.bcf'.format(chrom)
				vcffilter_b = '{} filter -e "MQ < 40 || QUAL < 20 || DP <10 " -s LOWQUAL -Ou {} '\
				'|{} view -f PASS -Ov -o filter.{}.vcf'.format(bcftools_path,out1file_b,bcftools_path,chrom)
				out2file_b = 'filter.{}.vcf'.format(chrom)
				vcftools_b = '{}  --vcf {} '\
				'--min-alleles 2.0 --max-alleles 2.0 --max-missing 0.95 --non-ref-af 0.05 --max-non-ref-af 0.95 --recode -c >finnal.{}.vcf'.format(vcftools_path,out2file_b,chrom) # maf 0.05
				out3file_b = 'finnal.{}.vcf'.format(chrom)
				
				final_scrb = snpcalling_b+"\n"+vcffilter_b+"\n"+vcftools_b
				OUT2.write(final_scrb+"\n")

				if imputation in ['true','t','T']:
					impute_b = 'java -jar {} gt={} ne=1000 out=imputation.{}.vcf'.format(beagle_path,out3file_b,chrom)
					OUT2.write(impute_b+"\n")
				else:
					pass

		os.chdir(os.path.abspath(os.path.dirname(os.getcwd())))
		if not os.path.exists("Gatk_call"):
			os.mkdir('Gatk_call')
			os.chdir('Gatk_call')
		else:
			print("the Folder of gatk_call already exists")
			os._exit()

		if not os.path.exists('java_tmp'):
			os.makedirs('java_tmp')
		else:
			print("the Folder of java_tmp already exists")
			os._exit()
		for i in all3:
			os.mkdir(i[0].stem)
			os.chdir(i[0].stem)
			chrom = re.findall(r'chr[0-9]+',i[0].name,re.IGNORECASE)[0]
			with open(outfile_g,'w') as OUT3,open(i[1]) as in1:
				tmp_scr = []
				tmp_gvcf = []
				IN = in1.readlines()
				for chrlist in IN:
					sample = Path(chrlist).stem
					HaplotypeCaller = '{} -Xmx10G -Djava.io.tmpdir=java_tmp -jar {} -R {} -T HaplotypeCaller '\
					'--emitRefConfidence GVCF --variant_index_type LINEAR --variant_index_parameter 128000 -nct 8 -I {} -o {}.g.vcf'.format(java7_path,jar_gatk,i[0],chrlist.rstrip('\n'),sample)
					tmp_scr.append(HaplotypeCaller)
					tmp_gvcf.append(f'-V {sample}.g.vcf')

				final_scr = '\n'.join(tmp_scr)
				final_gvcf = ' '.join(tmp_gvcf)
				GenotypeGVCFs = '{} -Xmx10G -Djava.io.tmpdir=java_tmp -jar {}  -T GenotypeGVCFs -R {} {} -o {}.raw.snps.indels.vcf'.format(java7_path,jar_gatk,i[0],final_gvcf,chrom)
				VariantRecalibrator_snp = '{} -Xmx10G -jar {}  -T SelectVariants -selectType SNP  -R {} -V {}.raw.snps.indels.vcf -o {}.snp.vcf'.format(java7_path,jar_gatk,i[0],chrom,chrom)
				VariantFiltration_snp = '{} -Xmx5G -jar {}  -T VariantFiltration  -R {} --filterExpression "QD < 2.0 || MQ < 40.0 || ReadPosRankSum < -8.0 || FS > 60.0|| HaplotypeScore > 13.0 || MQRankSum < -12.5" '\
								' --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant {}.snp.vcf -o {}.snpfiltered.vcf'.format(java7_path,jar_gatk,i[0],chrom,chrom)
				grep_snp = 'grep -E "PASS|#" {}.snpfiltered.vcf  > {}.snp.final.vcf'.format(chrom,chrom)
				vcf2genotype = '{} Formtools VCF2Genotype  -InPut {}.snp.final.vcf -OutPut {}.snp.genotype -WithHeader'.format(iTools_path,chrom,chrom)
				VariantRecalibrator_indel = '{} -Xmx10G -jar {}  -T INDEL -selectType SNP  -R {} -V {}.raw.snps.indels.vcf -o {}.indel.vcf'.format(java7_path,jar_gatk,i[0],chrom,chrom)
				VariantFiltration_indel = '{} -Xmx5G -jar {}  -T VariantFiltration  -R {} --filterExpression "QD < 2.0 || MQ < 40.0 || ReadPosRankSum < -20.0 || FS > 200|| MQRankSum < -12.5||SOR > 10.0" '\
								' --filterName LowQualFilter --missingValuesInExpressionsShouldEvaluateAsFailing --logging_level ERROR --variant {}.indel.vcf -o {}.indelfiltered.vcf'.format(java7_path,jar_gatk,i[0],chrom,chrom)
				grep_indel = 'grep -E "PASS|#" {}.indelfiltered.vcf  > {}.indel.final.vcf'.format(chrom,chrom)
				final_scrg = final_scr+"\n"+GenotypeGVCFs +"\n"+VariantRecalibrator_snp +"\n"+VariantFiltration_snp +"\n"+grep_snp +"\n"+vcf2genotype+"\n"+VariantRecalibrator_indel+"\n"+VariantFiltration_indel+"\n"+grep_indel
				OUT3.write(final_scrg+"\n")
			os.chdir(os.path.abspath(os.path.dirname(os.getcwd())))
			
	except TypeError :
			print("File not found,please input existing file or use option --help")
	except FileNotFoundError :
			print("File not found,please input existing file or use option --help")
	else:
		print("successfully! The script has been saved in 05.Freebayes_call/{},06.Bcftools_call/{} and 07.Gatk_call/{}.".format(outfile_f,outfile_b,outfile_g))
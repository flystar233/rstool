#!/usr/bin/env python
# coding=utf-8
import os
import click
from pathlib import Path

@click.command("vcf_merge", help="Take the intersection of vcf files")
@click.option(
    "-a", "--listfile1", type=click.Path(exists=True), help="vcf file list from freebayes,the list must be compressioned by bgzip and create index"
    )
@click.option(
    "-b", "--listfile2", type=click.Path(exists=True), help="vcf file list from bcftools,the list must be compressioned by bgzip and create index"
    )
@click.option(
    "-o", "--outfile", default="work_vcf_merge.sh", type=str, help="Output file", show_default=True)
def vcf_merge(listfile1,listfile2,outfile):
	listfile1 = os.path.abspath(listfile1)
	listfile2 = os.path.abspath(listfile2)
	try:
		if not os.path.exists("08.Vcf_merge"):
			os.mkdir('08.Vcf_merge')
			os.chdir('08.Vcf_merge')
		else:
			print("the Folder of vcf_merge already exists")
			os._exit()
		with open(outfile,'w') as OUT:
			vcf_concat_f = '{} concat -f {} -o freebayes_concat.bcf -O b'.format(bcftools_path,listfile1)
			vcf_concat_b = '{} concat -f {} -o bcftools_concat.bcf -O b'.format(bcftools_path,listfile2)
			vcf_sort_f = '{} sort  freebayes_concat.bcf -o freebayes_sort.bcf -O b'.format(bcftools_path)
			vcf_sort_b = '{} sort  bcftools_concat.bcf -o bcftools_sort.bcf -O b'.format(bcftools_path)
			vcf_merge =  '{} merge  freebayes_sort.bcf  bcftools_sort.bcf -o merge.vcf.gz -O z'.format(bcftools_path)
			vcf2ped = '{}KIT/plink --vcf merge.vcf.gz --maf 0.05 --recode --out myplink --noweb'.format(plink_path)
			vcf2tped = '{}KIT/plink --vcf merge.vcf.gz --maf 0.05 --recode transpose --out myplink --noweb'.format(plink_path)
			vcf2bed = '{}KIT/plink --vcf merge.vcf.gz --maf 0.05 --make-bed --out myplink --noweb'.format(plink_path)
			
			final_scr = vcf_concat_f+"\n"+vcf_concat_b+"\n"+vcf_sort_f+"\n"+vcf_sort_b+"\n"+vcf_merge+"\n"+vcf2ped+"\n"+vcf2tped+"\n"+vcf2bed
			OUT.write(final_scr+"\n")

	except TypeError:
		print("File not found,please input existing file or use option --help")
	except FileNotFoundError :
		print("File not found,please input existing file or use option --help")
	else:
		print("successfully! The script has been saved in 08.Vcf_merge/{}.".format(outfile))
#!/usr/bin/env python
# coding=utf-8

import os
import click
import re
import csv
from pathlib import Path
from  collections import defaultdict
import toml
from scr.data_filter import data_filter
from scr.index import index
from scr.bwa import bwa
from scr.realign import realign
from scr.snpcall import snpcall
from scr.cnvcall import cnvcall
from scr.vcf_merge import vcf_merge
from scr.find_gene import find_gene

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS) # turn on subcommand options
@click.version_option(version='1.6')
def entrance(): # defind the entrance of method
	""" This tool can process raw fq data to clean vcf data.Each function is independent and does not affect each other.
	To prevent file accidents, please use the absolute path of files.
	"""
	"""                       
                                                                      $$$$              
                                                                      $$$$              
                                $$$                                   $$$$              
        $$$$$$$$   $$$$$$$$   $$$$$$$$    ;$$$$$$$       ;$$$$$$$     $$$$              
        $$$$$$*  $$$$    $$$  $$$$$$$$  +$$$$   $$$$   -$$$$   $$$$   $$$$              
        $$$$     $$$$$%         $$$     $$$       $$$  $$$       $$$  $$$$              
        $$$;        $$$$$$$$    $$$     $$$       $$$  $$$       $$$  $$$$              
        $$$;     $$      $$$$   $$$     $$$^     $$$$  $$$^     $$$$  $$$$              
        $$$;     $$$$$$$$$$$    $$$$$$$  $$$$$$$$$$~    $$$$$$$$$$~   $$$$              
        $$$;          $$$^         $$$$       $$             $$        $$$$
	"""
	pass
entrance.add_command(data_filter)
entrance.add_command(index)
entrance.add_command(bwa)
entrance.add_command(realign)
entrance.add_command(snpcall)
entrance.add_command(cnvcall)
entrance.add_command(vcf_merge)
entrance.add_command(find_gene)
if __name__ == "__main__":
    if  os.path.exists('./input.toml'):
        Configuration = toml.load("./input.toml")
        print("Find input.toml, will use software path in input.toml.")
        cut_adapter_path = Configuration['configuration']['cut_adapter_path']
        sickle_path = Configuration['configuration']['sickle_path']
        bwa_path = Configuration['configuration']['bwa_path']
        samtools_path = Configuration['configuration']['samtools_path']
        iTools_path = Configuration['configuration']['iTools_path']
        vcffilter_path = Configuration['configuration']['vcffilter_path']
        vcftools_path = Configuration['configuration']['vcftools_path']
        freebayes_path = Configuration['configuration']['freebayes_path']
        bcftools_path = Configuration['configuration']['bcftools_path']
        beagle_path = Configuration['configuration']['beagle_path']
        plink_path = Configuration['configuration']['plink_path']
        java7_path = Configuration['configuration']['java7_path']
        java8_path = Configuration['configuration']['java8_path']
        jar_gatk = Configuration['configuration']['jar_gatk']
        jar_picard = Configuration['configuration']['jar_picard']
        sambamba_path = Configuration['configuration']['sambamba_path']
        qualimap_path = Configuration['configuration']['qualimap_path']
        cnvcaller_path = Configuration['configuration']['cnvcaller_path']
        blasr_path = Configuration['configuration']['blasr_path']
        sawriter_path = Configuration['configuration']['sawriter_path']
    #software path
    else:
        print("Not find input.toml, will use defind software path.")
        cut_adapter_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/cutadapt'
        sickle_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/sickle-1.33/sickle'
        bwa_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/bwa-0.7.17/bwa'
        samtools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/samtools'
        qualimap_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/qualimap_v2.2.1/qualimap'
        iTools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/iTools_Code/iTools'
        vcffilter_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/freebayes/vcflib/bin/vcffilter'
        vcftools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/vcftools_0.1.13/bin/vcftools'
        freebayes_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/freebayes/bin/freebayes'
        bcftools_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/bcftools'
        beagle_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/beagle.28Sep18.793.jar'
        plink_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/plink'
        java7_path = '/ifs4/BC_PUB/biosoft/pipeline/Package/jre1.7.0_55/bin/java'
        java8_path ='/ifs4/BC_PUB/biosoft/pipeline/newblc/03.Soft_ALL/jdk1.8.0_131/bin/java'
        jar_gatk = '/ifs4/BC_PUB/biosoft/pipeline/Package/GATK-3.3.0/GenomeAnalysisTK.jar'
        jar_picard = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/picard.jar'
        sambamba_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/KIT/sambamba-0.6.9-linux-static'
        cnvcaller_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/cnvcaller'
        blasr_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/blasr'
        sawriter_path = '/zfssz3/NASCT_BACKUP/MS_PMO2017/xutengfei1/software/miniconda3/bin/sawriter'

    entrance()
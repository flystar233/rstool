#!/usr/bin/env python
# coding=utf-8

import os
import click
from configparser import ConfigParser
from scr.data_filter import data_filter
from scr.index import index
from scr.bwa import bwa
from scr.realign import realign
from scr.snpcall import snpcall
from scr.cnvcall import cnvcall
from scr.vcf_merge import vcf_merge
from scr.find_gene import find_gene
from scr.vcf_stat import vcf_stat
from scr.svcall import svcall
from scr.sv_merge import sv_merge

CONTEXT_SETTINGS = dict(help_option_names=['-h', '--help'])
@click.group(context_settings=CONTEXT_SETTINGS) # turn on subcommand options
@click.version_option(version='1.6.4')
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
entrance.add_command(svcall)
entrance.add_command(vcf_merge)
entrance.add_command(find_gene)
entrance.add_command(vcf_stat)
entrance.add_command(sv_merge)
if __name__ == "__main__":
    if  os.path.exists('config.ini'):
        print(">>> Find config.ini, will use software path in config.ini.")
    else:
        print(">>> Not find config.ini, will use defind software path.")

    entrance()

#!/usr/bin/env python
# coding=utf-8

import os
import click
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
    else:
        print("Not find input.toml, will use defind software path.")

    entrance()

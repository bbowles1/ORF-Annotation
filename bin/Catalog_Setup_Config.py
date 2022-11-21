#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 13:12:05 2022

Perform uORF Catalog setup by creating the following files:
    uorf_catalog.tsv (an adapted version of Mcgillivray et al's uORF catalog)
    uorf.bed (a bed file representing Mcgillivray et al's uORF catalog')
    range_df.pkl (a file contianing FASTA sequence and exon boundary positions for all 5'UTR regions)
    TITER_CDS.tsv (a file containing the last 100nt of sequence for all 5'UTR regions, needed to predict TITER mean ribosome load)

@author: Bradley Bowles
"""

import sys
import os
import argparse
import configparser

# parse inputs, if config path is provided
parser = argparse.ArgumentParser(description='Set up uORF Catalog.')
parser.add_argument("-c", "--config", nargs='?', type=str, default=None, const=None,
                    help='Path to config.txt. Defaults to searching in main for this file.')
parser.add_argument("-v", "--verbose", nargs='?', type=str, default='N', const='N',
                    help='Display python warnings? (Y/N)')
parser.add_argument("-t", "--titer", nargs='?', type=str, default='Y', const='Y',
                    help='Create a reference file for TITER analysis (Y/N).')
parser.add_argument("-u", "--update", nargs='?', type=str, default='N', const='N',
                    help='Update config file with paths to uorfs bed and catalog files (Y/N).')
args = parser.parse_args()

if args.verbose.lower() == 'n':
    import warnings
    warnings.filterwarnings('ignore') # setting ignore as a parameter


# set up module import
SCRIPT_DIR = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.abspath(os.path.sep.join([SCRIPT_DIR, "..", "modules"])))
import uorf_setup

# parse config file
conf = configparser.ConfigParser(inline_comment_prefixes="#")
try:
    if args.config == None:
        confpath = os.path.abspath(os.path.sep.join([SCRIPT_DIR, "..", "config.txt"]))
        conf.read(confpath)
    else:
        confpath = os.path.abspath(args.config)
        conf.read(confpath)
except:
    confpath = os.path.abspath(os.path.sep.join([SCRIPT_DIR, "..", "config.txt"]))
    conf.read(confpath)

# set variables
outdir = conf.get('all', 'outdir')
McGillivray_path = conf.get('all', 'McGillivray_path')
gff3_path = conf.get('all', 'gff3_path')
FASTA_path = conf.get('all', 'FASTA_path')
refseq_path = conf.get('all', 'refseq_path')
refseq_only = conf.get('all', 'refseq_only')
if refseq_only.lower()=='yes':
    refseq_only=True
else:
    refseq_only=False

# create uorf_catalog.tsv and uorf.bed
print("Creating uORF catalog and bed files.")
catalogpath, bedpath = uorf_setup.create_catalog(McGillivray_path, gff3_path, FASTA_path, outdir, refseq_only, refseq_path)

# update config with catalog paths
if args.update.lower()!='n':
    with open(confpath, "a") as file: # append mode
        file.write('catalogpath='+catalogpath+'\n')
        file.write('bedpath='+bedpath+'\n')


# create range_df.pkl
print("Creating backend file: range_df.pkl. This will take some time.")
uorf_setup.create_range_df(gff3_path, FASTA_path, outdir)

if args.titer.lower() == 'y':
    # create TITER_CDS.tsv
    print("Creating backend file: TITER_CDS.tsv. This will take some time.")
    titerpath = uorf_setup.create_TITER_df(gff3_path, FASTA_path, outdir)

    # update config with catalog paths
    if args.update.lower()!='n':
        with open(confpath, "a") as file: # append mode
            file.write('titerpath='+titerpath+'\n')

    print('Done! Created uorf_catalog.tsv, uorf.bed, range_df.pkl, and TITER_CDS.tsv.')
else:
    print('Done! Created uorf_catalog.tsv, uorf.bed, and range_df.pkl.')

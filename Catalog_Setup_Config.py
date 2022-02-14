#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 13:12:05 2022

Perform uORF Catalog setup by creating the following files:
    uorf_catalog.tsv (an adapted version of Mcgillivray et al's uORF catalog)
    uorf.bed (a bed file representing Mcgillivray et al's uORF catalog')
    range_df.pkl (a file contianing FASTA sequence and exon boundary positions for all 5'UTR regions)
    TITER_CDS.tsv (a file containing the last 100nt of sequence for all 5'UTR regions, needed to predict TITER mean ribosome load)

@author: m181414
"""

import uorf_setup

# output directory path for all created files
outdir = '/home/mayo/m181414/research/uorf_project/bior_script/'

# Define a path to your McGillivray et al. 2018 published catalog (doi: 10.1093/nar/gky188, Git: http://github.gersteinlab.org/uORFs/)
McGillivray_path = '/home/mayo/m181414/research/uorf_project/bior_script/complete_uORF_predictions_hg19.tsv'

# Define a path to an Ensembl .gff3 file (Ensembl release 87 or higher preferred)
# can be downloaded from http://ftp.ensembl.org/pub/grch37/current/gff3/homo_sapiens/Homo_sapiens.GRCh37.87.chr.gff3.gz
gff3_path = '/home/mayo/m181414/research/uorf_project/bior_script/Ensembl_GRCH37.87.gff3'

# Define a path to a reference FASTA (.fa) file - to check that site identities match a FASTA reference build
# this is required for create_range_df and create_range_df
FASTA_path = '/research/bsi/data/refdata/app/cava/human/2019/cgsl_wes/downloaded/2019_05_24/hg19_93.fa'

# Subset to RefSeq definitions only (using Ensembl:RefSeq matched transcripts download from Biomart)
refseq_only = True

# path to an Ensembl:RefSeq GRCh37 biomart file - can use the example file in the Git directory
# must have 'Transcript stable ID' and 'RefSeq mRNA ID' columns if you download your own biomart file
refseq_path = 'GRCh37_biomart_refseq_export.txt'



# create uorf_catalog.tsv and uorf.bed
uorf_setup.create_catalog(McGillivray_path, gff3_path, FASTA_path, outdir, True, refseq_path)

# create range_df.pkl
uorf_setup.create_range_df(gff3_path, FASTA_path, outdir)

# create TITER_CDS.tsv
uorf_setup.create_TITER_df(gff3_path, FASTA_path, outdir)

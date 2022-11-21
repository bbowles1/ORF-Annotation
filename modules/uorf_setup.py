#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Created on Wed Jun 24 15:09:28 2020

Perform uORF Catalog setup using the following scripts:
    create_catalog (writes uorf_catalog.tsv and uorf.bed)
    create_range_df (creates range_df.pkl)
    create_TITER_df (creates TITER_CDS.tsv)

@author: m181414
"""


import pandas as pd
import numpy as np
import pybedtools
import os


def complement_function(input_FASTA):  # This function translates negative strand nucleotides into their complements, but does not reverse the reading frame - must do this manually

    nucleotide_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A',
                       'a':'T', 'c':'G', 'g':'C', 't':'A'}

    input_FASTA = [nucleotide_dict[k] for k in input_FASTA]

    new_codon = ''.join(input_FASTA)

    return new_codon        # output new codons



def get_seq(BED_df, FASTA_path):

    # create a BED object and attempt to get sequence
    Bed_obj = pybedtools.BedTool.from_dataframe(BED_df) # create a BedTools object for pybedtools
    fasta = pybedtools.example_filename(FASTA_path)
    Bed_obj = Bed_obj.sequence(fi=fasta)
    FASTA_list = ((open(Bed_obj.seqfn).read())).split('\n') # opens FASTA text file and splits it to list
    FASTA_list = [ x for x in FASTA_list if ">" not in x ] # remove FASTA line header
    FASTA_list = [i for i in FASTA_list if i] # remove empty elements from list

    # if pybedtools return is empty, search FASTA for RefSeq identifiers instead
    if not bool(FASTA_list):
        # map FASTA names to correct identifiers
        try:
            FASTA_id_path = r'/Users/bbowles/Documents/PhD/uORF Work/uORF_Code/ORF-Annotation-main/build/FASTA_chrom_identifiers.txt'
            FASTA_ids = pd.read_csv(FASTA_id_path, sep='\t')[['Abbreviation', 'RefSeq sequence']].set_index(
                'Abbreviation').to_dict()['RefSeq sequence']
        except:
            print(
                'Could not locate a list of FASTA IDs. Is a FASTA_chrom_identifiers.txt included in the build folder?')

        # map FASTA IDs
        BED_df.loc[:, 'chrom'] = BED_df.chrom.map(FASTA_ids)

        # retry search
        Bed_obj = pybedtools.BedTool.from_dataframe(BED_df)  # create a BedTools object for pybedtools
        fasta = pybedtools.example_filename(FASTA_path)
        Bed_obj = Bed_obj.sequence(fi=fasta)
        FASTA_list = ((open(Bed_obj.seqfn).read())).split('\n')  # opens FASTA text file and splits it to list
        FASTA_list = [x for x in FASTA_list if ">" not in x]  # remove FASTA line header
        FASTA_list = [i for i in FASTA_list if i]  # remove empty elements from list

    return FASTA_list



def create_catalog(McGillivray_path, gff3_path, FASTA_path, outdir, refseq_only, refseq_path):

    # write catalog resources: a .tsv catalog and bed file
    # Pass None to FASTA_path to skip verification step

    # output catalog path
    outpath = os.path.join(outdir, 'uorf_catalog.tsv')

    # output bedpath
    bedpath = os.path.join(outdir, 'uorf.bed')

    # read in McGillivray et al dataframe
    ref_df = pd.read_csv(McGillivray_path, sep = '\t', header = 0, names = ['uORF_ID', 'uORF_score', 'start_codon', 'chromosome', 'strand', 'START', 'STOP', 'gene', 'additional_transcripts','type', 'peptide_score'],
                         dtype={'uORF_ID':str, 'uORF_score':float, 'start_codon':str, 'chromosome':str, 'strand':str, 'START':int, 'STOP':int, 'gene':str, 'additional_transcripts':str, 'type':str, 'peptide_score':str})



    # Perform formatting:

    # standardize 'None' instances across the dataframe
    ref_df.replace('None', np.NaN, inplace=True)

    # Drop -inf scoring uORFs from the McGillivray catalog - unlikely to be biologically relevant, and this improves pipeline efficiency
    ref_df = ref_df.loc[ref_df.uORF_score.astype(float) != np.NINF]

    # reformat negative strand columns to BED format (STOP > START)
    pos_df = ref_df.loc[ref_df.strand == '+']
    neg_df = ref_df.loc[ref_df.strand == '-']
    neg_df.rename(columns = {'START':'STOP','STOP':'START'}, inplace=True)
    ref_df = pd.concat([pos_df,neg_df], sort = True)


    # ADD ANNOTATION: Count associated uORFs for each gene
    count_dict = ref_df.gene.value_counts().to_dict()
    ref_df.loc[:, 'associated_uorfs'] = ref_df.gene.map(count_dict)


    # ADD ANNOTATION: start position of the CDS for the transcript in question
    ref_df['gene_start'] = np.NaN

    # Annotate gene_start: the start of coding at the canonical transcript ATG. This is the CDS start position in the Ensembl .gff3.
    ensg_df = pd.read_csv(gff3_path, header=0, sep='\t', comment='#', names=['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes'])

    # create column to define CDS start_sites
    ensg_df = ensg_df.loc[ensg_df.type == 'CDS'] # subset ENSG dataframe to include only CDS regions
    ensg_df['transcript'] = ensg_df.attributes.str.split(';').str[1]
    ensg_df.loc[:,'transcript'] = ensg_df.transcript.str.split(':').str[1]

    CDS_dict = {}
    for row in ensg_df.loc[ensg_df.strand == '+'].sort_values(by = 'start', ascending = True).drop_duplicates(subset = 'transcript', keep = 'first').itertuples():
        CDS_dict.update({row.transcript : row.start})  # create a dict of all ENSG:CDS_start pairs

    for row in ensg_df.loc[ensg_df.strand == '-'].sort_values(by = 'stop', ascending = False).drop_duplicates(subset = 'transcript', keep = 'first').itertuples():
        CDS_dict.update({row.transcript : row.stop})  # create a dict of all ENSG:CDS_start pairs

    ref_df.loc[:,'gene_start'] = ref_df.uORF_ID.str.split('.').str[0].map(CDS_dict)



    # check start/stop definitions - start site POS must always be less than stop site POS
    if not ref_df.loc[(ref_df.strand == '+') & (ref_df.gene_start < ref_df.START)].empty:
        raise Exception('Some + strand sites are in the wrong orientation (start POS > stop POS)')
    if not ref_df.loc[(ref_df.strand == '-') & (ref_df.gene_start > ref_df.STOP)].empty:
        raise Exception('Some - strand sites are in the wrong orientation (start POS > stop POS)')



    # Drop entries with start errors (the FASTA reference does not match the putative uORF start codon sequence)
    if bool(FASTA_path):


        pos_df = ref_df.loc[ref_df.strand == '+']
        neg_df = ref_df.loc[ref_df.strand == '-']

        # locate incorrect FASTA entries
        pos_df['chromStart'] = pos_df.START.astype(int) - 1
        pos_df['chromEnd'] = pos_df.chromStart.astype(int) + 3 # grabs the end of the new FASTA sequence
        BED_df = pos_df[['chromosome','chromStart','chromEnd']]
        BED_df.rename(columns={'chromosome':'chrom'}, inplace=True)
        
        # get FASTA sequences
        FASTA_list = get_seq(BED_df, FASTA_path)

        pos_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the pos_df column
        # pos_df now contains a corrected FASTA sequence at the site, which needs to be converted to the rev. complement for negative strands
        # drop incorrect FASTA entries

        pos_df = pos_df.loc[~(pos_df.FASTA != pos_df.start_codon)] # removes 670 rows

        # repeat process for neg_df

        neg_df['chromStart'] = neg_df.STOP.astype(int) - 3
        neg_df['chromEnd'] = neg_df.STOP.astype(int) # grabs the end of the new FASTA sequence
        BED_df = neg_df[['chromosome','chromStart','chromEnd']]
        BED_df.rename(columns={'chromosome':'chrom'}, inplace=True)

        # get FASTA sequences
        FASTA_list = get_seq(BED_df, FASTA_path)

        neg_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the neg_df column
        # neg_df now contains a corrected FASTA sequence at the site, which needs to be converted to the rev. complement for negative strands
        # drop incorrect FASTA entries
        neg_df = neg_df.loc[~neg_df.FASTA.str.contains('N')]
        neg_df.loc[:,'FASTA'] = neg_df.FASTA.apply(lambda x: x[::-1])
        neg_df.loc[:,'FASTA'] = neg_df.FASTA.apply(list).apply(complement_function)

        neg_df = neg_df.loc[~(neg_df.FASTA != neg_df.start_codon)] # removes 708 rows

        # Concat pos and neg df
        ref_df = pd.concat([pos_df, neg_df], sort=True)
        ref_df.sort_index(inplace=True)
        ref_df.drop(labels = ['FASTA', 'chromEnd', 'chromStart'], axis = 1, inplace = True)

        # drop 'chr' prefix from ref_df chromosome entries
        ref_df.loc[:, 'chromosome'] = ref_df.chromosome.str.strip('chr')



    # Add 5'UTR sites

    # Create a comprehensive ENST:ENSG dict
    ensg_names = ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes']
    ensg_dtypes = {'seqid':str, 'source':str, 'type':str, 'start':int, 'stop':int, 'score':str, 'strand':str, 'phase':str, 'attributes':str}
    ensg_df = pd.read_csv(gff3_path, header=0, sep='\t', comment='#', names=ensg_names, dtype = ensg_dtypes)

    ensg_df['ENST'] = np.NaN
    ensg_df['ENSG'] = np.NaN
    ensg_df.loc[(ensg_df.attributes.str.contains('ENSG')) & (ensg_df.attributes.str.contains('ENST')), 'ENST'] = ensg_df.loc[(ensg_df.attributes.str.contains('ENSG')) & (ensg_df.attributes.str.contains('ENST'))].attributes.str.split(';').str[0].str.split(':').str[1]
    ensg_df.loc[(ensg_df.attributes.str.contains('ENSG')) & (ensg_df.attributes.str.contains('ENST')), 'ENSG'] = ensg_df.loc[(ensg_df.attributes.str.contains('ENSG')) & (ensg_df.attributes.str.contains('ENST'))].attributes.str.split(';').str[1].str.split(':').str[1]

    ENST_ENSG_dict = {}
    for row in ensg_df.loc[(ensg_df.attributes.str.contains('ENSG')) & (ensg_df.attributes.str.contains('ENST'))].itertuples():
        ENST_ENSG_dict.update({row.ENST : row.ENSG})

    # create ENSG : gene dictionary
    ensg_df['gene'] = np.NaN
    ensg_df.loc[(ensg_df.type == 'gene'), 'gene'] = ensg_df.attributes.str.split(';').str[1].str.split('=').str[1]
    ensg_df.loc[(ensg_df.type == 'gene'), 'ENSG'] = ensg_df.loc[(ensg_df.type == 'gene')].attributes.str.split(';').str[0].str.split(':').str[1]

    ENSG_gene_dict = {}
    for row in ensg_df.loc[(ensg_df.type == 'gene')].itertuples():
        ENSG_gene_dict.update({row.ENSG : row.gene})

    # create UTR_df
    UTR_df = ensg_df.loc[ensg_df.type == 'five_prime_UTR'][['seqid', 'start', 'stop', 'strand', 'attributes']]
    UTR_df.rename(columns = {'start' : 'START', 'stop' : 'STOP', 'seqid' : 'chromosome'}, inplace = True)

    # create a list of all 23 chromosomes
    chrom_list = [str(i) for i in (list(range(1,23)))]
    chrom_list.append('X')
    chrom_list.append('Y')

    # Map gene names
    UTR_df['gene'] = UTR_df.attributes.str.split(':').str[1].map(ENST_ENSG_dict).map(ENSG_gene_dict) # get ENST info for UTR
    UTR_df = UTR_df.loc[UTR_df.gene.notna()] # drops 185 rows with gene IDs that could not be mapped

    # Identify regions with faulty chromosome annotations
    chrom_dict = {}
    for row in ensg_df.loc[ensg_df.type == 'gene'].itertuples():
        chrom_dict.update({row.gene : row.seqid})
    UTR_df.loc[~(UTR_df.chromosome.isin(chrom_list)), 'chromosome'] = UTR_df.loc[~(UTR_df.chromosome.isin(chrom_list))].gene.map(chrom_dict)
    UTR_df = UTR_df.loc[UTR_df.chromosome.isin(chrom_list)]
    # drops ~ 40 entries from the UTR_df

    # add additional UTR annotations to match the McGillivray catalog
    UTR_df['additional_transcripts'] = np.NaN
    UTR_df['peptide_score'] = np.NaN
    UTR_df['start_codon'] = np.NaN
    # strand already mapped
    UTR_df['type'] = 'gene_UTR'
    UTR_df['uORF_ID'] = UTR_df.attributes.str.split(':').str[1] # maps in ENST#
    UTR_df['uORF_score'] = np.NaN
    UTR_df['associated_uorfs'] = UTR_df.gene.map(count_dict)
    UTR_df['gene_start'] = UTR_df.uORF_ID.map(CDS_dict)
    UTR_df.drop(['attributes'], axis = 1, inplace = True)

    # Concat UTR regions to ref_df
    ref_df = pd.concat([ref_df, UTR_df], sort=True)
    ref_df.reset_index(drop=True, inplace=True)



    # Double check gene_start annotations across the dataframe for the 5'UTR sites we just added
    if bool(FASTA_path):

        start_codons = ['ATG', 'TTG', 'GTG', 'CTG','AGG', 'ACG', 'ATA', 'ATT', 'ATC'] # all possible start codons, AAG is not being considered due to very low ribosome loading in literature

        pos_df = ref_df.loc[(ref_df.strand == '+') & (ref_df.gene_start.notna())]
        pos_df.loc[:, 'gene_start'] = pos_df.gene_start.astype(float).apply(int)
        neg_df = ref_df.loc[(ref_df.strand == '-') & (ref_df.gene_start.notna())]
        neg_df.loc[:, 'gene_start'] = neg_df.gene_start.astype(float).apply(int)

        # locate incorrect FASTA entries
        pos_df['chromStart'] = pos_df.gene_start.astype(int) - 1
        pos_df['chromEnd'] = pos_df.chromStart.astype(int) + 3 # grabs the end of the new FASTA sequence
        BED_df = pos_df[['chromosome','chromStart','chromEnd']]
        BED_df.rename(columns={'chromosome':'chrom'}, inplace=True)
        BED_df.loc[:, 'chrom'] = 'chr' +  BED_df.chrom.astype(str)

        # get FASTA sequences
        FASTA_list = get_seq(BED_df, FASTA_path)

        pos_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the pos_df column
        # pos_df now contains a corrected FASTA sequence at the site, which needs to be converted to the rev. complement for negative strands
        # drop incorrect FASTA entries

        # adjust gene_starts that do not match input FASTA
        pos_df.loc[(~pos_df.FASTA.isin(start_codons)), 'gene_start'] = np.NaN

        # repeat process for neg_df

        neg_df['chromStart'] = (neg_df.gene_start.apply(int)) - 3
        neg_df['chromEnd'] = neg_df.gene_start.astype(int) # grabs the end of the new FASTA sequence
        BED_df = neg_df[['chromosome','chromStart','chromEnd']]
        BED_df.rename(columns={'chromosome':'chrom'}, inplace=True)
        BED_df.loc[:, 'chrom'] = 'chr' +  BED_df.chrom.astype(str)
        
        # get FASTA sequences
        FASTA_list = get_seq(BED_df, FASTA_path)

        neg_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the neg_df column
        # neg_df now contains a corrected FASTA sequence at the site, which needs to be converted to the rev. complement for negative strands
        # drop incorrect FASTA entries
        neg_df = neg_df.loc[~neg_df.FASTA.str.contains('N')]
        neg_df.loc[:,'FASTA'] = neg_df.FASTA.apply(lambda x: x[::-1])
        neg_df.loc[:,'FASTA'] = neg_df.FASTA.apply(list).apply(complement_function)

        # adjust gene_starts that do not match input FASTA
        neg_df.loc[(~neg_df.FASTA.isin(start_codons)), 'gene_start'] = np.NaN

        # Concat pos and neg df
        ref_df = pd.concat([pos_df, neg_df], sort=True)
        ref_df.sort_index(inplace=True)
        ref_df.drop(labels = ['FASTA', 'chromEnd', 'chromStart'], axis = 1, inplace = True)

        # drop 'chr' prefix from ref_df chromosome entries
        ref_df.loc[:, 'chromosome'] = ref_df.chromosome.str.strip('chr')



    # read in all ensembl transcripts that have a curated refseq protein coding transcript
    if refseq_only == True:
        refseq_df = pd.read_csv(refseq_path, sep = '\t', comment = '#')
        refseq_df.rename(columns = {'Gene stable ID':'ENSG','Transcript stable ID':'ENST', 'RefSeq mRNA ID':'refseq'}, inplace=True)
        ref_list = refseq_df.loc[refseq_df.refseq.notna()].ENST.unique()
        ref_df = ref_df.loc[ref_df.uORF_ID.str.split('.').str[0].isin(ref_list)]



    # get error stats for dataframe
    if not ref_df.loc[ref_df.gene.isna()].empty:
        print('Empty gene entries!')
    if not ref_df.loc[ref_df.uORF_ID.isna()].empty:
        print('Empty uORF_ID entries!')()
    if not ref_df.loc[(ref_df.STOP < ref_df.START)].empty:
        print('STOP precedes START for some rows!')
    if (ref_df.loc[(ref_df.strand == '+') & (ref_df.gene_start < ref_df.START)].shape[0]) + (ref_df.loc[(ref_df.strand == '-') & (ref_df.gene_start > ref_df.STOP)].shape[0]) != 0:
        print('gene_start entries are incorrect!')
    CDS_error = ref_df.loc[(ref_df.strand == '+') & (ref_df.gene_start < ref_df.STOP)].shape[0] + ref_df.loc[(ref_df.strand == '-') & (ref_df.gene_start > ref_df.START)].shape[0]
    print(ref_df.shape[0], ' rows.')
    print(CDS_error, ' entries with an intronic uORF STOP site.')
    if not ref_df.loc[ref_df.chromosome.str.contains('chr')].empty:
        print('Dataframe contains chromosome entries with "chr" prefix.')
    if not ref_df.loc[~ref_df.chromosome.isin(chrom_list)].empty:
        print('Some entries in the chromosome column are not chromosomes.')



    # Final reformatting steps
    # reorder and rename columns
    ref_df = ref_df[[ 'chromosome', 'START', 'STOP', 'strand', 'uORF_ID', 'type', 'gene', 'gene_start', 'associated_uorfs',
     'additional_transcripts', 'peptide_score', 'start_codon', 'uORF_score']]
    ref_df.rename(columns={'chromosome':'#CHROM'}, inplace=True)

    # replace NaN instances with '.'
    ref_df.replace(np.NaN, '.', inplace=True)


    # create a list of all 23 chromosomes
    chrom_list = [str(i) for i in (list(range(1,23)))]
    chrom_list.append('X')
    chrom_list.append('Y')

    # sort by CHROM
    ref_df['#CHROM'] = pd.Categorical(ref_df['#CHROM'], chrom_list)
    ref_df.sort_values(by=['#CHROM','START'], ascending=True, inplace=True)

    # write .tsv catalog to outpath
    ref_df.to_csv(outpath, sep='\t', index=False)
    print('Completed:', outpath)


    # write .bed file to bedpath
    bed_df = ref_df[['#CHROM', 'START', 'STOP']]
    bed_df.rename(columns={'#CHROM':'chrom', 'START':'chromStart', 'STOP':'chromEnd'}, inplace=True)
    # convert BED file to 0-based chromosome locations
    bed_df.loc[:, 'chromStart'] = bed_df.chromStart.astype(int) - 1
    bed_df.loc[:, 'chromEnd'] = bed_df.chromEnd.astype(int) - 1
    bed_df.sort_values(by=['chrom', 'chromStart'], inplace=True)
    bed_df.loc[:, 'chrom'] = 'chr' + bed_df.chrom.astype(str)
    bed_df.to_csv(bedpath, sep='\t', header=False, index=False)
    print('Completed:', bedpath)


    # return file paths to update config
    return outpath, bedpath



def create_range_df(gff3_path, FASTA_path, outdir):

    # Goal: get UTR sequence and UTR exon bounds for each ENST in the ensembl GRCh37 release
    # This is a .pkl file provided to the uORF interpreter script


    # Begin gathering FASTA sequences

    ensg_names = ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes']
    ensg_dtypes = {'seqid':str, 'source':str, 'type':str, 'start':int, 'stop':int, 'score':str, 'strand':str, 'phase':str, 'attributes':str}
    ensg_df = pd.read_csv(gff3_path, header=0, sep='\t', comment='#', names=ensg_names, dtype = ensg_dtypes)

    # create a list of all 23 chromosomes
    chrom_list = [str(i) for i in (list(range(1,23)))]
    chrom_list.append('X')
    chrom_list.append('Y')

    # remove non-chrom entries from seqid column
    ensg_df = ensg_df.loc[ensg_df.seqid.isin(chrom_list)]

    # subset to 5'UTR regions
    ensg_df = ensg_df.loc[ensg_df.type == 'five_prime_UTR']

    # define a transcript annotation (without version number)
    ensg_df['transcript'] = ensg_df.attributes.str.split(':').str[1]

    # get FASTA sequence for each 5'UTR exon
    BED_df = ensg_df[['seqid','start','stop']]
    BED_df.rename(columns = {'seqid':'chrom', 'start':'chromStart', 'stop':'chromEnd'}, inplace=True)
    BED_df.loc[:, 'chrom'] = 'chr'+ BED_df.chrom
    BED_df.loc[:, 'chromStart'] = BED_df.chromStart - 1 # adjusts to zero-based BEDtools nucleotide assignment

    # get FASTA sequences
    print("Obtaining FASTA sequences for all 5'UTR regions.")
    FASTA_list = get_seq(BED_df, FASTA_path)

    ensg_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the pos_df column


    # cut out duplicate entries in df, preferentially keeping entries in ensembl_havana, then havana, then ensembl
    ID_list = pd.Series(ensg_df.loc[(ensg_df.source == 'ensembl_havana')].transcript.unique())
    ensg_df.drop(ensg_df.loc[(ensg_df.source == 'havana') & (ensg_df.transcript.isin(ID_list))].index, axis=0, inplace=True)
    ID_list = ID_list.append((pd.Series(ensg_df.loc[(ensg_df.source == 'havana')].transcript.unique())))
    ensg_df.drop(ensg_df.loc[(ensg_df.source == 'ensembl') & (ensg_df.transcript.isin(ID_list))].index, axis=0, inplace=True)
    ID_list = ID_list.append((pd.Series(ensg_df.loc[(ensg_df.source == 'ensembl')].transcript.unique())))

    # rank values by start site in ascending order
    ensg_df = ensg_df.sort_values(by=['start'], ascending=True)

    # append FASTA sequences to the ensg_df - also appends range, not sure this is necessary
    ENST_FASTA_dict = {}
    range_dict = {}
    for ID in ensg_df.transcript.unique():
        FASTA_list = []
        range_list = []
        for row in ensg_df.loc[ensg_df.transcript == ID].itertuples():
            FASTA_list.append(row.FASTA)
            range_list.append(range((row.start), (row.stop + 1)))
        ENST_FASTA_dict.update({row.transcript:FASTA_list})
        range_dict.update({row.transcript:range_list})

    # create a range_df with transcript, start, stop, stranf, FASTA, and exon junctions ('FASTA_range' column)
    range_df = ensg_df[['transcript', 'start', 'stop', 'strand']].drop_duplicates(subset = 'transcript', keep = 'first')
    range_df.reset_index(drop=True, inplace=True)
    range_df['FASTA'] = range_df.transcript.map(ENST_FASTA_dict)
    range_df.loc[:, 'FASTA'] = range_df.apply(lambda row : ''.join(row.FASTA), axis = 1)
    range_df['FASTA_range'] = range_df.transcript.map(range_dict)

    # must write output as a .pkl
    rangepath = os.path.join(outdir, 'range_df.pkl')
    range_df.to_pickle( rangepath )
    print('Completed', rangepath)




def create_TITER_df(gff3_path, FASTA_path, outdir):

    # get first 100 nucleotides of CDS to use in a calculation with titer_df
    # outputs 'titer_CDS.tsv' in the outdir

    # import ensembl .gff3 file
    ensg_names = ['seqid', 'source', 'type', 'start', 'stop', 'score', 'strand', 'phase', 'attributes']
    ensg_dtypes = {'seqid':str, 'source':str, 'type':str, 'start':int, 'stop':int, 'score':str, 'strand':str, 'phase':str, 'attributes':str}
    ensg_df = pd.read_csv(gff3_path, header=0, sep='\t', comment='#', names=ensg_names, dtype = ensg_dtypes)

    # create a list of all 23 chromosomes
    chrom_list = [str(i) for i in (list(range(1,23)))]
    chrom_list.append('X')
    chrom_list.append('Y')

    # remove non-chrom entries from seqid column
    ensg_df = ensg_df.loc[ensg_df.seqid.isin(chrom_list)]

    # subset to coding sequence
    CDS_df = ensg_df.loc[ensg_df.type == 'CDS']

    # define transcript (without version number)
    CDS_df['transcript'] = CDS_df.attributes.str.split(';').str[1].str.split(':').str[1]

    # cut out duplicate entries in df, preferentially keeping entries in ensembl_havana, then havana, then ensembl (df agreement > manual annotation > Havana alone)
    ID_list = pd.Series(CDS_df.loc[(CDS_df.source == 'ensembl_havana')].transcript.unique())
    CDS_df.drop(CDS_df.loc[(CDS_df.source == 'havana') & (CDS_df.transcript.isin(ID_list))].index, axis=0, inplace=True)
    ID_list = ID_list.append((pd.Series(CDS_df.loc[(CDS_df.source == 'havana')].transcript.unique())))
    CDS_df.drop(CDS_df.loc[(CDS_df.source == 'ensembl') & (CDS_df.transcript.isin(ID_list))].index, axis=0, inplace=True)
    ID_list = ID_list.append((pd.Series(CDS_df.loc[(CDS_df.source == 'ensembl')].transcript.unique())))

    # use pybedtools to get the UTR FASTA sequence
    BED_df = CDS_df[['seqid','start','stop']]
    BED_df.rename(columns = {'seqid':'chrom', 'start':'chromStart', 'stop':'chromEnd'}, inplace=True)
    BED_df.loc[:, 'chrom'] = 'chr'+ BED_df.chrom
    BED_df.loc[:, 'chromStart'] = BED_df.chromStart - 1 # adjusts to zero-based BEDtools nucleotide assignment
    
    # get FASTA sequences
    print('Obtaining FASTA sequences for TITER analysis.')
    FASTA_list = get_seq(BED_df, FASTA_path)
        
    CDS_df['FASTA'] = np.asarray(FASTA_list) # read FASTA elements back to the pos_df column

    # read all regions like they are positive start codons, and then take the reverse complement for negative strands
    CDS_df = CDS_df.sort_values(by=['start'], ascending=True)

    # append FASTA sequences to the ensg_df
    ENST_CDS_dict = {}
    for ID in ID_list:
        FASTA_list = []
        for row in CDS_df.loc[CDS_df.transcript == ID].itertuples():
            FASTA_list.append(row.FASTA)
        CDS_sequence = ''.join(FASTA_list)
        ENST_CDS_dict.update({row.transcript:CDS_sequence})
    CDS_df['FASTA'] = CDS_df.transcript.map(ENST_CDS_dict)

    # convert negative strand FASTA to reverse complement
    CDS_df.loc[CDS_df.strand == '-', 'FASTA'] = CDS_df.loc[CDS_df.strand == '-'].FASTA.apply(complement_function).str[::-1]

    # subset FASTA sequence for all entries to start codon + last 100 nucleotides for FASTA sequence
    CDS_df['titer_str'] = CDS_df.FASTA.str[0:103]

    # slim dataframe and write to output directory
    CDS_df = CDS_df[['transcript', 'titer_str']]
    titerpath = os.path.join(outdir, 'TITER_CDS.tsv')
    CDS_df.to_csv(os.path.join(titerpath), sep = '\t', index=False)
    print('Completed:', titerpath)
    
    return titerpath

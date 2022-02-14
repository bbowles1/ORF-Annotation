# ORF-Annotation
An annotation pipeline to identify disruptive variation in 5'UTR regions of the human genome.

# Background
This builds upon work by McGillivray et al. to programmatically identify all predicted upstream open reading frames in the human genome: 
McGillivray P, Ault R, Pawashe M, Kitchen R, Balasubramanian S, Gerstein M. A comprehensive catalog of predicted functional upstream open reading frames in humans. Nucleic Acids Res. 2018 Apr 20;46(7):3326-3338. doi: 10.1093/nar/gky188. PMID: 29562350; PMCID: PMC6283423.
http://github.gersteinlab.org/uORFs/

We have used the hg19 reference build for the all ORF-Annotation work.

# Catalog Setup
To begin this work, the user will need to edit paths in Catalog_Setup_Config.py. We have provided a path to matched transcripts between Ensembl and RefSeq, but the user will need to provide their own paths to the following files:
 - A download of the McGillivray et al catalog (available in the manuscript supplementary data).
 - A path to an Ensembl .gff3 file (we used hg19 version 87).
 - A path to an hg19 FASTA file.
 - A path of an Ensembl biomart download containing Ensembl transcripts and their matched RefSeq equivalents. We have provided the "GRCh37_biomart_refseq_export.txt" file as an example.

The catalog setup config will create the following files:
 - A .tsv file containing McGillivray et al annotations and a .bed file for the same regions.
 - A range_df.tsv file, which contains Ensembl 5'UTR exon regions and their associated FASTA sequence.


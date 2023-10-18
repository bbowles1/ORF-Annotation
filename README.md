# ORFA Variant Interpreter
The Open Reading Frame Annotation (ORFA) variant interpreter assesses a 5'UTR variant's impact on upstream open reading frame (uORF) presence. It is intended to take a VCF file and label each variant with all possible effects on overlapping transcripts (it produces a one-to-many output where a single variant may have different effect on multiple overlapping transcripts). It labels all variants with a categorical prediction of effect on uORF regions, such as "new AUG" or "AUG removal."

# Dependencies
- Bedtools
- Tabix
- Pandas (1.4+ recommended)
- Pybedtools (0.8+ recommended)

# Background
This builds upon work by McGillivray et al. to programmatically identify all predicted upstream open reading frames in the human genome: 
McGillivray P, Ault R, Pawashe M, Kitchen R, Balasubramanian S, Gerstein M. A comprehensive catalog of predicted functional upstream open reading frames in humans. Nucleic Acids Res. 2018 Apr 20;46(7):3326-3338. doi: 10.1093/nar/gky188. PMID: 29562350; PMCID: PMC6283423.
http://github.gersteinlab.org/uORFs/

We have used the GRCh37/hg19 reference build for the all ORF-Annotation work.

# Required Input Files
Running all ORFA scripts requires several input files:
 - A download of the McGillivray et al. catalog (available in the manuscript supplementary data as "complete_uORF_predictions_hg19.tsv").
 - A path to an Ensembl .gff3 file (we used hg19 version 87).
 - A path to a hg19 FASTA file (build/FASTA_chrom_identifiers.txt has been included to handle RefSeq chromosome naming conventions in FASTA reference files).
 - A path to an Ensembl Biomart download containing Ensembl transcripts and their matched RefSeq equivalents. We have provided the "build/GRCh37_biomart_refseq_export.txt" file as an example.

# Building the Catalog
Users will need to run a first-time setup to generate several additional reference files.
1. To begin this work, the user will need to download the above files and provide their absolute paths in the "config.txt" file.
2. After updating config.txt, the user will need to run `bin/Catalog_Setup_Config.py`. This will create the following files:
 - uorf_catalog.tsv: Contains a filtered set of McGillivray et al. annotations and their start locations within the gneome.
 - uorf_catalog.bed: A .bed file for the same regions as above.
 - Range_df.pkl: A data structure mapping Ensembl 5'UTR exon regions to their associated FASTA sequence.
 - TITER_CDS.tsv: A data structure mapping transcript IDs to a 200bp sequenced use for input in the TITER algorithm to predict ribosome laod (Zhang et al. 2017, https://doi.org/10.1093/bioinformatics/btx247).

# Variant Annotation
We will be working to provide our variant labelling and interpreter scripts as we prepare this work for publication. Please check back soon.


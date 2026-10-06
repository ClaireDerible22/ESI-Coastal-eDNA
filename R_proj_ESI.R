# Code to merge ASV tables and taxonomy from eDNA data processed in Qiime2, using dada2 and BLAST or RDP Classifier for taxonomy.
##Using RDP Classifier for Taxonomy 
###Trying to read in Taxonomy and ASV files
asv_data <- read.delim("C:/Users/clair/Documents/GitHub/ESI-Coastal-eDNA/data/2023/COI-LERAYXT/ESI2023_featuretable_export.tsv")
taxrdp <- read.delim("C:/Users/clair/Documents/GitHub/ESI-Coastal-eDNA/data/2023/COI-LERAYXT/rdp.output")

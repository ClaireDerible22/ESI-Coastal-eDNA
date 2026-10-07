# Code to merge ASV tables and taxonomy from eDNA data processed in Qiime2, using dada2 and BLAST or RDP Classifier for taxonomy.
##Using RDP Classifier for Taxonomy 
###Trying to read in Taxonomy and ASV files
asv_data <- read.delim("C:/Users/clair/Documents/GitHub/ESI-Coastal-eDNA/data/2023/COI-LERAYXT/ESI2023_featuretable_export.tsv")
taxrdp <- read.delim("C:/Users/clair/Documents/GitHub/ESI-Coastal-eDNA/data/2023/COI-LERAYXT/rdp.output")
#testing to see if R project is syncing with GitHub by hitting save and then checking desktop GitHub

###Now we can begin
# Load libraries ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(tidyverse)

#simple function to filter reads in an ASV file - set to 0.01% now but can be changed to 0.1% etc
filter_low_reads <- function(asv_data) {
  numeric_asv_data <-asv_data[sapply(asv_data,is.numeric)]
  total_reads <- sum(numeric_asv_data, na.rm = TRUE)  # Total reads across species and sites
  species_sums <- rowSums(numeric_asv_data, na.rm = TRUE)  # Sum of reads per species
  asv_data_filtered <- asv_data[species_sums >= 0.0001 * total_reads, ]  # Keep species with at least 0.01% of total reads
  return(asv_data_filtered)
}

#at this point, I can't actually see the new filtered data, I am going to try and save the function as a new dataset
filtered_asv_data <- filter_low_reads(asv_data)

##ESI
### ---Merging data---
#ASV
esi12s<-read.table(file = "data/2021Data/NEW/12S/ESI2021_12S_feature_table_export.tsv", header = T, sep = "\t")
#taxonomy
esi12s.taxa <-read.table(file = "data/2021Data/NEW/12S/ESI_12Sblast_1results.tsv",header = F,sep="\t")
colnames(esi12s.taxa)<-c("ASV","NCBI","percentID", "evalue","length","species","group","commonname")

esi21.12s.merge <- left_join(esi12s, esi12s.taxa, by =c("OTU.ID"="ASV"))  %>% filter(group %in% c("bony fishes","whales & dolphins", "sharks & rays", "birds"), percentID > 90)


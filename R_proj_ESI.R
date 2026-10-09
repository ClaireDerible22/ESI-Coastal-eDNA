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
##ESI COI seining data



#Now the 2023 COI seining data (ASV data and rdp.ouput as taxonomy)
esi23.coi.asv <-read.table("filtered_asv_data.tsv, row.names=FALSE", header = T, sep="\t") %>% glimpse()
esi23.coi.taxa <-read.table("C:/Users/clair/Documents/GitHub/ESI-Coastal-eDNA/data/2023/COI-LERAYXT/rdp.output", header = F, sep="\t") %>% select(c("V1","V12","V15", "V24","V26","V27","V29"))

#MERGE

esi23.coi.merge <-left_join(esi23.coi.asv,esi23.coi.taxa, by=c("OTU.ID"="V1"))

#filter
esi23.coi.filt <- esi23.coi.merge %>% filter(V26>0.94, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida","Mollusca","Nematoda","Rhodophyta","Gastrotricha","Chlorophyta","Echinodermata","Brachiopoda","Porifera","Cnidaria","Nemertea","Haptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies","Tardigrada","Rotifera", "Chaetognatha")) %>% select(!starts_with(c("ENEG","EXT","PCRB"))) %>%
  rename(Phylum=V12, Class=V15, Species=V27, Confidence=V29)

esi23.coi.filt2 <- filter_low_reads(esi23.coi.filt) %>% filter(!Species %in% c("Sus_scrofa","Homo_sapiens")) %>%
  select(-c(OTU.ID,Phylum,Class,V24,V26)) %>%
  relocate(Species) %>%
  relocate(Confidence, .after=Species)

##write csv

write.csv(x = esi23.coi.filt2, file = "C:/Users/clair/Documents/GitHub/ESI-coastal-eDNA/data/2023/COI-LERAYXT/ESI2023_Coastal_COI_GOTeDNA.csv", quote = F, row.names = F)


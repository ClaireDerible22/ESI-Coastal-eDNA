# Code to merge ASV tables and taxonomy from eDNA data processed in Qiime2, using dada2 and BLAST or RDP Classifier for taxonomy.
##Using RDP Classifier for Taxonomy 
###Trying to read in Taxonomy and ASV files
asv_data <- read.delim("C:/Users/clair/Documents/GitHub/ESI-Coastal-eDNA/data/2023/COI-LERAYXT/ESI2023_featuretable_export.tsv")
taxrdp <- read.delim("C:/Users/clair/Documents/GitHub/ESI-Coastal-eDNA/data/2023/COI-LERAYXT/rdp.output")
#testing to see if R project is syncing with GitHub by hitting save and then checking desktop GitHub

###Now we can begin (Working with 2023 ESI Seining COI data)
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


#####Filtering and Merging 2023 ESI 12S data################################################################

#Imported both data sets required (ASV and taxonomy)

##We're going to skip the filtering by low reads step for now...

###Trying to read in (and rename for R use) 12S Taxonomy and ASV files ##Header = false so that first data column isn't deleted (sep=, is used for csv to keep variables from merging into one)
esi23.12S.asv <- read.csv("data/2023/12S/Seining_ASV_TaxonTable_Filtered_FishOnly.csv", header = T, sep=",")
esi23.12S.taxa<- read.delim("data/2023/12S/seining23_blast.tsv", header = F)


###MERGE###

esi23.12S.merge <-left_join(esi23.12S.asv,esi23.12S.taxa, by=c("ID"="V1"))

###filter 12S for fish, dolphin, whales, sharks

esi23.12S.filt.merge <- esi23.12S.merge %>% filter(PercentID>97.99 & V7>97.99 %in% c("bony fishes","whales & dolphins","sharks & rays"))

################################### MAKING MAP #######################################################

##USING SCRIPT OUTLINE FROM NICK##


#Make a map of eDNA sites in the ESI
install.packages("sf")
install.packages("raster")
install.packages("ggplot2")
install.packages("dplyr")
library(sf)
library(raster)
library(ggplot2)
library(dplyr)

##reading in master data file with coordinates (saved as a csv file)

sites <- read.csv("C:\\Users\\clair\\Desktop\\HONOURS\\SeiningMasterData2023_mapping.csv")
##following two lines filters out only the site, long, and lat columns
sites_clean <- sites%>%
  select(Site, Lat, Long)

## removing Wolfes Island samples
cleaned_sites <- sites_clean%>%
  filter(Site != "WolfesIsland")


target_sites <- c("TaylorHead", "Moosehead", "Goldboro", "Little Harbour", "Conrod Beach")

cleaned_sites <- cleaned_sites %>% 
  filter(Site %in% target_sites) %>%
  group_by(Site) %>%
  slice_sample(n = 1) %>%
  ungroup()




latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"


novsco <-read_sf("data/mapping_shapefiles/EasternShoreIslands_networksite.shp") %>% 
  st_transform(latlong) %>% mutate(name="Nova Scotia")%>% dplyr::select(name, geometry)

esi <-read_sf("data/mapping_shapefiles/NS_coastline_project_Erase1.shp")%>%
  st_transform(latlong)%>%
  mutate(name="ESI")%>%
  dplyr::select(name,geometry)

#perley.sites <-read.csv("~/GitHub/easternshoreislands_aoi/data/SamplingDesign_65by20km.csv", header = T) %>% glimpse()

focal_bound <- esi%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_transform(utm)%>%
  st_buffer(4)%>% 
  st_transform(latlong)%>%
  st_bbox()

ggplot()+
  geom_sf(data=novsco, fill="red", alpha=0.5 )+
  geom_sf(data=esi, fill="grey40", linewidth=1.1, colour="grey40")+
  coord_sf(expand=0, xlim=focal_bound[c(1,3)], ylim=focal_bound[c(2,4)])+
  xlab(label = "Longitude")+
  ylab(label="Latitude")+
  geom_point(data=cleaned_sites, aes(x=Long, y=Lat, color=Site, shape=Site), size=4)+
  scale_color_manual(name="Site",values = c("Moosehead"="blue4", "Conrod Beach"="yellow", "Goldboro"="green", "Little Harbour"="blue3", "TaylorHead"="blue"))+
  scale_shape_manual(name="Site", values= c("Moosehead"=15, "Conrod Beach"=16, "Goldboro"=16, "Little Harbour"=17, "TaylorHead"=18))+
  #geom_point(data=perley.sites, aes(x=lon, y=lat), shape=21, colour='black', fill="cornflowerblue", size=4)+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 14, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14),
        text=element_text(size=18))


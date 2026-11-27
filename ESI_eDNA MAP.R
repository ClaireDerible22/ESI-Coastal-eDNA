################################### MAKING MAP #######################################################

##USING SCRIPT OUTLINE FROM NICK##
## From Nick - I would make this map script a separate script Claire!

#Make a map of eDNA sites in the ESI
#install.packages("sf")
#install.packages("raster")
#install.packages("ggplot2")
#install.packages("dplyr")
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

##this is to rename a specific cell (TaylorHead to Taylor Head)
cleaned_sites$Site[5] <- "Taylor Head"


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
  scale_color_manual(name="Site",values = c("Moosehead"="blue4", "Conrod Beach"="green3", "Goldboro"="green3", "Little Harbour"="blue4", "Taylor Head"="blue4"))+
  scale_shape_manual(name="Site", values= c("Moosehead"=18, "Conrod Beach"=17, "Goldboro"=16, "Little Harbour"=20, "Taylor Head"=15))+
  #geom_point(data=perley.sites, aes(x=lon, y=lat), shape=21, colour='black', fill="cornflowerblue", size=4)+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 14, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14),
        text=element_text(size=18))

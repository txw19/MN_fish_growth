library(ggplot2)
library(dplyr)
library(reshape)

library(mapproj)
all_states <- map_data("state")
states <- subset(all_states, region %in% c(  "minnesota") )

#read in location
all.lakes<-read.csv("D:/Users/grhansen/Documents/MN niche modeling/Data/lakelist.csv", header=TRUE, sep=",")

#count species from gillnets and trapnets
GN<- read.csv("~/MN fish trends/Data/3.10.2014 GN CPUE DATA_BEST.csv", header=T)
TN <- read.csv("~/MN fish trends/Data/3.10.2014 TN CPUE DATA_BEST.csv", header=T)

#simplest way - list all species ever caught by lake, count unique species.
GN.spp=select(GN, DOW, Spp)
TN.spp=select(TN, DOW, Spp)
both.spp=rbind(GN.spp, TN.spp)
Spp.unique=unique(both.spp)
spp.counts=summarise(group_by(Spp.unique, DOW), Species=length(Spp))

spp.counts=merge(spp.counts, all.lakes, by="DOW")

windows()
p <- ggplot()
p <- p + geom_polygon( data=states, aes(x=long, y=lat, group = group),colour="black", fill="grey90" )
p=p+geom_point( data=spp.counts, aes(LAKE_CENTER_LONG_DD5, LAKE_CENTER_LAT_DD5,  colour=Species), size=3)
p=p+scale_colour_distiller(name="Spp. Count", palette="YlOrRd", trans = "reverse")
p=p+theme(
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.background = element_blank(),
  axis.text.x= element_text(size=12),
  axis.text.y=element_text(size=12),
  legend.position=c(.75,.45) , legend.key = element_blank(),
  #legend.title=element_text(size=16, hjust=-0.3),legend.text=element_text(size=16),
  axis.title.y=element_blank(),
  axis.title.x=element_blank() 
)
p

ggsave(file="GN_TN_species_count_MN_map.png", width=6, height=6, units="in", dpi=300)


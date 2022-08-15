#script 4/4/20
#B subtype paper

#####################
#Seq Prevalence Maps#
#####################

library(sqldf)
library(dplyr)
library(ggplot2)
library(smacpod)
library(spdep)
library(rsatscan)
library(rgdal)
library(spatstat)
library(rgeos)
library(maptools)
library(RColorBrewer)
library(DCluster)
library(sp)
library(gstat)
library(tigris)
library(dismo)
library(reshape2)
library(Hmisc)
library(INLA) #not available
library(readxl)
library(ggthemes)
library(stringr)

#import data
phylo <- read.csv("insert data file")
head(phylo)
countyfips<-read.csv("insert data file",stringsAsFactors = T,header = T)
countyfips$fips<-sub(".", "", countyfips$fips)
phylo$county_name <- as.character(phylo$county_name)
phylo$county = substr(phylo$county_name,1,nchar(phylo$county_name)-4)


#load datasets for clustered and singleton seqs, clean up counties
clusters_only$county_name <- as.character(clusters_only$county_name)
clusters_only$county = substr(clusters_only$county_name,1,nchar(clusters_only$county_name)-4)
singletons_only$county_name <- as.character(singletons_only$county_name)
singletons_only$county = substr(singletons_only$county_name,1,nchar(singletons_only$county_name)-4)

phylo$county[phylo$county == "*UNK" | phylo$county == ""] <- NA
phylo$county[phylo$county == "*DE SOTO"] <- "DESOTO"
phylo$county[phylo$county == "ST LUCIE"] <- "ST. LUCIE"

clusters_only$county[clusters_only$county == "*UNK"] <- ""
clusters_only$county[clusters_only$county == "*DE SOTO"] <- "DESOTO"
clusters_only$county[clusters_only$county == "ST LUCIE"] <- "ST. LUCIE"
singletons_only$county[singletons_only$county == "*UNK"] <- ""
singletons_only$county[singletons_only$county == "*DE SOTO"] <- "DESOTO"
singletons_only$county[singletons_only$county == "ST LUCIE"] <- "ST. LUCIE"

#sanity checks
table(phylo$county)
table(countyfips$county)
table(clusters_only$county)
table(singletons_only$county)
table(countyfips$fips)

#subset for seqs after 2012, subtype B
table(phylo$subtype)
phylo_sub <- subset(phylo, genotype_yr >=2012 & subtype == "B")

phylo_fp<-merge(phylo_sub,countyfips,by="county",all.x = T,all.y =T)
phylo_fp<-merge(clusters_only,countyfips,by="county",all.x = T,all.y=T)
phylo_fp<-merge(singletons_only,countyfips,by="county",all.x = T,all.y=T)

phylo_sub <- clusters_only

phylo_sub <- phylo_fp


#seqs by county
q<-"
SELECT fips, county, count(pat_idx) as n_seq
FROM phylo_sub
GROUP BY fips
;
"

cty_seq<-sqldf(q)

#remove obs with cty unknown
cty_seq_notNA <- subset(cty_seq, fips != "")
cty_seq_notNA <- cty_seq

#categorize
fivenum(cty_seq_notNA$n_seq)
cty_seq_notNA$seq_cat[cty_seq_notNA$n_seq < 10] <- "<10"
cty_seq_notNA$seq_cat[cty_seq_notNA$n_seq >= 10 & cty_seq_notNA$n_seq <100] <- "10-99"
cty_seq_notNA$seq_cat[cty_seq_notNA$n_seq >= 100 & cty_seq_notNA$n_seq <1000] <- "100-999"
cty_seq_notNA$seq_cat[cty_seq_notNA$n_seq >= 1000] <- "1000+"
table(cty_seq_notNA$seq_cat)

#download FL shp data from Census
countyfl2010<-readOGR(dsn="insert data file",encoding="ESRI Shapefile")
countyfl2010@data$fips <- countyfl2010@data$COUNTYFP10
head(countyfl2010@data)
table(countyfl2010@data$fips)

#merge
temp<-countyfl2010@data
temp$id2<-c(1:nrow(temp))
temp<-merge(temp,cty_seq_notNA,by="fips",all.x=T)
#temp$seq_cat[is.na(temp$seq_cat)] <- "<10"
temp<-temp[order(temp$id2),]
table(temp$fips, temp$seq_cat)
countyfl2010@data<-temp

#create map of number of seqs available for the entire period 2012-2017
#gpclibPermit()
countyfl2010_f<-fortify(countyfl2010,region="fips")
countyfl2010_f_data<-countyfl2010@data
countyfl2010_f<-merge(countyfl2010_f,countyfl2010_f_data,by.x="id",by.y="fips",all.x=T)
countyfl2010_f$x<-countyfl2010_f$long
countyfl2010_f$y<-countyfl2010_f$lat


table(countyfl2010_f_data$fips, countyfl2010_f_data$seq_cat)

basemap<-ggplot(countyfl2010_f)+geom_polygon()+geom_path(color="white")+coord_equal()+theme(panel.background =element_blank(),axis.title=element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust=0.5,face="bold",size=16))

(mapseq<-basemap+scale_fill_brewer("Sequences",type="seq",palette = "PuRd")+aes(long,lat,group=group, fill=seq_cat)+labs(title="Frequency of HIV-1 pol Subtype B Sequences by FL County, 2012-2017") + geom_path(data=countyfl2010_f,size=1,aes(long,lat,group=group)))

(mapclusseq<-basemap+scale_fill_brewer("Sequences",type="seq",palette = "PuRd")+aes(long,lat,group=group, fill=seq_cat)+labs(title="Frequency of Clustered HIV-1 pol Subtype B Sequences by FL County, 2012-2017") + geom_path(data=countyfl2010_f,size=1,aes(long,lat,group=group)))

(mapsingleseq<-basemap+scale_fill_brewer("Sequences",type="seq",palette = "PuRd")+aes(long,lat,group=group, fill=seq_cat)+labs(title="Frequency of Unclustered HIV-1 pol Subtype B Sequences by FL County, 2012-2017") + geom_path(data=countyfl2010_f,size=1,aes(long,lat,group=group)))


ggsave("insert data path/mapclusseq.png",dpi=300)

################
#Prevalence map#
################
#create map of number of people living with HIV for the entire period 2012-2017
prev_by_county <- read.csv("insert data file")

prev_by_county$Cases<-as.numeric(sub(",", "", prev_by_county$Cases,fixed = TRUE))
fivenum(prev_by_county$Cases)
prev_by_county$prev_cat[prev_by_county$Cases < 100] <- "<100"
prev_by_county$prev_cat[prev_by_county$Cases >= 100 & prev_by_county$Cases <1000] <- "100-999"
prev_by_county$prev_cat[prev_by_county$Cases >= 1000 & prev_by_county$Cases <10000] <- "1000-9999"
prev_by_county$prev_cat[prev_by_county$Cases >= 10000] <- "10000+"
table(prev_by_county$prev_cat)

#merge
temp<-countyfl2010@data
temp$id2<-c(1:nrow(temp))
temp<-merge(temp,prev_by_county,by.x="county",by.y="County",all.x=T)
#temp$seq_cat[is.na(temp$seq_cat)] <- "<10"
temp<-temp[order(temp$id2),]
table(temp$fips, temp$prev_cat)
countyfl2010@data<-temp

countyfl2010_f<-fortify(countyfl2010,region="fips")
countyfl2010_f_data<-countyfl2010@data
countyfl2010_f<-merge(countyfl2010_f,countyfl2010_f_data,by.x="id",by.y="fips",all.x=T)
countyfl2010_f$x<-countyfl2010_f$long
countyfl2010_f$y<-countyfl2010_f$lat

basemap<-ggplot(countyfl2010_f)+geom_polygon()+geom_path(color="white")+coord_equal()+theme(panel.background =element_blank(),axis.title=element_blank(),axis.text = element_blank(),axis.ticks = element_blank(),plot.title = element_text(hjust=0.5,face="bold",size=16))

#saved as TIFF 850x600
(mapprev<-basemap+scale_fill_brewer("Cases",type="seq",palette = "PuRd")+aes(long,lat,group=group, fill=prev_cat)+labs(title="Prevalence of HIV by FL County, 2018") + geom_path(data=countyfl2010_f,size=1,aes(long,lat,group=group)))

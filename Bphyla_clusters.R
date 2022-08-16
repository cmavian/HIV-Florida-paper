#Analysis of B subtype clusters and singletons in Florida
#By: Shannan Rich

#load libraries
library("readxl")
library("sqldf")
library("dplyr")
library("car")
library("countrycode")
library("vegan")
library("RColorBrewer")

#import data
phylo_data <- read.csv("INSERT FILE")


#######################
#import super clusters#
#######################
cs2to4 <- read.csv("INSERT FILE")
cs5to10 <- read.csv("INSERT FILE")
cs11to30 <- read.csv("INSERT FILE")
cs31to73 <- read.csv("INSERT FILE")

#bind
clusteredseqs <- rbind(cs2to4, cs5to10, cs11to30, cs31to73)

#add column for cluster_yn
clusteredseqs$cluster_yn <- 1

#trim categorize phylo_data for only necessary data fields
phylo_data1 <- sqldf("SELECT geno_idx, pat_idx, genotype_yr, county_name, subtype, birth_sex,
                     dob_yr, expo_categ, race, birth_country_cd, hiv_dx_yr, dod_yr FROM phylo_data")

#create new categories
attach(phylo_data1)
phylo_data1$age_cat[age_dx <= 14] <- "0to14 years"
phylo_data1$age_cat[age_dx >= 15 & age_dx <= 19] <- "15to19 years"
phylo_data1$age_cat[age_dx >= 20 & age_dx <= 25] <- "20to25 years"
phylo_data1$age_cat[age_dx >= 26 & age_dx <= 35] <- "26to35 years"
phylo_data1$age_cat[age_dx >= 36 & age_dx <= 45] <- "36to45 years"
phylo_data1$age_cat[age_dx >= 46 & age_dx <= 55] <- "46to55 years"
phylo_data1$age_cat[age_dx >= 56] <- "56+ years"
table(phylo_data1$age_cat)
phylo_data1$expo <- recode(phylo_data1$expo_categ, "c(9, 10, 11) = 'Unknown';
                                                    c(03) = 'Heterosexual';
                                                    c(08) = 'Perinatal';
                                                    c(01, 06) = 'MSM';
                                                    c(02, 04, 05, 07) = 'IDU'")
table(phylo_data1$expo)
phylo_data1$race_eth <- recode(phylo_data1$race, "c(1) = 'Hispanic';
                                                  c(4) = 'Black';
                                                  c(6) = 'White';
                                                  c(2,3,5,7, 8) = 'Other'")
table(phylo_data1$race_eth)

#create primary dataset
fullphylaset <- sqldf("SELECT phylo_data1.*, clusteredseqs.cluster_size FROM phylo_data1 
              LEFT JOIN clusteredseqs USING(geno_idx)")

#subset for sequences from 2012-2017
fullphylaset_sub <- subset(fullphylaset, genotype_yr >= 2012)

#subset for B sequences only
Bphylaset <- fullphylaset[which(fullphylaset$subtype == 'B'),]

#recode country of birth
table(Bphylaset$birth_country_cd)
Bphylaset$birthregion <- countrycode(Bphylaset$birth_country_cd, origin="iso3c", destination="region")
table(Bphylaset$birthregion)
Bphylaset$region_CAT <- recode(Bphylaset$birthregion, "c('Australia and New Zealand', 'Central Asia', 'Eastern Asia', 'Melanesia', 'Micronesia', 'Polynesia', 
                                                             'South-Eastern Asia', 'Southern Asia', 'Western Asia') = 'AsiaPacific';
                                                       c('Central America', 'South America') = 'Latin America';
                                                       c('Eastern Africa', 'Middle Africa', 'Northern Africa', 'Southern Africa', 'Western Africa') = 'Africa';
                                                       c('Eastern Europe', 'Northern Europe', 'Southern Europe', 'Western Europe') = 'Europe'")
table(Bphylaset$region_CAT)


#export cleaned dataset
write.csv(Bphylaset, "INSERT PATH\\Bphylaset.csv")


############################################
#import dataset, compute descriptive stats!#
############################################
Bphylaset <- read.csv("INSERT PATH\\Bphylaset.csv", stringsAsFactors = FALSE)
attach(Bphylaset)


##################
#additional fixes#
##################
#add singleton option
Bphylaset$cluster_size[is.na(Bphylaset$cluster_size)] <- "Singleton"

#recode age
Bphylaset$new_age <- Bphylaset$age_cat
Bphylaset$new_age[age_cat == '0-14 years']  <- "0-19 years"
Bphylaset$new_age[age_cat == '15-19 years'] <- "0-19 years"

#recode sex
Bphylaset$birth_sex[birth_sex == 'M'] <- 'Male'
Bphylaset$birth_sex[birth_sex == 'F'] <- 'Female'

#recode region
Bphylaset$new_region <- Bphylaset$region_CAT
Bphylaset$new_region[region_CAT == 'Africa'] <- 'Foreign'
Bphylaset$new_region[region_CAT == 'AsiaPacific'] <- 'Foreign'
Bphylaset$new_region[region_CAT == 'Caribbean'] <- 'Foreign'
Bphylaset$new_region[region_CAT == 'Europe'] <- 'Foreign'
Bphylaset$new_region[region_CAT == 'Latin America'] <- 'Foreign'

#recode expo group
Bphylaset$expo_new <- Bphylaset$expo
Bphylaset$expo_new[expo == 'Perinatal'] <- 'Other'
Bphylaset$expo_new[expo == 'Unknown'] <- 'Other'

#recode county
Bphylaset$county_new <- Bphylaset$county_name
Bphylaset$county_new[Bphylaset$county_name == "*DE SOTO CO."] <- "DESOTO CO."
Bphylaset$county_new[Bphylaset$county_name == "ST LUCIE CO."] <- "ST. LUCIE CO."
Bphylaset$county_new[Bphylaset$county_name == "*UNKNOWN"] <- NA
length(unique(Bphylaset$county_new))

#create districts (based on FLGISA URL = https://www.flgisa.org/districts/)
Bphylaset$district <- Bphylaset$county_new

Bphylaset$district <- recode(Bphylaset$county_new, "
                      c('ESCAMBIA CO.','SANTA ROSA CO.','OKALOOSA CO.','WALTON CO.','HOLMES CO.',
                        'WASHINGTON CO.','BAY CO.','JACKSON CO.','CALHOUN CO.','GULF CO.',
                        'GADSDEN CO.','LIBERTY CO.','FRANKLIN CO.','LEON CO.','WAKULLA CO.',
                        'JEFFERSON CO.') = 'Northwest';
                      c('MADISON CO.','TAYLOR CO.','HAMILTON CO.','SUWANNEE CO.','LAFAYETTE CO.',
                        'DIXIE CO.','COLUMBIA CO.','GILCHRIST CO.','LEVY CO.','BAKER CO.','UNION CO.',
                        'ALACHUA CO.','BRADFORD CO.','MARION CO.','NASSAU CO.','DUVAL CO.','CLAY CO.',
                        'PUTNAM CO.','ST. JOHNS CO.') = 'Northeast';
                      c('CITRUS CO.','HERNANDO CO.','PASCO CO.','HILLSBOROUGH CO.','PINELLAS CO.',
                        'SUMTER CO.') = 'Central East';
                      c('LAKE CO.','POLK CO.','FLAGLER CO.','VOLUSIA CO.','SEMINOLE CO.','ORANGE CO.',
                        'OSCEOLA CO.','BREVARD CO.','INDIAN RIVER CO.') = 'Central West';
                      c('MANATEE CO.','SARASOTA CO.','CHARLOTTE CO.','LEE CO.','COLLIER CO.',
                        'HARDEE CO.','DESOTO CO.','HIGHLANDS CO.','GLADES CO.','HENDRY CO.') = 'Southwest';
                      c('OKEECHOBEE CO.','ST. LUCIE CO.','MARTIN CO.','PALM BEACH CO.','BROWARD CO.',
                        'MIAMI-DADE CO.','MONROE CO.') = 'Southeast'
                             ")
                   
#create urban/rural categ based on US Census 2010
Bphylaset$rural_urban <- "Urban"
Bphylaset$rural_urban <- recode(Bphylaset$county_new, "
                                c('WALTON CO.','HOLMES CO.','WASHINGTON CO.','JACKSON CO.','CALHOUN CO.',
                                  'GULF CO.','LIBERTY CO.','FRANKLIN CO.','GADSDEN CO.','WAKULLA CO.',
                                  'JEFFERSON CO.','MADISON CO.','TAYLOR CO.','HAMILTON CO.','SUWANNEE CO.',
                                  'LAFAYETTE CO.','DIXIE CO.','COLUMBIA CO.','GILCHRIST CO.','LEVY CO.',
                                  'BAKER CO.','UNION CO.','BRADFORD CO.','HARDEE CO.','DESOTO CO.',
                                  'HIGHLANDS CO.','GLADES CO.','HENDRY CO.','MONROE CO.','OKEECHOBEE CO.') = 'Rural'
                                ")
Bphylaset$rural_urban[Bphylaset$rural_urban != 'Rural'] <- 'Urban'
Bphylaset$rural_urban[Bphylaset$county_new == ''] <- NA

#create HIV dx year categ
Bphylaset$hiv_dx_yr_NEW <- Bphylaset$hiv_dx_yr
Bphylaset$hiv_dx_yr_NEW[Bphylaset$hiv_dx_yr < 2010] = '<2010'
Bphylaset$hiv_dx_yr_NEW[Bphylaset$hiv_dx_yr >= 2010 & Bphylaset$hiv_dx_yr <= 2011] <- '2010-11'
Bphylaset$hiv_dx_yr_NEW[Bphylaset$hiv_dx_yr >= 2012 & Bphylaset$hiv_dx_yr <= 2013] <- '2012-13'
Bphylaset$hiv_dx_yr_NEW[Bphylaset$hiv_dx_yr >= 2014 & Bphylaset$hiv_dx_yr <= 2015] <- '2014-15'
Bphylaset$hiv_dx_yr_NEW[Bphylaset$hiv_dx_yr >= 2016 & Bphylaset$hiv_dx_yr <= 2017] <- '2016-17'
table(Bphylaset$hiv_dx_yr_NEW)


#remove people with missing county (n=96)
Bphylaset1 <- subset(Bphylaset, county_new != '')
#remove people with genotypes before 2012
Bphylaset2 <- subset(Bphylaset1, genotype_yr >= 2012)

#specify EHE counties
Bphylaset2$EHEcounty<-0
Bphylaset2$EHEcounty[Bphylaset2$county_new %in% c("MIAMI-DADE CO.","BROWARD CO.","PINELLAS CO.","DUVAL CO.","HILLSBOROUGH CO.","PALM BEACH CO.","ORANGE CO.")]<-1
table(Bphylaset2$EHEcounty)

clusters_only <- subset(Bphylaset2, cluster_size != "Singleton")

#relabel cluster size
clusters_only$cluster_sizeFIX[clusters_only$cluster_size == "2to4_"] <- "a) 2-4"
clusters_only$cluster_sizeFIX[clusters_only$cluster_size == "5to10_"] <- "b) 5-10"
clusters_only$cluster_sizeFIX[clusters_only$cluster_size == "11to30_"] <- "c) 11-28"
clusters_only$cluster_sizeFIX[clusters_only$cluster_size == "31to73_"] <- "d) 29-70"
#changed to 2-4, 5-10, 11-28, 29-70 since data were restricted to 2012-2017

singletons_only <- subset(Bphylaset2, cluster_size == "Singleton") 

attach(clusters_only)

#tables
stat_age <- table(new_age, cluster_sizeFIX)
stat_sex <- table(birth_sex, cluster_sizeFIX)
stat_race_eth <- table(race_eth, cluster_sizeFIX)
stat_district <- table(district, cluster_sizeFIX)
stat_rural_urban <- table(rural_urban, cluster_sizeFIX)
stat_region <- table(new_region, cluster_sizeFIX)
stat_expo <- table(expo_new, cluster_sizeFIX)
stat_hiv_dx_yr <- table(hiv_dx_yr_NEW, cluster_sizeFIX)
stat_genotype_yr <- table(genotype_yr, cluster_sizeFIX)

##################
#Chi-square tests#
##################
stat_age
age.per<-round(prop.table(stat_age, 2)*100, digits = 1)
chisq.test(stat_age) 
barplot(stat_age, legend=TRUE,col=brewer.pal(6,name = "Spectral"))
barplot(age.per, legend=TRUE,col=brewer.pal(6,name = "Spectral"),xlab="Cluster size (# sequences)",ylab="Percent",main="Cluster Size by Age Category")

stat_sex
sex.per<-round(prop.table(stat_sex, 2)*100, digits = 1)
chisq.test(stat_sex) 
barplot(stat_sex, legend=TRUE,col=brewer.pal(2,name = "Spectral"))
barplot(sex.per, legend=TRUE,col=brewer.pal(2,name = "Spectral"))

stat_race_eth
race_eth.per<-round(prop.table(stat_race_eth, 2)*100, digits = 1)
chisq.test(stat_race_eth)
barplot(stat_race_eth, legend=TRUE,col=brewer.pal(4,name = "Spectral"))
barplot(race_eth.per, legend=TRUE,col=brewer.pal(4,name = "Spectral"))

stat_district
district.per<-round(prop.table(stat_district, 2)*100, digits = 1)
chisq.test(stat_district) 
barplot(stat_district, legend=TRUE,col=brewer.pal(6,name = "Spectral"))
barplot(district.per, legend=TRUE,col=brewer.pal(6,name = "Spectral"))

stat_rural_urban
round(prop.table(stat_rural_urban, 2)*100, digits = 1)
chisq.test(stat_rural_urban) 

stat_region
reg.per<-round(prop.table(stat_region, 2)*100, digits = 1)
chisq.test(stat_region) 
barplot(stat_region, legend=TRUE,col=brewer.pal(2,name = "Spectral"))
barplot(reg.per, legend=TRUE,col=brewer.pal(2,name = "Spectral"))

stat_expo
risk.per<-round(prop.table(stat_expo, 2)*100, digits = 1)
chisq.test(stat_expo) 
barplot(stat_expo, legend=TRUE,col=brewer.pal(4,name = "Spectral"))
barplot(risk.per, legend=TRUE,col=brewer.pal(4,name = "Spectral"))

stat_hiv_dx_yr
dx.per<-round(prop.table(stat_hiv_dx_yr, 2)*100, digits = 1)
chisq.test(stat_hiv_dx_yr)
barplot(stat_hiv_dx_yr, legend=TRUE,col=brewer.pal(5,name = "Spectral"))
barplot(dx.per, legend=TRUE,col=brewer.pal(5,name = "Spectral"))

stat_genotype_yr
round(prop.table(stat_genotype_yr, 2)*100, digits = 1)
chisq.test(stat_genotype_yr) 

#calculate genotype year mean & IQR by cluster_size & perform kruskal.test for 1 IV with 2+ groups and an ordinal DV 
aggregate(genotype_yr, list(cluster_size), median)
aggregate(genotype_yr, list(cluster_size), fivenum)
aggregate(genotype_yr, list(cluster_size), sd)
aggregate(genotype_yr, list(cluster_size), mean)

kruskal.test(genotype_yr, cluster_size)

aggregate(age_dx, list(cluster_size), mean)
aggregate(age_dx, list(cluster_size), sd)
aggregate(age_dx, list(cluster_size), fivenum)
kruskal.test(age_dx, cluster_size)


#####################################################################
#multiplot of significant findings, with column for unclustered seqs#
#####################################################################
Bphylaset2$cluster_sizeFIX<- "e) Unclustered"
Bphylaset2$cluster_sizeFIX[Bphylaset2$cluster_size == "2to4_"] <- "a) 2-4"
Bphylaset2$cluster_sizeFIX[Bphylaset2$cluster_size == "5to10_"] <- "b) 5-10"
Bphylaset2$cluster_sizeFIX[Bphylaset2$cluster_size == "11to30_"] <- "c) 11-28"
Bphylaset2$cluster_sizeFIX[Bphylaset2$cluster_size == "31to73_"] <- "d) 29-70"
table(Bphylaset2$cluster_sizeFIX)

attach(Bphylaset2)

#tables
stat_age <- table(new_age, cluster_sizeFIX)
stat_sex <- table(birth_sex, cluster_sizeFIX)
stat_race_eth <- table(race_eth, cluster_sizeFIX)
stat_district <- table(district, cluster_sizeFIX)
stat_expo <- table(expo_new, cluster_sizeFIX)
stat_hiv_dx_yr <- table(hiv_dx_yr_NEW, cluster_sizeFIX)

stat_age
age.per<-round(prop.table(stat_age, 2)*100, digits = 1)

stat_sex
sex.per<-round(prop.table(stat_sex, 2)*100, digits = 1)

stat_race_eth
race_eth.per<-round(prop.table(stat_race_eth, 2)*100, digits = 1)

stat_district
district.per<-round(prop.table(stat_district, 2)*100, digits = 1)

stat_expo
risk.per<-round(prop.table(stat_expo, 2)*100, digits = 1)

stat_hiv_dx_yr
dx.per<-round(prop.table(stat_hiv_dx_yr, 2)*100, digits = 1)

library(wesanderson) #color palettes
windowsFonts(times = windowsFont("Times New Roman")) 
par(family = "times",mfrow=c(3,2),mar=c(4, 5, 3, 9))
barplot(age.per, legend=TRUE,col=wes_palette("Zissou1", n = 6,type = "continuous"),ylab="Percent",main="Age",args.legend = list(x ='right', bty='n', inset=c(-0.25,0)))
barplot(sex.per, legend=TRUE,col=wes_palette("Zissou1", n = 2,type = "continuous"),main="Sex",args.legend = list(x ='right', bty='n', inset=c(-0.20,0)))
barplot(race_eth.per, legend=TRUE,col=wes_palette("Zissou1", n = 4,type = "continuous"),main="Race and Ethnicity",ylab="Percent",args.legend = list(x ='right', bty='n', inset=c(-0.2,0)))
barplot(district.per, legend=TRUE,col=wes_palette("Zissou1", n = 6,type = "continuous"),main="District",args.legend = list(x ='right', bty='n', inset=c(-0.3,0)))
barplot(risk.per, legend=TRUE,col=wes_palette("Zissou1", n = 4,type = "continuous"),xlab="Cluster size (# sequences)",ylab="Percent",main="Transmission Category",args.legend = list(x ='right', bty='n', inset=c(-0.3,0)))
barplot(dx.per, legend=TRUE,col=wes_palette("Zissou1", n = 5,type = "continuous"),xlab="Cluster size (# sequences)",main="Diagnosis Year",args.legend = list(x ='right', bty='n', inset=c(-0.2,0)))


###############################
#binary: cluster vs no cluster#
###############################
#Bphylaset2 <- Bphylaset1
Bphylaset2$cluster_status <- Bphylaset2$cluster_size
Bphylaset2$cluster_status[cluster_size == '2to4_'] <- 'Clustered'
Bphylaset2$cluster_status[cluster_size == '5to10_'] <- 'Clustered'
Bphylaset2$cluster_status[cluster_size == '11to30_'] <- 'Clustered'
Bphylaset2$cluster_status[cluster_size == '31to73_'] <- 'Clustered'

attach(Bphylaset2)

#tables
stat_age <- table(age_cat, cluster_status)
stat_sex <- table(birth_sex, cluster_status)
stat_race_eth <- table(race_eth, cluster_status)
stat_district <- table(district, cluster_status)
stat_rural_urban <- table(rural_urban, cluster_status)
stat_region <- table(region_CAT, cluster_status)
stat_expo <- table(expo, cluster_status)
stat_hiv_dx_yr <- table(hiv_dx_yr_NEW, cluster_status)
stat_genotype_yr <- table(genotype_yr, cluster_status)

##################
#Chi-square tests#
##################
stat_age
round(prop.table(stat_age, 2)*100, digits = 1)
chisq.test(stat_age)

stat_sex
round(prop.table(stat_sex, 2)*100, digits = 1)
chisq.test(stat_sex)

stat_race_eth
round(prop.table(stat_race_eth, 2)*100, digits = 1)
chisq.test(stat_race_eth)

stat_district
round(prop.table(stat_district, 2)*100, digits = 1)
chisq.test(stat_district)

stat_rural_urban
round(prop.table(stat_rural_urban, 2)*100, digits = 1)
chisq.test(stat_rural_urban)

stat_region
round(prop.table(stat_region, 2)*100, digits = 1)
chisq.test(stat_region)

stat_expo
round(prop.table(stat_expo, 2)*100, digits = 1)
chisq.test(stat_expo)

stat_hiv_dx_yr
round(prop.table(stat_hiv_dx_yr, 2)*100, digits = 1)
chisq.test(stat_hiv_dx_yr)

#calculate genotype year mean & IQR by cluster_status & perform wilcox 1 IV with 2 groups 
aggregate(genotype_yr, list(cluster_status), median)
aggregate(genotype_yr, list(cluster_status), fivenum)
wilcox.test(genotype_yr~cluster_status)

aggregate(age_dx, list(cluster_status), mean)
aggregate(age_dx, list(cluster_status), sd)
aggregate(age_dx, list(cluster_status), fivenum)

wilcox.test(age_dx~cluster_status)

#####################
#multivariable model#
#####################
#set clustered = 1
Bphylaset2$cluster_recode <- 0
Bphylaset2$cluster_recode[Bphylaset2$cluster_status == "Clustered"] <- 1

#change reference groups
Bphylaset2$region_CAT <- as.factor(Bphylaset2$region_CAT)
Bphylaset3 <- within(Bphylaset2, region_CAT <- relevel(region_CAT, ref = 'Northern America'))
Bphylaset3$hiv_dx_yr_NEW <- as.factor(Bphylaset3$hiv_dx_yr_NEW)
Bphylaset4 <- within(Bphylaset3, hiv_dx_yr_NEW <- relevel(hiv_dx_yr_NEW, ref = '2016-17'))

#run glm 
fullmodel <- glm(cluster_recode ~ age_dx + birth_sex + race_eth + district + rural_urban + region_CAT + expo + hiv_dx_yr_NEW + genotype_yr, 
                 family=binomial(link="logit"), data = Bphylaset4)
round(exp(cbind(coef(fullmodel), confint(fullmodel))), digits = 2)
summary(fullmodel)

#run boosted glm 
library("mboost")
Bphyla.gb <- glmboost(cluster_recode ~ age_dx + birth_sex + race_eth + district + rural_urban + region_CAT + expo + hiv_dx_yr_NEW + genotype_yr, data = Bphylaset4,
                    control = boost_control(mstop = 2000), center = FALSE)
exp(coef(Bphyla.gb), confint(Bphyla.gb))
confint.glmboost(Bphyla.gb)




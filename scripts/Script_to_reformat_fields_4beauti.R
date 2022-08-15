##################################################################
# Script to put seq location & time in correct format for BEAUTI #
##################################################################
#load libs
library("bios2mds")
library("tidyr")

#import fasta file with bios2mds; extract names into dataframe;  
dat <- import.fasta("insert data file", aa.to.upper = TRUE, gap.to.dash = TRUE)
dat_names <- as.data.frame(names(dat))
dat_names$name <- dat_names$`names(dat)`

#separate seq names into new columns with tidyr
dat_names1 <- separate(dat_names, name, c("id","expo","region","country","year"), sep = "_", remove = TRUE,
         convert = FALSE, extra = "warn", fill = "warn")

#reformat location: retain NA-FL for FDOH sequences, set rest to region
dat_names1$location = dat_names1$region
dat_names1$location[dat_names1$country == "UNITED-STATES-FL"] <- "UNITED-STATES-FL"

#concatenate name & location/year with '=' deliminiter 
#location
paste(dat_names1$name, dat_names1$location, sep = "=", collapse = ",")

#year
paste(dat_names1$name, dat_names1$year, sep = "=", collapse = ",")


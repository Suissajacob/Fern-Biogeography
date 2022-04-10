# This is the scrip used to clean the GBIF data, change the taxonomy, and create a working dataset for downstream analyses

#read packages
library(tidyverse)
library(CoordinateCleaner)
library(maps)
library(rbokeh)
library("widgetframe")
library(taxize)




####################
# maket he world plot for the rest of the script
####################
world <- map_data("world") 

ga <- ggplot()
gg <- ga + geom_map(data=world, map=world, aes(x=long, y=lat, map_id=region),color="white", fill="#7f7f7f", size=0.05, alpha=1/4)


####################
#read in data
####################
gbifRawFern.data<- read.csv("data", header = TRUE)



####################
# check number of types of observations
####################
head(gbifRawFern.data, n = 1)

recordType.code<- unique(gbifRawFern.data$taxonRank)


####################
# take the data and make the factors numeric and remove all 0, 0 coordinates and anything below -60 lat
# filter the data so that you remove all decimal lat and long that are 0 and below -60
####################
step1_cleaned.data<- gbifRawFern.data%>%
  dplyr::select(order, family, genus, species,scientificName,taxonRank, decimalLatitude, decimalLongitude, elevation, basisOfRecord, issue,  year, month, countryCode, recordNumber, recordedBy)%>%
  dplyr::filter(taxonRank== "SPECIES"| taxonRank== "SUBSPECIES" | taxonRank=="VARIETY", .preserve = T)%>% 
  dplyr::mutate(decimalLatitude = as.numeric(as.character(decimalLatitude)))%>%
  dplyr::mutate(decimalLongitude = as.numeric(as.character(decimalLongitude)))%>%
  dplyr::mutate(elevation = as.numeric(as.character(elevation)))%>%
  dplyr::filter(!decimalLongitude == 0 & !decimalLatitude == 0 & !decimalLatitude < -60, .preserve = T)

####################
# check number of distinct names
####################
check_dist_names<- step1_cleaned.data%>%
  distinct(step1_cleaned.data$species)  



###################
#Make a new column with just the variety
###################

step2_cleaned.data<-step1_cleaned.data %>% 
  separate(col =scientificName, sep = "var\\.", into = c(NA, "var"), remove = FALSE)%>%
  separate(col =scientificName, sep = "subsp.", into = c(NA, "subsp"), remove = FALSE)%>%
  separate(col =var,  sep = " ", into = c(NA, "var"))%>%
  separate(col =subsp,  sep = " ", into = c(NA, "subsp."))

gbif_test<- gbif_test%>%
  arrange(taxonRank)


####################
# check to see how many and what type of flaggs come out with the data 
# clean the dataset with clean coordinates
####################
cord_cleaner_output <- CleanCoordinates(step2_cleaned.data, lat = "decimalLatitude", lon = "decimalLongitude", species = step1_cleaned.data$species, inst.rad = 1000, centroids.rad = 5000, capitals.rad=5000)

####################
# check out the bad cords
####################
bad_cords.fern<- cord_cleaner_output%>%
  filter(.summary==FALSE, .preserve = T)

####################
# plot the bad cords
####################

gg+ geom_point(data= bad_cords.fern, aes(decimalLongitude, decimalLatitude), size=1, alpha = 0.5) +xlim(-200, 200)+ ggtitle("GBIF occurrence records removed by Coordinate Cleaner: 89216")


####################
#filter out all bad cords just keeping good ones
####################
fernCleanCords.data<- cord_cleaner_output%>%
  filter(!.summary==FALSE, .preserve = T)%>%
  select(order, family, genus, species,scientificName, subsp., var, taxonRank,decimalLatitude, decimalLongitude, elevation,month, year,countryCode, recordNumber, recordedBy)

####################
# check out how many species and genera  we have left
####################
dist <- fernCleanCords.data%>%
  distinct(species)%>%
  arrange(species) 

dist.genus <- fernCleanCords.data%>%
  distinct(genus)%>%
  arrange(genus) 

####################
# save the data
####################
setwd("path")

write.csv(fernCleanCords.data, "CleanCords_pre_taxize22Apr2020.csv")


###################
#split the periods and make the author names the same as that in the Hassler data
##################

dat<- read.csv("path to: CleanCords_pre_taxize22Apr2020.csv")

##################
# delete spaces and make the auth the same as hassler
##################

#Step 1 remove all spaces after periods to standardize everything
dat$GBIF_scientificName <- gsub("\\. ", ".", dat$GBIF_scientificName)

# Step 2 add 1 space after all periods
dat$GBIF_scientificName <- gsub("\\.", "\\. ", dat$GBIF_scientificName)

# Step 3 remove all spaces in between . and )
dat$GBIF_scientificName <- gsub("\\. \\)", "\\.\\)", dat$GBIF_scientificName)

# Step 4 remove all spaces at the end of the last word
dat$GBIF_scientificName <- gsub(" $","", dat$GBIF_scientificName)

# change subsp. to ssp. to corroborate the hassler list
dat$GBIF_scientificName <- gsub("subsp.", "ssp.", dat$GBIF_scientificName)


setwd("path")

write.csv(dat, "CleanCords_pre_taxize23Apr2020.csv")



################
#Get tautonyms
################

################
# seperate the gbif species column by spaces
#################
dat<-dat %>% 
  separate(col =GBIF_scientificName, sep = " ", into = c("genus", "spp"), remove = FALSE)

##################
#Find out how many tautonyms there are by asking how many columns spp = var/species column
##################

# filter out only the rows where the spp and var or spp. and subsp column are equal i.e. only the tautonyms 
tautCheck<- dat%>%
  filter(spp==var | spp==subsp., .preserve = TRUE)%>%
  mutate(GBIF_species=as.character(GBIF_species))

##############
#Make distinct dataset of tauts
###############

taut_distinct<-tautCheck%>%
  distinct(GBIF_scientificName,.keep_all = TRUE)

# save
setwd("path")

write.csv(taut_distinct, "tautolonyms_distinct.csv")


###############
#read in data
###############

OCC_dat<- read.csv("path to cleaned occ dat", header = TRUE, stringsAsFactors = FALSE)



##########
# filter out all taxa with NO_MATCH in the correct_GBIF_NAME column
#########

OCC_filtered<-OCC_dat%>%
  mutate(ACCEPTED_Names=FIXED_GBIF_names)%>%
  filter(!ACCEPTED_Names=="NO_MATCH", .preserve = TRUE)%>%
  select(!c(FIXED_GBIF_names))

# now with removing NO_MATCH we have 853950

#########
#Check unique taxa
#########
unique_SPECIES<- OCC_filtered%>%
  distinct(ACCEPTED_Names)

#########
#Check unique genera
#########
unique_genus<- OCC_filtered%>%
  distinct(genus)

setwd("path")

write.csv(OCC_filtered, "Final_Synon_TOT_OCC_27Apr2020.csv" )

###############
#Find most of the centroids and remove them  
###############


##############
#After all points that should be deleted are determined from the interactive maps, make a list of them and remove them
##############


##############
#Remove centroids
##############

OCC_filtered<- read.csv("path to clean occ ", header = TRUE, stringsAsFactors = FALSE)

##############
# all centroids of bolivia. This was done for many countries, cities, etc.
##############
centroid_bolivia<- OCC_filtered%>%
  filter(decimalLatitude== -17.77722, .preserve = TRUE)%>%
  filter(decimalLongitude== -64.72638, .preserve = TRUE)%>%
  select(genus, X)




##########
#merge the data points from the centroids to remove
##########
centroid_points_to_remove <- bind_rows(input all of the centroids you fond)

points.list<- as.integer(centroid_points_to_remove$X)

##########
#remove all the promblematic centroid points from the dataframe
##########
OCC_centroid_removed<- OCC_filtered%>%
  filter(., !X %in% points.list, .preserve = TRUE)

##########
#make a new column called Accepted_genus
##########

cleaned_occ_dat3May<-OCC_centroid_removed %>% 
  separate(col =ACCEPTED_Names, sep = " ", into = c("Accepted_genus", NA), remove = FALSE)



# Save the new cleaned data
setwd("path")
write.csv(cleaned_occ_dat3May, "Cleaned_OCC_tobeMapped_3.0_3May_2020.csv")



#########################################################################################################################################    Make genus maps to find the problem points     ######################################


##############
# make interactive maps of the genera and remove problematic points
##############

#Load the packages
library(tidyverse)
library(rgbif)
library(maps)
library(CoordinateCleaner)
library(RColorBrewer)
library(rbokeh)
library("widgetframe")

#################################
# read in the data
#################################

Fern_div.df<- read.csv("input the path to the data")


###################################################################
# Make the maps
###################################################################
setwd("path")

genus_list<- length(Fern_div.df$genus)


for(i in 1:length(genus_list)){
  
  genus <- levels(Fern_div.df$genus)# this takes the genus names and makes them a character list
  
  i=1 # this will make the genus of interest the value for the map to plot
  
  a <- Fern_div.df[which(Fern_div.df$genus==genus[i]),] # subset the dataframe to just have the data from the genus of interest
  
  plot <- suppressWarnings(figure(width = 800, height = 450, padding_factor = 0)%>%         ly_map("world", col = "gray") %>%
                             ly_points(decimalLongitude, decimalLatitude, data = a, size = 5,hover = c(species, decimalLongitude, decimalLatitude, X)))# make the plot
  
  htmlwidgets::saveWidget(as_widget(plot), "index.html")# this is how you save it as an html
  
  #widgetframe::frameWidget(plot,width=600,height=400)# make it interactive # this would be how you plot it in R
}




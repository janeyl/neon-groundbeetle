#__________________________________________
# library----
#__________________________________________

library(tidyverse)
library(car)       #Anova (not anova)
library(RColorBrewer)
library(neonUtilities)
library(vegan)
library(dplyr)
library(plotrix) #std.error function


#clean work area
rm(list=ls(all=TRUE))

# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

#__________________________________________
# Read in relevant raw data sets----
#__________________________________________

PlantHARV <-loadByProduct(dpID="DP1.10098.001", site="HARV", 
                          package="expanded", check.size=T)

BeetleHARV<- loadByProduct(dpID="DP1.10022.001", site="HARV", 
                           package="expanded", check.size=T) 

PlantBART <-loadByProduct(dpID="DP1.10098.001", site="BART", 
                          package="expanded", check.size=T)

BeetleBART<- loadByProduct(dpID="DP1.10022.001", site="BART", 
                           package="expanded", check.size=T) 

PlantTALL <-loadByProduct(dpID="DP1.10098.001", site="TALL", 
                          package="expanded", check.size=T)

BeetleTALL<- loadByProduct(dpID="DP1.10022.001", site="TALL", 
                           package="expanded", check.size=T) 

BeetleJERC<- loadByProduct(dpID="DP1.10022.001", site="JERC", 
                           package="expanded", check.size=T) 

Mycorrhizal.df<-read.csv ("plantSpecies.csv", header=TRUE)

#__________________________________________
# Create new tidy data sets HARV----
#__________________________________________
if(TRUE){
  Mycorrhizal.df<-Mycorrhizal.df%>%                                            #dataset with leaf habit and fungal association by taxonID
    select(taxonID,                                                            # download from:
           treeType, 
           mycorrhizal)                                                        #remove scientific name from taxonID
  
  plots_HARV.df<-BeetleHARV$bet_fielddata%>%
    distinct(plotID, .keep_all = FALSE)                                        #beetle pitfall trap sampling locations
  
  mappingandtagging_HARV.df<-PlantHARV$vst_mappingandtagging%>%
    filter(plotID %in% plots_HARV.df$plotID)%>%                                #mapped and tagged trees with the species info
    select(eventID,                                                            #filter by the plots that have beetle samples
           individualID,                                                       #remove unwanted columns from df
           taxonID, 
           scientificName)
  
  perplotperyear_HARV.df<-PlantHARV$vst_perplotperyear%>%
    filter(plotID %in% plots_HARV.df$plotID)%>%                                #filter by the plots that have beetle samples
    filter(eventID=="vst_HARV_2019")%>%                                        #CHECK EACH SITE: retain only latest fully sampled year. Harv:2019
    select(plotID,                                                             #Checking to see the plot area sampled
           nlcdClass, 
           totalSampledAreaTrees)                                              #Harv: consistent 400 square m per site
  
  apparentindividual_HARV.df<-PlantHARV$vst_apparentindividual%>%
    filter(plotID %in% plots_HARV.df$plotID)%>%
    filter(eventID=="vst_HARV_2019")%>%                                        #retain only year 2019 for Harv
    filter(growthForm=="multi-bole tree"|growthForm=="single bole tree")%>%    #keep only single & multi bole trees    
    filter(!stemDiameter < 10)%>%                                              #filter out trees <10 cm in Diameter
    filter(plantStatus=="Live"|
             plantStatus=="Live, insect damaged"|
             plantStatus=="Live,  other damage"|
             plantStatus=="Live, disease damaged"|
             plantStatus=="Live, physically damaged"|
             plantStatus=="Live, broken bole"|
             plantStatus=="other damage")%>%
    select(plotID,
           individualID,
           stemDiameter)%>%
    left_join(., perplotperyear_HARV.df, by = "plotID")%>%                     #join in sampleing area
    left_join(., mappingandtagging_HARV.df, by = "individualID")%>%            #add species ID and names
    left_join(., Mycorrhizal.df, by = "taxonID")%>%                           #fungal associations and leaf habit
    mutate(taxonID=dplyr::recode(taxonID, 'ACSA3' = 'ACSA2'),                         #check for same species
           taxonID=dplyr::recode(taxonID, 'ACSAS' = 'ACSA2'), 
           scientificName=dplyr::recode(scientificName, 'Acer saccharum Marshall' = 'Acer saccharinum L.'),
           scientificName=dplyr::recode(scientificName, 'Acer saccharum Marshall var. saccharum' = 'Acer saccharinum L.'),
           basalAreaCm2=(3.142*((stemDiameter/2)^2)))                          #calculate the basal area of each tree from the stem diameter
}
#__________________________________________
# Create new data frame based on basal area ----
#__________________________________________
if(TRUE){
  BasalAreaHARV.df<-apparentindividual_HARV.df%>%                              #make a df with calculations based on basal area
    spread(treeType, basalAreaCm2)%>%                                          #making a column for leaf habits, and the basal area of each tree to category
    mutate(basalAreaCm2=(3.142*((stemDiameter/2)^2)))%>%         
    spread(mycorrhizal, basalAreaCm2)%>%                                       #making a column for leaf habits, and the basal area of each tree to category
    mutate(Evergreen = replace_na(Evergreen, 0),                               #replace na's with 0 
           Deciduous = replace_na(Deciduous, 0), 
           ECM = replace_na(ECM, 0), 
           AM = replace_na(AM, 0), 
           totalBA=Deciduous+Evergreen)%>%                                     #calculate the totalBA for each plot
    group_by(plotID)%>%                                                        #calculate the total area of deciduous, evergreen, am, ecm, species at the plot level
    summarise(Deciduous = sum(Deciduous),
              Evergreen = sum(Evergreen),
              AM = sum(AM),
              ECM = sum(ECM),
              totalBA = sum(totalBA))%>%                                       #total species density at each plot by basal area of stems
    mutate(perEvergreen=Evergreen/totalBA*100,                                 #calculate the percent evergreen species at each plot
           perECM=ECM/totalBA*100)%>%                                          #calculate the percent ECM at each plot
    left_join(., apparentindividual_HARV.df, by = "plotID")%>%                 #re-join taxonID to create a diversity calculation
    select(plotID,  
           nlcdClass,
           taxonID,
           basalAreaCm2,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalBA,
           perEvergreen,
           perECM)%>%
    group_by(plotID)%>%                                                        #designate basal area of each species into columns, grouped by plotID
    pivot_wider(names_from = "taxonID", values_from = "basalAreaCm2", values_fn = list(basalAreaCm2 = sum))%>%
    replace(is.na(.),0)
  species.df<-BasalAreaHARV.df[,10:20]                                         #CHECK EACH SITE:for vegan package, subset species into columns for diversity
  BasalAreaHARV.df$ShannonPlant<-diversity(species.df, index = "shannon")      #Shannon Diversity Index
  BasalAreaHARV.df$SimpsonPlant<-diversity(species.df, index = "simpson")      #Simpson Diversity Index
  BasalAreaHARV.df<-BasalAreaHARV.df%>%         
    select(plotID,                                                             #retain necessary columns
           nlcdClass,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalBA,
           perEvergreen,
           perECM,
           ShannonPlant,
           SimpsonPlant)
  rm(species.df)                                                               #remove the species df
}

#__________________________________________
# Create new data frame based on number of stems----
#__________________________________________
if(TRUE){
  StemsHARV.df<-apparentindividual_HARV.df%>%                                  #make a df with calculations based on stems of each species
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%                   #count each stem 
    spread(treeType, stemCount)%>%                                             #count by leaf habit
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%
    spread(mycorrhizal, stemCount)%>%                                          #count by fungal association
    mutate(Evergreen = replace_na(Evergreen, 0),                               #replace NA's
           Deciduous = replace_na(Deciduous, 0),
           ECM = replace_na(ECM, 0),
           AM = replace_na(AM, 0),
           totalStems=Deciduous+Evergreen)%>%                                  #total count of each species
    group_by(plotID)%>%                                                        #summarize at the plot level
    summarise(Deciduous = sum(Deciduous),
              Evergreen = sum(Evergreen),
              AM = sum(AM),
              ECM = sum(ECM),
              totalStems = sum(totalStems))%>%
    mutate(perEvergreen=Evergreen/totalStems*100,                              #make a percent evergreen calculation
           perECM=ECM/totalStems*100)%>%                                       #make a percent ECM calculation
    left_join(., apparentindividual_HARV.df, by = "plotID")%>%                 #join in species data
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%
    select(plotID,                                                             #retain necessary columns
           nlcdClass,
           taxonID,
           stemCount,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalStems,
           perECM,
           perEvergreen)%>%
    group_by(plotID)%>%                                                        #regroup by plotID and place each species into columns
    pivot_wider(names_from = "taxonID", values_from = "stemCount", values_fn = list(stemCount = sum))%>%
    replace(is.na(.),0)
  speciesStems.df<-StemsHARV.df[,10:20]                                        #CHECK SITE: Harv has 10 species, so subset columns 10-20
  StemsHARV.df$ShannonPlant<-diversity(speciesStems.df, index = "shannon")     #Shannon Diversity Index
  StemsHARV.df$SimpsonPlant<-diversity(speciesStems.df, index = "simpson")     #Simpson Diversity Index
  StemsHARV.df<-StemsHARV.df%>%                                         
    select(plotID,                                                             #retain only necessary columns
           nlcdClass,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalStems,
           perEvergreen,
           perECM,
           ShannonPlant,
           SimpsonPlant)
  rm(speciesStems.df)
}
#__________________________________________
# Descriptive Statistics SITE level----
#__________________________________________
if(TRUE){
  SiteStatsHARV.df<-apparentindividual_HARV.df%>%                              #look at abundance of each species at HARV by stem and basal area
    select(.,plotID,
           taxonID,
           scientificName,
           basalAreaCm2)%>%
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%                   #create a species count
    group_by(taxonID)%>%                                                       #group by species ID
    count(taxonID, taxonID, name = "stems")%>%                                 #makes a column for the number of stems 
    left_join(., apparentindividual_HARV.df, by = "taxonID")%>%                #rejoin to add basal area column
    select(.,
           taxonID,
           stems,
           scientificName,
           basalAreaCm2)%>%
    group_by(taxonID, stems, scientificName)%>%
    summarise(basalAreaCm2 = (sum(basalAreaCm2)))                              #sum basal area by species ID
  SiteStatsHARV.df<-SiteStatsHARV.df%>%
    mutate(relAbunBA = basalAreaCm2/sum(SiteStatsHARV.df$basalAreaCm2)*100,    #calculate the relative adundance of a species at HARV
           relAbunSt = stems/sum(SiteStatsHARV.df$stems)*100)                  #relative abundance = 1 species/ total species at site
}
#__________________________________________
# Descriptive Statistics PLOT level----
#__________________________________________
if(TRUE){
  PlotStatsHARV.df<-BasalAreaHARV.df%>%                                        #create a df with descriptive statists at the plot level
    select(plotID,                                                             #select parameters based on basal area
           perECM,
           perEvergreen,
           totalBA)%>%
    left_join(., StemsHARV.df, by = "plotID")%>%
    select(., plotID,                                                          #select parameters based on number of stems
           totalStems,
           totalBA,
           perECM.x,
           perEvergreen.x,
           perECM.y,
           perEvergreen.y
    )
  PlotStatsHARV.df<-PlotStatsHARV.df%>%                                        #calculate MEAN basal area per plot
    summarise(.,                                                               #calculate STD ERROR OF MEAN basal area per plot
              #calculate MIN & MAX basal area per plot
              avgBA = mean(PlotStatsHARV.df$totalBA), 
              seBA = std.error(PlotStatsHARV.df$totalBA),
              minBA = min(PlotStatsHARV.df$totalBA),
              maxBA = max(PlotStatsHARV.df$totalBA),
              sumtotalBA = sum(PlotStatsHARV.df$totalBA),
              
              avgperECMba = mean(PlotStatsHARV.df$perECM.x),                   #calculate MEAN percent ECM fungal association based on basal area per plot
              seperECMba = std.error(PlotStatsHARV.df$perECM.x),               #calculate STD ERROR OF MEAN percent ECM fungal association based on basal area per plot
              minperECMba = min(PlotStatsHARV.df$perECM.x),                    #calculate MIN & MAX percent ECM fungal association based on basal area per plot
              maxperECMba = max(PlotStatsHARV.df$perECM.x),
              
              avgperEGba = mean(PlotStatsHARV.df$perEvergreen.x),              #calculate MEAN percent Evergreen leaf habib by basal area per plot
              seperEGba = std.error(PlotStatsHARV.df$perEvergreen.x),          #calculate STD ERROR OF MEAN percent Evergreen leaf habib by basal area per plot
              minperEGba = min(PlotStatsHARV.df$perEvergreen.x),               #calculate MIN & MAX percent Evergreen leaf habib by basal area per plot
              maxperEGba = max(PlotStatsHARV.df$perEvergreen.x),
              
              avgStems = mean(PlotStatsHARV.df$totalStems),                    #MEAN total species by stem and per plot
              seStems = std.error(PlotStatsHARV.df$totalStems),                #STD ERROR OF MEAN by total species by stem and per plot
              minStems = min(PlotStatsHARV.df$totalStems),                     #MIN & MAX total species by stem and per plot
              maxStems = max(PlotStatsHARV.df$totalStems),
              sumtotalStems = sum(PlotStatsHARV.df$totalStems),
              
              avgperECMstems = mean(PlotStatsHARV.df$perECM.y),                #MEAN percent ECM fungal association by stem and per plot
              seperECMstems = std.error(PlotStatsHARV.df$perECM.y),            #STD ERROR OF MEAN percent ECM fungal association by stem and per plot
              minperECMstems = min(PlotStatsHARV.df$perECM.y),                 #MIC & MAX percent ECM fungal association by stem and per plot
              maxperECMstems = max(PlotStatsHARV.df$perECM.y),
              
              avgperEGstems = mean(PlotStatsHARV.df$perEvergreen.y),           #MEAN percent evergreen leaf habit by stem and per plot
              seperEGstems = std.error(PlotStatsHARV.df$perEvergreen.y),       #STD ERROR OF MEAN percent evergreen leaf habit by stem and per plot
              minperEGstms = min(PlotStatsHARV.df$perEvergreen.y),             #MIN & MAX percent evergreen leaf habit by stem and per plot
              maxperEGstms = max(PlotStatsHARV.df$perEvergreen.y),
              
              avgShannonBA = mean(BasalAreaHARV.df$ShannonPlant),
              seShannonBA = std.error(BasalAreaHARV.df$ShannonPlant),
              avgSimpsonBA = mean(BasalAreaHARV.df$SimpsonPlant),
              seSimpsonBA = std.error(BasalAreaHARV.df$SimpsonPlant),
              
              avgShannonStems = mean(StemsHARV.df$ShannonPlant),
              seShannonStems = std.error(StemsHARV.df$ShannonPlant),
              avgSimpsonStems = mean(StemsHARV.df$SimpsonPlant),
              seSimpsonStems = std.error(StemsHARV.df$SimpsonPlant))%>%
    select(!plotID)%>%                                                  #remove plotID to make a final df
    distinct()%>%
    pivot_longer(., everything(), 
                 names_to = "HARV", 
                 values_to = "Value")
}
#__________________________________________
# Create new tidy data sets TALL----
#__________________________________________
if(TRUE){
  plots_TALL.df<-BeetleTALL$bet_fielddata%>%
    distinct(plotID, .keep_all = FALSE)            #beetle pitfall trap sampling locations
  # PLOTS: 1-9, 13
  
  mappingandtagging_TALL.df<-PlantTALL$vst_mappingandtagging%>%
    filter(plotID %in% plots_TALL.df$plotID)%>%    #mapped and tagged trees with the species info, filter by beetle plots sampled
    select(individualID,                           #remove unwanted columns from df
           taxonID, 
           scientificName)
  
  perplotperyear_TALL.df<-PlantTALL$vst_perplotperyear%>%
    filter(plotID %in% plots_TALL.df$plotID)%>%     #CHECK EACH SITE: retain only latest fully sampled year
    select(plotID,                                  #Checking to see the plot area sampled
           nlcdClass,                               #TALL: don't select eventID because >1 sampling event
           totalSampledAreaTrees)%>%                #TALL: consistent 400 square m per site
    distinct()
  
  apparentindividual_TALL.df<-PlantTALL$vst_apparentindividual%>%
    filter(plotID %in% plots_TALL.df$plotID)%>%                              #retain all years where these plots were sampled because a 2nd complete samplying event has not occured since the initial one
    filter(growthForm=="multi-bole tree"|growthForm=="single bole tree")%>%  #keep only single & multi bole trees    
    filter(!stemDiameter < 10)%>%                                            #filter out trees <10 cm in Diameter
    filter(plantStatus=="Live"|
             plantStatus=="Live, insect damaged"|
             plantStatus=="Live,  other damage"|
             plantStatus=="Live, disease damaged"|
             plantStatus=="Live, physically damaged"|
             plantStatus=="Live, broken bole"|
             plantStatus=="other damage")%>%
    filter(!plotID=="TALL_003")%>%
    filter(!plotID=="TALL_005")%>%                        #removed plots 003 and 005 because they only have 1 and 3 trees larger than 10cm DBH
    select(plotID,
           individualID,
           stemDiameter)%>%
    left_join(., mappingandtagging_TALL.df, by = "individualID")%>%
    left_join(., perplotperyear_TALL.df, by = "plotID")%>%
    filter(!taxonID=="2PLANT")%>%
    left_join(., Mycorrhizal.df, by = "taxonID")%>%
    mutate(., basalAreaCm2=(3.142*((stemDiameter/2)^2)))%>%
    distinct()
  
}
#__________________________________________
# Create new data frame based on basal area ----
#__________________________________________
if(TRUE){
  BasalAreaTALL.df<-apparentindividual_TALL.df%>%                    #make a df with calculations based on basal area
    spread(treeType, basalAreaCm2)%>%                           #making a column for leaf habits, and the basal area of each tree to category
    mutate(basalAreaCm2=(3.142*((stemDiameter/2)^2)))%>%         
    spread(mycorrhizal, basalAreaCm2)%>%                         #making a column for leaf habits, and the basal area of each tree to category
    mutate(Evergreen = replace_na(Evergreen, 0),                 #replace na's with 0 
           Deciduous = replace_na(Deciduous, 0), 
           ECM = replace_na(ECM, 0), 
           AM = replace_na(AM, 0), 
           totalBA=Deciduous+Evergreen)%>%                      #calculate the totalBA for each plot
    group_by(plotID)%>%                                         #calculate the total area of deciduous, evergreen, am, ecm, species at the plot level
    summarise(Deciduous = sum(Deciduous),
              Evergreen = sum(Evergreen),
              AM = sum(AM),
              ECM = sum(ECM),
              totalBA = sum(totalBA))%>%                        #total species density at each plot by basal area of stems
    mutate(perEvergreen=Evergreen/totalBA*100,                  #calculate the percent evergreen species at each plot
           perECM=ECM/totalBA*100)%>%                           #calculate the percent ECM at each plot
    left_join(., apparentindividual_TALL.df, by = "plotID")%>%  #re-join taxonID to create a diversity calculation
    select(plotID,  
           nlcdClass,
           taxonID,
           basalAreaCm2,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalBA,
           perEvergreen,
           perECM)%>%
    group_by(plotID)%>%                                          #designate basal area of each species into columns, grouped by plotID
    pivot_wider(names_from = "taxonID", values_from = "basalAreaCm2", values_fn = list(basalAreaCm2 = sum))%>%
    replace(is.na(.),0)
  species.df<-BasalAreaTALL.df[,10:25]                               #CHECK EACH SITE:for vegan package, subset species into columns for diversity
  BasalAreaTALL.df$ShannonPlant<-diversity(species.df, index = "shannon")      #Shannon Diversity Index
  BasalAreaTALL.df$SimpsonPlant<-diversity(species.df, index = "simpson")      #Simpson Diversity Index
  BasalAreaTALL.df<-BasalAreaTALL.df%>%         
    select(plotID,                                              #retain necessary columns
           nlcdClass,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalBA,
           perEvergreen,
           perECM,
           ShannonPlant,
           SimpsonPlant)
  rm(species.df)                                                #remove the species df
}
#__________________________________________
# Create new data frame based on number of stems----
#__________________________________________
if(TRUE){
  StemsTALL.df<-apparentindividual_TALL.df%>%                       #make a df with calculations based on stems of each species
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%    #count each stem 
    spread(treeType, stemCount)%>%                              #count by leaf habit
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%
    spread(mycorrhizal, stemCount)%>%                           #count by fungal association
    mutate(Evergreen = replace_na(Evergreen, 0),                #replace NA's
           Deciduous = replace_na(Deciduous, 0),
           ECM = replace_na(ECM, 0),
           AM = replace_na(AM, 0),
           totalStems=Deciduous+Evergreen)%>%                   #total count of each species
    group_by(plotID)%>%                                         #summarize at the plot level
    summarise(Deciduous = sum(Deciduous),
              Evergreen = sum(Evergreen),
              AM = sum(AM),
              ECM = sum(ECM),
              totalStems = sum(totalStems))%>%
    mutate(perEvergreen=Evergreen/totalStems*100,               #make a percent evergreen calculation
           perECM=ECM/totalStems*100)%>%                        #make a percent ECM calculation
    left_join(., apparentindividual_TALL.df, by = "plotID")%>%  #join in species data
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%
    select(plotID,
           nlcdClass,
           taxonID,
           stemCount,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalStems,
           perECM,
           perEvergreen)%>%
    group_by(plotID)%>%                                         #regroup by plotID and place each species into columns
    pivot_wider(names_from = "taxonID", values_from = "stemCount", values_fn = list(stemCount = sum))%>%
    replace(is.na(.),0)
  speciesStems.df<-StemsTALL.df[,10:25]                             #CHECK SITE: Harv has 10 species, so subset columns 10-20
  StemsTALL.df$ShannonPlant<-diversity(speciesStems.df, index = "shannon")        #Shannon Diversity Index
  StemsTALL.df$SimpsonPlant<-diversity(speciesStems.df, index = "simpson")        #Simpson Diversity Index
  StemsTALL.df<-StemsTALL.df%>%                                         
    select(plotID,
           Deciduous,
           Evergreen,
           AM,
           ECM,
           totalStems,
           perEvergreen,
           perECM,
           ShannonPlant,
           SimpsonPlant)
  rm(speciesStems.df)
}
#__________________________________________
# Descriptive Statistics SITE level----
#__________________________________________
if(TRUE){
  SiteStatsTALL.df<-apparentindividual_TALL.df%>%                      #look at abundance of each species at HARV by stem and basal area
    select(.,plotID,
           taxonID,
           scientificName,
           basalAreaCm2)%>%
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%             #create a species count
    group_by(taxonID)%>%                                                 #group by species ID
    count(taxonID, taxonID, name = "stems")%>%                           #makes a column for the number of stems 
    left_join(., apparentindividual_TALL.df, by = "taxonID")%>%          #rejoin to add basal area column
    select(.,
           taxonID,
           stems,
           scientificName,
           basalAreaCm2)%>%
    group_by(taxonID, stems, scientificName)%>%
    summarise(basalAreaCm2 = (sum(basalAreaCm2)))                        #sum basal area by species ID
  SiteStatsTALL.df<-SiteStatsTALL.df%>%
    mutate(relAbunBA = basalAreaCm2/sum(SiteStatsTALL.df$basalAreaCm2)*100,  #calculate the relative adundance of a species at HARV
           relAbunSt = stems/sum(SiteStatsTALL.df$stems)*100)                #relative abundance = 1 species/ total species at site
}
#__________________________________________
# Descriptive Statistics PLOT level----
#__________________________________________
library(plotrix)
if(TRUE){
  PlotStatsTALL.df<-BasalAreaTALL.df%>%                                          #create a df with descriptive statists at the plot level
    select(plotID,                                                       #select parameters based on basal area
           perECM,
           perEvergreen,
           totalBA)%>%
    left_join(., StemsTALL.df, by = "plotID")%>%
    select(., plotID,                                                    #select parameters based on number of stems
           totalStems,
           totalBA,
           perECM.x,
           perEvergreen.x,
           perECM.y,
           perEvergreen.y
    )
  PlotStatsTALL.df<-PlotStatsTALL.df%>%                                          #calculate MEAN basal area per plot
    summarise(.,                                                         #calculate STD ERROR OF MEAN basal area per plot
              #calculate MIN & MAX basal area per plot
              avgBA = mean(PlotStatsTALL.df$totalBA), 
              seBA = std.error(PlotStatsTALL.df$totalBA),
              minBA = min(PlotStatsTALL.df$totalBA),
              maxBA = max(PlotStatsTALL.df$totalBA),
              sumtotalBA = sum(PlotStatsTALL.df$totalBA),
              
              avgperECMba = mean(PlotStatsTALL.df$perECM.x),                 #calculate MEAN percent ECM fungal association based on basal area per plot
              seperECMba = std.error(PlotStatsTALL.df$perECM.x),             #calculate STD ERROR OF MEAN percent ECM fungal association based on basal area per plot
              minperECMba = min(PlotStatsTALL.df$perECM.x),                  #calculate MIN & MAX percent ECM fungal association based on basal area per plot
              maxperECMba = max(PlotStatsTALL.df$perECM.x),
              
              avgperEGba = mean(PlotStatsTALL.df$perEvergreen.x),           #calculate MEAN percent Evergreen leaf habib by basal area per plot
              seperEGba = std.error(PlotStatsTALL.df$perEvergreen.x),       #calculate STD ERROR OF MEAN percent Evergreen leaf habib by basal area per plot
              minperEGba = min(PlotStatsTALL.df$perEvergreen.x),            #calculate MIN & MAX percent Evergreen leaf habib by basal area per plot
              maxperEGba = max(PlotStatsTALL.df$perEvergreen.x),
              
              avgStems = mean(PlotStatsTALL.df$totalStems),                  #MEAN total species by stem and per plot
              seStems = std.error(PlotStatsTALL.df$totalStems),              #STD ERROR OF MEAN by total species by stem and per plot
              minStems = min(PlotStatsTALL.df$totalStems),                   #MIN & MAX total species by stem and per plot
              maxStems = max(PlotStatsTALL.df$totalStems),
              sumtotalStems = sum(PlotStatsTALL.df$totalStems),
              
              avgperECMstems = mean(PlotStatsTALL.df$perECM.y),                #MEAN percent ECM fungal association by stem and per plot
              seperECMstems = std.error(PlotStatsTALL.df$perECM.y),            #STD ERROR OF MEAN percent ECM fungal association by stem and per plot
              minperECMstems = min(PlotStatsTALL.df$perECM.y),                 #MIC & MAX percent ECM fungal association by stem and per plot
              maxperECMstems = max(PlotStatsTALL.df$perECM.y),
              
              avgperEGstems = mean(PlotStatsTALL.df$perEvergreen.y),           #MEAN percent evergreen leaf habit by stem and per plot
              seperEGstems = std.error(PlotStatsTALL.df$perEvergreen.y),       #STD ERROR OF MEAN percent evergreen leaf habit by stem and per plot
              minperEGstms = min(PlotStatsTALL.df$perEvergreen.y),            #MIN & MAX percent evergreen leaf habit by stem and per plot
              maxperEGstms = max(PlotStatsTALL.df$perEvergreen.y),
              
              avgShannonBA = mean(BasalAreaTALL.df$ShannonPlant),
              seShannonBA = std.error(BasalAreaTALL.df$ShannonPlant),
              avgSimpsonBA = mean(BasalAreaTALL.df$SimpsonPlant),
              seSimpsonBA = std.error(BasalAreaTALL.df$SimpsonPlant),
              
              avgShannonStems = mean(StemsTALL.df$ShannonPlant),
              seShannonStems = std.error(StemsTALL.df$ShannonPlant),
              avgSimpsonStems = mean(StemsTALL.df$SimpsonPlant),
              seSimpsonStems = std.error(StemsTALL.df$SimpsonPlant))%>%
    select(!plotID)%>%                                                  #remove plotID to make a final df
    distinct()%>%
    pivot_longer(., everything(), 
                 names_to = "TALL", 
                 values_to = "Value")
}
#__________________________________________
# Remove unwanted df
#__________________________________________
if(TRUE){
  rm(apparentindividual_TALL.df,
     mappingandtagging_TALL.df,
     perplotperyear_TALL.df,
     plots_TALL.df)}

#__________________________________________
# join stem and BA into one df----
#__________________________________________
if(TRUE){
  allHARV.df<-left_join(BasalAreaHARV.df, StemsHARV.df, by = "plotID")%>%
    rename(., DF_BA = Deciduous.x,
           EGBA = Evergreen.x,
           AM_BA = AM.x,
           ECM_BA = ECM.x,
           PerEG_BA = perEvergreen.x,
           PerECM_BA = perECM.x,
           Shan_BA = ShannonPlant.x,
           Sim_BA = SimpsonPlant.x,
           DF_ST = Deciduous.y,
           EF_ST = Evergreen.y,
           AM_ST = AM.y,
           ECM_ST = ECM.y,
           PerEG_ST = perEvergreen.y,
           PerECM_ST = perECM.y,
           Shan_ST = ShannonPlant.y,
           Sim_ST = SimpsonPlant.y,
           nlcdClass = nlcdClass.x)%>%
    select(., !nlcdClass.y)
  
  
  allTALL.df<-left_join(BasalAreaTALL.df, StemsTALL.df, by = "plotID")%>%
    rename(., DF_BA = Deciduous.x,
           EGBA = Evergreen.x,
           AM_BA = AM.x,
           ECM_BA = ECM.x,
           PerEG_BA = perEvergreen.x,
           PerECM_BA = perECM.x,
           Shan_BA = ShannonPlant.x,
           Sim_BA = SimpsonPlant.x,
           DF_ST = Deciduous.y,
           EF_ST = Evergreen.y,
           AM_ST = AM.y,
           ECM_ST = ECM.y,
           PerEG_ST = perEvergreen.y,
           PerECM_ST = perECM.y,
           Shan_ST = ShannonPlant.y,
           Sim_ST = SimpsonPlant.y)
}

allVeg.df<-bind_rows(allTALL.df,
                     allHARV.df)

#__________________________________________
# Remove unwanted df
#__________________________________________
if(TRUE){
  rm(BasalAreaTALL.df,
     StemsHARV.df,
     StemsTALL.df)}





#__________________________________________
# BEETLE----
#__________________________________________

#__________________________________________
# Read in relevant raw data sets HARV----
#__________________________________________

sort_HARV.df <- BeetleHARV$bet_sorting%>%
  filter(collectDate < "2019-12-31",                                           # incomple sampling year
         collectDate > "2013-12-31",                                           # incomplete sampling year
         !plotID=="HARV_025",                                                  # plot sampling incomplete (only 3 events)
         sampleType == "carabid" | sampleType =="other carabid",               # removing species that are not ground beetles
         !setDate=="2014-09-12")%>%                                            # incomplete sampling data event
  mutate(taxonID = replace(taxonID, taxonID == 'APELUC', 'APELUC2'),           # changing subspecies to species
         taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%
  select(subsampleID,
         taxonID,
         individualCount,
         sampleID)


expertTax_HARV <- BeetleHARV$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate < "2019-12-31",
         collectDate > "2013-12-31",
         !plotID=="HARV_025",
         !setDate=="2014-09-12",
         family=="Carabidae")%>%
  mutate(taxonID = replace(taxonID, taxonID == 'APELUC', 'APELUC2'),
         taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%
  transform(scientificName=paste(genus, specificEpithet))%>%                   # renaming all species to be only genus and specific Epithet
  select(plotID,
         individualID,
         taxonID,
         scientificName)%>%
  distinct()

field_HARV.df <- BeetleHARV$bet_fielddata%>%
  mutate(eventID = replace(eventID, eventID %in% "HARV201427", "HARV.2014.27"),# changing typos 
         eventID = replace(eventID, eventID %in% "HARV201433", "HARV.2014.33"),
         eventID = replace(eventID, eventID %in% "HARV201435", "HARV.2014.35"),
         eventID = replace(eventID, eventID %in% "HARV201439", "HARV.2014.39"),
         eventID = replace(eventID, eventID %in% "HARV201528", "HARV.2015.28"),
         eventID = replace(eventID, eventID %in% "HARV201534", "HARV.2015.34"))%>%
  filter(collectDate < "2019-12-31",
         collectDate > "2013-12-31",
         !plotID=="HARV_025",
         !eventID=="HARV.2014.39")%>%
  select(plotID,
         trapID,
         nlcdClass,
         eventID,
         sampleID,
         cupStatus,
         lidStatus)%>%
  transform(., plotIDeventID=paste(plotID, eventID))%>%                        # make column with plotID & eventID for later grouping 
  distinct()%>%
  filter(., !sampleID == "HARV_005.S.20150611",                                # duplicated sampling locations
         !sampleID == "HARV_005.N.20150611",
         !sampleID == "HARV_005.W.20150611",
         !sampleID == "HARV_005.E.20150611",)


if(TRUE){
  paraTax_HARV <- BeetleHARV$bet_parataxonomistID%>%
    filter(collectDate < "2019-12-31",
           collectDate > "2013-12-31",
           !plotID=="HARV_025",
           !setDate=="2014-09-12")%>%
    mutate(taxonID = replace(taxonID, taxonID == 'APELUC', 'APELUC2'),
           taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
           taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
           taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%
    select(individualID,
           collectDate, 
           subsampleID,
           taxonID)%>%
    left_join(., expertTax_HARV, by = "individualID")%>%                         # joining expertly identified ground beetles to the pinned ground beetles
    distinct()%>%
    dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%                 # changing the non-expert ID's to expert ID's
    left_join(., sort_HARV.df, by = "subsampleID")%>%
    mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%
    group_by(subsampleID)%>%
    add_count(subsampleID, name = "beetleCount")%>%                              # number of beetles per subsample
    group_by(subsampleID)%>%
    add_count(taxonID.y, name = "speciesCount")%>%                               # number of species per subsample
    distinct(., subsampleID, .keep_all = TRUE)%>%
    mutate(individualCount = coalesce(individualCount, beetleCount))%>%          # change the individual count from the sort df to calculated number
    select(!taxonID.x)%>%
    select(!taxonID)%>%
    select(!beetleCount)%>%
    select(!speciesCount)%>%
    select(!plotID)%>%
    select(!scientificName)                                                      # make a row for each individual species
  paraTax_HARV<-paraTax_HARV[rep(row.names(paraTax_HARV), paraTax_HARV$individualCount), ]  
  paraTax_HARV<-field_HARV.df%>%                                                 # join in field data
    select(plotIDeventID,
           plotID,
           sampleID,
           nlcdClass)%>%
    distinct()%>%
    right_join(paraTax_HARV, field_HARV.df, by = "sampleID")%>%                  # !!wrong individual count- each row is a species - keep for later
    na.omit()
}


#__________________________________________
# Density PlotID-EventID grouping----
#__________________________________________
if(TRUE){
  density_HARV.df<-paraTax_HARV%>%                                               # make a density calculation by plotID and eventID
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, name = "BeetleCount")%>%
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    select(plotIDeventID,
           plotID,
           collectDate,
           nlcdClass,
           BeetleCount)%>%
    left_join(., field_HARV.df, by = "plotIDeventID")%>%                         #join field data
    select(!plotID.y)%>%
    select(!nlcdClass.y)%>%
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, plotIDeventID, name = "cupNum")%>%                  #calculated # cups at each plotID-eventID
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    mutate(density = BeetleCount/cupNum)%>%                                      # density calculated by total beetles at a plotID-eventID, divided by total cups at a plotID-eventID
    select(plotIDeventID,
           plotID.x,
           collectDate,
           nlcdClass.x,
           BeetleCount,
           cupNum,
           density)%>%
    mutate(plotID = plotID.x,
           nlcdClass = nlcdClass.x)
}
#__________________________________________
# Diversity plotID-eventID grouping----
#__________________________________________
if(TRUE){
  diversityHARV.df<-paraTax_HARV%>%                                              # make a diversity df from the transformed pin df
    select(plotIDeventID,
           taxonID.y,
           individualCount)%>%
    group_by(plotIDeventID, taxonID.y)%>%                                        # group by plotID-eventID and each taxonID - make a column for each species with the number of individuals 
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)%>%
    left_join(., density_HARV.df, by = "plotIDeventID")                          # join in density info and other field data 
  species.df<-diversityHARV.df[,2:34]                                            # species only df for vegan diversity analysis - right amount of species 2388
  diversityHARV.df$ShannonBeetle<-diversity(species.df, index = "shannon")       # shannon diversity index
  diversityHARV.df$SimpsonBeetle<-diversity(species.df, index = "simpson")       # simpson diversity index
  diversityHARV.df<-diversityHARV.df%>%
    select(plotIDeventID,
           plotID,
           nlcdClass,
           collectDate,
           BeetleCount,
           cupNum,
           density,
           ShannonBeetle,
           SimpsonBeetle)                                                        # retain necessary columns
  rm(species.df)
}
#__________________________________________
#  Density and Diversity df at plotID level----
#__________________________________________

if(TRUE){
  densityPLOT_HARV.df<-density_HARV.df%>%                                        # density calculation
    group_by(plotID)%>%
    summarise(cupNum = sum(cupNum),
              BeetleCount = sum(BeetleCount))%>%
    mutate(density = BeetleCount/cupNum)
  
  diversityPLOT_HARV.df<-paraTax_HARV%>%                                         # diversity df
    select(plotID,
           taxonID.y,
           individualCount)%>%
    group_by(plotID)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  speciesPLOT.df<-diversityPLOT_HARV.df[,2:34]                                 # right amount of species 2388
  diversityPLOT_HARV.df$ShannonBeetle<-diversity(speciesPLOT.df, index = "shannon")        
  diversityPLOT_HARV.df$SimpsonBeetle<-diversity(speciesPLOT.df, index = "simpson")
  diversityPLOT_HARV.df<-diversityPLOT_HARV.df%>%
    select(plotID,
           ShannonBeetle,
           SimpsonBeetle)%>%
    left_join(., densityPLOT_HARV.df, by = "plotID")%>%
    left_join(., field_HARV.df, by = "plotID")%>%
    distinct(., plotID, .keep_all = TRUE)%>%
    select(plotID,
           nlcdClass, 
           density,
           ShannonBeetle,
           SimpsonBeetle,
           BeetleCount)
  rm(speciesPLOT.df)
}

#__________________________________________
#  Density and Diversity at plotID-Year level----
#__________________________________________
if(TRUE){
  densityPlotYr_HARV.df<-density_HARV.df%>%                                      # density df 
    mutate(year = format(collectDate, format="%Y"))%>%                           # make a year column
    group_by(year, plotID)%>%                                                    # group by plot-year
    summarise(cupNum = sum(cupNum),
              BeetleCount = sum(BeetleCount))%>%
    mutate(density = BeetleCount/cupNum)%>%
    transform(., plotIDyear=paste(plotID, year))
  
  diversityPlotYr_HARV.df<-paraTax_HARV%>%                                       # diversity df
    mutate(year = format(collectDate, format="%Y"))%>%
    transform(., plotIDyear=paste(plotID, year))%>%
    select(plotID,
           taxonID.y,
           individualCount,
           year,
           plotIDyear)%>%
    group_by(year)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  speciesPlotYR.df<-diversityPlotYr_HARV.df[,4:36]                              # right amount of species 2388
  diversityPlotYr_HARV.df$ShannonBeetle<-diversity(speciesPlotYR.df, index = "shannon")        
  diversityPlotYr_HARV.df$SimpsonBeetle<-diversity(speciesPlotYR.df, index = "simpson")
  diversityPlotYr_HARV.df<-diversityPlotYr_HARV.df%>%
    select(plotIDyear,
           ShannonBeetle,
           SimpsonBeetle,
           year)%>%
    left_join(., densityPlotYr_HARV.df, by = "plotIDyear") %>%                   # join with density
    left_join(., field_HARV.df, by = "plotID")%>%                                # add field data
    distinct(., plotIDyear, .keep_all = TRUE)%>%
    select(plotID,
           plotIDyear,
           nlcdClass, 
           density, 
           ShannonBeetle,
           SimpsonBeetle,
           BeetleCount,
           cupNum,
           year.x)
  diversityPlotYr_HARV.df$siteID<- "HARV"
  rm(speciesPlotYR.df)
}


#__________________________________________
#  density of most abundant beetle species at plotID-Year level----
#__________________________________________
speciesDensityPlotYr_HARV.df<-density_HARV.df%>%                                       
  mutate(year = format(collectDate, format="%Y"))%>%                            #grab year from collectDate
  group_by(year, plotID)%>%
  summarise(cupNum = sum(cupNum),                                               #sum number of cups and beetles collected at each plot for each year
            BeetleCount = sum(BeetleCount))%>%
  mutate(density = BeetleCount/cupNum)%>%
  transform(., plotIDyear=paste(plotID, year))

speciesPlotYr_HARV.df<-paraTax_HARV%>%
  mutate(year = format(collectDate, format="%Y"))%>%
  transform(., plotIDyear=paste(plotID, year))%>%
  select(taxonID.y,
         year,
         individualCount,
         plotIDyear)%>%
  group_by(year)%>%
  pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)

speciesDensityPlotYr_HARV.df<-speciesDensityPlotYr_HARV.df%>%
  left_join(., speciesPlotYr_HARV.df, by = "plotIDyear")%>%
  mutate(CARGORdensity = CARGOR/cupNum)%>%
  mutate(SYNIMPdensity = SYNIMP/cupNum)%>%
  mutate(SPHSTE3density = SPHSTE3/cupNum)%>%
  mutate(PTETRI3density = PTETRI3/cupNum)%>%
  mutate(SPHCANdensity = SPHCAN/cupNum)%>%
  select(plotID,
         year.x,
         plotIDyear, 
         CARGORdensity, 
         SYNIMPdensity,
         SPHSTE3density, 
         PTETRI3density, 
         SPHCANdensity)

speciesDensityPlotYr_HARV.df<-speciesDensityPlotYr_HARV.df%>%
  left_join(., allVeg.df, by = "plotID")%>%                                        
  select(plotID,
         year.x,
         plotIDyear, 
         CARGORdensity, 
         SYNIMPdensity,
         SPHSTE3density, 
         PTETRI3density, 
         SPHCANdensity, 
         totalBA,
         PerEG_BA,
         PerECM_BA)


#__________________________________________
# TALL                              ----
#__________________________________________


sort_TALL.df <- BeetleTALL$bet_sorting%>%
  filter(collectDate < "2019-12-31",   
         collectDate > "2015-12-31",                                           #remove 2020 data: incomplete data 
         sampleType == "carabid" | sampleType =="other carabid")%>%            #remove samples that are not ground beetles
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),          #aggregate taxa to the species-level
         taxonID = replace(taxonID, taxonID == 'DICFUR2', 'DICFUR'),
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'),
         taxonID = replace(taxonID, taxonID == 'CICPUN1', 'CICPUN'))%>%
  select(subsampleID,
         taxonID,
         individualCount,
         sampleID)
sum(sort_TALL.df$individualCount)                                              #770 individuals in total 2014-2019


expertTax_TALL <- BeetleTALL$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2015-12-31",
         family=="Carabidae")%>%                                               #remove samples that are not ground beetles
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),          #change subspecies to species
         taxonID = replace(taxonID, taxonID == 'DICFUR2', 'DICFUR'),
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'),
         taxonID = replace(taxonID, taxonID == 'CICPUN1', 'CICPUN'))%>%
  transform(scientificName=paste(genus, specificEpithet))%>%                   #this code will convert subspecies names to the species name 
  select(plotID,
         individualID,
         taxonID,
         scientificName)%>%
  distinct()

field_TALL.df <- BeetleTALL$bet_fielddata%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2015-12-31",
         !sampleCollected=="N") %>%                                            #removing missing samples 
  select(plotID, 
         trapID,
         nlcdClass,
         eventID,
         sampleID,
         cupStatus,
         lidStatus)%>%
  transform(., plotIDeventID=paste(plotID, eventID))%>%
  distinct()

if(TRUE){
  paraTax_TALL <- BeetleTALL$bet_parataxonomistID%>%
    filter(collectDate < "2019-12-31",                                         #remove 2020 data: incomplete data 
           collectDate > "2015-12-31")%>%
    mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),        #change subspecies to species
           taxonID = replace(taxonID, taxonID == 'DICFUR2', 'DICFUR'),
           taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'),
           taxonID = replace(taxonID, taxonID == 'CICPUN1', 'CICPUN'))%>%
    select(individualID,
           collectDate, 
           subsampleID,
           taxonID)%>%
    left_join(., expertTax_TALL, by = "individualID")%>%                       #add expert ID data to paratax data via individual ID
    distinct()%>%
    dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%               #use expert IDs if available, otherwise use parataxonomist IDs
    left_join(., sort_TALL.df, by = "subsampleID")%>%                          #add sorting data to paratax data via subsample ID
    mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%                        #use expert/paratax IDs if available, otherwise sorting IDs
    group_by(subsampleID)%>%
    add_count(subsampleID, name = "beetleCount")%>%                            #sum the number of beeltes in each subsample
    group_by(subsampleID)%>%
    add_count(taxonID.y, name = "speciesCount")%>%                             #sum the number of species in each subsample
    distinct(., subsampleID, .keep_all = TRUE)%>%
    mutate(individualCount = coalesce(individualCount, beetleCount))%>%        #use individual count data if available, otherwise number of beetles
    select(!taxonID.x)%>%                                                      #remove unneeded columns
    select(!taxonID)%>%
    select(!beetleCount)%>%
    select(!speciesCount)%>%
    select(!plotID)%>%
    select(!scientificName)                                                    #755 individuals 
  paraTax_TALL<-paraTax_TALL[rep(row.names(paraTax_TALL), paraTax_TALL$individualCount), ]
  paraTax_TALL<-field_TALL.df%>%
    select(plotIDeventID,                                                      #select field data to add to paratax data 
           plotID,
           sampleID,
           nlcdClass)%>%
    distinct()%>%
    right_join(paraTax_TALL, field_TALL.df, by = "sampleID")%>%                #add field data to paratax data via sample ID
    select(plotID,                                                             #remove unneeded columns
           nlcdClass,
           collectDate, 
           plotIDeventID,
           sampleID, 
           subsampleID, 
           taxonID.y,
           individualCount)
  #total 755 individuals with clean data
}

#__________________________________________
# Density at plotIDeventID level   ----                                             
#__________________________________________
if(TRUE){
  densityTALL.df<-paraTax_TALL%>%  
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, name = "BeetleCount")%>%                           #sum the number of beetles collected from each plot at each sampling date
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    select(plotIDeventID,
           plotID,
           collectDate,
           nlcdClass,
           BeetleCount)%>%
    left_join(., field_TALL.df, by = "plotIDeventID")%>%                        #add field data to density data via plotIDeventID
    select(!plotID.y)%>%                                                        #remove unneeded columns
    select(!nlcdClass.y)%>%
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, plotIDeventID, name = "cupNum")%>%                 #use total "events" at each plot at each collection time to calcualte number of cups used for collection at each collection time
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    mutate(density = BeetleCount/cupNum)%>%                                     #calculate beetle density at each plot per event
    select(plotIDeventID,                                                       #remove unneeded columns
           plotID.x,
           collectDate,
           nlcdClass.x,
           BeetleCount,
           cupNum,
           density)%>%
    mutate(plotID = plotID.x,
           nlcdClass = nlcdClass.x)
  sum(densityTALL.df$BeetleCount) #correct!
}

#__________________________________________
#  density and diversity at plotID-Year level----
#__________________________________________
if(TRUE){
  densityPlotYr_TALL.df<-densityTALL.df%>%                                       
    mutate(year = format(collectDate, format="%Y"))%>%                            #grab year from collectDate
    group_by(year, plotID)%>%
    summarise(cupNum = sum(cupNum),                                               #sum number of cups and beetles collected at each plot for each year
              BeetleCount = sum(BeetleCount))%>%
    mutate(density = BeetleCount/cupNum)%>%
    transform(., plotIDyear=paste(plotID, year))
  
  diversityPlotYr_TALL.df<-paraTax_TALL%>%
    mutate(year = format(collectDate, format="%Y"))%>%
    transform(., plotIDyear=paste(plotID, year))%>%
    select(plotID,
           taxonID.y,
           individualCount,
           year,
           plotIDyear)%>%
    group_by(year)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  speciesPlotYR.df<-diversityPlotYr_TALL.df[,4:50]# right amount of species 750
  diversityPlotYr_TALL.df$ShannonBeetle<-diversity(speciesPlotYR.df, index = "shannon") #don't be alarmed it doesn't show up right away       
  diversityPlotYr_TALL.df$SimpsonBeetle<-diversity(speciesPlotYR.df, index = "simpson")
  diversityPlotYr_TALL.df<-diversityPlotYr_TALL.df%>%
    select(plotIDyear,
           ShannonBeetle,
           SimpsonBeetle)%>%
    left_join(., densityPlotYr_TALL.df, by = "plotIDyear")%>%
    left_join(., field_TALL.df, by = "plotID")%>%
    distinct(., plotIDyear, .keep_all = TRUE)%>%
    select(plotID,
           plotIDyear,
           nlcdClass, 
           density, 
           ShannonBeetle,
           SimpsonBeetle,
           BeetleCount,
           cupNum)
  diversityPlotYr_TALL.df$siteID<- "TALL"
  rm(speciesPlotYR.df)
}

#__________________________________________
#  density of most abundant beetle species at plotID-Year level----
#__________________________________________
speciesDensityPlotYr_TALL.df<-densityTALL.df%>%                                       
  mutate(year = format(collectDate, format="%Y"))%>%                            #grab year from collectDate
  group_by(year, plotID)%>%
  summarise(cupNum = sum(cupNum),                                               #sum number of cups and beetles collected at each plot for each year
            BeetleCount = sum(BeetleCount))%>%
  mutate(density = BeetleCount/cupNum)%>%
  transform(., plotIDyear=paste(plotID, year))

speciesPlotYr_TALL.df<-paraTax_TALL%>%
  mutate(year = format(collectDate, format="%Y"))%>%
  transform(., plotIDyear=paste(plotID, year))%>%
  select(taxonID.y,
         year,
         individualCount,
         plotIDyear)%>%
  group_by(year)%>%
  pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)

speciesDensityPlotYr_TALL.df<-speciesDensityPlotYr_TALL.df%>%
  left_join(., speciesPlotYr_TALL.df, by = "plotIDyear")%>%
  mutate(CYCCON2density = CYCCON2/cupNum)%>%
  mutate(ANIMER2density = ANIMER2/cupNum)%>%
  mutate(SELOPAdensity = SELOPA/cupNum)%>%
  mutate(ANIHAPdensity = ANIHAP/cupNum)%>%
  mutate(DICDIL5density = DICDIL5/cupNum)%>%
  select(plotID,
         year.x,
         plotIDyear, 
         CYCCON2density, 
         ANIMER2density,
         SELOPAdensity, 
         ANIHAPdensity, 
         DICDIL5density)

speciesDensityPlotYr_TALL.df<-speciesDensityPlotYr_TALL.df%>%
  left_join(., allVeg.df, by = "plotID")%>%
  filter(., !plotID=="TALL_003",
         !plotID=="TALL_005")%>%                                          #plots without veg data 
  select(plotID,
         year.x,
         plotIDyear, 
         CYCCON2density, 
         ANIMER2density,
         SELOPAdensity, 
         ANIHAPdensity, 
         DICDIL5density, 
         totalBA,
         PerEG_BA,
         PerECM_BA)



#__________________________________________
# Read in relevant raw data sets JERC ----
#__________________________________________

sort_JERC.df <- BeetleJERC$bet_sorting%>%
  filter(collectDate > "2013-12-31",                                           #remove 2013 data: incomplete data 
         !plotID=="JERC_022",                                                  #remove plot 022 (only sample for 2 events), 023 (not eg, df, or mf)
         !plotID=="JERC_023",                       
         sampleType == "carabid" | sampleType =="other carabid",               #remove samples that are not ground beetles - leave in carabid identified to family, maybe expertly ID'ed
         !setDate=="2014-09-12")%>%
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'))%>%       #samplying event not at JERC
  select(subsampleID,
         taxonID,
         individualCount,
         sampleID)
sum(sort_JERC.df$individualCount)                                              #455 species in total 2014-2020


expertTax_JERC <- BeetleJERC$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate > "2013-12-31",
         !plotID=="JERC_022",
         !plotID=="JERC_023",
         !setDate=="2014-09-12",
         family=="Carabidae")%>%
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'))%>%                                              #there are three subspecies: NOT changing taxon ID becuase species not present DICDIL1, DICFUR1, TETCAR1
  transform(scientificName=paste(genus, specificEpithet))%>%                   #this code will convert their subspecies to the species name 
  select(plotID,
         individualID,
         taxonID,
         scientificName)%>%
  distinct()

field_JERC.df <- BeetleJERC$bet_fielddata%>%
  mutate(eventID = replace(eventID, eventID %in% "JERC201528", "JERC.2015.28"), 
         eventID = replace(eventID, eventID %in% "JERC201437", "JERC.2014.37"),
         eventID = replace(eventID, eventID %in% "JERC201534", "JERC.2015.34"),
         eventID = replace(eventID, eventID %in% "JERC201530", "JERC.2015.30"))%>% 
  filter(collectDate > "2013-12-31",
         !plotID=="JERC_022",
         !plotID=="JERC_023",
         !eventID=="JERC.2014.39",
         !sampleCollected=="N")%>%                                             #removing missing sample | Meghan should I remove disturbed cups? - some still appear in sort / expert df
  select(siteID,
         plotID, 
         trapID,
         nlcdClass,
         eventID,
         sampleID,
         cupStatus,
         lidStatus)%>%
  transform(., plotIDeventID=paste(plotID, eventID))%>%
  distinct()


if(TRUE){
  paraTax_JERC <- BeetleJERC$bet_parataxonomistID%>%
    filter(collectDate > "2013-12-31",
           !plotID=="JERC_022",
           !plotID=="JERC_023")%>%
    mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),
           taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'))%>%  
    select(individualID,
           collectDate, 
           subsampleID,
           taxonID)%>% #411
    left_join(., expertTax_JERC, by = "individualID")%>%
    distinct()%>%
    dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%
    left_join(., sort_JERC.df, by = "subsampleID")%>%
    mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%
    group_by(subsampleID)%>%
    add_count(subsampleID, name = "beetleCount")%>%
    group_by(subsampleID)%>%
    add_count(taxonID.y, name = "speciesCount")%>%
    distinct(., subsampleID, .keep_all = TRUE)%>%
    mutate(individualCount = coalesce(individualCount, beetleCount))%>%
    select(!taxonID.x)%>%
    select(!taxonID)%>%
    select(!beetleCount)%>%
    select(!speciesCount)%>%
    select(!plotID)%>%
    select(!scientificName)
  paraTax_JERC<-paraTax_JERC[rep(row.names(paraTax_JERC), paraTax_JERC$individualCount), ]
  paraTax_JERC<-field_JERC.df%>%
    select(siteID,
           plotIDeventID,
           plotID,
           sampleID,
           nlcdClass)%>%
    distinct()%>%
    right_join(paraTax_JERC, field_JERC.df, by = "sampleID")                   #don't worry about wrong individual count- each row is a species
  #total 417 species with clean data
}

#__________________________________________
# Density----
#__________________________________________
if(TRUE){
  density.df<-paraTax_JERC%>%  
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, name = "BeetleCount")%>%
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    select(plotIDeventID,
           plotID,
           collectDate,
           nlcdClass,
           BeetleCount)%>%
    left_join(., field_JERC.df, by = "plotIDeventID")%>%
    select(!plotID.y)%>%
    select(!nlcdClass.y)%>%
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, plotIDeventID, name = "cupNum")%>%
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    mutate(density = BeetleCount/cupNum)%>%
    select(plotIDeventID,
           plotID.x,
           collectDate,
           nlcdClass.x,
           BeetleCount,
           cupNum,
           density)%>%
    mutate(plotID = plotID.x,
           nlcdClass = nlcdClass.x)
  sum(density.df$BeetleCount) #correcto!
}
#__________________________________________
# Diversity plotIDeventID level----
#__________________________________________
if(TRUE){
  diversityJERC.df<-paraTax_JERC%>%
    select(plotIDeventID,
           taxonID.y,
           individualCount)%>%
    group_by(plotIDeventID, taxonID.y)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)%>%
    left_join(., density.df, by = "plotIDeventID")
  species.df<-diversityJERC.df[,2:44]
  diversityJERC.df$ShannonBeetle<-diversity(species.df, index = "shannon")        
  diversityJERC.df$SimpsonBeetle<-diversity(species.df, index = "simpson")
  diversityJERC.df<-diversityJERC.df%>%
    select(plotIDeventID,
           plotID,
           nlcdClass,
           collectDate,
           BeetleCount,
           cupNum,
           density,
           ShannonBeetle,
           SimpsonBeetle)
  rm(species.df)
}
#__________________________________________
#  plotID level denstiy ----
#__________________________________________

if(TRUE){
  densityPLOT.df<-density.df%>%
    group_by(plotID)%>%
    summarise(cupNum = sum(cupNum),
              BeetleCount = sum(BeetleCount))%>%
    mutate(density = BeetleCount/cupNum)
  
  diversityPLOT_JERC.df<-paraTax_JERC%>%
    select(plotID,
           taxonID.y,
           individualCount)%>%
    group_by(plotID)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  speciesPLOT.df<-diversityPLOT_JERC.df[,2:44]
  diversityPLOT_JERC.df$ShannonBeetle<-diversity(speciesPLOT.df, index = "shannon")        
  diversityPLOT_JERC.df$SimpsonBeetle<-diversity(speciesPLOT.df, index = "simpson")
  diversityPLOT_JERC.df<-diversityPLOT_JERC.df%>%
    select(plotID,
           ShannonBeetle,
           SimpsonBeetle)%>%
    left_join(., densityPLOT.df, by = "plotID")%>%
    left_join(., field_JERC.df, by = "plotID")%>%
    distinct(., plotID, .keep_all = TRUE)%>%
    select(plotID,
           nlcdClass, 
           density,
           ShannonBeetle,
           SimpsonBeetle,
           BeetleCount)
  rm(speciesPLOT.df)
}
#__________________________________________
#  diversity & density plotID-Year level----
#__________________________________________
if(TRUE){
  densityPlotYr_JERC.df<-density.df%>%
    mutate(year = format(collectDate, format="%Y"))%>%
    group_by(year, plotID)%>%
    summarise(cupNum = sum(cupNum),
              BeetleCount = sum(BeetleCount))%>%
    mutate(density = BeetleCount/cupNum)%>%
    transform(., plotIDyear=paste(plotID, year))
  
  diversityPlotYr_JERC.df<-paraTax_JERC%>%
    mutate(year = format(collectDate, format="%Y"))%>%
    transform(., plotIDyear=paste(plotID, year))%>%
    select(plotID,
           taxonID.y,
           individualCount,
           year,
           plotIDyear)%>%
    group_by(year)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  speciesPlotYR.df<-diversityPlotYr_JERC.df[,4:46]# right amount of species 417
  diversityPlotYr_JERC.df$ShannonBeetle<-diversity(speciesPlotYR.df, index = "shannon")        
  diversityPlotYr_JERC.df$SimpsonBeetle<-diversity(speciesPlotYR.df, index = "simpson")
  diversityPlotYr_JERC.df<-diversityPlotYr_JERC.df%>%
    select(plotIDyear,
           ShannonBeetle,
           SimpsonBeetle)%>%
    left_join(., densityPlotYr_JERC.df, by = "plotIDyear")%>% #need to make a density by PlotYr
    left_join(., field_JERC.df, by = "plotID")%>%
    distinct(., plotIDyear, .keep_all = TRUE)%>%
    select(plotID,
           plotIDyear,
           nlcdClass, 
           density, #need to make a density by PlotYr
           ShannonBeetle,
           SimpsonBeetle,
           BeetleCount,
           cupNum)
  diversityPlotYr_JERC.df$siteID<- "JERC"
  rm(speciesPlotYR.df)
}

#__________________________________________
# Read in relevant data BART----
#__________________________________________

# BART 
sort_BART.df <- BeetleBART$bet_sorting%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         !plotID=="BART_005",                                                  #remove plot 005,006 (only sample for each) 
         !plotID=="BART_006",
         !plotID=="BART_008",
         sampleType == "carabid" | sampleType =="other carabid")%>%            #remove samples that are not ground beetles - leave in carabid identified to family, maybe expertly ID'ed 
  mutate(taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%         #changing subspecies taxonID to species taxonID
  select(subsampleID,
         taxonID,
         individualCount,
         sampleID)
sum(sort_BART.df$individualCount)                                              #9938 species in total 2014-2029


expertTax_BART <- BeetleBART$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate < "2019-12-31",                                           
         !plotID=="BART_005",                                                 
         !plotID=="BART_006",
         !plotID=="BART_008",
         family=="Carabidae")%>%                                               #there are 3 subspecies: only CYMPLA3 has a species level present - changin scientific name but leaving taxonID
  mutate(taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%
  transform(scientificName=paste(genus, specificEpithet))%>%                   #this code will convert their subspecies to the species name 
  select(plotID,
         individualID,
         taxonID,
         scientificName)%>%
  distinct()


field_BART.df <- BeetleBART$bet_fielddata%>%
  mutate(eventID = replace(eventID, eventID %in% "BART201426", "BART.2014.26"), 
         eventID = replace(eventID, eventID %in% "BART201428", "BART.2014.28"),
         eventID = replace(eventID, eventID %in% "BART201436", "BART.2014.36"))%>%
  filter(collectDate < "2019-12-31",                                           
         !plotID=="BART_005",                                                 
         !plotID=="BART_006",
         !plotID=="BART_008",
         !sampleCollected=="N")%>%                                             #removing missing sample | Meghan should I remove disturbed cups? - some still appear in sort / expert df
  select(plotID, 
         trapID,
         nlcdClass,
         eventID,
         sampleID,
         cupStatus,
         lidStatus)%>%
  transform(., plotIDeventID=paste(plotID, eventID))%>%
  distinct()


if(TRUE){
  paraTax_BART <- BeetleBART$bet_parataxonomistID%>%
    filter(collectDate < "2019-12-31",                                           
           !plotID=="BART_005",                                                 
           !plotID=="BART_006",
           !plotID=="BART_008")%>%
    mutate(taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
           taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
           taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%
    select(individualID,
           collectDate, 
           subsampleID,
           taxonID)%>%
    left_join(., expertTax_BART, by = "individualID")%>%
    distinct()%>%
    dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%
    left_join(., sort_BART.df, by = "subsampleID")%>%
    mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%
    group_by(subsampleID)%>%
    add_count(subsampleID, name = "beetleCount")%>%
    group_by(subsampleID)%>%
    add_count(taxonID.y, name = "speciesCount")%>%
    distinct(., subsampleID, .keep_all = TRUE)%>%
    mutate(individualCount = coalesce(individualCount, beetleCount))%>%
    select(!taxonID.x)%>%
    select(!taxonID)%>%
    select(!beetleCount)%>%
    select(!speciesCount)%>%
    select(!plotID)%>%
    select(!scientificName)                                                    #3771 species 
  paraTax_BART<-paraTax_BART[rep(row.names(paraTax_BART), paraTax_BART$individualCount), ]
  paraTax_BART<-field_BART.df%>%
    select(plotIDeventID,
           plotID,
           sampleID,
           nlcdClass)%>%
    distinct()%>%
    right_join(paraTax_BART, field_BART.df, by = "sampleID")                   #don't worry about wrong individual count- each row is a species
  #total 3771 species with clean data
}

#__________________________________________
# Density----
#__________________________________________
if(TRUE){
  densityBART.df<-paraTax_BART%>%  
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, name = "BeetleCount")%>%
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    select(plotIDeventID,
           plotID,
           collectDate,
           nlcdClass,
           BeetleCount)%>%
    left_join(., field_BART.df, by = "plotIDeventID")%>%
    select(!plotID.y)%>%
    select(!nlcdClass.y)%>%
    group_by(plotIDeventID)%>%
    add_count(plotIDeventID, plotIDeventID, name = "cupNum")%>%
    distinct(., plotIDeventID, .keep_all = TRUE)%>%
    mutate(density = BeetleCount/cupNum)%>%
    select(plotIDeventID,
           plotID.x,
           collectDate,
           nlcdClass.x,
           BeetleCount,
           cupNum,
           density)%>%
    mutate(plotID = plotID.x,
           nlcdClass = nlcdClass.x)
  sum(densityBART.df$BeetleCount) #correcto!
}
#__________________________________________
# Diversity plotIDeventID level----
#__________________________________________
if(TRUE){
  diversityBART.df<-paraTax_BART%>%
    select(plotIDeventID,
           taxonID.y,
           individualCount)%>%
    group_by(plotIDeventID, taxonID.y)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)%>%
    left_join(., densityBART.df, by = "plotIDeventID")
  species.df<-diversityBART.df[,2:25]                                          #correct 3771 species
  diversityBART.df$ShannonBeetle<-diversity(species.df, index = "shannon")         
  diversityBART.df$SimpsonBeetle<-diversity(species.df, index = "simpson")
  diversityBART.df<-diversityBART.df%>%
    select(plotIDeventID,
           plotID,
           nlcdClass,
           collectDate,
           BeetleCount,
           cupNum,
           density,
           ShannonBeetle,
           SimpsonBeetle)
  rm(species.df)
}
#__________________________________________
#  plotID level denstiy ----
#__________________________________________

if(TRUE){
  densityPLOT_BART.df<-densityBART.df%>%
    group_by(plotID)%>%
    summarise(cupNum = sum(cupNum),
              BeetleCount = sum(BeetleCount))%>%
    mutate(density = BeetleCount/cupNum)
  
  diversityPLOT_BART.df<-paraTax_BART%>%
    select(plotID,
           taxonID.y,
           individualCount)%>%
    group_by(plotID)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  speciesPLOT.df<-diversityPLOT_BART.df[,2:25]
  diversityPLOT_BART.df$ShannonBeetle<-diversity(speciesPLOT.df, index = "shannon")        
  diversityPLOT_BART.df$SimpsonBeetle<-diversity(speciesPLOT.df, index = "simpson")
  diversityPLOT_BART.df<-diversityPLOT_BART.df%>%
    select(plotID,
           ShannonBeetle,
           SimpsonBeetle)%>%
    left_join(., densityPLOT_BART.df, by = "plotID")%>%
    left_join(., field_BART.df, by = "plotID")%>%
    distinct(., plotID, .keep_all = TRUE)%>%
    select(plotID,
           nlcdClass, 
           density,
           ShannonBeetle,
           SimpsonBeetle,
           BeetleCount)
  rm(speciesPLOT.df)
}
#__________________________________________
#  diversity plotID-Year level----
#__________________________________________
if(TRUE){
  densityPlotYr_BART.df<-densityBART.df%>%
    mutate(year = format(collectDate, format="%Y"))%>%
    group_by(year, plotID)%>%
    summarise(cupNum = sum(cupNum),
              BeetleCount = sum(BeetleCount))%>%
    mutate(density = BeetleCount/cupNum)%>%
    transform(., plotIDyear=paste(plotID, year))
  
  diversityPlotYr_BART.df<-paraTax_BART%>%
    mutate(year = format(collectDate, format="%Y"))%>%
    transform(., plotIDyear=paste(plotID, year))%>%
    select(plotID,
           taxonID.y,
           individualCount,
           year,
           plotIDyear)%>%
    group_by(year)%>%
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  speciesPlotYR.df<-diversityPlotYr_BART.df[,4:27]# right amount of species 3771
  diversityPlotYr_BART.df$ShannonBeetle<-diversity(speciesPlotYR.df, index = "shannon")        
  diversityPlotYr_BART.df$SimpsonBeetle<-diversity(speciesPlotYR.df, index = "simpson")
  diversityPlotYr_BART.df<-diversityPlotYr_BART.df%>%
    select(plotIDyear,
           ShannonBeetle,
           SimpsonBeetle)%>%
    left_join(., densityPlotYr_BART.df, by = "plotIDyear")%>% 
    left_join(., field_BART.df, by = "plotID")%>%
    distinct(., plotIDyear, .keep_all = TRUE)%>%
    select(plotID,
           plotIDyear,
           nlcdClass, 
           density, 
           ShannonBeetle,
           SimpsonBeetle,
           BeetleCount,
           cupNum)
  diversityPlotYr_BART.df$siteID<- "BART"
  rm(speciesPlotYR.df)
}

#__________________________________________
#Descriptive Stats----
#__________________________________________
if(TRUE){

BART.df<-diversityPlotYr_BART.df
TALL.df<-diversityPlotYr_TALL.df
HARV.df<-diversityPlotYr_HARV.df
JERC.df<-diversityPlotYr_JERC.df

#__________________________________________
# BART stats----
#__________________________________________

statsBART.df<-BART.df%>%                                        #create a df with descriptive statists at the plot level
  select(plotID,                                                             #select parameters based on basal area
         density,
         ShannonBeetle,
         SimpsonBeetle,
         BeetleCount)                                   #calculate MEAN basal area per plot

statsBART.df<-statsBART.df%>%     
  summarise(.,   
            avgDensity = mean(statsBART.df$density), 
            seDensity = std.error(statsBART.df$density),
            minDensity = min(statsBART.df$density),
            maxDensity = max(statsBART.df$density),
            
            avgShanB = mean(statsBART.df$ShannonBeetle),                   #calculate MEAN percent ECM fungal association based on basal area per plot
            seShanB = std.error(statsBART.df$ShannonBeetle),               #calculate STD ERROR OF MEAN percent ECM fungal association based on basal area per plot
            minShanB = min(statsBART.df$ShannonBeetle),                    #calculate MIN & MAX percent ECM fungal association based on basal area per plot
            maxShanB = max(statsBART.df$ShannonBeetle),
            
            avgSimB = mean(statsBART.df$SimpsonBeetle),              #calculate MEAN percent Evergreen leaf habib by basal area per plot
            seSimB = std.error(statsBART.df$SimpsonBeetle),          #calculate STD ERROR OF MEAN percent Evergreen leaf habib by basal area per plot
            minSimB = min(statsBART.df$SimpsonBeetle),               #calculate MIN & MAX percent Evergreen leaf habib by basal area per plot
            maxSimB = max(statsBART.df$SimpsonBeetle),
            
            avgBeetle = mean(statsBART.df$BeetleCount),                    #MEAN total species by stem and per plot
            seBeetle = std.error(statsBART.df$BeetleCount),                #STD ERROR OF MEAN by total species by stem and per plot
            minBeetle = min(statsBART.df$BeetleCount),                     #MIN & MAX total species by stem and per plot
            maxBeetle = max(statsBART.df$BeetleCount),
            sumtotalBeetle = sum(statsBART.df$BeetleCount))%>%        #remove plotID to make a final df
  distinct()%>%
  pivot_longer(., everything(), 
               names_to = "BART", 
               values_to = "Value")

#__________________________________________
#HARV stats----
#__________________________________________

statsHARV.df<-HARV.df%>%                                        #create a df with descriptive statists at the plot level
  select(plotID,                                                             #select parameters based on basal area
         density,
         ShannonBeetle,
         SimpsonBeetle,
         BeetleCount)  

statsHARV.df<-statsHARV.df%>%     
  summarise(.,   
            avgDensity = mean(statsHARV.df$density), 
            seDensity = std.error(statsHARV.df$density),
            minDensity = min(statsHARV.df$density),
            maxDensity = max(statsHARV.df$density),
            
            avgShanB = mean(statsHARV.df$ShannonBeetle),                   #calculate MEAN percent ECM fungal association based on basal area per plot
            seShanB = std.error(statsHARV.df$ShannonBeetle),               #calculate STD ERROR OF MEAN percent ECM fungal association based on basal area per plot
            minShanB = min(statsHARV.df$ShannonBeetle),                    #calculate MIN & MAX percent ECM fungal association based on basal area per plot
            maxShanB = max(statsHARV.df$ShannonBeetle),
            
            avgSimB = mean(statsHARV.df$SimpsonBeetle),              #calculate MEAN percent Evergreen leaf habib by basal area per plot
            seSimB = std.error(statsHARV.df$SimpsonBeetle),          #calculate STD ERROR OF MEAN percent Evergreen leaf habib by basal area per plot
            minSimB = min(statsHARV.df$SimpsonBeetle),               #calculate MIN & MAX percent Evergreen leaf habib by basal area per plot
            maxSimB = max(statsHARV.df$SimpsonBeetle),
            
            avgBeetle = mean(statsHARV.df$BeetleCount),                    #MEAN total species by stem and per plot
            seBeetle = std.error(statsHARV.df$BeetleCount),                #STD ERROR OF MEAN by total species by stem and per plot
            minBeetle = min(statsHARV.df$BeetleCount),                     #MIN & MAX total species by stem and per plot
            maxBeetle = max(statsHARV.df$BeetleCount),
            sumtotalBeetle = sum(statsHARV.df$BeetleCount))%>%        #remove plotID to make a final df
  distinct()%>%
  pivot_longer(., everything(), 
               names_to = "HARV", 
               values_to = "Value")

#__________________________________________
# JERC stats----
#__________________________________________

statsJERC.df<-JERC.df%>%                                        #create a df with descriptive statists at the plot level
  select(plotID,                                                             #select parameters based on basal area
         density,
         ShannonBeetle,
         SimpsonBeetle,
         BeetleCount)

statsJERC.df<-statsJERC.df%>% 
  summarise(.,   
            avgDensity = mean(statsJERC.df$density), 
            seDensity = std.error(statsJERC.df$density),
            minDensity = min(statsJERC.df$density),
            maxDensity = max(statsJERC.df$density),
            
            avgShanB = mean(statsJERC.df$ShannonBeetle),                   #calculate MEAN percent ECM fungal association based on basal area per plot
            seShanB = std.error(statsJERC.df$ShannonBeetle),               #calculate STD ERROR OF MEAN percent ECM fungal association based on basal area per plot
            minShanB = min(statsJERC.df$ShannonBeetle),                    #calculate MIN & MAX percent ECM fungal association based on basal area per plot
            maxShanB = max(statsJERC.df$ShannonBeetle),
            
            avgSimB = mean(statsJERC.df$SimpsonBeetle),              #calculate MEAN percent Evergreen leaf habib by basal area per plot
            seSimB = std.error(statsJERC.df$SimpsonBeetle),          #calculate STD ERROR OF MEAN percent Evergreen leaf habib by basal area per plot
            minSimB = min(statsJERC.df$SimpsonBeetle),               #calculate MIN & MAX percent Evergreen leaf habib by basal area per plot
            maxSimB = max(statsJERC.df$SimpsonBeetle),
            
            avgBeetle = mean(statsJERC.df$BeetleCount),                    #MEAN total species by stem and per plot
            seBeetle = std.error(statsJERC.df$BeetleCount),                #STD ERROR OF MEAN by total species by stem and per plot
            minBeetle = min(statsJERC.df$BeetleCount),                     #MIN & MAX total species by stem and per plot
            maxBeetle = max(statsJERC.df$BeetleCount),
            sumtotalBeetle = sum(statsJERC.df$BeetleCount))%>%        #remove plotID to make a final df
  distinct()%>%
  pivot_longer(., everything(), 
               names_to = "JERC", 
               values_to = "Value")

#__________________________________________
# TALL stats----
#__________________________________________


statsTALL.df<-TALL.df%>%                                        #create a df with descriptive statists at the plot level
  select(plotID,                                                             #select parameters based on basal area
         density,
         ShannonBeetle,
         SimpsonBeetle,
         BeetleCount)                                       #calculate MEAN basal area per plot

statsTALL.df<-statsTALL.df%>%
  summarise(.,   
            avgDensity = mean(statsTALL.df$density), 
            seDensity = std.error(statsTALL.df$density),
            minDensity = min(statsTALL.df$density),
            maxDensity = max(statsTALL.df$density),
            
            avgShanB = mean(statsTALL.df$ShannonBeetle),                   #calculate MEAN percent ECM fungal association based on basal area per plot
            seShanB = std.error(statsTALL.df$ShannonBeetle),               #calculate STD ERROR OF MEAN percent ECM fungal association based on basal area per plot
            minShanB = min(statsTALL.df$ShannonBeetle),                    #calculate MIN & MAX percent ECM fungal association based on basal area per plot
            maxShanB = max(statsTALL.df$ShannonBeetle),
            
            avgSimB = mean(statsTALL.df$SimpsonBeetle),              #calculate MEAN percent Evergreen leaf habib by basal area per plot
            seSimB = std.error(statsTALL.df$SimpsonBeetle),          #calculate STD ERROR OF MEAN percent Evergreen leaf habib by basal area per plot
            minSimB = min(statsTALL.df$SimpsonBeetle),               #calculate MIN & MAX percent Evergreen leaf habib by basal area per plot
            maxSimB = max(statsTALL.df$SimpsonBeetle),
            
            avgBeetle = mean(statsTALL.df$BeetleCount),                    #MEAN total species by stem and per plot
            seBeetle = std.error(statsTALL.df$BeetleCount),                #STD ERROR OF MEAN by total species by stem and per plot
            minBeetle = min(statsTALL.df$BeetleCount),                     #MIN & MAX total species by stem and per plot
            maxBeetle = max(statsTALL.df$BeetleCount),
            sumtotalBeetle = sum(statsTALL.df$BeetleCount))%>%        #remove plotID to make a final df
  distinct()%>%
  pivot_longer(., everything(), 
               names_to = "TALL", 
               values_to = "Value")


}

#__________________________________________
# Figure 1: Shan/Sim/ D BEETLE Y ~ X nlcdClass----
#__________________________________________
if(TRUE){
  

vegB.df<-bind_rows(HARV.df, TALL.df, JERC.df, BART.df)
  
p1 <- ggplot(vegB.df, aes(x=nlcdClass, y=ShannonBeetle, fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = vegB.df, aes(shape = siteID), size = 2, alpha = 0.5, width = 0.25, show.legend = T)+
  scale_shape_manual(values=c(15, 16,18,17))+
  scale_fill_brewer(palette="Dark2", guide=FALSE)+
  labs(x = NULL, y = 'Ground Beetle Diversity\n(Shannon Index)')+
  scale_x_discrete(name =NULL, 
                   limits = c("deciduousForest","mixedForest","evergreenForest"), 
                   labels = NULL)+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="top",
    legend.direction="horizontal")+
  scale_y_continuous(breaks=c(0,.25,1,1.75,2))
p1


p2 <- ggplot(vegB.df, aes(x=nlcdClass, y=SimpsonBeetle, fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = vegB.df, aes(shape = siteID), size = 2, alpha = 0.5, width = 0.25, show.legend = FALSE)+
  scale_shape_manual(values=c(15, 16,18,17))+
  scale_fill_brewer(palette="Dark2", guide=FALSE)+
  labs(x = NULL, y = 'Ground Beetle Diversity\n(Simpson Index)')+
  scale_x_discrete(name =NULL, 
                   limits = c("deciduousForest","mixedForest","evergreenForest"), 
                   labels = NULL)+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1))
p2


p3 <- ggplot(vegB.df, aes(x=nlcdClass, y=log(density), fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = vegB.df, aes(shape = siteID), size = 2, alpha = 0.5, width = 0.25, show.legend = F)+
  scale_shape_manual(values=c(15, 16,18,17))+
  scale_fill_brewer(palette="Dark2", guide=FALSE)+
  labs(x = 'Forest Cover Class', y = 'Ground Beetle Density (log)\n(Individuals/Cup)')+
  scale_x_discrete(name ="NLCD Forest Cover Class", 
                   limits = c("deciduousForest","mixedForest","evergreenForest"), 
                   labels = c("Deciduous", "Mixed","Evergreen"))+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  annotate(geom="text", x=1, y=1.78, label="a",
           color="black",alpha = 0.70, size=8)+
  annotate(geom="text", x=2, y=1.78, label="ab",
           color="black",alpha = 0.70, size=8)+
  annotate(geom="text", x=3, y=1.78, label="b",
           color="black",alpha = 0.70, size=8)+
  scale_y_continuous(breaks=c(-1.0,-.5,0,.5,1.0,1.5))
p3


#cow plot to stack them 

nlcdgrid<-plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 20, ncol = 1)
nlcdgrid
pdf("/Users/yourfilepath/Desktop/Model_Graphs/nlcd_grid.pdf", width = 7, height = 21)
plot(nlcdgrid)
dev.off()
}


#__________________________________________
# stats and graphs----
#__________________________________________

if(TRUE){
veg.df<-allVeg.df
beetle.df<-bind_rows(HARV.df,TALL.df)
  
  
#clean data frams 
vegB.df<-beetle.df%>%
  left_join(., veg.df, by = "plotID")%>%
  filter(., !siteID=="JERC",
         !siteID=="BART",
         !plotID=="TALL_003",
         !plotID=="TALL_005")%>%  #plots without veg data 
  rename(nlcdClass = nlcdClass.x)%>%
  mutate(PerAM_BA=AM_BA/totalBA*100,
         PerAM_ST=AM_ST/totalStems*100)

#__________________________________________
# 3 variables ~ nlcdClass
#__________________________________________

m <- nlme::lme(ShannonBeetle ~ nlcdClass, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m));Anova(m)
m <- nlme::lme(SimpsonBeetle ~ nlcdClass,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m));Anova(m)
m <- nlme::lme(density ~ nlcdClass, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m));Anova(m);lsmeans(m, pairwise~nlcdClass, adjust="tukey")   

#__________________________________________
# 3 variables ~ Tdiv, %EG, %ECM
#__________________________________________

m <- nlme::lme(ShannonBeetle ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(SimpsonBeetle ~ Sim_BA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

m <- nlme::lme(ShannonBeetle ~ Shan_ST + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(SimpsonBeetle ~ Sim_ST + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))


#__________________________________________
# 3 variables ~ Tden, %EG, %ECM
#__________________________________________

m <- nlme::lme(ShannonBeetle ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(SimpsonBeetle ~ totalBA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

#__________________________________________
#diversity begets diversity, density begets density
#__________________________________________

m <- nlme::lme(ShannonBeetle ~Shan_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(SimpsonBeetle ~ Sim_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ totalBA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

#__________________________________________
#diversity predicted by tree density
#__________________________________________
m <- nlme::lme(density ~ totalStems,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ totalBA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(ShannonBeetle ~totalBA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(ShannonBeetle ~ totalStems,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

m <- nlme::lme(SimpsonBeetle ~totalBA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(SimpsonBeetle ~ totalStems,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
#__________________________________________
#Top 5 
#m<-nlme::lme(Species1 ~ PerEG_BA, random = ~1| siteID/plotID,data= ?);summary(m);shapiro.test(resid(m))

#Species2~ PerEG_BA
#Species3~ PerEG_BA
#Species4~ PerEG_BA
#Species5~ PerEG_BA
#__________________________________________


#__________________________________________
#LOG versions----
#__________________________________________
#
#
#
#
#
#
#
#
#
#

#__________________________________________
#Log Transformed BeetleDiversity ~ Percent Evergreen
#__________________________________________

vegB.df$log<-log(vegB.df$ShannonBeetle)
vegB.df<-vegB.df%>%
  filter(!log == "-Inf") #sometimes one data point comes back -Inf

#__________________________________________
# 3 variables ~ nlcdClass
#__________________________________________

m <- nlme::lme(log(ShannonBeetle) ~ nlcdClass, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(SimpsonBeetle) ~ nlcdClass,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(density) ~ nlcdClass, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

#__________________________________________
# 3 variables ~ Tdiv, %EG, %ECM
#__________________________________________

m <- nlme::lme(log(ShannonBeetle) ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(SimpsonBeetle) ~ Sim_BA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(density) ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))


#__________________________________________
# 3 variables ~ Tden, %EG, %ECM
#__________________________________________

m <- nlme::lme(log(ShannonBeetle) ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(SimpsonBeetle) ~ totalBA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(density) ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

m <- nlme::lme(log(density) ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

#__________________________________________
#diversity begets diversity, density begets density
#__________________________________________

m <- nlme::lme(log(ShannonBeetle) ~Shan_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(SimpsonBeetle) ~ Sim_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(density) ~ totalBA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
#__________________________________________
#log version diversity predicted by tree density
#__________________________________________
m <- nlme::lme(log(density) ~ totalStems,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(density) ~ totalBA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(ShannonBeetle) ~totalBA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(ShannonBeetle) ~ totalStems,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

m <- nlme::lme(log(SimpsonBeetle) ~totalBA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(log(SimpsonBeetle) ~ totalStems,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
#__________________________________________


#__________________________________________
# Figure 2: ShanBEETLE Y ~ X %EG----
#__________________________________________

p<-ggplot(vegB.df,aes(x=PerEG_BA, y=ShannonBeetle))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape = siteID))+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance\n of Evergreen Trees (%)', y = 'Ground Beetle Diversity\n(Shannon Index)', color='Site')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p

pdf("/Users/yourfilepath/Desktop/Model_Graphs/ShanB_PerEF.pdf", width = 7, height = 7)
plot(p)
dev.off()

#__________________________________________
# Figure 3: DenBEETLE Y ~ X %EG----
#__________________________________________

p<-ggplot(vegB.df,aes(x=PerEG_BA, y=density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape = siteID))+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance of Evergreen Trees (%)', y = 'Ground Beetle Density\n(Individuals/Cup)', color='Site')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  guides(color=FALSE)
p

pdf("/Users/yourfilepath/Desktop/Model_Graphs/Den_PerEG.pdf", width = 7, height = 7)
plot(p)
dev.off()

#testing
p<-ggplot(vegB.df,aes(x=PerEG_BA, y=density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape = siteID))+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance\n of Evergreen Trees (%)', y = 'Ground Beetle Density (log)\n(Individuals/Cup)', color='Site')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p
pdf("/Users/yourfilepath/Desktop/Model_Graphs/Den_PerEG.pdf", width = 7, height = 7)
plot(p)
dev.off()

#__________________________________________
# Plot - perEG ~ Beetle diversity----
#__________________________________________
#simpsonB ~ PerEG_BA
p<-ggplot(vegB.df,aes(x=PerEG_BA, y=ShannonBeetle))+
  geom_smooth(method = 'lm', formula = 'y ~ x') + 
  geom_point(aes(color = plotID))+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance of Evergreen Trees (%)', y = 'Ground Beetle Diversity\n(Shannon Index)', color='Plot')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=16), 
        axis.title.y=element_text(size=16), 
        axis.text.x=element_text(size=12), 
        axis.text.y=element_text(size=12))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())
p

pdf("/Users/yourfilepath/Desktop/Model_Graphs/Div_PerEGallplot.pdf", width = 7, height = 7)
plot(p)
dev.off()
}

#__________________________________________
#  species abundance figures and stats----
#__________________________________________
if(TRUE){

  #harv
  speciesDensityPlotYr_HARV.df$sqrtCARGORdensity<-sqrt(speciesDensityPlotYr_HARV.df$CARGORdensity)
  hist(speciesDensityPlotYr_HARV.df$sqrtCARGORdensity)  
  
  speciesDensityPlotYr_HARV.df$sqrtSYNIMPdensity<-sqrt(speciesDensityPlotYr_HARV.df$SYNIMPdensity)
  hist(speciesDensityPlotYr_HARV.df$SPHCANdensity) 
  
  
  m <- nlme::lme(sqrtCARGORdensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_HARV.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(SYNIMPdensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_HARV.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(SPHSTE3density ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_HARV.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(PTETRI3density ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_HARV.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(SPHCANdensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_HARV.df);summary(m); shapiro.test(resid(m))
#tall
  speciesDensityPlotYr_TALL.df$sqrtCYCCON2density<-sqrt(speciesDensityPlotYr_TALL.df$CYCCON2density)
  hist(speciesDensityPlotYr_TALL.df$sqrtCYCCON2density)  
  
  speciesDensityPlotYr_TALL.df$logDICDIL5density<-log(speciesDensityPlotYr_TALL.df$DICDIL5density+1)
  hist(speciesDensityPlotYr_TALL.df$logDICDIL5density)  
  
  m <- nlme::lme(sqrtCYCCON2density ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(ANIMER2density ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(SELOPAdensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(ANIHAPdensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(logDICDIL5density ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  
  #just to see if EG alone had any effects
  m <- nlme::lme(sqrtCYCCON2density ~ PerEG_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(ANIMER2density ~ PerEG_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(SELOPAdensity ~ PerEG_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(ANIHAPdensity ~ PerEG_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))
  m <- nlme::lme(logDICDIL5density ~ PerEG_BA, random = ~1|plotID, 
                 data = speciesDensityPlotYr_TALL.df);summary(m); shapiro.test(resid(m))

  
p1<-ggplot(speciesDensityPlotYr_HARV.df,aes(x=PerEG_BA, y=CARGORdensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Carabus goryi', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)+
  ggtitle("HARV")+
  theme(title=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))
p1

p2<-ggplot(speciesDensityPlotYr_HARV.df,aes(x=PerEG_BA, y=SYNIMPdensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Synuchus impunctatus', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p2


p3<-ggplot(speciesDensityPlotYr_HARV.df,aes(x=PerEG_BA, y=SPHSTE3density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Sphaeroderus stenostomus', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p3

p4<-ggplot(speciesDensityPlotYr_HARV.df,aes(x=PerEG_BA, y=PTETRI3density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Pterostichus tristis', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p4

p5<-ggplot(speciesDensityPlotYr_HARV.df,aes(x=PerEG_BA, y=SPHCANdensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance\n of Evergreen Trees (%)', y = 'Sphaeroderus canadensis', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p5

p6<-ggplot(speciesDensityPlotYr_TALL.df,aes(x=PerEG_BA, y=CYCCON2density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Cyclotrachelus convivus', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)+
  ggtitle("TALL")+
  theme(title=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))
p6

p7<-ggplot(speciesDensityPlotYr_TALL.df,aes(x=PerEG_BA, y=ANIMER2density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Anisodactylus merula', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p7

p8<-ggplot(speciesDensityPlotYr_TALL.df,aes(x=PerEG_BA, y=SELOPAdensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Selenophorus opalinus', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p8

p9<-ggplot(speciesDensityPlotYr_TALL.df,aes(x=PerEG_BA, y=ANIHAPdensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Anisodactylus haplomus', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p9

p10<-ggplot(speciesDensityPlotYr_TALL.df,aes(x=PerEG_BA, y=DICDIL5density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance\n of Evergreen Trees (%)', y = 'Dicaelus dilatatus', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=20), 
        axis.title.y=element_text(size=20), 
        axis.text.x=element_text(size=18), 
        axis.text.y=element_text(size=18))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  labs(shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="right",
    legend.direction="vertical")+
  guides(color=FALSE)
p10
#cow plot to stack them 

#most abundant species plot----

Harv.plot<-plot_grid(p1,p2,p3,p4,p5, ncol = 1)
Harv.plot
Tall.plot<-plot_grid(p6,p7,p8,p9,p10, ncol = 1)
speciesgrid<-plot_grid(Harv.plot, Tall.plot, ncol = 2)
speciesgrid
pdf("/Users/yourfilepath/Desktop/Model_Graphs/species_grid.pdf", width = 10, height = 21)
plot(speciesgrid)
dev.off()


}

vegB.df
write_csv(vegB.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/vegB.df.csv")

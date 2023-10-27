#rm(list=ls(all=TRUE))
#clean work area

#__________________________________________
# library----
#__________________________________________

library(car)       #Anova (not anova)
library(RColorBrewer)
library(neonUtilities)
library(vegan)
library(plotrix) #std.error function
library(lubridate)
library(cowplot)
library(emmeans) #lsmeans to compare 
library(tidyverse)
library(dplyr)

# Set global option to NOT convert all character variables to factors
options(stringsAsFactors=F)

#__________________________________________
# Read in relevant raw data sets----
#__________________________________________

#PlantHARV <-loadByProduct(dpID="DP1.10098.001", site="HARV", #This is the file path to download from NEON directly, but I've saved the files in our GitHib to make it easier
 #                         package="expanded", check.size=T)

#BeetleHARV<- loadByProduct(dpID="DP1.10022.001", site="HARV", 
 #                          package="expanded", check.size=T) 

#PlantBART <-loadByProduct(dpID="DP1.10098.001", site="BART", 
 #                         package="expanded", check.size=T)

#BeetleBART<- loadByProduct(dpID="DP1.10022.001", site="BART", 
 #                          package="expanded", check.size=T) 

#PlantTALL <-loadByProduct(dpID="DP1.10098.001", site="TALL", 
 #                         package="expanded", check.size=T)

#BeetleTALL<- loadByProduct(dpID="DP1.10022.001", site="TALL", 
 #                          package="expanded", check.size=T) 

#BeetleJERC<- loadByProduct(dpID="DP1.10022.001", site="JERC", 
 #                          package="expanded", check.size=T) 
#no plant for JERC

#Load NEONdata.RData directly from GitHub Repository into your Environment

Mycorrhizal.df<-read.csv ("plantSpecies.csv", header=TRUE) #this is also already in the GitHub Repository

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
  
  mappingandtagging_HARV.df<-PlantHARV$vst_mappingandtagging%>%                #mapped and tagged trees with the species info
    filter(plotID %in% plots_HARV.df$plotID)%>%                                #filter by the plots that have beetle samples
    select(eventID,                                                            #remove unwanted columns from df
           individualID,                                                       
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
    filter(!plotID=="HARV_025")%>%                        #removed plot 025 because there were only 3 cups collected
    select(plotID,
           individualID,
           stemDiameter)%>%
    left_join(., perplotperyear_HARV.df, by = "plotID")%>%                     #join in sampling area
    left_join(., mappingandtagging_HARV.df, by = "individualID")%>%            #add species ID and names
    left_join(., Mycorrhizal.df, by = "taxonID")%>%                            #fungal associations and leaf habit
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
  species.df<-BasalAreaHARV.df[,10:18]                                         #CHECK EACH SITE:for vegan package, subset species into columns for diversity
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
  speciesStems.df<-StemsHARV.df[,10:18]                                        #CHECK SITE: Harv has 10 species, so subset columns 10-20
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
    mutate(stemCount=(replace(basalAreaCm2, values = 1)))%>%                   #create an individual count
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
    mutate(relAbunBA = basalAreaCm2/sum(SiteStatsHARV.df$basalAreaCm2)*100,    #calculate the relative abundance of a species at HARV
           relAbunSt = stems/sum(SiteStatsHARV.df$stems)*100)                  #relative abundance = 1 species/ total species at site
}

#__________________________________________
# Descriptive Statistics PLOT level----
#__________________________________________
if(TRUE){
  PlotStatsHARV.df<-BasalAreaHARV.df%>%                                        #create a df with descriptive statistics at the plot level
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
    filter(eventID=="vst_TALL_2015")%>%                                        #retain only year 2015 for Tall
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

expertTax_HARV <- BeetleHARV$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate < "2019-12-31",                                           # incomplete sampling year
         collectDate > "2013-12-31",                                           # incomplete sampling year
         !plotID=="HARV_025",                                                  # incomplete sampling year
         !setDate=="2014-09-12",                        
         family=="Carabidae")%>%                                               # removing species that are not ground beetles
  mutate(taxonID = replace(taxonID, taxonID == 'APELUC', 'APELUC2'),           # changing subspecies to species
         taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%
  transform(scientificName=paste(genus, specificEpithet))%>%                   # renaming all species to be only genus and specific Epithet (helpful later)
  select(individualID,
         taxonID)
#write_csv(expertTax_HARV, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/expertTax_HARV.csv")
#quickly write a version with scientific name

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
           subsampleID,
           taxonID)%>%
    left_join(., expertTax_HARV, by = "individualID")%>%                         # joining expertly identified ground beetles to the pinned ground beetles
    dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%                 # changing the non-expert ID's to expert ID'
    select(subsampleID,
           taxonID.y)%>%
    distinct()

#make a list of para and expert IDed subsamples that contain only 1 species
paraTaxUnique<-paraTax_HARV%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n==1)

#filter subsamples containing more than 1 species out of paraTax table
paraTax_HARV_unique<-paraTax_HARV %>%
  filter(subsampleID %in% paraTaxUnique$subsampleID)

sort_HARV.df <- BeetleHARV$bet_sorting%>%
  filter(collectDate < "2019-12-31",                                           # incomplete sampling year
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
sum(sort_HARV.df$individualCount) #3823

#join subsamples with only 1 species to sort table
sort_HARV.df<-sort_HARV.df%>%
  left_join(paraTax_HARV_unique, by="subsampleID")

#make a list of para and expert IDed subsamples that contain more than 1 species
paraTaxMulti<-paraTax_HARV%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n>1)

#keep subsamples containing more than 1 species in paraTax table
paraTax_HARV_multi<-paraTax_HARV %>%
  filter(subsampleID %in% paraTaxMulti$subsampleID)

#join sort data to paraTax multi
paraTax_HARV_multi<-paraTax_HARV_multi%>%
  left_join(sort_HARV.df, by="subsampleID")%>%
  select(subsampleID, 
         taxonID.y.x, 
         individualCount, 
         sampleID)%>%
  left_join(paraTaxMulti,by="subsampleID")%>%
  mutate (ind=individualCount/n)       #assume expertly IDed species composition is reflective of whole subsample (other fewer individuals expertly IDed than there are individuals in a subsample)

data_wide <- spread(paraTax_HARV_multi, taxonID.y.x, ind)%>%
  select(!individualCount)%>%
  select(!n)

#remove subsamples with multiple species from sort table
sort_HARV.df<-sort_HARV.df %>%
  filter(!subsampleID %in% paraTaxMulti$subsampleID)

#if there are no para or expert ID's, assume sort ID are correct
sort_HARV.df<-sort_HARV.df %>%
  dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%
  select(!taxonID)
  
data_wide1<-spread(sort_HARV.df, taxonID.y, individualCount)
    
#put two dataframes of species and counts in subsamples together

Example <- full_join(data_wide, data_wide1, by = "subsampleID")%>% #once again 2028 subsamples
  dplyr::mutate(sampleID.x = coalesce(sampleID.y, sampleID.x))%>% 
  replace(is.na(.), 0)%>%
  mutate(APELUC2=APELUC2.x+APELUC2.y,                             #must be a better way to do this
         CARGOR=CARGOR.x+CARGOR.y, 
         CARSER=CARSER.x+CARSER.y, 
         CYMCRI=CYMCRI.x+CYMCRI.y, 
         CYMLIM=CYMLIM.x+CYMLIM.y, 
         CYMPLA2=CYMPLA2.x+CYMPLA2.y, 
         DICPOL=DICPOL.x+DICPOL.y,
         MYACYA=MYACYA.x+MYACYA.y, 
         NOTAEN=NOTAEN.x+NOTAEN.y, 
         PLADEC=PLADEC.x+PLADEC.y, 
         PTEADO=PTEADO.x+PTEADO.y, 
         PTECOR1=PTECOR1.x+PTECOR1.y, 
         PTEMUT2=PTEMUT2.x+PTEMUT2.y, 
         PTEPEN=PTEPEN.x+PTEPEN.y, 
         PTEROS=PTEROS.x+PTEROS.y, 
         PTETRI3=PTETRI3.x+PTETRI3.y, 
         SPHCAN=SPHCAN.x+SPHCAN.y, 
         SPHSTE3=SPHSTE3.x+SPHSTE3.y, 
         SYNIMP=SYNIMP.x+SYNIMP.y, 
         TRIAUT=TRIAUT.x+TRIAUT.y)%>%
  select(!APELUC2.x)%>%
  select(!APELUC2.y)%>%
  select(!CARGOR.y)%>%
  select(!CARGOR.x)%>%
  select(!CARSER.x)%>%
  select(!CARSER.y)%>%
  select(!CYMCRI.x)%>%
  select(!CYMCRI.y)%>%
  select(!CYMLIM.x)%>%
  select(!CYMLIM.y)%>%
  select(!CYMPLA2.x)%>%
  select(!CYMPLA2.y)%>%
  select(!DICPOL.x)%>%
  select(!DICPOL.y)%>%
  select(!MYACYA.x)%>%
  select(!MYACYA.y)%>%
  select(!NOTAEN.x)%>%
  select(!NOTAEN.y)%>%
  select(!PLADEC.x)%>%
  select(!PLADEC.y)%>%
  select(!PTEADO.x)%>%
  select(!PTEADO.y)%>%
  select(!PTECOR1.x)%>%
  select(!PTECOR1.y)%>%
  select(!PTEMUT2.x)%>%
  select(!PTEMUT2.y)%>%
  select(!PTEPEN.x)%>%
  select(!PTEPEN.y)%>%
  select(!PTEROS.x)%>%
  select(!PTEROS.y)%>%
  select(!PTETRI3.x)%>%
  select(!PTETRI3.y)%>%
  select(!SPHCAN.x)%>%
  select(!SPHCAN.y)%>%
  select(!SPHSTE3.x)%>%
  select(!SPHSTE3.y)%>%
  select(!SYNIMP.x)%>%
  select(!SYNIMP.y)%>%
  select(!TRIAUT.x)%>%
  select(!TRIAUT.y)%>%
  select(!sampleID.y)%>%
  rename(sampleID=sampleID.x)

df<-Example %>% 
  select(!subsampleID)%>%
  group_by(sampleID) %>% 
  summarise(across(everything(), sum))
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
         collectDate,
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

field_try.df<-field_HARV.df%>%
  left_join(df, by="sampleID")%>%
  mutate(cupStatus = replace_na("Ok"))%>%  #recode NA cup and lid status to Ok
  mutate(lidStatus = replace_na("Ok"))%>%
  filter (cupStatus!="Disturbed")%>%
  filter (cupStatus!="Missing")%>%
  filter(lidStatus!="Disturbed")%>%
  filter(lidStatus!="Missing")%>%          #remove cups that were disturbed or missing or had disturbed or missing lids; assume NA OK
  replace(is.na(.), 0)%>%                  #replace NAs with 0
  mutate(collectDate = ymd(collectDate),   #define collect date as a date
         collectDate = as_date(collectDate))%>%
  mutate (year=year(collectDate))%>%       #extract year
  transform(., plotIDyear=paste(plotID, year))%>% #create plot_year variable
  mutate(cups=1)%>%
  select(plotID, ACUHYD:cups)

field_try2<-field_try.df%>%    #summarize # of cups sampled and number of individuals of each species by plot_year
  select(!plotID)%>%
  select(!year)%>%
  group_by(plotIDyear)%>%
  summarise(across(everything(), sum))%>%
  rowwise() %>%                #sum total individuals by plot_year
  mutate(totalIndividuals = sum(c_across(ACUHYD:TRIAUT)))%>%   
  mutate(density=totalIndividuals/cups)  #calculate yearly density for each plot

#calculate diversity for each year for each plot
species.df<-field_try2[,2:34]                                            # species only df for vegan diversity analysis - right amount of species 2388
field_try2$ShannonBeetle<-diversity(species.df, index = "shannon")       # shannon diversity index
field_try2$SimpsonBeetle<-diversity(species.df, index = "simpson")       # simpson diversity index
field_try2$spNum<-specnumber(species.df)

#should diversity indices be weighted somehow by the number of cups of individuals caught?

sum(field_try2$totalIndividuals) #3823 individuals collected; 33 species
ave(field_try2$density) #1.90923 ind/cup on average
std.error(field_try2$density) #SE: 0.1888564
ave(field_try2$ShannonBeetle) #1.617419
std.error(field_try2$ShannonBeetle) #SE: 0.03778079
ave(field_try2$SimpsonBeetle) #0.7248093
std.error(field_try2$SimpsonBeetle) #SE: 0.01216914

as.data.frame(colSums(field_try2[,c(2:34)]))           #find most abundant species

sum(field_try2$CARGOR)/sum(field_try2$totalIndividuals)*100 #22.21%; 849 ind
sum(field_try2$SYNIMP)/sum(field_try2$totalIndividuals)*100 #34.01%; 1300 ind
sum(field_try2$SPHSTE3)/sum(field_try2$totalIndividuals)*100 #11.85%; 453 ind
sum(field_try2$PTETRI3)/sum(field_try2$totalIndividuals)*100 #6.34%; 242 ind
sum(field_try2$PTEPEN)/sum(field_try2$totalIndividuals)*100 #6.00%; 230 ind
sum(field_try2$SPHCAN)/sum(field_try2$totalIndividuals)*100 #4.64%; 178 ind

#species with <5 ind: ACUHYD, AMPINT, ANIRUN, CALFRI, CICSEX, GASHON, HARRUB, OLIPAR, POELUC (9 of 33 at HARV)
  
#calculate density and relative abundance of dominant species
field_try2<-field_try2%>%
  mutate(CARGOR_den=CARGOR/cups*100)%>%
  mutate(SYNIMP_den=SYNIMP/cups*100)%>%
  mutate(SPHSTE3_den=SPHSTE3/cups*100)%>%
  mutate(PTETRI3_den=PTETRI3/cups*100)%>%
  mutate(PTEPEN_den=PTEPEN/cups*100)%>%
  mutate(CARGOR_relabun=CARGOR/totalIndividuals*100)%>%
  mutate(SYNIMP_relabun=SYNIMP/totalIndividuals*100)%>%
  mutate(SPHSTE3_relabun=SPHSTE3/totalIndividuals*100)%>%
  mutate(PTETRI3_relabun=PTETRI3/totalIndividuals*100)%>%
  mutate(PTEPEN_relabun=PTEPEN/totalIndividuals*100)
  
#join tree data to beetle data
DF2 <- field_try.df %>%
  select(plotIDyear, plotID)%>%
  distinct()
HARV_data.df<-field_try2%>%
  left_join(DF2, by="plotIDyear")%>%
  left_join(allHARV.df, by="plotID")

#quick and dirty graphs and stats to check that same patterns hold as before
hist(HARV_data.df$density)
HARV_data.df$logDen<-log(HARV_data.df$density)
hist(HARV_data.df$logDen)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logDen)
abline(lm(HARV_data.df$logDen~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logDen~HARV_data.df$PerEG_BA)) #still sig

boxplot(HARV_data.df$logDen~HARV_data.df$nlcdClass)
summary(aov(HARV_data.df$logDen~HARV_data.df$nlcdClass)) #still sig

hist(HARV_data.df$ShannonBeetle)
plot(HARV_data.df$PerEG_BA, HARV_data.df$ShannonBeetle)
abline(lm(HARV_data.df$ShannonBeetle~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$ShannonBeetle~HARV_data.df$PerEG_BA)) #not sig

boxplot(HARV_data.df$ShannonBeetle~HARV_data.df$nlcdClass)
summary(aov(HARV_data.df$ShannonBeetle~HARV_data.df$nlcdClass)) #not sig

hist(HARV_data.df$SimpsonBeetle)
plot(HARV_data.df$PerEG_BA, HARV_data.df$SimpsonBeetle)
abline(lm(HARV_data.df$SimpsonBeetle~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$SimpsonBeetle~HARV_data.df$PerEG_BA)) #not sig

boxplot(HARV_data.df$SimpsonBeetle~HARV_data.df$nlcdClass)
summary(aov(HARV_data.df$SimpsonBeetle~HARV_data.df$nlcdClass)) #not sig

#checking out evergreen effects on abundant sp density and relative abundance
#SYNIMP
hist(HARV_data.df$SYNIMP_den)
HARV_data.df$logSYN<-log(HARV_data.df$SYNIMP_den+1)
hist(HARV_data.df$logSYN)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logSYN)
abline(lm(HARV_data.df$logSYN~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logSYN~HARV_data.df$PerEG_BA)) #marg sig

hist(HARV_data.df$SYNIMP_relabun)
HARV_data.df$sqrtSYN<-sqrt(HARV_data.df$SYNIMP_relabun)
hist(HARV_data.df$sqrtSYN)

plot(HARV_data.df$PerEG_BA, HARV_data.df$sqrtSYN)
abline(lm(HARV_data.df$sqrtSYN~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$sqrtSYN~HARV_data.df$PerEG_BA)) # sig

#CAROGR
hist(HARV_data.df$CARGOR_den)
HARV_data.df$logcCARGOR<-log(HARV_data.df$CARGOR_den+1)
hist(HARV_data.df$logcCARGOR)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logcCARGOR)
abline(lm(HARV_data.df$logcCARGOR~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logcCARGOR~HARV_data.df$PerEG_BA)) #sig

hist(HARV_data.df$CARGOR_relabun)
HARV_data.df$logCARGOR_relabun<-log(HARV_data.df$CARGOR_relabun+1)
hist(HARV_data.df$logCARGOR_relabun)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logCARGOR_relabun)
abline(lm(HARV_data.df$logCARGOR_relabun~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logCARGOR_relabun~HARV_data.df$PerEG_BA)) # marg sig

#PTETRI
hist(HARV_data.df$PTETRI3_den)
HARV_data.df$logPTETRI3<-log(HARV_data.df$PTETRI3_den+1)
hist(HARV_data.df$logPTETRI3)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logPTETRI3)
abline(lm(HARV_data.df$logPTETRI3~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logPTETRI3~HARV_data.df$PerEG_BA)) #not sig

hist(HARV_data.df$PTETRI3_relabun)
HARV_data.df$logPTETRI3_relabun<-log(HARV_data.df$PTETRI3_relabun+1)
hist(HARV_data.df$logPTETRI3_relabun)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logPTETRI3_relabun)
abline(lm(HARV_data.df$logPTETRI3_relabun~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logPTETRI3_relabun~HARV_data.df$PerEG_BA)) # marg sig

#PTEPEN
hist(HARV_data.df$PTEPEN_den)
HARV_data.df$logPTEPEN_den<-log(HARV_data.df$PTEPEN_den+1)
hist(HARV_data.df$logPTEPEN_den)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logPTEPEN_den)
abline(lm(HARV_data.df$logPTEPEN_den~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logPTEPEN_den~HARV_data.df$PerEG_BA)) #not sig

hist(HARV_data.df$PTEPEN_relabun)
HARV_data.df$logPTEPEN_relabun<-log(HARV_data.df$PTEPEN_relabun+1)
hist(HARV_data.df$logPTEPEN_relabun)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logPTEPEN_relabun)
abline(lm(HARV_data.df$logPTEPEN_relabun~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logPTEPEN_relabun~HARV_data.df$PerEG_BA)) # not sig

#SPHSTE
hist(HARV_data.df$SPHSTE3_den)
HARV_data.df$logSPHSTE3_den<-log(HARV_data.df$SPHSTE3_den+1)
hist(HARV_data.df$logSPHSTE3_den)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logSPHSTE3_den)
abline(lm(HARV_data.df$logSPHSTE3_den~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logSPHSTE3_den~HARV_data.df$PerEG_BA)) #not sig

hist(HARV_data.df$SPHSTE3_relabun)
HARV_data.df$logSPHSTE3_relabun<-log(HARV_data.df$SPHSTE3_relabun+1)
hist(HARV_data.df$logSPHSTE3_relabun)

plot(HARV_data.df$PerEG_BA, HARV_data.df$logSPHSTE3_relabun)
abline(lm(HARV_data.df$logSPHSTE3_relabun~HARV_data.df$PerEG_BA))
summary(lm(HARV_data.df$logSPHSTE3_relabun~HARV_data.df$PerEG_BA)) # not sig

#__________________________________________
# TALL                              ----
#__________________________________________

expertTax_TALL <- BeetleTALL$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2015-12-31",                                           #remove 2014 and 2015: incomplete data (started sampling mid-May, but late March in other years)  
         family=="Carabidae")%>%                                               #remove samples that are not ground beetles
  filter(!plotID=="TALL_003")%>%
  filter(!plotID=="TALL_005")%>%                                               #removed plots 003 and 005 because they only have 1 and 3 trees larger than 10cm DBH
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),          #change subspecies to species
         taxonID = replace(taxonID, taxonID == 'DICFUR2', 'DICFUR'),
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'),
         taxonID = replace(taxonID, taxonID == 'CICPUN1', 'CICPUN'))%>%
  transform(scientificName=paste(genus, specificEpithet))%>%                   #this code will convert subspecies names to the species name 
  select(individualID,
         taxonID)
#write_csv(expertTax_TALL, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/expertTax_TALL.csv")
#saved a version with the scientific name

paraTax_TALL <- BeetleTALL$bet_parataxonomistID%>%
  filter(collectDate < "2019-12-31",                                         #remove 2020 data: incomplete data 
         collectDate > "2015-12-31")%>%
  filter(!plotID=="TALL_003")%>%
  filter(!plotID=="TALL_005")%>% 
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),        #change subspecies to species
         taxonID = replace(taxonID, taxonID == 'DICFUR2', 'DICFUR'),
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'),
         taxonID = replace(taxonID, taxonID == 'CICPUN1', 'CICPUN'))%>%
  select(individualID,
         subsampleID,
         taxonID)%>%
  left_join(., expertTax_TALL, by = "individualID")%>%                         # joining expertly identified ground beetles to the pinned ground beetles
  dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%                 # changing the non-expert ID's to expert ID'
  select(subsampleID,
         taxonID.y)%>%
  distinct()

#make a list of para and expert IDed subsamples that contain only 1 species
paraTaxUnique<-paraTax_TALL%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n==1)

#filter subsamples containing more than 1 species out of paraTax table
paraTax_TALL_unique<-paraTax_TALL %>%
  filter(subsampleID %in% paraTaxUnique$subsampleID)

sort_TALL.df <- BeetleTALL$bet_sorting%>%
  filter(collectDate < "2019-12-31",   
         collectDate > "2015-12-31",                                           #remove 2020 data: incomplete data 
         sampleType == "carabid" | sampleType =="other carabid")%>%            #remove samples that are not ground beetles
  filter(!plotID=="TALL_003")%>%
  filter(!plotID=="TALL_005")%>% 
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),          #aggregate taxa to the species-level
         taxonID = replace(taxonID, taxonID == 'DICFUR2', 'DICFUR'),
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'),
         taxonID = replace(taxonID, taxonID == 'CICPUN1', 'CICPUN'))%>%
  select(subsampleID,
         taxonID,
         individualCount,
         sampleID)
sum(sort_TALL.df$individualCount)                                              #470 individuals in total 2016-2019

#join subsamples with only 1 species to sort table
sort_TALL.df<-sort_TALL.df%>%
  left_join(paraTax_TALL_unique, by="subsampleID")

#make a list of para and expert IDed subsamples that contain more than 1 species
paraTaxMulti<-paraTax_TALL%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n>1)

#keep subsamples containing more than 1 species in paraTax table
paraTax_TALL_multi<-paraTax_TALL %>%
  filter(subsampleID %in% paraTaxMulti$subsampleID)

#join sort data to paraTax multi
paraTax_TALL_multi<-paraTax_TALL_multi%>%
  left_join(sort_TALL.df, by="subsampleID")%>%
  select(subsampleID, 
         taxonID.y.x, 
         individualCount, 
         sampleID)%>%
  left_join(paraTaxMulti,by="subsampleID")%>%
  mutate (ind=individualCount/n)       #assume expertly IDed species composition is reflective of whole subsample (other fewer individuals expertly IDed than there are individuals in a subsample)

data_wide <- spread(paraTax_TALL_multi, taxonID.y.x, ind)%>%
  select(!individualCount)%>%
  select(!n)

#remove subsamples with multiple species from sort table
sort_TALL.df<-sort_TALL.df %>%
  filter(!subsampleID %in% paraTaxMulti$subsampleID)

#if there are no para or expert ID's, assume sort ID are correct
sort_TALL.df<-sort_TALL.df %>%
  dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%
  select(!taxonID)

data_wide1<-spread(sort_TALL.df, taxonID.y, individualCount)

#put two dataframes of species and counts in subsamples together

Example <- full_join(data_wide, data_wide1, by = "subsampleID")%>% #once again 2028 subsamples
  dplyr::mutate(sampleID.x = coalesce(sampleID.y, sampleID.x))%>% 
  replace(is.na(.), 0)%>%
  mutate(ANIHAP=ANIHAP.x+ANIHAP.y,                             #must be a better way to do this
         ANIMER2=ANIMER2.x+ANIMER2.y, 
         PASDEP=PASDEP.x+PASDEP.y, 
         PASPUN=PASPUN.x+PASPUN.y)%>%
  select(!ANIHAP.x)%>%
  select(!ANIHAP.y)%>%
  select(!ANIMER2.x)%>%
  select(!ANIMER2.y)%>%
  select(!PASDEP.x)%>%
  select(!PASDEP.y)%>%
  select(!PASPUN.x)%>%
  select(!PASPUN.y)%>%
  select(!sampleID.y)%>%
  rename(sampleID=sampleID.x)

df<-Example %>% 
  select(!subsampleID)%>%
  group_by(sampleID) %>% 
  summarise(across(everything(), sum))

field_TALL.df <- BeetleTALL$bet_fielddata%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2015-12-31",
         !sampleCollected=="N") %>%                                            #removing missing samples 
  filter(!plotID=="TALL_003")%>%
  filter(!plotID=="TALL_005")%>% 
  select(plotID,
         trapID,
         nlcdClass,
         collectDate,
         eventID,
         sampleID,
         cupStatus,
         lidStatus)%>%
  transform(., plotIDeventID=paste(plotID, eventID))%>%
  distinct()

field_try.df<-field_TALL.df%>%
  left_join(df, by="sampleID")%>%
  mutate(cupStatus = replace_na("Ok"))%>%  #recode NA cup and lid status to Ok
  mutate(lidStatus = replace_na("Ok"))%>%
  filter (cupStatus!="Disturbed")%>%
  filter (cupStatus!="Missing")%>%
  filter(lidStatus!="Disturbed")%>%
  filter(lidStatus!="Missing")%>%          #remove cups that were disturbed or missing or had disturbed or missing lids; assume NA OK
  replace(is.na(.), 0)%>%                  #replace NAs with 0
  mutate(collectDate = ymd(collectDate),   #define collect date as a date
         collectDate = as_date(collectDate))%>%
  mutate (year=year(collectDate))%>%       #extract year
  transform(., plotIDyear=paste(plotID, year))%>% #create plot_year variable
  mutate(cups=1)%>%
  select(plotID, AGOCON:cups)

field_try2<-field_try.df%>%    #summarize # of cups sampled and number of individuals of each species by plot_year
  select(!plotID)%>%
  select(!year)%>%
  group_by(plotIDyear)%>%
  summarise(across(everything(), sum))%>%
  rowwise() %>%                #sum total individuals by plot_year
  mutate(totalIndividuals = sum(c_across(AGOCON:PASPUN)))%>%   
  mutate(density=totalIndividuals/cups)  #calculate yearly density for each plot

#calculate diversity for each year for each plot
species.df<-field_try2[,2:34]                                            # species only df for vegan diversity analysis - right amount of species 2388
field_try2$ShannonBeetle<-diversity(species.df, index = "shannon")       # shannon diversity index
field_try2$SimpsonBeetle<-diversity(species.df, index = "simpson")       # simpson diversity index
field_try2$spNum<-specnumber(species.df)

#should diversity indices be weighted somehow by the number of cups of individuals caught?

sum(field_try2$totalIndividuals) #470 individuals collected; 33 species
ave(field_try2$density) #0.336069 ind/cup on average
std.error(field_try2$density) #SE: 0.04528136
ave(field_try2$ShannonBeetle) #1.283471 
std.error(field_try2$ShannonBeetle) #SE: 0.08528649
ave(field_try2$SimpsonBeetle) #0.6249106 
std.error(field_try2$SimpsonBeetle) #SE: 0.03445304

as.data.frame(colSums(field_try2[,c(2:34)]))           #find most abundant species

sum(field_try2$CYCCON2)/sum(field_try2$totalIndividuals)*100 #35.32%; 166 ind
sum(field_try2$DICDIL5)/sum(field_try2$totalIndividuals)*100 #8.51%; 40 ind
sum(field_try2$CYCFRE)/sum(field_try2$totalIndividuals)*100 #6.17%; 29 ind
sum(field_try2$ANIHAP)/sum(field_try2$totalIndividuals)*100 #5.85%; 27.5 ind
sum(field_try2$PASDEP)/sum(field_try2$totalIndividuals)*100 #5.32%; 25 ind

#species with <5 ind: AGOCON, ANIRUS, APESIN, CHLEMA, DICELO, HARPEN, HELNIG, LOXCRE, NOTNOV, POLLAE, PTESP23, SCAUNI1, SELELL , SELGRA, SPHSTE1 (15 of 33 at HARV)

#calculate density and relative abundance of dominant species
field_try2<-field_try2%>%
  mutate(CYCCON2_den=CYCCON2/cups*100)%>%
  mutate(DICDIL5_den=DICDIL5/cups*100)%>%
  mutate(CYCFRE_den=CYCFRE/cups*100)%>%
  mutate(ANIHAP_den=ANIHAP/cups*100)%>%
  mutate(PASDEP_den=PASDEP/cups*100)%>%
  mutate(CYCCON2_relabun=CYCCON2/totalIndividuals*100)%>%
  mutate(DICDIL5_relabun=DICDIL5/totalIndividuals*100)%>%
  mutate(CYCFRE_relabun=CYCFRE/totalIndividuals*100)%>%
  mutate(ANIHAP_relabun=ANIHAP/totalIndividuals*100)%>%
  mutate(PASDEP_relabun=PASDEP/totalIndividuals*100)

#join tree data to beetle data
DF2 <- field_try.df %>%
  select(plotIDyear, plotID)%>%
  distinct()
TALL_data.df<-field_try2%>%
  left_join(DF2, by="plotIDyear")%>%
  left_join(allTALL.df, by="plotID")

#quick and dirty graphs and stats to check that same patterns hold as before
hist(TALL_data.df$density)
TALL_data.df$sqrtDen<-sqrt(TALL_data.df$density)
hist(TALL_data.df$sqrtDen)

plot(TALL_data.df$PerEG_BA, TALL_data.df$sqrtDen)
abline(lm(TALL_data.df$sqrtDen~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$sqrtDen~TALL_data.df$PerEG_BA)) #still sig

boxplot(TALL_data.df$sqrtDen~TALL_data.df$nlcdClass)
summary(aov(TALL_data.df$sqrtDen~TALL_data.df$nlcdClass)) #not sig, but clear pattern

hist(TALL_data.df$ShannonBeetle)
plot(TALL_data.df$PerEG_BA, TALL_data.df$ShannonBeetle)
abline(lm(TALL_data.df$ShannonBeetle~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$ShannonBeetle~TALL_data.df$PerEG_BA)) #sig

boxplot(TALL_data.df$ShannonBeetle~TALL_data.df$nlcdClass)
summary(aov(TALL_data.df$ShannonBeetle~TALL_data.df$nlcdClass)) #not sig

hist(TALL_data.df$SimpsonBeetle)
plot(TALL_data.df$PerEG_BA, TALL_data.df$SimpsonBeetle)
abline(lm(TALL_data.df$SimpsonBeetle~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$SimpsonBeetle~TALL_data.df$PerEG_BA)) #not sig

boxplot(TALL_data.df$SimpsonBeetle~TALL_data.df$nlcdClass)
summary(aov(TALL_data.df$SimpsonBeetle~TALL_data.df$nlcdClass)) #not sig

#checking out evergreen effects on abundant sp density and relative abundance
#CYCCON2
hist(TALL_data.df$CYCCON2_den)
TALL_data.df$logCYCCON2<-log(TALL_data.df$CYCCON2_den+1)
hist(TALL_data.df$logCYCCON2)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logCYCCON2)
abline(lm(TALL_data.df$logCYCCON2~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logCYCCON2~TALL_data.df$PerEG_BA)) #sig

hist(TALL_data.df$CYCCON2_relabun)
TALL_data.df$sqrtCYCCON2<-sqrt(TALL_data.df$CYCCON2_relabun)
hist(TALL_data.df$sqrtCYCCON2)

plot(TALL_data.df$PerEG_BA, TALL_data.df$sqrtCYCCON2)
abline(lm(TALL_data.df$sqrtCYCCON2~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$sqrtCYCCON2~TALL_data.df$PerEG_BA)) # margsig

#DICDIL5
hist(TALL_data.df$DICDIL5_den)
TALL_data.df$logDICDIL5<-log(TALL_data.df$DICDIL5_den+1)
hist(TALL_data.df$logDICDIL5)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logDICDIL5)
abline(lm(TALL_data.df$logDICDIL5~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logDICDIL5~TALL_data.df$PerEG_BA)) #sig

hist(TALL_data.df$DICDIL5_relabun)
TALL_data.df$logDICDIL5_relabun<-log(TALL_data.df$DICDIL5_relabun+1)
hist(TALL_data.df$logDICDIL5_relabun)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logDICDIL5_relabun)
abline(lm(TALL_data.df$logDICDIL5_relabun~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logDICDIL5_relabun~TALL_data.df$PerEG_BA)) # not sig

#CYCFRE
hist(TALL_data.df$CYCFRE_den)
TALL_data.df$logCYCFRE<-log(TALL_data.df$CYCFRE_den+1)
hist(TALL_data.df$logCYCFRE)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logCYCFRE)
abline(lm(TALL_data.df$logCYCFRE~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logCYCFRE~TALL_data.df$PerEG_BA)) #not sig

hist(TALL_data.df$CYCFRE_relabun)
TALL_data.df$logCYCFRE_relabun<-log(TALL_data.df$CYCFRE_relabun+1)
hist(TALL_data.df$logCYCFRE_relabun)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logCYCFRE_relabun)
abline(lm(TALL_data.df$logCYCFRE_relabun~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logCYCFRE_relabun~TALL_data.df$PerEG_BA)) # not sig

#ANIHAP
hist(TALL_data.df$ANIHAP_den)
TALL_data.df$logANIHAP<-log(TALL_data.df$ANIHAP_den+1)
hist(TALL_data.df$logANIHAP)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logANIHAP)
abline(lm(TALL_data.df$logANIHAP~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logANIHAP~TALL_data.df$PerEG_BA)) #not sig

hist(TALL_data.df$ANIHAP_relabun)
TALL_data.df$logANIHAP_relabun<-log(TALL_data.df$ANIHAP_relabun+1)
hist(TALL_data.df$logANIHAP_relabun)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logANIHAP_relabun)
abline(lm(TALL_data.df$logANIHAP_relabun~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logANIHAP_relabun~TALL_data.df$PerEG_BA)) # not sig

#PASDEP
hist(TALL_data.df$PASDEP_den)
TALL_data.df$logPASDEP<-log(TALL_data.df$PASDEP_den+1)
hist(TALL_data.df$logPASDEP)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logPASDEP)
abline(lm(TALL_data.df$logPASDEP~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logPASDEP~TALL_data.df$PerEG_BA)) #not sig

hist(TALL_data.df$PASDEP_relabun)
TALL_data.df$logPASDEP_relabun<-log(TALL_data.df$PASDEP_relabun+1)
hist(TALL_data.df$logPASDEP_relabun)

plot(TALL_data.df$PerEG_BA, TALL_data.df$logPASDEP_relabun)
abline(lm(TALL_data.df$logPASDEP_relabun~TALL_data.df$PerEG_BA))
summary(lm(TALL_data.df$logPASDEP_relabun~TALL_data.df$PerEG_BA)) # not sig


#__________________________________________
# Read in relevant raw data sets JERC ----
#__________________________________________

expertTax_JERC <- BeetleJERC$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2014-12-31",                                           #2014 sampling began in mid-May
         family=="Carabidae")%>%                                               #remove samples that are not ground beetles
  filter(!plotID=="JERC_022")%>%
  filter(!plotID=="JERC_023")%>%                                               #remove plot 022 (only sample for 2 events), 023 (not eg, df, or mf)
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),          #change subspecies to species
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'))%>% 
  transform(scientificName=paste(genus, specificEpithet))%>%                   #this code will convert subspecies names to the species name 
  select(individualID,
         scientificName,
         taxonID)

paraTax_JERC <- BeetleJERC$bet_parataxonomistID%>%
  filter(collectDate < "2019-12-31",                                         #remove 2020 data: incomplete data 
         collectDate > "2014-12-31")%>%
  filter(!plotID=="JERC_022")%>%
  filter(!plotID=="JERC_023")%>%  
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),          #change subspecies to species
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'))%>% 
  select(individualID,
         subsampleID,
         taxonID)%>%
  left_join(., expertTax_JERC, by = "individualID")%>%                         # joining expertly identified ground beetles to the pinned ground beetles
  dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%                 # changing the non-expert ID's to expert ID'
  select(subsampleID,
         taxonID.y)%>%
  distinct()

#make a list of para and expert IDed subsamples that contain only 1 species
paraTaxUnique<-paraTax_JERC%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n==1)

#filter subsamples containing more than 1 species out of paraTax table
paraTax_JERC_unique<-paraTax_JERC %>%
  filter(subsampleID %in% paraTaxUnique$subsampleID)

sort_JERC.df <- BeetleJERC$bet_sorting%>%
  filter(collectDate < "2019-12-31",   
         collectDate > "2014-12-31",                                           #remove 2020 data: incomplete data 
         sampleType == "carabid" | sampleType =="other carabid")%>%            #remove samples that are not ground beetles
  filter(!plotID=="JERC_022")%>%
  filter(!plotID=="JERC_023")%>%  
  mutate(taxonID = replace(taxonID, taxonID == 'DICDIL1', 'DICDIL5'),          #change subspecies to species
         taxonID = replace(taxonID, taxonID == 'TETCAR1', 'TETCAR2'))%>% 
  select(subsampleID,
         taxonID,
         individualCount,
         sampleID)
sum(sort_JERC.df$individualCount)                                              #410 individuals in total 2015-2019

#join subsamples with only 1 species to sort table
sort_JERC.df<-sort_JERC.df%>%
  left_join(paraTax_JERC_unique, by="subsampleID")

#make a list of para and expert IDed subsamples that contain more than 1 species
paraTaxMulti<-paraTax_JERC%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n>1)

#keep subsamples containing more than 1 species in paraTax table
paraTax_JERC_multi<-paraTax_JERC %>%
  filter(subsampleID %in% paraTaxMulti$subsampleID)

#join sort data to paraTax multi
paraTax_JERC_multi<-paraTax_JERC_multi%>%
  left_join(sort_JERC.df, by="subsampleID")%>%
  select(subsampleID, 
         taxonID.y.x, 
         individualCount, 
         sampleID)%>%
  left_join(paraTaxMulti,by="subsampleID")%>%
  mutate (ind=individualCount/n)       #assume expertly IDed species composition is reflective of whole subsample (other fewer individuals expertly IDed than there are individuals in a subsample)

data_wide <- spread(paraTax_JERC_multi, taxonID.y.x, ind)%>%
  select(!individualCount)%>%
  select(!n)

#remove subsamples with multiple species from sort table
sort_JERC.df<-sort_JERC.df %>%
  filter(!subsampleID %in% paraTaxMulti$subsampleID)

#if there are no para or expert ID's, assume sort ID are correct
sort_JERC.df<-sort_JERC.df %>%
  dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%
  select(!taxonID)

data_wide1<-spread(sort_JERC.df, taxonID.y, individualCount)

#put two dataframes of species and counts in subsamples together

Example <- full_join(data_wide, data_wide1, by = "subsampleID")%>% #once again 2028 subsamples
  dplyr::mutate(sampleID.x = coalesce(sampleID.y, sampleID.x))%>% 
  replace(is.na(.), 0)%>%
  mutate(SELELL=SELELL.x+SELELL.y, 
         SELGRA=SELGRA.x+SELGRA.y)%>%
  select(!SELELL.x)%>%
  select(!SELELL.y)%>%
  select(!SELGRA.x)%>%
  select(!SELGRA.y)%>%
  select(!sampleID.y)%>%
  rename(sampleID=sampleID.x)

df<-Example %>% 
  select(!subsampleID)%>%
  group_by(sampleID) %>% 
  summarise(across(everything(), sum))

field_JERC.df <- BeetleJERC$bet_fielddata%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2014-12-31",
         !sampleCollected=="N") %>%                                            #removing missing samples 
  filter(!plotID=="JERC_022")%>%
  filter(!plotID=="JERC_023")%>%  
  select(plotID,
         trapID,
         nlcdClass,
         collectDate,
         eventID,
         sampleID,
         cupStatus,
         lidStatus)%>%
  transform(., plotIDeventID=paste(plotID, eventID))%>%
  distinct()

field_try.df<-field_JERC.df%>%
  left_join(df, by="sampleID")%>%
  mutate(cupStatus = replace_na("Ok"))%>%  #recode NA cup and lid status to Ok
  mutate(lidStatus = replace_na("Ok"))%>%
  filter (cupStatus!="Disturbed")%>%
  filter (cupStatus!="Missing")%>%
  filter(lidStatus!="Disturbed")%>%
  filter(lidStatus!="Missing")%>%          #remove cups that were disturbed or missing or had disturbed or missing lids; assume NA OK
  replace(is.na(.), 0)%>%                  #replace NAs with 0
  mutate(collectDate = ymd(collectDate),   #define collect date as a date
         collectDate = as_date(collectDate))%>%
  mutate (year=year(collectDate))%>%       #extract year
  transform(., plotIDyear=paste(plotID, year))%>% #create plot_year variable
  mutate(cups=1)%>%
  select(plotID, AGOCON:cups)

field_try2<-field_try.df%>%    #summarize # of cups sampled and number of individuals of each species by plot_year
  select(!plotID)%>%
  select(!year)%>%
  group_by(plotIDyear)%>%
  summarise(across(everything(), sum))%>%
  rowwise() %>%                #sum total individuals by plot_year
  mutate(totalIndividuals = sum(c_across(AGOCON:SELGRA)))%>%   
  mutate(density=totalIndividuals/cups)  #calculate yearly density for each plot

#calculate diversity for each year for each plot
species.df<-field_try2[,2:49]                                            # species only df for vegan diversity analysis - right amount of species 2388
field_try2$ShannonBeetle<-diversity(species.df, index = "shannon")       # shannon diversity index
field_try2$SimpsonBeetle<-diversity(species.df, index = "simpson")       # simpson diversity index
field_try2$spNum<-specnumber(species.df)

#should diversity indices be weighted somehow by the number of cups of individuals caught?

sum(field_try2$totalIndividuals) #410 individuals collected; 48 species
ave(field_try2$density) #0.1856806   ind/cup on average
std.error(field_try2$density) #SE: 0.03749157
ave(field_try2$ShannonBeetle) #0.8112845   
std.error(field_try2$ShannonBeetle) #SE: 0.09638695 
ave(field_try2$SimpsonBeetle) #0.5854742  
std.error(field_try2$SimpsonBeetle) #SE: 0.04615976

as.data.frame(colSums(field_try2[,c(2:49)]))           #find most abundant species

sum(field_try2$PASSUB2)/sum(field_try2$totalIndividuals)*100 #23.41%; 96 ind
sum(field_try2$ANIMER2)/sum(field_try2$totalIndividuals)*100 #10.49%; 43 ind
sum(field_try2$CYCSIG)/sum(field_try2$totalIndividuals)*100 #9.76%; 40 ind
sum(field_try2$CYCLAE)/sum(field_try2$totalIndividuals)*100 #7.56%; 31 ind
sum(field_try2$CYCOVU)/sum(field_try2$totalIndividuals)*100 #7.56%; 31 ind

#31 species with <5 ind

#join tree data to beetle data
DF2 <- field_try.df %>%
  select(plotIDyear, plotID)%>%
  distinct()
JERC_data.df<-field_try2%>%
  left_join(DF2, by="plotIDyear")%>%
  left_join(field_JERC.df %>% select(plotID, nlcdClass), by="plotID")%>%
  distinct()

#quick and dirty graphs and stats to check that same patterns hold as before
hist(JERC_data.df$density)
JERC_data.df$logDen<-log(JERC_data.df$density+1)
hist(JERC_data.df$logDen)

boxplot(JERC_data.df$logDen~JERC_data.df$nlcdClass)
summary(aov(JERC_data.df$logDen~JERC_data.df$nlcdClass)) #sig

hist(JERC_data.df$ShannonBeetle)
boxplot(JERC_data.df$ShannonBeetle~JERC_data.df$nlcdClass)
summary(aov(JERC_data.df$ShannonBeetle~JERC_data.df$nlcdClass)) #marg sig

hist(JERC_data.df$SimpsonBeetle)
boxplot(JERC_data.df$SimpsonBeetle~JERC_data.df$nlcdClass)
summary(aov(JERC_data.df$SimpsonBeetle~JERC_data.df$nlcdClass)) #not sig

#__________________________________________
# Read in relevant data BART----
#__________________________________________

# BART
expertTax_BART <- BeetleBART$bet_expertTaxonomistIDProcessed%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2014-12-31",                                           #2014 sampling began in mid-May
         family=="Carabidae")%>%                                               #remove samples that are not ground beetles
  mutate(taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%         #changing subspecies taxonID to species taxonID
  transform(scientificName=paste(genus, specificEpithet))%>%                   #this code will convert subspecies names to the species name 
  select(individualID,
         scientificName,
         taxonID)

paraTax_BART <- BeetleBART$bet_parataxonomistID%>%
  filter(collectDate < "2019-12-31",                                         #remove 2020 data: incomplete data 
         collectDate > "2014-12-31")%>%
  mutate(taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%         #changing subspecies taxonID to species taxonID
  select(individualID,
         subsampleID,
         taxonID)%>%
  left_join(., expertTax_JERC, by = "individualID")%>%                         # joining expertly identified ground beetles to the pinned ground beetles
  dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID.x))%>%                 # changing the non-expert ID's to expert ID'
  select(subsampleID,
         taxonID.y)%>%
  distinct()

#make a list of para and expert IDed subsamples that contain only 1 species
paraTaxUnique<-paraTax_BART%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n==1)

#filter subsamples containing more than 1 species out of paraTax table
paraTax_BART_unique<-paraTax_BART %>%
  filter(subsampleID %in% paraTaxUnique$subsampleID)

sort_BART.df <- BeetleBART$bet_sorting%>%
  filter(collectDate < "2019-12-31",   
         collectDate > "2014-12-31",                                           #remove 2020 data: incomplete data 
         sampleType == "carabid" | sampleType =="other carabid")%>%            #remove samples that are not ground beetles
  mutate(taxonID = replace(taxonID, taxonID == 'CYMPLA3', 'CYMPLA2'),
         taxonID = replace(taxonID, taxonID == 'SPHCAN1', 'SPHCAN'),
         taxonID = replace(taxonID, taxonID == 'SPHSTE1', 'SPHSTE3'))%>%         #changing subspecies taxonID to species taxonID
  select(subsampleID,
         taxonID,
         individualCount,
         sampleID)%>%
  replace(is.na(.), 1)                                                        #one subsample has no individualCount. assume is 1
sum(sort_BART.df$individualCount)                                              #9373 individuals in total 2015-2019

#join subsamples with only 1 species to sort table
sort_BART.df<-sort_BART.df%>%
  left_join(paraTax_BART_unique, by="subsampleID")

#all para and expert IDed subsamples contain only 1 species
paraTaxMulti<-paraTax_BART%>%
  select(subsampleID)%>%
  count(subsampleID)%>%
  filter(n>1)

#if there are no para or expert ID's, assume sort ID are correct
sort_BART.df<-sort_BART.df %>%
  dplyr::mutate(taxonID.y = coalesce(taxonID.y, taxonID))%>%
  select(!taxonID)

data_wide1<-spread(sort_BART.df, taxonID.y, individualCount)

#put two dataframes of species and counts in subsamples together

Example <- data_wide1%>% 
  replace(is.na(.), 0)

df<-Example %>% 
  select(!subsampleID)%>%
  group_by(sampleID) %>% 
  summarise(across(everything(), sum))

field_BART.df <- BeetleBART$bet_fielddata%>%
  filter(collectDate < "2019-12-31",                                           #remove 2020 data: incomplete data 
         collectDate > "2014-12-31",
         !sampleCollected=="N") %>%                                            #removing missing samples 
  select(plotID,
         trapID,
         nlcdClass,
         collectDate,
         eventID,
         sampleID,
         cupStatus,
         lidStatus)%>%
  transform(., plotIDeventID=paste(plotID, eventID))%>%
  distinct()

field_try.df<-field_BART.df%>%
  left_join(df, by="sampleID")%>%
  mutate(cupStatus = replace_na("Ok"))%>%  #recode NA cup and lid status to Ok
  mutate(lidStatus = replace_na("Ok"))%>%
  filter (cupStatus!="Disturbed")%>%
  filter (cupStatus!="Missing")%>%
  filter(lidStatus!="Disturbed")%>%
  filter(lidStatus!="Missing")%>%          #remove cups that were disturbed or missing or had disturbed or missing lids; assume NA OK
  replace(is.na(.), 0)%>%                  #replace NAs with 0
  mutate(collectDate = ymd(collectDate),   #define collect date as a date
         collectDate = as_date(collectDate))%>%
  mutate (year=year(collectDate))%>%       #extract year
  transform(., plotIDyear=paste(plotID, year))%>% #create plot_year variable
  mutate(cups=1)%>%
  select(plotID, AGORET:cups)

field_try2<-field_try.df%>%    #summarize # of cups sampled and number of individuals of each species by plot_year
  select(!plotID)%>%
  select(!year)%>%
  group_by(plotIDyear)%>%
  summarise(across(everything(), sum))%>%
  rowwise() %>%                #sum total individuals by plot_year
  mutate(totalIndividuals = sum(c_across(AGORET:SYNIMP)))%>%   
  mutate(density=totalIndividuals/cups)  #calculate yearly density for each plot

#calculate diversity for each year for each plot
species.df<-field_try2[,2:29]                                            # species only df for vegan diversity analysis - right amount of species 2388
field_try2$ShannonBeetle<-diversity(species.df, index = "shannon")       # shannon diversity index
field_try2$SimpsonBeetle<-diversity(species.df, index = "simpson")       # simpson diversity index
field_try2$spNum<-specnumber(species.df)

#should diversity indices be weighted somehow by the number of cups of individuals caught?

sum(field_try2$totalIndividuals) #9373 individuals collected; 28 species
ave(field_try2$density) #5.838311    ind/cup on average
std.error(field_try2$density) #SE: 0.5707725
ave(field_try2$ShannonBeetle) #1.722614    
std.error(field_try2$ShannonBeetle) #SE: 0.03798594 
ave(field_try2$SimpsonBeetle) #0.744534   
std.error(field_try2$SimpsonBeetle) #SE: 0.01186412

as.data.frame(colSums(field_try2[,c(2:29)]))           #find most abundant species

sum(field_try2$SYNIMP)/sum(field_try2$totalIndividuals)*100 #40.99%; 3842 ind
sum(field_try2$PTEPEN)/sum(field_try2$totalIndividuals)*100 #11.50%; 1072 ind
sum(field_try2$PTETRI3)/sum(field_try2$totalIndividuals)*100 #11.21%; 1051 ind
sum(field_try2$CYMNEG)/sum(field_try2$totalIndividuals)*100 #7.07%; 663 ind
sum(field_try2$PTELAC2)/sum(field_try2$totalIndividuals)*100 #6.89%; 646 ind

#7 species with <5 ind

#join tree data to beetle data
DF2 <- field_try.df %>%
  select(plotIDyear, plotID)%>%
  distinct()
BART_data.df<-field_try2%>%
  left_join(DF2, by="plotIDyear")%>%
  left_join(field_BART.df %>% select(plotID, nlcdClass), by="plotID")%>%
  distinct()

#quick and dirty graphs and stats to check that same patterns hold as before
hist(BART_data.df$density)
BART_data.df$logDen<-log(BART_data.df$density+1)
hist(BART_data.df$logDen)

boxplot(BART_data.df$logDen~BART_data.df$nlcdClass)
summary(aov(BART_data.df$logDen~BART_data.df$nlcdClass)) #not sig

hist(BART_data.df$ShannonBeetle)
boxplot(BART_data.df$ShannonBeetle~BART_data.df$nlcdClass)
summary(aov(BART_data.df$ShannonBeetle~BART_data.df$nlcdClass)) #not sig

hist(BART_data.df$SimpsonBeetle)
boxplot(BART_data.df$SimpsonBeetle~BART_data.df$nlcdClass)
summary(aov(BART_data.df$SimpsonBeetle~BART_data.df$nlcdClass)) #not sig



#__________________________________________
#Make DFs needed for stats and graphs
#__________________________________________
#add site ID to each "_data.df" for HARV, TALL, JERC, and BART
HARV_data.df$siteID <- "HARV"
TALL_data.df$siteID <- "TALL"
JERC_data.df$siteID <- "JERC"
BART_data.df$siteID <- "BART"
rm(vegB.df)
vegB.df<-bind_rows(HARV_data.df, TALL_data.df, JERC_data.df, BART_data.df)

#df for stats of all sites and nlcd class
vegB1.df <- vegB.df%>%
  select(ShannonBeetle,SimpsonBeetle, density, nlcdClass, siteID, plotID) #only have columns needed

#df for stats and figures of beetle + veg at Harv and Tall
vegB2.df <- vegB.df%>%
  filter(., !siteID=="JERC",
         !siteID=="BART")%>%  
  mutate(PerAM_BA=AM_BA/totalBA*100,
         PerAM_ST=AM_ST/totalStems*100)

#__________________________________________
# stats
#__________________________________________
  
#__________________________________________
# 3 variables ~ nlcdClass (all sites: HARV, TALL, BART, JERC)
#__________________________________________
#removing no observations 
subset_df <- vegB1.df
subset_df <- subset(subset_df, density != 0) #remove no observations

m1 <- nlme::lme(ShannonBeetle ~ nlcdClass, random = ~1| siteID/plotID, 
               data = vegB1.df);summary(m1);anova(m1); shapiro.test(resid(m1));
m2 <- nlme::lme(SimpsonBeetle ~ nlcdClass,random = ~1| siteID/plotID, 
               data = vegB1.df);summary(m2);anova(m2);shapiro.test(resid(m2));
m3 <- nlme::lme(log(density) ~ nlcdClass, random = ~1| siteID/plotID, 
               data = subset_df);summary(m3); shapiro.test(resid(m3));anova(m3);lsmeans(m3, pairwise~nlcdClass, adjust="tukey")   

#__________________________________________
# 3 variables ~ Tree div, %EG, %ECM (TALL & HARV)
#__________________________________________
#Basal area models
m4 <- nlme::lme(ShannonBeetle ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB2.df);summary(m4); shapiro.test(resid(m4))
m5 <- nlme::lme(SimpsonBeetle ~ Sim_BA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB2.df);summary(m5); shapiro.test(resid(m5));anova(m5)
m6 <- nlme::lme(log(density) ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB2.df);summary(m6); shapiro.test(resid(m6));anova(m6)

#removing no observations 
subset2_df <- vegB2.df
subset2_df <- subset(subset2_df, density != 0) #remove no observations
#total stems models
m7 <- nlme::lme(ShannonBeetle ~ Shan_ST + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB2.df);summary(m7); shapiro.test(resid(m7));anova(m7)
m8 <- nlme::lme(SimpsonBeetle ~ Sim_ST + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB2.df);summary(m8); shapiro.test(resid(m8));anova(m8)
m9 <- nlme::lme(log(density) ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = subset2_df);summary(m9); shapiro.test(resid(m9));anova(m9)


#__________________________________________
#Top 5 
#__________________________________________
#HARV - density models
  #BA models
  
  m10 <- nlme::lme(logcCARGOR ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m10); shapiro.test(resid(m10));anova(m10)
  m11 <- nlme::lme(logSYN ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m11); shapiro.test(resid(m11));anova(m11)
  m12 <- nlme::lme(logSPHSTE3_den ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m12); shapiro.test(resid(m12));anova(m12)
  m13 <- nlme::lme(logPTETRI3 ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m13); shapiro.test(resid(m13));anova(m13)
  m14 <- nlme::lme(logPTEPEN_den ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m14); shapiro.test(resid(m14));anova(m14)

#TALL - density models
  #BA models
  m15 <- nlme::lme(logCYCCON2 ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m15); shapiro.test(resid(m15));anova(m15)
  m16 <- nlme::lme(logDICDIL5 ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID,  #PerEG_BA p value is 0.04428 *
                 data = TALL_data.df);summary(m16); shapiro.test(resid(m16));anova(m16)
  m17 <- nlme::lme(logCYCFRE ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m17); shapiro.test(resid(m17));anova(m17)
  m18 <- nlme::lme(logANIHAP ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m18); shapiro.test(resid(m18));anova(m18)
  m19 <- nlme::lme(logPASDEP ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m19); shapiro.test(resid(m19));anova(m19)
 

#HARV - relative abundance models
  
  #BA models
  m20 <- nlme::lme(logCARGOR_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m20); shapiro.test(resid(m20));anova(m20)
  m21 <- nlme::lme(SYNIMP_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, #PerEG_BA p value on the edge 0.0573 .
                 data = HARV_data.df);summary(m21); shapiro.test(resid(m21));anova(m21)
  m22 <- nlme::lme(logSPHSTE3_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m22); shapiro.test(resid(m22));anova(m22)
  m23 <- nlme::lme(logPTETRI3_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m23); shapiro.test(resid(m23));anova(m23)
  m24 <- nlme::lme(logPTEPEN_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = HARV_data.df);summary(m24); shapiro.test(resid(m24));anova(m24)
  

#TALL - relative abundance
  #BA models
  m25 <- nlme::lme(CYCCON2_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, #Total BA NOT EG!! predicts abundnace (p= 0.04359 *)
                 data = TALL_data.df);summary(m25); shapiro.test(resid(m25));anova(m25)
  m26 <- nlme::lme(logDICDIL5_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m26); shapiro.test(resid(m26));anova(m26)
  m27 <- nlme::lme(logCYCFRE_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m27); shapiro.test(resid(m27));anova(m27)
  m28 <- nlme::lme(logANIHAP_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m28); shapiro.test(resid(m28));anova(m28)
  m29 <- nlme::lme(PASDEP_relabun ~ totalBA + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                 data = TALL_data.df);summary(m29); shapiro.test(resid(m29));anova(m29)
#_____________________________________
#TOP FIVE Stem models
  #HARV - density models
  #Stem models
  
  m30 <- nlme::lme(logcCARGOR ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary( m30); shapiro.test(resid( m30));anova(m30)
  m31 <- nlme::lme(logSYN ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m31); shapiro.test(resid(m31));anova(m31)
  m32 <- nlme::lme(logSPHSTE3_den ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m32); shapiro.test(resid(m32));anova(m32)
  m33 <- nlme::lme(logPTETRI3 ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m33); shapiro.test(resid(m33));anova(m33)
  m34 <- nlme::lme(logPTEPEN_den ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m34); shapiro.test(resid(m34));anova(m34)
  
  #TALL - density models
  #Stem models
  m35 <- nlme::lme(logCYCCON2 ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m35); shapiro.test(resid(m35));anova(m35)
  m36 <- nlme::lme(logDICDIL5 ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID,  #PerEG_BA p value is 0.04428 *
                   data = TALL_data.df);summary(m36); shapiro.test(resid(m36));anova(m36)
  m37 <- nlme::lme(logCYCFRE ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m37); shapiro.test(resid(m37));anova(m37)
  m38 <- nlme::lme(logANIHAP ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m38); shapiro.test(resid(m38));anova(m38)
  m39 <- nlme::lme(logPASDEP ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m39); shapiro.test(resid(m39));anova(m39)
  
  
  #HARV - relative abundance models
  
  #Stem models
  m40 <- nlme::lme(logCARGOR_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m40); shapiro.test(resid(m40));anova(m40)
  m41 <- nlme::lme(SYNIMP_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, #PerEG_BA p value on the edge 0.0573 .
                   data = HARV_data.df);summary(m41); shapiro.test(resid(m41));anova(m41)
  m42 <- nlme::lme(logSPHSTE3_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m42); shapiro.test(resid(m42));anova(m42)
  m43 <- nlme::lme(logPTETRI3_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m43); shapiro.test(resid(m43));anova(m43)
  m44 <- nlme::lme(logPTEPEN_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = HARV_data.df);summary(m44); shapiro.test(resid(m44));anova(m44)
  
  
  #TALL - relative abundance
  #Stem models
  m45 <- nlme::lme(CYCCON2_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, #Total BA NOT EG!! predicts abundnace (p= 0.04359 *)
                   data = TALL_data.df);summary(m45); shapiro.test(resid(m45));anova(m45)
  m46 <- nlme::lme(logDICDIL5_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m46); shapiro.test(resid(m46));anova(m46)
  m47 <- nlme::lme(logCYCFRE_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m47); shapiro.test(resid(m47));anova(m47)
  m48 <- nlme::lme(logANIHAP_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m48); shapiro.test(resid(m48));anova(m48)
  m49 <- nlme::lme(PASDEP_relabun ~ totalStems + PerEG_BA+PerECM_BA, random = ~1|plotID, 
                   data = TALL_data.df);summary(m49); shapiro.test(resid(m49));anova(m49)  


#__________________________________________
# Figure 1: Three panel of Beetle Shan/Sim/ Density ~ forest cover (nlcdClass) ----
#__________________________________________


p1 <- ggplot(vegB1.df, aes(x=nlcdClass, y=ShannonBeetle, fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = vegB1.df, aes(shape = siteID), size = 4, alpha = 0.5, width = 0.25, show.legend = T)+
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
  theme(axis.title.x=element_text(size=24), 
        axis.title.y=element_text(size=24), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))+
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
  ylim(0, 2.5)+
    scale_y_continuous(breaks=c(0,.25,.5,.75,1,1.25,1.5,1.75,2,2.25))
p1


p2 <- ggplot(vegB1.df, aes(x=nlcdClass, y=SimpsonBeetle, fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = vegB1.df, aes(shape = siteID), size = 4, alpha = 0.5, width = 0.25, show.legend = FALSE)+
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
  theme(axis.title.x=element_text(size=24), 
        axis.title.y=element_text(size=24), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  scale_y_continuous(breaks=c(0,.25,.5,.75,1))
p2


p3 <- ggplot(vegB1.df, aes(x=nlcdClass, y=log(density), fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = vegB1.df, aes(shape = siteID), size = 4, alpha = 0.5, width = 0.25, show.legend = F)+
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
  theme(axis.title.x=element_text(size=24), 
        axis.title.y=element_text(size=24), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))+
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
  annotate(geom="text", x=1, y=3.5, label="a",
           color="black",alpha = 0.70, size=8)+
  annotate(geom="text", x=2, y=3.5, label="ab",
           color="black",alpha = 0.70, size=8)+
  annotate(geom="text", x=3, y=3.5, label="b",
           color="black",alpha = 0.70, size=8)+
  scale_y_continuous(labels = scales::label_number(accuracy = 0.1), breaks=c(-4.0,-2.0, 0.0, 2.0, 4.0), limits=c(-4.0,4.0))
p3

#cow plot to stack them 
library(cowplot)
fig1<-plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 20, ncol = 1)
fig1
pdf("/Users/JaneyLienau/Desktop/Model_Graphs/fig1.pdf", width = 7, height = 21)
plot(fig1)
dev.off()

#__________________________________________
# Figure 2: Three panel of Beetle diversity/density ~ %EG----
#__________________________________________

p1<-ggplot(vegB2.df,aes(x=PerEG_BA, y=ShannonBeetle))+
 # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(shape = siteID),alpha = 0.5, size = 4)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Ground Beetle Diversity\n(Shannon Index)')+
  scale_x_discrete(name =NULL,labels = NULL)+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=24), 
        axis.title.y=element_text(size=24), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))+
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
  guides(color=FALSE)+
  scale_y_continuous(breaks=c(0.0,0.5,1,1.5,2,2.25)) #0.0,.25,0.5,1,1.25,1.5,2
p1


p2<-ggplot(vegB2.df,aes(x=PerEG_BA, y=SimpsonBeetle))+
  #geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(shape = siteID),alpha = 0.5, size = 4, show.legend = F)+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Ground Beetle Diversity\n(Simpson Index)')+
  scale_x_discrete(name =NULL,labels = NULL)+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=24), 
        axis.title.y=element_text(size=24), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))+
  theme(axis.ticks.length=unit(-0.25, "cm"))+
  theme(axis.ticks = element_line(colour = "black", size = 0.4))+
  theme(axis.ticks.x = element_blank())+
  guides(color=FALSE)
p2


p3<-ggplot(vegB2.df,aes(x=PerEG_BA, y=log(density)))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(shape = siteID),alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Evergreen Trees\nRelative Abundance (%)', y = 'Ground Beetle Density (log)\n(Individuals/Cup)', color='Site')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=24), 
        axis.title.y=element_text(size=24), 
        axis.text.x=element_text(size=20), 
        axis.text.y=element_text(size=20))+
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
  scale_y_continuous(breaks=c(-2.5,-1.0,0,1.0,2.5))
p3


fig2<-plot_grid(p1, p2, p3, labels = c('A', 'B', 'C'), label_size = 20, ncol = 1)
fig2
pdf("/Users/JaneyLienau/Desktop/Model_Graphs/fig2.pdf", width = 7, height = 21)
plot(fig2)
dev.off()

#__________________________________________
#  Supplemental Figures 2: density of top 10 species ~ evergreen tree abundance)
#__________________________________________

p1<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=CARGOR_den))+
  geom_smooth(method = 'lm', formula = 'y ~ x', color="gray") +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'CARGOR\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)+
  ggtitle("HARV")+
  theme(title=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))
p1

p2<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=SYNIMP_den))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'SYNIMP\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p2


p3<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=SPHSTE3_den))+
  #geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'SPHSTE3\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p3

p4<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=PTETRI3_den))+
  #  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'PTETRI3\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p4

p5<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=PTEPEN_den))+
  # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Evergreen Trees\nRelative Abundance (%)', y = 'PTEPEN\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p5
#tall
names(TALL_data.df)
p6<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=CYCCON2_den))+
  #geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'CYCCON2\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)+
  ggtitle("TALL")+
  theme(title=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))
p6

p7<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=DICDIL5_den))+
  geom_smooth(method = 'lm', formula = 'y ~ x', color="gray") +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'DICDIL5\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p7

p8<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=CYCFRE_den))+
  # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'CYCFRE\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p8

p9<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=ANIHAP_den))+
  # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'ANIHAP\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p9

p10<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=PASDEP_den))+
  #geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Evergreen Trees\nRelative Abundance (%)', y = 'PASDEP\n(Individuals/Cup)', color='plotID')+
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
  guides(color=none)
p10
#cow plot to stack them 

#most abundant species plot----relative abundance currently in graphs, change as needed

Harv.plot<-plot_grid(p1,p2,p3,p4,p5, ncol = 1, labels = c('a', 'b', 'c', 'd','e'), label_size = 16)
Harv.plot
Tall.plot<-plot_grid(p6,p7,p8,p9,p10, ncol = 1, labels = c('f', 'g', 'h', 'i','j'), label_size = 16)
speciesgrid<-plot_grid(Harv.plot, Tall.plot, ncol = 2)
speciesgrid
pdf("/Users/JaneyLienau/Desktop/Model_Graphs/Supp2_species_grid_density.pdf", width = 10, height = 21)
plot(speciesgrid)
dev.off()

#__________________________________________
#  Supplemental Figures 3 -- relative abundance of top 10 species ~ evergreen tree abundance
#__________________________________________

names(HARV_data.df)
#plot 1
p1<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=CARGOR_relabun))+
  geom_smooth(method = 'lm', formula = 'y ~ x', color ="gray") +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'CARGOR\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)+
  ggtitle("HARV")+
  theme(title=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))
p1

p2<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=SYNIMP_relabun))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'SYNIMP\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)
p2


p3<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=SPHSTE3_relabun))+
  #geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'SPHSTE3\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)
p3

p4<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=PTETRI3_relabun))+
#  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'PTETRI3\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)
p4

p5<-ggplot(HARV_data.df,aes(x=PerEG_BA, y=PTEPEN_relabun))+
 # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Evergreen Trees\nRelative Abundance (%)', y = 'PTEPEN\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)
p5
#tall
names(TALL_data.df)
p6<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=CYCCON2_relabun))+
 # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'CYCCON2\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)+
  ggtitle("TALL")+
  theme(title=element_text(size=20))+
  theme(plot.title = element_text(hjust = 0.5))
p6

p7<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=DICDIL5_relabun))+
  #geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'DICDIL5\nRelative Abundance (%)', color='plotID')+
  theme(axis.title.x = element_text(margin = margin(t = 5, b=5)), 
        axis.title.y = element_text(margin = margin(l = 5, r=5)), 
        axis.text.x=element_text(margin = margin(t=10)), 
        axis.text.y=element_text(margin = margin(r = 10)))+
  theme(axis.title.x=element_text(size=18), 
        axis.title.y=element_text(size=18), 
        axis.text.x=element_text(size=16), 
        axis.text.y=element_text(size=16))+
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
  guides(color=none)
p7

p8<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=CYCFRE_relabun))+
 # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'CYCFRE\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)
p8

p9<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=ANIHAP_relabun))+
 # geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'ANIHAP\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)
p9

p10<-ggplot(TALL_data.df,aes(x=PerEG_BA, y=PASDEP_relabun))+
  #geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(alpha = 0.5, size = 4, show.legend = F)+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Evergreen Trees\nRelative Abundance (%)', y = 'PASDEP\nRelative Abundance (%)', color='plotID')+
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
  guides(color=none)
p10
#cow plot to stack them 

#most abundant species plot----relative abundance currently in graphs, change as needed

Harv.plot<-plot_grid(p1,p2,p3,p4,p5, ncol = 1, labels = c('a', 'b', 'c', 'd','e'), label_size = 16)
Harv.plot
Tall.plot<-plot_grid(p6,p7,p8,p9,p10, ncol = 1, labels = c('f', 'g', 'h', 'i','j'), label_size = 16)
speciesgrid<-plot_grid(Harv.plot, Tall.plot, ncol = 2)
speciesgrid
pdf("/Users/JaneyLienau/Desktop/Model_Graphs/Supp3_species_grid_relAbun.pdf", width = 10, height = 21)
plot(speciesgrid)
dev.off()

#saving key data frames for functional analysis and ordination
write_csv(HARV_data.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/HARV_data.csv")
write_csv(TALL_data.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/TALL_data.csv")
write_csv(JERC_data.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/JERC_data.csv")
write_csv(BART_data.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/BART_data.csv")

write_csv(allHARV.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/allHARV.csv")
write_csv(allTALL.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/allTALL.csv")

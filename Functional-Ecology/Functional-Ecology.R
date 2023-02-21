
# -----------------------------------------------------------
#  functional ecology
# -----------------------------------------------------------

library(dplyr)
library(nlme)

write.csv(diversityPlotYr_TALL.df, "/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/SpeciesPlotYr_TALL.csv")

# -----------------------------------------------------------
#  read in data
# -----------------------------------------------------------
HARV.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/Functional-Ecology/diversityPlotYr_HARV.csv")

TALL.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/Functional-Ecology/SpeciesPlotYr_TALL.csv")

vegB.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/vegB.csv")

SpeciesCodes.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/Functional-Ecology/SpeciesCodes.df.csv")

function.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/Functional-Ecology/gb-functional-groups.csv")


# -----------------------------------------------------------
#  Functional Data
# -----------------------------------------------------------

#test.df <- left_join(SpeciesCodes.df, function.df, by = "scientificName")


test <- colSums(HARV.df[,5:37])
test

test1 <- colSums(TALL.df[,5:50])
test1



# -----------------------------------------------------------
#  lme models
# -----------------------------------------------------------

m <- nlme::lme("functional group" ~ nlcdClass, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m));Anova(m)


# -----------------------------------------------------------
# 3 variables ~ nlcdClass
# -----------------------------------------------------------
#make density of each functional group predicted by nlcdClass
m <- nlme::lme(density of feeding mode and other groups~ nlcdClass, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m));Anova(m)


# -----------------------------------------------------------
# 3 variables ~ Tdiv, %EG, %ECM
# -----------------------------------------------------------

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


# -----------------------------------------------------------
# 3 variables ~ Tden, %EG, %ECM
# -----------------------------------------------------------

m <- nlme::lme(ShannonBeetle ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(SimpsonBeetle ~ totalBA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

# -----------------------------------------------------------
#diversity begets diversity, density begets density
# -----------------------------------------------------------

m <- nlme::lme(ShannonBeetle ~Shan_BA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(SimpsonBeetle ~ Sim_BA,random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))
m <- nlme::lme(density ~ totalBA, random = ~1| siteID/plotID, 
               data = vegB.df);summary(m); shapiro.test(resid(m))

# -----------------------------------------------------------
#diversity predicted by tree density
# -----------------------------------------------------------
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
# -----------------------------------------------------------

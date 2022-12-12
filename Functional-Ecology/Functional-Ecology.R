
# -----------------------------------------------------------
#  functional ecology
# -----------------------------------------------------------

library(dplyr)
library(nlme)
library(tidyverse)
library(car)   
library(RColorBrewer)
library(cowplot)

# -----------------------------------------------------------
#  read in data
# -----------------------------------------------------------
HARV.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/Functional-Ecology/diversityPlotYr_HARV.csv")

TALL.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/Functional-Ecology/SpeciesPlotYr_TALL.csv")

vegB.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/vegB.csv")

function.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/Evergreen-abundance-drives-ground-beetle-diversity-and-density-in-eastern-temperate-forests/Functional-Ecology/function-table-dec.csv")
# -----------------------------------------------------------
#  Data cleaning
# -----------------------------------------------------------
function.df <- function.df%>% #cleaning data frame a bit
  select(-scientificName)%>%
  select(-DP.filght)%>%
  na.omit()


#removed species observations that were less than 5 observations total at a site
TALL.df <- TALL.df%>%
  filter(!plotIDyear == "TALL_004 2016") #removing a year that didn't have any observations


HARV.df <- pivot_longer(HARV.df, 5:37, values_to = "count")%>% #reformatting data
 rename(., taxonID = "name")%>%
  select(-X)
sum(HARV.df$count) #2389

TALL.df <- pivot_longer(TALL.df, 5:52, values_to = "count")%>% #reformatting data
  rename(., taxonID = "name")%>%
  select(-X)
sum(TALL.df$count) #755

harvtall.df <- rbind(TALL.df, HARV.df)#joining the data sets from both sites

#remove species that were observed less than 5 times
test <- harvtall.df%>%
  select(taxonID, count)%>%
  group_by(taxonID)%>%
  summarise(count = sum(count)) #there are 12 species
sum(test$count)#total 3144
tostay <- test%>%
  filter(., count > 4)
sum(tostay$count)#3085
toremove <- test%>%
  filter(., count < 4)
sum(toremove$count)#47 observations total
toremove$taxonID

harvtall.df <- harvtall.df%>% 
  filter(!taxonID == "ACUHYD",
         !taxonID == "AMBMEX",
         !taxonID == "AMPINT",
         !taxonID == "ANIRUS",
         !taxonID == "APESIN",
         !taxonID == "CALFRI",
         !taxonID == "CICSEX",
         !taxonID == "CRADUB",
         !taxonID == "DICELO",
         !taxonID == "GASHON",
         !taxonID == "HARKAT",
         !taxonID == "HARPEN",
         !taxonID == "HARRUB",
         !taxonID == "HELNIG",
         !taxonID == "LEBPUL",
         !taxonID == "NOTNOV",
         !taxonID == "NOTSAY",
         !taxonID == "NOTTER",
         !taxonID == "OLIPAR",
         !taxonID == "POELUC",
         !taxonID == "POLLAE",
         !taxonID == "PTESP23",
         !taxonID == "SELFOS",
         !taxonID == "SELGRA",
         !taxonID == "SELPAL",
         !taxonID == "SPHSTE1",
         !taxonID == "STEPLE")
sum(harvtall.df$count)#3097
#left join functional data into main df
harvtall.df <- left_join(harvtall.df, 
                         function.df, 
                         by = "taxonID")


harvtall.df <- na.omit(harvtall.df) #remove a few obervations that didn't have functional data (these species also were >10 observations(low abundance) at each site)
sum(harvtall.df$count)#2857 - lost 240 that we don't have food data for

#joing with veg data
harvtall.df <-  pivot_wider(harvtall.df, 
                            names_from = "Food", 
                            values_from = "count")%>%
  replace(is.na(.),0)

sum(harvtall.df[,5:6])##2857
#summarize table by plot year for each feeding mode
#calculate percent omnivore and predator
harvtall.df <- harvtall.df%>%
  group_by(plotIDyear)%>%                                                        #summarize at the plot level
  summarise(Omnivore = sum(Omnivore),
            Predator = sum(Predator))%>%
  mutate(., PerOmnivore = Omnivore/(Omnivore+Predator)*100,
         PerPredator = Predator/(Omnivore+Predator)*100)


harvtall.df <- left_join(harvtall.df, vegB.df, by = "plotIDyear")
#calculate density count/cups

harvtall.df <- harvtall.df%>%
  mutate(., OmnDensity = Omnivore/cupNum, 
         PreDensity = Predator/cupNum)

# -----------------------------------------------------------
#  lme models - Feeding group by nlcdClass
# -----------------------------------------------------------
#
m <- nlme::lme(OmnDensity ~ nlcdClass, random = ~1| siteID/plotID, 
               data = harvtall.df);summary(m); shapiro.test(resid(m));Anova(m)
m <- nlme::lme(PreDensity ~ nlcdClass, random = ~1| siteID/plotID, 
               data = harvtall.df);summary(m); shapiro.test(resid(m));Anova(m)
#-2.263075  0.0535 less in eg forest


#omnivore count plot
p1 <- ggplot(harvtall.df, aes(x=nlcdClass, y=OmnDensity, fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = harvtall.df, aes(shape = siteID), size = 2, alpha = 0.5, width = 0.25, show.legend = T)+
  scale_shape_manual(values=c(15, 16,18,17))+
  scale_fill_brewer(palette="Dark2", guide=FALSE)+
  labs(x = NULL, y = 'Omivorous Ground Beetle Density')+
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
    legend.direction="horizontal")
p1

#predator count plot
p2 <- ggplot(harvtall.df, aes(x=nlcdClass, y=PreDensity, fill=nlcdClass)) + 
  geom_boxplot(alpha=0.42)+
  geom_jitter(data = harvtall.df, aes(shape = siteID), size = 2, alpha = 0.5, width = 0.25, show.legend = F)+
  scale_shape_manual(values=c(15, 16,18,17))+
  scale_fill_brewer(palette="Dark2", guide=FALSE)+
  labs(x = NULL, y = 'Predatory Ground Beetle Density')+
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
    legend.position="top",
    legend.direction="horizontal")
p2

#cow plot to stack omnivore/predator count plots
p1p2<-plot_grid(p1, p2, labels = c('A', 'B'), label_size = 20, ncol = 1)
p1p2
pdf("/Users/JaneyLienau/Desktop/p1p2.pdf", width = 7, height = 14)
plot(p1p2)
dev.off()


# -----------------------------------------------------------
# 3 variables ~ Tdiv, %EG, %ECM EG sig predictor of all version predator + omnivore except "Predator"
# -----------------------------------------------------------

m <- nlme::lme(OmnDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

p3<-ggplot(harvtall.df,aes(x=PerEG_BA, y=OmnDensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Omnivore Species Density', color='plotID')+
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

p4<-ggplot(harvtall.df,aes(x=PerEG_BA, y=PreDensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point()+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance\n of Evergreen Trees (%)', y = 'Predator Species Density', color='plotID')+
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

#cow plot to stack omnivore/predator percent plots
p3p4<-plot_grid(p3, p4, labels = c('A', 'B'), label_size = 20, ncol = 1)
p3p4
pdf("/Users/JaneyLienau/Desktop/p3p4.pdf", width = 7, height = 14)
plot(p3p4)
dev.off()

# -----------------------------------------------------------
# 3 variables ~ Tden, %EG, %ECM - EG sig predictor of predator + omnivore
# -----------------------------------------------------------

m <- nlme::lme(OmnDensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~ totalBA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

#samre results as total basal area
m <- nlme::lme(OmnDensity ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~ totalStems + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))


# -----------------------------------------------------------
#diversity predicted by tree density - NS
# -----------------------------------------------------------
m <- nlme::lme(OmnDensity ~ totalBA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~totalBA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(OmnDensity ~ totalStems,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~totalStems, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))



# -----------------------------------------------------------

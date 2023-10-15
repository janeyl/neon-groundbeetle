#________________________________________________________________________
#Packages---------------------------
#________________________________________________________________________
library(dplyr)
library(nlme)
library(tidyverse)
library(car)   
library(RColorBrewer)
library(cowplot)

#________________________________________________________________________
#Raw data-------------------
#________________________________________________________________________
HARV.df <- read.csv("HARV_data.csv")

TALL.df <- read.csv("TALL_data.csv")

BART.df <- read.csv("BART_data.csv")

JERC.df <- read.csv("JERC_data.csv")

vegB.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/vegB.csv")

function.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/Functional-Ecology/function-table-dec.csv")

#________________________________________________________________________
#Data Cleaning -----------------------------
#________________________________________________________________________

funspecies.df <- function.df%>%
  select(-DP.filght,
         -scientificName)

function.df <- function.df%>% #cleaning data frame a bit
  select(-scientificName)%>%
  select(-DP.filght)%>%
  na.omit()


sum(HARV.df[,2:34])#3823
HARV.df <- pivot_longer(HARV.df, 2:34, values_to = "count")%>% #reformatting data
 rename(., taxonID = "name")%>%
  select(plotIDyear,
         plotID,
         count,
         taxonID,
         siteID)
sum(HARV.df$count) #3823

sum(TALL.df[,2:34])#470
TALL.df <- pivot_longer(TALL.df, 2:34, values_to = "count")%>% #reformatting data
  rename(., taxonID = "name")%>%
  select(plotIDyear,
         plotID,
         count,
         taxonID,
         siteID)
sum(TALL.df$count) #470


harvtall.df <- rbind(TALL.df, HARV.df)#joining the data sets from both sites

harvtall.df <- harvtall.df%>%
  select(-siteID)

#remove species that were observed less than 5 times
test <- harvtall.df%>%
  select(taxonID, count)%>%
  group_by(taxonID)%>%
  summarise(count = sum(count)) #there are 12 species
sum(test$count)#total 4293
tostay <- test%>%
  filter(., count > 5)
sum(tostay$count)#4240
toremove <- test%>%
  filter(., count < 5)
sum(toremove$count)#43 observations total
toremove$taxonID

harvtall.df <- harvtall.df%>% 
  filter(!taxonID == "ACUHYD",
         !taxonID == "AGOCON",
         !taxonID == "AMPINT",
         !taxonID == "ANIRUS",
         !taxonID == "APESIN",
         !taxonID == "CALFRI",
         !taxonID == "CICSEX",
   
         !taxonID == "DICELO",
         !taxonID == "GASHON",

         !taxonID == "HARPEN",
         !taxonID == "HARRUB",
         !taxonID == "HELNIG",
         !taxonID == "LOXCRE",
         !taxonID == "NOTNOV",

         !taxonID == "OLIPAR",
         !taxonID == "POELUC",
         !taxonID == "POLLAE",
         !taxonID == "PTESP23",
         !taxonID == "SCAUNI1",
         !taxonID == "SELELL",
         !taxonID == "SELGRA",

         !taxonID == "SPHSTE1")
sum(harvtall.df$count)#4250
#left join functional data into main df
harvtall.df <- left_join(harvtall.df, 
                         function.df, 
                         by = "taxonID")

#species that we need info:
infospp <- harvtall.df%>%
  select(count,
         taxonID)%>%
  group_by(taxonID)%>%
  summarise(across(everything(), sum))%>%
  rowwise()%>%
  left_join(., function.df)
sum(infospp$count)#4250 #6 species we don't have feeding data for 
#90.16667 don't have info 90.16667/4250*100 =2.121569%

#TRIAUT 31.16667 indi no info
#ANIHAP 27.50000 indi no info :(

harvtall.df <- na.omit(harvtall.df) #remove a few obervations that didn't have functional data (these species also were <10 observations(low abundance) at each site)
sum(harvtall.df$count)#2857 - lost 240 that we don't have food data for

#joing with veg data
harvtall.df <-  pivot_wider(harvtall.df, 
                            names_from = "Food", 
                            values_from = "count")%>%
  replace(is.na(.),0)

sum(harvtall.df[,4:5])##4154.833
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

sum(harvtall.df$PreDensity)#90.26837
sum(harvtall.df$OmnDensity)#46.9515 - predators were  times more dense than omnivore speces

##proportion density
harvtall.df <- harvtall.df%>%
  mutate(., PerOmnivoreDensity = OmnDensity/(OmnDensity+PreDensity)*100,
         PerPredatorDensity = PreDensity/(OmnDensity+PreDensity)*100)

#________________________________________________________________________
#lnlcd Class lme-----------------------------------------------------------
#________________________________________________________________________
#green out models that we don't use- 

#m <- nlme::lme(OmnDensity ~ nlcdClass, random = ~1| siteID/plotID, 
#               data = harvtall.df);summary(m); shapiro.test(resid(m));Anova(m)
#m <- nlme::lme(PreDensity ~ nlcdClass, random = ~1| siteID/plotID, 
#               data = harvtall.df);summary(m); shapiro.test(resid(m));Anova(m)

#m <- nlme::lme(LOGOmnDensity ~ nlcdClass, random = ~1| siteID/plotID, 
#               data = harvtall.df,na.action = na.omit);summary(m); shapiro.test(resid(m));Anova(m)
#m <- nlme::lme(LOGPreDensity ~ nlcdClass, random = ~1| siteID/plotID, 
#               data = harvtall.df,na.action = na.omit);summary(m); shapiro.test(resid(m));Anova(m)

harvtall.df$LOGOmnDensity <- log(harvtall.df$OmnDensity)
harvtall.df <- harvtall.df%>%
  mutate(LOGOmnDensity = na_if(LOGOmnDensity, "-Inf"))

harvtall.df$LOGPreDensity <- log(harvtall.df$PreDensity)
harvtall.df <- harvtall.df%>%
  mutate(LOGPreDensity = na_if(LOGPreDensity, "-Inf"))
#________________________________________________________________________
##nlcd class plots----
#________________________________________________________________________
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
pdf("/Users/JaneyLienau/Desktop/Model_Graphs/p1p2.pdf", width = 7, height = 14)
plot(p1p2)
dev.off()
#________________________________________________________________________
#Tree diversity lme----
#3 variables ~ Tdiv, %EG, %ECM EG sig predictor of all version predator + omnivore except "Predator"
#________________________________________________________________________
m <- nlme::lme(OmnDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(LOGOmnDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGPreDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
#________________________________________________________________________
##%EG plot~----
#________________________________________________________________________
p3<-ggplot(harvtall.df,aes(x=PerEG_BA, y=OmnDensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape =siteID))+
  theme(legend.position = "top")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Omnivore Species Density\n(Individuals/Cup)', color='Site')+
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
  labs(shape = "Site ID", color = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="top",
    legend.direction="horizontal")+
  scale_y_continuous(breaks = seq(0, 4, by = .5))
p3

p4<-ggplot(harvtall.df,aes(x=PerEG_BA, y=PreDensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape = siteID))+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Evergreen Trees\nRelative Abundance (%)', y = 'Predator Species Density\n(Individuals/Cup)', color='Site')+
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
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="none",
    legend.direction="vertical")+
  scale_y_continuous(breaks = seq(0, 4, by = .5))
p4

#cow plot to stack omnivore/predator percent plots
p3p4<-plot_grid(p3, p4, labels = c('A', 'B'), label_size = 20, ncol = 1)
p3p4
pdf("/Users/JaneyLienau/Desktop/Model_Graphs/Supp4_PredatorOmnivore_density.pdf", width = 7, height = 14)
plot(p3p4)
dev.off()

#________________________________________________________________________
#Tree density------
#3 variables ~ Tden, %EG, %ECM - EG sig predictor of predator + omnivore
#________________________________________________________________________
##  Tden, %EG, %ECM ----
#________________________________________________________________________
m <- nlme::lme(OmnDensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGOmnDensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~ totalBA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGPreDensity ~ totalBA + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

#samre results as total basal area
m <- nlme::lme(OmnDensity ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGOmnDensity ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(PreDensity ~ totalStems + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(LOGPreDensity ~ totalStems + PerEG_BA+PerECM_BA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

#________________________________________________________________________
##  Tden----
#________________________________________________________________________
m <- nlme::lme(OmnDensity ~ totalBA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGOmnDensity ~ totalBA,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PreDensity ~totalBA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGPreDensity ~totalBA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(OmnDensity ~ totalStems,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGOmnDensity ~ totalStems,random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(PreDensity ~totalStems, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(LOGPreDensity ~totalStems, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))


##BEFORE YOU RUN THIS CODE!!! Make sure the harvtall.df you use is a version that doesn't exclude the species for the functional analysis
m <- nlme::lme(density ~totalStems, random = ~1| siteID/plotID.x, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m)) #ns
m <- nlme::lme(log(density) ~totalStems, random = ~1| siteID/plotID.x, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))#ns

m <- nlme::lme(density ~totalBA, random = ~1| siteID/plotID.x, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m)) #ns
m <- nlme::lme(log(density) ~totalBA, random = ~1| siteID/plotID.x, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m)) #ns
#________________________________________________________________________
##Tden plot--------
#________________________________________________________________________
p5<-ggplot(harvtall.df,aes(x=totalBA, y=density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape = siteID))+
  theme(legend.position = "top")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance of Trees\n (Basal Area)', y = 'Ground Beetle Density', color='Site')+
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
    labs(shape = "Site ID", color = "Site ID")
p5

pdf("/Users/JaneyLienau/Desktop/Model_Graphs/p5.pdf", width = 7, height = 7)
plot(p5)
dev.off()


p6<-ggplot(harvtall.df,aes(x=totalStems, y=density))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape= siteID))+
  theme(legend.position = "top")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Relative Abundance of Trees\n (Total Stems)', y = 'Ground Beetle Density', color='Site')+
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
  labs(shape = "Site ID", color = "Site ID")
p6

pdf("/Users/JaneyLienau/Desktop/p6.pdf", width = 7, height = 7)
plot(p6)
dev.off()

#________________________________________________________________________
# Prop. by %EG plot----
#________________________________________________________________________
hist(log(harvtall.df$PerOmnivoreDensity))
m <- nlme::lme(PerOmnivoreDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PerPredatorDensity ~ Shan_BA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(PerOmnivoreDensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PerPredatorDensity ~ totalBA + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))

m <- nlme::lme(PerOmnivoreDensity ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))
m <- nlme::lme(PerPredatorDensity ~ totalStems + PerEG_BA+PerECM_BA, random = ~1| siteID/plotID, 
               data = harvtall.df, na.action = na.omit);summary(m); shapiro.test(resid(m))


p7<-ggplot(harvtall.df,aes(x=PerEG_BA, y=PerOmnivoreDensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape = siteID))+
  theme(legend.position = "right")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = NULL, y = 'Omnivore\nRelative Abundance (%)', color='Site')+
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
  labs(shape = "Site ID", color = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="top",
    legend.direction="horizontal")+
  scale_y_continuous(breaks = seq(0, 100, by = 10))
p7

#________________________________________________________________________
# Figure 4: Prop. of Predators ~ %EG 
#________________________________________________________________________

p8<-ggplot(harvtall.df,aes(x=PerEG_BA, y=PerPredatorDensity))+
  geom_smooth(method = 'lm', formula = 'y ~ x') +
  geom_point(aes(color = siteID, shape = siteID))+
  theme(legend.position = "top")+
  scale_fill_brewer(palette="Dark2")+
  labs(x = 'Evergreen Trees\nRelative Abundance (%)', y = 'Predator\nRelative Abundance (%)')+
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
  labs(color = "Site ID",shape = "Site ID")+
  theme(
    axis.text=element_blank(),
    title=element_text(size=rel(1.5)),
    legend.text=element_text(size=rel(1.5)),
    legend.position="top",
    legend.direction="horizontal")+
  scale_y_continuous(breaks = seq(0, 99, by = 10))
p8

pdf("/Users/JaneyLienau/Desktop/Model_Graphs/fig4.pdf", width = 7, height = 7)
plot(p8)
dev.off()

#fig with both in a panel------------------

#cow plot to stack omnivore/predator count plots
p7p8<-plot_grid(p7, p8, labels = c('A', 'B'), label_size = 20, ncol = 1)
p7p8
pdf("/Users/JaneyLienau/Desktop/Model_Graphs/p7p8.pdf", width = 7, height = 14)
plot(p7p8)
dev.off()

##----------------------------------------------------
##-------------Summary Beetle Data Table----------##
##----------------------------------------------------

#BART
sum(BART.df[,2:29])#9373
BART.df <- pivot_longer(BART.df, 2:29, values_to = "count")%>% #reformatting data
  rename(., taxonID = "name")%>%
  select(plotIDyear,
         plotID,
         count,
         taxonID,
         siteID)
sum(BART.df$count) #9373

#JERC
sum(JERC.df[,2:49])#410
JERC.df <- pivot_longer(JERC.df, 2:49, values_to = "count")%>% #reformatting data
  rename(., taxonID = "name")%>%
  select(plotIDyear,
         plotID,
         count,
         taxonID,
         siteID)
sum(JERC.df$count) #410

speciestable.df <- rbind(HARV.df, TALL.df, BART.df, JERC.df)

testspecies <- speciestable.df%>%
  select(-plotIDyear, -plotID)
testspecies$count <-   trunc(testspecies$count)

testspecies <- testspecies%>%
  pivot_wider(names_from = "siteID", values_from = "count", values_fn = list(count = sum))%>%
  replace(is.na(.),0)

sppTable.df <-left_join(testspecies, funspecies.df,by = "taxonID") 

#collecting scientific names from Beetle data
if(TRUE){
  sci.df <- BeetleHARV$bet_expertTaxonomistIDProcessed%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci2.df <- BeetleHARV$bet_expertTaxonomistIDRaw%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci3.df <- BeetleHARV$bet_parataxonomistID%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci4.df <- BeetleHARV$bet_sorting%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sciCombind <- rbind(sci.df, sci2.df, sci3.df, sci4.df)%>%
    distinct()#%>%
  #filter(., !grepl("sp.", scientificName))
  
  #try to find names TALL
  sci5.df <- BeetleTALL$bet_expertTaxonomistIDProcessed%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci6.df <- BeetleTALL$bet_expertTaxonomistIDRaw%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci7.df <- BeetleTALL$bet_parataxonomistID%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci8.df <- BeetleTALL$bet_sorting%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sciCombind2 <- rbind(sci5.df, sci6.df, sci7.df, sci8.df)%>%
    distinct()#%>%
  #filter(., !grepl("sp.", scientificName))
  
  #try to find names JERC
  sci9.df <- BeetleJERC$bet_expertTaxonomistIDProcessed%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci10.df <- BeetleJERC$bet_expertTaxonomistIDRaw%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci11.df <- BeetleJERC$bet_parataxonomistID%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci12.df <- BeetleJERC$bet_sorting%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sciCombind3 <- rbind(sci9.df, sci10.df, sci11.df, sci12.df)%>%
    distinct()#%>%
  #filter(., !grepl("sp.", scientificName))
  
  #try to find names BART
  sci13.df <- BeetleBART$bet_expertTaxonomistIDProcessed%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci14.df <- BeetleBART$bet_expertTaxonomistIDRaw%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci15.df <- BeetleBART$bet_parataxonomistID%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sci16.df <- BeetleBART$bet_sorting%>%
    select(taxonID, scientificName)%>%
    distinct()
  
  sciCombind4 <- rbind(sci13.df, sci14.df, sci15.df, sci16.df)%>%
    distinct()#%>%
  # filter(., !grepl("sp.", scientificName))
  
  names.df <- rbind(sciCombind,sciCombind2, sciCombind3, sciCombind4)%>%
    distinct()
}

sppTable.df <- left_join(sppTable.df, names.df, by = "taxonID")%>%
  select(scientificName,
         taxonID,
         Food,
         HARV,
         TALL,
         BART,
         JERC)%>%
  filter(!taxonID == "HARPEN.HARTEX")%>%
  mutate(total = rowSums(.[4:7]))

#write_csv(sppTable.df, "/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/sppTable.csv")
###-------------------------------------END Table-----------------------------------##


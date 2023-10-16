
#Packages-------------------------------------------------
library(tidyverse)
library(vegan)
library(cowplot)
library(ggrepel)                   # makes data labels not overlap
#Load data-------------------------------------------------
HARV_data.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/HARV_data.csv")
TALL_data.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/TALL_data.csv")
JERC_data.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/JERC_data.csv")
BART_data.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/BART_data.csv")

allTALL.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/allTALL.csv")
allHARV.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/allHARV.csv")
BasalAreaHARV.df <- read.csv("/Users/JaneyLienau/Desktop/GitHubRepository/neon-groundbeetle/Ordination/data/BasalAreaHARV.df.csv")

#TALL-------------------------------------------------
  #ord stats and figures 
ordinationTALL.df <- TALL_data.df%>%
  select(plotID,AGOCON:PASPUN)
ordinationTALL.df <- ordinationTALL.df%>%
  group_by(plotID)%>%
  summarise(across(everything(), sum))%>%
  rowwise()

veg.df<-allTALL.df

ordinationTALL.df<-allTALL.df%>%
  select(plotID,
         nlcdClass,
         PerEG_BA,
         Shan_BA)%>%
  left_join(ordinationTALL.df, by = "plotID")

species.df<-ordinationTALL.df[,5:37]                       
  
  # Remove singletons:
  species.df = species.df[,colSums(species.df) > 1]
  
  my_nmds_result <- species.df%>% 
    vegan::metaMDS()
  
  # plot stress
  my_nmds_result$stress ## [1] 0.08367846 great (>0.05 is excellen <.2 is poor)
  
  data.scores = as.data.frame(scores(my_nmds_result, display = "sites", "species"))
  data.scores$plotID = ordinationTALL.df$plotID #need for plot
  
  ordinationTALL.df<-ordinationTALL.df%>%
    select(!nlcdClass)%>%
    select(!PerEG_BA)%>%
    select(!Shan_BA)
  
  en = envfit(my_nmds_result, ordinationTALL.df, permutations = 999, na.rm = TRUE)
  
  en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #ordiarrowmul = for arrows
  en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en) #plotID - need to keep nlcdClass off
  
  data.scores<-veg.df%>%
    select(plotID,
           nlcdClass,
           PerEG_BA,
           totalBA)%>%
    left_join(., data.scores, by = "plotID")
  
  
  #with o lines
  gg1 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=nlcdClass),alpha=0.30) + # add the convex hulls
    geom_point(data = data.scores, aes(shape = nlcdClass), size = 5)+ 
    geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2, max.overlaps = 20) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(shape = "NLCD Forest Cover Class", 
         fill = "NLCD Forest Cover Class", 
         title = "Talledega Stress = 0.084",#0.08367846
         )+
    theme(plot.title = element_text(hjust = 0.5))
  gg1
  #pdf("/Users/JaneyLienau/Desktop/Model_Graphs/TALL_Ord.pdf", width = 7, height = 7)
  #plot(gg1)
  #dev.off()
  
  # PERMANOVA 
  
  per1<- adonis2(species.df ~ EG, data = data.scores %>% mutate(EG = nlcdClass == "evergreenForest"), permutations = 999, method="bray")
  per1 #sig EG 0.026
  TALL<-as.data.frame(per1)
  TALL$siteID<-"TALL"
  #print
  ordinationTALL.df$siteID<-"TALL"

  #what does the data need to look like?
  #are the variances conparable? variances of what?
  
       
#removing mixed forest
  d <- dist(species.df[-2,])
 mod <-  betadisper(d, data.scores$nlcdClass[-2])
mod
anova(mod)#null should be variances are same, supports
#HARV-------------------------------------------------
  ordinationHARV.df <- HARV_data.df%>%
    select(plotID,ACUHYD:TRIAUT)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    group_by(plotID)%>%
    summarise(across(everything(), sum))%>%
    rowwise()
  
  veg.df<-allHARV.df
  
  ordinationHARV.df<-allHARV.df%>%
    select(plotID,
           nlcdClass,
           PerEG_BA,
           Shan_BA)%>%
    left_join(ordinationHARV.df, by = "plotID")
  
  species.df<-ordinationHARV.df[,5:37]                                         
  sum(species.df)
  # Remove singletons:
  species.df = species.df[,colSums(species.df) > 1]
  
  
  veg.df<-BasalAreaHARV.df
  
  ordinationHARV.df<-ordinationHARV.df%>%
    left_join(veg.df, by = "plotID")
  
  my_nmds_result <- species.df%>% 
    vegan::metaMDS()
  
  # plot stress
  my_nmds_result$stress ## [1] 0.1295496 ok (>0.05 is excellen <.2 is poor)
  goodness(my_nmds_result)
  #inertcomp(my_nmds_result)
  #
  ordinationHARV.df<-ordinationHARV.df%>%
    select(!perECM)%>%
    select(!ShannonPlant)%>%
    select(!SimpsonPlant)%>%
    select(!totalBA)
  
  ordinationHARV.df<-ordinationHARV.df%>%
    select(!Deciduous)%>%
    select(!Evergreen)%>%
    select(!AM)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!nlcdClass.y)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!ECM)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!perEvergreen)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!PerEG_BA)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!Shan_BA)
  
  
  data.scores = as.data.frame(scores(my_nmds_result, display = "sites", "species"))
  data.scores$plotID = ordinationHARV.df$plotID #need for plot
  
  en = envfit(my_nmds_result, ordinationHARV.df, permutations = 999, na.rm = TRUE)
  
  en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #ordiarrowmul = for arrows
  en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en) #plotID - need to keep nlcdClass off
  
  data.scores<-veg.df%>%
    select(plotID,
           nlcdClass,
           perEvergreen)%>%
    left_join(., data.scores)%>%
    filter(!plotID == "HARV_025")
  

  #with o lines
  gg2 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=nlcdClass),alpha=0.30) + # add the convex hulls
    geom_point(data = data.scores, aes(shape = nlcdClass), size = 5)+ 
    geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2, check_overlap= FALSE) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(shape = "NLCD Forest Cover Class", 
         fill = "NLCD Forest Cover Class", 
         title = "Harvard Stress = 0.130")+ #0.1295496
    theme(plot.title = element_text(hjust = 0.5))
  gg2
  #pdf("/Users/JaneyLienau/Desktop/Model_Graphs/HARV_Ord.pdf", width = 7, height = 7)
  #plot(gg2)
  #dev.off()
  
  # PERMANOVA 
  per2<- adonis2(species.df ~ EG, data = data.scores %>% mutate(EG = nlcdClass == "evergreenForest"), permutations = 999, method="bray")
  per2 #sig EG 0.042

  HARV<-as.data.frame(per2)
  HARV$siteID<-"HARV"
  
  #print
  ordinationHARV.df$siteID<-"HARV"
  
  
  #species driving site distribution
  spp.fit <- envfit(my_nmds_result, ordinationHARV.df, permutations = 999)
  head(spp.fit)
  
  d <- dist(species.df)
  mod <-  betadisper(d, data.scores$nlcdClass)
  mod
  anova(mod)#
  permutest(mod, permutations = 99)
  TukeyHSD(mod)
  #showing that there is no evidence that they are unequal, there isn't evidence for different variances, so the tests are robust, 
  #for the TALL forests, one only had 1 point, so we couldn't assess. In Harv, EG was a little less but test says ok

#JERC-------------------------------------------------
  #ord stats and figures 
  ordinationJERC.df <- JERC_data.df%>%
    select(plotID,AGOCON:SELGRA)
  ordinationJERC.df <- ordinationJERC.df%>%
    group_by(plotID)%>%
    summarise(across(everything(), sum))%>%
    rowwise()
  nlcdClass <- JERC_data.df%>%
    select(plotID, nlcdClass)%>%
    distinct()
  
  ordinationJERC.df<-ordinationJERC.df%>%
    left_join(nlcdClass, by = "plotID")
  
  species.df<-ordinationJERC.df[,2:49]                       
  
  # Remove singletons:
  species.df = species.df[,colSums(species.df) > 1]
  
  my_nmds_result <- species.df%>% 
    vegan::metaMDS()
  
  # plot stress
  my_nmds_result$stress ## [1] 0.08736128 great (>0.05 is excellen <.2 is poor)
  
  data.scores = as.data.frame(scores(my_nmds_result, display = "sites", "species"))
  data.scores$plotID = ordinationJERC.df$plotID #need for plot
  
  
  en = envfit(my_nmds_result, ordinationJERC.df, permutations = 999, na.rm = TRUE)
  
  en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #ordiarrowmul = for arrows
  en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en) #plotID - need to keep nlcdClass off
  
  data.scores<-nlcdClass%>%
    left_join(., data.scores, by = "plotID")
  
  
  #with o lines
  gg3 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=nlcdClass),alpha=0.30) + # add the convex hulls
    geom_point(data = data.scores, aes(shape = nlcdClass), size = 5)+ 
    geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2, max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(shape = "NLCD Forest Cover Class", fill = "NLCD Forest Cover Class", title = "The Jones Center At Ichauway Stress = 0.087")+ #0.08736128
    theme(plot.title = element_text(hjust = 0.5))
  gg3
  #pdf("/Users/JaneyLienau/Desktop/Model_Graphs/JERC_Ord.pdf", width = 7, height = 7)
  #plot(gg3)
  #dev.off()
  

  
  per3<- adonis2(species.df ~ EG, data = data.scores %>% mutate(EG = nlcdClass == "evergreenForest"), permutations = 999, method="bray")
  per3 #not sig EG 0.64
  
  JERC<-as.data.frame(per3)
  JERC$siteID<-"JERC"
  #print
  ordinationJERC.df$siteID<-"JERC"


  
  #BART-------------------------------------------------
  #ord stats and figures 
  ordinationBART.df <- BART_data.df%>%
    select(plotID,AGORET:SYNIMP)
  ordinationBART.df <- ordinationBART.df%>%
    group_by(plotID)%>%
    summarise(across(everything(), sum))%>%
    rowwise()
  nlcdClass <- BART_data.df%>%
    select(plotID, nlcdClass)%>%
    distinct()
  
  ordinationBART.df<-ordinationBART.df%>%
    left_join(nlcdClass, by = "plotID")
  
  species.df<-ordinationBART.df[,2:29]                       
  
  # Remove singletons:
  species.df = species.df[,colSums(species.df) > 1]
  
  my_nmds_result <- species.df%>% 
    vegan::metaMDS()
  
  # plot stress
  my_nmds_result$stress ## [1] 0.1024312 great (>0.05 is excellen <.2 is poor)
  
  data.scores = as.data.frame(scores(my_nmds_result, display = "sites", "species"))
  data.scores$plotID = ordinationBART.df$plotID #need for plot
  
  
  en = envfit(my_nmds_result, ordinationBART.df, permutations = 999, na.rm = TRUE)
  
  en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #ordiarrowmul = for arrows
  en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en) #plotID - need to keep nlcdClass off
  
  data.scores<-nlcdClass%>%
    left_join(., data.scores, by = "plotID")
  
  
  #with o lines
  gg4 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=nlcdClass),alpha=0.30) + # add the convex hulls
    geom_point(data = data.scores, aes(shape = nlcdClass), size = 5)+ 
    geom_text_repel(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(shape = "NLCD Forest Cover Class", fill = "NLCD Forest Cover Class", title = "Bartlett Stress = 0.102")+ #0.1024312
    theme(plot.title = element_text(hjust = 0.5))
  gg4
  #pdf("/Users/JaneyLienau/Desktop/Model_Graphs/BART_Ord.pdf", width = 7, height = 7)
  #plot(gg4)
  #dev.off()
  

  per4<- adonis2(species.df ~ EG, data = data.scores %>% mutate(EG = nlcdClass == "evergreenForest"), permutations = 999, method="bray")
  per4 #not sig EG 0.116
  
  BART<-as.data.frame(per4)
  BART$siteID<-"BART"
  
  #print
  ordinationBART.df$siteID<-"BART"
  
  #Supp Figure-----------------------------------------------------
  
  
  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    gg3 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  prow <- plot_grid(
    gg1 + theme(legend.position="none"),
    gg2 + theme(legend.position="none"),
    gg3 + theme(legend.position="none"),
    gg4 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B", "C","D"),
    hjust = -1,
    nrow = 2
  )
  prow
  
  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  p <-  plot_grid(prow, legend, rel_widths = c(2, .4))
  p
  #PDF-------------------------------------------------------------
  pdf("/Users/JaneyLienau/Desktop/Model_Graphs/allNMDS-notFig.pdf", width = 14, height = 10)
  plot(p)
  dev.off()
  
  
  
  #------Figure 3: NMDS of HARV and TALl 
  
  prow1 <- plot_grid(
    gg1 + theme(legend.position="none"),
    gg2 + theme(legend.position="none"),
    align = 'hv',
    labels = c("A", "B"),
    hjust = -1,
    ncol  = 2
  )
  prow1
  
  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  p1 <-  plot_grid(prow1, legend, rel_widths = c(2, .4))
  p1
  
  pdf("/Users/JaneyLienau/Desktop/Model_Graphs/fig3.pdf", width = 12, height = 7)
  plot(p1)
  dev.off()
  
  
  #------Supplemental Figure 1: NMDS of BART and JERC
  
  prow2 <- plot_grid(
    gg3+ theme(legend.position="none"),
    gg4 + theme(legend.position="none"),
    align = 'hv',
    labels = c("A", "B"),
    hjust = -1,
    ncol  = 2
  )
  prow2
  
  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
  p2 <-  plot_grid(prow2, legend, rel_widths = c(2, .4))
  p2
  
  pdf("/Users/JaneyLienau/Desktop/Model_Graphs/Supp1_BARTandJERC_NMDS.pdf", width = 12, height = 7)
  plot(p2)
  dev.off()
  
  
  #---table with summary statistics
  
SummaryStats<- rbind(TALL, BART, JERC, HARV)





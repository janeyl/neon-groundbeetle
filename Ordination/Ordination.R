
# -----------------------------------------------------------
#  ordination
# -----------------------------------------------------------

library(tidyverse)
library(vegan)

paraTax_TALL = read_csv("Ordination/data/paraTax_TALL.csv")
allTALL.df = read_csv("Ordination/data/allTALL.df.csv")

paraTax_HARV = read_csv("Ordination/data/paraTax_HARV.csv")
BasalAreaHARV.df = read_csv("Ordination/data/BasalAreaHARV.df.csv")

if(TRUE){
  #ord stats and figures 
  ordinationTALL.df<-paraTax_TALL%>%                                               
    select(taxonID.y,
           individualCount,
           plotID)%>%
    group_by(plotID, taxonID.y)%>%                                        
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  
  veg.df<-allTALL.df
  
  ordinationTALL.df<-allTALL.df%>%
    select(plotID,
           nlcdClass,
           PerEG_BA,
           Shan_BA,
           totalBA)%>%
    left_join(ordinationTALL.df, by = "plotID")
  #remove species that are 0? because in some cases veg is missing?
  
  ordinationTALL.df<-ordinationTALL.df%>%
    select(!STEPLE)%>%
    select(!CICPUN)%>%
    select(!TETCAR2)%>%
    select(!CRADUB)%>%
    select(!AMBMEX)%>%
    select(!TETVIR)%>%
    select(!HARKAT)%>%
    select(!ACUTES)%>%
    select(!NOTSAY)%>%
    select(!NOTTER)%>%
    select(!LEBPUL)%>%
    select(!SELFOS)%>%
    select(!SELFAT)%>%
    select(!SELPAL)%>%
    select(!SELCON2)
  species.df<-ordinationTALL.df[,6:38]                        
  
  # Remove singletons:
  species.df = species.df[,colSums(species.df) > 1]
  
  my_nmds_result <- species.df%>% 
    vegan::metaMDS()
  
  # plot stress
  my_nmds_result$stress ## [1] 0.0773343 great (>0.05 is excellen <.2 is poor)
  
  
  data.scores = as.data.frame(scores(my_nmds_result))
  data.scores$plotID = ordinationTALL.df$plotID #need for plot
  
  ordinationTALL.df<-ordinationTALL.df%>%
    select(!nlcdClass)%>%
    select(!PerEG_BA)%>%
    select(!Shan_BA) %>%
    select(!totalBA)
  
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
  gg3 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(data = data.scores, aes(colour = PerEG_BA, shape = nlcdClass), size = 3)+ #per EG
    scale_fill_gradient2()+ 
    #geom_point(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), size = .5, colour = "red")+
    geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "red", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2, check_overlap= TRUE) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(colour = "Percent Evergreen", shape = "Site ID")+
    guides(colour = guide_colourbar(order = 1),
           shape = guide_legend(order = 2))
  gg3
  pdf("/Users/yourfilepath/Desktop/Model_Graphs/TALL_Ord3.pdf", width = 7, height = 7)
  plot(gg3)
  dev.off()
  
  # PERMANOVA 
  per1<- adonis2(species.df ~ PerEG_BA, data = data.scores, permutations = 999, method="bray")
  per1 # sig
  per2<- adonis2(species.df ~ plotID, data = data.scores, permutations = 999, method="bray")
  per2
  
  per3<- adonis2(species.df ~ EG, data = data.scores %>% mutate(EG = nlcdClass == "evergreenForest"), permutations = 999, method="bray")
  per3
  
  #print
  ordinationTALL.df$siteID<-"TALL"
  
  
  ordinationHARV.df<-paraTax_HARV%>%                                               
    select(taxonID.y,
           individualCount,
           plotID)%>%
    group_by(plotID, taxonID.y)%>%                                        
    pivot_wider(names_from = "taxonID.y", values_from = "individualCount", values_fn = length, values_fill = 0)
  species.df<-ordinationHARV.df[,2:34]                                            # 2388 species
  
  # Remove singletons:
  species.df = species.df[,colSums(species.df) > 1]
  
  
  veg.df<-BasalAreaHARV.df
  
  ordinationHARV.df<-ordinationHARV.df%>%
    left_join(veg.df, by = "plotID")
  
  my_nmds_result <- species.df%>% 
    vegan::metaMDS()
  
  # plot stress
  my_nmds_result$stress ## [1] 0.1066321 ok (>0.05 is excellen <.2 is poor)
  goodness(my_nmds_result)
  inertcomp(my_nmds_result)
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
    select(!nlcdClass)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!ECM)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!perEvergreen)
  
  
  data.scores = as.data.frame(scores(my_nmds_result))
  data.scores$plotID = ordinationHARV.df$plotID #need for plot
  
  en = envfit(my_nmds_result, ordinationHARV.df, permutations = 999, na.rm = TRUE)
  
  en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en) #ordiarrowmul = for arrows
  en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en) #plotID - need to keep nlcdClass off
  
  data.scores<-veg.df%>%
    select(plotID,
           nlcdClass,
           perEvergreen)%>%
    left_join(., data.scores)%>%
    filter(!plotID=="HARV_025")
  
  
  #with o lines
  gg3 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_point(data = data.scores, aes(colour = perEvergreen, shape = nlcdClass), size = 3)+ #per EG
    scale_fill_gradient2()+ 
    #geom_point(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), size = .5, colour = "red")+
    geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "red", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2, check_overlap= TRUE) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(colour = "Percent Evergreen", shape = "Site ID")+
    guides(colour = guide_colourbar(order = 1),
           shape = guide_legend(order = 2))
  gg3
  pdf("/Users/yourfilepath/Desktop/Model_Graphs/HARV_Ord3.pdf", width = 7, height = 7)
  plot(gg3)
  dev.off()
  
  # PERMANOVA 
  per1<- adonis2(species.df ~ EG, data = data.scores %>% mutate(EG = nlcdClass == "evergreenForest"), permutations = 999, method="bray")
  per1
  #no sig
  per2<- adonis2(species.df ~ perEvergreen, data = data.scores, permutations = 999, method="bray")
  per2
  #no sig
  
  
  #print
  ordinationHARV.df$siteID<-"HARV"
  
  
  #species driving site distribution
  spp.fit <- envfit(my_nmds_result, ordinationHARV.df, permutations = 999)
  head(spp.fit)
}



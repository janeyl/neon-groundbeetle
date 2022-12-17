
#Packages-------------------------------------------------
library(tidyverse)
library(vegan)
library(cowplot)
#Load data-------------------------------------------------
paraTax_TALL = read_csv("Ordination/data/paraTax_TALL.csv")
allTALL.df = read_csv("Ordination/data/allTALL.df.csv")

paraTax_HARV = read_csv("Ordination/data/paraTax_HARV.csv")
BasalAreaHARV.df = read_csv("Ordination/data/BasalAreaHARV.df.csv")


#TALL-------------------------------------------------
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
  my_nmds_result$stress ## [1] 0.07596668 great (>0.05 is excellen <.2 is poor)
  
  data.scores = as.data.frame(scores(my_nmds_result, display = "sites", "species"))
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
    geom_polygon(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=nlcdClass),alpha=0.30) + # add the convex hulls
    geom_point(data = data.scores, aes(shape = nlcdClass), size = 3)+ 
    geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "red", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2, check_overlap= TRUE) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(shape = "Site ID", fill = "Site ID", title = "Talledega Stress = 0.076")+
    theme(plot.title = element_text(hjust = 0.5))
  gg3
  pdf("/Users/JaneyLienau/Desktop/TALL_Ord3.pdf", width = 7, height = 7)
  plot(gg3)
  dev.off()
  
  # PERMANOVA 
  per1<- adonis2(species.df ~ PerEG_BA, data = data.scores, permutations = 999, method="bray")
  per1 # sig
  per2<- adonis2(species.df ~ plotID, data = data.scores, permutations = 999, method="bray")
  per2
  
  per3<- adonis2(species.df ~ EG, data = data.scores %>% mutate(EG = nlcdClass == "evergreenForest"), permutations = 999, method="bray")
  per3 #sig EG 0.036
  
  #print
  ordinationTALL.df$siteID<-"TALL"

#HARV-------------------------------------------------
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
  my_nmds_result$stress ## [1] 0.09850331 ok (>0.05 is excellen <.2 is poor)
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
    select(!nlcdClass)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!ECM)
  
  ordinationHARV.df <- ordinationHARV.df%>%
    select(!perEvergreen)
  
  
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
    filter(!plotID=="HARV_025")
  
  
  #with o lines
  gg4 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=data.scores,aes(x=NMDS1,y=NMDS2,fill=nlcdClass),alpha=0.30) + # add the convex hulls
    geom_point(data = data.scores, aes(shape = nlcdClass), size = 3)+ 
    geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "red", #species name
              fontface = "bold", label = row.names(en_coord_cont), size = 2, check_overlap= TRUE) + 
    theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
          axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
          legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
          legend.text = element_text(size = 9, colour = "grey30")) +
    labs(shape = "Site ID", fill = "Site ID", title = "Harvard Stress = 0.099")+
    theme(plot.title = element_text(hjust = 0.5))
  
  gg4
  pdf("/Users/JaneyLienau/Desktop/HARV_Ord3.pdf", width = 7, height = 7)
  plot(gg4)
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
  
#Supp Figure-----------------------------------------------------
  
  
  # extract the legend from one of the plots
  legend <- get_legend(
    # create some space to the left of the legend
    gg3 + theme(legend.box.margin = margin(0, 0, 0, 12))
  )
  
  prow <- plot_grid(
    gg3 + theme(legend.position="none"),
    gg4 + theme(legend.position="none"),
    align = 'vh',
    labels = c("A", "B"),
    hjust = -1,
    nrow = 1
  )
  prow
 
  # add the legend to the row we made earlier. Give it one-third of 
  # the width of one plot (via rel_widths).
 p <-  plot_grid(prow, legend, rel_widths = c(2, .4))
  p
#PDF-------------------------------------------------------------
  pdf("/Users/JaneyLienau/Desktop/supp.pdf", width = 14, height = 7)
  plot(p)
  dev.off()
  




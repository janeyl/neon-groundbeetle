##making a stats table
#rename log(density), and then manually make the names for each model
### VEG Practice
#library(broom)# use this next time you do stats!

library(kableExtra)
library(tidyr)
library(dplyr)

#________________________________________________________________________
##  MAKE SURE TO LOAD m1-49 in v2JLienau_RCode_GB and m50-61 from v2Functional-Ecology
#________________________________________________________________________

# Create a list to store ANOVA results
anova_results_list <- list()

# Loop through the models (adjust the range as needed)
for (i in 1:61) {
  # Calculate ANOVA for each model (a1 through a49)
  anova_result <- anova(get(paste0("m", i)))
  anova_result$Response <- as.character(formula(get(paste0("m", i)))[[2]])
  
  # Add the ANOVA result to the list
  anova_results_list[[i]] <- data.frame(anova_result)
}
#rename columns to be more specific 


# Combine ANOVA results into a single data frame
combined_df <- do.call(rbind, anova_results_list)

# Format the F value and p-value
format_p_value <- function(f_value, df1, df2, p_value) {
  formatted_f <- sprintf("%.3f", f_value)
  formatted_p <- sprintf("%.3f", p_value)
  return(paste0("F(", df1, ",", df2, ")=", formatted_f, ", p=", formatted_p))
}

combined_df$F_value <- format_p_value(combined_df$`F.value`, combined_df$`numDF`, combined_df$`denDF`, combined_df$`p.value`)

# Create a data frame with terms from row names without incrementing the numbers
combined_df$term <- rownames(combined_df)

new_data <- combined_df %>%
  mutate(term = case_when(
    term == "(Intercept)" ~ "Intercept",
    term == "(Intercept)1" ~ "Intercept",
    # Add more conditions as needed
    TRUE ~ term  # Keep other values unchanged
  ))

combined_df <- combined_df %>%
  mutate(term = sub("\\(Intercept\\)\\d+", "(Intercept)", term))%>%
  mutate(term = sub("PerECM_BA\\d+", "PerECM_BA", term))%>%
  mutate(term = sub("nlcdClass\\d+", "nlcdClass", term))%>%
  mutate(term = sub("PerEG_BA\\d+", "PerEG_BA", term))%>%
  mutate(term = sub("Shan_BA\\d+", "Shan_BA", term))%>%
  mutate(term = sub("Shan_ST\\d+", "Shan_ST", term))%>%
  mutate(term = sub("	Sim_BA\\d+", "	Sim_BA", term))%>%
  mutate(term = sub("Sim_ST\\d+", "Sim_ST", term))%>%
  mutate(term = sub("totalBA\\d+", "totalBA", term))%>%
  mutate(term = sub("totalStems\\d+", "totalStems", term))%>%
  mutate(Order = 1:238)

rm(table_df)
# Create a table using kableExtra
table_df <- combined_df %>%
  select(Response, term, F_value, Order) %>%
  pivot_wider(names_from = "term", values_from = "F_value")%>%
  mutate(Response = ifelse(Response == "log", "density", Response))%>%
  mutate(Response = ifelse(Response == "density", "logDensity", Response))

table_1 <- table_df[1:30,]%>%
  mutate(SiteID = "many")
table_2 <- table_df[31:50,]%>%
  mutate(SiteID = "HARV_BA")
table_3 <- table_df[51:70,]%>%
  mutate(SiteID = "TALL_BA")
table_4 <- table_df[71:90,]%>%
  mutate(SiteID = "HARV_BA")
table_5 <- table_df[91:110,]%>%
  mutate(SiteID = "TALL_BA")
table_6 <- table_df[111:130,]%>%
  mutate(SiteID = "HARV_ST")
table_7 <- table_df[131:150,]%>%
  mutate(SiteID = "TALL_ST")
table_8 <- table_df[151:170,]%>%
  mutate(SiteID = "HARV_ST")
table_9 <- table_df[171:190,]%>%
  mutate(SiteID = "TALL_ST")
table_10 <- table_df[191:238,]%>%
  mutate(SiteID = "Function")

### table NCLD class
table_nlcd <- tabel_first %>% select(Response, Vars_b[2], SiteID) %>% na.omit()->summary_nlcd



##beetle diversity ~ tree diversity 
tabel_first <- table_1%>%
  select(-"Order", -"(Intercept)")

tabel_first <- tabel_first%>%
  mutate(SiteID = ifelse(row_number() <= 6, "nlcdclass", SiteID))%>%
  mutate(SiteID = ifelse(row_number() >= 7 & row_number() <= 18, "treedivBA", SiteID))%>%
  mutate(SiteID = ifelse(row_number() >= 19 & row_number() <= 30, "treedivST", SiteID))

Vars_b <- colnames(tabel_first)

tabel_first %>% select(Response, Vars_b[4], SiteID) %>% na.omit() %>%
  left_join(tabel_first %>% select(Response, Vars_b[2], SiteID) %>% na.omit()) %>%
  left_join(tabel_first %>% select(Response, Vars_b[3], SiteID) %>% na.omit())%>%
  left_join(tabel_first %>% select(Response, Vars_b[5], SiteID) %>% na.omit())%>%
  left_join(tabel_first %>% select(Response, Vars_b[6], SiteID) %>% na.omit())%>%
  left_join(tabel_first %>% select(Response, Vars_b[7], SiteID) %>% na.omit())%>%
  left_join(tabel_first %>% select(Response, Vars_b[8], SiteID) %>% na.omit())%>%
  left_join(tabel_first %>% select(Response, Vars_b[9], SiteID) %>% na.omit())%>%
  left_join(tabel_first %>% select(Response, Vars_b[10], SiteID) %>% na.omit())-> summary_first

#### table most
tabel_most <- rbind(table_2, table_3, table_4, table_5, table_6, table_7, table_8, table_9)%>%
  select(-"Order", -"(Intercept)")
Vars_a <- colnames(tabel_most)

tabel_most %>% select(Response, Vars_a[4], SiteID) %>% na.omit() %>%
  left_join(tabel_most %>% select(Response, Vars_a[5], SiteID) %>% na.omit()) %>%
  left_join(tabel_most %>% select(Response, Vars_a[7], SiteID) %>% na.omit())%>%
  left_join(tabel_most %>% select(Response, Vars_a[10], SiteID) %>% na.omit()) -> summary_most


###table function
tabel_function <- table_10%>%
  select(-"Order", -"(Intercept)",-"nlcdClass",-"Sim_BA",-"Shan_ST",-"Sim_ST")
Vars_c <- colnames(tabel_function)
rm(tabel_function)

tabel_function <- tabel_function%>%
  mutate(SiteID = ifelse(row_number() <= 8, "InduvidualCup", SiteID))%>%
  mutate(SiteID = ifelse(row_number() >= 9 & row_number() <= 16, "normdensity", SiteID))%>%
  mutate(SiteID = ifelse(row_number() >= 17 & row_number() <= 24, "logdensity", SiteID))%>%
  mutate(SiteID = ifelse(row_number() >= 25 & row_number() <= 32, "Stem", SiteID))%>%
  mutate(SiteID = ifelse(row_number() >= 33 & row_number() <= 40, "BA", SiteID))%>%
  mutate(SiteID = ifelse(row_number() >= 41 & row_number() <= 48, "Stem2", SiteID))


######______________NEED TO FIX THIS, NOT RIGHT< 12 models total->need to rename different columsn
tabel_function  %>% select(Response, Vars_c[2], SiteID) %>% na.omit() %>%
  full_join(tabel_function %>% select(Response, Vars_c[3], SiteID) %>% na.omit()) %>%
  full_join(tabel_function %>% select(Response, Vars_c[4], SiteID) %>% na.omit())%>%
  full_join(tabel_function %>% select(Response, Vars_c[5], SiteID) %>% na.omit())%>%
  full_join(tabel_function %>% select(Response, Vars_c[6], SiteID) %>% na.omit())  -> summary_function

names(summary_function)

rm(allmodels)
allmodels <- full_join(summary_nlcd,summary_first)%>%
  full_join(summary_most)
allmodels <- allmodels%>%
  full_join(summary_function)

allmodels <- allmodels%>%
  select(Response, SiteID, nlcdClass, Shan_BA,Sim_BA,Shan_ST,Sim_ST,PerEG_BA,PerECM_BA,totalBA,totalStems)

allmodels <- allmodels%>%
  mutate(Response = ifelse(Response == "ShannonBeetle", "Ground Beetle Diversity (Shannon Index)", Response))%>%
  mutate(Response = ifelse(Response == "SimpsonBeetle", "Ground Beetle Diversity (Simpson Index)", Response))%>%
  mutate(Response = ifelse(Response == "logDensity", "Ground Beetle Density (log) (Individuals/Cup)", Response))%>%
  
  mutate(Response = ifelse(Response == "SYNIMP_relabun", "SYNIMP Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "PASDEP_relabun", "PASDEP Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "logSPHSTE3_relabun", "SPHSTE3 (log) Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "logPTETRI3_relabun", "PTETRI3 (log) Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "logPTEPEN_relabun", "PTEPEN (log) Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "logDICDIL5_relabun", "DICDIL5 (log) Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "logCYCFRE_relabun", "CYCFRE (log) Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "logCARGOR_relabun", "CARGOR (log) Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "logANIHAP_relabun", "ANIHAP (log) Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "CYCCON2_relabun", "CYCCON2 Relative Abundance (%)", Response))%>%
  
  mutate(Response = ifelse(Response == "logANIHAP", "ANIHAP (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logSYN", "SYNIMP (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logSPHSTE3_den", "SPHSTE3 (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logPTETRI3", "PTETRI3 (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logPTEPEN_den", "PTEPEN (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logPASDEP", "PASDEP (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logDICDIL5", "DICDIL5 (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logCYCFRE", "CYCFRE (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logCYCCON2", "CYCCON2 (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "logcCARGOR", "CARGOR (log) (Individuals/Cup)", Response))%>%
  
  mutate(Response = ifelse(Response == "LOGOmnDensity", "Omnivore Density (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "LOGPreDensity", "Predator Density (log) (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "OmnDensity", "Omnivore Density (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "PreDensity", "Predator Density (Individuals/Cup)", Response))%>%
  mutate(Response = ifelse(Response == "PerOmnivoreDensity", "Omnivore Relative Abundance (%)", Response))%>%
  mutate(Response = ifelse(Response == "PerPredatorDensity", "Predator Relative Abundance (%)", Response))

table <- allmodels%>%
  kbl() %>%
  kable_styling()

# Print the table
table

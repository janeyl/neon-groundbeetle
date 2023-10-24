##making a stats table


library(kableExtra)
library(tidyr)


#----
# Create a list to store ANOVA results
anova_results_list <- list()

# Loop through the models (adjust the range as needed)
for (i in 1:49) {
  # Calculate ANOVA for each model (a1 through a49)
  anova_result <- anova(get(paste0("m", i)))
  anova_result$Response <- as.character(formula(get(paste0("m", i)))[[2]])
  
  # Add the ANOVA result to the list
  anova_results_list[[i]] <- data.frame(anova_result)
  
  #add column with
}




# Combine ANOVA results into a single data frame
combined_df <- do.call(rbind, anova_results_list)

# Format the F value and p-value
format_p_value <- function(f_value, df1, df2, p_value) {
  formatted_f <- sprintf("%.3f", f_value)
  formatted_p <- sprintf("%.3f", p_value)
  return(paste0("F(", df1, ",", df2, ")=", formatted_f, ", p =", formatted_p))
}

combined_df$F_value <- format_p_value(combined_df$`F-value`, combined_df$`numDF`, combined_df$`denDF`, combined_df$`p.value`)

# Create a data frame with terms from row names without incrementing the numbers
combined_df$term <- rownames(combined_df)

# Add model names (you can replace these with the actual response variable names)
combined_df$model <- paste("response_variables", 0:48)


# Create a table using kableExtra
table <- combined_df %>%
  select(model, term, F_value) %>%
  pivot_wider(names_from = "model", values_from = "F_value") %>%
  kbl() %>%
  kable_styling()

# Print the table
table


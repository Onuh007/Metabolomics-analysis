library(tidyverse)
library(FactoMineR)
library(factoextra)
library(corrplot)


diabetes_data <- data.frame(
  sample_id  = paste0("P", 1:12),
  group      = c(rep("Healthy", 6), rep("T2D", 6)),
  glucose    = c(4.9,5.2,4.7,5.1,4.8,5.0, 8.3,9.1,7.8,9.8,8.6,9.3),
  lactate    = c(1.0,1.2,0.9,1.1,1.0,1.1, 2.6,2.9,2.4,3.2,2.7,3.0),
  pyruvate   = c(0.4,0.5,0.4,0.5,0.4,0.5, 1.1,1.3,1.0,1.4,1.2,1.3),
  glutamine  = c(2.2,2.4,2.1,2.3,2.2,2.3, 4.1,4.6,3.8,5.0,4.3,4.7),
  alanine    = c(1.3,1.5,1.2,1.4,1.3,1.4, 2.8,3.2,2.6,3.5,2.9,3.1),
  citrate    = c(0.8,1.0,0.7,0.9,0.8,0.9, 1.7,2.0,1.5,2.2,1.8,1.9),
  leucine    = c(1.1,1.3,1.0,1.2,1.1,1.2, 2.3,2.6,2.1,2.8,2.4,2.5),
  valine     = c(1.2,1.4,1.1,1.3,1.2,1.3, 2.5,2.8,2.3,3.0,2.6,2.7),
  tyrosine   = c(0.6,0.7,0.5,0.6,0.6,0.7, 1.4,1.6,1.3,1.8,1.5,1.6),
  isoleucine = c(NA, 1.2,1.0,1.1,1.0,1.1, 2.2,NA, 2.0,2.5,2.2,2.3)
)

diabetes_data          # print the whole table
nrow(diabetes_data)    # number of rows (samples)
ncol(diabetes_data)    # number of columns (variables)
dim(diabetes_data)     # both at once: rows, columns
head(diabetes_data)
tail(diabetes_data)
str(diabetes_data)
summary(diabetes_data)

# from the data, we observe that disease group are usually more than double 
# he value of control. 


# check missing values and Impute missing values
cat("Missing values per metabolite:\n")  # missing data
print(colMeans(is.na(diabetes_data)) * 100)

min_isoleucine <- min(diabetes_data$leucine, na.rm = TRUE)
min_isoleucine

half_min <- min_isoleucine / 2
half_min

diabetes_data$isoleucine[is.na(diabetes_data$isoleucine)] <- half_min

sum(is.na(diabetes_data$isoleucine))   # Confirm no NAs remain
diabetes_data$isoleucine

# Calculate RSD for all metabolite columns
diabetes_matrix <- diabetes_data %>%         # Extract the metabolite matrix 
  dplyr::select(-sample_id, -group)          # (remove sample_id and group)

rsd_values <- apply(diabetes_matrix, 2, function(x) {
  (sd(x) / mean(x)) * 100
})
print(rsd_values)

# No metabolite will be removed at RSD < 30%

log_transformed <- log2(diabetes_matrix + 1)
scaled_matrix <- scale(log_transformed)
scaled_matrix[1:3, ]

# Univariate statistics: t-test and fold change
metabolite_cols <- colnames(diabetes_matrix)
results_df <- data.frame(
  metabolite   = metabolite_cols,
  p_value      = NA,
  mean_control = NA,
  mean_disease = NA
)

for (i in seq_along(metabolite_cols)) {
  met  <- metabolite_cols[i]
  test <- t.test(diabetes_data[[met]] ~ diabetes_data$group)
  results_df$p_value[i]      <- test$p.value
  results_df$mean_control[i] <- test$estimate[1]
  results_df$mean_disease[i] <- test$estimate[2]
}

results_df$fold_change <- results_df$mean_disease / results_df$mean_control
results_df$log2FC      <- log2(results_df$fold_change)
results_df$p_adjusted  <- p.adjust(results_df$p_value, method = "BH")
results_df$significant <- results_df$p_adjusted < 0.05
results_df$neg_log10_p <- -log10(results_df$p_value)

cat("\nStatistical results:\n")
print(results_df)


# Figure 1: Volcano plot
fig1 <- ggplot(results_df,
               aes(x = log2FC, y = neg_log10_p,
                   colour = significant, label = metabolite)) +
  geom_point(size = 4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1),     linetype = "dashed") +
  geom_text(vjust = -0.8, size = 3.5) +
  scale_colour_manual(values = c("grey60", "#D55E00")) +
  labs(title  = "Figure 1: Volcano Plot",
       x      = "Log2 Fold Change",
       y      = "-log10(p-value)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Figure 2: Faceted box plots
long_data <- diabetes_data %>%
  pivot_longer(cols      = all_of(metabolite_cols),
               names_to  = "metabolite",
               values_to = "concentration")

fig2 <- ggplot(long_data,
               aes(x = group, y = concentration, fill = group)) +
  geom_boxplot(alpha = 0.6) +
  geom_jitter(width = 0.1, size = 2) +
  facet_wrap(~ metabolite, scales = "free_y") +
  scale_fill_manual(values = c("Control" = "#0072B2",
                               "Disease" = "#D55E00")) +
  labs(title = "Figure 2: Metabolite Distributions",
       x = "Group", y = "Concentration (mmol/L)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

# Figure 3: Correlation heatmap using corrplot

cor_matrix <- cor(diabetes_matrix)
fig3 <- corrplot(cor_matrix, method = "color", type = "upper",
         tl.col = "black", addCoef.col = "black", 
         number.cex = 0.7, title = "Figure 3: Metabolite Distributions")

print(fig1)
print(fig2)
print(fig3)

# Save all figures
ggsave("figure1_volcano.pdf",      fig1, width = 8, height = 6, dpi = 300)
ggsave("figure2_Faceted box.pdf",  fig2, width = 8, height = 6, dpi = 300)
ggsave("figure3_Correlation heatmap.pdf", fig3, width = 10, height = 8, dpi = 300)


# TASK 6: Multivariate analysis
pca_result <- PCA(scaled_matrix, graph = FALSE)
cat("\nVariance explained:\n")
print(pca_result$eig[1:3, ])

# PCA score plot
fig4 <- fviz_pca_ind(pca_result,
                     col.ind      = diabetes_data$group,
                     palette      = c("steelblue", "tomato"),
                     addEllipses  = TRUE,
                     legend.title = "Group",
                     title        = "Figure 4: PCA Scores Plot")

# PCA loading plot
fig5 <- fviz_pca_var(pca_result,
             col.var = "contrib",
             gradient.cols = c("steelblue", "white", "tomato"),
             repel = TRUE,
             title = "Figure 5:PCA Loadings Plot")

print(fig4)
print(fig5)

# 93.9%


# TASK 7: Biological interpretation
# a) all the tested metabolites are statistically significant
# b) lactate,  pyruvate, alanine,  citrate, leucine, valine, and tyrosine
# c) All metabolic pathways involving all tested metabolites are represented
# d) the shows that is clear difference between Healthy and T2D groups
# it also shows how closely related each member of the group are with each other
# e) the T2D showed higher valuation in their biosynthetic pathway such as 
# amino acid conversion, TCA cycle, and glycolysis. This shows an alteration in 
# disease metabolism of nitrogen leading to reprogramming of metabolic pathway
# where cells shift accordingly leading to higher biosynthesis as observed in 
# the disease group.







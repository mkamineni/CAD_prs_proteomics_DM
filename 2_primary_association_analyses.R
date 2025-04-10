### 1 - Import Libraries and data frames ####
  # 1a - Libraries ####

library(data.table)
library(boot)
library(broom)
library(dplyr)
library(survival)

  # 1b - Load dataframe and make new dataframes with time-varying covariates

outdir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/output/"

df_inc <- fread("/medpop/esp/mkaminen/ukb_proteomics_cvd/input/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")
df_inc$Sex_numeric <- factor(df_inc$Sex_numeric)
df_inc$mergedrace <- factor(df_inc$mergedrace)
df_inc$ever_smoked <- factor(df_inc$ever_smoked)
df_inc$antihtnbase <- factor(df_inc$antihtnbase)
df_inc$cholmed <- factor(df_inc$cholmed)

print(colnames(df_inc))
print(unique(df_inc$alc))

# Make dataframes of all CAD, diabetes, and non-diabetes
#df_cad_dm <- df_inc[df_inc$dm2_prev==1, ]
#df_cad_no_dm <- df_inc[df_inc$dm2_prev==0, ]
df_cad_dm <- df_inc[df_inc$a1c>=48, ]
df_cad_no_dm <- df_inc[df_inc$a1c<48, ]
df_cad_all <- df_inc

print(nrow(df_cad_dm))
print(nrow(df_cad_no_dm))

#proteins <- df_cad_all[, which(colnames(df_cad_all) == "CLIP2"): which(colnames(df_cad_all) == "SCARB2")]
proteins <- colnames(df_cad_all)[which(colnames(df_cad_all) == "CLIP2"): which(colnames(df_cad_all) == "SCARB2")]

#general_feat <- c("cad_prs", "age", "Sex_numeric", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "ever_smoked", "BMI_final", "SBP_final", "antihtnbase", "tchol_final", "hdl_final", "cholmed", "tdi_log_final", "creat_final")
general_feat <- c("cad_prs", "age", "Sex_numeric", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

for (protein in proteins) {
	model_data <- cbind(df_cad_all[, ..general_feat], proteins)

	cad_model <- lm(formula = cad_prs ~ ., data = model_data)

	# Print the coefficients
	cat("Coefficients:\n")	
	coef_cad <- coef(cad_model)
	ci_cad <- confint(cad_model)
	print(ci_cad)
	cad_feat <- names(coef_cad)
}



model_data <- cbind(df_cad_all[, ..general_feat], proteins)

cad_model <- lm(formula = cad_prs ~ ., data = model_data)

# Print the coefficients
cat("Coefficients:\n")
coef_cad <- coef(cad_model)
ci_cad <- confint(cad_model)
print(ci_cad)
cad_feat <- names(coef_cad)

# Print the residuals
cat("\nResiduals:\n")
resid_cad <- residuals(cad_model)

print("finished cad model")

model_data <- cbind(df_cad_dm[, ..general_feat], proteins)

cad_dm_model <- lm(formula = cad_prs ~ ., data = model_data)

# Print the coefficients
cat("Coefficients:\n")
coef_cad_dm <- coef(cad_dm_model)
ci_cad_dm <- confint(cad_dm_model)
print(ci_cad_dm)
cad_dm_feat <- names(coef_cad_dm)

# Print the residuals
cat("\nResiduals:\n")
resid_cad_dm <- residuals(cad_dm_model)

print("finished cad dm model")


proteins_no_dm <- df_cad_no_dm[, which(colnames(df_cad_no_dm) == "CLIP2"): which(colnames(df_cad_no_dm) == "SCARB2")]

general_feat <- c("cad_prs", "age", "Sex_numeric", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "ever_smoked", "BMI_final", "SBP_final", "antihtnbase", "tchol_final", "hdl_final", "cholmed", "tdi_log_final", "creat_final")

model_data_no_dm <- cbind(df_cad_no_dm[, ..general_feat], proteins_no_dm)

cad_no_dm_model <- lm(formula = cad_prs ~ ., data = model_data_no_dm)

# Print the coefficients
cat("Coefficients:\n")
coef_cad_no_dm <- coef(cad_no_dm_model)
ci_cad_no_dm <- confint(cad_no_dm_model)
cad_no_dm_feat <- names(coef_cad_no_dm)


# Print the residuals
cat("\nResiduals:\n")
resid_cad_no_dm <- residuals(cad_no_dm_model)
print("finished cad no dm model")



# Visualize into plots
# Save data into csv
p_values = summary(cad_model)$coefficients[,4]
cad_csv <- data.frame(
  Feat = names(coef_cad),
  Coef = coef_cad,
  CI_Lower = ci_cad[, 1],
  CI_Upper = ci_cad[, 2],
  P_value = p_values,
  FDR = p.adjust(p_values, method = "BH")
)

cad_csv <- cad_csv[order(cad_csv$FDR, decreasing=FALSE),]
write.csv(cad_csv, paste0(outdir, "cad_coef_limited.csv"))

p_values = summary(cad_dm_model)$coefficients[,4]
cad_dm_csv <- data.frame(
  Feat = names(coef_cad_dm),
  Coef = coef_cad_dm,
  CI_Lower = ci_cad_dm[, 1],
  CI_Upper = ci_cad_dm[, 2],
  P_value = p_values,
  FDR = p.adjust(p_values, method = "BH")
)

cad_dm_csv <- cad_dm_csv[order(cad_dm_csv$FDR, decreasing=FALSE),]  
write.csv(cad_dm_csv, paste0(outdir, "cad_dm_coef.csv"))

p_values = summary(cad_no_dm_model)$coefficients[,4]
cad_no_dm_csv <- data.frame(
  Feat = names(coef_cad_no_dm),
  Coef = coef_cad_no_dm,
  CI_Lower = ci_cad_no_dm[, 1],
  CI_Upper = ci_cad_no_dm[, 2],
  P_value = p_values, 
  FDR = p.adjust(p_values, method = "BH")
)

cad_no_dm_csv <- cad_no_dm_csv[order(cad_no_dm_csv$FDR, decreasing=FALSE),]
write.csv(cad_no_dm_csv, paste0(outdir, "cad_no_dm_coef.csv"))

# create a plot to visualize coefficients for features and confidence intervals
png(paste0(outdir, "cad_no_dm_coef.png"))
lower_bound <- ci_cad_no_dm[, 1]
upper_bound <- ci_cad_no_dm[, 2]
plot_cad <- cad_no_dm_csv[cad_no_dm_csv$FDR<0.05, ]

plot(plot_cad, ylim = range(c(plot_cad$CI_Lower, plot_cad$CI_Upper)),
     main = "Coefficients (95% CI)", 
     col = "blue", pch = 19, ylab = "Coefficient Value", xaxt = "n", xlab = "Coefficient Names")

axis(1, at = 1:length(plot_cad$Feat), labels = plot_cad$Feat, las = 2)

# Add arrows for the confidence interval
arrows(1:length(plot_cad), plot_cad$CI_Lower, 1:length(plot_cad), plot_cad$Upper, angle = 90, code = 3, length = 0.1, col = "red")
dev.off()



### 6 - Merge data frames, restructure and write to csv ####
  # 6a - Combine summary stats for base models of all outcomes to one data frame ####

#adj_df <- select(adj_df, -c('std.error', 'statistic'))
#colnames(adj_df) <- c('Protein', 'HR', 'P_Value', 'CI_Lower', 'CI_Upper', 'Outcome')
#adj_df <- adj_df[,c('Protein', 'Outcome', 'HR', 'CI_Lower', 'CI_Upper', 'P_Value')]

  # 6c - Write to csv ####

#write.csv(adj_df, '/medpop/esp/mkaminen/ukb_proteomics_cvd/output/1_adj_no_prev_timevar_model_hr.csv')


### 1 - Import Libraries and data frames ####
# 1a - Libraries ####

library(data.table)
library(boot)
library(broom)
library(dplyr)
library(survival)

# 1b - Load dataframe and make new dataframes with time-varying covariates

outdir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/output/"

df_inc <- fread("/medpop/esp/mkaminen/ukb_proteomics_cvd/input/ukb_proteomics_baseline_excl_and_imput_withprevcvd.tsv.gz")
df_inc$Sex_numeric <- factor(df_inc$Sex_numeric)
df_inc$mergedrace <- factor(df_inc$mergedrace)
df_inc$ever_smoked <- factor(df_inc$ever_smoked)
df_inc$antihtnbase <- factor(df_inc$antihtnbase)
df_inc$cholmed <- factor(df_inc$cholmed)

print(nrow(df_inc))
print(colnames(df_inc))
print(unique(df_inc$alc))

# age, cad_prev, cad_prs
# want to find associatons between CAD ~ CAD PRS, CAD ~ age, CAD ~ sex

# make final data frame
columns = c("covariate", "coef", "lower", "upper", "pval")
all_coefs <- data.frame(matrix(nrow=0, ncol=length(columns)))
colnames(all_coefs) = columns

# find association between CAD and CAD PRS
cad_prs_model <- lm(as.formula(paste("cad_prev ~ cad_prs")), data = df_inc)
summary_model <- summary(cad_prs_model)
print(summary_model)
all_coef <- summary_model$coefficients["cad_prs", "Estimate"]
all_pval <- summary_model$coefficients["cad_prs", "Pr(>|t|)"]
all_confint <- confint(cad_prs_model)
all_lower <- all_confint["cad_prs", 1]
all_upper <- all_confint["cad_prs", 2]
all_coefs[nrow(all_coefs)+1, ] <- c("cad_prs", all_coef, all_lower, all_upper, all_pval)

# find association between CAD and age
cad_age_model <- lm(as.formula(paste("cad_prev ~ age")), data = df_inc)
summary_model <- summary(cad_age_model)
print(summary_model)
all_coef <- summary_model$coefficients["age", "Estimate"]
all_pval <- summary_model$coefficients["age", "Pr(>|t|)"]
all_confint <- confint(cad_age_model)
all_lower <- all_confint["age", 1]
all_upper <- all_confint["age", 2]
all_coefs[nrow(all_coefs)+1, ] <- c("age", all_coef, all_lower, all_upper, all_pval)

# find association between CAD and sex
cad_sex_model <- lm(as.formula(paste("cad_prev ~ Sex_numeric")), data = df_inc)
summary_model <- summary(cad_sex_model)
print(summary_model)
all_coef <- summary_model$coefficients["Sex_numeric1", "Estimate"]
all_pval <- summary_model$coefficients["Sex_numeric1", "Pr(>|t|)"]
all_confint <- confint(cad_sex_model)
all_lower <- all_confint["Sex_numeric1", 1]
all_upper <- all_confint["Sex_numeric1", 2]
all_coefs[nrow(all_coefs)+1, ] <- c("sex", all_coef, all_lower, all_upper, all_pval)

write.csv(all_coefs, paste0(outdir, "baseline_assoc_sens.csv"))


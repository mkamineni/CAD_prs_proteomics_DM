# 1a - Libraries ####
library(data.table)
library(boot)
library(broom)
library(dplyr)
library(survival)
library(splines)


# read in cohort
df_inc <- fread("/medpop/esp/mkaminen/ukb_proteomics_cvd/input/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")
df_inc$Sex_numeric <- factor(df_inc$Sex_numeric)
df_inc$mergedrace <- factor(df_inc$mergedrace)
df_inc$ever_smoked <- factor(df_inc$ever_smoked)
df_inc$antihtnbase <- factor(df_inc$antihtnbase)
df_inc$cholmed <- factor(df_inc$cholmed)

# limit to diabetic cohort
df_inc$dm <- ifelse(df_inc$a1c >= 48 | df_inc$dm2_prev==1, 1, 0)
df_inc <- df_inc[df_inc$dm == 1, ]
print(nrow(df_inc))
#print(colnames(df_inc))

# print out distribution of L1CAM
print(summary(df_inc$L1CAM))

# print out distribution of CAD PRS
print(summary(df_inc$cad_prs))

# print out distribution of CAD
print(summary(df_inc$cad_inc))

# print out correlation of L1CAM and CAD PRS

# make linear regression with L1CAM and CAD PRS as covariates
cad_model <- glm(as.formula("cad_inc ~ cad_prs + L1CAM + age + Sex_numeric + cad_prs:L1CAM"), family = binomial, data = df_inc)
summary_model <- summary(cad_model)
ci <- confint.default(cad_model)

features = c("cad_prs:L1CAM", "cad_prs", "L1CAM")

for (feature in features) {
	print(feature)
	est <- summary_model$coefficients[feature, "Estimate"]
	se <- summary_model$coefficients[feature, "Std. Error"]
	pval <- summary_model$coefficients[feature, "Pr(>|z|)"]
	ci_lower <- ci[feature, 1]
	ci_upper <- ci[feature, 2]
	print(paste("coeff estimate", as.character(est)))
	print(paste("coeff lower", as.character(ci_lower)))
        print(paste("coeff upper", as.character(ci_upper)))
        print(paste("coeff pval", as.character(pval)))

}

# so this means L1CAM is not modifying the effect of CAD PRS on CAD but is it mediating it

# Step 2: Outcome model (PRS + L1CAM → CAD)
model <- glm(cad_inc ~ cad_prs +  age + Sex_numeric, family = binomial, data = df_inc)

model_L1CAM <- glm(cad_inc ~ cad_prs + L1CAM + age + Sex_numeric, family = binomial, data = df_inc)

print(coef(summary(model))["cad_prs", ])
print(coef(summary(model_L1CAM))["cad_prs", ])



# In diabetics only
df_inc$high_cad_prs  <- ifelse(df_inc$cad_prs > median(df_inc$cad_prs), 1, 0)
df_inc$low_L1CAM <- ifelse(df_inc$L1CAM < median(df_inc$L1CAM), 1, 0)
df_inc$group <- factor(interaction(df_inc$high_cad_prs, df_inc$low_L1CAM),
			levels = c("0.0", "0.1", "1.0", "1.1"),
			labels = c("LowPRS_HighL1CAM", "LowPRS_LowL1CAM", "HighPRS_HighL1CAM", "HighPRS_LowL1CAM"))
model <- glm(cad_inc ~ group + age + Sex_numeric, family = binomial, data = df_inc)
print(exp(cbind(OR = coef(model), confint(model))))


#:q## 1 - Import Libraries and data frames ####
  # 1a - Libraries ####

library(data.table)
library(boot)
library(broom)
library(dplyr)
library(survival)
library(splines)

# Helper functions
# Function to convert cad_prs to spline in any dataset
add_spline_formula <- function(df, outcome, spline_var = "cad_prs", interact_var = NULL) {
	# Get all predictors except outcome
	predictors <- setdiff(names(df), outcome)
      
        # Replace spline_var with ns(spline_var, df=4)
	predictors <- sapply(predictors, function(x) {
		if(x == spline_var) {
			paste0("ns(", spline_var, ", df=4)")
		} else {
			x
		}
	})
									  
	# Combine predictors
	formula_text <- paste(outcome, "~", paste(predictors, collapse = " + "))
	
	# add interaction term separately if needed
	if (!is.null(interact_var)) {
		formula_text <- paste(formula_text, "+ ns(", spline_var, ", df=4):", interact_var)
	}

	return(as.formula(formula_text))
}

# Helper function to extract "cumulative" effect
extract_effect <- function(model, coef_pattern, spline=FALSE){
	summary_model <- summary(model)
  
	if(spline){
		# Find coefficients matching the pattern
		idx <- grep(coef_pattern, rownames(summary_model$coefficients))
    
		# cumulative coefficient
		coef_vec <- coef(model)[idx]
		vcov_mat <- vcov(model)[idx, idx]
		est <- sum(coef_vec)
    
		# approximate SE and CI
		se <- if(length(idx) == 1) sqrt(vcov_mat) else sqrt(sum(vcov_mat))
		ci_lower <- est - 1.96 * se
		ci_upper <- est + 1.96 * se
    
		# approximate p-value using chi-squared (Wald test)
		chisq <- if(length(idx) == 1){
			(coef_vec^2 / vcov_mat)[1]
		} else {
			(coef_vec %*% solve(vcov_mat) %*% coef_vec)[1,1]
		}
		pval <- pchisq(chisq, df=length(idx), lower.tail=FALSE)
    	} else {
		est <- summary_model$coefficients[coef_pattern, "Estimate"]
		se <- summary_model$coefficients[coef_pattern, "Std. Error"]
		pval <- summary_model$coefficients[coef_pattern, "Pr(>|t|)"]
		ci <- confint(model)
		ci_lower <- ci[coef_pattern, 1]
		ci_upper <- ci[coef_pattern, 2]
	}
  	return(list(estimate=est, se=se, ci_lower=ci_lower, ci_upper=ci_upper, pval=pval))
}

# 1b - Load dataframe and make new dataframes with time-varying covariates
spline <- TRUE

outdir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/output/"
if (spline) {
	outdir <- paste0(outdir, "spline_")
}

df_inc <- fread("/medpop/esp/mkaminen/ukb_proteomics_cvd/input/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")
df_inc$Sex_numeric <- factor(df_inc$Sex_numeric)
df_inc$mergedrace <- factor(df_inc$mergedrace)
df_inc$ever_smoked <- factor(df_inc$ever_smoked)
df_inc$antihtnbase <- factor(df_inc$antihtnbase)
df_inc$cholmed <- factor(df_inc$cholmed)

print(nrow(df_inc))
print(colnames(df_inc))
print(unique(df_inc$alc))

# Make dataframes of all CAD, diabetes, and non-diabetes
#df_cad_dm <- df_inc[df_inc$dm2_prev==1, ]
#df_cad_no_dm <- df_inc[df_inc$dm2_prev==0, ]
df_cad_all <- df_inc
df_cad_all$dm <- ifelse(df_cad_all$a1c >= 48 | df_cad_all$dm2_prev==1, 1, 0)
df_cad_dm <- df_cad_all[df_cad_all$dm == 1, ]
df_cad_no_dm <- df_cad_all[df_cad_all$dm == 0, ]

print(nrow(df_cad_dm))
print(nrow(df_cad_no_dm))

#proteins <- df_cad_all[, which(colnames(df_cad_all) == "CLIP2"): which(colnames(df_cad_all) == "SCARB2")]
proteins <- colnames(df_cad_all)[which(colnames(df_cad_all) == "CLIP2"): which(colnames(df_cad_all) == "SCARB2")]

#general_feat <- c("cad_prs", "age", "Sex_numeric", "mergedrace", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10", "ever_smoked", "BMI_final", "SBP_final", "antihtnbase", "tchol_final", "hdl_final", "cholmed", "tdi_log_final", "creat_final")
general_feat <- c("cad_prs", "age", "Sex_numeric", "PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

# make final data frame
columns = c("protein", "all_coef", "all_lower", "all_upper", "all_pval", "all_inter_coef", "all_inter_lower", "all_inter_upper", "all_inter_pval", "dm_coef", "dm_lower", "dm_upper", "dm_pval", "nodm_coef", "nodm_lower", "nodm_upper", "nodm_pval")
all_prot_coefs <- data.frame(matrix(nrow=0, ncol=length(columns)))
colnames(all_prot_coefs) = columns
num_proteins <- length(proteins)
print(num_proteins)
for (protein in proteins) {
        feat <- append(general_feat, protein)
        protein_safe <- make.names(protein)
        
	if (protein_safe %in% colnames(df_cad_all)){
		
		# all cad model
		cad_feat <- append(feat, "dm")
		model_data <- df_cad_all[, ..cad_feat]
	       
		cad_model <- if(spline){lm(add_spline_formula(model_data, outcome = protein_safe, interact_var = "dm"), data = model_data)} else {lm(as.formula(paste(protein_safe, "~ . + cad_prs:dm")), data = model_data)}
    
		res <- extract_effect(cad_model, coef_pattern = if(spline) "^ns\\(cad_prs" else "cad_prs", spline=spline)
		all_coef <- res$estimate
		all_pval <- res$pval
		all_lower <- res$ci_lower
		all_upper <- res$ci_upper
		    
		print(rownames(summary(cad_model)$coefficients))
		res_inter <- extract_effect(cad_model, coef_pattern = if(spline) ":dm$" else "cad_prs:dm", spline=spline)
		all_inter_coef <- res_inter$estimate
		all_inter_pval <- res_inter$pval
		all_inter_lower <- res_inter$ci_lower
		all_inter_upper <- res_inter$ci_upper

		dm_model_data <- df_cad_dm[, ..feat]
        	dm_model <- if(spline){lm(add_spline_formula(dm_model_data, outcome = protein_safe), data = dm_model_data)} else {lm(as.formula(paste(protein_safe, "~ .")), data = dm_model_data)}
		res_dm <- extract_effect(dm_model, coef_pattern = if(spline) "^ns\\(cad_prs" else "cad_prs", spline=spline)
		dm_coef <- res_dm$estimate
		dm_pval <- res_dm$pval
		dm_lower <- res_dm$ci_lower
		dm_upper <- res_dm$ci_upper

		# no dm model
		nodm_model_data <- df_cad_no_dm[, ..feat]
		nodm_model <- if(spline){lm(add_spline_formula(nodm_model_data, outcome = protein_safe), data = nodm_model_data)} else {lm(as.formula(paste(protein_safe, "~ .")), data = nodm_model_data)}
		res_nodm <- extract_effect(nodm_model, coef_pattern = if(spline) "^ns\\(cad_prs" else "cad_prs", spline=spline)
		nodm_coef <- res_nodm$estimate
		nodm_pval <- res_nodm$pval
		nodm_lower <- res_nodm$ci_lower
		nodm_upper <- res_nodm$ci_upper

        	# add new row to final df
		all_prot_coefs[nrow(all_prot_coefs) + 1,] = list(protein_safe, all_coef, all_lower, all_upper, all_pval, all_inter_coef, all_inter_lower, all_inter_upper, all_inter_pval, dm_coef, dm_lower, dm_upper, dm_pval, nodm_coef, nodm_lower, nodm_upper, nodm_pval)
	} else {
		print(protein_safe)
	}
}

all_prot_coefs <- all_prot_coefs[order(all_prot_coefs$dm_pval, decreasing=FALSE),]

write.csv(all_prot_coefs, paste0(outdir, "all_dm_cohort_cad.csv"))

all_prot_coefs <- all_prot_coefs[order(all_prot_coefs$nodm_pval, decreasing=FALSE),]

write.csv(all_prot_coefs, paste0(outdir, "all_nodm_cohort_cad.csv"))




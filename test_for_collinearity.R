library(splines)
library(data.table)

outdir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/output/"

df_cad_all <- fread("/medpop/esp/mkaminen/ukb_proteomics_cvd/input/ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz")
df_cad_all$dm <- ifelse(df_cad_all$a1c >= 48 | df_cad_all$dm2_prev==1, 1, 0)

proteins <- c("CD34", "EPCAM", "GPA33", "ITGA11", "L1CAM", "MCAM", "CTSB", "F3", "GCG", "CNTN5", "NTRK3", "PLTP", "CD58", "GDF15", "OMD", "PCSK9", "NPY")

pc_terms <- paste(paste0("PC", 1:10), collapse = " + ")

results <- data.frame(protein=character(), anova_pval=numeric(), recommendation=character(), stringsAsFactors=FALSE)

for (protein in proteins) {
	protein_safe <- make.names(protein)
	if (protein_safe %in% colnames(df_cad_all)) {
		formula_linear <- as.formula(paste(protein_safe, "~ cad_prs + dm + cad_prs:dm + age + Sex_numeric +", pc_terms))
			        
		formula_spline <- as.formula(paste(protein_safe, "~ ns(cad_prs, df=3) + dm + ns(cad_prs, df=3):dm + age + Sex_numeric +", pc_terms))

		linear_mod <- lm(formula_linear, data = df_cad_all)
		spline_mod <- lm(formula_spline, data = df_cad_all)

		anova_res <- anova(linear_mod, spline_mod)
		pval <- anova_res$`Pr(>F)`[2]
		rec <- ifelse(pval < 0.05, "spline", "linear")

		results[nrow(results)+1, ] <- list(protein_safe, round(pval, 6), rec)
		cat(protein_safe, "- anova p:", round(pval, 4), "->", rec, "\n")
	}
}

write.csv(results, paste0(outdir, "linearity_test.csv"), row.names=FALSE)


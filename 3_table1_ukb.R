# Table 1 (UK Biobank column) + Results "Study Population" numbers
# Uses measured (non-imputed) covariates; reports n missing for each.
library(data.table)

indir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/input/"
outdir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/output/"

df <- fread(paste0(indir, "ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz"))
df$dm <- ifelse(df$a1c >= 48 | df$dm2_prev == 1, 1, 0)

cat("\n### Sample ####\n")
cat("Rows in cohort file:", nrow(df), "\n")
print(table(dm = df$dm, useNA = "ifany"))
cat("-> Analytic sample = non-missing dm (same as the regressions in script 2)\n")
df <- df[!is.na(df$dm), ]
cat("Analytic N:", nrow(df), " | diabetes:", sum(df$dm == 1),
	sprintf("(%.1f%%)\n", 100 * mean(df$dm == 1)))

cat("\n### Codings to check before using the table ####\n")
cat("Sex_numeric (UKB field 31: 0 = female, 1 = male):\n"); print(table(df$Sex_numeric, useNA = "ifany"))
cat("race:\n"); print(table(df$race, useNA = "ifany"))
cat("mergedrace:\n"); print(table(df$mergedrace, useNA = "ifany"))
cat("ever_smoked:\n"); print(table(df$ever_smoked, useNA = "ifany"))
cat("Other candidate columns:\n")
print(grep("htn|hypert|lipid|statin|chol|smok|white|ethnic", colnames(df), value = TRUE, ignore.case = TRUE))

# Derived definitions (check these match what you want to report)
df$female <- as.integer(df$Sex_numeric == 0)
df$white <- as.integer(df$mergedrace == names(which.max(table(df$mergedrace))))  # assumes most common category = White; confirm above
df$htn <- as.integer(df$SBP >= 140 | df$DBP >= 90 | df$antihtnbase == 1)
df$hld <- as.integer(df$cholmed == 1 | df$ldl >= 4.9)  # LDL >= 190 mg/dL

cont <- c(Age = "age", BMI = "BMI", SBP = "SBP", DBP = "DBP", TotalChol = "tchol",
	LDL = "ldl", HDL = "hdl", TG = "tg", CAD_PRS = "cad_prs")
bin <- c(Female = "female", White = "white", EverSmoked = "ever_smoked",
	Hypertension = "htn", Hyperlipidemia = "hld", LipidMed = "cholmed",
	AntiHTNMed = "antihtnbase", Diabetes = "dm")

summ <- function(d) {
	out <- c(N = format(nrow(d), big.mark = ","))
	for (v in names(cont)) {
		x <- as.numeric(d[[cont[v]]])
		out[v] <- sprintf("%.1f (%.1f) [missing %d]", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE), sum(is.na(x)))
	}
	for (v in names(bin)) {
		x <- as.numeric(as.character(d[[bin[v]]]))
		out[v] <- sprintf("%s (%.1f%%) [missing %d]", format(sum(x == 1, na.rm = TRUE), big.mark = ","),
			100 * mean(x == 1, na.rm = TRUE), sum(is.na(x)))
	}
	out
}

tab <- data.frame(Characteristic = c("N", names(cont), names(bin)),
	All = summ(df), Diabetes = summ(df[df$dm == 1, ]), NoDiabetes = summ(df[df$dm == 0, ]),
	row.names = NULL)

cat("\n### Table 1: mean (SD) or n (%) ####\n")
print(tab, right = FALSE)
cat("Units: lipids in mmol/L (UKB native); BMI kg/m2; BP mmHg.\n")
write.csv(tab, paste0(outdir, "table1_ukb.csv"), row.names = FALSE)

cat("\n### Exact p-values the manuscript still needs ####\n")
res <- fread(paste0(outdir, "all_dm_cohort_cad.csv"))
cat("Output file date:", format(file.info(paste0(outdir, "all_dm_cohort_cad.csv"))$mtime), "\n")
print(res[protein %in% c("GDF15", "GCG", "L1CAM"),
	.(protein, dm_coef, dm_pval, all_inter_coef, all_inter_lower, all_inter_upper, all_inter_pval)])

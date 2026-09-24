# All remaining UK Biobank numbers for the manuscript, in one run.
# Printed output is also saved to output/manuscript_numbers.txt
#   1 - Sample flow (for the exclusion flowchart / Figure 1)
#   2 - Analytic sample and diabetes counts
#   3 - Table 1 (all, diabetes, no diabetes)
#   4 - CAD PRS by diabetes status
#   5 - Protein missingness before imputation
#   6 - Checks on the saved association output (exact p-values, skipped proteins)
#   7 - Incident CAD and L1CAM mediation (events, interaction, mediation, 4-group ORs)
# Uses measured (non-imputed) covariates for Table 1; reports n missing for each.

library(data.table)
library(dplyr)

indir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/input/"
outdir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/output/"
run_flow <- TRUE     # section 1 re-reads the raw file and redoes the kinship step (slower)
n_boot <- 1000       # bootstrap draws for the mediation CI
set.seed(1)

sink(paste0(outdir, "manuscript_numbers.txt"), split = TRUE)
pcs <- paste0("PC", 1:10)

### 1 - Sample flow ####
cat("\n### 1 - Sample flow ####\n")
if (run_flow) {
	library(ukbtools)
	raw <- fread(paste0(indir, "ukb_proteomics_baseline_noexcl_noimput.tsv.gz"))
	prs <- fread(paste0(indir, "CADprs.tsv.gz"))
	prs <- prs[prs[[1]] > 0]
	setnames(prs, c("IND_ID", "MixRawPRS_CAD2ndLayer"), c("id", "cad_prs"), skip_absent = TRUE)
	cat("Olink baseline file:", nrow(raw), "\n")
	raw <- merge(raw, prs[, .(id, cad_prs)], by = "id")
	cat("With CAD PRS:", nrow(raw), "\n")
	idx <- c(which(colnames(raw) == "id"), which(colnames(raw) == "CLIP2"):which(colnames(raw) == "NPM1"))
	b <- raw[, ..idx]
	cat("Proteins on platform:", ncol(b) - 1, "\n")
	b <- b %>% dplyr::select(where(~mean(is.na(.)) <= 0.1))
	cat("Proteins after >10% missingness exclusion:", ncol(b) - 1, "\n")
	b <- b[rowMeans(is.na(b[, -1])) <= 0.1, ]
	raw <- raw[raw$id %in% b$id, ]
	cat("After excluding people with >10% protein missingness:", nrow(raw), "\n")
	raw <- raw[!is.na(raw$race), ]
	cat("After excluding missing self-reported race:", nrow(raw), "\n")
	raw <- raw[!is.na(raw$PC1), ]
	cat("After excluding missing genetic PCs:", nrow(raw), "\n")
	rel <- fread(paste0(indir, "ukb7089_rel_s488363.dat"))
	set.seed(1)
	rm_ids <- ukb_gen_samples_to_remove(rel, raw$id, cutoff = 0.0884)
	raw <- raw[!(raw$id %in% rm_ids), ]
	cat("After excluding relatives:", nrow(raw), "\n")
	cat("Prevalent CAD excluded:", sum(raw$cad_prev == 1, na.rm = TRUE), "\n")
	raw <- raw[raw$cad_prev == 0, ]
	cat("After excluding prevalent CAD:", nrow(raw), "\n")
	rm(raw, b, rel, prs); gc()
}

### 2 - Analytic sample ####
df <- fread(paste0(indir, "ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz"))
df$dm <- ifelse(df$a1c >= 48 | df$dm2_prev == 1, 1, 0)

cat("\n### 2 - Analytic sample ####\n")
cat("Rows in cohort file (after A1c merge):", nrow(df), "\n")
print(table(dm = df$dm, useNA = "ifany"))
cat("A1c missing:", sum(is.na(df$a1c)), " | dm2_prev == 1:", sum(df$dm2_prev == 1, na.rm = TRUE),
	" | A1c >= 48 only (no dm2_prev):", sum(df$a1c >= 48 & df$dm2_prev != 1, na.rm = TRUE), "\n")
df <- df[!is.na(df$dm), ]
cat("Analytic N (non-missing diabetes status, as in script 2):", nrow(df),
	" | diabetes:", sum(df$dm == 1), sprintf("(%.1f%%)\n", 100 * mean(df$dm == 1)))

### 3 - Table 1 ####
cat("\n### 3 - Table 1 ####\n")
cat("Codings to check:\n")
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
	LDL = "ldl", HDL = "hdl", TG = "tg", CAD_PRS_raw = "cad_prs")
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
cat("\nTable 1: mean (SD) or n (%)\n")
print(tab, right = FALSE)
cat("Units: lipids in mmol/L (UKB native); BMI kg/m2; BP mmHg.\n")
write.csv(tab, paste0(outdir, "table1_ukb.csv"), row.names = FALSE)

### 4 - CAD PRS by diabetes status ####
cat("\n### 4 - CAD PRS by diabetes status ####\n")
cat("Raw CAD PRS: mean", mean(df$cad_prs), "SD", sd(df$cad_prs), "\n")
df$cad_prs_std <- as.numeric(scale(df$cad_prs))
print(df[, .(n = .N, mean_prs_std = mean(cad_prs_std), sd_prs_std = sd(cad_prs_std)), by = dm])
print(t.test(cad_prs_std ~ dm, data = df))
fit <- lm(as.formula(paste("cad_prs_std ~ dm + age + factor(Sex_numeric) +", paste(pcs, collapse = " + "))), data = df)
cat("PRS difference (SD units), diabetes vs not, adjusted for age, sex, PCs:\n")
print(cbind(coef(summary(fit))["dm", , drop = FALSE], confint(fit)["dm", , drop = FALSE]))

### 5 - Protein missingness before imputation ####
cat("\n### 5 - Protein missingness before imputation ####\n")
mf <- paste0("/medpop/esp/mkaminen/ukb_proteomics_cvd/", "missingness_basefile.csv")
if (file.exists(mf)) {
	miss <- fread(mf)
	prot_cols <- colnames(df)[which(colnames(df) == "CLIP2"):which(colnames(df) == "SCARB2")]
	m <- miss[colname %in% prot_cols]
	cat("File date:", format(file.info(mf)$mtime), "| proteins matched:", nrow(m), "\n")
	cat(sprintf("Overall %% of protein values imputed: %.2f%%\n", 100 * sum(m$missing_count) / (nrow(m) * nrow(df))))
	cat(sprintf("Mean per-protein missingness: %.2f%% (max %.2f%%)\n", 100 * mean(m$missing_proportion), 100 * max(m$missing_proportion)))
} else cat("missingness_basefile.csv not found\n")

### 6 - Checks on the saved association output ####
cat("\n### 6 - Saved association output ####\n")
of <- paste0(outdir, "all_dm_cohort_cad.csv")
res <- fread(of)
cat("File:", of, "| date:", format(file.info(of)$mtime), "| rows:", nrow(res), "\n")
prot_cols <- colnames(df)[which(colnames(df) == "CLIP2"):which(colnames(df) == "SCARB2")]
cat("Protein columns in cohort file:", length(prot_cols), "\n")
cat("Proteins missing from output (hyphenated names skipped by make.names):",
	paste(setdiff(make.names(prot_cols), res$protein), collapse = ", "), "\n")
thr <- 0.05 / length(prot_cols)
hits <- res[dm_pval < thr][order(dm_pval)]
cat("Hits at 0.05/", length(prot_cols), ": ", nrow(hits), "\n", sep = "")
print(hits[, .(protein, dm_coef, dm_lower, dm_upper, dm_pval, all_inter_coef, all_inter_lower, all_inter_upper, all_inter_pval)])
cat("Interaction significant at 0.05/", nrow(hits), ": ",
	paste(hits[all_inter_pval < 0.05 / nrow(hits), protein], collapse = ", "), "\n", sep = "")

### 7 - Incident CAD and L1CAM mediation (participants with diabetes) ####
cat("\n### 7 - Incident CAD and L1CAM mediation ####\n")
inc_col <- grep("^cad_inc$|^cad_incident$|^inc_cad$|^incident_cad$", colnames(df), value = TRUE, ignore.case = TRUE)[1]
cat("Candidate outcome / follow-up columns:\n")
print(grep("inc|time|fu_|follow|censor|event", colnames(df), value = TRUE, ignore.case = TRUE))

if (is.na(inc_col)) {
	cat("No incident CAD column found; set inc_col by hand from the list above and rerun section 7.\n")
} else {
	cat("Using incident CAD column:", inc_col, "\n")
	d <- df[df$dm == 1 & !is.na(df[[inc_col]]), ]
	d$inc <- as.integer(d[[inc_col]] == 1)
	d$sex <- factor(d$Sex_numeric)
	cat("Participants with diabetes:", nrow(d), "| incident CAD events:", sum(d$inc),
		sprintf("(%.1f%%)\n", 100 * mean(d$inc)))
	tcol <- grep(paste0("^", inc_col, ".*(time|fu|years|days)"), colnames(df), value = TRUE, ignore.case = TRUE)
	if (length(tcol)) { cat("Follow-up (", tcol[1], "): median / IQR\n"); print(quantile(d[[tcol[1]]], c(.25, .5, .75), na.rm = TRUE)) }

	or_row <- function(m, term) {
		ci <- confint.default(m)[term, ]
		sprintf("beta %.3f | OR %.2f (%.2f, %.2f) | p %.2e", coef(m)[term], exp(coef(m)[term]),
			exp(ci[1]), exp(ci[2]), coef(summary(m))[term, 4])
	}

	cat("\nRaw PRS scale (compare to manuscript 0.417 -> 0.409):\n")
	m0r <- glm(inc ~ cad_prs + age + sex, family = binomial, data = d)
	m1r <- glm(inc ~ cad_prs + L1CAM + age + sex, family = binomial, data = d)
	cat("  base:     ", or_row(m0r, "cad_prs"), "\n")
	cat("  + L1CAM:  ", or_row(m1r, "cad_prs"), "\n")

	for (adj in c("age + sex", paste("age + sex +", paste(pcs, collapse = " + ")))) {
		cat("\nPer SD of CAD PRS, adjusted for:", ifelse(grepl("PC1", adj), "age, sex, PC1-10", "age, sex"), "\n")
		f0 <- as.formula(paste("inc ~ cad_prs_std +", adj))
		f1 <- as.formula(paste("inc ~ cad_prs_std + L1CAM +", adj))
		fi <- as.formula(paste("inc ~ cad_prs_std * L1CAM +", adj))
		fa <- as.formula(paste("L1CAM ~ cad_prs_std +", adj))
		m0 <- glm(f0, family = binomial, data = d); m1 <- glm(f1, family = binomial, data = d)
		cat("  PRS, base model:      ", or_row(m0, "cad_prs_std"), "\n")
		cat("  PRS, + L1CAM:         ", or_row(m1, "cad_prs_std"), "\n")
		cat("  L1CAM (per SD):       ", or_row(m1, "L1CAM"), "\n")
		cat("  PRS x L1CAM:          ", or_row(glm(fi, family = binomial, data = d), "cad_prs_std:L1CAM"), "\n")
		a <- lm(fa, data = d)
		cat(sprintf("  a-path, L1CAM ~ PRS:   beta %.3f, p %.2e\n", coef(a)["cad_prs_std"], coef(summary(a))["cad_prs_std", 4]))
		pm <- function(dd) { b0 <- coef(glm(f0, family = binomial, data = dd))["cad_prs_std"]
			b1 <- coef(glm(f1, family = binomial, data = dd))["cad_prs_std"]; (b0 - b1) / b0 }
		est <- pm(d)
		boots <- replicate(n_boot, pm(d[sample(nrow(d), replace = TRUE), ]))
		cat(sprintf("  Proportion mediated (difference method): %.1f%% (bootstrap 95%% CI %.1f%%, %.1f%%)\n",
			100 * est, 100 * quantile(boots, .025), 100 * quantile(boots, .975)))
	}

	cat("\n4-group analysis (median splits within participants with diabetes)\n")
	d$prs_grp <- ifelse(d$cad_prs_std > median(d$cad_prs_std), "HighPRS", "LowPRS")
	d$l1_grp <- ifelse(d$L1CAM > median(d$L1CAM), "HighL1CAM", "LowL1CAM")
	d$grp <- relevel(factor(paste(d$prs_grp, d$l1_grp, sep = "_")), ref = "LowPRS_HighL1CAM")
	print(d[, .(n = .N, events = sum(inc), pct = round(100 * mean(inc), 1)), by = grp][order(grp)])
	m4 <- glm(inc ~ grp + age + sex, family = binomial, data = d)
	for (g in grep("^grp", names(coef(m4)), value = TRUE)) cat(" ", sub("^grp", "", g), "vs LowPRS_HighL1CAM:", or_row(m4, g), "\n")
}

sink()

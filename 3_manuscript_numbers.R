# Manuscript numbers not produced by other scripts:
#   1 - Table 1 (all, diabetes, no diabetes)
#   2 - CAD PRS by diabetes status

library(data.table)

indir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/input/"
outdir <- "/medpop/esp/mkaminen/ukb_proteomics_cvd/output/"

df <- fread(paste0(indir, "ukb_proteomics_baseline_excl_and_imput_noprevcvd.tsv.gz"))
df$dm <- ifelse(df$a1c >= 48 | df$dm2_prev == 1, 1, 0)
df <- df[!is.na(df$dm), ]
df$cad_prs_std <- as.numeric(scale(df$cad_prs))

### 1 - Table 1 ####
print(table(Sex_numeric = df$Sex_numeric))
print(table(mergedrace = df$mergedrace))

df$female <- as.integer(df$Sex_numeric == 0)
df$white <- as.integer(df$mergedrace == names(which.max(table(df$mergedrace))))  # most common category = White; confirm above
df$htn <- as.integer(df$SBP >= 140 | df$DBP >= 90 | df$antihtnbase == 1)
df$hld <- as.integer(df$cholmed == 1 | df$ldl >= 4.9)  # LDL >= 190 mg/dL

cont <- c(Age = "age", BMI = "BMI", SBP = "SBP", DBP = "DBP", TotalChol = "tchol",
	LDL = "ldl", HDL = "hdl", TG = "tg")
bin <- c(Female = "female", White = "white", EverSmoked = "ever_smoked",
	Hypertension = "htn", Hyperlipidemia = "hld", Diabetes = "dm")

summ <- function(d) {
	out <- c(N = format(nrow(d), big.mark = ","))
	for (v in names(cont)) {
		x <- as.numeric(d[[cont[v]]])
		out[v] <- sprintf("%.1f (%.1f)", mean(x, na.rm = TRUE), sd(x, na.rm = TRUE))
	}
	for (v in names(bin)) {
		x <- as.numeric(as.character(d[[bin[v]]]))
		out[v] <- sprintf("%s (%.1f%%)", format(sum(x == 1, na.rm = TRUE), big.mark = ","), 100 * mean(x == 1, na.rm = TRUE))
	}
	out
}

tab <- data.frame(Characteristic = c("N", names(cont), names(bin)),
	All = summ(df), Diabetes = summ(df[df$dm == 1, ]), NoDiabetes = summ(df[df$dm == 0, ]))
print(tab, right = FALSE, row.names = FALSE)
write.csv(tab, paste0(outdir, "table1_ukb.csv"), row.names = FALSE)

### 2 - CAD PRS by diabetes status ####
print(df[, .(n = .N, mean_prs_std = mean(cad_prs_std), sd_prs_std = sd(cad_prs_std)), by = dm])
fit <- lm(as.formula(paste("cad_prs_std ~ dm + age + factor(Sex_numeric) +", paste0("PC", 1:10, collapse = " + "))), data = df)
print(cbind(coef(summary(fit))["dm", , drop = FALSE], confint(fit)["dm", , drop = FALSE]))

rm(list = ls())

library(openxlsx)
library(impute)
library(lme4)
library(lumi)
library(biomaRt)
library(org.Hs.eg.db)
library(lumiHumanIDMapping)


### Loading the original data

load("Data/Original/DNAmethylation_Geneexprs_Covariates_Lung_cancer_data_set.RData")
covars <- readRDS("Data/Original/full_covariates_nowac.rds")
hm450 <- readRDS("Data/Original/hm450.rds")


### Methylation

# Extracting smoking-related CpG sites
london_nc <- read.xlsx("Data/Meta_analyses/001506_supplemental_tables.xlsx", sheet = 3, startRow = 3)
london <- read.xlsx("Data/Meta_analyses/001506_supplemental_tables.xlsx", sheet = 6, startRow = 3)
london_bonf <- london[which(london$`P-value` < (0.05 / nrow(beta.LC))), ]
ids <- intersect(london_bonf$Probe.ID, rownames(beta.LC))
dnam <- beta.LC[ids, ]
dnam <- t(dnam)

# Removing CpG sites with more than 30% of missing data
dnam <- dnam[, apply(dnam, 2, FUN = function(x) {
  sum(is.na(x)) / length(x)
}) < 0.3]

# Logit2-transformation
M <- log2(dnam / (1 - dnam))

# Imputation
M_imputed <- impute.knn(M)$data

# Linkage with covariate data
LCsamples$Row.names.1 <- paste0(LCsamples$chip, "_", LCsamples$chip.pos)
rownames(M_imputed) <- LCsamples[rownames(M_imputed), "Row.names.1"]
covars <- covars[!is.na(covars$sampleID), ]
rownames(covars) <- covars$sampleID
ids <- intersect(rownames(M_imputed), rownames(covars))
M_imputed <- M_imputed[ids, ]
covars <- covars[ids, ]

# Denoising on technical effects and lung cancer
M_denoised <- matrix(NA, nrow = nrow(M_imputed), ncol = ncol(M_imputed))
rownames(M_denoised) <- rownames(M_imputed)
colnames(M_denoised) <- colnames(M_imputed)
pvalues <- rep(NA, ncol(M_imputed))
for (k in 1:ncol(M_imputed)) {
  print(k)
  y <- covars$caco
  z <- as.character(factor(covars$smoking.status, levels = c("Never", "Former", "Current"), labels = c(0, 1, 2)))
  z1 <- ifelse(z == 1, yes = 1, no = 0)
  z2 <- ifelse(z == 2, yes = 1, no = 0)
  model0 <- lmer(M_imputed[, k] ~ (1 | covars$chip) + (1 | covars$chip.pos) + z1 + z2)
  model <- lmer(M_imputed[, k] ~ (1 | covars$chip) + (1 | covars$chip.pos) + y + z1 + z2)
  pvalues[k] <- anova(model0, model)$`Pr(>Chisq)`[2]
  M_denoised[, k] <- residuals(model) + fixef(model)["z1"] * z1 + fixef(model)["z2"] * z2
}
M_denoised <- scale(M_denoised)

# Annotation
annot <- hm450[colnames(M_denoised), ]

# Saving prepared dataset and annotation file
saveRDS(M_denoised, paste0("Data/NOWAC_MTT_smoking_", ncol(M_denoised), ".rds"))
saveRDS(annot, paste0("Data/NOWAC_MTT_annot_", ncol(M_denoised), ".rds"))


### Gene expression

# Annotation
nuID <- rownames(exprs)
annot <- getNuIDMappingInfo(nuID, lib.mapping = "lumiHumanIDMapping")
x <- org.Hs.egCHR
chr_genes <- mappedkeys(x)
annot <- cbind(annot, chr = unlist(as.list(x[chr_genes]))[annot[, 2]])

# Extracting smoking-related transcripts
huan <- read.xlsx("Data/Meta_analyses/HMG-2016-JH-00187_Huan_Fig1SupplTables.xlsx", sheet = 2, startRow = 2)
entrez_smoking <- huan$Entrez.Gene.ID[which((huan$Meta.P.Value < 0.05 / 20000))] # number of genes
entrez_smoking <- entrez_smoking[entrez_smoking %in% annot[, 2]]
annot <- annot[annot[, 2] %in% entrez_smoking, ]
tr <- exprs[rownames(annot), ]
tr <- t(tr)
tr <- log(tr)
tr <- scale(tr)

# Renaming the rows
rownames(covars) <- paste(covars$case_ctrl, covars$pair.no)
ids <- intersect(rownames(covars), rownames(tr))
tr <- tr[ids, ]
covars <- covars[ids, ]
rownames(tr) <- covars$sampleID

# Denoising on lung cancer
T_denoised <- matrix(NA, nrow = nrow(tr), ncol = ncol(tr))
rownames(T_denoised) <- rownames(tr)
colnames(T_denoised) <- colnames(tr)
pvalues <- rep(NA, ncol(tr))
for (k in 1:ncol(tr)) {
  print(k)
  y <- covars$caco
  z <- as.character(factor(covars$smoking.status, levels = c("Never", "Former", "Current"), labels = c(0, 1, 2)))
  z1 <- ifelse(z == 1, yes = 1, no = 0)
  z2 <- ifelse(z == 2, yes = 1, no = 0)
  model <- lm(tr[, k] ~ y + z1 + z2)
  pvalues[k] <- summary(model)$coefficients["y", 4]
  T_denoised[, k] <- residuals(model) + coef(model)["z1"] * z1 + coef(model)["z2"] * z2
}
T_denoised <- scale(T_denoised)

# Saving prepared dataset and annotation file
saveRDS(tr, paste0("Data/NOWAC_TTX_smoking_", ncol(tr), ".rds"))
saveRDS(annot, paste0("Data/NOWAC_TTX_annot_", ncol(tr), ".rds"))

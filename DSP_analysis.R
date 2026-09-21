# GeoMx DSP spatial proteomics analysis 
# Hip versus knee osteoarthritis synovium
library(readxl)
library(standR)
library(SpatialExperiment)
library(SummarizedExperiment)
library(edgeR)
library(limma)
library(scater)
library(dplyr)
library(tibble)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(openxlsx)


# ==========================================================
# 1. Import data

counts <- read_xlsx("data/counts.xlsx") |> as.data.frame()
featureanno <- read_xlsx("data/featureanno.xlsx") |> as.data.frame()
sampleanno <- read_xlsx("data/sampleanno.xlsx") |> as.data.frame()

spe0 <- readGeoMx(
  countFile = counts,
  sampleAnnoFile = sampleanno,
  featureAnnoFile = featureanno,
  colnames.as.rownames = c("TargetName"),
  coord.colnames = c("ROICoordinateX", "ROICoordinateY"),
  rmNegProbe = FALSE
)

# ==========================================================
# 2. Protein and ROI quality control

spe <- addPerROIQC(spe0, rm_genes = TRUE)

colData(spe)$part <- factor(colData(spe)$part)
colData(spe)$region <- factor(colData(spe)$region)
colData(spe)$CD45 <- factor(colData(spe)$CD45)

plotROIQC(
  spe,
  x_axis = "AOINucleiCount", y_axis = "lib_size",
  x_lab = "Nuclei count", y_lab = "Library size",
  x_threshold = 50, y_threshold = 10000,
  col = part,
  threshold_col = "red",
  threshold_linetype = "dashed"
)

keep_roi <- colData(spe)$AOINucleiCount >= 50 &
  colData(spe)$lib_size >= 10000

spe <- spe[, keep_roi]
spe_raw <- spe

# ==========================================================
# 3. Candidate normalization factor assessment

candidate_controls <- c("GAPDH", "Calreticulin", "RPS6", "Histone H3", "TOMM20")

control_expr <- assay(spe, "counts")[candidate_controls, , drop = FALSE]
control_expr <- log2(control_expr + 1)

control_cor <- cor(
  t(control_expr),
  use = "pairwise.complete.obs",
  method = "pearson"
)

pheatmap(
  control_cor,
  display_numbers = round(control_cor, 2),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  color = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
  main = "Correlation of candidate normalization factors"
)

# GAPDH and Calreticulin showed poor concordance and were excluded.
NCGs <- c("RPS6", "Histone H3", "TOMM20")

# ==========================================================
# 4. TMM normalization and RUVg correction

spe_tmm <- geomxNorm(spe, method = "TMM")

spe_norm <- geomxBatchCorrection(
  spe_tmm,
  factors = "part",
  NCGs = NCGs,
  k = 2,
  method = "RUVg"
)

# ==========================================================
# 5. Normalization assessment

plotRLExpr(spe_raw, ordannots = "part", assay = 2, col = part) +
  ggtitle("Raw")

plotRLExpr(spe_norm, ordannots = "part", assay = 2, col = part) +
  ggtitle("Normalization")

set.seed(100)
spe_norm <- scater::runPCA(spe_norm)
pca_results <- reducedDim(spe_norm, "PCA")

drawPCA(
  spe_norm,
  precomputed = pca_results,
  col = part
)

# ==========================================================
# 6. Differential protein analysis

colData(spe_norm)$group <- interaction(
  colData(spe_norm)$part,
  colData(spe_norm)$region,
  colData(spe_norm)$CD45,
  sep = "_",
  drop = TRUE
)

design <- model.matrix(
  ~ 0 + group + ruv_W1 + ruv_W2,
  data = as.data.frame(colData(spe_norm))
)

dge <- spe2dge(spe_norm)

# Filter low-abundance protein targets
keep <- filterByExpr(dge, design = design)
dge <- dge[keep, , keep.lib.sizes = FALSE]

# limma-voom
v <- voom(dge, design, plot = TRUE)
fit <- lmFit(v, design)


# ==========================================================
# 7. Hip versus knee contrasts

contrast_matrix <- makeContrasts(
  Intima_CD45P = groupHip_Intima_CD45P - groupKnee_Intima_CD45P,
  Intima_CD45N = groupHip_Intima_CD45N - groupKnee_Intima_CD45N,
  Subintima_CD45P = groupHip_Subintima_CD45P - groupKnee_Subintima_CD45P,
  Subintima_CD45N = groupHip_Subintima_CD45N - groupKnee_Subintima_CD45N,
  levels = design
)

fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# ==========================================================
# 8. Export differential abundance results
contrast_names <- colnames(contrast_matrix)

results <- lapply(contrast_names, function(x) {
  topTable(
    fit2,
    coef = x,
    number = Inf,
    adjust.method = "BH"
  ) |>
    rownames_to_column("TargetName")
})

names(results) <- contrast_names

wb <- createWorkbook()

for (nm in names(results)) {
  addWorksheet(wb, nm)
  writeDataTable(wb, nm, results[[nm]], tableStyle = "TableStyleMedium2")
  setColWidths(wb, nm, cols = seq_len(ncol(results[[nm]])), widths = "auto")
}

saveWorkbook(
  wb,
  "Supplementary_Table_S3_DSP_DE.xlsx",
  overwrite = TRUE
)


# ==========================================================
# 9. Volcano plots

dir.create("results/volcano", recursive = TRUE, showWarnings = FALSE)

make_volcano <- function(df, contrast_name) {
  
  plot_df <- df |>
    mutate(
      neg_log10_p = -log10(P.Value),
      significance = ifelse(
        P.Value < 0.05 & abs(logFC) > 0.5,
        "Significant",
        "Not significant"
      )
    )
  
  p <- ggplot(plot_df, aes(logFC, neg_log10_p)) +
    geom_point(aes(color = significance), alpha = 0.8, size = 2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed") +
    labs(
      title = contrast_name,
      x = "log2 fold change (Hip - Knee)",
      y = "-log10(p-value)"
    ) +
    theme_classic()
  
  ggsave(
    file.path("results/volcano", paste0(contrast_name, "_volcano.pdf")),
    p, width = 6, height = 5
  )
  
  p
}

for (nm in names(results)) {
  make_volcano(results[[nm]], nm)
}


# ==========================================================
# 10. Export significant proteins for IPA

dir.create("results/IPA", recursive = TRUE, showWarnings = FALSE)

for (nm in names(results)) {
  
  ipa_input <- results[[nm]] |>
    filter(P.Value < 0.05, abs(logFC) > 0.5) |>
    select(TargetName, logFC, P.Value, adj.P.Val)
  
  write.csv(
    ipa_input,
    file.path("results/IPA", paste0(nm, "_IPA_input.csv")),
    row.names = FALSE
  )
}


# ==========================================================
# 11. Session information

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

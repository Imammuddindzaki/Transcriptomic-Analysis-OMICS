#Modul: Analisis Ekspresi Gen Respon Sepsis
#Dataset: GSE65682 (CAP vs No CAP) Pneumonia Diagnoses 
#Platform: 	Expression profiling by array ([HG-U219] Affymetrix Human Genome U219 Array - GPL13667)
#Tujuan: Mengidentifikasi Differentially Expressed Genes (DEG)
#PART A. PENGANTAR KONSEP
#Analisis ekspresi gen bertujuan untuk membandingkan tingkat ekspresi gen
#antara dua kondisi biologis (Community-Acquired Pneumonia vs No CAP)

#PERSIAPAN LINGKUNGAN KERJA (INSTALL & LOAD PACKAGE)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("GEOquery", "limma"), ask = FALSE, update = FALSE)
install.packages(c("pheatmap", "ggplot2", "dplyr"))
if (!requireNamespace("umap", quietly = TRUE)) {
  install.packages("umap")
}
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(dplyr)
library(umap)

# PENGAMBILAN DATA DARI GEO
gset <- getGEO("GSE65682", GSEMatrix = TRUE, AnnotGPL = TRUE)[[1]] 
# kalau ngga bisa download, bisa manual dengan cara
# Buka link ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65682/matrix/ dan klik filenya
# Kemudian gset <- getGEO(filename = "GSE65682_series_matrix.txt.gz")

# PRE-PROCESSING DATA EKSPRESI
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
LogTransform <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (LogTransform) {
  ex[ex <= 0] <- NA
  ex <- log2(ex)
}

# DEFINISI KELOMPOK SAMPEL
#Kita mengambil kolom ‘Pneumonia diagnoses’ yang berisi info CAP vs No CAP
group_info <- pData(gset)[["pneumonia diagnoses:ch1"]]
groups <- make.names(group_info)
gset$group <- factor(groups)
nama_grup <- levels(gset$group)
print(nama_grup)

#DESIGN MATRIX (KERANGKA STATISTIK)
design <- model.matrix(~0 + gset$group)
colnames(design) <- levels(gset$group)
grup_CAP <- "cap"
grup_No_CAP <- "no.cap"
contrast_formula <- paste(grup_CAP, "-", grup_No_CAP)
print(paste("Kontras yang dianalisis:", contrast_formula))

# ANALISIS DIFFERENTIAL EXPRESSION (LIMMA)
fit <- lmFit(ex, design)
contrast_matrix <- makeContrasts(contrasts = contrast_formula, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
topTableResults <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "B",
  number = Inf,
  p.value = 0.05
)
head(topTableResults)

# 1 BOXPLOT DISTRIBUSI NILAI EKSPRESI
group_colors <- as.numeric(gset$group)
boxplot(
  ex,
  col = group_colors,
  las = 2,
  outline = FALSE,
  main = "Boxplot Distribusi Nilai Ekspresi per Sampel",
  ylab = "Expression Value (log2)"
)
legend(
  "topright",
  legend = levels(gset$group),
  fill = unique(group_colors),
  cex = 0.8
)

# 2 DISTRIBUSI NILAI EKSPRESI (DENSITY PLOT)
expr_long <- data.frame(
  Expression = as.vector(ex),
  Group = rep(gset$group, each = nrow(ex))
)
ggplot(expr_long, aes(x = Expression, color = Group)) +
  geom_density(linewidth = 1) +
  theme_minimal() +
  labs(
    title = "Distribusi Nilai Ekspresi Gen (CAP & No_CAP)",
    x = "Expression Value (log2)",
    y = "Density"
  )

# 3 UMAP (VISUALISASI DIMENSI RENDAH)
umap_input <- t(ex)
umap_result <- umap(umap_input)
umap_df <- data.frame(
  UMAP1 = umap_result$layout[, 1],
  UMAP2 = umap_result$layout[, 2],
  Group = gset$group
)
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  theme_minimal() +
  labs(
    title = "UMAP Plot: CAP vs No_CAP",
    x = "UMAP 1",
    y = "UMAP 2"
  )

# VISUALISASI VOLCANO PLOT
topTableResults <- cbind(
  PROBEID = rownames(topTableResults),
  topTableResults
  )
length(topTableResults$PROBEID)
idx <- match(topTableResults$PROBEID, rownames(fData(gset)))
topTableResults$SYMBOL <- fData(gset)[idx, "Gene Symbol"]
topTableResults$ENTREZID <- fData(gset)[idx, "Entrez Gene"]
volcano_data <- data.frame(
  logFC = topTableResults$logFC,
  adj.P.Val = topTableResults$adj.P.Val,
  Gene = topTableResults$SYMBOL
)
volcano_data$status <- "NO"
volcano_data$status[volcano_data$logFC > 1 & volcano_data$adj.P.Val < 0.05] <- "UP"
volcano_data$status[volcano_data$logFC < -1 & volcano_data$adj.P.Val < 0.05] <- "DOWN"
ggplot(volcano_data, aes(x = logFC, y = -log10(adj.P.Val), color = status)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_minimal() +
  ggtitle("Volcano Plot DEG CAP")

# VISUALISASI HEATMAP
topTableResults <- topTableResults[
  order(topTableResults$adj.P.Val),
]
top50 <- head(topTableResults, 50)
mat_heatmap <- ex[top50$PROBEID, ]
gene_label <- ifelse(
  is.na(top50$SYMBOL) | top50$SYMBOL == "",
  top50$PROBEID, # jika SYMBOL kosong → probe ID
  top50$SYMBOL # jika ada → gene symbol
)
rownames(mat_heatmap) <- gene_label
mat_heatmap <- mat_heatmap[
  rowSums(is.na(mat_heatmap)) == 0,
]
gene_variance <- apply(mat_heatmap, 1, var)
mat_heatmap <- mat_heatmap[gene_variance > 0, ]
annotation_col <- data.frame(
  Group = gset$group
)
rownames(annotation_col) <- colnames(mat_heatmap)
pheatmap(
  mat_heatmap,
  scale = "row", 
  annotation_col = annotation_col,
  show_colnames = FALSE, 
  show_rownames = TRUE,
  fontsize_row = 7,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  main = "Top 50 Differentially Expressed Genes"
)


# Kita mencoba ENRICHMENT disini
# Kita melihat gen apa saja yang UP dan DOWN sehingga dapat kita kaitkan pathway terkait
# Melihat berapa banyak yang DOWN dan UP kemudian meng-groupkannya
table(volcano_data$status)
up_genes <- topTableResults[
  topTableResults$logFC >= 1 & topTableResults$adj.P.Val < 0.05,
]
down_genes <- topTableResults[
  topTableResults$logFC <= -1 & topTableResults$adj.P.Val < 0.05,
]

# Kita cek kolom ENTREZID dan pastikan bersih, pastikan tidak ada "1234 /// 5678"
head(up_genes$ENTREZID)
up_entrez <- unique(na.omit(up_genes$ENTREZID))
down_entrez <- unique(na.omit(down_genes$ENTREZID))
length(up_entrez)
length(down_entrez)

# Kita akan pakai clusterProfiler.
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
library(clusterProfiler)
library(org.Hs.eg.db)
ego_up <- enrichGO(
  gene          = up_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Kita baca yang "UP" dulu
dotplot(ego_up, showCategory = 10) +
  ggtitle("GO Enrichment - UP genes (CAP)")
# Bisa lihat TOP 20 pakai ini
head(as.data.frame(ego_up), 20)

# Kemudian baca yang "DOWN"
ego_down <- enrichGO(
  gene          = down_entrez,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "BP",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05,
  readable      = TRUE
)
dotplot(ego_down, showCategory = 10) +
  ggtitle("GO Enrichment - DOWN genes (CAP)")
head(as.data.frame(ego_down), 20)

# KEGG untuk mempertajam, Kali ini untuk UP Genes
kegg_up <- enrichKEGG(
  gene         = up_entrez,
  organism     = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)
dotplot(kegg_up, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment - UP genes (CAP)")

# KEGG untuk Down Genes
kegg_down <- enrichKEGG(
  gene         = down_entrez,
  organism     = "hsa",
  pAdjustMethod = "BH",
  qvalueCutoff  = 0.05
)

dotplot(kegg_down, showCategory = 10) +
  ggtitle("KEGG Pathway Enrichment - DOWN genes (CAP)")

# Metode GSEA dapat dilakukan jika hasil tidak muncul
# Pertama, kita siapkan Gene list
all_results <- topTable(
  fit2,
  adjust = "fdr",
  sort.by = "none",
  number = Inf
)
all_results <- cbind(
  PROBEID = rownames(all_results),
  all_results
)
idx_all <- match(all_results$PROBEID, rownames(fData(gset)))
all_results$ENTREZID <- fData(gset)[idx_all, "Entrez Gene"]
all_results <- all_results[!is.na(all_results$ENTREZID), ]
all_results <- all_results[order(abs(all_results$logFC), decreasing = TRUE), ]
all_results <- all_results[!duplicated(all_results$ENTREZID), ]
gene_list <- all_results$logFC
names(gene_list) <- all_results$ENTREZID
gene_list <- sort(gene_list, decreasing = TRUE)

# Jalankan GSEA
gsea_kegg <- gseKEGG(
  geneList = gene_list,
  organism = "hsa",
  pAdjustMethod = "BH",
  verbose = FALSE
)
dotplot(gsea_kegg, showCategory = 10) +
  ggtitle("GSEA KEGG - CAP vs No_CAP")
head(as.data.frame(gsea_kegg)) 

# Kita coba interprestasi arah NES, innate vs adaptive
gsea_df <- as.data.frame(gsea_kegg)
innate_keywords <- c("Toll", "NOD", "Cytokine", "Complement", 
                     "bacterial", "LPS", "NF-kappa")
adaptive_keywords <- c("Th", "T cell", "B cell", "IgA", 
                       "antigen", "graft", "arthritis", 
                       "diabetes", "immune network")
gsea_df$category <- "Other"
gsea_df$category[grepl(paste(innate_keywords, collapse="|"), 
                       gsea_df$Description, ignore.case=TRUE)] <- "Innate"
gsea_df$category[grepl(paste(adaptive_keywords, collapse="|"), 
                       gsea_df$Description, ignore.case=TRUE)] <- "Adaptive"
aggregate(NES ~ category, data=gsea_df, mean)
table(gsea_df$NES > 0)
sig <- subset(gsea_df, p.adjust < 0.05)
table(sig$category)

# Kita periksa lebih dalam bagian "Adaptive"
adaptive_only <- subset(gsea_df, category == "Adaptive")
table(adaptive_only$NES > 0)
table(
  Direction = ifelse(adaptive_only$NES > 0, "Positive (CAP)", "Negative (CAP)")
)
mean(subset(gsea_df, category == "Adaptive" & p.adjust < 0.05)$NES)

# Bar Perbandingan Rerata NES: Adaptive vs Innate Pathways
sig_df <- subset(gsea_df, p.adjust < 0.05)
mean_nes <- aggregate(NES ~ category, data = sig_df, mean)
mean_nes
library(ggplot2)
ggplot(mean_nes, aes(x = category, y = NES, fill = category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme_minimal(base_size = 14) +
  scale_fill_manual(values = c(
    "Adaptive" = "#2C7BB6",
    "Innate" = "#D7191C",
    "Other" = "gray70"
  )) +
  labs(
    title = "Mean Normalized Enrichment Score (NES) by Immune Category",
    x = "",
    y = "Mean NES (Significant Pathways)"
  ) +
  theme(legend.position = "none")

# Visualisasi GSEA Enrichment Plot (Adaptive Pathway)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("enrichplot")
gsea_df[gsea_df$Description == "Th1 and Th2 cell differentiation", ]
library(enrichplot)
gseaplot2(
  gsea_kegg,
  geneSetID = "hsa04658", # Th1/Th2 differentiation
  title = "GSEA Enrichment: Th1 and Th2 Cell Differentiation"
)

# Tambahan Density Plot NES Distribution
ggplot(sig_df, aes(x = NES, fill = category)) +
  geom_density(alpha = 0.5) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Distribution of NES Across Immune Categories",
    x = "Normalized Enrichment Score",
    y = "Density"
  )

# MENYIMPAN HASIL
write.csv(topTableResults, "Hasil_GSE65682_DEG.csv")
message("Analisis selesai. File hasil telah disimpan.")

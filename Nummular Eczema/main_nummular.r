# Load required libraries
library(tidyverse)
library(DESeq2)
library(edgeR)
library(EnhancedVolcano)
library(VennDiagram)
library(grid)


data <- read.csv("C:\\Users\\pksha\\Desktop\\autoimmune\\Nummular Eczema\\counts.csv", header = TRUE, row.names = 1)
info <- read.csv("C:\\Users\\pksha\\Desktop\\autoimmune\\Nummular Eczema\\meta.csv", header = TRUE, sep = ",")
meta <- info
count <- data

# Remove rows with NA values from count matrix
count <- na.omit(count)

# Remove samples with NA values from metadata
meta <- meta %>% filter(!is.na(Sample_id))

# Ensure 'Sample_id' in metadata matches columns in count data
meta <- meta %>% filter(Sample_id %in% colnames(count)) %>% arrange(match(Sample_id, colnames(count)))

# Set row names of meta to Sample_id for consistency
rownames(meta) <- meta$Sample_id

# Ensure count data and metadata dimensions match
if (ncol(count) != nrow(meta)) {
  stop("Number of columns in count matrix does not match number of rows in metadata after removing NAs.")
}

# Ensure 'Site' column in metadata is a factor
meta$Label <- as.factor(meta$Label)

# Renaming the levels of the Label column
meta$Label <- gsub(" ", "_", meta$Label)

# Perform DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = ~Label)
dds <- DESeq(dds)
res_de <- results(dds)
de_DEGs <- rownames(res_de[which(res_de$padj < 0.05 & abs(res_de$log2FoldChange) > 1), ])

# Perform edgeR analysis
group <- meta$Label
dge <- DGEList(counts = count, group = group)
dge <- calcNormFactors(dge)
dge <- estimateDisp(dge)
fit <- glmFit(dge)
lrt <- glmLRT(fit)
res_edger <- topTags(lrt, n = nrow(dge))
edger_DEGs <- rownames(res_edger[which(res_edger$table$PValue < 0.05 & abs(res_edger$table$logFC) > 1), ])

# Create Venn diagram
venn.plot <- venn.diagram(
  x = list(DESeq2 = de_DEGs, edgeR = edger_DEGs),
  category.names = c("DESeq2", "edgeR"),
  fill = c("red", "blue"),
  alpha = 0.5,
  cex = 2,
  cat.cex = 2,
  cat.pos = 0,
  cat.dist = 0.05,
  cat.fontface = "bold",
  main = "Venn Diagram of DEGs: DESeq2 vs edgeR",
  filename = NULL # Prevents the function from saving the plot to a file automatically
)

# Save the Venn diagram to a file
venn.plot.file <- "C:\\Users\\pksha\\Desktop\\autoimmune\\Nummular Eczema\\venn_diagram_DEGs.png"
png(filename = venn.plot.file)
grid.draw(venn.plot)

# Add legend to the Venn diagram
legend_grob <- legendGrob(
  labels = c("DESeq2", "edgeR"),
  pch = 15,
  gp = gpar(col = c("red", "blue"), fontsize = 12, fontface = "bold")
)
pushViewport(viewport(x = 0.85, y = 0.85, just = c("right", "top")))
grid.draw(legend_grob)
popViewport()

dev.off()

# Create volcano plot for DESeq2 results
pdf("C:\\Users\\pksha\\Desktop\\autoimmune\\Nummular Eczema\\volcano_DESeq2.pdf", width = 8, height = 6)
EnhancedVolcano(res_de,
  lab = rownames(res_de),
  x = "log2FoldChange",
  y = "pvalue",
  title = "Volcano Plot: DESeq2",
  pCutoff = 0.05,
  FCcutoff = 1,
  labSize = 3.0,
  pointSize = 3.0,
  colAlpha = 1,
  legendPosition = "right"
)
dev.off()

# Create volcano plot for edgeR results
pdf("C:\\Users\\pksha\\Desktop\\autoimmune\\Nummular Eczema\\volcano_edgeR.pdf", width = 8, height = 6)
EnhancedVolcano(res_edger$table,
  lab = rownames(res_edger$table),
  x = "logFC",
  y = "PValue",
  title = "Volcano Plot: edgeR",
  pCutoff = 0.05,
  FCcutoff = 1,
  labSize = 3.0,
  pointSize = 3.0,
  colAlpha = 1,
  legendPosition = "right"
)
dev.off()


# prerequisits
library(readxl)
library(tidyverse)
library(pheatmap)

# import data
RNA <- read_excel("data/GFPvIRF4.xlsx")

# Step 1: Generate a set of genes, that includes only the top 30 genes by padj
# at the same time making sure that padj is below < 0.01
RNA_hm <- RNA %>%
  filter(padj < 0.01) %>%
  slice_min(padj, n = 50, with_ties = FALSE) %>%
  arrange(log2FoldChange)

# Step 2: Create a row-names vector for plotting
Row_names <- as_vector(RNA_hm$geneSymbol)

# Step 3: rename columns to more sensical names
RNA_hm_2 <- RNA_hm %>%
  select(X_alignedData.Sample25_lisa_T11_GFP_r1_Aligned.out.bam,
         X_alignedData.Sample26_lisa_T11_GFP_r2_Aligned.out.bam,
         X_alignedData.Sample27_lisa_T11_IRF4_r1_Aligned.out.bam,
         X_alignedData.Sample28_lisa_T11_IRF4_r2_Aligned.out.bam
  ) %>%
  rename(
    GFP1 = X_alignedData.Sample25_lisa_T11_GFP_r1_Aligned.out.bam,
    GFP2 = X_alignedData.Sample26_lisa_T11_GFP_r2_Aligned.out.bam,
    IRF41 = X_alignedData.Sample27_lisa_T11_IRF4_r1_Aligned.out.bam,
    IRF42 = X_alignedData.Sample28_lisa_T11_IRF4_r2_Aligned.out.bam
  )

# Step 4: plot using pheatmap, euclidian clustering normalization to rows
dat <- as.matrix(RNA_hm_2)
rownames(dat) <- Row_names
pheatmap(dat, cutree_rows = 5,
         scale = "row",
         cluster_cols = FALSE,
         fontsize_row = 6
)


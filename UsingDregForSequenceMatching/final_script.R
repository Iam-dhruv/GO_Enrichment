# 1. Load the libraries (Mamba already installed them!)
library(clusterProfiler)
library(org.Hs.eg.db)

# 2. Read your specific file
genes = read.table("10genes_with_nrf1.txt", sep = "\t", header = F, stringsAsFactors = F)$V1

# 3. Run the Gene Ontology Enrichment
ego <- enrichGO(gene          = genes,
                OrgDb         = org.Hs.eg.db,
                keyType       = "SYMBOL",
                ont           = "BP",
                pAdjustMethod = "BH",
                qvalueCutoff  = 0.01)

# 4. Save the results
write.csv(as.data.frame(ego), "11GO_annotation_results.csv")

# 5. Generate the plot
pdf("12nrf1_dotplot.pdf", height = 8, width = 8)
dotplot(ego, showCategory = 20, font.size = 6)
dev.off()
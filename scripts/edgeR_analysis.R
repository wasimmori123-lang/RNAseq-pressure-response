library(edgeR)
library(ggplot2)

# ----------------------------
# 1. Load data
# ----------------------------

counts <- read.table("../data/tagcount.out",
                     header=TRUE,
                     row.names=1,
                     sep="\t",
                     quote="")

counts <- as.matrix(counts)

# ----------------------------
# 2. Define groups
# ----------------------------

group <- factor(c("LP","LP","HP","HP"))

# ----------------------------
# 3. Create DGEList
# ----------------------------

d <- DGEList(counts=counts, group=group)

# ----------------------------
# 4. Filter low-expression genes
# ----------------------------

keep <- filterByExpr(d)

d <- d[keep,,keep.lib.sizes=FALSE]

# ----------------------------
# 5. Normalize
# ----------------------------

d <- calcNormFactors(d)

# ----------------------------
# 6. Design matrix
# ----------------------------

design <- model.matrix(~group)

# ----------------------------
# 7. Estimate dispersion
# ----------------------------

d <- estimateDisp(d, design)

# ----------------------------
# 8. Fit model
# ----------------------------

fit <- glmQLFit(d, design)

# ----------------------------
# 9. Differential expression test
# ----------------------------

out <- glmQLFTest(fit, coef=2)

# ----------------------------
# 10. Extract p-values
# ----------------------------

p.value <- out$table$PValue

q.value <- p.adjust(p.value,
                    method="BH")

# ----------------------------
# 11. DEG filtering
# ----------------------------

deg <- which(q.value < 0.05 &
             abs(out$table$logFC) > 1.584)

# Upregulated genes
up <- which(q.value < 0.05 &
            out$table$logFC > 1.584)

# Downregulated genes
down <- which(q.value < 0.05 &
              out$table$logFC < -1.584)

# ----------------------------
# 12. Summary
# ----------------------------

cat("Total DEG:", length(deg), "\n")

cat("Upregulated:", length(up), "\n")

cat("Downregulated:", length(down), "\n")

# ----------------------------
# 13. Save DEG table
# ----------------------------

deg_table <- out$table[deg, ]

write.table(deg_table,
            "../results/DEG_filtered.txt",
            sep="\t",
            quote=FALSE)

# ----------------------------
# 14. Volcano plot
# ----------------------------

volcano <- out$table

volcano$qvalue <- q.value

volcano$negLogFDR <- -log10(volcano$qvalue)

volcano$group <- "Not Significant"

volcano$group[
  volcano$qvalue < 0.05 &
  volcano$logFC > 1.584
] <- "Upregulated"

volcano$group[
  volcano$qvalue < 0.05 &
  volcano$logFC < -1.584
] <- "Downregulated"

p <- ggplot(volcano,
            aes(x=logFC,
                y=negLogFDR,
                color=group)) +

  geom_point(size=2,
             alpha=0.7) +

  scale_color_manual(values=c(
    "Upregulated"="#d73027",
    "Downregulated"="#4575b4",
    "Not Significant"="grey70"
  )) +

  geom_vline(xintercept=c(-1.584,1.584),
             linetype="dashed") +

  geom_hline(yintercept=-log10(0.05),
             linetype="dashed") +

  labs(
    x="log2 Fold Change",
    y="-log10(FDR)"
  ) +

  theme_classic(base_size=14) +

  theme(
    legend.title=element_blank(),
    legend.position="top"
  )

ggsave("../figures/Volcano_plot.png",
       p,
       width=6,
       height=5,
       dpi=300)

# ----------------------------
# 15. MA plot
# ----------------------------

png("../figures/tagcount_FDR.png",
    width=800,
    height=600)

plotSmear(d,
          de.tags=rownames(deg_table))

dev.off()

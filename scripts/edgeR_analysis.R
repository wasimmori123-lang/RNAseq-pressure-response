library(edgeR)

# Load count matrix
counts <- read.table("../data/tagcount.out",
                     header=TRUE,
                     row.names=1,
                     sep="\t",
                     quote="")

counts <- as.matrix(counts)

# Define groups
group <- factor(c("LP","LP","HP","HP"))

# Create DGEList
d <- DGEList(counts=counts, group=group)

# Filter low-expression genes
keep <- filterByExpr(d)
d <- d[keep,,keep.lib.sizes=FALSE]

# Normalize
d <- calcNormFactors(d)

# Design matrix
design <- model.matrix(~group)

# Estimate dispersion
d <- estimateDisp(d, design)

# Fit model
fit <- glmQLFit(d, design)

# Differential expression test
out <- glmQLFTest(fit, coef=2)

# Extract p-values
p.value <- out$table$PValue
q.value <- p.adjust(p.value, method="BH")

# DEG criteria
deg <- which(q.value < 0.05 & abs(out$table$logFC) > 1.584)

# Upregulated genes
up <- which(q.value < 0.05 & out$table$logFC > 1.584)

# Downregulated genes
down <- which(q.value < 0.05 & out$table$logFC < -1.584)

# Print summary
cat("Total DEG:", length(deg), "\n")
cat("Upregulated:", length(up), "\n")
cat("Downregulated:", length(down), "\n")

# Save DEG table
deg_table <- out$table[deg, ]

write.table(deg_table,
            "../results/DEG_filtered.txt",
            sep="\t",
            quote=FALSE)

# MA plot
png("../figures/MA_plot.png",
    width=800,
    height=600)

plotSmear(d, de.tags=rownames(deg_table))

dev.off()

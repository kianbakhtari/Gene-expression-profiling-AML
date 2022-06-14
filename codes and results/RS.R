## STUDENT INFO: Kian Bakhtari - 97110025

## Set working directory
setwd("/Users/kian/Desktop/bio-project/")

## Imports
library(GEOquery)
library(limma)
library(umap)
library(ggplot2)
library(reshape2)
library(plyr)
library(pheatmap)
library(RColorBrewer)
library(readr)
library(Rtsne)

#### Data download and preparation
series <- "GSE48558"
platform <- "GPL6244"
gset_all <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "/Users/kian/Desktop/bio-project/")

if (length(gset_all) > 1) idx <- grep(platform, attr(gset_all, "names")) else idx <- 1
gset_all <- gset_all[[idx]]

# Make proper column names to match toptable
fvarLabels(gset_all) <- make.names(fvarLabels(gset_all))


# Group membership for all samples, with "X" being excluded, 1 being
# normal and 0 being AML patient
groups <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
                 "XXXXXXXXXXXXXXXXXX1X1XXX1X1111X1XX11XX11X1X1X1X1X1",
                 "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
                 "11111111111111111111")



#### Plot correlation heatmap for all samples (172 samples)
groups <- strsplit(groups, split="")[[1]]
ex_all <- exprs(gset_all)
pdf("initial_corr_all.pdf")
pheatmap(cor(ex_all), labels_row = groups, labels_col = groups, width = 82, height = 82,
         fontsize = 3, border_color = NA,
         main = 'correlation heatmap of all samples')
dev.off()


#### Exclude unwanted samples
groups <- paste0("0000000000000XXXXXXXXXXXXXXXXXXXXXXXXXXX1XXX1XXXXX",
                 "XXXXXXXXXXXXXXXXXX1X1XXX1X1111X1XX11XX11X1X1X1X1X1",
                 "XXX1XXX1XXXXXXXXXXXXXXXXXXXXXXXXXXXXX1111111001000",
                 "11111111111111111111")
gr <- strsplit(groups, split="")[[1]]
inc_gr <- which(gr != "X")
inc_lables = gr[inc_gr]
gset <- gset_all[ ,inc_gr]


#### Plot initial box plot to see if normalizing is required.
ex <- exprs(gset)
pdf("initial-boxplot.pdf", width = 16)
# png("initial-boxplot.png", width = 64)
boxplot(ex,
        main = "Boxplot of gene expression values",
        xlab = "Samples",
        ylab = "Gene expression values",
        col = "orange",
        border = "brown",
        )
dev.off()


#### Plot correlation heatmap
pdf("initial_corr.pdf")
pheatmap(cor(ex), labels_row = inc_lables, labels_col = inc_lables, width = 32, height = 32,
         fontsize = 5, border_color = NA,
         main = 'correlation heatmap of samples')
dev.off()


#### Dimension reduction without scaling
pc <- prcomp(ex)
pdf("inital_pca.pdf")
plot(pc, 
     main = "Variance preserved by each PC",
     xlab = "Principal components",
     # ylab = "Variance",
     col = "orange",
     )

plot(pc$x[,1:2],
     main = "Genes mapped to reduced space",
     xlab = "PC1",
     ylab = "PC2",
     )
dev.off()

#### Dimension reduction with scaling
ex.scale <- t(scale(t(ex), scale = F))
pc <- prcomp(ex.scale)
pdf("scaled-pca.pdf")
plot(pc, 
     main = "Variance preserved by each PC (after scaling)",
     xlab = "Principal components",
     # ylab = "Variance",
     col = "orange",
)

plot(pc$x[,1:2],
     main = "Genes mapped to reduced space (after scaling)",
     xlab = "PC1",
     ylab = "PC2",
)
dev.off()


#### Plot samples in PCA space
pcr <- data.frame(pc$r[,1:3], Group=inc_lables)
pdf("PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3) + theme_bw() + 
  ggtitle("Samples in reduced space") + theme(plot.title = element_text(hjust = 0.5))
dev.off()


#### Plot heatmap after PCA
pdf("pca_corr.pdf")
pheatmap(cor(t(data.matrix(pcr[,1:3]))), labels_row = inc_lables,
           labels_col = inc_lables, width = 32, height = 32,
           fontsize = 5, border_color = NA,
           main = 'correlation heatmap of samples (after PCA)')
dev.off()



#### t-SNE
ex_t <- t(ex)
ex_t.labels <- as.factor(inc_lables)
tsne <- Rtsne(ex_t[,-1], dims = 2, perplexity=20, verbose=TRUE, max_iter = 1000)

pdf("tSNE-3.pdf")
colors = rainbow(length(unique(ex_t.labels)))
names(colors) = unique(ex_t.labels)
par(mgp=c(2.5,1,0))
plot(tsne$Y, t='n', main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2",
     "cex.main"=2, "cex.lab"=1.5)
text(tsne$Y, labels=ex_t.labels, col=colors[ex_t.labels])
dev.off()


tsne2 <- Rtsne(ex_t[,-1], dims = 3, perplexity=10, verbose=TRUE, max_iter = 2000)
pdf("tsne_corr.pdf")
pheatmap(cor(t(tsne2$Y)), labels_row = inc_lables,
         labels_col = inc_lables, width = 32, height = 32,
         fontsize = 5, border_color = NA,
         main = 'correlation heatmap of samples (after t-SNE)')
dev.off()



############################################## Gene expression analysis
gs <- factor(inc_lables)
groups <- make.names(c("AML","Healthy"))
levels(gs) <- groups
gset$group <- gs
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gs)
fit <- lmFit(gset, design) 

cts <- paste(groups[1], groups[2], sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)

tT <- subset(tT, select=c("Gene.symbol","ID","adj.P.Val","logFC"))
write.table(tT, file="/Users/kian/Desktop/bio-project/results.txt",
            row.names=F, sep="\t", quote = F)



#### Pathway and gene anthology
aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <- unique(as.character(strsplit2(aml.up$Gene.symbol, "///")))
write.table(aml.up.genes, file = "/Users/kian/Desktop/bio-project/aml-up-genes.txt", 
            quote = F, row.names = F, col.names = F)

aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <- unique(as.character(strsplit2(aml.down$Gene.symbol, "///")))
write.table(aml.down.genes, file = "/Users/kian/Desktop/bio-project/aml-down-genes.txt",
            quote = F, row.names = F, col.names = F)




# First line
setwd("D:/MASTER/term_1_Problems/Intro_To_Bioinformatics/R project/")
library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(gplots)
library(ggplot2)
library(reshape2)
library(plyr)
series <- "GSE48558"
platform <- "GPL6244"

#### Load data  
gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
#length(gset)
#class(gset)
#names(gset)

#gset <- gset[[1]]
if (length(gset) > 1) idx <- grep(platform, attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
#length(gset)
#### Grouping samples

# group membership for all samples

gsms <- paste0("1111111111111XXXXXXXXXXXXXXXXXXXXXXXXXXX0XXX0XXXXX",
               "XXXXXXXXXXXXXXXXXX0X0XXX0X0000X0XX00XX00X0X0X0X0X0",
               "XXX0XXX0XXXXXXXXXXXXXXXXXXXXXXXXXXXXX0000000110111",
               "00000000000000000000")
sml <- strsplit(gsms, split="")[[1]]

# filter out excluded samples (marked as "X")
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# grouping
gr <- c(rep("AML",13),
        rep("Granulocytes",2),
        rep("B", 1),
        rep("T", 1),
        rep("Granulocytes",2),
        rep("Monocytes",2),
        rep("B", 1),
        rep("T", 5),
        rep("B", 1),
        rep("T", 1),
        rep("B", 1),
        rep("T", 1),
        rep("CD34",3),
        rep("Granulocytes",7),
        rep("AML",2),
        rep("T", 1),
        rep("AML",3),
        rep("B", 7),
        rep("T", 1),
        rep("Monocytes",4),
        rep("Granulocytes",1),
        rep("T", 7))

#### expression matrix
#### Log2 scale, if required
ex <- exprs(gset)
#dim(ex)

#### find out that data is normal or not
# ex <- log2(ex + 1) # it is not needed hear
# exprs(gset) <- ex # it is not needed hear
#max(ex)
#min(ex)

#### quality control ####
pdf("Results/boxplot.pdf", width = 64)
boxplot(ex)
dev.off() 

#### Normalize, if required
#### Normalizing
# ex <- normalizeQuantiles(ex) # it is not needed hear
# exprs(gset) <- ex # it is not needed hear

#?scale
ex.scale <- t(scale(t(ex), scale = F))

#### Principal Components Analysis
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

pc <- prcomp(ex.scale)
pdf("Results/PC_scaled.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()

# for scanning samples, samples are in rotation of pc
pcr <- data.frame(pc$rotation[,1:3], Group = gr)
pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(x = PC1, y = PC2, color = Group)) + geom_point(size = 3) + theme_bw()
dev.off()

#### Correlation Heatmap
ex.cor <- cor(ex.scale)
rownames(ex.cor) <- gr
colnames(ex.cor) <- gr 
ex.cor.gr <- ex.cor[, !colnames(ex.cor) %in% c("AML")]
ex.cor.gr <- ex.cor.gr[rownames(ex.cor) %in% c("AML"), ]
max(ex.cor.gr)

pdf("Results/CorHeatmap.pdf", width = 15,height = 15)
pheatmap(ex.cor, labels_row = gr, labels_col = gr,color = bluered(256), border_color = NA)
dev.off()

#### Differential Expression Analysis ####
gr <- factor(gr)
gset$group <-  (gr)
design <- model.matrix(~group + 0, gset)
colnames(design) <- levels(gr)

fit <- lmFit(gset, design)  # fit linear model

# set up contrasts of interest and recalculate model coefficients
cts <- paste("Monocytes", "AML", sep="-")
cont.matrix <- makeContrasts(contrasts=cts, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)

# compute statistics and table of top significant genes
fit2 <- eBayes(fit2, 0.01)
tT <- topTable(fit2, adjust="fdr", sort.by="P", number=Inf)

tT <- subset(tT, select=c("Gene.symbol","Gene.ID", "adj.P.Val", "logFC"))
write.table(tT, "Results/AML_Monocytes.txt", row.names=F, sep="\t", quote = F)

aml.up <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.up.genes <- unique(aml.up$Gene.symbol)
#aml.up.genes <- sub("///.*", "", aml.up.genes)
aml.up.genes <- unique(as.character(strsplit2(aml.up.genes, "///")))
write.table(aml.up.genes, file = "Results/Monocytes_AML_up.txt", quote = F, row.names = F, col.names = F)

aml.down <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.down.genes <- unique(aml.down$Gene.symbol)
#aml.down.genes <- sub("///.*", "", aml.down.genes)
aml.down.genes <- unique(as.character(strsplit2(aml.down.genes, "///")))
write.table(aml.down.genes, file = "Results/Monocytes_AML_down.txt", quote = F, row.names = F, col.names = F)

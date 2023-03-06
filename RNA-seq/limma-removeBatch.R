library(limma)
library(FactoMineR)
library(factoextra)
library(patchwork)

########
tpm<-read.table("for-normalized-four-cohort-tpm.txt", header = T, row.names = 1, sep = "\t")
batch <- c(rep('CAGA',48),rep('CPGEA',134),rep('TCGA',497),rep('WCDT',99))
batch <- as.factor(batch)
design <- model.matrix(~0 + batch)
########

condition <- data.frame(sample = colnames(tpm),
                        group = c(rep('CAGA',48),rep('CPGEA',134),rep('TCGA',497),rep('WCDT',99)),
                        batch = c(rep('batch1',48),rep('batch2',134),rep('batch3',497),rep('batch4',99)))

log_tpm<-log10(tpm+1)
pre.pca <- PCA(t(log_tpm),graph = FALSE)
pre<-fviz_pca_ind(pre.pca,
                  geom= c("point", "point"),
                  col.ind = condition$group,
                  addEllipses = TRUE,
                  legend.title="Group") + theme_bw() + 
  ggtitle("All gene\nPCA (before batch correction)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

########
tpm=normalizeBetweenArrays(tpm)
new_tpm <- removeBatchEffect(tpm, batch = batch)
########

#log_new_tpm<-log10(new_tpm+1)
limma.pca <- PCA(t(new_tpm),graph = FALSE)
limma<-fviz_pca_ind(limma.pca,
                     col.ind = condition$group,
                     geom = c("point", "point"),
                     addEllipses = TRUE,
                     legend.title="Group") + theme_bw() +
  ggtitle("All gene\nPCA (after batch correction)") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5)
  
write.table(new_tpm, file="removeBatch.txt", sep="\t", quote=F, row.names=T, col.names=T)

pdf("ALLgene-PCA.pdf", width = 9.5, height = 4)
pre+limma
dev.off()

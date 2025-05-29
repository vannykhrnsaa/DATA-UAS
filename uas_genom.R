if(!requireNamespace("BiocManager", quitely = TRUE)){install.packages("BiocManager")}
BiocManager::install(c("Biobase", "GEOquery"))
BiocManager::install("genefilter")
require(genefilter)
library(genefilter)
library(Biobase)
library(GEOquery)

#microarray
BiocManager::install("hgu133plus2.db")
library(hgu133plus2.db)
library(ggplot2)

##GSE
dtgeo <-getGEO('GDS3610')
##Melihat Data
dtgeo

#Convert to Eset
eset <- GDS2eSet(dtgeo, do.log2 = TRUE)
eset

#Phenotype Data
phdtgeo <- pData(eset)
head(phdtgeo)

#Expression Data
expdtgeo <- exprs(eset)
dim(expdtgeo)

#Anotasi
Meta(dtgeo)$platform
annotation(eset) <- "hgu133plus2.db"

#Filterasi
esetFilt <- nsFilter(eset)
print(esetFilt)

#Visual sebelum dan sesudah filterisasi
expdgeoFilt <- exprs(esetFilt$eset)
par(mfrow=c(1,2))
hist(expdtgeo, main='original')
hist(expdgeoFilt, main='filtered')
dim(exprs(esetFilt$eset))

#labelisasi
vargrp <- phdtgeo[,2]
table(vargrp)
group <- ifelse(vargrp == "control", 0,
                ifelse(vargrp == "nasopharyngeal carcinoma", 1, NA))

#save data
data_filter_t <-t(expdgeoFilt)
write.csv(data_filter_t, "uas_clean_data.csv", row.names = TRUE)
group_df=data.frame(SampleID=colnames(expdgeoFilt), group=group)
write.csv(group_df, "uas_label.csv", row.names = FALSE)

#DGE
#Limma
BiocManager::install("limma")
library(limma)
design <- model.matrix(~vargrp)
fit <- eBayes(lmFit(expdgeoFilt, design))
fit
#Top Genes
topResult <- topTable(fit, coef=2, number=50)
rownames(topResult)
#Extract top Genes
selected <- rownames(expdgeoFilt) %in% row.names(topResult)
expdtgeosel <- expdgeoFilt[selected,]
topgene <- t(expdtgeosel)

#download data
write.csv(topResult, "topgene50.csv", row.names = TRUE)
write.csv(selected, "topgene50extract.csv", row.names = TRUE)
write.csv(topgene, "uas__top50.csv", row.names = TRUE)
write.csv(topgene, "C:/Users/vanny/Documents/UI/KULIAH/S6/Genom/uas_topgene.csv", row.names = TRUE)
"C:\Users\vanny\Documents\UI\KULIAH\S6\Genom\UTS\topgene50extract.csv"

#Heatmap
heatmap(expdtgeosel)

#Visualisasi p-value
boxplot(group, expdgeoFilt[1,])

#Volcano
names <- rownames(expdgeoFilt)
result <- data.frame(Gene = names, Rheumatoid = NA, Osteoarth = NA,
                     t = NA, p.value = NA, diff = NA)

for (i in names) {
  Rheumatoid_vals <- expdgeoFilt[i, group == 1]
  Osteoarth_vals <- expdgeoFilt[i, group == 2]
  
  # T-test antara Rheumatoid dan Osteoarth
  ttest <- t.test(Rheumatoid_vals, Osteoarth_vals)
  
  result[result$Gene == i, ] <- c(
    i,
    mean(Rheumatoid_vals),
    mean(Osteoarth_vals),
    ttest$statistic,
    ttest$p.value,
    mean(Rheumatoid_vals) - mean(Osteoarth_vals)  # diff = Rheuma - Osteo
  )
}

result$log10_pval <- -log10(as.numeric(result$p.value))
result$diff <- as.numeric(result$diff)

# Volcano plot
ggplot(result, aes(x = diff, y = log10_pval)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  xlab("Difference (Rheumatoid - Osteoarthritis)") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot Rheumatoid vs Osteoarthritis") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red")  # Garis threshold p-value 0.

# Volcano Lima
allResults <- topTable(fit, coef=2, number=Inf)

# Tambah kolom -log10(p-value)
allResults$log10_pval <- -log10(allResults$P.Value)

# Library untuk plotting
library(ggplot2)

# Volcano plot
ggplot(allResults, aes(x = logFC, y = log10_pval)) +
  geom_point(alpha = 0.8) +
  theme_minimal() +
  xlab("log2 Fold Change (Rheumatoid vs Osteoarthritis)") +
  ylab("-log10(p-value)") +
  ggtitle("Volcano Plot - limma result") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "gray") +  # logFC threshold
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") # p-value threshold

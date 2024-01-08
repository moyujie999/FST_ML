rm(list = ls())
library(xlsx) 
library(edgeR)
library(limma)
library("readxl")
library(indicspecies)
library(igraph)
library(Hmisc)
library(sciplot)
library(reshape2)
library(ggpmisc)


##### import data #####
`ps-orig_16S` <- readRDS("D:/Rprojet/soil_V2/dataset/ps16S.RDS")
`ps-orig_ITS` <- readRDS("D:/Rprojet/soil_V2/dataset/psITS.RDS")

otu <- `ps-orig_16S`@otu_table
otu_ITS <- `ps-orig_ITS`@otu_table
otu_all <- cbind(otu, otu_ITS)
taxa <- `ps-orig_16S`@tax_table
taxa_ITS <- `ps-orig_ITS`@tax_table
design16 <- `ps-orig_16S`@sam_data
design_mat <- read_excel("D:/Rprojet/Data/design.xlsx")

# remove sample 78P
design <- design16[design16$Sample_Number != 78,]
otu<- otu[design16$Sample_Number != 78,]
otu_ITS<- otu_ITS[design16$Sample_Number != 78,]
design_mat <- design_mat[design16$Sample_Number != 78,]
otu_all <-otu_all[design16$Sample_Number != 78,]

##### preprocess  #####
ra <- otu/rowSums(otu)
ra_ITS <- otu_ITS/rowSums(otu_ITS)
ra <- cbind(ra, ra_ITS)
taxa <- rbind(taxa, taxa_ITS)
dim(ra)


depth_range <- which(design$Sample_Depth == '0-10 cm')
ra <- ra[depth_range, ]
design <-design[depth_range, ]
design_mat <-design_mat[depth_range, ]
otu_all <-otu_all[depth_range, ]


ra <- t(ra)
ra <- ra[which(rowSums(ra > 0) >= 12),]
taxa <- taxa[rownames(ra),]
otu_all <- t(otu_all)

otu_all <- otu_all[rownames(ra),]

#####indicator species#####
edgeR_Fertility  <- DGEList(counts=ra,
                            group=design$Fertility_Source,
                            genes=taxa)

edgeR_Tillage  <- DGEList(counts=ra,
                          group=design$Tillage,
                          genes=taxa)
edgeR_Cover  <- DGEList(counts=ra,
                        group=design$Cover_Crop,
                        genes=taxa)

otu_norm <- cpm(edgeR_Fertility, normalized.lib.sizes=T, log=F)
indic <- as.data.frame(t(otu_norm))

indic_design_Fertility <- design$Fertility_Source
indic_design_Tillage <- design$Tillage
indic_design_Cover <- design$Cover_Crop


set.seed(8046)
indicatorsp_Fertility  <- multipatt(indic,indic_design_Fertility,func = "r.g",control=how(nperm=1000))
indic_df_Fertility<- indicatorsp_Fertility$sign
indicatorsp_Tillage  <- multipatt(indic,indic_design_Tillage,func = "r.g",control=how(nperm=1000))
indic_df_Tillage<- indicatorsp_Tillage$sign
indicatorsp_Cover  <- multipatt(indic,indic_design_Cover,func = "r.g",control=how(nperm=1000))
indic_df_Cover<- indicatorsp_Cover$sign
# p.value < 0.05
Fertility <- as.matrix(indic_df_Fertility[which(indic_df_Fertility$p.value < 0.05),])
Tillage <- as.matrix(indic_df_Tillage[which(indic_df_Tillage$p.value < 0.05),])
Cover <- as.matrix(indic_df_Cover[which(indic_df_Cover$p.value < 0.05),])

# EDGER
colnames(design)
model_matt_Tillage <- model.matrix(~Tillage, data=design_mat)
model_matt_Fertility <- model.matrix(~Fertility_Source, data=design_mat)
model_matt_Cover <- model.matrix(~Cover_Crop, data=design_mat)


edgeR_Tillage <- DGEList(counts=otu_all, group=design$Tillage, genes=taxa)
edgeR_Fertility <- DGEList(counts=otu_all, group=design$Fertility_Source, genes=taxa)
edgeR_Cover <- DGEList(counts=otu_all, group=design$Cover_Crop, genes=taxa)

edgeR_Tillage <- calcNormFactors(edgeR_Tillage)
edgeR_Fertility <- calcNormFactors(edgeR_Fertility)
edgeR_Cover <- calcNormFactors(edgeR_Cover)

dge_Tillage <- estimateGLMRobustDisp(edgeR_Tillage, design=model_matt_Tillage)
dge_Fertility <- estimateGLMRobustDisp(edgeR_Fertility, design=model_matt_Fertility)
dge_Cover <- estimateGLMRobustDisp(edgeR_Cover, design=model_matt_Cover)



fit_Tillage <- glmFit(dge_Tillage, design=model_matt_Tillage)
fit_Fertility <- glmFit(dge_Fertility, design=model_matt_Fertility)
fit_Cover <- glmFit(dge_Cover, design=model_matt_Cover)

lrt_Tillage <- glmLRT(fit_Tillage)
lrt_Fertility <- glmLRT(fit_Fertility)
lrt_Cover <- glmLRT(fit_Cover)

Tillage_by_edgeR <- topTags(lrt_Tillage, n=Inf, p.value=0.05)
Fertility_by_edgeR <- topTags(lrt_Fertility, n=Inf, p.value=0.05)

Cover_by_edgeR <- topTags(lrt_Cover, n=Inf, p.value=0.05)

edgeR_list_Tillage <- rownames(Tillage_by_edgeR$table)
edgeR_list_Fertility <- rownames(Fertility_by_edgeR$table)
edgeR_list_Cover <- rownames(Cover_by_edgeR$table)

Fertility_list <- intersect(rownames(Fertility), edgeR_list_Fertility)
Tillage_list <- intersect(rownames(Tillage), edgeR_list_Tillage)
Cover_list <- intersect(rownames(Cover), edgeR_list_Cover)



ra_cor <- rcorr(t(ra),type=c("spearman"))
CorrDF <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    from = rownames(cormat)[col(cormat)[ut]],
    to = rownames(cormat)[row(cormat)[ut]],
    cor =(cormat)[ut],
    p = pmat[ut]
  )
}

cor_df <- CorrDF(ra_cor$r,ra_cor$P)
min(cor_df$cor)

cor_df$padj <- p.adjust(cor_df$p, method = "fdr") 




# cor_df_padj <- cor_df[which(abs(cor_df$cor) > 0.76),]
cor_df_padj <- cor_df[which(cor_df$cor  > 0.76),]
cor_df_padj <- cor_df_padj[which(cor_df_padj$padj < 0.005),]


nodeattrib <- data.frame(node = union(cor_df_padj$from,cor_df_padj$to))


nodeattrib$Fertility <- 0
nodeattrib$Tillage <- 0
nodeattrib$Cover <- 0
nodeattrib$times <- 0
for (i in as.character(nodeattrib$node))
{
  if (i %in% Fertility_list == TRUE)
  {nodeattrib[nodeattrib$node==i,"Fertility"] <- "1"
  nodeattrib[nodeattrib$node==i,"color"] <- "Fertility"
  nodeattrib[nodeattrib$node==i,"times"] <- nodeattrib[nodeattrib$node==i,"times"]+1}
  else
  { nodeattrib[nodeattrib$node==i,"color"] <- "NA"}
  if (i %in% Cover_list == TRUE)
  {nodeattrib[nodeattrib$node==i,"Cover"] <- "1"
  nodeattrib[nodeattrib$node==i,"color"] <- "Cover"
  nodeattrib[nodeattrib$node==i,"times"] <- nodeattrib[nodeattrib$node==i,"times"]+1}
  if (i %in% Tillage_list == TRUE)
  {nodeattrib[nodeattrib$node==i,"Tillage"] <- "1"
  nodeattrib[nodeattrib$node==i,"color"] <- "Tillage"
  nodeattrib[nodeattrib$node==i,"times"] <- nodeattrib[nodeattrib$node==i,"times"]+1}
}


nodeattrib <- cbind(nodeattrib,taxa[as.character(nodeattrib$node),])

net <- graph_from_data_frame(cor_df_padj,direct=F, vertices = nodeattrib)
save(nodeattrib,file="D:/Rprojet/soil_V2/network_data/node_0cm_pos.Rdata")

coords_its <- layout_(net,with_fr(niter=9999, grid="nogrid"))

save(coords_its,file="D:/Rprojet/soil_V2/network_data/coords_0cm_pos.Rdata")
save(net,file="D:/Rprojet/soil_V2/network_data/net_0cm_pos.Rdata")

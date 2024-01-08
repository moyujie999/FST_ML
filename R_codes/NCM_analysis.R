rm(list = ls())
library(xlsx) 
library("readxl")

##### import data #####
`ps-orig_16S` <- readRDS("D:/Rprojet/soil_V2/dataset/ps16S.RDS")
`ps-orig_ITS` <- readRDS("D:/Rprojet/soil_V2/dataset/psITS.RDS")

otu <- `ps-orig_16S`@otu_table
otu_ITS <- `ps-orig_ITS`@otu_table
taxa <- `ps-orig_16S`@tax_table
design16 <- `ps-orig_16S`@sam_data
design_mat <- read_excel("D:/Rprojet/Data/design.xlsx")

# remove sample 78P
design <- design16[design16$Sample_Number != 78,]
otu<- otu[design16$Sample_Number != 78,]
otu_ITS<- otu_ITS[design16$Sample_Number != 78,]

design_mat <- design_mat[design16$Sample_Number != 78,]

##### preprocess  #####

ra <- otu_ITS
ra <- t(ra)
dim(ra)


library(Hmisc)
library(minpack.lm)
library(stats4)
spp <- ra
spp<-t(spp)
N <- mean(apply(spp, 1, sum))

#get the mean of species relative abundance in metacommmunity

p.m <- apply(spp, 2, mean)
p.m <- p.m[p.m != 0]
p <- p.m/N

pp <- as.data.frame(p)




spp.bi <- 1*(spp>0)
freq <- apply(spp.bi, 2, mean)
freq <- freq[freq != 0]
#get a table record species percentage and occurrence frequency in metacommunity
C <- merge(p, freq, by=0)
#sort the table according to occurence frquency of each species
C <- C[order(C[,2]),]
C <- as.data.frame(C)
#delete rows containning zero
C.0 <- C[!(apply(C, 1, function(y) any(y == 0))),]
p <- C.0[,2]
freq <- C.0[,3]
names(p) <- C.0[,1]
names(freq) <- C.0[,1]
d = 1/N



##Fit model parameter m (or Nm) using Non-linear least squares (NLS)
m.fit <- nlsLM(freq ~ pbeta(d, N*m*p, N*m*(1 -p), lower.tail=FALSE),start=list(m=0.008))

m.fit #get the m value
m.ci <- confint(m.fit, 'm', level=0.95)
freq.pred <- pbeta(d, N*coef(m.fit)*p, N*coef(m.fit)*(1 -p), lower.tail=FALSE)
pred.ci <- binconf(freq.pred*nrow(spp), nrow(spp), alpha=0.05, method="wilson", return.df=TRUE)
Rsqr <- 1 - (sum((freq - freq.pred)^2))/(sum((freq - mean(freq))^2))
Rsqr
nm <- coef(m.fit)*N
nm



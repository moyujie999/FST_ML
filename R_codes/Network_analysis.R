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

load(file="D:/Rprojet/soil_V2/network_data/net_0cm_pos.Rdata")

load(file="D:/Rprojet/soil_V2/network_data/node_0cm_pos.Rdata")

load(file="D:/Rprojet/soil_V2/network_data/coords_0cm_pos.Rdata")

sensitive_nodes <- rownames(nodeattrib[nodeattrib$color != "NA",])



# triangle vertex shape
mytriangle <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x=coords[,1], y=coords[,2], bg=vertex.color,
          stars=cbind(vertex.size, vertex.size, vertex.size),
          add=TRUE, inches=FALSE)
}
# clips as a circle
add_shape("triangle", clip=shapes("circle")$clip,
          plot=mytriangle)

################################################################# vertex shape
# generic star vertex shape, with a parameter for number of rays
mystar <- function(coords, v=NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size  <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  norays <- params("vertex", "norays")
  if (length(norays) != 1 && !is.null(v)) {
    norays <- norays[v]
  }
  
  mapply(coords[,1], coords[,2], vertex.color, vertex.size, norays,
         FUN=function(x, y, bg, size, nor) {
           symbols(x=x, y=y, bg=bg,
                   stars=matrix(c(size,size/2), nrow=1, ncol=nor*2),
                   add=TRUE, inches=FALSE)
         })
}
# no clipping, edges will be below the vertices anyway
add_shape("star", clip=shape_noclip,
          plot=mystar, parameters=list(vertex.norays=5))


Fertility_nodes1 <- rownames(nodeattrib[nodeattrib$color=="Fertility",])
Cover_nodes1 <- rownames(nodeattrib[nodeattrib$color=="Cover",])
Tillage_nodes1 <- rownames(nodeattrib[nodeattrib$color=="Tillage",])

Fertility_nodes <- rownames(nodeattrib[nodeattrib$Fertility=="1",])
Management_nodes <- rownames(nodeattrib[nodeattrib$Management=="1",])
Cover_nodes <- rownames(nodeattrib[nodeattrib$Cover=="1",])
Tillage_nodes <- rownames(nodeattrib[nodeattrib$Tillage=="1",])

cal <- nodeattrib[Cover_nodes,]
sum(cal$Kingdom == 'Archaea')
sum(cal$Kingdom == 'Bacteria')
sum(cal$Kingdom == 'k__Fungi')

Fertility_Cover <- Fertility_nodes[which(Fertility_nodes %in% Cover_nodes ) ]
Fertility_Tillage <- Fertility_nodes[which(Fertility_nodes %in% Tillage_nodes ) ]
Cover_Tillage <- Cover_nodes1[which(Cover_nodes %in% Tillage_nodes ) ]
Fertility_Cover_Tillage <- Fertility_nodes[which(Fertility_nodes %in% Cover_Tillage ) ]

# color
V(net)$color <- "#666666"

V(net)$color[V(net)$name %in% Cover_nodes] <- "#FFFF33"
V(net)$color[V(net)$name %in% Fertility_nodes] <- "#FF3300"
V(net)$color[V(net)$name %in% Tillage_nodes] <- "#0099FF"

V(net)$color[V(net)$name %in% Fertility_Cover] <- "#FF9900"
V(net)$color[V(net)$name %in% Fertility_Tillage] <- "#9933FF"
V(net)$color[V(net)$name %in% Cover_Tillage] <- "#66CC00"

V(net)$color[V(net)$name %in% Fertility_Cover_Tillage] <- "#ff0099"

sensitive_nodes <- c(Cover_nodes,Fertility_nodes,Tillage_nodes,Fertility_Cover,Fertility_Tillage,Cover_Tillage, Fertility_Cover_Tillage)


V(net)$frame.color <- V(net)$color

V(net)$frame.width <- 2

## set the size
V(net)$size <- V(net)$name
V(net)$size[!V(net)$size %in% sensitive_nodes] <- 1.5
V(net)$size[V(net)$size %in% sensitive_nodes] <- 3
# V(net)$size[V(net)$name %in% key_stone_list] <- 10
nodesizes <- as.numeric(V(net)$size)




# set the shape
V(net)$shape <- V(net)$Kingdom
V(net)$shape[V(net)$shape == 'Bacteria'] <- "circle"
V(net)$shape[V(net)$shape == 'k__Fungi'] <- "triangle"
V(net)$shape[V(net)$shape == 'Archaea'] <- "pie"
# V(net)$shape[V(net)$name %in% key_stone_list] <- "star"
# key_stone_list <- names(key_stone)
nodesizes[V(net)$Kingdom =='k__Fungi' ] <- nodesizes[V(net)$Kingdom =='k__Fungi' ] *2

ecol <- rep("gray80", ecount(net))



edgelist <- get.edgelist(net)
edgelist <- as.data.frame(edgelist)




pdf(paste0("D:/Rprojet/soil_V2/network_figure/network_pos_ensem.pdf"),width=7,height=5)
# par(mfrow=c(2,2), mar=c(1,1,1,1))
cols <- c("#FF3300","#FFFF33","#0099FF","#FF9900","#9933FF", "#66CC00","#ff0099","#FFFFCC", "#FFFFCC","#FFFFCC")
col_group <-c("#FFF68F","#66FFFF")

plot(net,vertex.label=NA,vertex.size=nodesizes, layout=coords_its,edge.color=ecol,
     mark.groups=list(),main = "a 0-10cm",
     mark.col=col_group, mark.border=col_group)
# plot(net,vertex.label=NA,vertex.size=nodesizes, layout=coords_its,edge.color=ecol, main = "b 10-20cm")

# legend("topright",legend=c("F-s","C-s","T-s" ,"FC-s" ,"FT-s" ,"CT-s" ,"FCT-s"),col=cols,
#        bty="n",fill=cols,border=cols)
# 
# legendList <- c("Archaea","Bacteria","fungi")
# legend("topleft",40,legendList,col=c("dodgerblue3","limegreen"),
#        bty="n",pch = c(11,12),lty = c(NA,NA) ,text.font=1.0,cex = 1.5,lwd=2, ncol = 2);
# legend("topleft",legend=c("Module 1","Module 2"),col=cols,
#        bty="n",fill=col_group,border=col_group)
colss <- c("#666666","#666666","#666666","#666666","#666666","#666666")
# legend("bottomright", legend=legendList  , col = colss , bty = "n", pch=c(1,16,17) ,border=colss, text.col=colss )

dev.off()
print(1)





vertex_connectivity(net)

###################################edge
edgelist <- get.edgelist(net)
edgelist <- as.data.frame(edgelist)


edgelist$type <- 0
dim(edgelist)
data <- data.frame(
  from = edgelist[,1],
  to = edgelist[,2]
)
data$type <- 0
data$positive <- 0
for( i in 1:dim(data)[1]){
  
  edge.infos1 <- nodeattrib[nodeattrib$node == data[i,1],]
  edge.infos2 <- nodeattrib[nodeattrib$node == data[i,2],]
  
  if (edge.infos1$Kingdom == "Bacteria" & edge.infos2$Kingdom == "Bacteria"){data[i,3] <- "b2b"}
  if (edge.infos1$Kingdom == "Bacteria" & edge.infos2$Kingdom == "k__Fungi"){data[i,3] <- "b2f"}
  if (edge.infos1$Kingdom == "k__Fungi" & edge.infos2$Kingdom == "Bacteria"){data[i,3] <- "b2f"}
  if (edge.infos1$Kingdom == "k__Fungi" & edge.infos2$Kingdom == "k__Fungi"){data[i,3] <- "f2f"}
  if (edge.infos1$Kingdom == "Archaea" & edge.infos2$Kingdom == "Archaea"){data[i,3] <- "a2a"}
  if (edge.infos1$Kingdom == "k__Fungi" & edge.infos2$Kingdom == "Archaea"){data[i,3] <- "a2f"}
  if (edge.infos1$Kingdom == "Archaea" & edge.infos2$Kingdom == "k__Fungi"){data[i,3] <- "a2f"}
  if (edge.infos1$Kingdom == "Archaea" & edge.infos2$Kingdom == "Bacteria"){data[i,3] <- "a2b"}
  if (edge.infos1$Kingdom == "Bacteria" & edge.infos2$Kingdom == "Archaea"){data[i,3] <- "a2b"}
  # data[i,4]<- cor(ra[data[i,1],],ra[data[i,2],], method = "spearman")
  
  print(i)
}

min(data[,4])

# ra[data[i,1],]

edge_a2a <- data[data$type == "a2a",]
edge_a2b <- data[data$type == "a2b",]
edge_b2b <- data[data$type == "b2b",]
edge_f2f <- data[data$type == "f2f",]
edge_b2f <- data[data$type == "b2f",]
edge_a2f <- data[data$type == "a2f",]

#################################################degree 
deg <- degree(net )
deg.dit <- degree.distribution(net)
mean(deg)
deg.dit


# deg_management <- deg[Management_nodes]
deg_cover <- deg[Cover_nodes]
deg_till <- deg[Tillage_nodes]
deg_fertility <- deg[Fertility_nodes]
mean(deg_cover)


################################edge_density
edge_density(net, loops=F)
ecount(net)/(vcount(net)*(vcount(net)-1)) #for a directed network

################################network diameter

diameter(net, directed=F, weights=NA)
diam <- get_diameter(net, directed=T)

#################################distance
mean_distance(net, directed=F)

################################subgroup
sapply(cliques(net.sym), length) # clique sizes
########################################## strength 
stren <- strength(net)
mean(stren)
#########################
transitivity(net) 
############################
average.path.length(net)

########################modularity
wtc <- cluster_walktrap(net)
modularity(wtc)
modularity(net, membership(wtc))

####################################Vertex and edge betweenness centrality

btw <- betweenness(net)
mean(btw)
btw_edge <- edge_betweenness(net)
mean(btw_edge)
centr_betw_ <- centr_betw(net)


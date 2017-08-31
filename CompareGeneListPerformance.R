#install.packages(  "minerva")
require(ggplot2)
require(ComplexHeatmap)
require(circlize)

require(clValid)
require(ConsensusClusterPlus)
library(reshape2)
library(igraph)

require(doMC)
registerDoMC(32)

load("../DATA/colorRanges.Rdata")
dtMeta=read.table(file="../DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
dtMeta$HER.2[dtMeta$HER.2 == "0"] = NA


load("../DATA/Prot.Rdata")
rownames(dtProt) = dtProt$Gene.Symbol
dtProt$Gene.Symbol = NULL

load("../DATA/geneListsHJ.RData")
geneLists[readLines(con = "Data_and_Figures/CoreNet-Unfiltered.lst"), "CoreNetAll"] = T
geneLists[readLines(con = "Data_and_Figures/CoreNet.lst"), "CoreNetFilt"] = T
geneLists[rownames(dtProt),"All"] =T
geneLists[is.na(geneLists)] =F


samps  = readLines(con = "Data_and_Figures/Outliers.lst")




plotGeneListPerformance = function(lstname, gns,ctClusters){
  # Doesn't work in loop ran manually
  # "All","FDATargets","CoreNetAll","CoreNetFilt","PAM50"
  # lstname="CoreNetFilt"
  # gns=rownames(geneLists)[  geneLists[,lstname] ]
  # ctClusters=4
  gns=gns[gns %in% rownames(dtProt)]
  dtx=log2(dtProt[gns ,samps])
  
  y=ConsensusClusterPlus(data.matrix(dtx), maxK=8, reps=500, pFeature = 0.8, pItem=1, clusterAlg = "hc",innerLinkage =  "ward.D2", finalLinkage = "ward.D2",  distance = "pearson",plot="pdf", title=paste0("byGeneList/",lstname)  )
  ConsPlot = y[[ ctClusters]]
  
  
  haL = rowAnnotation(df = data.frame( dtMeta[colnames(dtx),c("PAM50", "Her2_PAM50_Cor" ,  "LumB_PAM50_Cor" ,  "LumA_PAM50_Cor"  ,"Normal_PAM50_Cor","Basal_PAM50_Cor"  )]),  
                      col=list(PAM50 = xcolorRanges$PAM50,
                               Basal_PAM50_Cor = colorRamp2(c(-1, 1), c("black", xcolorRanges$PAM50["Basal"])),
                               Her2_PAM50_Cor = colorRamp2(c(-1, 1), c("black", xcolorRanges$PAM50["Her2"])),
                               LumA_PAM50_Cor = colorRamp2(c(-1, 1), c("black", xcolorRanges$PAM50["LumA"])),
                               LumB_PAM50_Cor = colorRamp2(c(-1, 1), c("black", xcolorRanges$PAM50["LumB"])),
                               Normal_PAM50_Cor = colorRamp2(c(-1, 1), c("black", xcolorRanges$PAM50["Normal"]))
                      )
  )
  
  
  
  xcol=rainbow(ctClusters)
  names(xcol) = paste0("CSNS",1:ctClusters)
  
  
  xcolNP= c("#d8b365", "#5ab4ac" )
  names(xcolNP) =c("neg", "pos")
  
  xcolH= c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99")
  names(xcolH) =c("Ductal" ,  "DCIS"  ,   "Mucinous" ,"Lobular" , "Mixed" )
  
  dtxhaT = data.frame ( Histology = dtMeta[colnames(dtx),"Histology"], stringsAsFactors = F)
  dtxhaT = cbind(dtxhaT, dtMeta[colnames(dtx),c("HER.2","ER","PR") ])
  dtxhaT = cbind(dtxhaT, data.frame(Consensus=paste0("CSNS",ConsPlot$consensusClass)))
  
  haT =  HeatmapAnnotation(df = dtxhaT,  
                           col=list('Consensus'=xcol, 'ER'=xcolNP, 'PR'=xcolNP, 'HER.2'=xcolNP, 'Histology'=xcolH ))
  
  
  ctClusters=8
  X=y[[ctClusters]][["consensusMatrix"]]
  lines(ecdf(X))
  heatmap()
  rownames(X) = colnames(dtx)
  
  # 
  # pdf(file=paste0("byGeneList/",lstname ,".Annot.pdf"), width=20, height=15)
  # Heatmap(X,  top_annotation = haT, row_names_side = "left" , show_row_dend = FALSE ) + haL
  # dev.off()

  a = silhouette(cutree( ConsPlot$consensusTree, ctClusters), 1- ConsPlot$consensusMatrix  )[,3]
  b = ConsPlot$consensusClass
  id = names(ConsPlot$consensusClass)
  names(b) = NULL
  data.frame(id, a,lstname,b)
  
}







#install.packages("coloramp2")

#haT = HeatmapAnnotation(df = data.frame(ProtSub9506= dtMeta[colnames(dtx),"ProtSub9506" ]),  col=list(ProtSub9506 = xcolorRanges$ProtSub9506))


x=foreach (lst = c("All","FDATargets","CoreNetAll","CoreNetFilt","PAM50"), .combine = rbind) %dopar% {
  gns=rownames(geneLists)[ geneLists[,lst]]
  plotGeneListPerformance(lst,gns,4)
}


write.table(x,file="Data_and_Figures/ConsensusClases.tab",sep="\t", quote = T)
pdf("Data_and_Figures/SilhouetteGeneListCompare.pdf")
ggplot(x, aes(factor(lstname), a)) + geom_boxplot() + ylim(0.7,1)
dev.off()




calcClusterCompare = function(cluA,cluB){
  grps=list()
  for (grp in unique(cluA)  ){
    grps[[grp]] = names(cluA)[cluA == grp]
  }
  
  for (grp in unique(cluB)  ){
    grps[[grp]] = names(cluB)[cluB == grp]
  }
  
  nodes = data.frame(name=names(grps),idx=1:length(grps) )
  nodes$idx = nodes$idx -1
  rownames(nodes) = nodes$name
  
  grpAll = intersect(names(cluA),names(cluB))
  
  edges = foreach (a = unique(cluA), .combine = rbind) %:%
    foreach (b = unique(cluB), .combine = rbind) %do% {
      grpA = grps[[a]]
      grpB = grps[[b]]
      grpA = intersect(grpA,grpAll)
      grpB = intersect(grpB,grpAll)
      weight = sum(grpA %in% grpB)
      c(nodes[a,"idx"],nodes[b,"idx"],weight)
    }
  edges=data.frame(edges, stringsAsFactors = F, row.names = 1:nrow(edges))
  colnames(edges) = c("source","target","weight")
  edges=edges[edges$weight > 0,]
  list(nodes=nodes,edges=edges)
}


idx=x$lstname == "CoreNetFilt"
cA = paste0("CSNS",x$b[idx],"_CoreNetFilt")
names(cA) = x$id[idx]
  
cB = dtMeta$PAM50[dtMeta$id %in% names(cA)]
names(cB) = dtMeta$id[dtMeta$id %in% names(cA)]

res=calcClusterCompare(cA,cB)
library(networkD3)

sankeyNetwork(Links = res$edges, Nodes = res$nodes, Source = "source",
              Target = "target", Value = "weight", NodeID = "name",fontSize = 14)


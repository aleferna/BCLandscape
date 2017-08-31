require(data.table)
require(stringr)

require(reshape2)
require(igraph)
require(doMC)
require(ggplot2)
registerDoMC(35)
#install.packages("beanplot")
require(beanplot)


getMTX = function(tpCor, tpSRC){
  if(tpSRC == "RNA"){
    load("../DATA/RNA.Rdata")
    rownames(dtRNA) = dtRNA$geneName
    dtRNA$geneName = NULL
    mtx=t(dtRNA)
  }
  
  if(tpSRC == "RNA45"){ 
    load("../DATA/Prot.Rdata")
    dtProt$Gene.Symbol = NULL
    load("../DATA/RNA.Rdata")
    rownames(dtRNA) = dtRNA$geneName
    dtRNA$geneName = NULL
    dtRNA = dtRNA[,colnames(dtProt)[]]
    mtx=t(dtRNA)
  }
  
  if(tpSRC == "prot"){  
    load("../DATA/Prot.Rdata")
    rownames(dtProt) = dtProt$Gene.Symbol
    dtProt$Gene.Symbol = NULL
    mtx=t(log2(dtProt))
  }
  mtx
}

getCorum = function(gns){
  
  dtCorum = read.table("../DATA/CORUM_30102016_coreComplexes.txt",sep="\t", 
                       comment.char = '?',quote = '?', header = T, stringsAsFactors = F)
  dtx = foreach (x = dtCorum$subunits.Gene.name., .combine=rbind) %do% {
    x=unlist(str_split(x,";"))
    x=x[x %in% gns]
    if (length(x) > 1){
      t(combn(x,2))  
    }else{
      NULL
    }
  }
  colnames(dtx) = c("node1","node2")
  dtx=data.frame(dtx,stringsAsFactors = F)
  idx = dtx$node1 > dtx$node2
  x = dtx$node1[idx]
  dtx$node1[idx] = dtx$node2[idx]
  dtx$node2[idx] = x
  dtx$X = T
  unique(dtx)
  
}

getBioGrid = function(gns){
  dtBioGrid = read.table("/Z/aleferna/GeneLists/BioGrid/BIOGRID-ALL-3.4.137.GenePairs",header = T,comment.char = "#",sep = "\t",stringsAsFactors = F)
  idx = dtBioGrid$Entrez %in% gns & dtBioGrid$Gene %in% gns
  dtBioGrid = dtBioGrid[idx,]
  colnames(dtBioGrid) = c("node1","node2")
  #remove duplicates
  dtBioGrid = dtBioGrid[!duplicated(dtBioGrid),]
  #sort the rows
  idx = dtBioGrid$node1 > dtBioGrid$node2
  x = dtBioGrid$node1[idx]
  dtBioGrid$node1[idx] = dtBioGrid$node2[idx]
  dtBioGrid$node2[idx] = x
  #remove self loops
  idx = dtBioGrid$node1 != dtBioGrid$node2
  dtBioGrid = dtBioGrid[idx,]
  dtBioGrid = data.table(dtBioGrid)
  dtBioGrid$X = T
  dtBioGrid
}


############################################################################################################################################################
################################################               Genome Wide Correlation          ############################################################
############################################################################################################################################################
doViolinALL = function(tpSRC,tpCor, mtx, dtKnown, srcKnown){
  
  mtxCor=cor( mtx, method = tpCor)
  edges=melt(mtxCor)
  colnames(edges) = c("node1", "node2","cor")
  edges$node1 = paste(edges$node1)
  edges$node2 = paste(edges$node2)
  idx = edges$node1 < edges$node2
  edges=edges[idx,]
  edges=data.table(edges)
  idx=NULL
  
  edges=merge(edges,dtKnown, all=T, by=c("node1","node2"))
  edges[is.na(edges$X), X := F]
  
  xedges=data.frame(edges)
  
  pdf(file=paste0("Figures/CorrelationVsBioGrid/",tpSRC,".",tpCor,".",srcKnown,".Violin-ALL.pdf"), height = 8, width=8)
  
  
  
  beanplot(xedges$cor ~ xedges$X, ll = 0.04,
           main = paste("Correlation of All Genes: ",srcKnown) , side = "both", xlab=paste(tpSRC, tpCor),
           col = list("#be1200", c("#006eff", "#FFFF00")),
           axes=T,bw="nrd",what = c(1,1,1,0), ylim = c(-0.5,1.2))
  
  
  dev.off()
}



############################################################################################################################################################
################################################               SubGroup Correlation          ###############################################################
############################################################################################################################################################

doSubGroup = function(tpSRC,tpCor, mtx, dtKnown, srcKnown){
  load("../DATA/geneListsHJ.RData")
  
  geneLists$Nik.Zainal_Drivers =NULL
  geneLists$FDATargets = NULL
  geneLists$CancerSensus = NULL
  geneLists$PAM50 = NULL
  geneLists$DrugableGenome = NULL
  rownames(geneLists) = geneLists$GeneName
  geneLists$GeneName = NULL
  geneLists$No.enrichment = NULL
  geneLists$FredlundNet = NULL
  geneLists$TumorGenes9506 = NULL 
  
  xedges=foreach (ListName = colnames(geneLists), .combine = rbind ) %do% {
    #ListName = colnames(geneLists)[1]
    gns = rownames( geneLists)[geneLists[,ListName]]
    
    # 
    
    gns = gns[gns %in% colnames(mtx)]
    print(c(ListName, length(gns)))
    
    mtxCor=cor(mtx[,gns],method = tpCor)
    edges=melt(mtxCor)
    colnames(edges) = c("node1", "node2","cor")
    edges$node1 = paste(edges$node1)
    edges$node2 = paste(edges$node2)
    idx = edges$node1 < edges$node2
    edges=edges[idx,]
    edges=data.table(edges)
    idx=NULL
    edges=merge(edges,dtKnown, all.x=T, by=c("node1","node2"))
    edges[is.na(edges$X), X := F]
    edges$ListName = ListName
    edges  
  }
  
  
  
  pdf(file=paste0("Figures/CorrelationVsBioGrid/",tpSRC,".",tpCor,".",srcKnown,".Violin.pdf"), height = 8, width=12)
  idx = xedges$ListName %in% c("Erythrocyte"  , "Plasma",        "Adipose"     ,  "ECM.1"      ,   "ECM.2"     ,    "Luminal"     )
  
  beanplot(xedges$cor[idx] ~ xedges$X[idx]*xedges$ListName[idx], ll = 0.04,
           main = paste0(srcKnown," Correlation by Gene Cluster"), side = "both", xlab=paste(tpSRC, tpCor),
           col = list("#be1200", c("#006eff", "#FFFF00")),
           axes=T,bw="nrd",what = c(1,1,1,0), ylim = c(-0.5,1.2))
  
  
  idx = xedges$ListName %in% c("Immune"      ,  "Golgi"   ,      "Basal"   ,      "Mitochondrion" ,"Proliferation"  , "Transcription") 
  
  beanplot(xedges$cor[idx] ~ xedges$X[idx]*xedges$ListName[idx], ll = 0.04,
           main = paste0(srcKnown," Correlation by Gene Cluster"), side = "both", xlab=paste(tpSRC, tpCor),
           col = list("#be1200", c("#006eff", "#FFFF00")),
           axes=T,bw="nrd",what = c(1,1,1,0), ylim = c(-0.5,1.2))
  
  dev.off()
  
}


############################################################################################################################################################
################################################              RNA vs Protein @ protein Complx         ###############################################################
############################################################################################################################################################

doKnown.RNAvPROT.ByComplex = function(mtxR, mtxP, tpCor){
  
  
  geneLists=list()
  geneLists$GINS = c("GINS1","GINS2","GINS3","GINS4")
  geneLists$MCM = c("MCM2","MCM3","MCM4","MCM5","MCM6","MCM7")
  geneLists$Condensing_I = c("NCAPD2","NCAPG","NCAPH","SMC2","SMC4")
  geneLists$Condensing_II = c("NCAPD3","NCAPG2","NCAPH2","SMC2","SMC4")
  geneLists$Mitotic_14s_cohesin_1 = c("RAD21","SMC1A","SMC3","STAG1")
  geneLists$DNA_PolAlpha = c("CHEK1","PRIM1","PRIM2","POLA1","POLA2")
  
  xedges=foreach (ListName = names(geneLists), .combine = rbind ) %do% {
    #ListName = names(geneLists)[1]
    gns = geneLists[[ListName]]
    gns = gns[gns %in% colnames(mtxP)]
    gns = gns[gns %in% colnames(mtxR)]
    
    mtxCor=cor(mtxR[,gns],method = tpCor)
    edgesR=melt(mtxCor)
    edgesR$Src = "RNA"
    
    mtxCor=cor(mtxP[,gns],method = tpCor)
    edgesP=melt(mtxCor)
    edgesP$Src = "Protein"
    edges = rbind(edgesR,edgesP)
    
    colnames(edges) = c("node1", "node2","cor","Src")
    edges$node1 = paste(edges$node1)
    edges$node2 = paste(edges$node2)
    idx = edges$node1 < edges$node2
    edges=edges[idx,]
    edges=data.table(edges)
    idx=NULL
    edges$ListName = ListName
    edges  
  }
  
  pdf(file=paste0("Figures/CorrelationVsBioGrid/Complex.",tpCor,".RNAvsProtein.pdf"), height = 8, width=12)
  
  beanplot(show.names = T, xedges$cor ~ xedges$Src * xedges$ListName, ll = 0.04,
           main = "RNA vs Protein Correlation of Complexes", side = "both", xlab=paste( tpCor),
           col = list("red", "orange"),
           axes=T,bw="nrd",what = c(1,1,1,0), ylim = c(0,1.2) )
  
  dev.off()
  
  
}

############################################################################################################################################################
################################################        Correlation of BioGrid Genes RNA vs Prot      ###############################################################
############################################################################################################################################################

doKnown.RNAvsProt.All = function(mtxR, mtxP, tpCor, dtKnown, srcKnown){
  
  
  
  mtxCor=cor(mtxR,method = tpCor)
  edgesR=melt(mtxCor)
  edgesR$Src = "RNA"
  
  mtxCor=cor(mtxP,method = tpCor)
  edgesP=melt(mtxCor)
  edgesP$Src = "Protein"
  
  edges = rbind(edgesR,edgesP)
  
  colnames(edges) = c("node1", "node2","cor","Src")
  edges$node1 = paste(edges$node1)
  edges$node2 = paste(edges$node2)
  idx = edges$node1 < edges$node2
  edges=edges[idx,]
  edges=data.table(edges)
  idx=NULL
  edges=merge(edges,dtKnown, all.x=T, by=c("node1","node2"))
  edges[is.na(edges$X), X := F]
  
  pdf(file=paste0("Figures/CorrelationVsBioGrid/",srcKnown,".Known.",tpCor,".RNAvsProtein.pdf"), height = 8, width=12)
  
  beanplot(show.names = T, edges$cor ~ edges$X * edges$Src, ll = 0.04,
           main = paste("RNA vs Protein Correlation:",srcKnown), side = "both", xlab=paste(srcKnown, tpCor),
           col = list("red", "orange", "blue", "lightblue"),
           axes=T,bw="nrd",what = c(1,1,1,0), ylim = c(-0.5,1.2) )
  
  dev.off()
}











############################################################################################################################################################
################################################               p-value calculation           ###############################################################
############################################################################################################################################################



doPValCalc = function(tpSRC, tpCor, mtx, dtKnown, srcKnown){
  
  mtxCor=cor( mtx, method = tpCor)
  edges=melt(mtxCor)
  colnames(edges) = c("node1", "node2", "cor")
  edges$node1 = paste(edges$node1)
  edges$node2 = paste(edges$node2)
  idx = edges$node1 < edges$node2
  edges=edges[idx,]
  edges=data.table(edges)
  idx=NULL
  
  edges=merge(edges, dtKnown, all=T, by=c("node1","node2"))
  edges[is.na(edges$X), X := F]
  
  xedges=data.frame(edges)
  
  
  load("../DATA/geneListsHJ.RData")
  rownames(geneLists) = geneLists$GeneName
  geneLists$GeneName = NULL
  gns=colnames(mtx)
  geneLists$All =F 
  geneLists[gns,"All"] = T
  geneLists[is.na(geneLists)] = F 
  
  
  xstats=foreach (ListName = colnames(geneLists), .combine=rbind, .errorhandling = "remove") %dopar% {
    gns = rownames( geneLists)[geneLists[,ListName]]
    gns = gns[gns %in% colnames(mtx)]
    mtxCor=cor(mtx[,gns],method = tpCor)
    edges=melt(mtxCor)
    colnames(edges) = c("node1", "node2","cor")
    edges$node1 = paste(edges$node1)
    edges$node2 = paste(edges$node2)
    idx = edges$node1 < edges$node2
    edges=edges[idx,]
    edges=data.table(edges)
    idx=NULL
    edges=merge(edges,dtKnown, all.x=T, by=c("node1","node2"))
    edges[is.na(edges$X), X := F]
    edges$ListName = ListName
    
    fg=edges$cor[edges$X]
    bg=sample(edges$cor[!edges$X], sum(edges$X))
    p=wilcox.test(fg,bg, alternative = "greater")
    print(ListName)
    data.frame(correlType=tpCor, dataSource=tpSRC, knownSrc=srcKnown, test=ListName, 
               pval=p$p.value, tstat=p$statistic, 
               alternative=p$alternative , 
               length(fg), fgAvg = mean(fg), fgStd = sd(fg),
               length(bg), bgAvg = mean(bg), bgStd = sd(bg),
               delta = mean(fg) - mean(bg)
    )
  }
  xstats
}


############################################################################################################################################################
################################################               Main Flow           ###############################################################
############################################################################################################################################################
tpCor="pearson"

load("../DATA/Prot.Rdata")
rownames(dtProt) = dtProt$Gene.Symbol
dtProt$Gene.Symbol = NULL
dtProt = log2(dtProt)

load("../DATA/RNA.Rdata")
rownames(dtRNA) = dtRNA$geneName
dtRNA$geneName = NULL

gns=intersect(rownames(dtRNA),rownames(dtProt) )

dtBioGrid = getBioGrid(gns)
dtCorum = getCorum(gns)

gns=intersect(gns, unique(c(dtBioGrid$node1, dtBioGrid$node2))  )
smps=intersect(colnames(dtProt),colnames(dtRNA))

mtxP = t(dtProt[gns,smps])
mtxR = t( dtRNA[gns,smps])

doKnown.RNAvPROT.ByComplex(mtxR, mtxP, tpCor)

doKnown.RNAvsProt.All(mtxR, mtxP, tpCor, dtCorum, "CORUM")
doKnown.RNAvsProt.All(mtxR, mtxP, tpCor, dtBioGrid, "BioGrid" )

doSubGroup("RNA", tpCor, mtxR, dtBioGrid, "BioGrid" )
doSubGroup("RNA", tpCor, mtxR, dtCorum, "CORUM" )
doSubGroup("Prot", tpCor, mtxP, dtBioGrid, "BioGrid" )
doSubGroup("Prot", tpCor, mtxP, dtCorum, "CORUM" )



stats = foreach (srcKnown = c("BioGrid","CORUM"), .combine=rbind) %dopar% {
  if (srcKnown == "BioGrid"){
    dtKnown = dtBioGrid
  }else{
    dtKnown = dtCorum
  }
  foreach (tpSRC = c("RNA","Prot"), .combine=rbind) %dopar% {
    if (tpSRC == "RNA"){
      doPValCalc(tpSRC, tpCor, mtxR, dtKnown, srcKnown )    
    }else{
      doPValCalc(tpSRC, tpCor, mtxP, dtKnown, srcKnown )  
    }
  }
}
write.table(stats, file=paste0("Figures/CorrelationVsBioGrid/Stats.tab"), sep="\t", row.names = F, quote = F)








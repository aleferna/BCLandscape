#install.packages(  "clValid")
#install.packages(  "kohonen")
#install.packages(  "minerva")
require(amap)
require(data.table)
require(doMC)
require(fitdistrplus)
require(logspline)
require(ComplexHeatmap)
require(curl)
require(networkD3)
require(circlize) 
require(clValid)
require(ConsensusClusterPlus)
require(reshape2)
require(igraph)
require(mixtools)
#install.packages("networkD3")


registerDoMC(20)

reportHighExpVarGenes = function(){
  #1. Assume that some genes have the "same" expression while others have differentially regulated expression
  #2. Assume that both of these populations are gaussian
  #3. Filter all values that are at least 3 std devs away from the "low" population
  
  sdQ=function(x){ sd(sort(x)[2:(length(x)-1)]) }
  x=log2(apply(dtProt,1,sdQ))
  plot(density(x))
  
  gauMix=normalmixEM(x)
  if (gauMix$mu[1] > gauMix$mu[2]){
    idxfg  = 1  
    idxbg  = 2
  }else{
    idxfg  = 2  
    idxbg  = 1
  }
  gmix=data.frame(fgMu = gauMix$mu[idxfg], 
                  fgSigma = gauMix$sigma[idxfg], 
                  fgLambda = gauMix$lambda[idxfg], 
                  bgMu = gauMix$mu[idxbg], 
                  bgSigma = gauMix$sigma[idxbg], 
                  bgLambda = gauMix$lambda[idxbg])
  
  fg=rnorm(1E5, gmix$fgMu,   gmix$fgSigma)
  bg=rnorm(1E5, gmix$bgMu,   gmix$bgSigma)
  
  pdf(file=paste0("Figures/",SetName,"/VarianceGaussianMix.pdf"))
  d=density(x)
  plot(d$x,d$y/max(d$y),  type='l', lwd=5)
  d=density(fg)
  lines(d$x, gmix$fgLambda*d$y/max(d$y),  col='blue', lwd=2)
  
  d=density(bg)
  lines(d$x, gmix$bgLambda*d$y/max(d$y),  col='red', lwd=2)
  
  ###Simulate the data
  bg=rnorm(1E5*gmix$bgLambda, mean =  gmix$bgMu, sd = gmix$bgSigma)
  fg=rnorm(1E5*gmix$fgLambda, mean =  gmix$fgMu, sd = gmix$fgSigma)
  d=density(c(fg,bg))
  lines(d$x,d$y/max(d$y), col='green')
  
  ###Estimate ROC analysis
  roc=foreach (qt = 1:99 , .combine=rbind ) %dopar% {
    ## generate multiple iterations 
    
    xroc=foreach ( i = 1:10 , .combine=rbind ) %do% {
      bg=rnorm(1E5*gmix$bgLambda, mean =  gmix$bgMu, sd = gmix$bgSigma)
      fg=rnorm(1E5*gmix$fgLambda, mean =  gmix$fgMu, sd = gmix$fgSigma)
      t=quantile(c(fg,bg),qt/100.0)
      TP = sum(fg>t)
      FP = sum(bg>t)
      TN = sum(bg<t)
      FN = sum(fg<t)
      Sensitivity =  TP / (TP + FN)
      Specificty = TN / (TN + FP)
      TPr = TP / (TP + FP)
      score = sum(fg>t) - sum(bg>t)
      data.frame(qt, t, Sensitivity,Specificty, TP, FP, TN , FN, TPr, score)
    }
    colMeans(xroc)
  }
  
  roc = data.frame(roc)
  idx = roc$score == max(roc$score)
  t=roc$t[idx]
  t=ceiling(t*2)/2
  idx2 = which.min(abs(roc$t - t))
  
  
  plot(roc$t, roc$score,   type = 'l', lwd=2, xlab='Variance Threshold', ylab='Score (#TP - #FP)' ,  main=paste('Threshold=',2^roc$t[idx]))
  points(roc$t[idx],  roc$score[idx], col='orange', cex=2, xlim=c(0,1), pch=16)
  points(roc$t[idx2], roc$score[idx2], col='red', cex=2, xlim=c(0,1), pch=16)
  
  plot( 1-roc$Specificty, roc$Sensitivity, type='l', lwd=3, ylab="True Positive Rate", xlab='False Positive Rate', main='ROC Analysis') 
  points(1-roc$Specificty[idx],roc$Sensitivity[idx], col='orange', pch=16, cex=2)
  points(1-roc$Specificty[idx2],roc$Sensitivity[idx2], col='red', pch=16, cex=2)
  
  
  gns=names(x[x>t])
  text(0.5, 0.5, labels = paste0("Selected Genes: ", length(gns)) )
  text(0.5, 0.47, labels = paste0("Estimated True Positives: ~",  as.integer( length(gns)*roc$TPr[idx2]   ) ))
  text(0.5, 0.44, labels = paste0("Estimated False Positives: ~", as.integer( length(gns)* (1 - roc$TPr[idx2] )   ) ))
  
  dev.off()
  names(x[x>t])
  
}

###4.2 Adipose and Immune genes that may be from contaminations. 
reportTumorGenes = function(dtx){
  gnsNT = readLines("../DATA/None-tumor_proteins_161116.txt")
  idx = !(rownames(dtx) %in% gnsNT)
  rownames(dtx)[idx]
}

generateGraph = function(dtx){
  mtx = cor(t(dtx),method = "pearson")
  xpairs=subset(melt(mtx), value!=0)
  xpairs=data.table::data.table(xpairs)
  colnames(xpairs) = c("Source","Target","pearson")
  xpairs$Source = paste(xpairs$Source)
  xpairs$Target = paste(xpairs$Target)
  xpairs=xpairs[xpairs$Source > xpairs$Target,] 
  graph_from_data_frame(xpairs,directed=F)
}

getGraph = function(dtx, sdMult){
  #dtx=log2(dtProt[rownames(dtProt) %in% gnsHighVar,])
  mtx = cor(t(dtx),method = "pearson")
  xpairs=subset(melt(mtx), value!=0)
  xpairs=data.table::data.table(xpairs)
  colnames(xpairs) = c("Source","Target","pearson")
  xpairs$Source = paste(xpairs$Source)
  xpairs$Target = paste(xpairs$Target)
  xpairs=xpairs[xpairs$Source > xpairs$Target,]
  
  
  pdf(file=paste0("Figures/",SetName,"/PearsonNetCorrelation.pdf"))
  
  plot(density(xpairs$pearson), main='Distribution of pearson correlations', xlim=c(-1,1))
  xx=rnorm(1E5, mean=mean(xpairs$pearson), sd = sd(xpairs$pearson))
  lines(density(xx), col='red')
  tpear = mean(xpairs$pearson) + sdMult * sd(xpairs$pearson)
  abline(v=tpear)
  dev.off()
  
  idx = abs(xpairs$pearson) > tpear
  sum(idx)
  
  xgraph = graph_from_data_frame(xpairs[idx,],directed=F)
  
  xgraph
  
}





reportNonOutliers = function(dtx, ctCluster){
  pdf(file=paste0("Figures/",SetName,"/Outlier_Detection.pdf"))
  t=0.5
  M=data.matrix(cor(dtx, method = 'pearson')  )
  diag(M) = 0
  M=apply(M,1,max) 
  M=sort(M)
  plot(M, col='green',pch=16, cex=1.5, main='Max correlation to a different sample' )
  abline(h=t)
  points(M[M<t],labels = names(M[M<t]), col='red', pch=16, cex=1.5)
  
  
  outliers=names(M[M<t])
  notoutliers=names(M[M>t])
  dtx$RAND = apply(dtx,1,function(x){sample(x,1)})
  M=foreach (i = 1:10, .combine = rbind, .errorhandling = 'remove' ) %dopar% {
    Z=foreach (ii = 1:100, .combine = c, .errorhandling = 'remove' ) %do% {
      dtx$RAND = apply(dtx,1,function(x){sample(x,1)})
      X=foreach (ii = 1:100, .combine = rbind, .errorhandling = 'remove' ) %do% {
        idxX = sample(nrow(dtx), nrow(dtx)*runif(1))
        idxY = sample(ncol(dtx), ncol(dtx)*runif(1))
        X=dtx[idxX,idxY ]
        melt(cor(X, method = 'pearson')  )
      }
      X=data.table(X)
      Y=X[,list(RPear=mean(value)),by=list(Var1,Var2)]
      Y=reshape2::acast(Y,Var1~Var2,value.var = 'RPear')
      diag(Y) = 0
      Y=apply(Y,1,max)
      names(Y)[Y<t]
    }
    table(Z)[colnames(dtx)]
  }
  M[is.na(M)] = 0
  colnames(M) = colnames(dtx)
  colnames(M)[ colnames(M) %in% outliers] = paste0(">>",colnames(M)[ colnames(M) %in% outliers])
  a=Heatmap(t(M), name = "Random Iterations")
  draw(a)
  
  dev.off()
  notoutliers
}


generateConsClustering = function(dtx, ctClusters){
  
  y=ConsensusClusterPlus(d =  data.matrix(dtx), maxK=15, reps=2000, pFeature = 0.8, pItem=0.8, clusterAlg = "hc",
                         innerLinkage =  "average", finalLinkage = "average",  distance = "pearson",plot="pdf", 
                         title=paste0("Figures/",SetName,"/"))
  
  
  ConsPlot = y[[ ctClusters]]
  X=y[[ctClusters]][["consensusMatrix"]]
  Xclu=y[[ctClusters]]$consensusClass
  
  xx = names(table(Xclu))[table(Xclu) < 3]
  for (i in 1:length(Xclu)){
    
    Xclu[i] = paste0("CONS",Xclu[i])
    
  }
  
  
  rownames(X) = colnames(dtx)
  colnames(X) = colnames(dtx)
  
  xcolorRanges$ConsensusClust = rainbow(length(unique(Xclu)))
  names(xcolorRanges$ConsensusClust) = sort(unique(Xclu))
  xcolorRanges <<- xcolorRanges
  
  dtMeta$ConsensusClust = ""
  dtMeta[names(Xclu ), "ConsensusClust"] = Xclu
  dtMeta <<- dtMeta
  
  
  
  
  
  
  
  haT = HeatmapAnnotation(df = dtMeta[colnames(dtx), c( "PAM50", "ConsensusClust","HER.2","ER","PR","Histology")]  ,
                          col=list(
                            PAM50 = xcolorRanges$PAM50,
                            ConsensusClust = xcolorRanges$ConsensusClust,
                            'ER'= xcolorRanges$xcolNPE, 'PR'= xcolorRanges$xcolNPR, 'HER.2'= xcolorRanges$xcolNPH, 'Histology'= xcolorRanges$xcolH
                          ))
  
  dtMeta=dtMeta[colnames(dtx),]
  Y=apply(dtMeta[,c("Her2_PAM50_Cor","LumB_PAM50_Cor","LumA_PAM50_Cor","Normal_PAM50_Cor","Basal_PAM50_Cor")], 1, 
          function(xx){
            xx=unlist(as.double( xx))
            xx=xx[order(xx,decreasing = T)]
            min(xx[1] - xx[2],0.5)
          } 
  )
  dtMeta$PAM50_Conf = Y
  
  xcolPAM50x = colorRamp2(c(-0.75, 0, 0.75), c("red", "black","green"))
  haL = rowAnnotation(df =  dtMeta[colnames(dtx),c("PAM50_Conf", "PAM50", "Her2_PAM50_Cor" ,  "LumB_PAM50_Cor" ,  "LumA_PAM50_Cor"  ,"Normal_PAM50_Cor","Basal_PAM50_Cor"  )],  
                      col=list(PAM50 = xcolorRanges$PAM50,
                               Basal_PAM50_Cor  = xcolPAM50x,
                               Her2_PAM50_Cor   = xcolPAM50x,
                               LumA_PAM50_Cor   = xcolPAM50x,
                               LumB_PAM50_Cor   = xcolPAM50x,
                               Normal_PAM50_Cor = xcolPAM50x,
                               PAM50_Conf = colorRamp2(c(0, 0.5), c("black", "purple"))
                      )
  )  
  
  haR = rowAnnotation(df =  dtMeta[colnames(dtx),c("Invasive.Tumour","DCIS","Stroma","Lymphocytes","Adipose","Normal.Epith")] ,
                      col=list( Invasive.Tumour  = colorRamp2(c(0, 100), c("black","brown")),
                                DCIS             = colorRamp2(c(0, 100), c("black", "cyan")),
                                Stroma   = colorRamp2(c(0, 100), c("black", "green")),
                                Lymphocytes   = colorRamp2(c(0, 100), c("black", "orange")),
                                Adipose = colorRamp2(c(0, 100), c("black", "blue")),
                                Normal.Epith = colorRamp2(c(0, 100), c("black", "purple"))
                      )
                      
  )
  
  
  
  
  
  
  x=haR + Heatmap(
    matrix = data.matrix(X), 
    top_annotation = haT,
    clustering_method_columns = "average" , 
    clustering_distance_columns  = "pearson"    ,
    clustering_method_rows = "average" , 
    clustering_distance_rows  = "pearson" ,name = "Consensus Measure"
  )  + haL
  
  pdf(paste0("Figures/",SetName,"/ConsensusClust.pdf"), height = 12, width = 18)
  draw(x,heatmap_legend_side = "bottom")
  dev.off()
  
  
  
  
  #dev.off()
  
  writeLines(text = paste(names(ConsPlot$consensusClass ), paste0("GRP",ConsPlot$consensusClass) , sep="\t"), paste0("Figures/",SetName,"/ConsensusClustering.lst"))
  
  
}


generateSampleGraph = function(){
  #dtx=log2(dtProt[idx,])
  xx=cor(dtx,method = "pearson")
  xx=graph_from_adjacency_matrix(xx,weighted = T,mode = "undirected")
  xx= set_vertex_attr(xx,"type", V(xx) , dtMeta[ names(V(xx)),"PAM50"])
  igraph::write_graph(xx, file = "Temp/samps.gml", format = "gml")
  plot(xx)
  
}


generateGeneListHeatmaps = function(){
  dtx=log2(dtProt)
  load("../DATA/dtGeneCategories.Rdata")
  load("../DATA/RNA.Rdata")
  
  gnsPAM50=c("UBE2T","BIRC5","NUF2","CDC6","CCNB1","TYMS","MYBL2","CEP55","MELK","NDC80","RRM2","UBE2C","CENPF","PTTG1","EXO1","ORC6L","ANLN","CCNE1",
             "CDC20","MKI67","KIF2C","ACTR3B","MYC","EGFR","KRT5","PHGDH","CDH3","MIA","KRT17","FOXC1","SFRP1","KRT14","ESR1","SLC39A6","BAG1","MAPT",
             "PGR","CXXC5","MLPH","BCL2","MDM2","NAT1","FOXA1","BLVRA","MMP11","GPR160","FGFR4","GRB7","TMEM45B","ERBB2")
  
  
  
  xcolNPE= c("grey", "blue", "white" )
  names(xcolNPE) =c("neg", "pos","NA")
  xcolNPR= c("grey", "orange", "white" )
  names(xcolNPR) =c("neg", "pos","NA")
  xcolNPH= c("grey", xcolorRanges$PAM50["Her2"], "white" )
  names(xcolNPH) =c("neg", "pos","NA")
  xcolH= c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99")
  names(xcolH) =c("Ductal" ,  "DCIS"  ,   "Mucinous" ,"Lobular" , "Mixed" )
  
  dtMeta$ConsensusClust[dtMeta$ConsensusClust == ""] = "NA"
  xcolorRanges$ConsensusClust["NA"] = "grey"
  
  xcolPAM50x = colorRamp2(c(-0.75, 0, 0.75), c("red", "black","green"))
  
  dtMeta$HClustFull = NA
  Y=apply(dtMeta[,c("Her2_PAM50_Cor","LumB_PAM50_Cor","LumA_PAM50_Cor","Normal_PAM50_Cor","Basal_PAM50_Cor")], 1, 
          function(xx){
            xx=unlist(as.double( xx))
            xx=xx[order(xx,decreasing = T)]
            min(xx[1] - xx[2],0.5)
          } 
  )
  dtMeta$PAM50_Conf = Y
  
  
  
  haB = columnAnnotation(df =  dtMeta[colnames(dtx),c("PAM50_Conf",  "Her2_PAM50_Cor" ,  "LumB_PAM50_Cor" ,  "LumA_PAM50_Cor"  ,"Normal_PAM50_Cor","Basal_PAM50_Cor" )],  
                         col=list(Basal_PAM50_Cor  = xcolPAM50x,
                                  Her2_PAM50_Cor   = xcolPAM50x,
                                  LumA_PAM50_Cor   = xcolPAM50x,
                                  LumB_PAM50_Cor   = xcolPAM50x,
                                  Normal_PAM50_Cor = xcolPAM50x,
                                  PAM50_Conf = colorRamp2(c(0, 0.5), c("black", "purple"))
                         )
  )  
  
  
  pdf(paste0("Figures/",SetName,"/GeneLists.pdf"),width = 15, height = 20)
  
  ########################################################################################################
  ####  ProtPAM37
  ########################################################################################################
  idx=rownames(dtx) %in% gnsPAM50
  d=as.dist(1-cor(dtx[idx,], method="pearson"))
  d=hclust(d,method="average")
  d=cutree(d,6)
  dtMeta$HClustFull = "NA"
  dtMeta[names(d),"HClustFull"] = paste0("ProtPAM37.",as.integer(d))
  
  xHClustFull = topo.colors(length(unique(dtMeta$HClustFull)))
  names(xHClustFull) = unique(dtMeta$HClustFull)
  haT = HeatmapAnnotation(df = dtMeta[colnames(dtx), c( "PAM50","HClustFull", "ConsensusClust","HER.2","ER","PR","Histology")]  ,
                          col=list(
                            PAM50 = xcolorRanges$PAM50,
                            ConsensusClust = xcolorRanges$ConsensusClust,
                            ER=xcolNPE, 
                            PR=xcolNPR, 
                            'HER.2'=xcolNPH, 
                            Histology=xcolH,
                            HClustFull = xHClustFull
                          ))
  
  
  x=Heatmap(matrix = data.matrix(dtx[idx,]),  
            top_annotation = haT, 
            clustering_method_columns = "average" , 
            clustering_distance_columns  = "pearson" , 
            name = "Protein PAM37" , bottom_annotation = haB  ) 
  draw(x )
  
  ########################################################################################################
  ####  RNA PAM50
  ########################################################################################################
  
  
  rownames(dtRNA) = dtRNA$geneName
  dtRNA$geneName = NULL
  idx = rownames(dtRNA) %in% gnsPAM50
  dtxR=dtRNA[idx,colnames(dtProt)]
  
  d=as.dist(1-cor(dtxR, method="pearson"))
  d=hclust(d,method="average")
  d=cutree(d,6)
  
  dtMeta$HClustFull = "NA"
  dtMeta[names(d),"HClustFull"] = paste0("RNA.PAM50",as.integer(d))
  
  xHClustFull = topo.colors(length(unique(dtMeta$HClustFull)))
  names(xHClustFull) = unique(dtMeta$HClustFull)
  haT = HeatmapAnnotation(df = dtMeta[colnames(dtx), c( "PAM50","HClustFull", "ConsensusClust","HER.2","ER","PR","Histology")]  ,
                          col=list(
                            PAM50 = xcolorRanges$PAM50,
                            ConsensusClust = xcolorRanges$ConsensusClust,
                            ER=xcolNPE, 
                            PR=xcolNPR, 
                            'HER.2'=xcolNPH, 
                            Histology=xcolH,
                            HClustFull = xHClustFull
                          ))
  
  
  x=Heatmap(matrix = data.matrix(dtxR),  
            top_annotation = haT, 
            clustering_method_columns = "average" , 
            clustering_distance_columns  = "pearson" , 
            name = "mRNA PAM50" , bottom_annotation = haB ) 
  draw(x )
  
  ########################################################################################################
  ####  RNA PAM37
  ########################################################################################################
  
  idx = rownames(dtRNA) %in% gnsPAM50 & rownames(dtRNA) %in% rownames(dtProt)
  dtxR=dtRNA[idx,colnames(dtProt)]
  
  
  d=as.dist(1-cor(dtxR, method="pearson"))
  d=hclust(d,method="average")
  d=cutree(d,6)
  
  dtMeta$HClustFull = "NA"
  dtMeta[names(d),"HClustFull"] = paste0("RNA.PAM37",as.integer(d))
  
  xHClustFull = topo.colors(length(unique(dtMeta$HClustFull)))
  names(xHClustFull) = unique(dtMeta$HClustFull)
  haT = HeatmapAnnotation(df = dtMeta[colnames(dtx), c( "PAM50","HClustFull", "ConsensusClust","HER.2","ER","PR","Histology")]  ,
                          col=list(
                            PAM50 = xcolorRanges$PAM50,
                            ConsensusClust = xcolorRanges$ConsensusClust,
                            ER=xcolNPE, 
                            PR=xcolNPR, 
                            'HER.2'=xcolNPH, 
                            Histology=xcolH,
                            HClustFull = xHClustFull
                          ))
  
  
  x=Heatmap(matrix = data.matrix(dtRNA[idx,colnames(dtProt)]),  
            top_annotation = haT, 
            clustering_method_columns = "average" , 
            clustering_distance_columns  = "pearson" , 
            name = "mRNA PAM37"  , bottom_annotation = haB )  
  draw(x )  
  
  dev.off()
  
  
}

generateHClust = function (dtx, ctCluster){
  Xclu=hclust(as.dist(1 - cor(dtx,method = "pearson")), method = "average")
  Xclu=cutree(Xclu,k=ctCluster)
  xx = names(table(Xclu))[table(Xclu) < 3]
  for (i in 1:length(Xclu)){
    Xclu[i] = paste0("HClust",Xclu[i])
  }
  
  xcolorRanges$HClust = rainbow(length(unique(Xclu)))
  names(xcolorRanges$HClust) = sort(unique(Xclu))
  
  
  dtMeta$HClust = ""
  dtMeta[names(Xclu ), "HClust"] = Xclu
  
  cts=read.table("../DATA/geneCats.tab",sep="\t", header = T, stringsAsFactors = F)
  cts = cts[!is.na(cts$Id),]
  rownames(cts) = cts$Id 
  
  
  pdf(paste0("Figures/",SetName,"/Hclust_Long.pdf"),width = 15, height = as.integer(nrow(dtx)/8))
  haT = HeatmapAnnotation(df = dtMeta[colnames(dtx), c( "PAM50", "HClust", "ConsensusClust")]  ,
                          col=list(
                            PAM50 = xcolorRanges$PAM50,
                            HClust = xcolorRanges$HClust,
                            ConsensusClust = xcolorRanges$ConsensusClust
                          ))
  xdtx = cts[rownames(dtx),c("BCLandscape_Protein_Clustering","Fredlund_Classifier" ) ]
  xdtx$Fredlund_Classifier = NULL
  
  xdtx$BCLandscape_Protein_Clustering[is.na(xdtx$BCLandscape_Protein_Clustering)] = "No.enrichment"
  haR = rowAnnotation(df =  xdtx, 
                      col=list(
                        BCLandscape_Protein_Clustering=xcolorRanges$ProtColRange
                      )
                      
  )
  
  
  x=Heatmap(matrix = data.matrix(dtx),  top_annotation = haT, clustering_method_columns = "average" , 
            clustering_distance_columns  = "pearson" , name = "Protein"  ) + haR
  draw(x )
  dev.off()
  
  pdf(paste0("Figures/",SetName,"/Hclust_Short.pdf"),width = 15, height = 8)
  idxLabels=rownames(dtx) %in% c("TYMS","SOX10","CBS","CASC3","JMJD4","BRCA1","FOXA1","GATA3","GREB1","CCND1","CDK12",
                                 "GRB7","STARD3","RARA","IGF1R","KERA","AR","ERBB2","MIEN1","PGAP3","ESR1","ERBB4",
                                 "GINS1","CCNB1","MKI67","TOP2A","CDC20","TMEM26","NTRK2","NRIP2","TFPI","NRIP1","MX1",
                                 "PDK1","MET","SFRP1") 
  
  
  idxLabels= which(idxLabels)
  lbls=rownames(dtx)[idxLabels]
  
  
  
  ra = rowAnnotation(link = row_anno_link(at = idxLabels, labels = lbls),
                     width = unit(1, "cm") + max_text_width(lbls))
  
  haT = HeatmapAnnotation(df = dtMeta[colnames(dtx), c( "PAM50", "ConsensusClust","HER.2","ER","PR","Histology")]  ,
                          col=list(
                            PAM50 = xcolorRanges$PAM50,
                            ConsensusClust = xcolorRanges$ConsensusClust,
                            'ER'= xcolorRanges$xcolNPE, 'PR'= xcolorRanges$xcolNPR, 'HER.2'= xcolorRanges$xcolNPH, 'Histology'= xcolorRanges$xcolH
                          ))
  
  x=Heatmap(matrix = dtx,  top_annotation = haT, clustering_method_columns = "average" ,cluster_rows = T,
            clustering_distance_columns  = "pearson", show_row_names = F,name = "Protein"  ) + ra + haR 
  
  draw(x )
  dev.off()
  
  
  
  dtMeta <<- dtMeta
}


getTCGASubtype = function(fl){
  dtx=read.table(fl,sep="\t",header  = T,stringsAsFactors = F, row.names = 1)
  dtx[is.na(dtx)] = 0
  dtx=data.table(dtx)
  colnames(dtx)[1] = 'SubType'
  dtx=dtx[, lapply(.SD, mean), by=SubType]
  nms = dtx$SubType
  nms[nms == ""] = "NA"
  nms=str_replace_all(nms," ",".")
  nms=str_replace_all(nms,",","")
  dtx$SubType = NULL
  dtx=t(dtx)
  
  colnames(dtx) = paste0(str_replace( basename(fl),".tsv",""),".",nms)
  dtx
  
  
}


writeGML = function(dtx){
  #dtx=log2(dtProt[idx,SampsNonOutliers] )
  xgraph = getGraph(dtx, tPearsonNetSDMult)
  
  ## Generate average by ConsensusClass 
  for (cns in unique(dtMeta$ConsensusClust)){
    if (cns != ""){
      nms = rownames(dtMeta)[dtMeta$ConsensusClust == cns]
      if (length(nms) == 1){
        v= dtx[names(V(xgraph)), nms]
      }else{
        v= rowMeans( dtx[names(V(xgraph)), nms]) 
      }
      xgraph = set_vertex_attr(xgraph,  cns, V(xgraph),  v)
      
    }
  }
  
  ## Generate average by ConsensusClass 
  for (cns in unique(dtMeta$PAM50)){
    nms = dtMeta$id[dtMeta$PAM50 == cns]
    nms = nms[nms %in% colnames(dtx)]
    v= rowMeans( dtx[names(V(xgraph)), nms]) 
    xgraph = set_vertex_attr(xgraph,  paste0("PAM50_",cns), V(xgraph),  v)
  }
  
  cts=read.table("../DATA/geneCats.tab",sep="\t", header = T, stringsAsFactors = F)
  cts = cts[!is.na(cts$Id),]
  rownames(cts) = cts$Id 
  V(xgraph)$BCLandscape_Protein_Clustering = cts[names(V(xgraph)), "BCLandscape_Protein_Clustering"]
  V(xgraph)$Fredlund_Classifier = cts[names(V(xgraph)), "Fredlund_Classifier"]
  
  for (fl in c("../DATA/TCGA/BRCAUS_CNV_byPAM50.tsv","../DATA/TCGA/BRCAUS_Mut_byPAM50.tsv","../DATA/TCGA/OV_CNV_ByStage.tsv","../DATA/TCGA/PRAD_CNV.tsv")){
    dtx=getTCGASubtype(fl)
    for (x in colnames(dtx)){
      nms=names(V(xgraph))
      nms=nms[nms %in% rownames(dtx)]
      xgraph = set_vertex_attr(xgraph,  x, V(xgraph),  0)
      v=dtx[nms,x]
      xgraph = set_vertex_attr(xgraph,  x, V(xgraph)[nms],  v)
    }
  }
  
  
  selectedGns  = readLines(con=paste0("Figures/",SetName,"/SelectedGeneList.lst"))
  V(xgraph)$SelectedGene = "NO"
  V(xgraph)$SelectedGene[names(V(xgraph)) %in% selectedGns] = "YES"
  
  
  flgml=paste0("Figures/",SetName,"/CoreNet.gml")
  igraph::write_graph(xgraph, file = flgml, format = "gml")
  
}

### Not working
getGeneGroupSummarizedData = function(dtx,flGMT){
  flGMT="/Z/aleferna/GeneLists/MSigDB/c2.cp.reactome.v5.2.symbols.gmt"
  gmt=GSA::GSA.read.gmt(flGMT)
  
  dtx=abs(log2(dtProt))
  res=foreach (gs = 1:length(gmt$geneset.names), .combine=rbind) %dopar% {
    
    gsName = gmt$geneset.names[gs]
    gns = gmt$genesets[[gs]]
    idx = gns %in% rownames(dtx)
    if (sum(idx)>5){
      gns=gns[idx]
      x=data.frame(t(colMeans(dtx[gns,])))
      x$GeneSetName=gsName
      x  
    }else{
      NULL
    }
  }
  rownames(res) = res$GeneSetName
  res$GeneSetName  = NULL
  x= apply(res, 1, sd)
  
  res = res[order(x, decreasing = T) < 51, ] #top 50
  
  
  y=scale(res)
  
  haT = HeatmapAnnotation(df = dtMeta[colnames(res), c("ProtSub9506", "PAM50")]  ,  
                          col=list(PAM50 = xcolorRanges$PAM50, ProtSub9506 = xcolorRanges$ProtSub9506 ))
  res[res>1] = 1
  Heatmap(res, top_annotation = haT)
  res
}


AnnotateGeneClusters = function(dtx){
  
  dir.create(paste0( "Figures/",SetName,"/GeneEnrich/") , showWarnings = F)
  
  x=foreach (i = 3:20, .combine = rbind, .errorhandling = 'remove' ) %do% {
    d = as.dist(1 - cor(t(dtx),method = "pearson"))
    Xclu=hclust(d, method = "average")
    Xclu=cutree(Xclu,k=i)
    a=cluster::silhouette(Xclu , d)
    mean(a[,3])
  }
  plot(x)
  Xclu=hclust(as.dist(1 - cor(t(dtx),method = "pearson")), method = "average")
  
  ### Educated guess
  Xclu=cutree(Xclu,k=11)
  writeLines( paste(names(Xclu), paste0("GRP",Xclu), sep="\t") , con=paste0("Figures/",SetName,"/GeneEnrich/GeneGroups.tab" ))
  
  xx = names(table(Xclu))[table(Xclu) < 3]
  for (i in 1:length(Xclu)){
    if (Xclu[i] %in% xx){
      Xclu[i] = "NA"
    }else{
      Xclu[i] = paste0("GRP",Xclu[i])
    }
  }
  
  foreach  (flGmt = list.files("/Z/aleferna/GeneLists/MSigDB/", full.names = T)  ) %dopar%{
    gmt = GSA::GSA.read.gmt(flGmt)
    res=foreach (lstName = unique(Xclu), .combine=rbind) %do% {
      lst = names(Xclu)[Xclu == lstName]
      if (length(lst) > 5){
        foreach (i = 1:length(gmt$geneset.names), .combine=rbind) %do% {
          gslst = gmt$genesets[[i]]
          idx = gslst %in% rownames(dtx)
          if (sum(idx) > 5){
            gsn=gmt$geneset.names[i]
            gslst = gslst[idx]
            ctHits = sum(lst %in% gslst)
            lstHits = paste(lst[lst %in% gslst],collapse=",")
            ctRef = length(gslst)
            p = -10 * log10 ( 1.0-phyper(ctHits, ctRef, nrow(dtx)-ctRef, length(lst))) 
            data.frame(lstName,gsn,p,ctHits,lstHits, stringsAsFactors = F)  
          }
        }
      }
    }
    if (!is.null(res)){
      res$p[is.infinite(res$p)] = 200
      res$p = as.integer(res$p)
      res = reshape::cast(res, gsn ~ lstName, value = "p" )
      write.table(res,file=paste0( "Figures/",SetName,"/GeneEnrich/" ,basename(flGmt),".tab"), sep="\t", row.names = F, quote = F)
      rownames(res) = res$gsn
      res$gsn = NULL
      idx = apply(res,1,max) > 30
      if (sum(idx) > 2){
        res = data.frame(res[idx,])
        pdf(paste0( "Figures/",SetName,"/GeneEnrich/" ,basename(flGmt),".pdf") , height=max(15,as.integer(nrow(res)/4)) , width=15 )
        x=Heatmap(data.matrix(res), row_names_max_width = unit(100, "mm")   )
        draw(x, heatmap_legend_side = "bottom" )
        dev.off()
      }
    }
    0
  }
}

plot2DDendrogram = function(flOut, dtx){
  
  smps = colnames(dtx)
  d=as.dist( 1-cor( data.matrix(dtx),method = "pearson" )  )
  dtx = hclust(d, method = "average")
  dtx=as.dendrogram(dtx)
  
  ul <- function(x, parent ) {
    if (is.list(x)) {
      i=0
      foreach (y = x, .combine=rbind) %do% {
        i = i + 1 
        child=paste0(parent,".",i)
        rbind(ul(y, child ),c(parent,child,attr(y,"height")))
      }
    }else{
      name = attr(x, "label")
      c(parent,name,0)
    }
    
  }
  dtx=ul(dtx,"root")
  
  edges=data.frame(dtx,stringsAsFactors = F)
  colnames(edges)=c("source","target","weight")
  idxBad = grepl(edges$target, pattern = "OSL")
  Bad=edges[idxBad,]
  edges=edges[!idxBad,]
  for (i in 1:nrow(Bad)){
    idx = edges$target == Bad$source[i] 
    edges$target[idx] = Bad$target[i] 
  } 
  
  
  nodes=data.frame(name=sort(unique(c(edges$source,edges$target))),subtype="" , stringsAsFactors = F)
  rownames(nodes) = nodes$name
  nodes[smps,"subtype"] = dtMeta[smps,"ConsensusClust"]
  nodes$subtype[nodes$subtype == ""] = "Weak"
  idx=!nodes$name %in% smps
  nodes$subtype[idx] = "Group"
  nodes[!idx,"sz"] = 100
  nodes[idx,"sz"] = 2
  
  nodes$id = 0:(nrow(nodes)-1)
  rownames(nodes) = nodes$name
  edges$source = nodes[edges$source,"id"]
  edges$target = nodes[edges$target,"id"]
  edges$weight=10*as.double( edges$weight)+0.1
  
  
  cls=xcolorRanges[["ConsensusClust"]]
  xvals=paste0(paste0('"',c(labels(cls),"Group","Weak"),'"'), collapse=',')
  xcols=paste0(paste0('"',c(cls,"#808080","#808000"),'"'), collapse=',')
  xscale = paste0('d3.scaleOrdinal().domain([',xvals,']).range([',xcols,"])")
  xscale = 'd3.scaleOrdinal()
  .domain(["CONS1","CONS2","CONS3","CONS4","CONS5","CONS6","Group","Weak"])
  .range(["#FF0000","#FFFF00","#00FF00","#00FFFF","#0000FF","#FF00FF","#808080","#404040"]);'
  
  a=forceNetwork(Links=edges, Nodes=nodes, NodeID="name", Group="subtype", 
                 legend = T, 
                 zoom = T,  
                 Value = "weight",
                 linkDistance = JS('function(d){return d.value }') ,
                 Nodesize = "sz",
                 fontSize = 20,
                 opacity = 1 ,
                 colourScale  = xscale
  )
  
  saveNetwork(a, flOut, selfcontained = TRUE)
  
  
}




ProcessData = function(tPearsonNetSDMult,NTFilt){
  dir.create(paste0( "Figures/",SetName),showWarnings = T)
  file.copy(from = "Generate_Clustering_and_Network.R", to=paste0("Figures/",SetName))
  
  gnsHighVar = reportHighExpVarGenes()
  
  SampsNonOutliers = reportNonOutliers(log2(dtProt[gnsHighVar,]), ctClusters)
  writeLines(text =  SampsNonOutliers , paste0("Figures/",SetName,"/SampsNonOutliers.lst"))
  
  #Remove genes that no longer have std dev
  X=log2( apply( dtProt[gnsHighVar,], 1, sd))
  Y=log2( apply( dtProt[gnsHighVar,SampsNonOutliers], 1, sd))
  gnsHighVar = names(Y)[Y > min(X)]
  writeLines(text =  gnsHighVar , paste0("Figures/",SetName,"/gnsHighVar.lst"))
  
  
  
  
  
  gnsTumor = reportTumorGenes(log2(dtProt))
  writeLines(text =  gnsTumor , paste0("Figures/",SetName,"/gnsTumor.lst"))
  
  
  
  
  idx = rownames(dtProt) %in% gnsHighVar &
    rownames(dtProt) %in% gnsTumor 
  
  writeLines(text =  rownames(dtProt)[idx] , paste0("Figures/",SetName,"/SelectedGeneList.lst"))
  
  
}

GenerateFigures = function(){
  
  gnsHighVar = readLines(paste0("Figures/",SetName,"/gnsHighVar.lst"))
  gnsSel = readLines(paste0("Figures/",SetName,"/SelectedGeneList.lst"))
  idx = rownames(dtProt) %in% gnsSel
  SampsNonOutliers = readLines(paste0("Figures/",SetName,"/SampsNonOutliers.lst"))
  
  generateConsClustering(log2(dtProt[idx,SampsNonOutliers]), ctClusters)
  generateHClust(log2(dtProt[idx,SampsNonOutliers]), ctClusters)
  
  write.table(dtMeta[dtMeta$ProtSub9506 != "",], file=paste0("Figures/",SetName,"/ClusteringAndMeta.tab"), sep="\t", quote = F, row.names = F)
  idx = rownames(dtProt) %in% gnsHighVar
  
  generateGeneListHeatmaps()
  
  idx = rownames(dtProt) %in% gnsSel
  plot2DDendrogram(paste0("Figures/",SetName,"/Dendro2D.VIP.html"),log2(dtProt[idx,SampsNonOutliers]))
  plot2DDendrogram(paste0("Figures/",SetName,"/Dendro2D.ALL.html"),log2(dtProt[idx,]))
  
  ## Write GML 
  idx = rownames(dtProt) %in% gnsSel
  
  writeGML(log2(dtProt[idx,SampsNonOutliers] ))
  
  AnnotateGeneClusters(log2(dtProt[idx,SampsNonOutliers]))
}




load("../DATA/Prot.Rdata")
rownames(dtProt) = dtProt$Gene.Symbol
dtProt$Gene.Symbol = NULL
dtProt = dtProt[sample(1:nrow(dtProt) ),]
dtProt = dtProt[,sample(1:45)]


SetName = ""
dtMeta=read.table(file="../DATA/meta.tab",sep = "\t",header = T,stringsAsFactors = F)
rownames(dtMeta) = dtMeta$id
dtMeta$HER.2[dtMeta$HER.2 == "0"] = NA
dtMeta$HER.2[is.na(dtMeta$HER.2)] = "NA"
dtMeta$ER[is.na(dtMeta$ER)] = "NA"
dtMeta$PR[is.na(dtMeta$PR)] = "NA"


load("../DATA/colorRanges.Rdata")
xcolorRanges$ProtColRange= c("Erythrocyte"="#e31a1c","Plasma"="#ff7f00","Adipose"="#b15928",
                             "ECM.1"="#6a3d9a","ECM.2"="#cab2d6","Luminal"="#1f78b4",
                             "Transcription"="#a6cee3","Immune"="#b2df8a","No enrichment"="#d8b365",
                             "Golgi"="#fdbf6f","Basal"="#bc2025","Mitochondrion"="#33a02c",
                             "Proliferation"="#ffff99","HER2.amplicon"="#fb9a99","NA"="#e5e1e8",
                             "No.enrichment"="#e5e1e8")


xcolorRanges$xcolNPE= c("grey", "blue", "white" )
names(xcolorRanges$xcolNPE) =c("neg", "pos","NA")

xcolorRanges$xcolNPR= c("grey", "orange", "white" )
names(xcolorRanges$xcolNPR) =c("neg", "pos","NA")

xcolorRanges$xcolNPH= c("grey", xcolorRanges$PAM50["Her2"], "white" )
names(xcolorRanges$xcolNPH) =c("neg", "pos","NA")

xcolorRanges$xcolH= c("#e66101","#fdb863","#f7f7f7","#b2abd2","#5e3c99")
names(xcolorRanges$xcolH) =c("Ductal" ,  "DCIS"  ,   "Mucinous" ,"Lobular" , "Mixed" )




NTFilt=T
tPearsonNetSDMult=1.96
ctClusters=6

SetName = "HighVar_NT_Outlier_ct6_MAR3"
#ProcessData(tPearsonNetSDMult,NTFilt)  
GenerateFigures()



#  Sent to anna lise
# load("/Z/aleferna/Breast_Cancer_Landscape/DATA/geneListsHJ.RData")
# x=readLines("/Z/aleferna/Breast_Cancer_Landscape/CoreNet_and_Clustering/Figures/HighVar_NT_Outlier_ct6_FEB27/SelectedGeneList.lst")
# geneLists$CoreNet[!geneLists$GeneName %in% x] = F
# write.table(geneLists , file="/Z/aleferna/Breast_Cancer_Landscape/CoreNet_and_Clustering/Figures/HighVar_NT_Outlier_ct6_FEB27/GeneMtx.tab", sep="\t",quote = F, row.names = F)
# geneLists$GeneName[geneLists$ECM.2]

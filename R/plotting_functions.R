
require(igraph)
require(scales)
require(RColorBrewer)
require(MASS)

plot_feat <- function(x,red="umap",feat=NULL,label=NULL,assay="RNA",pch=16,
                      bg=NA,font.labels=1,cex.labels=1,cex=.3,dims=c(1,2),mins=0,
                      add_graph=NULL,percent_connections=1,nbin=400,n=10,main=NULL,
                      col=c("grey90","grey80","grey70","navy","black"),...){
  fn <- feat

  if(feat %in% rownames(x@assays[[assay]]@data) ){
    feat <- x@assays[[assay]]@data[feat,]
  } else if(feat %in% colnames(x@meta.data) ) { feat <- as.numeric(x@meta.data[,feat])
  } else { message("Feature or metadata not found!!")
    feat <- rep(0,ncol(x))
  }

  if(is.null(mins)){mins <- min(feat,na.rm = T)}
  if( sum(is.na(feat)) > 0 ){ feat[is.na(feat)] <- 0 }
  
  if(max(feat) != 0){
    feat <- (feat - mins)/ ( sort(feat,T,na.last = T)[ min(10,sum(feat!=0,na.rm = T))  ] - mins)
    feat[feat > 1] <- 1}
  o <- order(feat,na.last = T)

  pal <- c( col[1],colorRampPalette(col[-1])(10))[round(feat*9)+1][o]
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)

  #creates plot
  plot( x@reductions[[red]]@cell.embeddings[o,dims] , type="n", xaxt="n",yaxt="n", main=ifelse( !is.null(main), main, paste0(assay,"_",fn)),...)

  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(x@graphs)  ){
    add_graph(x,red=x@reductions[[red]]@cell.embeddings[o,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}

  #compute density and centroids
  if(!is.null(label) ){
    labels <- factor(as.character(x@meta.data[[as.character(label)]]))
    if( !is.na(sum(as.numeric(levels(feat)))) ){
      labels <- factor(as.numeric(as.character(x@meta.data[[as.character(label)]])))}
    centroids <-  sapply( as.character(levels(labels)) ,
                          red=x@reductions[[red]]@cell.embeddings[,dims],
                          cl1=labels,
                          function(jj,red,cl1) { pmean(red[cl1==jj,])  })
  }

  #adds points
  points(x@reductions[[red]]@cell.embeddings[o,dims],pch=pch,cex=cex,bg=bg, col=paste0(pal) )
  options(warn=0)

  #adds labels
  if(!is.null(label)){
    # points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = levels(labels),cex=cex.labels,font=font.labels,xpd=T)
    par(xpd=F)
  }
}




plot_meta <- function(x,red="umap",feat=NULL,pch=16,cex=.3,label=F,dims=c(1,2),font.labels=1,cex.labels=1, col = c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Pastel2") ),
                      add_graph=NULL,percent_connections=1,nbin=400,add_lines=F,main=NULL,...){
  fn <- feat
  feat <- factor(as.character(x@meta.data[[feat]]))
  if( !is.na(sum(as.numeric(levels(feat)))) ){
    feat <- factor(as.numeric(as.character(x@meta.data[[fn]])))}
  try(col <- colorRampPalette(col)( max( length(col) , length(unique(feat))) ))
  col <- col[feat]
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)

  #creates plot
  plot( x@reductions[[red]]@cell.embeddings[,dims], type="n", xaxt="n",yaxt="n",main=ifelse( !is.null(main), main, paste0(fn)),...)

  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(x@graphs)  ){
    add_graph(x,red=x@reductions[[red]]@cell.embeddings[,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}

  #compute density and centroids
  if(label | add_lines){
    centroids <-  sapply( as.character(levels(feat)) , 
                          reds=as.data.frame(x@reductions[[red]]@cell.embeddings[,dims]), 
                          cl1=feat, function(jj,reds,cl1) { pmean(reds[cl1==jj,])  })
    }

  if(add_lines){
    add_centroid_lines(x@reductions[[red]]@cell.embeddings[,dims],feat,col,centroids)}

  #adds points
  points(x@reductions[[red]]@cell.embeddings[,dims],pch=pch,cex=cex,bg=paste0(col,90), col=col )
  options(warn=0)

  #adds labels
  if(label){
    # points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = levels(feat),cex=cex.labels,font=font.labels,xpd=T)
    par(xpd=F)
    }
}





plot_meta2 <- function(data,red="umap",feat=NULL,pch=16,cex=.3,label=F,dims=c(1,2),
                       font.labels=1,cex.labels=1, col = c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Pastel2") ),
                      add_graph=NULL,percent_connections=1,nbin=400,add_lines=F,main=NULL,gapx=.2,gapy=.2,
                      ncol=1,...){
  
  fn <- feat
  n <- length(fn)
  
  
  #creates plot
  plot( c( 0 , ncol ), c( 0 , ceiling(n/ncol) ), type="n", 
        frame=F,
        xaxs="i",yaxs="i",
        xaxt="n",yaxt="n",
        xlab="",ylab="",
        main="",...)
  
  rowmin <- ceiling(n/ncol) - 1
  colmin <- 0
  for(i in 1:n){
    feat <- fn[i]
    if( colmin >= ncol ){rowmin <- rowmin - 1}
    if( colmin >= ncol ){colmin <- 0}
    
    
    feat <- factor(as.character(data@meta.data[[fn[i]]]))
    if( !is.na(sum(as.numeric(levels(feat)))) ){
      feat <- factor(as.numeric(as.character(data@meta.data[[fn[i]]])))}
    try(col <- colorRampPalette(col)( max( length(col) , length(unique(feat))) ))
    col2 <- col[feat]
    #par(mar=c(1.5,1.5,1.5,1.5))
    options(warn=-1)
    
    
    x <- data@reductions[[red]]@cell.embeddings[,dims[1]]
    x <- (x - min(x,na.rm = T)) / (max(x,na.rm = T) - min(x,na.rm = T) )*(1-gapx)+gapx/2
    x <- x + colmin
    
    y <- data@reductions[[red]]@cell.embeddings[,dims[2]]
    y <- (y - min(y,na.rm = T)) / (max(y,na.rm = T) - min(y,na.rm = T))*(1-gapy)+gapy/2
    y <- y + rowmin
    
    #adds underlying graph
    if(!is.null(add_graph) ){ if( add_graph %in% names(data@graphs)  ){
      add_graph(data,red=cbind(x,y),graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}
    
    #compute density and centroids
    if(label | add_lines){
      centroids <-  sapply( as.character(levels(feat)) , 
                            reds=as.data.frame(cbind(x,y)), 
                            cl1=feat, function(jj,reds,cl1) { pmean(reds[cl1==jj,])  })
    }
    
    if(add_lines){
      add_centroid_lines(cbind(x,y),feat,col2,centroids)}
    
    #adds points
    points(cbind(x,y),pch=pch,cex=cex,bg=paste0(col2,90), col=col2 )
    options(warn=0)
    
    text(1+colmin,1+rowmin,labels = fn[i],adj=c(1,1),cex=cex.labels,font=font.labels,xpd=T)
    
    #adds labels
    if(label){
      # points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
      text(centroids[1,],centroids[2,],labels = levels(feat),cex=cex.labels,font=font.labels,xpd=T)
      par(xpd=F)
    }
    colmin <- colmin + 1
  }
}




add_graph <- function(a,red,graph="RNA_snn",percent_connections=1){
  g <- graph_from_adjacency_matrix(a@graphs[[graph]],weighted = T,diag = F)
  g <- simplify(g)
  eee <- as.data.frame(as_edgelist(g))
  sel <- (1:nrow(eee))[ 1:nrow(eee) %in% round(seq(1,nrow(eee),length.out = round(nrow(eee)*percent_connections)))]
  eee <- eee[sel,]
  eee$x0 <- red[ match(eee[,1],rownames(red)), 1 ]
  eee$x1 <- red[ match(eee[,2],rownames(red)), 1 ]
  eee$y0 <- red[ match(eee[,1],rownames(red)), 2 ]
  eee$y1 <- red[ match(eee[,2],rownames(red)), 2 ]
  apply(eee,1,function(x){ lines(x[3:4],x[5:6],col="grey80",lwd=.1) })
}




compute_density <- function(red,nbin=200){
  map <- MASS::kde2d(x = red[,1], y = red[,2],h = c(diff(apply(red, 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)) / 7.5),
                      n = c(nbin,nbin))
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  xbin <- cut(red[, 1], mkBreaks(map$x), labels = FALSE)
  ybin <- cut(red[, 2], mkBreaks(map$y), labels = FALSE)
  dens <- map$z[cbind(xbin, ybin)]
  dens[is.na(dens)] <- 0
  dens <- (dens-min(dens))/(max(dens)-min(dens))
  return(dens)
}


add_centroid_lines <- function(red,feat,pal,centroids){
  for(i in unique(feat)){
    df <- red[feat==i,]
    this_col <- pal[feat==i][1]
    c1 <- centroids[1,i]
    c2 <- centroids[2,i]
    apply(df,1,function(x){ lines(c(c1,x[1]),c(c2,x[2]),col=this_col,lwd=.2) })
  }

}



plot_gene_cloud <- function(TOM, gene_module, mm, main=NULL,...){
  gtemp <- igraph::graph_from_adjacency_matrix(TOM[names(gene_module[ gene_module %in% mm ]),names(gene_module[ gene_module %in% mm ])],weighted = T,diag = F)
  l <- layout_nicely(gtemp)
  plot(l,type="n",axes=F,frame=F,xlab="",ylab = "",main=ifelse(is.null(main), mm, main))
  text(l , labels=names(gene_module[ gene_module %in% mm ]) , ...)
}



pmean <- function(x, steps=30, minf=.6){
  f <- 1
  if(is.vector(x)){
    mm <- mean(x)
    for(i in 1:steps){
      dx <-  x - mm
      ind <- order( abs(dx) ) [ 1:ceiling(length(x) * max(minf,f) )]
      mm <- mean(x[ind])
      f <- f - .01 }
  } else {
    mm <- colMeans(x)
    for(i in 1:steps){
      dx <- ( colSums( (t(x)-mm )^2 ) ) ^ (1/2)
      ind <- order( abs(dx) ) [ 1:ceiling(nrow(x) * max(minf,f) )]
      mm <- colMeans( x[ ind , ] )
      f <- f - .01 }
   }
  return( mm ) }


lower_sd <- function(x, mm=mean(x)){
  sds <- sqrt( sum((x[ x < mm ] - mm)^2)/length(x) )
  return(sds)
}




plot_tree <- function( data, slingshot_curves, gene, rotate90=F, assay="RNA",edge.weights=F,
                       pal=c("grey90","grey70","blue3","navy"),minsize=.5,sizefactor=2,...){
  g <- igraph::graph_from_adjacency_matrix( curves@adjacency, mode = "undirected")
  
  temp <- factor(colSums(t(curves@clusterLabels) * as.numeric(colnames(curves@clusterLabels))))
  if( !is.na(sum(as.numeric(levels(temp)))) ){
    temp <- factor(as.numeric(as.character(temp))) }
  
  x1 <- rowsum(as.matrix(a@assays[[assay]]@data[gene,]), temp)
  x1 <- t(x1 / c(table(temp)))
  x1 <- t(apply(t(x1) , 2,function(i) (i-0)/(max(i)-0) ) )
  
  x2 <- rowsum(( as.matrix(data@assays[[assay]]@data[gene,]!=0) *1), temp)
  x2 <- t(x2 / c(table(temp)))
  
  l <-  -igraph::layout_as_tree(g, root = curves@slingParams$start.clus)[,2:1]
  if(rotate90){l <-  igraph::layout_as_tree(g, root = curves@slingParams$start.clus)}
  
  plot( l , type="n",frame=F,axes=F,xlab="",ylab="",main=gene,font.main=1 ,cex.main=1)
  
  eee <- as.data.frame(as_edgelist(g))
  eee$x0 <- l[ match(eee[,1],1:nrow(l)), 1 ]
  eee$x1 <- l[ match(eee[,2],1:nrow(l)), 1 ]
  eee$y0 <- l[ match(eee[,1],1:nrow(l)), 2 ]
  eee$y1 <- l[ match(eee[,2],1:nrow(l)), 2 ]
  if(edge.weights){
    eee$weight <- apply(eee,1,function(x){ curves@slingParams$dist[ x[1], x[2] ] })
    eee$weight <- 1 - (eee$weight - 0) / ( max(eee$weight) - 0)
  } else {eee$weight <- .5}
  
  apply(eee,1,function(x){ lines(x[3:4],x[5:6],col="black",lwd=(as.numeric(x[7])*2+.5) ) })
  
  points(l , pch=21,
         bg=c( "grey95",colorRampPalette(pal)(19))[c(t(x1) )*18+1 ],
         cex=(x2+minsize)*sizefactor )
}






graph_abstraction <- function( data , red="umap" , clustering , graph="SNN", cutoff=0.01){
  clustering_use <- factor(data@meta.data[,clustering])
  g <- graph_from_adjacency_matrix(data@graphs[[graph]],weighted = T,diag = F)
  g <- simplify(g)
  eee <- as.data.frame(as_edgelist(g))
  eee$g1 <- clustering_use[ match(eee[,1],colnames(data)) ]
  eee$g2 <- clustering_use[ match(eee[,2],colnames(data)) ]
  
  res <- data.frame()
  for(k in levels(clustering_use)){
    temp <- table(eee[eee$g1 == k,"g2"])
    # temp <- temp / table(clustering_use)
    temp <- temp/sum(temp)
    temp[temp < cutoff] <- 0
    res <- rbind(res, data.frame(s=k,p=temp))
  }
  
  centroids <-  t(sapply( as.character(unique(clustering_use)) ,
                          red=data@reductions[[red]]@cell.embeddings,
                          cl1=as.character(clustering_use),
                          function(jj,red,cl1) { pmean(red[cl1==jj,])  }))
  
  res$x0 <- centroids[ match(res[,1],rownames(centroids)), 1 ]
  res$x1 <- centroids[ match(res[,2],rownames(centroids)), 1 ]
  res$y0 <- centroids[ match(res[,1],rownames(centroids)), 2 ]
  res$y1 <- centroids[ match(res[,2],rownames(centroids)), 2 ]
  
  return(res)
}




empty_plot <- function(){
  plot( c(0,1),c(0,1) , type="n", axes=F,
        main="",frame=F,xlab="",ylab="")
}









require(igraph)
require(scales)
require(RColorBrewer)



#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname plot_feat
plot_feat <- function(x,red="umap",feat=NULL,label=NULL,assay="RNA",pch=16,
                      bg=NA,font.labels=1,cex.labels=1,cex=.3,dims=c(1,2),mins=NULL,
                      add_graph=NULL,percent_connections=1,nbin=400,n=10,main=NULL,maxs=NULL,
                      col=c("grey90","grey80","grey60","navy","black"),add=F,frame=F,font.main=1,add_legend=T,cex.main=1,...){
  fn <- feat
  if(is(x,"Seurat")){
    red <- x@reductions[[red]]@cell.embeddings
    if(feat %in% rownames(x@assays[[assay]]@data) ){
      feat <- x@assays[[assay]]@data[feat,]
    } else if(feat %in% colnames(x@meta.data) ) { feat <- as.numeric(x@meta.data[,feat])
    } else { message("Feature or metadata not found!!")
      feat <- rep(0,ncol(x))
    }
  } else {
    feat <- x[feat,]
  }

  if(is.null(mins)){mins <- min(c(feat,0),na.rm = T)}
  if(is.null(maxs)){
    maxs <- quantile(feat,0.99,na.rm = T)
    if(maxs==0){maxs <- max(feat,na.rm = T)}
  }
  
  if( sum(is.na(feat)) > 0 ){ feat[is.na(feat)] <- 0 }
  
  if(max(feat,na.rm = T) != 0){
    feat <- (feat - mins) / ( maxs - mins)
    feat[feat > 1] <- 1}
  o <- order(feat,na.last = T)

  pal <- c( col[1],colorRampPalette(col[-1])(99))[round(feat*98)+1][o]
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)

  #creates plot
  if(!add){
  plot( red[o,dims] ,xlab="",ylab="", type="n", xaxt="n",yaxt="n", frame=frame, font.main=font.main,cex.main=cex.main, main=ifelse( !is.null(main), main, paste0(fn)),...)
}
  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(x@graphs)  ){
    add_graph(x,red=red[o,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}

  #compute density and centroids
  if(!is.null(label) ){
    labels <- factor(as.character(x@meta.data[[as.character(label)]]))
    if( !is.na(sum(as.numeric(levels(feat)))) ){
      labels <- factor(as.numeric(as.character(x@meta.data[[as.character(label)]])))}
    centroids <-  sapply( as.character(levels(labels)) ,
                          red=red[,dims],
                          cl1=labels,
                          function(jj,red,cl1) { pmean(red[cl1==jj,])  })
  }

  #adds points
  points(red[o,dims],pch=pch,cex=cex,bg=bg, col=paste0(pal) )
  options(warn=0)

  #adds labels
  if(!is.null(label)){
    # points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = levels(labels),cex=cex.labels,font=font.labels,xpd=T)
    par(xpd=F)
  }
  
  add_corner_axis(xlab=colnames(red[,dims])[1],ylab=colnames(red[,dims])[2])
  
  if(add_legend){
    add_scale_legend(labels = c("min","max"),
                     pal = c( col[1],colorRampPalette(col[-1])(98)))
  }
}




#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname plot_meta
plot_meta <- function(x,red="umap",feat=NULL,pch=16,cex=.3,label=F,dims=c(1,2),font.labels=1,cex.labels=1, 
                      col = c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Pastel2") ),
                      add_graph=NULL,percent_connections=1,nbin=400,add_lines=F,main=NULL,add=F,frame=F,font.main=1,cex.main=1,...){
  fn <- feat
    if(is(x,"Seurat")){
      red <- x@reductions[[red]]@cell.embeddings
    feat <- factor(as.character(x@meta.data[[feat]]))
    if( !is.na(sum(as.numeric(levels(feat)))) ){
      feat <- factor(as.numeric(as.character(x@meta.data[[fn]])))}
    try(col <- colorRampPalette(col)( max( length(col) , length(unique(feat))) ))
    cols <- col[feat]
  } else {
    feat <- x[feat,]
    try(col <- colorRampPalette(col)( max( length(col) , length(unique(feat))) ))
    cols <- col[feat]
  }
    
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)

  #creates plot
  if(!add){
  plot( red[,dims],xlab="",ylab="", type="n", xaxt="n",yaxt="n",font.main=font.main,cex.main=cex.main,main=ifelse( !is.null(main), main, paste0(fn)),frame=frame,...)
  }

  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(x@graphs)  ){
    add_graph(x,red=red[,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}

  #compute density and centroids
  if(label | add_lines){
    centroids <-  sapply( as.character(levels(feat)) , 
                          reds=as.data.frame(red[,dims]), 
                          cl1=feat, function(jj,reds,cl1) { pmean(reds[cl1==jj,])  })
    }

  if(add_lines){
    add_centroid_lines(red[,dims],feat,cols,centroids)}

  #adds points
  points(red[,dims],pch=pch,cex=cex,bg=paste0(col,90), col=cols )
  options(warn=0)

  #adds labels
  if(label){
    # points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = levels(feat),cex=cex.labels,font=font.labels,xpd=T)
    par(xpd=F)
  }
  
  add_corner_axis(xlab=colnames(red[,dims])[1],ylab=colnames(red[,dims])[2])
  
  legend(par("usr")[2],par("usr")[4],legend = levels(feat),
         pch=16,col=col, bty = "n",
         cex = 1, pt.cex = 1,xjust = 0,yjust = 1,
         title.adj = 1,xpd=T,y.intersp = .7)
  
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





#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname add_graph
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




# compute_density <- function(red,nbin=200){
#   map <- MASS::kde2d(x = red[,1], y = red[,2],h = c(diff(apply(red, 2, quantile, probs=c(0.05, 0.95), na.rm=TRUE)) / 7.5),
#                       n = c(nbin,nbin))
#   mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
#   xbin <- cut(red[, 1], mkBreaks(map$x), labels = FALSE)
#   ybin <- cut(red[, 2], mkBreaks(map$y), labels = FALSE)
#   dens <- map$z[cbind(xbin, ybin)]
#   dens[is.na(dens)] <- 0
#   dens <- (dens-min(dens))/(max(dens)-min(dens))
#   return(dens)
# }



#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname add_centroid_lines
add_centroid_lines <- function(red,feat,pal,centroids){
  for(i in unique(feat)){
    df <- red[feat==i,]
    this_col <- pal[feat==i][1]
    c1 <- centroids[1,i]
    c2 <- centroids[2,i]
    apply(df,1,function(x){ lines(c(c1,x[1]),c(c2,x[2]),col=this_col,lwd=.2) })
  }

}



#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname plot_gene_cloud
plot_gene_cloud <- function(TOM, gene_module, mm, main=NULL,...){
  gtemp <- igraph::graph_from_adjacency_matrix(TOM[names(gene_module[ gene_module %in% mm ]),names(gene_module[ gene_module %in% mm ])],weighted = T,diag = F)
  l <- layout_nicely(gtemp)
  plot(l,type="n",axes=F,frame=F,xlab="",ylab = "",main=ifelse(is.null(main), mm, main))
  text(l , labels=names(gene_module[ gene_module %in% mm ]) , ...)
}



#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname pmean
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




#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname plot_tree
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







#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname graph_abstraction
graph_abstraction <- function( data , red="umap" , clustering , graph="SNN", cutoff=0, nn=0){
  
  clustering_use <- factor(data@meta.data[,clustering])
  mm <- model.matrix( ~ 0 + clustering_use )
  colnames(mm) <- levels(clustering_use)
  
  # res <- data@graphs[[graph]]
  # # res <- res %*% Matrix::t(res)
  # res <- res %*% mm
  # res <- as.matrix(Matrix::t(res) %*% mm)
  # res <- res / c(table(clustering_use))
  # res <- t(res) / c(table(clustering_use))
  # # res[res < cutoff ] <- 0
  # 
  # res <- data.frame(s=c(sapply(colnames(res),res=res,function(x,res) rep(x,nrow(res)))),
  #                   p.Var=rep( rownames(res) , ncol(res)),
  #                   p.Freq= c(res) )
  g <- graph_from_adjacency_matrix(data@graphs[[graph]],weighted = T,diag = F)
  g <- simplify(g)
  eee <- as.data.frame(as_edgelist(g))
  eee$g1 <- clustering_use[ match(eee[,1],colnames(data)) ]
  eee$g2 <- clustering_use[ match(eee[,2],colnames(data)) ]
  # 
  res <- data.frame()
  for(k in levels(clustering_use)){
    temp <- table(eee[eee$g1 == k,"g2"])
    temp <- temp / table(clustering_use)
    temp <- temp / sum(temp)
    if(nn > 0){
      nn2 <- sort(temp,decreasing = T)[nn]
      temp[temp < nn2] <- 0
    }
    temp[temp < cutoff] <- 0
    res <- rbind(res, data.frame(s=k,p=temp))
  }
  
  centroids <-  t(sapply( as.character(unique(clustering_use)) ,
                          red=data@reductions[[red]]@cell.embeddings,
                          cl1=as.character(clustering_use),
                          function(jj,red,cl1) { pmean(red[cl1==jj,])  }))
  
  res$x0 <- as.numeric(centroids[ match(res[,1],rownames(centroids)), 1 ])
  res$x1 <- as.numeric(centroids[ match(res[,2],rownames(centroids)), 1 ])
  res$y0 <- as.numeric(centroids[ match(res[,1],rownames(centroids)), 2 ])
  res$y1 <- as.numeric(centroids[ match(res[,2],rownames(centroids)), 2 ])
  
  return(res)
}




#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname empty_plot
empty_plot <- function(...,main="",frame=F,xlab="",ylab="",cex.main=1,font.main=1){
  plot( c(0,1),c(0,1) , type="n", axes=F,
        main=main,frame=frame,xlab=xlab,ylab=ylab,cex.main=cex.main,font.main=font.main,...)
}






#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname plot_spatial_feat
plot_spatial_feat <- function(x,red="slice1",feat=NULL,res="lowres",label=NULL,assay="Spatial",plot_tissue=T,rescale=T,transparency="",pch=16,
                              bg=NA,font.labels=1,cex.labels=1,cex=1,mins=NULL,
                              add_graph=NULL,percent_connections=1,nbin=400,n=10,main=NULL,maxs=NULL,
                              col=c("grey95","grey70","navy","black"),...){
  fn <- feat
  
  if(feat %in% rownames(x@assays[[assay]]@data) ){
    feat <- x@assays[[assay]]@data[feat,]
  } else if(feat %in% colnames(x@meta.data) ) { feat <- as.numeric(x@meta.data[,feat])
  } else { message("Feature or metadata not found!!")
    feat <- rep(0,ncol(x))
  }
  
  if(is.null(mins)){mins <- 0}
  if(is.null(maxs)){
    maxs <- sort(feat,decreasing = T)[min(5,length(feat[feat!=0]))]
    
    if(maxs==0){ maxs <- max(feat,na.rm = T) }
    }
  
  if( sum(is.na(feat)) > 0 ){ feat[is.na(feat)] <- mins }
  
  if(max(feat,na.rm = T) != 0){
    if(rescale){
      feat <- (feat - mins) / ( maxs - mins)
      feat[is.na(feat)] <- 0
      feat[feat > 1] <- 1
    } else {
      feat <- feat / max(x@assays[[assay]]@data)
    }
  }
  o <- order(feat,na.last = T)
  
  pal <- paste0(c( colorRampPalette(col[1])(1),colorRampPalette(col[-1])(99))[round(feat*99)+1][o], transparency )
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)
  
  #creates plot
  coo <- x@images$slice1@coordinates
  coo$imagecol <- coo$imagecol*x@images$slice1@scale.factors$hires
  coo$imagerow <- dim(x@images$slice1@image)[1] - coo$imagerow*x@images$slice1@scale.factors$hires
  spot_col_lims <- range(coo$imagecol)
  spot_row_lims <- range(coo$imagerow)
  empty_plot(xlim=spot_col_lims+c(-5,5), ylim=spot_row_lims+c(-5,5),frame=F,yaxs="i",xaxs="i",
             main=ifelse( !is.null(main), main, paste0(assay,"_",fn)),...)
  if(plot_tissue){
    rasterImage(x@images$slice1@image,1,1,dim(x@images$slice1@image)[2],dim(x@images$slice1@image)[1])
  }
  # plot( x@reductions[[red]]@cell.embeddings[o,dims] , type="n", xaxt="n",yaxt="n", main=ifelse( !is.null(main), main, paste0(assay,"_",fn)),...)
  
  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(x@graphs)  ){
    add_graph(x,red=cbind(coo$imagecol , coo$imagerow)[o,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}
  
  #compute density and centroids
  if(!is.null(label) ){
    labels <- factor(as.character(x@meta.data[[as.character(label)]]))
    if( !is.na(sum(as.numeric(levels(feat)))) ){
      labels <- factor(as.numeric(as.character(x@meta.data[[as.character(label)]])))}
    centroids <-  sapply( as.character(levels(labels)) ,
                          red=cbind(coo$imagecol , coo$imagerow)[,dims],
                          cl1=labels,
                          function(jj,red,cl1) { pmean(red[cl1==jj,])  })
  }
  
  #adds points
  points(cbind(coo$imagecol , coo$imagerow)[o,],pch=pch,cex=cex,bg=bg, col=pal )
  options(warn=0)
  
  #adds labels
  if(!is.null(label)){
    # points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = levels(labels),cex=cex.labels,font=font.labels,xpd=T)
    par(xpd=F)
  }
}




#' @title AAA
#' @description AAA
#' @details AAA
#' @export
#' @rdname plot_spatial_meta
plot_spatial_meta <- function(x,red="slice1",feat=NULL,assay="Spatial",rescale=T,plot_tissue=T,pch=16,cex=1,label=F,dims=c(1,2),font.labels=1,cex.labels=1, bg=NA,
                              col = c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2"),RColorBrewer::brewer.pal(8,"Accent"),RColorBrewer::brewer.pal(9,"Pastel1"),RColorBrewer::brewer.pal(8,"Pastel2") ),
                      add_graph=NULL,percent_connections=1,nbin=400,add_lines=F,main=NULL,transparency="",...){
  fn <- feat
  feat <- factor(as.character(x@meta.data[[feat]]))
  if( !is.na(sum(as.numeric(levels(feat)))) ){
    feat <- factor(as.numeric(as.character(x@meta.data[[fn]])))}
  try(col <- paste0(colorRampPalette(col)( max( length(col) , length(levels(feat))) ),transparency) )
  col <- col[feat]
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)
  
  #creates plot
  coo <- x@images$slice1@coordinates
  coo$imagecol <- coo$imagecol*x@images$slice1@scale.factors$lowres
  coo$imagerow <- dim(x@images$slice1@image)[1] - coo$imagerow*x@images$slice1@scale.factors$lowres
  spot_col_lims <- range(coo$imagecol)
  spot_row_lims <- range(coo$imagerow)
  empty_plot(xlim=spot_col_lims+c(-5,5), ylim=spot_row_lims+c(-5,5),frame=F,yaxs="i",xaxs="i",
             main=ifelse( !is.null(main), main, paste0(assay,"_",fn)),...)
  if(plot_tissue){
    rasterImage(x@images$slice1@image,1,1,dim(x@images$slice1@image)[2],dim(x@images$slice1@image)[1])
  }
  
  
  
  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(x@graphs)  ){
    add_graph(x,red=cbind(coo$imagecol , coo$imagerow)[,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}
  
  #compute density and centroids
  if(label | add_lines){
    centroids <-  sapply( as.character(levels(feat)) , 
                          reds=as.data.frame(cbind(coo$imagecol , coo$imagerow)[,dims]), 
                          cl1=feat, function(jj,reds,cl1) { pmean(reds[cl1==jj,])  })
  }
  
  if(add_lines){
    add_centroid_lines(cbind(coo$imagecol , coo$imagerow)[,dims],feat,col,centroids)}
 
  #adds points
  points(cbind(coo$imagecol , coo$imagerow)[,dims],pch=pch,cex=cex,bg=bg, col=col )
  options(warn=0)
  
  #adds labels
  if(label){
    # points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = levels(feat),cex=cex.labels,font=font.labels,xpd=T)
    par(xpd=F)
  }
}


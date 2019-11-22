
require(igraph)
require(scales)
require(RColorBrewer)
require(MASS)
require(spatstat)

plot_feat <- function(x,red="umap",feat=NULL,label=NULL,assay="RNA",pch=16,bg=NA,cex=.3,dims=c(1,2),mins=0,add_graph=NULL,percent_connections=1,nbin=400,main=NULL,...){
  fn <- feat
  
  if(feat %in% rownames(x@assays[[assay]]@data) ){
    feat <- x@assays[[assay]]@data[feat,]
  } else if(feat %in% colnames(x@meta.data) ) { feat <- x@meta.data[,feat]  
  } else { message("Feature or metadata not found!!") }
  
  if(is.null(mins)){mins <- min(feat)}
  
  feat <- (feat - mins)/ ( sort(feat,T)[ min(10,sum(feat!=0))  ] - mins)
  feat[feat > 1] <- 1
  o <- order(feat)
  
  pal <- colorRampPalette(c("grey90","navy","blue3"))(10)[round(feat*9)+1][o]
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)
  
  #creates plot
  plot( x@reductions[[red]]@cell.embeddings[o,dims] , xlab="",ylab="",type="n", xaxt="n",yaxt="n", main=ifelse(main!=NULL, main, paste0(fn,"_",assay)),...)
  
  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(a@graphs)  ){
    add_graph(a,red=x@reductions[[red]]@cell.embeddings[o,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}
  
  #compute density and centroids
  if(!is.null(label) ){
    dens <- compute_density(x@reductions[[red]]@cell.embeddings[,dims],nbin=nbin)
    centroids <-  sapply( unique(label) , red=x@reductions[[red]]@cell.embeddings[,dims], cl1=as.character(label), function(jj,red,cl1) { apply( red[cl1==jj,],2,function(ii) {ndens <- dens[cl1==jj];weighted.median(x = ii,w = 1/(1.1-ndens/(max(ndens))) )} )  })
    colnames(centroids) <- unique(label)}
  
  #adds points
  points(x@reductions[[red]]@cell.embeddings[o,dims],pch=pch,cex=cex,bg=bg, col=paste0(pal) )
  options(warn=0)
  
  #adds labels
  if(!is.null(label)){
    label <- factor(x@meta.data[[label]])
    points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = unique(label))
  }
}




plot_meta <- function(x,red="umap",feat=NULL,pch=16,cex=.3,label=F,dims=c(1,2),
                      add_graph=NULL,percent_connections=1,nbin=400,add_lines=F,main=NULL,...){
  fn <- feat
  feat <- factor(x@meta.data[[feat]])
  if( !is.na(sum(as.numeric(levels(feat)))) ){
    levels(feat) <- levels(feat) [ order(as.numeric(levels(feat))) ] }
  
  pal <- c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") )
  pal <- pal[feat]
  #par(mar=c(1.5,1.5,1.5,1.5))
  options(warn=-1)
  
  #creates plot
  plot( x@reductions[[red]]@cell.embeddings[,dims], type="n", xaxt="n",yaxt="n", main=ifelse(main!=NULL, main, fn),...)

  #adds underlying graph
  if(!is.null(add_graph) ){ if( add_graph %in% names(a@graphs)  ){
    add_graph(a,red=x@reductions[[red]]@cell.embeddings[,dims],graph=add_graph,percent_connections=percent_connections) }else{message("Graph not found!")}}
  
  #compute density and centroids
  if(label | add_lines){
    dens <- compute_density(x@reductions[[red]]@cell.embeddings[,dims],nbin=nbin)
    centroids <-  sapply( as.character(unique(feat)) , reds=as.data.frame(x@reductions[[red]]@cell.embeddings[,dims]), cl1=feat, function(jj,reds,cl1) { apply( reds[cl1==jj,],2,function(ii) {ndens <- dens[cl1==jj];weighted.median(x = ii,w = 1/(1.1-ndens/(max(ndens))) )} )  })
    colnames(centroids) <- unique(feat)}
  
  if(add_lines){
    add_centroid_lines(x@reductions[[red]]@cell.embeddings[,dims],feat,pal,centroids)}
  
  #adds points
  points(x@reductions[[red]]@cell.embeddings[,dims],pch=pch,cex=cex,bg=paste0(pal,90), col=pal )
  options(warn=0)
  
  #adds labels
  if(label){
    points(centroids[1,],centroids[2,],col="#ffffff90",pch=16,cex=2)
    text(centroids[1,],centroids[2,],labels = unique(feat))
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
  apply(eee,1,function(x){ lines(x[3:4],x[5:6],col="#b3b3b3",lwd=.05) })
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






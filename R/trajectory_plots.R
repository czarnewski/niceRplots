plot_smoothers <- function(seuratObject,
                           SlingshotCurves,
                           gene,
                           clustering_use,
                           plotOrder=T,
                           col=pal,
                           colorByPseudotime=T,
                           spar=.9,
                           pt.size=.3,
                           factorScale=1.3){
  
  pt <- slingshot::slingPseudotime(curves)
  pt <- apply(pt,1,function(x) min(x,na.rm = T))
  ooo <- order(pt)[order(pt)]
  
  if(plotOrder){
    x <- 1:length(pt)
    y <- seuratObject@assays$RNA@data[gene,order(pt)]
    col <- col[seuratObject@meta.data[order(pt),clustering_use]]
  } else {
    x <- pt
    y <- seuratObject@assays$RNA@data[gene,]
    col <- col[seuratObject@meta.data[,clustering_use]]
  }
  
  plot(x,y,col=col,cex=pt.size,pch=16,xlab="Pseudotime",ylab="",las=1,main=gene)
  sm <- smooth.spline( x,y , spar = spar)
  lines(sm$x,sm$y*factorScale,lwd=2)
}




















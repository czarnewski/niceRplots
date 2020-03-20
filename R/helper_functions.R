require(scales)
require(raster)
#MAIN Violin plot function
#---------------
violins <- function(data, gene, clustering, plot_points=T,plot_y_axis=T,plot_x_axis=T,smooth=2,method="log",points_method="uniform",col="default",
                    pt.col="grey",pt.cex=.5,pt.pch=16,bw=.7,max_points=200,assay="RNA",ylab="expression",cex.main=1,main=gene,cex.axis=1,...){
  if(gene %in% rownames(data@assays[[assay]]@data) ){
    feat <- data@assays[[assay]]@data[gene,]
  } else if(gene %in% colnames(data@meta.data) ) {
    feat <- data@meta.data[, gene]
  } else { message("Feature or metadata not found!!") }

  #par(mar=c(2,3,2,1))
  n <- length(unique(data@meta.data[,clustering]))
  my_max <- max(max(feat,na.rm = T),.00000001,na.rm = T)*1.1

  plot(c(.4,n+.6),c(-1,-1), ylim=c(-.1,my_max),...,ylab="",type="n" ,frame.plot = F,yaxs="i",xaxs="i",las=1,xlab="",main=main,xaxt = "n",cex.main=cex.main)
  lines(c(0,n+.6),c(my_max,my_max),col="white",lwd=6,xpd=T)
  lines(c(n+.6,n+.6),c(0,my_max),col="white",lwd=6,xpd=T)
  mtext(side = 2, text = ylab, line = 2,las=3)
  #col_pal <- hue_pal()(length(unique(data@meta.data[,clustering])))
  col_pal <- c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") )
  col_pal <- col_pal[sort(unique(data@meta.data[,clustering]))]

  if(col=="default"){col <- paste0(col_pal,95)} else { col <- rep(col, length(unique(data@meta.data[,clustering])) )}

  for(i in 1:length(unique(data@meta.data[,clustering]))){
    cl <- sort(unique(data@meta.data[,clustering]))[i]
    #if(plot_points){
      #points(runif(sum(DATA@meta.data[,clustering] == cl),min = i-.4,max = i+.4),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey40")
      #points(rnorm(sum(DATA@meta.data[,clustering] == cl),mean = i,sd = .12),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey60")
    #}
    x <- na.omit(feat[data@meta.data[,clustering] == cl])
    suppressWarnings(suppressMessages( try(draw_violin(x, at = i,base=0.0,col = col[i], smooth=smooth,plot_points=plot_points,method=method,points_method=points_method,
               bw = bw,border =  "grey20",max_points=max_points,pt.col=pt.col,pt.cex=pt.cex,pt.pch=pt.pch)) ))
    #paste0(col_pal[i])
    #vioplot(x,at = i,add=T,col = paste0(col_pal[i],95),
            #drawRect = F,wex = 1,h = .01, border =  paste0(col_pal[i]))
  }

  abline(h=-.1,v=.4,xpd=F,lwd=2)
  if(plot_x_axis){
    axis(1, at=1:n, labels=sort(unique(data@meta.data[,clustering])),cex.axis=cex.axis)
  }
  # if(plot_y_axis){
  #     axis(2, at=seq(-100,100,by = 1), labels=seq(-100,100,by = 1),cex.axis=cex.axis,las=1)
  # }
}
#---------------





#MAIN Violin plot function
#---------------
violist <- function(data, genes, clustering, plot_points=T,plot_y_axis=T,plot_x_axis=T,smooth=2,method="log",points_method="uniform",col="default",
                    pt.col="grey",pt.cex=.5,pt.pch=16,bw=.7,max_points=200,assay="RNA",ylab="expression",cex.main=1,main=gene,cex.axis=1,...){
  n <- length(unique(data@meta.data[,clustering]))
  plot(c(.4,n+.6),c(-1,-1), ylim=c(0,length(genes)),ylab="",type="n" ,frame.plot = F,yaxs="i",xaxs="i",las=1,xlab="",main="",xaxt = "n",yaxt = "n",cex.main=cex.main)

  #col_pal <- hue_pal()(length(unique(data@meta.data[,clustering])))
  col_pal <- c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") )
  col_pal <- col_pal[sort(unique(data@meta.data[,clustering]))]
  if(col=="default"){col <- paste0(col_pal,95)} else { col <- rep(col, length(unique(data@meta.data[,clustering])) )}

  panel_row <- length(genes)
  for(gene in genes){
    panel_row <- panel_row - 1

    if(gene %in% rownames(data@assays[[assay]]@data) ){
      feat <- data@assays[[assay]]@data[gene,]
    } else if(gene %in% colnames(data@meta.data) ) {
      feat <- data@meta.data[, gene]
    } else { message("Feature or metadata not found!!") }

    #par(mar=c(2,3,2,1))
    my_max <- max(max(feat,na.rm = T),.00000001,na.rm = T)*1.1

    for(i in 1:length(unique(data@meta.data[,clustering]))){
      cl <- sort(unique(data@meta.data[,clustering]))[i]
      #if(plot_points){
      #points(runif(sum(DATA@meta.data[,clustering] == cl),min = i-.4,max = i+.4),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey40")
      #points(rnorm(sum(DATA@meta.data[,clustering] == cl),mean = i,sd = .12),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey60")
      #}
      x <- na.omit(feat[data@meta.data[,clustering] == cl]) / my_max
      suppressWarnings(suppressMessages( try(draw_violin(x, base=panel_row+.02,at = i,col = col[i], smooth=smooth,plot_points=plot_points,method=method,points_method=points_method,
                                                         bw = bw,border =  "grey20",max_points=max_points,pt.col=pt.col,pt.cex=pt.cex,pt.pch=pt.pch)) ))
      #paste0(col_pal[i])
      #vioplot(x,at = i,add=T,col = paste0(col_pal[i],95),
      #drawRect = F,wex = 1,h = .01, border =  paste0(col_pal[i]))
    }

  # if(plot_y_axis){
  #     axis(2, at=seq(-100,100,by = 1), labels=seq(-100,100,by = 1),cex.axis=cex.axis,las=1)
  # }
    lines(c(0,n+.6),c(panel_row,panel_row),col="black",lwd=1,xpd=F)
  }


  lines(c(0,n+.6),c(length(genes),length(genes)),col="white",lwd=6,xpd=T)
  lines(c(n+.6,n+.6),c(0,length(genes)),col="white",lwd=6,xpd=T)
  abline(h=-.1,v=.4,xpd=F,lwd=2)
  # if(plot_x_axis){
  #   axis(1, at=1:n, labels=sort(unique(data@meta.data[,clustering])),cex.axis=cex.axis)
  # }

  mtext(at = (length(genes):1)-.5 , side = 2, text = genes, las=1 )

  text( 1:n, par("usr")[3], labels = sort(unique(data@meta.data[,clustering])), srt = 45, adj = c(1,1.5), xpd = T, cex=1)

}
#---------------




#Function to calculate violin density
#---------------
draw_violin <- function(x,base=0,method="log",plot_points=F,points_method="proportional",smooth=2,col="grey",
                        border="grey",at=1,pt.col="grey",pt.cex=0.5,pt.pch=16,bw=0.45,max_points=200){
  r <- sum(x!=0)/length(x)
  if(plot_points){
    if(points_method == "proportional"){points(rnorm(length(x),mean = at,r/5),x+base,cex=pt.cex,col=pt.col,pch=pt.pch)}
    if(points_method == "uniform"){points(rnorm(length(x),mean = at,sd = .12),x+base,cex=pt.cex,col=pt.col,pch=pt.pch)}
  }

  if(method == "uniform"){
    a <- density(x,bw = .25*smooth,n = 200)
    ys <- a$x
    xs <- a$y/max(a$y)*bw/2
  }

  if(method == "log"){
    x2 <- log2(log2(x+1)+1)
    a <- density(x2,bw = .05*smooth,n = 200)
    ys <- 2^(2^a$x-1)-1
    xs <- a$y/max(a$y)*bw/2
  }

  if(method == "mixed"){
    d1 <- density(x[x==0],bw = .01*smooth,n = 200)
    d1$y <- d1$y/max(d1$y)*(1-r)
    d2 <- density(c(x[x!=0]),bw = .2,n = 200)
    d2$y <- d2$y/max(d2$y)*(r)
    ys <- c(d1$x,d2$x)
    xs <- c(d1$y,d2$y)*bw
  }
  ulim <- max(quantile(x,1),0.015)
  llim <- max(min(x),-0.1)
  polygon( c(xs[ys<ulim & ys>llim] , -rev(xs[ys<ulim& ys>llim]) )+at,
           c(ys[ys<ulim & ys>llim]+base, rev(ys[ys<ulim& ys>llim])+base) ,col = col,border = border)
}
#---------------




#Function to plot dot gene averages
#---------------
plot_dots <- function(data, gene_list, clustering, pal=c("grey80","blue3"),main="",srt=90,cex.row=1,cex.col=1,show_grid=T,min_size=.5,show_axis=F){
  x1 <- rowsum(t(as.matrix(data@assays$RNA@data[unique(gene_list),])), data@meta.data[,clustering])
  x1 <- t(x1 / c(table(factor(data@meta.data[,clustering]))))
  x1 <- t(apply(t(x1) , 2,function(i) (i-0)/(max(i)-0) ) )

  x2 <- rowsum(( t (as.matrix(data@assays$RNA@data[unique(gene_list),]!=0)) *1), data@meta.data[,clustering])
  x2 <- t(x2 / c(table(factor(data@meta.data[,clustering]))))


  plot(0,0,type="n",las=1,xlim=c(1,ncol(x1)),ylim=c(1,nrow(x1)),axes=F,ylab="",xlab="",main=main)

  if(show_grid){
    for(i in 1:length(unique(data@meta.data[,clustering]) )){
      lines(c(i,i),c(length(unique(gene_list))+.5,0),col="grey95",lwd=.5,xpd=F)
    }
    for(i in 1:length(unique(gene_list)) ){
      lines(c(0,length(unique(data@meta.data[,clustering]) )+0.5),c(i,i),col="grey95",lwd=.5,xpd=F)
    }
  }

  points(rep(1:ncol(x1),nrow(x1)),sort(rep(1:(nrow(x1)),ncol(x1))), cex=c(t(x2) )*2+min_size,
         pch=16, col=c( "grey95",colorRampPalette(pal)(19))[c(t(x1) )*18+1 ],xpd=T)

  text(1:ncol(x1), par("usr")[3]-.3, labels = colnames(x1), srt = srt, adj = c(1,0.5), xpd = TRUE, cex=cex.row)
  text(0.4, 1:nrow(x1) , labels = rownames(x1), srt = 0, adj = c(1,0.5), xpd = TRUE, cex=cex.col)
  if(show_axis){
    lines(.5+c(ncol(x1),0,0),c(0,0,nrow(x1))+.5)
    lines(.5+c(0,length(unique(data@meta.data[,clustering]) )),c(par("usr")[3],par("usr")[3]),col="black",lwd=1,xpd=T)
  }
}
#---------------




#Heatmap plot
#---------------
plot_heat <- function(data, genes, order_metadata=NULL, annot=NULL, cut_max=2, row.cex=1, main="", heat_color=colorRampPalette(c("purple4","black","black","yellow4","yellow1") )(90),...){
  plot(0,type="n",ylim=c(0,1.105),xlim=c(0,1.25),col="white",axes=F,asp=F,main=main,xlab="",ylab="")

  if(!is.null(annot)){
    increment <- min( 0.1 / length(annot), .02)
    begin <- 1.01
    for(i in annot){
      end <- begin + increment
      if(!is.null(order_metadata)){
        annnn <- raster(matrix(as.numeric(data@meta.data[,i][order(data@meta.data[,order_metadata ])]),nrow = 1))
      } else {
        annnn <- raster(matrix(as.numeric(data@meta.data[,i]),nrow = 1))
      }
      annnn@extent <- extent(c(0,1),c(begin,end))
      plot(annnn,col=hue_pal()(length(unique(data@meta.data[,i]) )),axes=F,asp=F,ylim=c(0,4),legend=F,add=T,interpolate=F)
      text( 1.01 , (begin+end)/2 , labels = i, srt = 0, adj = c(0,0.5), xpd = TRUE, cex=row.cex,font=2)
      begin <- begin + increment
    }
    #text( 1.01 , seq(1,length(increment))*increment+1 , labels = annot, srt = 0, adj = c(0,0.5), xpd = TRUE, cex=row.cex,font=2)
  }
  # annnn <- raster(matrix(as.numeric(data$cell_tissue_donor_plate[order(data$cell_tissue_donor_plate)]),nrow = 1))
  # annnn@extent <- extent(c(0,1),c(1.01,1.025))
  # plot(annnn,col=hue_pal()(length(levels(data$cell_tissue_donor_plate) )),axes=F,asp=F,ylim=c(0,4),legend=F,add=T,interpolate=F)
  #
  # annnn <- raster(matrix(as.numeric(data$Donor[order(data$cell_tissue_donor_plate)]),nrow = 1))
  # annnn@extent <- extent(c(0,1),c(1.03,1.045))
  # plot(annnn,col=hue_pal()(length(levels(data$Donor) )),axes=F,asp=F,ylim=c(0,4),legend=F,add=T,interpolate=F)
  #
  # annnn <- raster(matrix(as.numeric(data$cell_tissue[order(data$cell_tissue_donor_plate)]),nrow = 1))
  # annnn@extent <- extent(c(0,1),c(1.05,1.065))
  # plot(annnn,col=hue_pal()(length(levels(data$cell_tissue) )),axes=F,asp=F,ylim=c(0,4),legend=F,add=T,interpolate=F)
  #
  # annnn <- raster(matrix(as.numeric(data$Celltype[order(data$cell_tissue_donor_plate)]),nrow = 1))
  # annnn@extent <- extent(c(0,1),c(1.07,1.085))
  # plot(annnn,col=hue_pal()(length(levels(data$Celltype) )),axes=F,asp=F,ylim=c(0,4),legend=F,add=T,interpolate=F)
  #
  # annnn <- raster(matrix(as.numeric(data$Tissue[order(data$cell_tissue_donor_plate)]),nrow = 1))
  # annnn@extent <- extent(c(0,1),c(1.09,1.105))
  # plot(annnn,col=hue_pal()(length(levels(data$Tissue) )),axes=F,asp=F,ylim=c(0,4),legend=F,add=T,interpolate=F)
  #
  if(!is.null(order_metadata)){
    teeeeest <- t(apply(data@assays$RNA@data[genes,] , 1, function(x) scale(x,T,T)))[,order(data@meta.data[,order_metadata])]
  } else {
    teeeeest <- t(apply(data@assays$RNA@data[genes,] , 1, function(x) scale(x,T,T)))[,order(data$cell_tissue_donor_plate)]
  }
  teeeeest[teeeeest > 2]  <- cut_max
  teeeeest[teeeeest < -2] <- cut_max
  #teeeeest <- t(apply(teeeeest , 1, function(x) (x - min(x)) / (max(x)-min(x)) ))
  #teeeeest[,1000] <- NA
  plot(raster(teeeeest),col= heat_color,asp=F,axes=F,add=T,interpolate=F,...)
  text( 1.01 , (1:nrow(teeeeest)-.5)/nrow(teeeeest) , labels = rev(rownames(teeeeest)), srt = 0, adj = c(0,0.5), xpd = TRUE, cex=row.cex,add=T)

}
#---------------






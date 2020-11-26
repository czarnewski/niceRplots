require(scales)


#MAIN Violin plot function
#---------------
violins <- function(data, gene=NULL, clustering=NULL, plot_points=T,plot_y_axis=T,transparency=NULL,plot_x_axis=T,smooth=.1,method="log",points_method="proportional",col=c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") ),
                    pt.col="grey",pt.cex=.5,pt.pch=16,bw=.7,max_points=200,assay="RNA",srt=0,ylab="expression",cex.main=1,main=gene,cex.axis=1,...){


  if( is(object = data,"Seurat") ){
    if(gene %in% rownames(data@assays[[assay]]@data) ){
      feat <- data@assays[[assay]]@data[gene,]
    } else if(gene %in% colnames(data@meta.data) ) {
      feat <- data@meta.data[, gene]
    } else { message("Feature or metadata not found!!") }
    temp <- factor(as.character(data@meta.data[,clustering]))
    N <- ncol(data)
  } else {
    feat <- c(as.matrix(data))
    temp <- factor(c(sapply(colnames(data),function(x){rep(x,nrow(data))})),levels = colnames(data))
    N <- nrow(data)
  }

  pt.col=rep( pt.col,N)[1:N]
  pt.cex=rep(pt.cex,N)[1:N]
  pt.pch=rep(pt.pch,N)[1:N]

  if( !is.na(sum(as.numeric(levels(temp)))) ){
    temp <- factor(as.numeric(as.character(temp))) }

  #par(mar=c(2,3,2,1))
  n <- length(levels(temp))
  print(levels(temp))
  my_max <- max(max(feat,na.rm = T),.00000001,na.rm = T)*1.1

  plot(c(.4,n+.6),c(-1,-1), ylim=c(-.1,my_max),...,ylab="",type="n" ,frame.plot = F,
       yaxs="i",xaxs="i",las=1,xlab="",main=main,xaxt = "n",cex.main=cex.main)
  lines(c(0,n+.6),c(my_max,my_max),col="white",lwd=6,xpd=T)
  lines(c(n+.6,n+.6),c(0,my_max),col="white",lwd=6,xpd=T)
  mtext(side = 2, text = ylab, line = 2,las=3)
  #col_pal <- hue_pal()(length(unique(data@meta.data[,clustering])))
  col <- col[1:n]

  if(col[1]=="default"){col <- paste0(col,95)} else { col <- rep(col, n )}

  for(i in 1:n){
    cl <- levels(temp)[i]
    #if(plot_points){
      #points(runif(sum(DATA@meta.data[,clustering] == cl),min = i-.4,max = i+.4),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey40")
      #points(rnorm(sum(DATA@meta.data[,clustering] == cl),mean = i,sd = .12),DATA@data[gene,DATA@meta.data[,clustering] == cl],cex=.3,pch=16,col="grey60")
    #}
    x <- na.omit(feat[temp == cl])
    suppressWarnings(suppressMessages( try(draw_violin(x, at = i,base=0.0,col = col[i],
                                                       smooth=smooth,
                                                       plot_points=plot_points,
                                                       method=method,
                                                       points_method=points_method,
                                                       transparency=transparency,
               bw = bw,border =  "grey20",max_points=max_points,pt.col=pt.col,pt.cex=pt.cex,pt.pch=pt.pch)) ))
    #paste0(col_pal[i])
    #vioplot(x,at = i,add=T,col = paste0(col_pal[i],95),
            #drawRect = F,wex = 1,h = .01, border =  paste0(col_pal[i]))
  }

  abline(h=-.1,v=.4,xpd=F,lwd=2)
  if(plot_x_axis){
    # axis(1, at=1:n, labels=sort(unique(data@meta.data[,clustering])),cex.axis=cex.axis)
    text( 1:n, par("usr")[3] - (par("usr")[4])/200, labels = levels(temp), srt = srt, adj = c(ifelse(srt==0,.5,ifelse(srt==90,1,1)),ifelse(srt==0,1,ifelse(srt==90,.5,1))),, xpd = T, cex=cex.axis)
  }
  # if(plot_y_axis){
  #     axis(2, at=seq(-100,100,by = 1), labels=seq(-100,100,by = 1),cex.axis=cex.axis,las=1)
  # }
}
#---------------





#MAIN Violin plot function
#---------------
violist <- function(data, genes, clustering, plot_points=T,plot_y_axis=T,plot_x_axis=T,smooth=.5,method="log",points_method="proportional",srt=0,transparency=NULL,
                    pt.col="grey",
                    pt.cex=.5,
                    pt.pch=16,
                    bw=.7,max_points=200,assay="RNA",ylab="expression",cex.main=1,main=gene,cex.axis=1,col = c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") ),...){

  if(is(data,"Seurat")){
    data <- data@assays[[assay]]@data[genes,]
    grouping <- data@meta.data[,clustering]
  } else {grouping <- factor(clustering)}


  n <- length(unique(grouping))
  N <- ncol(data)
  pt.col=rep( pt.col,N)[1:N]
  pt.cex=rep(pt.cex,N)[1:N]
  pt.pch=rep(pt.pch,N)[1:N]



  plot(c(.4,n+.6),c(-1,-1), ylim=c(0,length(genes)),ylab="",type="n" ,frame.plot = F,yaxs="i",xaxs="i",las=1,xlab="",main="",xaxt = "n",yaxt = "n",cex.main=cex.main)

  col <- colorRampPalette(col)(length(col))
  col <- col[as.factor(sort(unique(grouping)))]

  temp <- factor(as.character(grouping))
  if( !is.na(sum(as.numeric(levels(temp)))) ){
    temp <- factor(as.numeric(as.character(grouping))) }

  panel_row <- length(genes)
  for(gene in genes){
    panel_row <- panel_row - 1

    if(gene %in% rownames(data) ){
      feat <- data[gene,]
    } else if(gene %in% colnames(data@meta.data) ) {
      feat <- data@meta.data[, gene]
    } else { message("Feature or metadata not found!!") }

    #par(mar=c(2,3,2,1))
    my_max <- max(max(feat,na.rm = T),.00000001,na.rm = T)*1.1
    # message(paste0("feat:",gene,"\tmax:",my_max,"\tmin:",panel_row+.02,"\tdatamax:",max(feat),"\tdatamin:",min(feat)))

    for(i in 1:length(levels(temp))){
      cl <- levels(temp)[i]
      #if(plot_points){
      #points(runif(sum(grouping == cl),min = i-.4,max = i+.4),DATA@data[gene,grouping == cl],cex=.3,pch=16,col="grey40")
      #points(rnorm(sum(grouping == cl),mean = i,sd = .12),DATA@data[gene,grouping == cl],cex=.3,pch=16,col="grey60")
      #}
      x <- na.omit(feat[temp == cl]) / my_max

      suppressWarnings(suppressMessages( try(draw_violin(x, base=panel_row+.02,at = i,col = col[i],
                                                         smooth=smooth,
                                                         plot_points=plot_points,
                                                         method=method,
                                                         points_method=points_method,
                                                         bw = bw,
                                                         border =  "grey20",
                                                         max_points=max_points,
                                                         pt.col=pt.col[temp == cl],
                                                         pt.cex=pt.cex[temp == cl],
                                                         pt.pch=pt.pch[temp == cl],
                                                         transparency=transparency)) ))
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
  #   axis(1, at=1:n, labels=sort(unique(grouping)),cex.axis=cex.axis)
  # }

  # mtext(at = (length(genes):1)-.5 , side = 2, text = genes, las=1 )
  text( par("usr")[1] - (par("usr")[2])/200 , (length(genes):1)-.5, labels = genes,
        srt = 0, adj = c(1,.5), xpd = T, cex=cex.axis)
  text( 1:n, par("usr")[3] - (par("usr")[4])/200, labels = levels(temp), srt = srt, adj = c(ifelse(srt==0,.5,ifelse(srt==90,1,1)),ifelse(srt==0,1,ifelse(srt==90,.5,1))),, xpd = T, cex=cex.axis)

}
#---------------




#Function to calculate violin density
#---------------
draw_violin <- function(x,base=0,method="log",plot_points=F,points_method="proportional",smooth=2,col="grey",
                        border="grey",at=1,pt.col="grey",pt.cex=0.5,pt.pch=16,bw=0.45,max_points=200,transparency=NULL){
  r <- sum(x!=0)/length(x)
  if(plot_points){
    if(points_method == "proportional"){
      set.seed(1)
      points(rnorm(length(x),mean = at,sd = .12),x+base,cex=pt.cex,col=pt.col,pch=pt.pch)}
    if(points_method == "uniform"){
      set.seed(1)
      points(runif(length(x),min = at-.4,max = at+.4),x+base,cex=pt.cex,col=pt.col,pch=pt.pch)}
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
  ulim <- max(max(x),0.015)
  llim <- max(min(x),-0.1)
  polygon( c(xs[ys<ulim & ys>llim] , -rev(xs[ys<ulim& ys>llim]) )+at,
           c(ys[ys<ulim & ys>llim]+base, rev(ys[ys<ulim& ys>llim])+base) ,
           col = ifelse(is.null(transparency) , col, paste0(col,transparency) ),
           border = border)
}
#---------------




#Function to plot dot gene averages
#---------------
plot_dots <- function(data, genes, clustering, pal=c("grey90","grey70","navy"),main="",
                      srt=0,cex.row=1,cex.col=1,show_grid=T,min_size=.5,show_axis=T,...){

  if(is(data,"Seurat")){
    data <- data@assays[[assay]]@data[genes,]
    grouping <- data@meta.data[,clustering]
  } else {grouping <- factor(clustering)}


  temp <- factor(as.character(grouping))
  if( !is.na(sum(as.numeric(levels(temp)))) ){
    temp <- factor(as.numeric(as.character(grouping))) }

  x1 <- rowsum(t(as.matrix(data[rev(genes),])), temp)
  x1 <- t(x1 / c(table(temp)))
  x1 <- t(apply(t(x1) , 2,function(i) (i-0)/(max(i)-0) ) )

  x2 <- rowsum(( t (as.matrix(data[rev(genes),]!=0)) *1), temp)
  x2 <- t(x2 / c(table(temp)))


  plot(0,0,type="n",las=1,xlim=c(.5,ncol(x1)+.5),ylim=c(.5,nrow(x1)+.5),
       axes=F,ylab="",xlab="",main=main,xaxs="i",yaxs="i",...)

  if(show_grid){
    for(i in 1:length(levels(temp) )){
      lines(c(i,i),c(length(unique(genes))+.5,0),col="grey95",lwd=.5,xpd=F)
    }
    for(i in 1:length(unique(genes)) ){
      lines(c(0,length(levels(temp) )+0.5),c(i,i),col="grey95",lwd=.5,xpd=F)
    }
  }

  points(rep(1:ncol(x1),nrow(x1)),sort(rep(1:(nrow(x1)),ncol(x1))), cex=c(t(x2) )*2+min_size,
         pch=16, col=c( "grey95",colorRampPalette(pal)(19))[c(t(x1) )*18+1 ],xpd=T)

  text(1:ncol(x1), par("usr")[3] - (par("usr")[4])/200, labels = colnames(x1), srt = srt,
       adj = c(ifelse(srt==0,.5,ifelse(srt==90,1,1)),ifelse(srt==0,1,ifelse(srt==90,.5,1))),
       xpd = TRUE, cex=cex.col)
  text(par("usr")[1] - (par("usr")[2])/200, 1:nrow(x1) , labels = rownames(x1), srt = 0, adj = c(1,0.5), xpd = TRUE, cex=cex.row)

  if(show_axis){
    lines(c(.5,.5),c(.5,nrow(x1)+.5),col="black",lwd=1,xpd=T)
    lines(.5+c(0, length(levels(temp) )),c(.5,.5),col="black",lwd=1,xpd=T)
    par(xpd=F)
  }
}
#---------------




#Heatmap plot
#---------------
plot_heat <- function(data, genes, order_metadata=NULL, annot=NULL, cut_max=2, row.cex=1, main="", heat_color=c("grey90",colorRampPalette(c("grey80","navy","navy") )(90)),...){
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
    teeeeest <- t(apply(data@assays$RNA@data[genes,] , 1, function(x) (x - min(x)) / (max(x)-min(x)) ))[,order(data@meta.data[,order_metadata])]
  } else {
    teeeeest <- t(apply(data@assays$RNA@data[genes,] , 1, function(x) (x - min(x)) / (max(x)-min(x)) ))
  }
  teeeeest[teeeeest > 2]  <- cut_max
  teeeeest[teeeeest < -2] <- cut_max
  #teeeeest <- t(apply(teeeeest , 1, function(x) (x - min(x)) / (max(x)-min(x)) ))
  #teeeeest[,1000] <- NA
  plot(raster(teeeeest),col= heat_color,asp=F,axes=F,add=T,interpolate=F,...)
  text( 1.01 , (1:nrow(teeeeest)-.5)/nrow(teeeeest) , labels = rev(rownames(teeeeest)), srt = 0, adj = c(0,0.5), xpd = TRUE, cex=row.cex)

}
#---------------



###Gene expression barplots
#---------------
plot_bars <- function(data, gene, clustering,assay="RNA",
                      col=c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") ),
                      ...){
  if(gene %in% rownames(data@assays[[assay]]@data) ){
    feat <- data@assays[[assay]]@data[gene,]
  } else if(gene %in% colnames(data@meta.data) ) {
    feat <- data@meta.data[, gene]
  } else { message("Feature or metadata not found!!") }

  temp <- factor(as.character(data@meta.data[,clustering]))
  if( !is.na(sum(as.numeric(levels(temp)))) ){
    temp <- factor(as.numeric(as.character(data@meta.data[,clustering]))) }

  o <- order(temp)
  barplot( feat[o] , col=col [ temp[o] ],
           border= col[ temp[o]],
           names.arg="" , las=1, xaxs="i", yaxs="i",
           ylim=c( min(c(feat,0))*1.2, max(c(feat,0))*1.2 ),...)
  lines( c(0,length(feat)*1.2-.5), c(min(c(feat,0))*1.2,min(c(feat,0))*1.2),xpd=T ) ; par(xpd=F)
  mmm <- sapply(as.character(levels(temp)),
                function(i) mean( feat [temp == i] ))
  lines((1:length(feat))*1.2-.5,mmm[match(temp[o], as.character(unique(temp[o])))],lwd=2)
}
#---------------



###Gene expression barplots
#---------------
barlist <- function(data, genes, clustering, plot_y_axis=T,plot_x_axis=T,labels=NULL,srt=0,orderby = NULL,
                    assay="RNA",ylab="expression",font.main=1,cex.main=1,main="",cex.axis=1,
                    col = "default",draw_mean_lines=T,...){

  if(is(data,"Seurat")){
    data <- data@assays[[assay]]@data[genes,]
    grouping <- data@meta.data[,clustering]
  } else {grouping <- factor(clustering)}

  n <- length(unique(grouping))
  plot(c(0,ncol(data)*1.2),c(-1,-1),
       ylim=c(0,length(genes)),ylab="",type="n" ,frame.plot = F,yaxs="i",
       xaxs="i",las=1,xlab="",main=main,xaxt = "n",yaxt = "n",
       cex.main=cex.main,font.main=font.main)



  temp <- factor(as.character(grouping))
  if( !is.na(sum(as.numeric(levels(temp)))) ){
    temp <- factor(as.numeric(as.character(grouping))) }
  o <- order(temp)


  if(!is.null(orderby)){
    if(orderby %in% rownames(data) ){
      feat <- data[orderby,]
    } else if( is(data,"Seurat") ) {
      if(orderby %in% colnames(data@meta.data) ) {
        feat <- data@meta.data[, orderby]} else { message("Feature or metadata not found!!") }
    } else { message("Feature or metadata not found!!") }

    for(i in levels(temp)){
      o[temp[o] == i] <- o[temp[o] == i][order(feat[o][temp[o] == i],decreasing = F)]
    }
  }


  if(col[1]=="default"){
    col <- c(scales::hue_pal()(8),RColorBrewer::brewer.pal(9,"Set1"),RColorBrewer::brewer.pal(8,"Set2") )
    col <- col [ temp[o] ]
  } else if( length(col) == ncol(data) ){
    col <- rep(col, ncol(data))[1:ncol(data)]
    col <- col[o]
  } else {
    col <- col [ temp[o] ]
  }

  panel_row <- length(genes)
  for(gene in genes){
    panel_row <- panel_row - 1

    if(gene %in% rownames(data) ){
      feat <- data[gene,]
    } else if( is(data,"Seurat") ) {
      if(gene %in% colnames(data@meta.data) ) {
      feat <- data@meta.data[, gene]} else { message("Feature or metadata not found!!") }
    } else { message("Feature or metadata not found!!") }

    #par(mar=c(2,3,2,1))
    my_max <- max(c(feat,.00000001),na.rm = T)*1.1
    x <- (na.omit(feat) / my_max)
    barplot( x[o] , col=col, axes=F,
             border= col,
             names.arg="" , las=1, xaxs="i", yaxs="i",offset = panel_row,
             ylim=c( 0 , length(genes) ),
             add=T,xpd=F,...)
    if(draw_mean_lines){
      mmm <- sapply(as.character(unique(temp[o])),
                    function(i) mean( x [temp == i]))
      lines((1:length(x))*1.2-.5,
            mmm[match(temp[o],
                      as.character(unique(temp[o])))]+panel_row,lwd=2)
    }

    lines(c(0,0),c(panel_row,panel_row+.9),xpd=F,lwd=2)
    lines(c(0, length(x)*1.2 ),c(panel_row,panel_row),xpd=F,lwd=2)

    # if(plot_y_axis){
    #     axis(2, at=seq(-100,100,by = 1), labels=seq(-100,100,by = 1),cex.axis=cex.axis,las=1)
    # }
    # lines( c(0,length(x)*1.2-.5), c(panel_row,panel_row),xpd=T ) ; par(xpd=F)
  }


  # lines(c(0,n+.6),c(length(genes),length(genes)),col="white",lwd=6,xpd=T)
  # lines(c(n+.6,n+.6),c(0,length(genes)),col="white",lwd=6,xpd=T)
  # abline(h=0,v=0,xpd=F,lwd=2)
  # if(plot_x_axis){
  #   axis(1, at=1:n, labels=sort(unique(grouping)),cex.axis=cex.axis)
  # }

  # mtext(at = (length(genes):1)-.5 , side = 2, text = genes, las=1 , cex = cex.axis/2)
  # mtext(at = (cumsum(table(data@meta.data[o,clustering])) - (table(data@meta.data[o,clustering])/2))*1.2 - .5 , side = 2, text = genes, las=1 )
  if(!is.null(labels)){ genes <- labels  }
  text( par("usr")[1] - (par("usr")[2])/200 , (length(genes):1)-.5, labels = genes,
        srt = 0, adj = c(1,.5), xpd = T, cex=cex.axis)

  text( (cumsum(table(temp))*1.2 - .5)  - ((table(temp)/2)*1.2 - .5 ), par("usr")[3],
        labels = sort(unique(temp)), srt = srt, adj = c(ifelse(srt==0,.5,ifelse(srt==90,1,1)),ifelse(srt==0,1,ifelse(srt==90,.5,1))),, xpd = T, cex=cex.axis)

}
#---------------




getcluster <- function(data, gene, clustering, lowest=F){

  if(is(data,"Seurat")){
    data <- data@assays[[assay]]@data[genes,]
    grouping <- data@meta.data[,clustering]
  } else {
    grouping <- clustering
  }

  data <- sapply(unique(as.character(grouping)),function(x){
    mean( data[ gene , as.character(grouping) == x ])
  })
  if(lowest){
    return( unique(as.character(grouping))[which.min(data)] )
  } else {
    return( unique(as.character(grouping))[which.max(data)] )
  }
}




pointcluster <- function(a, gene, clustering, red="umap", cex=3,col="black",lowest=F){
  clust_use <- getcluster(a,gene,clustering,lowest)
  points(t(pmean(a@reductions[[red]]@cell.embeddings[a@meta.data[,clustering] == clust_use,])) ,
         pch=1, cex=cex, xpd=T,col=col)
}




add_annotation_percentages <- function(data,annotation_file,
                                       annotation_column_use="external_gene_name",
                                       annotation_column_target="gene_biotype",assay="RNA"){
  item <- annot[match(rownames(data@assays[[assay]]@counts), annot[, annotation_column_use]), annotation_column_target]
  item[is.na(item)] <- "unknown"

  data@meta.data <- data@meta.data[, !(colnames(data@meta.data) %in% annotation_column_target)]
  data@assays[[assay]]@meta.features <-

  # Calculate the percentage of each gene biotype
  perc <- rowsum(as.matrix(data@assays[[assay]]@counts), group = item)
  perc <- ( t(perc) / Matrix::colSums(data@assays[[assay]]@counts) )
  o <- order(apply(perc, 2, mean), decreasing = FALSE)
  perc <- perc[, o]
  print(dim(perc))

  # Add table to the object
  annot_table <-
    setNames(as.data.frame((perc*100)[, names(sort(table(item), decreasing = TRUE))]),
             paste0("percent_", names(sort(table(item), decreasing = TRUE))))
  data@meta.data <- data@meta.data[, !(colnames(data@meta.data) %in% colnames(annot_table))]

  data@meta.data <- cbind(
    data@meta.data,
    annot_table)

  return(data)
}


plot_sankey <- function(df, plot_labels=T, plot_weights=T, color_by=1,
                        order_1_by="NULL", order_2_by="NULL",
                        xlim=c(0,1), ylim1=c(0,1), ylim2=c(0,1),
                        use_w1 = T, use_w2 = T, gapv = 0.03, gap2v = 0.005,
                        smoothing = 0.15, nbreaks = 100, add=F, pal=NULL,...){

  if( order_1_by == "total"){
    df[,1] <- factor(df[,1],levels = names(sort(table(df[,1]),decreasing = T))) }
  if( order_2_by == "total"){
    df[,2] <- factor(df[,2],levels = names(sort(table(df[,2]),decreasing = T))) }

  if(is.null(pal)){
    pal <- scales::hue_pal()(length(levels(df[,1])))
  } else {
    pal <- colorRampPalette(pal)(length(pal))
  }
  weights <- c(t(table(df)))
  nzw <- weights!=0

  eee <- data.frame(x=c(sapply(levels(df[,1]),function(x)rep(x,ncol(table(df))))),
                    y=rep(levels(df[,2]),nrow(table(df))),
                    w=c(t(table(df))) )
  eee <- eee[eee$w>0,]

  eee$u1 <- rev(cumsum(rev(eee$w)))

  eee$x_orig <- factor(eee$x,levels=levels(df[,1]))
  eee$y_orig <- factor(eee$y,levels=levels(df[,2]))

  eee$x <- sprintf("%03d", as.numeric(eee$x_orig))
  eee$y <- sprintf("%03d", as.numeric(eee$y_orig))

  eee$order1 <- order(paste0(eee$x,"_",eee$y))
  eee$order2 <- order(paste0(eee$y,"_",eee$x))

  rownames(eee) <- as.character(1:nrow(eee))


  # use_w1 <- T
  maxv <- ylim1[2]
  minv <- ylim1[1]
  # gapv <- .03
  # gap2v <- .005
  # smoothing <- 0.15
  # nbreaks <- 100
  ngaps <- length(unique(eee$x))-1
  sp <- (maxv-gapv*ngaps-gap2v*nrow(eee)-minv)/(ifelse(use_w1,sum(eee$w),nrow(eee)))

  tempu <- c()
  tempd <- c()
  eee <- eee[as.character(eee$order1),]
  for(i in 1:nrow(eee)){
    tempu <- c(tempu,maxv)
    maxv <- maxv - sp*ifelse(use_w1,eee$w[i],1)
    tempd <- c(tempd,maxv)
    if(i != nrow(eee)){
      if( eee$x[i] != eee$x[i+1] ){
        # message("gap!")
        maxv <- maxv - gapv
      }
    }
    maxv <- maxv - gap2v
  }
  eee$y1 <- tempu
  eee$y2 <- tempd


  # use_w2 <- T
  maxv <- ylim2[2]
  minv <- ylim2[1]
  ngaps <- length(unique(eee$y))-1
  sp <- (maxv-gapv*ngaps-gap2v*nrow(eee)-minv)/(ifelse(use_w2,sum(eee$w),nrow(eee)))
  tempu <- c()
  tempd <- c()
  eee <- eee[as.character(eee$order2),]
  for(i in 1:nrow(eee) ){
    tempu <- c(tempu,maxv)
    maxv <- maxv - sp*ifelse(use_w2,eee$w[i],1)
    tempd <- c(tempd,maxv)
    if(i != nrow(eee)){
      if( eee$y[i] != eee$y[i+1] ){
        # message("gap!")
        maxv <- maxv - gapv
      }
    }
    maxv <- maxv - gap2v
  }
  eee$y3 <- tempu
  eee$y4 <- tempd


  eee <- eee[order(eee$order1),]


  if(order_2_by == "disentangled"){

    temp <- apply( eee, 1, function(x){
      return(rep( (as.numeric(x["y1"])+as.numeric(x["y2"]) )/2 , as.numeric(x["w"]) )) })
    temp2 <- sapply( as.character(unique(eee$y_orig)), function(x){
      return(mean( as.numeric(unlist(temp[eee$y_orig==x])) )) })
    sort(temp2,T)

    # eee$y_orig <- factor(eee$y_orig,levels=names(sort(temp2,T)))
    eee$y <- sprintf("%03d", as.numeric(factor(eee$y_orig,levels=names(sort(temp2,T)))))

    eee <- eee[as.character(1:nrow(eee)),]

    eee$order1 <- order(paste0(eee$x,"_",eee$y))
    eee$order2 <- order(paste0(eee$y,"_",eee$x))


    # use_w1 <- T
    maxv <- ylim1[2]
    minv <- ylim1[1]
    # gapv <- .03
    # gap2v <- .005
    # smoothing <- 0.15
    # nbreaks <- 100
    ngaps <- length(unique(eee$x))-1
    sp <- (maxv-gapv*ngaps-gap2v*nrow(eee)-minv)/(ifelse(use_w1,sum(eee$w),nrow(eee)))

    tempu <- c()
    tempd <- c()
    # eee <- eee[as.character(eee$order1),]
    for(i in 1:nrow(eee)){
      tempu <- c(tempu,maxv)
      maxv <- maxv - sp*ifelse(use_w1,eee$w[i],1)
      tempd <- c(tempd,maxv)
      if(i != nrow(eee)){
        if( eee$x[i] != eee$x[i+1] ){
          # message("gap!")
          maxv <- maxv - gapv
        }
      }
      maxv <- maxv - gap2v
    }
    eee$y1 <- tempu
    eee$y2 <- tempd

    # use_w2 <- T
    maxv <- ylim2[2]
    minv <- ylim2[1]
    ngaps <- length(unique(eee$y))-1
    sp <- (maxv-gapv*ngaps-gap2v*nrow(eee)-minv)/(ifelse(use_w2,sum(eee$w),nrow(eee)))
    tempu <- c()
    tempd <- c()
    eee <- eee[as.character(eee$order2),]
    for(i in 1:nrow(eee) ){
      tempu <- c(tempu,maxv)
      maxv <- maxv - sp*ifelse(use_w2,eee$w[i],1)
      tempd <- c(tempd,maxv)
      if(i != nrow(eee)){
        if( eee$y[i] != eee$y[i+1] ){
          # message("gap!")
          maxv <- maxv - gapv
        }
      }
      maxv <- maxv - gap2v
    }
    eee$y3 <- tempu
    eee$y4 <- tempd


    eee <- eee[order(eee$order1),]

  }

  x <- seq(xlim[1],xlim[2],length.out = nbreaks)

  if( !add ){
    plot( 0 , 0,type="n", ylim=c(min(ylim1,ylim2), max(ylim1,ylim2)),xlim=xlim,frame=F,axes=F,ylab="",xlab="",...)
  }

  for(i in nrow(eee):1){
    ys <- SSlogis(x ,1, mean(xlim), smoothing*diff(range(xlim)))
    ys <- (ys - min(ys)) / (max(ys) - min(ys))*(eee$y3[i]-eee$y1[i])
    ys <- ys + ifelse( eee$y3[i]-eee$y1[i] < 0 , max(eee$y1[i],eee$y3[i]), min(eee$y1[i],eee$y3[i]))

    ys2 <- SSlogis(x ,1, mean(xlim), smoothing*diff(range(xlim)))
    ys2 <- (ys2 - min(ys2)) / (max(ys2) - min(ys2))*(eee$y4[i]-eee$y2[i])
    ys2 <- ys2 + ifelse( eee$y4[i]-eee$y2[i] < 0 , max(eee$y2[i],eee$y4[i]), min(eee$y2[i],eee$y4[i]))

    polygon(x = c(x, rev(x)),col=paste0(pal[factor(eee[,c("x","y")[color_by] ])[i]],"90"),border = F,
            y = c(ys , rev(ys2)  ))
  }

  if(plot_labels){
    for(i in levels(eee$x_orig)){
      text( x=xlim[1],y = mean(range(eee[ eee$x_orig == i, c("y1","y2") ])),labels = i,pos = 2,xpd=T )
      lines(c(xlim[1],xlim[1]),range(eee[ eee$x_orig == i, c("y1","y2") ]))}
    for(i in levels(eee$y_orig)){
      text( x=xlim[2],y = mean(range(eee[ eee$y_orig == i, c("y3","y4") ])),labels = i,pos = 4,xpd=T )
      lines(c(xlim[2],xlim[2]),range(eee[ eee$y_orig == i, c("y3","y4") ]))}
  }

  if(plot_weights){
    for(i in 1:nrow(eee)){
      text( x=xlim[1]+diff(xlim)*.01, y = mean(range(eee[ i, c("y1","y2") ])),
            labels = eee$w[i],adj = 0,xpd=T,cex= .3+.7*(eee$w[i]/max(eee$w)) )}
    for(i in 1:nrow(eee)){
      text( x=xlim[2]-diff(xlim)*.01, y = mean(range(eee[ i, c("y3","y4") ])),
            labels = eee$w[i],adj = 1.5,xpd=T,cex= .3+.7*(eee$w[i]/max(eee$w)) )}
  }
  return(eee)
}



plot_enrich <- function(pathway_name,gmt,stats,enrichment_table=NULL,
                        frame=F,axes=F,xlab="",ylab="", main=NULL,cex.main=1,font.main=1,...){
  rnk <- rank(-stats)
  ord <- order(rnk)
  statsAdj <- stats[ord]
  statsAdj <- sign(statsAdj) * (abs(statsAdj)^1)
  statsAdj <- statsAdj/max(abs(statsAdj))

  pathway <- gmt[[pathway_name]]
  pathway <- sort(unname(as.vector(na.omit(match(pathway, names(statsAdj))))))
  pathway <- sort(pathway)
  gseaRes <- calcGseaStat(statsAdj, selectedStats = pathway,
                          returnAllExtremes = TRUE)

  bottoms <- gseaRes$bottoms
  tops <- gseaRes$tops
  n <- length(statsAdj)
  xs <- as.vector(rbind(pathway - 1, pathway))
  ys <- as.vector(rbind(bottoms, tops))
  toPlot <- data.frame(x = c(0, xs, n + 1), y = c(0, ys, 0))
  diff <- (max(tops) - min(bottoms))/8

  plot(toPlot,type="l",col="darkgreen",lwd=2, las=1,
       main=ifelse(is.null(main),pathway_name,main),
       xlim=c(0,length(stats)),xaxs="i",yaxs="i",
       ylim=c(min(toPlot$y),
              max(toPlot$y)),
       frame=frame,axes=axes,xlab=xlab,ylab=ylab,...)
  points(x=pathway,y=rep(0,length(pathway)),pch=73,xpd=T)
  lines(c(pathway[which.max(gseaRes$tops)],
          pathway[which.max(gseaRes$tops)]),
        c(0,gseaRes$res),lty=2,col="grey")

  if(!is.null(enrichment_table)){
    text(length(stats),gseaRes$res,adj=c(1,1),
         labels = paste0("p=",round(enrichment_table$pval[enrichment_table$pathway == pathway_name],5),
                         "\nNES=",round(enrichment_table$NES[enrichment_table$pathway == pathway_name],3),
                         "\nES=",round(enrichment_table$ES[enrichment_table$pathway == pathway_name],3) ) )
  }

  lines( c(0,length(stats)), c(0,0),xpd=T,lwd=1 )
  lines( c(0,0), c(-max(abs(range(toPlot$y))),max(abs(range(toPlot$y)))),xpd=T,lwd=1 )
  text( length(stats)/2, min(toPlot$y) ,labels="gene rank",adj=c(.5,1.5),xpd=T)
  text( 0, mean(range(toPlot$y)) ,labels="ES",srt=90,adj=c(.5,-.5),xpd=T)

}



# CLUSTER DENDROGRAM
# mypar(4,4)
# cluster_means <- as.matrix(sapply(unique(cl$membership),function(x){
#   rowMeans(microbiome[,cl$membership == x]) }))
# colnames(cluster_means) <- unique(cl$membership)
# d <- cor(cluster_means)
# d [ !upper.tri(d) ] <- 0
# image(d[nrow(d):1,] ,col=c("black",colorRampPalette(c("white", "grey80", "firebrick"))(90)),
#       breaks=c(-1,-.0001,seq(0,1,length.out = 90)),frame=F,axes=F)
#
#
#
# clusters_used <- c()
#
# d <- cor(cluster_means)
# d [ !upper.tri(d) ] <- 0
# image(d[nrow(d):1,] ,col=c("black",colorRampPalette(c("white", "grey80", "firebrick"))(90)),
#       breaks=c(-1,-.0001,seq(0,1,length.out = 90)),frame=F,axes=F)
#
# d[which.max(d)] <- -1
# image(d[nrow(d):1,] ,col=c("black",colorRampPalette(c("white", "grey80", "firebrick"))(90)),
#       breaks=c(-1,-.0001,seq(0,1,length.out = 90)),frame=F,axes=F)
# colN <- ceiling(which.max(d) / nrow(d))
# rowN <- ifelse( (which.max(d) %% nrow(d)) == 0 , nrow(d), (which.max(d) %% nrow(d)) )
# clusters_used <- c(clusters_used , colnames(cluster_means)[rowN] , colnames(cluster_means)[colN])
#
# cluster_means <- cbind( cluster_means [,-c(colN,rowN)], rowMeans(cluster_means[,c(colN,rowN)]) )
# colnames(cluster_means)[ncol(cluster_means)] <- paste0(rowN,"_",colN)
#
#

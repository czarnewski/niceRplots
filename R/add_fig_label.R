

add_letter <- function(label="", cex=2, font=2, plot=T,...) {
  totx <- (par("fin")[1] - par("mai")[2] - par("mai")[4])
  totx <- par("usr")[1] - ( (par("usr")[2] - par("usr")[1]) * par("mai")[2] / totx )
  
  toty <- (par("fin")[2] - par("mai")[1] - par("mai")[3])
  toty <- par("usr")[4] + ( (par("usr")[4] - par("usr")[3]) * par("mai")[3] / toty )
  
  if(plot){
    sw <- strwidth(label[1], cex=cex) * 60/100
    sh <- strheight(LETTERS[1], cex=cex) * 60/100
    
    text(label[1], x=totx+sw, y=toty-sh, xpd=T, cex=cex, font=font,...)
    # mtext( text = label[1], at = c(totx+sw, toty-sh), cex=cex, font=font, ... )
  } else{
    return(c(totx,toty))
  }
}




add_arrow <- function(side=1, w=.05, l=.05, invert=F, plot=T,col = "black",...) {
  totx <- (par("fin")[1] - par("mai")[2] - par("mai")[4])
  xmin <- par("usr")[1] - ( (par("usr")[2] - par("usr")[1]) * par("mai")[2] / totx )
  xmax <- par("usr")[2] + ( (par("usr")[2] - par("usr")[1]) * par("mai")[4] / totx )
  xmean <- mean(par("usr")[1:2])
  
  toty <- (par("fin")[2] - par("mai")[1] - par("mai")[3])
  ymin <- par("usr")[3] - ( (par("usr")[4] - par("usr")[3]) * par("mai")[1] / toty )
  ymax <- par("usr")[4] + ( (par("usr")[4] - par("usr")[3]) * par("mai")[3] / toty )
  ymean <- mean(par("usr")[3:4])
  
  if(plot){
    xf <- diff(c(xmin,xmax))
    yf <- diff(c(ymin,ymax))
    if(side == 1){
      xs <- c(xmean-w*xf,xmean+w*xf,xmean)
      ys <- c(ymin,ymin,ymin+l*yf)
    } else if(side == 2){
      ys <- c(ymean-w*yf,ymean+w*yf,ymean)
      xs <- c(xmin,xmin,xmin+l*xf)
    } else if(side == 3){
      xs <- c(xmean-w*xf,xmean+w*xf,xmean)
      ys <- c(ymax,ymax,ymax-l*yf)
    } else if(side == 4){
      ys <- c(ymean-w*yf,ymean+w*yf,ymean)
      xs <- c(xmax,xmax,xmax-l*xf)
    }
    polygon(xs,ys,col = col,xpd=T)
  } 
}











fig_label <- function(text, region="figure", pos="topleft", cex=NULL, ...) {
  
  region <- match.arg(region, c("figure", "plot", "device"))
  pos <- match.arg(pos, c("topleft", "top", "topright", 
                          "left", "center", "right", 
                          "bottomleft", "bottom", "bottomright"))
  
  if(region %in% c("figure", "device")) {
    ds <- dev.size("in")
    # xy coordinates of device corners in user coordinates
    x <- grconvertX(c(0, ds[1]), from="in", to="user")
    y <- grconvertY(c(0, ds[2]), from="in", to="user")
    
    # fragment of the device we use to plot
    if(region == "figure") {
      # account for the fragment of the device that 
      # the figure is using
      fig <- par("fig")
      dx <- (x[2] - x[1])
      dy <- (y[2] - y[1])
      x <- x[1] + dx * fig[1:2]
      y <- y[1] + dy * fig[3:4]
    } 
  }
  
  # much simpler if in plotting region
  if(region == "plot") {
    u <- par("usr")
    x <- u[1:2]
    y <- u[3:4]
  }
  
  sw <- strwidth(text, cex=cex) * 60/100
  sh <- strheight(text, cex=cex) * 60/100
  
  x1 <- switch(pos,
               topleft     =x[1] + sw, 
               left        =x[1] + sw,
               bottomleft  =x[1] + sw,
               top         =(x[1] + x[2])/2,
               center      =(x[1] + x[2])/2,
               bottom      =(x[1] + x[2])/2,
               topright    =x[2] - sw,
               right       =x[2] - sw,
               bottomright =x[2] - sw)
  
  y1 <- switch(pos,
               topleft     =y[2] - sh,
               top         =y[2] - sh,
               topright    =y[2] - sh,
               left        =(y[1] + y[2])/2,
               center      =(y[1] + y[2])/2,
               right       =(y[1] + y[2])/2,
               bottomleft  =y[1] + sh,
               bottom      =y[1] + sh,
               bottomright =y[1] + sh)
  
  old.par <- par(xpd=NA)
  on.exit(par(old.par))
  
  text(x1, y1, text, cex=cex, ...)
  return(invisible(c(x,y)))
}


#https://waterprogramming.wordpress.com/2015/12/02/easy-labels-for-multi-panel-plots-in-r/
#https://bitbucket.org/ggg121/r_figure_letter.git
put.fig.letter <- function(label, location="topleft", x=NULL, y=NULL, 
                           offset=c(0, 0), cex=cex,...) {
  if(length(label) > 1) {
    warning("length(label) > 1, using label[1]")
  }
  
  sw <- strwidth(label[1], cex=cex) * 60/100
  sh <- strheight(label[1], cex=cex) * 60/100
  
  if(is.null(x) | is.null(y)) {
    coords <- switch(location,
                     topleft = c(0,1)            + c(sw,-sh),
                     topcenter = c(0.5,1)        + c(-sw/2,-sh),
                     topright = c(1, 1)          + c(-sw,-sh),
                     bottomleft = c(0, 0)        + c(sw,+sh), 
                     bottomcenter = c(0.5, 0)    + c(-sw/2,+sh), 
                     bottomright = c(1, 0)       + c(-sw,+sh),
                     c(0, 1)                     + c(sw,-sh) )
  } else {
    coords <- c(x,y)
  }
  
  
  this.x <- grconvertX(coords[1] + offset[1], from="nfc", to="user")
  this.y <- grconvertY(coords[2] + offset[2], from="nfc", to="user")
  
  text(labels=label[1], x=this.x, y=this.y, xpd=T,cex=cex, ...)
}
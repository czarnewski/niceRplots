

#' Sets better default values for the base R grafics
#' @export
mypar <- function(
  R = 1,
  C = 1,
  cex.lab = 1,
  cex.axis = 1,
  cex.main = 1,
  font.main = 1,
  mar = c(2.5, 2.5, 1.6, 1.1),
  mgp = c(1.1, 0.3, 0),
  tcl = -0.2,
  las = 1, ...  ){
    par(cex.lab = cex.lab, cex.axis = cex.axis, cex.main = cex.main,
  font.main = font.main, mar = mar, mgp = mgp, tcl = tcl, las = las)
    par(mfrow = c(R, C), ...)
}



#' Prints joint-matrices to the console.
#' @export
LM <- function(
  X = NA,
  R = 1,
  C = 1,
  B = T,
  W = 1,
  H = 1,
  DN = NULL){
  if (is.object(X) || !is.atomic(X))
        X <- as.vector(X)
  X <- as.character(X)
  M <- base::.Internal(matrix(X, nrow = R, ncol = C, B, DN, missing(R),missing(C)))
  attr(M,"W") <- rep(W,ncol(M))[1:ncol(M)]
  attr(M,"H") <- rep(H,nrow(M))[1:nrow(M)]
  return(M)
}


#' Prints joint-matrices to the console.
#' @export
#' @param
#' @example
#' # Create a layout matrix and visualise it
#' A <- LM( X=1:6 , C=2 , W=c(3,4) , H=c(3,2,7) )
#' vis(A)

vis <- function(X){

  attr(X,"H") <- round(attr(X,"H"),2)
  attr(X,"W") <- round(attr(X,"W"),2)

  sprintf("% s", c("1","2") )

  if(nrow(X)==1){
    x <- paste0( c('\n    ',"[",X[1,],"]",'(',attr(X,"H")[1],')' ))
    cat(x)
  } else {
    x <- paste0( c('\n    ',"⎡",X[1,],"⎤",'⎛',attr(X,"H")[1],'⎞' ))
    cat(x)
    if(nrow(X)>2){
      for(i in 2:(nrow(X)-1) ){
        x <- paste0( c('\n    ',"⎥",X[i,],"⎥",'⎥',attr(X,"H")[i],'⎥' ))
        cat(x)
      }
    }
    x <- paste0( c('\n    ',"⎣",X[nrow(X),],"⎦",'⎝',attr(X,"H")[nrow(X)],'⎠' ))
    cat(x)
  }
  x <- paste0( c('\n    ',"(",attr(X,"W"),")\n" ))
  cat(x)
  cat('\n')
}



#' Subsets row, columns and cells matching desired filters. The extracted matrix
#' is expected to be a rectangular matrix. single-row/column matrices are allowed.
#' @param X A layout-matrix object created with `LM()` function.
#' @param pattern A pattern to extract from the layout-matrix `X`. This parameter
#' can be used in conjuction with `rows` and `cols` for more flexible subsetting.
#' @param rows A vector containing the row indexe(s) to extract from the layout-
#' matrix `X`. This parameter #' can be used in conjuction with `pattern` and
#' `cols` for more flexible subsetting.
#' @param cols A vector containing the column indexe(s) to extract from the
#' layout-matrix `X`. This parameter can be used in conjuction with `rows` and
#'  `pattern` for more flexible subsetting.
#' @return This function returns a layout matrix containing the selected pattern,
#' row and/or column.
#' @examples
#' # Create a layout matrix
#' C <- LM( X=1:6 , C=2 , W=c(3,4) , H=c(3,2,7) )
#'
#' # Subset layout matrix with cells containing the value "2"
#' select(X = C, pattern = "2" )
#'
#' # Subset layout matrix with cells containing the value "2" or "4"
#' select(X = C, pattern = c("2","4") )
#'
#' # Subset layout matrix to contain the 2nd row
#' select(X = C, rows = 2:3 )
#'
#' # Subset layout matrix to contain the 1st column
#' select(X = C, cols = 1 )
#'
#' #' # Subset layout matrix to contain the 1st column AND 2nd row
#' select(X = C, cols = 1 , rows = 2)
select <- function(X, rows=NULL, cols=NULL, pattern=NULL, inverse=FALSE){

  M <- X

  if(!is.null(pattern)){
    sel <- matrix(M %in% pattern,nrow = nrow(M),ncol = ncol(M))
    if( inverse ){ sel <- !sel }
    cc <- colSums( (sel)*1 )>0
    rr <- rowSums( (sel)*1 )>0

    M <- matrix( X[sel] , nrow = sum(rr) , ncol = sum(cc) )
    attr(M,"W") <- attr(X,"W")[cc]
    attr(M,"H") <- attr(X,"H")[rr]
  }

  if( !is.null(rows) ) {
    tmp <- M
    M <- matrix( M[rows,] , nrow = length(rows), ncol = ncol(M) )
    attr(M,"W") <- attr(tmp,"W")
    attr(M,"H") <- attr(tmp,"H")[rows]}

  if( !is.null(cols) ) {
    tmp <- M
    M <- matrix( M[,cols] , nrow = nrow(M) , ncol = length(cols) )
    attr(M,"W") <- attr(tmp,"W")[cols]
    attr(M,"H") <- attr(tmp,"H")

  }
  return( M )
}

viewer <- getOption("viewer")
viewer("https://export.uppmax.uu.se/snic2022-23-113/hdca_webdev/03_interactive_svg/")


#' Combines two layout matricies on the matching the pattern value.
#' @export
#' @param M1 The target layout-matrix containing a string matching the pattern
#' specified.
#' @param M2 The layout-matrix to be inserted in the target layout-matrix.
#' @param pattern A string present in the layout-matrix M1, that specifies where
#' the layout-matrix M2 will be inserted.
#' @examples
#' # Create a layout matrix
#' M1 <- LM( X=c(1,"A",2,"A") , C=2 , H=c(1,2) )
#' M2 <- LM( X=1:6 , C=2 , W=c(3,4) , H=c(3,2,6) )
#' vis(M1)
#' vis(M2)
#'
#' # Combine the two layout matrix on index 'A'
#' M3 <- combine( M1 = M1, M2 = C, pattern = "A")
#' vis(M3)
combine <- function( M1 , M2 , pattern = NULL ){
  A <- select(M1, pattern = pattern)

  # Define variables to iterate over
  q0 <- M2
  q1 <- A

  # Expand on the rows (height)
  b0 <- attr(A,"H")
  e0 <- attr(M2,"H")
  bp <- c(0,cumsum(b0)/sum(b0))
  ep <- c(0,cumsum(e0)/sum(e0))
  uu <- sort(unique(c(bp,ep)))[-1]
  o1 <- (1:(length(bp[-1])))[as.numeric(cut(x = uu, breaks = bp))]
  o2 <- (1:(length(ep[-1])))[as.numeric(cut(x = uu, breaks = ep))]
  q0 <- select(X = q0, rows = o2)
  q1 <- select(X = q1, rows = o1)

  # Expand on the columns (widths)
  b1 <- attr(A,"W")
  e1 <- attr(M2,"W")
  bp <- c(0,cumsum(b1)/sum(b1))
  ep <- c(0,cumsum(e1)/sum(e1))
  uu <- sort(unique(c(bp,ep)))[-1]
  o3 <- (1:(length(bp[-1])))[as.numeric(cut(x = uu, breaks = bp))]
  o4 <- (1:(length(ep[-1])))[as.numeric(cut(x = uu, breaks = ep))]
  q2 <- select(X = q0, cols = o4)
  q3 <- select(X = q1, cols = o3)

  M <- B
  sel <- matrix(M %in% pattern,nrow = nrow(M),ncol = ncol(M))
  cc <- colSums( (sel)*1 ) > 0
  rr <- rowSums( (sel)*1 ) > 0
  ci <- 1:length(cc)
  ri <- 1:length(rr)

  ec <- sort(c(ci[!cc],(min(ci[cc])+o3-1)))
  er <- sort(c(ri[!rr],(min(ri[rr])+o1-1)))
  wo <- c(attr(M,"W")[!cc], attr(q2,"W") )
  ho <- c(attr(M,"H")[!rr], attr(q2,"H") )

  # Expand rows
  M <- select(X = B, cols = ec[order(ec)], rows = er[order(er)])
  attr(M,"W") <- wo
  attr(M,"H") <- ho
  M[M == pattern] <- q2

  return(M)
}










B <- LM( X=c(1,"A",2,"A") , C=2 , H=c(1,3) )
vis(B)
C <- LM( X=1:6 , C=2 , W=c(3,4) , H=c(3,2,7) )
vis(C)
combine(M1 = B,M2 = C,pattern = "A")





X <- B
# Expand rows
tmp <- rep(1:nrow(X),rows)
X <- X[tmp,]

# Expand columns
tmp <- rep(1:ncol(X),cols)
X <- X[,tmp]



A


sel <- !grepl("[0-9]",X)
if( sum(sel) > 0){
  for(i in X[sel] ){


}}




CLM <- function(){

}


new_layout <- function(X){
  layout(
    mat     = X,
    widths  = attr(A,"W"),
    heights = attr(A,"H")
  )
}


ELM <- function(X,R,C){
  # Expand rows
  tmp <- rep(1:nrow(X),rows)
  X <- X[tmp,]

  # Expand columns
  tmp <- rep(1:ncol(X),cols)
  X <- X[,tmp]

  attr(X,"W") <- rep(1,ncol(X))
  attr(X,"H") <- rep(1,nrow(X))

  return(X)
}


new_layout( )
layout.show()



A <- LM( X=1:4 , C=2 , W=c(0.2,.35) , H=c(.5,2) )
A
A <- ELM(A)
A
B <- LM( X=c(1,"A",10,11) , C=2 , W=c(.5,1) , H=c(1,1) )
B


X <- B
rows <- round( (attr(X,"H"))/min(attr(X,"H")) )
cols <- round( (attr(X,"W"))/min(attr(X,"W")) )


sel <- !grepl("[0-9]",X)
if( sum(sel) > 0){
  for(i in X[sel] ){
    o <- unique(c(t(X)))
    o <- setNames(1:length(o),o)
    Ns <- as.numeric(grep("[0-9]",names(o),value = T))

    EI <- ELM(get(i))
    rows <- round( (attr(X,"H"))/min(attr(X,"H")) )
    cols <- round( (attr(X,"W"))/min(attr(X,"W")) )
    XX <- (X == i)*1

    rows <- rows * sum(attr(EI,"H"))
    cols <- cols * sum(attr(EI,"W"))
    attr(EI,"W") <- round(attr(EI,"W") * (cols[colSums(XX)>0] / min(attr(EI,"W"))))
    attr(EI,"H") <- round(attr(EI,"H") * (rows[rowSums(XX)>0] / min(attr(EI,"H"))))
    EI <- LM( EI ,R = nrow(EI) ,W = attr(EI,"W"),H = attr(EI,"H"))
    EI <- ELM(EI)

    # rows[rowSums(XX)>0] <- rows[rowSums(XX)>0] * sum(attr(EI,"H"))
    # cols[colSums(XX)>0] <- cols[colSums(XX)>0] * sum(attr(EI,"W"))

    # Expand rows
    tmp <- rep(1:nrow(X),rows)
    X <- X[tmp,]

    # Expand columns
    tmp <- rep(1:ncol(X),cols)
    X <- X[,tmp]

    oi <- unique(c(t(get(i))))
    oi <- setNames(1:length(oi),oi)
    ff <- grep(i,names(o))
    oi <- oi + (ff-1)
    sel <- (o > ff) & (grepl('[0-9]',names(o)))
    if( sum(sel)>0 ){
      nn <- as.numeric( factor(names(o)[sel]) )
      nn <- nn + max(oi)
      o[sel] <- nn
    }
    o[i] <- i
    X
    for(j in rev(names(o)) ){
      X[X == j] <- o[j]
    }
    for(j in rev(names(oi)) ){
      EI[EI == j] <- oi[j]
    }


    X[X==i] <- t(EI)
    X

    X[X==i] <- get(i)
    o <- unique(c(t(X)))
    Ns <- as.numeric(grep("[0-9]",o,value = T))


  }
}














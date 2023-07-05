fix_dates <- function(input){
  tmp <- (as.data.frame(strsplit(gsub( ".* on ","",input)," ")))
  monthz <- setNames( sprintf("%02d",1:12) , month.abb )
  tmp <- paste(tmp[3,],"-",monthz[unlist(tmp[1,])],"-",tmp[2,], sep="")
  return(tmp)
}



p <- function( ... , sep = "", collapse = NULL, recycle0 = FALSE ){
  return(paste(..., sep = "", collapse = collapse, recycle0 = recycle0))
}


p0 <- function( ... , collapse = NULL, recycle0 = FALSE ){
  return(paste0(..., collapse = collapse, recycle0 = recycle0))
}



GEO_read_metadata <- function( GEO_id = NULL , out_dir = "." ){
  
  if(!is.null(GEO_id)){
    # BUILD URL PATH
    url_path <- paste0("ftp://ftp.ncbi.nlm.nih.gov/geo/series/",
                       gsub("...$","nnn",GEO_id),"/",GEO_id,"/matrix/")
    
    # QUERY WHICH FILES ARE AVAILABLE
    fl <- system(paste0("curl --list-only ",url_path), intern = T,ignore.stderr = T )
    
    if(length(fl) > 1){
      message(paste0("GEO acession ",GEO_id," contains ",length(fl)," files.") )
    }
    # FOR EACH FILE
    res <- lapply(fl ,out_dir=out_dir, function(i,out_dir){
      
      # DOWNLOAD
      out_file <- paste0(out_dir,"/",i)
      if(!file.exists( out_file )){ 
        download.file( paste0(url_path,"/",i) , out_file, quiet=T ) } 
      
      # READ SAMPLE INFORMATION
      ii <- readLines( gzfile(out_file,"rt") )
      tmp <- t(read.delim( gzfile(out_file,"rt") , skip = (1:length(ii))[ii==""],header = F))
      colnames(tmp) <- gsub("^[!]","",tmp[1,])
      colnames(tmp) <- gsub("^Sample_","",colnames(tmp))
      rownames(tmp) <- tmp[,"geo_accession"]
      tmp <- tmp[-1,]
      
      # READ PROJECT INFORMATION
      tmp2 <- t(read.delim( gzfile(out_file,"rt") ,header = F, nrows = (1:length(ii))[ii==""]-1 ,sep = "\t"))
      tmp2 <- t(t(setNames(tmp2[2,],tmp2[1,])))
      rownames(tmp2) <- gsub("^!Series_","",rownames(tmp2))
      
      return( list(metadata = data.frame(tmp), projdata = (tmp2) ) )
    }); names(res) <- gsub("_series_matrix.*","",fl)
    
    return( res )
  } else {
    message("Usage:\nmetadata=GEO_read_metadata('GSE134809')")
  }
}




GEO_extract_raw_data <- function( path_to_tar ){
  
  path <- gsub("[.]tar","",path_to_tar)
  if(!dir.exists(path)){ dir.create(path) }
  system(paste0("tar -xvf ",path,".tar -C ",path) )
  
  # Move files into separate directories
  fl <- list.files(path)
  fl <- unique( gsub("_.*","",fl) )
  fl <- paste0(path,"/",fl)
  for(i in fl){  if(!dir.exists( i ) ){ dir.create( i )} }
  
  # Move files into separate directories
  dl <- list.dirs(path)[-1]
  for(i in dl){ system( paste0("mv ",i,"*.* ",i) ) }
}


rename_10x <- function( dl ){
  # Extract files
  for(i in dl){
    fl2 <- list.files(i)
    fl3 <- gsub(".*_","",fl2)
    file.rename(from = paste0(i,"/",fl2),
                to   = paste0(i,"/",fl3))
  }
}


GEO_download_files <- function( meta , columns , dir ){
  future_lapply( rownames(meta), dir_use=dir, meta=meta, columns=columns, function(i,dir_use,meta,columns){
    print(i)
    for(j in c(as.character( meta[i,columns] )) ) {
      dir.create(paste0( dir_use,"/",i), recursive = T )
      download.file( url = paste0(gsub("ftp:","https:",gsub(paste0(i,"_.*"),"",j)),
                                  gsub("_","%5F",gsub("[.]","%2E",gsub(".*[/]","",j)))) ,
                     destfile = paste0( dir_use,"/",i,"/",gsub(".*[/]","",j)) )
    }
  })
}




# 
# 
# GEO_read_metadata <- function( GEO_id = NULL , file_path = NULL, url_path = NULL ){
#   
#   # system("curl --list-only ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE134nnn/GSE134355/matrix/ > aaa.txt")
#   # a <- readLines("aaa.txt")
#   
#   if(!is.null(GEO_id)){
#     url_path <- paste0("https://ftp.ncbi.nlm.nih.gov/geo/series/",
#                        gsub("...$","nnn",GEO_id),"/",GEO_id,"/matrix/",
#                        GEO_id,"_series_matrix.txt.gz")
#     if(!is.null(file_path)) {
#       download.file( url_path , file_path )
#     } else {
#       download.file( url_path , destfile = paste0(GEO_id,"_series_matrix.txt.gz"))
#       file_path <- paste0(GEO_id,"_series_matrix.txt.gz")
#     }
#     ii <- readLines( file_path )
#   } else if(!is.null(file_path)) {
#     
#   } else {
#     file_path <- gzcon(url( url_path ))
#   }
#   # Read sample information
#   ii <- readLines( file_path )
#   tmp <- t(read.delim( file_path , skip = (1:length(ii))[ii==""],header = F))
#   colnames(tmp) <- gsub("^[!]","",tmp[1,])
#   colnames(tmp) <- gsub("^Sample_","",colnames(tmp))
#   rownames(tmp) <- tmp[,"geo_accession"]
#   tmp <- tmp[-1,]
#   
#   # Read project information
#   tmp2 <- t(read.delim( file_path ,header = F, nrows = (1:length(ii))[ii==""]-1 ,sep = "\t"))
#   tmp2 <- t(t(setNames(tmp2[2,],tmp2[1,])))
#   rownames(tmp2) <- gsub("^!Series_","",rownames(tmp2))
#   
#   
#   return( list(metadata = data.frame(tmp), projdata = (tmp2) ) )
# }

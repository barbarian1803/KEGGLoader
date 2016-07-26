trim <- function (x){
  #trim unnecessary blank space
  gsub("^\\s+|\\s+$", "", x)
}

getKEGGRelationCode <- function(relType){
  #kegg relation name to its numeric code
  as.numeric(kegg.relationship.db[kegg.relationship.db$name==relType,"code"])
}

convertGeneID <- function(query,from,to){
  if(from=="ensembl"){
    check <- unique(gene.db[gene.db$ensembl==query,])
    if (nrow(check)==0){
      message(paste("Gene ",query," can not be found in the database, convertion can not be done"))
      return(query)
    }
    return(unique(gene.db[gene.db$ensembl==query,to]))
  }
  if(from=="symbol"){
    check <- unique(gene.db[gene.db$symbol==query,])
    if (nrow(check)==0){
      message(paste("Gene ",query," can not be found in the database, convertion can not be done"))
      return(query)
    }
    return(unique(gene.db[gene.db$symbol==query,to]))
  }
  if(from=="entrez"){
    check <- unique(gene.db[gene.db$entrez==query,])
    if (nrow(check)==0){
      message(paste("Gene ",query," can not be found in the database, convertion can not be done"))
      return(query)
    }
    return(unique(gene.db[gene.db$entrez==query,to]))
  }
}

getTarget <- function(gene,matrix){
  #get gene target for the queried gene 
  output <- c()
  for (i in 1:ncol(matrix)){
    if (matrix[gene,i]>0){
      output <- c(output,colnames(matrix)[i])
    }
  }
  output
}

getTargetSimple <- function(gene,matrix){
  output <- sapply(getTarget(convertGeneID(gene,"symbol","entrez"),matrix),function(x){
    convertGeneID(x,"entrez","symbol")
  },simplify=TRUE)
  output
}

getTF <- function(gene,matrix){
  #get gene TF for the queried gene
  output <- c()
  for (i in 1:nrow(matrix)){
    if (matrix[i,gene]>0){
      output <- c(output,rownames(matrix)[i])
    }
  }
  output
}

getTFSimple <- function(gene,matrix){
  output <- sapply(getTF(convertGeneID(gene,"symbol","entrez"),matrix),function(x){
    convertGeneID(x,"entrez","symbol")
  },simplify=TRUE)
  output
}
source("initLibrary.R")
source("utilityFunction.R")

getPathwayGenes <- function(hsa_code){
  #get list of genes by its symbol for a specified pathway
  data <- getURL(paste("http://togows.dbcls.jp/entry/pathway/",hsa_code,"/genes.json",sep=""))
  data <- t(fromJSON(data))
  entrezID <- row.names(data)
  gene.db[gene.db$EntrezGene.ID%in%entrezID,]
}

loadKGML <- function(hsa_code){
  #download KGML file from the server
  message("Get KGML file from the server")
  dataXML <- xmlParse(paste("http://rest.kegg.jp/get/",hsa_code,"/kgml",sep=""))
  dataXML <- xmlToList(dataXML)
  dataXML
}

KGMLParser <- function(KGML){
  #parse KGML to get the data, only select entry and relationship type from the KGML
  #there are several tag in KGML data format. The network information is stored in entry and relationship tag.
  #entry defined vertices and relationship defined edges
  #entry in original KEGG network doesn't represent a single gene. It represents a group or a fimily of gene
  #each entry is given an entry ID and this entry ID is used in relationship data to connect 2 entries in network
  message("get the nodes from the pathway")
  
  entry <- list()
  idxEntry<-1
  
  relation <- list()
  idxRelation<-1
  
  for(i in 1:length(KGML)){
    if(names(KGML[i])=="entry"){
      entry[[idxEntry]] <- KGML[[i]]
      idxEntry <- idxEntry+1
    }else if(names(KGML[i])=="relation"){
      relation[[idxRelation]] <- KGML[[i]]
      idxRelation <- idxRelation+1
    }
  }
  
  
  entry <- lapply(entry,function(x)x[".attrs"])
  entry <- lapply(entry,function(x) as.data.frame(t(x[[1]]),stringsAsFactors=FALSE))
  entry <- lapply(entry,function(x){
    if(ncol(x)==4){
      x
    }
  })
  
  entry <- do.call(rbind,entry)
  entry <- entry[entry$type=="gene",]
  list("entry"=entry[entry$type=="gene",],"relationship"=relationshipDataFormatter(relation))
}

relationshipDataFormatter <- function(relation){
  #format the relationship data extracted from KGML as a table
  #data extracted from KGML is not formatted. It is just parsed string by R function
  #this function reformat it to table format
  message("formatting relationship data")
  relationship <- c()
  
  for(xx in 1:length(relation)){
    innerList <- relation[[xx]]
    #innerList consist of subtype and .attrs
    c <-c()
    attr <- NULL
    
    for(j in 1:length(innerList)){
      if(names(innerList[j])=="subtype"){
        c<-rbind(c,innerList[[j]])
      }else{
        attr <- as.data.frame(t(innerList[[j]]))
      }
    }
    relationship <- rbind(relationship,merge(c,attr))
  }
  relationship
}

createRelationshipTable <- function(entry,relationship){
  #create gene-gene relationnship
  #it extracts the gene information from entry table
  #use relationship data based on its entry ID to make gene-gene relationship
  relationshipTable <- data.frame()
  for (idx in 1:nrow(relationship)){
    row <- relationship[idx,]
    
    relType <- as.character(levels(relationship$name)[relationship[idx,"name"]])
    en1 <- as.numeric(levels(relationship$entry1)[relationship[idx,"entry1"]])
    en2 <- as.numeric(levels(relationship$entry2)[relationship[idx,"entry2"]])
    
    entry1 <- entry[entry$id==en1,"name"]
    entry2 <- entry[entry$id==en2,"name"]
    
    if ((length(entry1) == 0) && (typeof(entry1) == "character")){
      next
    }
    if ((length(entry2) == 0) && (typeof(entry2) == "character")){
      next
    }
    
    entry1genes <- strsplit(entry1, " ")[[1]]
    entry2genes <- strsplit(entry2, " ")[[1]]
    
    
    for (gene1 in entry1genes){
      gene1 <- substring(gene1,5)
      for (gene2 in entry2genes){
        gene2 <- substring(gene2,5)
        newRow <- data.frame("relType"=relType,"entry1"=en1,"entry2"=en2,stringsAsFactors = FALSE)
        newRow$gene1 <- gene1
        newRow$gene2 <- gene2
        relationshipTable <- rbind(relationshipTable,newRow)
      }
    }
  }
  relationshipTable
}

createAdjacencyMatrix <- function(relationshipTable,distinct=TRUE){
  #change the format of relationship table into adjacency matrix
  #because there is a possibilities for a gene has 2 entry ID, the entry ID is attached to gene name to identify it
  #it happens if in one KGML, there are several separated network, for example : alternative pathway
  #the matrix result is : column == target gene, row == TF gene
  if(distinct){
    element <- unique(c(paste(relationshipTable$gene1,relationshipTable$entry1,sep="."),paste(relationshipTable$gene2,relationshipTable$entry2,sep=".")))
  }else{
    element <- unique(c(relationshipTable$gene1,relationshipTable$gene2))
  }
  adjMat <- matrix(0, nrow = length(element), ncol = length(element))
  colnames(adjMat) <- element
  rownames(adjMat) <- element
  
  for (i in 1:nrow(relationshipTable)){
    r <- relationshipTable[i,]
    rel <- getKEGGRelationCode(r$relType)
    if(distinct){
      adjMat[paste(r$gene1,r$entry1,sep="."),paste(r$gene2,r$entry2,sep=".")] <- rel
    }else{
      adjMat[r$gene1,r$gene2] <- rel
    }
  }
  adjMat
}

tableToAdjmatrix <- function(table){
  element <- union(table$Target,table$TF)
  adjMat <- matrix(0, nrow = length(element), ncol = length(element))
  colnames(adjMat) <- element
  rownames(adjMat) <- element
  
  for (i in 1:nrow(table)){
    adjMat[table[i,"TF"],table[i,"Target"]] <- table[i,"relType"]
  }
  
  adjMat
}

setGeneSymbolAdjMatrix <- function(list){
  for (i in 1:length(list)){
    for(x in 1:(length(colnames(list[[i]])))){
      colnames(list[[i]])[x] <- convertGeneID(unlist(strsplit(colnames(list[[i]])[x],"[.]"))[1],"entrez","symbol")
    }
    
    for(x in 1:(length(rownames(list[[i]])))){
      rownames(list[[i]])[x] <- convertGeneID(unlist(strsplit(rownames(list[[i]])[x],"[.]"))[1],"entrez","symbol")
    }
  }
  return(list)
}

separateAdjMatrix <- function(matrix){
  #for adj matrix which contains more than one pathway,
  #separate it so that the adj matrix consist only 1 pathway
  #output is a list of matrix
  dfs <- function(gene,matrix){
    evaluated <<- c(evaluated,gene)
    TFs <- getTF(gene,matrix)
    if(!is.null(TFs)){
      for(tf in TFs){
        if(nrow(tesDFS[tesDFS$TF==tf & tesDFS$Target==gene,])==0){
          tesDFS <<- rbind(tesDFS,data.frame("TF"=tf,"Target"=gene,"relType"=matrix[tf,gene],stringsAsFactors = FALSE))
        }
      }
      for(tf in TFs){
        if(!tf%in%evaluated){
          dfs(tf,matrix)
        }
      }
    }
  }
  
  resList <- list()
  
  for(gene in getEndNode(matrix)){
    tesDFS <- data.frame()
    evaluated <- c()
    dfs(gene,matrix)
    resList[[gene]]<-tesDFS
  }
  
  resListBak <- resList
  mergedRes <- list()
  cur <- 1
  while(length(resListBak)>0){
    mergedRes[[cur]] <- resListBak[[1]]
    resListBak <- resListBak[-1]
    if(length(resListBak)==0){
      break
    }
    merged <- c()
    for(i in 1:length(resListBak)){
      temp <- union(mergedRes[[cur]],resListBak[[i]])
      if(nrow(temp)<nrow(mergedRes[[cur]])+nrow(resListBak[[i]])){
        mergedRes[[cur]] <- temp
        merged <- c(merged,i)
      }
    }
    
    if(!is.null(merged)){
      resListBak <- resListBak[-c(merged)]
    }
    
    cur <- cur+1
  }
  
  for(i in 1:length(mergedRes)){
    mergedRes[[i]]<-tableToAdjmatrix(mergedRes[[i]])
  }
  return(setGeneSymbolAdjMatrix(mergedRes))
}

getPathwayRelationshipTable<-function(hsa_code){
  #get KEGG pathway network in the form of relationship table
  KGML <- loadKGML(hsa_code)
  parsedKGML <- KGMLParser(KGML)
  entry <- parsedKGML[["entry"]]
  relationship <- parsedKGML[["relationship"]]
  relationshipTable <- createRelationshipTable(entry,relationship)
}

getEndNode <- function(matrix){
  #get the end node of the network represents in the adj matrix
  output <- c()
  for (g in colnames(matrix)){
    if (is.null(getTarget(g,matrix))){
      output <- c(output,g)
    }
  }
  output
}
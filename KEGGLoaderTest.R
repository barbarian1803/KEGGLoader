source("KEGGLoader.R")

wnt <- getPathwayRelationshipTable("hsa04310")
ajd.mat.wnt <- createAdjacencyMatrix(wnt,FALSE)
ajd.mat.wnt2 <- createAdjacencyMatrix(wnt,TRUE)


tes <- sapply(getTarget(convertGeneID("PORCN","symbol","entrez"),ajd.mat.wnt),function(x){
  convertGeneID(x,"entrez","symbol")
},simplify=TRUE)


resList <- separateAdjMatrix(ajd.mat.wnt2)


for(i in 1:length(resList)){
  print(nrow(resList[[i]]))
}




length(union(resList[[3]]$TF,resList[[3]]$Target))
for(i in 1:nrow(resList[[3]])){
  print(resList[[3]][1,])
  
  g1 <- unlist(strsplit(resList[[3]][1,"Target"],"[.]"))[1]
  g2 <- unlist(strsplit(resList[[3]][1,"TF"],"[.]"))[1]
  print(convertGeneID(g1,"entrez","symbol"))
  print(convertGeneID(g2,"entrez","symbol"))
  break
}


for(x in 1:(length(colnames(resList[[3]])))){
  colnames(resList[[3]])[x] <- convertGeneID(unlist(strsplit(colnames(resList[[3]])[x],"[.]"))[1],"entrez","symbol")
}

for(x in 1:(length(rownames(resList[[3]])))){
  rownames(resList[[3]])[x] <- convertGeneID(unlist(strsplit(rownames(resList[[3]])[x],"[.]"))[1],"entrez","symbol")
}

tes <- setGeneSymbolAdjMatrix(resList)

View(resList[[5]])

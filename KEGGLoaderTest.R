wnt <- getPathwayRelationshipTable("hsa04310")
ajd.mat.wnt <- createAdjacencyMatrix(wnt,FALSE)
ajd.mat.wnt2 <- createAdjacencyMatrix(wnt,TRUE)


tes <- sapply(getTarget(convertGeneID("PORCN","symbol","entrez"),ajd.mat.wnt),function(x){
  convertGeneID(x,"entrez","symbol")
},simplify=TRUE)


resList <- separateAdjMatrix(ajd.mat.wnt2)
for(x in names(resList)){
  print(paste(x,":",nrow(resList[[x]])))
}
xxx <- union(resList[["4609.36"]],resList[["3725.35"]])
xxx <- union(xxx,resList[["8061.34"]])

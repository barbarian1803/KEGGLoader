ProcessTF <- function(Target,graph){
  message(paste("Processing gene ",Target))
  allTF <- getTF(Target,graph)
  inhibitor <- c()
  activator <- c()
  unknown <- c()
  
  for(tf in allTF){
    tf.id <- convertGeneID(tf,"symbol","ensembl")
    if(identical(tf.id,character(0))){
      next
    }
    if(tf.id==tf){
      next
    }
    if(abs(logFCData[tf.id,"log2FoldChange"])<0.65){
      next
    }
    type <- getTFType(tf,Target,graph)
    
    if(type%in%c(3,5,9)){
      activator <- c(activator,tf)
    }else if(type%in%c(4,6,10)){
      inhibitor <- c(inhibitor,tf)
    }else{
      unknown <- c(unknown,tf)
    }
    
  }
  inhibitorAnalysis <- checkInhibitor(inhibitor,Target)
  activatorAnalysis <- checkActivator(activator,Target)
  output <- list(inhibitorAnalysis,activatorAnalysis)
  output
}

checkInhibitor <- function(inhibitor,target){
  Target <- convertGeneID(target,"symbol","ensembl")
  logTarget <- logFCData[Target,"log2FoldChange"][1]
  
  significant <- c()
  nonsignificant<-c()
  
  for(gene in inhibitor){
    Gene <- convertGeneID(gene,"symbol","ensembl")
    corr <- cor(AssayNormalCancer[Gene,],AssayNormalCancer[Target,])
    if(corr < (-0.6)){
      significant <- c(significant,gene)
    }else{
      nonsignificant <- c(nonsignificant,gene)      
    }
  }
  
  cause <- list()
  cause[[1]]<-significant
  cause[[2]]<-nonsignificant
  cause
}

checkActivator <- function(activator,target){
  Target <- convertGeneID(target,"symbol","ensembl")
  
  oriactivator <- activator
  activator <- filterActivator(activator,target)
  
  if(length(activator)<1){
    res <- list()
    res[[1]]<-c()
    res[[2]]<-oriactivator
    res[[3]]<-NULL
    return(res)
  }
  activatorEnsembl <- c()
  for (a in activator){
    activatorEnsembl <- c(activatorEnsembl,convertGeneID(a,"symbol","ensembl"))
  }
  activator <- activatorEnsembl
  max <- -1
  selected <- c()
  LMSummary <- NULL
  for(i in length(activator):1){
    list <- enum.choose(activator,i)
    for(j in 1:length(list)){
      formula <- paste(Target,paste(list[[j]],collapse="+"),sep="~")
      formula <- as.formula(formula)
      LM <- summary(lm(formula,as.data.frame(t(AssayNormalCancer))))
      rsqr <- LM$adj.r.squared
      coeff <- LM$coefficients
      if( rsqr>max & !AnyNegativeCoeff(coeff[-1,1]) & AnySignificant(coeff[-1,4]) ){
        max <- rsqr
        selected <- list[[j]]
        LMSummary <- LM
      }
    }
  }
  selectedSymbol <- c()
  for(n in selected){
    selectedSymbol <- c(selectedSymbol,convertGeneID(n,"ensembl","symbol"))
  }
  res <- list()
  res[[1]] <- selectedSymbol
  res[[2]]<-setdiff(oriactivator,selectedSymbol)
  res[[3]]<-LMSummary
  res
}

#remove activator which has negative correlation with its target
filterActivator <- function(activator,target){
  out <- c()
  Target <- convertGeneID(target,"symbol","ensembl")
  for( gene in activator){
    act <- convertGeneID(gene,"symbol","ensembl")
    corr <- cor(AssayNormalCancer[act,],AssayNormalCancer[Target,])
    if(corr > 0){
      out <- c(out,gene)
    }
  }
  out
}

enum.choose <- function(x, k) {
  if(k > length(x)) stop('k > length(x)')
  if(choose(length(x), k)==1){
    list(as.vector(combn(x, k)))
  } else {
    cbn <- combn(x, k)
    lapply(seq(ncol(cbn)), function(i) cbn[,i])
  }
}

AnyNegativeCoeff <- function(vector){
  out <- FALSE
  for(v in vector){
    if(v<0){
      out <- TRUE
      break
    }
  }
  out
}

AnySignificant <- function(vector){
  out <- FALSE
  for(v in vector){
    
    if(!is.nan(v)&v<=0.05){
      out <- TRUE
      break
    }
  }
  out
}
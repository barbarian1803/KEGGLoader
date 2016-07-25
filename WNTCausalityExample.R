load("WNTGraph.dat")
load("DESeq2resultBDC.dat")
source("KEGGLoader.R")
source("CausalityAnalysisWorkflow.R")

logFCData <- FC_BDC_Normal
ProcessTF("LEF1",canonWNTGraph)


all <- list()
for(g in colnames(canonWNTGraph)){
  ensembl <- convertGeneID(g,"symbol","ensembl") 
  if(ensembl==g || abs(logFCData[ensembl,"log2FoldChange"]<0.65)){
    message(paste("Gene ",g," is skipped"))
    next
  }
  all[[g]] <- ProcessTF(g,canonWNTGraph)
}

ProcessTF("WNT10B",canonWNTGraph)

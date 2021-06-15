search_rna_go <- function(RNA){
  #circRNA <- match.arg(RNA)
  load("data/rna2gene_GO_name.rda")
  if(! RNA %in% rna2gene_GO_name$RNA_ID){
    cat("The RNA has no corresponding GO function information ")
    return("NA")
  }
  infer_go <- rna2gene_GO_name[which(rna2gene_GO_name$RNA_ID == RNA),]
  return(infer_go)
}


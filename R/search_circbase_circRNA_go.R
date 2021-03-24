search_circbase_circRNA_go <- function(circRNA){
  #circRNA <- match.arg(circRNA)
  load("data/circ_gene_go.rda")
  if(! circRNA %in% circ_gene_go$circRNA_ID){
    cat("The circRNA has no corresponding GO function information ")
    return("NA")
  }
  infer_go <- circ_gene_go[which(circ_gene_go$circRNA_ID == circRNA),
                           c("circRNA_ID", "gene_symbol","Ensembl_ID","GO_term_accession",
                             "GO_term_name","GO_term_evidence_code", "GO_domain")]
  return(infer_go)
}


#' @title Run module search function
#'
#' @description runmodule constructs a node-weighted ncRNA network, performs module searching, generates simulation data from random networks,
#' @description normalizes module scores using simulation data, removes un-qualified modules, and orders resultant modules according to their significance.
#'
#' @param network A data frame containing a symbolic edge list of the ncRNA network in which the columns must contain "node_gene_ID", "type", "target_gene_Name"
#'
#' @param gene2weight  A weigth data frame containing three columns:"type","gene", "weight"
#'  the first "type" the type of the gene identifier; lncRNA, miRNA, circRNA
#'  the second gene is unique, gene identifier (should be coordinate with the node symbol used in ncRNA network);
#'  the third weight is gene-based  p-value or corrected p-value derived from differentially gene analysis or survival analysis
#' @param method	a character string indicating which the search method is to be computed
#'  . One of "global" (default, refer to Heinz method), "local ( refer to GS method)": can be abbreviated
#' @param d An integer used to define the order of neighbour genes to be
#'   searched  in  the method of  the method "local" . This parameter is default set up as 2
#' @param r A float indicating the cut-off for increment during module expanding
#'   process in  the method of  the method "local". Greater r will generate smaller module. Default is 0.1.
#' @param  maxsize An integer: the numbel of size of the module for user settings in the method of "global", default 15.
#' @param FDR Numeric value, from the false discovery rate a p-value threshold is calculated. P-values below this threshold are considered to be significant
#'  The FDR can be used to control the size of the maximum scoring module
#' @param seletN a vector: gene identifier IDs, or a gene identifier ID, for example "MIMAT0000461",c("MIMAT0000461", "ENSG00000250742")
#' @param issymbol Boolean value, whether to set the node attribute "symbol"(gene symbol) in the network.
#'  @return \code{runmodule} returns a list containing relevant data and results,
#'   including:
#'  \tabular{ll}{
#'    \code{GNCW} \tab the node-weighted network used for searching \cr
#'    \code{module} \tab list of genes comprising each module, named for the seed gene if the method is "local" or the igraph class of the module if the method is "global" \cr
#'    \code{module.score.matrix} \tab contains Zm, Zn \cr
#'  }
#'
#' @references Hongbo Shi, Jiayao Li, Qiong Song et al. (2019) Systematic identification and analysis of dysregulated miRNA and transcription factor feed-forward loops in hypertrophic cardiomyopathy
#' @references Peilin Jia, Siyuan Zheng, Jirong Rong, Wei Zheng, Zhongming Zhao. (2011) Bioformatics. dmGWAS: dense module searching for genome-wide association studies in protein-protein interaction networks.
#' @references Daniela Beisser, Gunnar W. Klau, Thomas Dandekar  et al. (2019)  BioNet: an R-Package for the functional analysis of biological networks
#' @examples
#' \dontrun{
#' data("dataN")
#' gene2weight <- combinp(dataN[,c("type","logFC","PValue")])
#' interac <- interStringency(type = "ncRNA",stringency = "strict")
#' interac <- interac[,c("node_gene_ID","type","target_gene_ID")]
#'  res.list_global <- runmodule(network = interac, gene2weight, method = "global",FDR = 1e-14)
#'  res.list_local <- runmodule(network = interac, gene2weight, method = "local",
#'  maxsize=15, seletN = "MIMAT0000461")
#' }
#'
#' @export



runmodule <- function (network, gene2weight, maxsize = 15, method = c("global","local"), d = 2, r = 0.1, seletN = NULL,FDR = 1e-14,issymbol = TRUE) #     c("Heinz","GS") infer to  c("global","local")
{
    #library(igraph)
    #library(stats)
    maxs <- maxsize
    if (min(gene2weight[, c("weight")]) <= 0 | max(gene2weight[, c("weight")]) >=
        1) {
        stop("P values out of range 0 < p < 1")
    }
   if(!missing(method) & length(method)>1) stop("Only one 'method' allowed.")
    method <- match.arg(method)
    cat("genes used: ", length(gene2weight[, c("gene")]), "\n", sep = "")
    rawG <- graph.data.frame(network[,c("node_gene_ID","target_gene_ID")], directed = F)
    rawG <- set_edge_attr(rawG, "types", value= as.character(network$type))
    rawG <- simplify(rawG,edge.attr.comb=toString)


  if(method == "local"){
    g.weight <- sapply(as.numeric(gene2weight[, c("weight")]), function(x) qnorm(1 -
        x))
  }
  if(method == "global"){
    pvals <- gene2weight$weight
    names(pvals) <- rownames(gene2weight)
    fb <- fitBumModel(pvals)
    g.weight <- scoreFunction(fb, fdr = FDR)
    }
    intG <- integGM(rawG, as.character(gene2weight[, c("gene")]), g.weight, as.character(gene2weight[, c("type")]), issymbol = issymbol )#, as.character(gene2weight[, c("symbol")])
    GNCW <- simplify(intG, edge.attr.comb=toString)
    cat("start searching at ", format(Sys.time(), "%H:%M, %b %d %Y"),
        " ...\n", sep = "")
	if( method == "global"){
		scores <- V(GNCW)$weight
		names(scores)<- V(GNCW)$name
		module <- FastHeinz(GNCW, scores)
		genes.idx <- seq(1, length(V(GNCW)$name))
		graph.g.weight = data.frame(GNCWene = V(GNCW)$name, gain.weight = V(GNCW)$weight)
		genes <- V(module)$name
		idx <- match(genes, graph.g.weight[, 1])
		idx <- idx[!is.na(idx)]
		ZM <- sum(graph.g.weight[idx, 2])/sqrt(length(idx))
		l.zperm <- c()
		for (j in 1:1e+05) {
            idx.pseudo = sample(genes.idx, size = length(V(module)))
            l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo,
                2])/sqrt(length(idx.pseudo)))
        }
		k.mean <- mean(l.zperm)
		k.sd <- sd(l.zperm)
		ZN = (ZM - k.mean)/k.sd
		res.list <- list()
		res.list[["GNCW"]] = GNCW
		res.list[["module"]] = module
		res.list[["module.score.matrix"]] = data.frame(Zm=ZM,Zn=ZN)
		cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"),
        " ...\n", sep = "")
		return(res.list)

	}
  if( method == "local"){
	sublist = list()
  nodes_user <- seletN
  if (is.null(nodes_user)){
    cat("please input Nodes "," ...\n", sep = "")
  }
	for (node in nodes_user) {#V(GNCW)$name
    if(!node %in% V(GNCW)$name){
     cat("Error: the Node  ",node," have no expression value or the interaction", sep = "")
      break
   }
		ng <- seedQuery_GS(GNCW, node, maxsize = maxs, search_r = 2, r = 0.1)
		if (vcount(ng) >= 5)
			sublist[[node]] <- ng
	}

	dm.result <- sublist
    cat("extracting modules...\n", sep = "")
    genesets <- list()
    for (k in 1:length(dm.result)) {
        node = names(dm.result[k])
        g = dm.result[[k]]
        genesets[[node]] <- V(g)$name
    }
    seed.genes <- names(genesets)
    #cat("removing identical modules...\n", sep = "")
    identical.idx <- list()
    for (k in 1:length(genesets)) {
        tmp.idx <- c(k)
        for (kt in 1:length(genesets)) {
            if (kt == k)
                (next)()
            genesk <- genesets[[k]]
            genest <- genesets[[kt]]
            if (length(genesk) != length(genest))
                (next)()
            overlap = intersect(genesk, genest)
            if (length(overlap) == length(genest)) {
                tmp.idx <- c(tmp.idx, kt)
            }
        }
        if (length(tmp.idx) > 1) {
            tmp.idx <- sort(tmp.idx)
            identical.idx[[seed.genes[k]]] <- tmp.idx
        }
    }
    toremove.idx <- c()
    for (k in 1:length(identical.idx)) {
        if(length(identical.idx) == 0)
          break
        tmp.idx <- identical.idx[[k]]
        toremove.idx <- c(toremove.idx, tmp.idx[-1])
    }
    toremove.idx <- unique(toremove.idx)
    if(is.null(toremove.idx)){
      genesets.clear <- genesets
    }else {
    genesets.clear <- genesets[-toremove.idx]
    }
    cat("permutation on random network...\n", sep = "")
    genesets.length <- c()
    for (k in 1:length(genesets.clear)) {
        genes <- genesets.clear[[k]]
        genesets.length <- c(genesets.length, length(genes))
    }

    genesets.length <- unique(genesets.length)
    genes.idx <- seq(1, length(V(GNCW)$name))
    graph.g.weight = data.frame(GNCWene = V(GNCW)$name, gain.weight = V(GNCW)$weight)
    genesets.length.null.dis <- list()
    length.max = max(genesets.length)
	length.min = min(genesets.length)
    for (k in length.min:length.max) {
        l.zperm <- c()
        for (j in 1:1e+05) {
            idx.pseudo = sample(genes.idx, size = k)
            l.zperm <- c(l.zperm, sum(graph.g.weight[idx.pseudo,
                2])/sqrt(length(idx.pseudo)))
        }
        genesets.length.null.dis[[as.character(k)]] = l.zperm
        cat(k, ".", sep = "")
    }
    genesets.length.null.stat <- list()
    for (k in length.min:length.max) {
        l.zperm <- genesets.length.null.dis[[as.character(k)]]
        k.mean <- mean(l.zperm)
        k.sd <- sd(l.zperm)
        genesets.length.null.stat[[as.character(k)]] = c(k.mean,
            k.sd)
    }
    zim <- data.frame(gene = names(genesets.clear), Zm = -9,
        Zn = -9, zcount = -9)
    for (k in 1:length(genesets.clear)) {
        genes <- genesets.clear[[k]]
        idx <- match(genes, graph.g.weight[, 1])
        idx <- idx[!is.na(idx)]
        zim[k, 2] <- sum(graph.g.weight[idx, 2])/sqrt(length(idx))
        tmp <- genesets.length.null.stat[[as.character(length(idx))]]
        zim[k, 3] = (zim[k, 2] - tmp[1])/tmp[2]
        zim[k, 4] = sum(genesets.length.null.dis[[as.character(length(idx))]] >=
            zim[k, 2])
    }
    zom = zim[order(zim[, 3], decreasing = T), ]
    res.list <- list()
    res.list[["GNCW"]] = GNCW
    #res.list[["graph.g.weight"]] = graph.g.weight
    res.list[["module"]] = genesets.clear
    #res.list[["genesets.length.null.dis"]] = genesets.length.null.dis
    #res.list[["genesets.length.null.stat"]] = genesets.length.null.stat
    #res.list[["zi.matrix"]] = zim
    res.list[["module.score.matrix"]] = zom[,c("gene", "Zm", "Zn")]
    save(res.list, file = "RESULT.list.RData")
    cat("finished at ", format(Sys.time(), "%H:%M, %b %d %Y"),
        " ...\n", sep = "")
    return(res.list)
	}
}

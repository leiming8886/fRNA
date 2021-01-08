#' @title Extract interactions according to stringency and interaction type
#'
#' @description interactions were extracted according to stringency and interaction type in the database of ENCORI

#' @param type	a character string indicating which interaction type is to be choosed,
#'  . One of "RBP" (RBP-circRNA,RBP-lncRNA,miRNA-circRNA,miRNA-lncRNA,miRNA-RBP), "ncRNA (miRNA-circRNA,miRNA-lncRNA)": can be abbreviated
#' @param stringency	a character string indicating which interaction stringency is to be choosed,
#'  . One of "low" ( number of supported experiments > = 1  ), "medium ( > = 2)","high ( > = 3)","strict ( > = 5)"
#'
#' @return interaction of setting
#' @examples
#' \dontrun{
#'   data("dataM2C")
#' data("dataM2L")
#' data("dataM2R")
#' data("dataR2C")
#' data("dataR2L")
#'  interac <- interStringency(type = "ncRNA",stringency = "strict")
#'  interac <- interac[,c("node_gene_ID","target_gene_ID")]
#' }
#'
#' @export
interStringency <- function(type = c("RBP","ncRNA"),stringency = c("low","medium","high","strict") ){
  #data(dataM2C)
  #data(dataM2L)
  #data(dataM2R)
  #data(dataR2C)
  #data(dataR2L)
  type <- match.arg(type)
  interaction <- NULL
  switch( type,
          "RBP"={
            interaction <-rbind(fRNC::dataR2C,fRNC::dataR2L)
            interaction <-rbind(interaction,fRNC::dataM2R)
            interaction <-rbind(interaction,fRNC::dataM2L)
            interaction <-rbind(interaction,fRNC::dataM2C)
          },
          "ncRNA"={
            interaction <-rbind(fRNC::dataM2C,fRNC::dataM2L)
            interaction <-rbind(interaction,fRNC::dataM2R)
          },

          stop("Enter something that switches me!")
  )


	if ( ! "clipExpNum" %in% colnames(interaction)){
		stop("clipExpNum must be in the col names!\n")
	}
	stringency <- match.arg(stringency)

	switch( stringency,
        "low"={
            interaction1 <- interaction[which(interaction$clipExpNum>=1),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
        },
        "medium"={ interaction1 <- interaction[which(interaction$clipExpNum>=2),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
        },
		     "high"={ interaction1 <- interaction[which(interaction$clipExpNum>=3),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
        },
		    "strict"={ interaction1 <- interaction[which(interaction$clipExpNum>=5),c("node_gene_ID", "node_gene_Name", "type", "target_gene_ID", "target_gene_Name")]
        },
        stop("Enter something that switches me!")
    )

	return(interaction1)

}

##fRNC: A R package to uncover RBP-ncRNA circuits from multi-omics data


##Description

The RNA binding protein (RBP) and non-coding RNA (ncRNA) interacting networks are increasingly recognized as main
mechanism in post-transcriptional regulation, and tightly associated with cellular malfunction and disease. Here,
we present fRNC, a systems biology tool to uncover dynamic spectrum of RBP-ncRNA circuits (RNC) by integrating 
transcriptomics, interactomics and clinical data. fRNC constructs the RBP-ncRNA network from experiment derived
CLIP-seq or PARE data. Given scoring on each node in the network, it finds a RNC containing global maximum significant 
genes. Alternatively, it can also search locally maximum RNCs according to user defined nodes. It enables users flexibly
to analyse and visualize the collective behaviors between a RBP and its interacting ncRNAs in a malfunctioned biological process.

##Getting Started
###Step 1. Install package dependencies 
     install_dependpackages <- function(){
         metr_pkgs <- c("survminer", "limma", "ggpubr", "XML", "igraph", "multtest","RBGL","edgeR")  
         list_installed <- installed.packages()
         new_pkgs <- subset(metr_pkgs, !(metr_pkgs %in% list_installed[, "Package"])) 
         if(length(new_pkgs)!=0){   
            if (!requireNamespace("BiocManager", quietly = TRUE))
                install.packages("BiocManager")
                BiocManager::install(new_pkgs)
                print(c(new_pkgs, " packages will be installed..."))
            }  
          if((length(new_pkgs)<1)){
                print("No new packageswill be installed...")
            }
        }
     install_dependpackages()

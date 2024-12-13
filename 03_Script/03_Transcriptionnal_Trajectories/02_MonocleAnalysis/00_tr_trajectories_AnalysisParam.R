# Enter the parameters of the analysis for tr_trajectories analysis


palette<-c(C1_ILCP = "#E41A1C", 
           "C2_Klrd1+" = "#377EB8", C3_LTiP = "#4DAF4A",
           C4_EILP = "#984EA3")

#table(NK_Seurat$seurat_clusters)


CLUSTERS_TO_REMOVE = c("C5_Apoe+", "C3_LTiP")

LIST_PARAM_CONTROL =list(ncenter=150)
root_pr_nodes = "Y_26"

#For graph test:
#Q_Value_Limit = 0.00005

#Coefficients
#k=9


#Heatmap pseudotime parameters
HEIGHT = unit(15, "cm")
WIDTH = unit(12, "cm")
POLICE= gpar(fontsize = 3)

NUMBER_CLUSTS = 6
NUMBER_GENES_BIG_HEATMAP = 150
NUMBER_GENES_SMALL_HEATMAP = 100

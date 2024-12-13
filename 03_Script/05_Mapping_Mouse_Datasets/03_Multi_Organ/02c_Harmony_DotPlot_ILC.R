## This script visualizes the expression of genes associated with different innate lymphoid cell (ILC) subsets across the integrated dataset.
## It begins by defining gene lists for each ILC subset (ILC1, ILC2, ILC3, ILCP, CLP, and all ILCs combined) and assigning colors to these groups.
## A DotPlot is created to display the expression levels of these genes, with row names color-coded according to their respective ILC groups.
## The plot highlights the ILC subset signatures within the dataset, providing a clear overview of their gene expression patterns.
## The resulting DotPlot is saved as both PNG and PDF files for further analysis or presentation.

ILC_All_List = c( "Id2" , "Ets1" , "Rora" , "Ahr" , "Cxcr6", "Arg1" , "Il7r" , "Il2rb" )
ILC1_List = c( "Tbx21" , "Eomes" , "Ccl5" , "Klrk1" , "Ncr1" , "Cxcr3" , "Xcl1" , "Ifng" )
ILC2_List = c( "Gata3" , "Bcl11b" , "Il1rl1" , "Il17rb" , "Icos" , "Hey1" , "Il9r" , "Il4" )
ILC3_List = c( "Rorc" , "Batf3" , "Batf", "Tox2" , "Il1r1", "Il23r", "Nrp1", "Lta" , "Il22" )
ILCP_List = c( "Zbtb16", "Tcf7" , "Tox", "Ikzf2" , "Runx3" , "Maf", "Bcl2" , "Cd7" , "Notch2" )
CLP_List = c( "Cd34" , "Bcl11a" , "Myb" , "Flt3" , "Tcf3" , "Nfil3" , "Notch1" , "Rag1", "Ebf1" )

gene_list = list.append(ILC_All_List , ILC1_List, ILC2_List, ILC3_List, ILCP_List, CLP_List)



# Create a named vector mapping genes to their groups
gene_groups = c(rep("ILC_All_List", length(ILC_All_List)),
                rep("ILC1_List", length(ILC1_List)),
                rep("ILC2_List", length(ILC2_List)),
                rep("ILC3_List", length(ILC3_List)),
                rep("ILCP_List", length(ILCP_List)),
                rep("CLP_List", length(CLP_List)))

names(gene_groups) = gene_list

# Define colors for each group
group_colors = c("ILC_All_List" = "#E41A1C",
                 "ILC1_List" = "#377EB8",
                 "ILC2_List" = "#4DAF4A",
                 "ILC3_List" = "#984EA3" ,
                 "ILCP_List" = "#FF7F00",
                 "CLP_List" = "#A65628")

# Function to get the color for each gene
gene_color = function(gene) {
  return(group_colors[gene_groups[gene]])
}

# Create the DotPlot
p11 = DotPlot(Merged_Seurat_Rescaled, features = gene_list, cols = "RdBu") + 
  NoLegend() + 
  ggtitle("ILC signatures") + 
  coord_flip()

# Customize the row names' colors
p11 = p11 + theme(axis.text.y = element_text(color = sapply(gene_list, gene_color)))

print(p11)


if (DO_SAVE_FIGURE == TRUE){
  
png(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/ILC_DotPlot.png"), width = 15, height = 30,  units = "cm", res=600 )
print(p11)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/Integrated_Analysis/ILC_DotPlot.pdf"), width = 15/2.54, height = 30/2.54 )
print(p11)
dev.off()

}









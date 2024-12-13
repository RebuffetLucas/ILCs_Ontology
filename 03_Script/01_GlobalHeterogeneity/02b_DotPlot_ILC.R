# Define gene lists for different ILC groups
ILC_All_List = c("Id2", "Ets1", "Rora", "Ahr", "Cxcr6", "Arg1", "Il7r", "Il2rb")
# Genes associated with all ILC populations.

ILC1_List = c("Tbx21", "Eomes", "Ccl5", "Klrk1", "Ncr1", "Cxcr3", "Xcl1", "Ifng")
# Genes specifically associated with ILC1 population.

ILC2_List = c("Gata3", "Bcl11b", "Il1rl1", "Il17rb", "Icos", "Hey1", "Il9r", "Il4")
# Genes specifically associated with ILC2 population.

ILC3_List = c("Rorc", "Batf3", "Batf", "Tox2", "Il1r1", "Il23r", "Nrp1", "Lta", "Il22")
# Genes specifically associated with ILC3 population.

ILCP_List = c("Zbtb16", "Tcf7", "Tox", "Ikzf2", "Runx3", "Maf", "Bcl2", "Cd7", "Notch2")
# Genes specifically associated with ILCP (ILC progenitor) population.

CLP_List = c("Cd34", "Bcl11a", "Myb", "Flt3", "Tcf3", "Nfil3", "Notch1", "Rag1", "Ebf1")
# Genes specifically associated with CLP (Common Lymphoid Progenitor) population.

# Combine all gene lists into a single list
gene_list = list.append(ILC_All_List, ILC1_List, ILC2_List, ILC3_List, ILCP_List, CLP_List)
# Appends all gene lists to create a unified list for plotting.

# Create a named vector mapping genes to their groups
gene_groups = c(
  rep("ILC_All_List", length(ILC_All_List)),
  rep("ILC1_List", length(ILC1_List)),
  rep("ILC2_List", length(ILC2_List)),
  rep("ILC3_List", length(ILC3_List)),
  rep("ILCP_List", length(ILCP_List)),
  rep("CLP_List", length(CLP_List))
)
# Assigns group labels to each gene based on the list they belong to.

names(gene_groups) = gene_list
# Sets the names of the gene_groups vector as the gene names.

# Define colors for each gene group
group_colors = c(
  "ILC_All_List" = "#E41A1C", # Red for all ILCs
  "ILC1_List" = "#377EB8",    # Blue for ILC1
  "ILC2_List" = "#4DAF4A",    # Green for ILC2
  "ILC3_List" = "#984EA3",    # Purple for ILC3
  "ILCP_List" = "#FF7F00",    # Orange for ILCP
  "CLP_List" = "#A65628"      # Brown for CLP
)

# Function to get the color for each gene
gene_color = function(gene) {
  return(group_colors[gene_groups[gene]])
}
# Defines a function to map a gene to its corresponding color based on the group it belongs to.

# Create the DotPlot
p11 = DotPlot(immune.combined, features = gene_list, cols = "RdBu") + 
  NoLegend() + 
  ggtitle("ILC signatures") + 
  coord_flip()
# Generates a DotPlot for the specified genes, removes the legend, adds a title, and flips the coordinate axes.

# Customize the row names' colors
p11 = p11 + theme(axis.text.y = element_text(color = sapply(gene_list, gene_color)))
# Modifies the y-axis text (gene names) to display their corresponding group colors.

print(p11)
# Displays the generated DotPlot.

# Save the DotPlot if saving is enabled
if (DO_SAVE_FIGURE == TRUE) {
  # Save as PNG
  png(
    file = paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures, "/DotPlot_ILCGenes_Vertical.png"), 
    width = 14, height = 30, units = "cm", res = 600
  )
  print(p11)
  dev.off()
  # Saves the plot as a high-resolution PNG file.
  
  # Save as PDF
  pdf(
    file = paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures, "/DotPlot_ILCGenes_Vertical.pdf"), 
    width = 14/2.54, height = 30/2.54
  )
  print(p11)
  dev.off()
  # Saves the plot as a PDF file with dimensions converted from cm to inches.
}

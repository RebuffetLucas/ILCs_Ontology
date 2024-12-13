# Scoring the klrd1 population with ILC1 vs ILC3 signature and vice versa

# Load the Seurat Object
NK_Seurat = readRDS(file.path(PATH_EXPERIMENT_OUTPUT_RDSFiles, "FL.rds"))

# Set cluster identities and visualize
NK_Seurat = SetIdent(NK_Seurat, value = "seurat_clusters")
DimPlot(NK_Seurat, cols = palette)

PBMC = NK_Seurat #Keep in memory for alternative tries
NK_Seurat = PBMC
DimPlot(NK_Seurat, cols = palette)

# Visualize without removing clusters
DimPlot(NK_Seurat, cols = palette)

# Scoring with ILC1 vs ILC3 signature
# Read signatures
Signatures_ILC1vsILC3 = read.csv(file.path(PATH_EXPERIMENT_REFERENCE, "Signatures_ILCs_Bubulgum_Seydina/ILC3_ILC1.csv"))

# Filter gene lists based on the Seurat object
List_ILC1_vs_ILC3 = list(intersect(Signatures_ILC1vsILC3$ILC1_vs_ILC3, rownames(NK_Seurat)))
List_ILC3_vs_ILC1 = list(intersect(Signatures_ILC1vsILC3$ILC3_vs_ILC1, rownames(NK_Seurat)))

# Add module scores for the signatures
NK_Seurat = AddModuleScore(NK_Seurat, features = List_ILC1_vs_ILC3, pool = NULL, name = "ILC1vsILC3", seed = 19)
NK_Seurat = AddModuleScore(NK_Seurat, features = List_ILC3_vs_ILC1, pool = NULL, name = "ILC3vsILC1", seed = 19)

# Remove the "1" at the end of metadata column names
colnames(NK_Seurat@meta.data) <- sapply(colnames(NK_Seurat@meta.data), function(x) {
  if (grepl("^ILC.*1$", x)) {
    return(gsub("1$", "", x))
  } else {
    return(x)
  }
})

# Feature plots for module scores
p4 = FeaturePlot(NK_Seurat, features = "ILC1vsILC3", pt.size = 0.6) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

p5 = FeaturePlot(NK_Seurat, features = "ILC3vsILC1", pt.size = 0.6) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

# Combine and visualize plots
p4 + p5

# Additional feature plot for "Rorc"
FeaturePlot(NK_Seurat, features = "Rorc", pt.size = 0.6) & 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu"))) & NoAxes()

# Install Libraries
# Install necessary libraries if not already installed
# Uncomment the lines below to install missing packages

# install.packages("sinaplot")
# BiocManager::install("destiny")
# install.packages("rgl")
# install.packages("magick")
# BiocManager::install("genefilter")

# Does not work
# remotes::install_github("Japrin/sscVis")
# devtools::install_github("Japrin/sscVis")

# library(Matrix.utils)


# Loading and preparing the data ###
PBMC = readRDS(file.path(PATH_EXPERIMENT_OUTPUT_RDSFiles, "FL.rds"))
DimPlot(PBMC)

# Remove APOE+ subset
object = subset(PBMC, idents = SUBSETS_TO_REMOVE, invert = TRUE)
object$seurat_clusters = droplevels(object$seurat_clusters)
object = SetIdent(object, value = "seurat_clusters")
DimPlot(object, cols = col_clust_traj)

# Identify variable features
object = FindVariableFeatures(object)
var.features = VariableFeatures(object)
set.seed(SEED)

# Extract data
log.norm.data = object@assays$RNA@data
log.norm.data.var = log.norm.data[var.features, ]
data = t(as.matrix(log.norm.data.var))

# Extract metadata
metadata = object@meta.data


# Diffusion map calculation
matrix = data %>% as.matrix()

# Perform Diffusion Map analysis
dm = DiffusionMap(data = matrix,
                  censor_val = 30,
                  censor_range = c(30, 40),
                  verbose = TRUE)

plot(dm)

# Save Diffusion Map
# saveRDS(dm, file.path(PATH_EXPERIMENT_OUTPUT_DESTINY_OBJECTS, "dm.rds"))

# Calculate diffusion pseudotime
dm_coor = dm@eigenvectors

group = factor(metadata[rownames(dm_coor), "seurat_clusters"], levels = c("C1_ILCP", "C2_Klrd1", "C3_ILC3P", "C4_EILP", "C5_Apoe"))
y = aggregate(dm_coor[, 1], by = list(group), median)
y2 = y[, 2] %>% setNames(y[, 1])
y2 = sort(y2)

if (names(y2)[1] == "C4_EILP") {
  root = which(dm_coor[, 1] == max(dm_coor[, 1]))
} else {
  root = which(dm_coor[, 1] == max(dm_coor[, 1]))
}

dpt_d = DPT(dm, tips = root)

pt = dpt_d[[paste0("DPT", root)]]
pt[pt > quantile(pt, 0.99, na.rm = TRUE)] = NA
pt[pt < quantile(pt, 0.01, na.rm = TRUE)] = NA

df = data.frame("pseudotime" = pt, cluster = group)
rownames(df) = names(dm$DC1)

# For colors
pt_color = dpt_d[[paste0("DPT", root)]]
pt_color[pt_color > quantile(pt_color, 0.99, na.rm = TRUE)] = quantile(pt_color, 0.99, na.rm = TRUE)
pt_color[pt_color < quantile(pt_color, 0.01, na.rm = TRUE)] = quantile(pt_color, 0.01, na.rm = TRUE)

clustering_d = metadata[rownames(dm_coor), "seurat_clusters"] %>%
  as.character() %>% setNames(rownames(dm_coor))

# 2D plots with ggplot
object_d = subset(object, cells = rownames(dm_coor))
object_d[["dm"]] = CreateDimReducObject(dm_coor, key = "DC")
object_d$pt = pt

df2 = data.frame(dm_coor[, 1:3], pseudotime = pt, pseudotime2 = pt_color, seurat_clusters = object_d$seurat_clusters)

# DC 1 & 2 plot
p1 = DimPlot(object_d, reduction = "dm", group.by = "seurat_clusters", pt.size = 1.4,
             cols = col_clust_traj) + theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())
p2 = ggplot(df2, aes(DC1, DC2, col = pseudotime2)) +
  geom_point() + theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_color_gradientn(colours = viridis::plasma(11, alpha = 1, begin = 0, end = 1, direction = 1))

p1 + p2

# Save outputs
png(file = file.path(PATH_EXPERIMENT_OUTPUT_DESTINY_FIGURES, "ClustersandPseudotimes_Normal_NoLTi.png"), width = 50, height = 15, units = "cm", res = 600)
print(p1 + p2)
dev.off()

# Additional pseudotime plots
p5 = ggplot(df, aes(pseudotime, cluster, col = cluster)) +  
  ggbeeswarm::geom_quasirandom(alpha = 0.7, cex = 1, show.legend = FALSE, groupOnX = FALSE) + ylab("")  + 
  ggtitle("Pseudotime in Dataset4") + theme_light() +
  stat_summary(
    aes(group = cluster), fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", color = "black", width = 0.7, lwd = 0.2,
    position = position_dodge(width = 0.75)
  ) +
  scale_color_manual(values = col_clust_traj)

p5

# Save additional plots
png(file = file.path(PATH_EXPERIMENT_OUTPUT_DESTINY_FIGURES, "Plot_PseudoTime_No_LTi.png"), width = 40, height = 15, units = "cm", res = 600)
print(p5)
dev.off()

pdf(file = file.path(PATH_EXPERIMENT_OUTPUT_DESTINY_FIGURES, "Plot_PseudoTime.pdf"), width = 40 / 2.54, height = 15 / 2.54)
print(p5)
dev.off()



#This script run the trajectory analysis for the entire subsets of Dim together

#install.packages("magick")
# This script runs the trajectory analysis for the entire subsets of Dim together

# Start from Seurat object
NK_Seurat = readRDS(file.path(PATH_EXPERIMENT_OUTPUT_RDSFiles, "FL.rds"))
NK_Seurat = SetIdent(NK_Seurat, value = "seurat_clusters")
DimPlot(NK_Seurat, cols = palette) # Plot Seurat clusters

PBMC = NK_Seurat
NK_Seurat = PBMC
DimPlot(NK_Seurat, cols = palette) # Re-plot after assignment

# Remove unwanted clusters
NK_Seurat = subset(PBMC, idents = CLUSTERS_TO_REMOVE, invert = TRUE)
cds = as.cell_data_set(NK_Seurat) # Convert Seurat object to Monocle cell dataset
cds = estimate_size_factors(cds) # Estimate size factors for normalization

# Assign partitions
recreate.partitions = c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) = cds@colData@rownames
recreate.partitions = as.factor(recreate.partitions)
cds@clusters@listData[["UMAP"]][["partitions"]] = recreate.partitions

#Assign cluster informations
list.cluster <- NK_Seurat@active.ident
cds@clusters@listData[["UMAP"]][["clusters"]] <- list.cluster


#Plot before trajectory
partition.before.traj <-plot_cells(cds, color_cells_by = "partition", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
partition.before.traj


cluster.before.traj <-plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F,  group_label_size = 5) + theme(legend.position = "right")
cluster.before.traj


#Learn Trajectory
cds <- learn_graph(cds,  use_partition = TRUE, close_loop = FALSE)


# Visualize trajectory graph
plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = T, trajectory_graph_segment_size = 2)


plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2)

p1 = plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
                label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2)

p1_bis = plot_cells(cds, color_cells_by = "cluster", label_groups_by_cluster = F ,
                label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
                group_label_size = 0, label_principal_points = F, trajectory_graph_segment_size = 2) + scale_color_manual(values = palette)




#Order cells in PseudoTime
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = root_pr_nodes)

plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T,
           label_branch_points = T, label_roots = F, label_leaves = F)

# Plot pseudotime
p2 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = F ,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5,
           graph_label_size = 5, label_principal_points = F, trajectory_graph_segment_size = 2)

p2

p3 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE , label_cell_groups= TRUE,
           label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5, 
           group_label_size =  15, label_principal_points = F, trajectory_graph_segment_size = 1.8, trajectory_graph_color= "blue")

p3

p4 = plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = TRUE , label_cell_groups= TRUE,
                label_branch_points = F, label_roots = F, label_leaves = F, cell_size = 0.5, 
                group_label_size =  15, label_principal_points = F, trajectory_graph_segment_size = 1.8, trajectory_graph_color= "red")

p4


#Save as figures
png(file=paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/TrajOntheUMAP.png"), width = 15, height =12.5, units= "cm" ,  res=600 )
p1
dev.off()

png(file=paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/TrajOntheUMAP_OtherColor.png"), width = 15, height =12.5, units= "cm" ,  res=600 )
p1_bis
dev.off()

png(file=paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/TrajPseudotime.png"), width = 15, height =12.5, units= "cm" ,  res=600 )
p2
dev.off()

png(file=paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/TrajPseudotimeblue.png"), width = 15, height =12.5, units= "cm" ,  res=600 )
p3
dev.off()

png(file=paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/TrajPseudotimered.png"), width = 15, height =12.5, units= "cm" ,  res=600 )
p4
dev.off()





#Filter the most significant genes

modulated_genes <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
genes <- row.names(subset(modulated_genes, q_value == 0 & morans_I > 0.25))
genes

  #Top50
modulated_genes %>% dplyr::arrange(q_value, desc(morans_I)) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(q_value < 0.05 ) %>%
  head(n=NUMBER_GENES_SMALL_HEATMAP) -> top50

row.names(top50)
genes=row.names(top50)


#Order them by pseudotime
pt.pseudotime = as.data.frame(pseudotime(cds)[order(pseudotime(cds))])
pt.pseudotime = t(pt.pseudotime)

#pt.pseudotime= as.matrix(pt.pseudotime)
#pt.pseudotime = rbind(rep(0,  length(pt.pseudotime)), pt.pseudotime)


pt.pseudotime = t(as.matrix(pseudotime(cds)[order(pseudotime(cds))]))


# Plot heatmap for pseudotime
htpseudo <- Heatmap(
  pt.pseudotime,
  name                         = "pseudotime",
  col                          = plasma(11, alpha=1, begin = 0, end=1, direction = 1) ,
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  km = 1,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =unit(0.5, "cm") )


print(htpseudo)


# Process gene expression data
pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]


#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes


# K-means clustering heatmap
#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  km = NUMBER_CLUSTS,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH,
  row_title = NULL)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH)

ht_list = htpseudo %v% htkm
draw(ht_list) #Tcheck that magick is working well if you have some ugly vertical white bars

ht_list2 = htpseudo %v% hthc
draw(ht_list2)


png(file= paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/FL_DynamicHeatmapTop", NUMBER_GENES_SMALL_HEATMAP,"_", NUMBER_CLUSTS ,"GeneClusters.png"), width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list)
dev.off()

png(file= paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/FL_DynamicHeatmapTop", NUMBER_GENES_SMALL_HEATMAP,"_NoGeneClusters.png"), width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list2)
dev.off()


#Top big heatmap
modulated_genes %>% dplyr::arrange(q_value<0.05, desc(morans_I)) %>%
  dplyr::filter(status == "OK") %>%
  dplyr::filter(!grepl("RPS", rownames(.))) %>%
  dplyr::filter(!grepl("RPL", rownames(.))) %>%
  dplyr::filter(!grepl("MT-", rownames(.))) %>%
  dplyr::filter(q_value < 0.05 ) %>%
  head(n=NUMBER_GENES_BIG_HEATMAP) -> top100

row.names(top100)
genes=row.names(top100)


#Order them by pseudotime
pt.pseudotime = as.data.frame(pseudotime(cds)[order(pseudotime(cds))])
pt.pseudotime = t(pt.pseudotime)

#pt.pseudotime= as.matrix(pt.pseudotime)
#pt.pseudotime = rbind(rep(0,  length(pt.pseudotime)), pt.pseudotime)


pt.pseudotime = t(as.matrix(pseudotime(cds)[order(pseudotime(cds))]))


htpseudo <- Heatmap(
  pt.pseudotime,
  name                         = "pseudotime",
  col                          = plasma(11, alpha=1, begin = 0, end=1, direction = 1) ,
  show_row_names               = FALSE,
  show_column_names            = FALSE,
  km = 1,
  row_title_rot                = 0,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =unit(0.5, "cm") )


print(htpseudo)

pt.matrix <- exprs(cds)[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]
pt.matrix <- normalized_counts(cds, norm_method = "log")[match(genes,rownames(rowData(cds))),order(pseudotime(cds))]


#Can also use "normalized_counts" instead of "exprs" to use various normalization methods, for example:
#normalized_counts(cds, norm_method = "log")

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=3)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
rownames(pt.matrix) <- genes



#K means with 6 groups
htkm <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  km = NUMBER_CLUSTS,
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH,
  row_title = NULL)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = POLICE,
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  cluster_rows                 = TRUE,
  show_row_dend = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE,
  heatmap_height =HEIGHT,
  heatmap_width = WIDTH)

ht_list = htpseudo %v% htkm
draw(ht_list)

ht_list2 = htpseudo %v% hthc
draw(ht_list2)




png(file= paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/FL_DynamicHeatmapTop", NUMBER_GENES_BIG_HEATMAP,"_", NUMBER_CLUSTS ,"GeneClusters.png"), width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list)
dev.off()

png(file= paste0(PATH_EXPERIMENT_OUTPUT_MONOCLE_FIGURES , "/FL_DynamicHeatmapTop", NUMBER_GENES_BIG_HEATMAP,"_NoGeneClusters.png"), width = 12, height = 16, units= "cm" ,  res=1200 )
draw(ht_list2)
dev.off()






#Create a Seurat Object from all the data from Liu et al
#Extract annotations
data_anno <- list.files(path=PATH_TO_Liu_Data, pattern = "anno.csv.gz", full.names = TRUE)

myfiles_anno <- lapply(data_anno, function(file) {
  read.csv(file, row.names = 1)  # Assuming the first column should be used as row names
})

# Extract the base names without the path, remove "GSE163587_" and the file extension
names(myfiles_anno) <- sapply(data_anno, function(x) {
  basename <- gsub(pattern = "GSE163587_", replacement = "", x = basename(x))
  sub(pattern = "\\.csv\\.gz$", replacement = "", basename)  # Remove the .csv.gz extension
})

# Optionally print the names to check
list_names_anno=names(myfiles_anno)

print(list_names_anno)

# Add a new column to each dataframe where each entry is the name of the dataframe
myfiles_anno <- lapply(names(myfiles_anno), function(name) {
  df <- myfiles_anno[[name]]
  # Create a new column 'DataFrameName' that stores the modified name (without '_anno')
  new_name <- sub("_anno$", "", name)  # Remove '_anno' from the name
  df <- transform(df, sample.name = new_name)
  return(df)
})

# Remove the '_anno' suffix from each element in the list
list_names_anno <- sub("_anno$", "", list_names_anno)
myfiles_anno <- setNames(myfiles_anno, list_names_anno)



#Extract annotations
data_rawdata<- list.files(path=PATH_TO_Liu_Data, pattern = "rawdata.csv.gz", full.names = TRUE)

myfiles_rawdata<- lapply(data_rawdata, function(file) {
  read.csv(file, row.names = 1)  # Assuming the first column should be used as row names
})

# Extract the base names without the path, remove "GSE163587_" and the file extension
names(myfiles_rawdata) <- sapply(data_rawdata, function(x) {
  basename <- gsub(pattern = "GSE163587_", replacement = "", x = basename(x))
  sub(pattern = "\\.csv\\.gz$", replacement = "", basename)  # Remove the .csv.gz extension
})

# Optionally print the names to check
print(names(myfiles_rawdata))

names(myfiles_rawdata) = sub("_rawdata$", "", list_names_anno)

Seurat_Object_List= list()


for (sample_name in names(myfiles_rawdata)){
  Seurat_Object_To_Add = CreateSeuratObject(myfiles_rawdata[[sample_name]] , meta.data = myfiles_anno[[sample_name]] )
  Seurat_Object_List = append( Seurat_Object_List, Seurat_Object_To_Add)
}

#test_seur= CreateSeuratObject(myfiles_rawdata[["w10_1_rawdata"]] , meta.data = myfiles_anno[["w10_1_anno"]] )

Merged_Seurat_Object = merge(x = Seurat_Object_List[[1]], y = Seurat_Object_List[-1] )


#Load the entire Metadata
All_MetaData = read.xlsx("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/01_Reference/Data/Liu_etal/MetadataTables/Metadata_AllCells_Liu_ESM.xlsx")


#Filter The messy metadata
# Find column names starting with 'integrated_snn'
cols_to_remove <- grep("^integrated_snn", names(All_MetaData), value = TRUE)

# Remove the columns from the dataframe
All_MetaData <- All_MetaData[, !names(All_MetaData) %in% cols_to_remove]
rownames(All_MetaData)= All_MetaData$X1

# Extract row names from Merged_Seurat_Object@meta.data
# Assuming this is a Seurat object, you may need to adjust how you access meta data depending on your Seurat version
row_names_to_match <- rownames(Merged_Seurat_Object@meta.data)

# Reorder All_MetaData to match the row order of Merged_Seurat_Object@meta.data
All_MetaData_reordered <- All_MetaData[match(row_names_to_match, rownames(All_MetaData)), ]

#Reinject the reordered Metadata
Merged_Seurat_Object@meta.data = All_MetaData_reordered

#Normalize the data
Merged_Seurat_Object = NormalizeData(Merged_Seurat_Object)

#Insert the data embeddings
Merged_Seurat_Object[['umap']] <- CreateDimReducObject(embeddings = as.matrix(All_MetaData_reordered[,c('UMAP_1','UMAP_2')]), key = 'UMAP_', assay = 'RNA') 
Merged_Seurat_Object[['tsne']] <- CreateDimReducObject(embeddings = as.matrix(All_MetaData_reordered[,c('tSNE_1','tSNE_2')]), key = 'TSNE_', assay = 'RNA') 


#Data Viz
DimPlot(Merged_Seurat_Object, reduction = "tsne", group.by = "cells_used")

#Save the object
saveRDS(Merged_Seurat_Object , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Seurat_Object/Liu_Merged_Integrated_Seurat_Object_AllCells.rds")




#Now take care of Lymphoid cells only (annotated as "used in their all data)
Merged_Seurat_Object2 = subset(Merged_Seurat_Object, subset= cells_used=="used")
MetaData_Lymphoid = read.xlsx("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/01_Reference/Data/Liu_etal/MetadataTables/Metadata_Lymphoidcells_Liu.xlsx")

#Filter The messy metadata
# Find column names starting with 'integrated_snn'
cols_to_remove2 <- grep("^integrated_snn", names(MetaData_Lymphoid), value = TRUE)

# Remove the columns from the dataframe
MetaData_Lymphoid <- MetaData_Lymphoid[, !names(MetaData_Lymphoid) %in% cols_to_remove2]
rownames(MetaData_Lymphoid)= MetaData_Lymphoid$X1

# Extract row names from Merged_Seurat_Object2@meta.data
# Assuming this is a Seurat object, you may need to adjust how you access meta data depending on your Seurat version
row_names_to_match <- rownames(Merged_Seurat_Object2@meta.data)

# Reorder MetaData_Lymphoid to match the row order of Merged_Seurat_Object2@meta.data
MetaData_Lymphoid_reordered <- MetaData_Lymphoid[match(row_names_to_match, rownames(MetaData_Lymphoid)), ]

#Reinject the reordered Metadata
Merged_Seurat_Object2@meta.data = MetaData_Lymphoid_reordered

#Normalize the data
Merged_Seurat_Object2 = NormalizeData(Merged_Seurat_Object2)

#Insert the data embeddings
Merged_Seurat_Object2[['umap']] <- CreateDimReducObject(embeddings = as.matrix(MetaData_Lymphoid_reordered[,c('UMAP_1','UMAP_2')]), key = 'UMAP_', assay = 'RNA') 
Merged_Seurat_Object2[['tsne']] <- CreateDimReducObject(embeddings = as.matrix(MetaData_Lymphoid_reordered[,c('tSNE_1','tSNE_2')]), key = 'TSNE_', assay = 'RNA') 


Merged_Seurat_Object2 = SetIdent(Merged_Seurat_Object2 ,  value = "Cluster")

#saveit
saveRDS(Merged_Seurat_Object2 , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Seurat_Object/Liu_Merged_Integrated_Seurat_Object_OnlyUsedCells.rds")








#Now take care of ILCP and Pre_ILC1, Pre_ILC2 and Pre_ILC3 cells only (corresponding to their Fig 7)
Metadata_ILCProge =read.xlsx("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/01_Reference/Data/Liu_etal/MetadataTables/Metadata_Pre_ILCandILCP_Liu.xlsx")


Merged_Seurat_Object3 = subset(Merged_Seurat_Object2, cells= Metadata_ILCProge$X1)

#Filter The messy metadata
# Find column names starting with 'integrated_snn'
cols_to_remove3 <- grep("^integrated_snn", names(Metadata_ILCProge), value = TRUE)

# Remove the columns from the dataframe
Metadata_ILCProge <- Metadata_ILCProge[, !names(Metadata_ILCProge) %in% cols_to_remove3]
rownames(Metadata_ILCProge)= Metadata_ILCProge$X1

# Extract row names from Merged_Seurat_Object2@meta.data
# Assuming this is a Seurat object, you may need to adjust how you access meta data depending on your Seurat version
row_names_to_match <- rownames(Merged_Seurat_Object3@meta.data)

# Reorder Metadata_ILCProge to match the row order of Merged_Seurat_Object3@meta.data
Metadata_ILCProge_reordered <- Metadata_ILCProge[match(row_names_to_match, rownames(Metadata_ILCProge)), ]

#Reinject the reordered Metadata
Merged_Seurat_Object3@meta.data = Metadata_ILCProge_reordered

#Normalize the data
Merged_Seurat_Object3 = NormalizeData(Merged_Seurat_Object3)

#Insert the data embeddings
Merged_Seurat_Object3[['umap']] <- CreateDimReducObject(embeddings = as.matrix(Metadata_ILCProge_reordered[,c('UMAP_1','UMAP_2')]), key = 'UMAP_', assay = 'RNA') 
Merged_Seurat_Object3[['tsne']] <- CreateDimReducObject(embeddings = as.matrix(Metadata_ILCProge_reordered[,c('tSNE_1','tSNE_2')]), key = 'TSNE_', assay = 'RNA') 


Merged_Seurat_Object3 = SetIdent(Merged_Seurat_Object3 ,  value = "cluster")

#saveit
saveRDS(Merged_Seurat_Object3 , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Seurat_Object/Liu_Merged_Integrated_Seurat_Object_Only_Progenitors.rds")



#Not used#
#Creating the Lymphoid alone (corresponding to their 


if (PREPARE_LIU_DATA ==TRUE){
  
  #Normalize and preprocess
  Merged_Seurat_Object = NormalizeData(Merged_Seurat_Object)
  all.genes = rownames(Merged_Seurat_Object)
  Merged_Seurat_Object = ScaleData(Merged_Seurat_Object, features = all.genes, do.scale = DATA_SCALE, do.center = DATA_CENTER)
  Merged_Seurat_Object= FindVariableFeatures(Merged_Seurat_Object)
  Merged_Seurat_Object = RunPCA(Merged_Seurat_Object)
  
  #BATCH CORRECTION
  Merged_Seurat_Object = RunHarmony(Merged_Seurat_Object, HARMONIZATION_CRITERIA)
  
  #Clust and Neighbors
  Merged_Seurat_Object <- FindNeighbors(Merged_Seurat_Object, reduction= "harmony", dims = FINDCLUSTERS_DIMS)
  Merged_Seurat_Object <- FindClusters(Merged_Seurat_Object, resolution = FINDCLUSTERS_RESOLUTION)
  
  #Run UMAP
  Merged_Seurat_Object <- RunUMAP(Merged_Seurat_Object, reduction = "harmony",dims = FINDCLUSTERS_DIMS )
  
}


#Have a quick Look / Diagnostic


All_Markers = FindAllMarkers(Merged_Seurat_Object , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10

p3 = DotPlot(Merged_Seurat_Object, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) &ggtitle("Spontaneous signatures")



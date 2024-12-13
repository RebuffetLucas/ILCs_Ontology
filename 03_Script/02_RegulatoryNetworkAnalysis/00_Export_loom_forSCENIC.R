### This script creates a loom file for BM and FL

#Import libraries
library(scater)
library(SeuratDisk)
library(loomR)
library(SCopeLoomR)

#Create a loom for BM:
BM_Object =readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/01_GlobalHeterogeneity/Seurat_Objects/BM.rds")

DimPlot(BM_Object, reduction = "tsne")

build_loom(file.name = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM.loom",
           dgem= BM_Object@assays$RNA@counts,
           title= "BM" ,
          default.embedding= BM_Object@reductions$tsne@cell.embeddings ,
           default.embedding.name="tsne"
           )

loom <- open_loom("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM.loom", mode = "r+")
close_loom(loom)

#Save the metadata
BM_Info= BM_Object@meta.data

saveRDS(BM_Info , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/BM/BM_Info.rds")



#Create a loom for FL:
FL_Object =readRDS("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/01_GlobalHeterogeneity/Seurat_Objects/FL.rds")

DimPlot(FL_Object, reduction = "umap")

build_loom(file.name = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/FL/FL.loom",
           dgem= FL_Object@assays$RNA@counts,
           title= "FL" ,
           default.embedding= FL_Object@reductions$umap@cell.embeddings ,
           default.embedding.name="umap"
)

loom <- open_loom("/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/FL/FL.loom", mode = "r+")
close_loom(loom)

#Save the metadata
FL_Info= FL_Object@meta.data
saveRDS(FL_Info , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/02_RegulatoryNetworkAnalysis/data/FL/FL_Info.rds")






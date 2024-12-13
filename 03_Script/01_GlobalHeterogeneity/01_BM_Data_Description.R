#This scripts gives a quick overview of the BM data
#Annotates the populations and save a new object
#Load the data sent by our colaborators
load(file.path(PATH_EXPERIMENT_REFERENCE_Data,SAMPLE_BM))
p0 = DimPlot(immune.combined, label = TRUE) #Have a quick look
p0


#Diagnostic with their granularity

All_Markers = FindAllMarkers(immune.combined , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)
#SaveRDS(All_Markers, file= "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/01_GlobalHeterogeneity/DEG/All_Markers.rds")

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10


p1 = DotPlot(immune.combined, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90))

p1

#Rename and save
levels(immune.combined$seurat_clusters) =  c("ILCP", "ILC2P", "EILP")
immune.combined = SetIdent(immune.combined, value = "seurat_clusters")
SAVE_RDS_FILE = paste0(PATH_EXPERIMENT_OUTPUT_RDSFiles, "/BM.rds")
saveRDS(immune.combined, file =  SAVE_RDS_FILE)


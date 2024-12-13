# Extract the signatures
#For BM
SAVE_RDS_FILE = paste0(PATH_EXPERIMENT_OUTPUT_RDSFiles, "/BM.rds")
immune.combined = readRDS(SAVE_RDS_FILE)
All_Markers = FindAllMarkers(immune.combined , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

saveRDS(All_Markers, file= "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/01_GlobalHeterogeneity/DEG/All_Markers_BM.rds")

#For FL
SAVE_RDS_FILE = paste0(PATH_EXPERIMENT_OUTPUT_RDSFiles, "/FL.rds")
immune.combined = readRDS(SAVE_RDS_FILE)
All_Markers = FindAllMarkers(immune.combined , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)
saveRDS(All_Markers, file= "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/01_GlobalHeterogeneity/DEG/All_Markers_FL.rds")

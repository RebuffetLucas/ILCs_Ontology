###############################################################################
# This file defines SAMPLE parameters as global variables that will be loaded
# before analysis starts. 

# Declare all the samples names
SAMPLE_BM = "BM.RData"
SAMPLE_FL= "FL.RData"


#Human samples


SAMPLE_FL_Science_YS_Atlas = "/01_Reference/Data/YS_paper/Seurat_Objects/Foetal_liver_data"
SAMPLE_YS_Science_YS_Atlas = "/01_Reference/Data/YS_paper/Seurat_Objects/YS_data"

SAMPLE_Liu_Data_All = "/05_Output/04_Mapping_Human_Datasets/Seurat_Object/Liu_Merged_Integrated_Seurat_Object_AllCells.rds"
SAMPLE_Liu_Data_Used = "/05_Output/04_Mapping_Human_Datasets/Seurat_Object/Liu_Merged_Integrated_Seurat_Object_OnlyUsedCells.rds"

SAMPLE_Liu_Data_Used_ILC_Proge ="/05_Output/04_Mapping_Human_Datasets/Seurat_Object/Liu_Merged_Integrated_Seurat_Object_Only_Progenitors.rds"


#Mice samples
SAMPLE_ROMAGNANI = "/01_Reference/Data/Chiara_R_Data/RORcGFP_CC_regressed.rds"

SAMPLE_Bai = "/01_Reference/Data/Bai_Etal/ILCIP_LSM_CCA.Rds"

SAMPLE_CHAP2_FL = "/05_Output/01_GlobalHeterogeneity/Seurat_Objects/FL.rds"
SAMPLE_CHAP2_BM = "/05_Output/01_GlobalHeterogeneity/Seurat_Objects/BM.rds"


# Group the samples in a vector (usefull for llops over the samples)
SAMPLE_SET = c( SAMPLE_BM, SAMPLE_FL)


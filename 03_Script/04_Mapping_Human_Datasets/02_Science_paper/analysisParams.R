###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis

#Directories etc...

ANALYSIS_STEP_NAME = "04_Mapping_Human_Datasets"
SUB_ANALYSIS_STEP_NAME = "02_Science_paper"



#Colors

#palette<-c('NK1C'='#F8766D','NKint'='#8494FF',
#           'NK1A'='#0CB702',
#           'NK1B'='#00BFC4','NK2'='#ED68ED',
#           'NK3'='#ABA300')

# Seed for pseudo-random numbers
SEED = 42;


# Scaling parameters (see Seurat::ScaleData())
DATA_CENTER       = TRUE
DATA_SCALE        = TRUE
DATA_VARS_REGRESS = NULL  # c("nCount_RNA") for UMIs (NULL to ignore)

#### Analysis parameters

# Maximum number of variable features to keep
VARIABLE_FEATURES_MAXNB   = 2000;  # For analysis (PCA)
VARIABLE_FEATURES_SHOWTOP = 200;   # For table in report


# Nearest-neighbor graph construction
FINDNEIGHBORS_K = 30

# Cluster identification parameters
FINDCLUSTERS_RESOLUTION     = 0.5;
FINDCLUSTERS_USE_PCA_NBDIMS = 30;  # Number of dimensions to use from PCA results
FINDCLUSTERS_ALGORITHM      = 1;   # 1 = Louvain; 2 = Louvain with multilevel refinement; 3 = SLM; 4 = Leiden
FINDCLUSTERS_DIMS = 1:20

# Parameters for identification of marker annotations for clusters (see Seurat::FindAllMarkers())
FINDMARKERS_METHOD    = "wilcox"  # Method used to identify markers
FINDMARKERS_ONLYPOS   = TRUE     # Only consider overexpressed annotations for markers ? (if FALSE downregulated genes can also be markers)

FINDMARKERS_MINPCT    = 0.2      # Only test genes that are detected in a minimum fraction of cells in either of the two populations. Speed up the function by not testing genes that are very infrequently expressed. Default is '0.1'.
FINDMARKERS_LOGFC_THR = 0.25     # Limit testing to genes which show, on average, at least X-fold difference (log-scale) between the two groups of cells. Default is '0.25'. Increasing logfc.threshold speeds up the function, but can miss weaker signals.
FINDMARKERS_PVAL_THR  = 0.001    # PValue threshold for identification of significative markers
FINDMARKERS_SHOWTOP   = 10      # Number of marker genes to show in report and tables (NULL for all)



# Scoring comparizon chapter 1
NUMBER_TOP_SCORING = 40
COR_pVALUE_THR = 5e-2
LOGFC_THR = 0.5

#SUBSETS
NOT_PROGENITORS = c("DC", "B_lymphoid", "Endothelium", "Erythroid", "Fibroblast", "Granulocyte_precursor", "Hepatocyte", "Kupffer_cell", "Mast_cell", "Monocyte", "MK", "NK")

DO_REMOVE_NK = TRUE
#Save_Figures

DO_SAVE_FIGURE = TRUE 
SAVE_FIG_PATH = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Figures"
SAVE_FIG_PATH = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Figures/Mapping_Science_Paper/Working_On_Subset_FL_YS"
#Resource
MOUSE_TO_HUM_FL_PATH = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Mouse_To_Human_Conversion/Conversion_MouseFLtoHuman.xlsx"


#SubClustering_Study
SUBSET_TO_SUBCLUSTER = c("ILC_precursor")
PREPARE_SUBSET_DATA = TRUE

RESOLUTION_SUBCLUST = 0.3

FINDMARKERS_SHOWTOP_SUBCLUST = 20


#Harmonization of the extracted cluster

HARMONIZATION_CRITERIA = "fetal.ids"





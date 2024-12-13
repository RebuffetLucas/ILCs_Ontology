###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis


#Colors
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
FINDMARKERS_SHOWTOP   = 20      # Number of marker genes to show in report and tables (NULL for all)



# Scoring comparizon chapter 1
NUMBER_TOP_SCORING = 40
COR_pVALUE_THR = 5e-2
LOGFC_THR = 0.5


#Source of data
FILE_SIGNATURES_PATH = file.path( "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/01_GlobalHeterogeneity/DEG/All_Markers_FL.rds")

FILTER_pvalue_THR = 0.05


ARRANGE_BY = "p_val_adj"

LENGTH_EXTRACT_BEST = 200 #Number of best mouse genes extracted in each group


#Manual inputs for specific signatures / genes that do not convert well
#df_Mouse_Human_Manual= data.frame(
#  gene = c("Dnaja1" ,  "H19", "Ldha", "Dbi", "Il4", "Klrd1", "Ifng", "Cd52", "Cd7", "Klrb1b", "Cd200r4", "Klrc1", "Tmsb10", "Ly6a", "Tmsb15a", "Prr29"),
#  human.gene.name = c("DNAJA1" , "H19", "LDHA", "DBI", "IL4" , "KLRD1", "IFNG" , "CD52", "CD7", "KLRB1", "CD200R1", "KLRC1", "TMSB10", "LY6H" , "TMSB15A", "PRR29" )
#)



PATH_SAVE_MOUSE_TO_HUMAN_SIGNATURES_ANDpvalue = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Mouse_To_Human_Conversion/Conversion_MouseFLtoHuman.xlsx"

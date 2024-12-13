###############################################################################
# This file defines ANALYSIS parameters as global variables that will be loaded
# before analysis starts. It should define common parameters used by the current
# analysis

#Directories etc...

ANALYSIS_STEP_NAME = "05_Mapping_Mouse_Datasets"
SUB_ANALYSIS_STEP_NAME = "03_Multi_Organ"


SIGNATURES_ROMAGNANI_ILC3 = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/01_Reference/Data/Chiara_R_Data/ExFig2a_Top30.csv"

SIGNATURES_ROMAGNANI_All = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/01_Reference/Data/Chiara_R_Data/ExFig1c_Top50.csv"

#Colors

#palette<-c('NK1C'='#F8766D','NKint'='#8494FF',
#           'NK1A'='#0CB702',
#           'NK1B'='#00BFC4','NK2'='#ED68ED',
#           'NK3'='#ABA300')


#When using public data
ANNOTATION_OF_INTEREST = "annotation"
REDUCTION_OF_INTEREST= "umap"

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
FINDCLUSTERS_RESOLUTION     = 0.4;
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
LOGFC_THR = 0.25

#Save_Figures

DO_SAVE_FIGURE = TRUE 

SAVE_FIG_PATH = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/05_Mapping_Mouse_Datasets/03_Multi_Organ/Figures"

#Resource
MOUSE_TO_HUM_FL_PATH = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/04_Mapping_Human_Datasets/Mouse_To_Human_Conversion/Conversion_MouseFLtoHuman.xlsx"

FILE_SIGNATURES_PATH = "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/01_GlobalHeterogeneity/DEG/All_Markers_FL.rds"



#SubClustering_Study

FINDMARKERS_SHOWTOP_SUBCLUST = 20


#Harmonization of the extracted cluster

HARMONIZATION_CRITERIA = "sample.name"

#PCA Analysis:
USE_MULTIBATCHNORM = TRUE


NUMBER_PCA_COMPONENT_psoleil = 10
NUMBER_VAR_BARPLOTS_PCA_COMPO = 20

PCA_PLOT_POINT_SIZE = 0.6
PCA_PLOT_TEXT_SIZE = 6
N_FEATURES_PCA = 500


#Harmony Integration
USE_MULTIBATCHNORM_BEFORE_HARMONY =  TRUE


VARIABLE_HARMONY = c("Organ")

N_INTEGRATION_FEATURES_HARMONY = 2000

FINDCLUSTERS_RESOLUTION_AFTER_HARMONY = 0.5


#After Harmony:
RENAME_CLUSTERS = c("ILC3P", "ILC3P" , "ILC3P", "ILC2P", "ILC2P", "Klrd1", "Klrd1", "EILP", "ILCP", "ILC2P", "Apoe", "ILC2P", "Klrd1_ILC3P", "EILP", "ILC3P", "ILC2P", "ILCP")


palette_3<-c('ILCP'='#F8766D',
             'Klrd1' = '#A3A500',
             "ILC3P" = '#00BF7D',
             'EILP' =  '#00B0F6' ,
             'Apoe' =  '#E76BF3',
             "Klrd1_ILC3P" ='#8494FF' ,
             "ILC2P"= '#FF61CC' )





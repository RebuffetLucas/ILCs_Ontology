# ##################################################
# Global declarations and libraries for the analysis
# ##################################################


######## R Libraries
library(Matrix)
library( digest)        # digest (hashing)
library( DT)            # datatable
library( forcats)       # fct_inorder (factors)
library( fs)            # path_sanitize
library( future)        # plan (multicore)
library( ggplot2)
#library( ggpubr)
library( ggrepel)
#library( grid)
#library( gridExtra)
library( htmltools)     # browsable
library( htmlwidgets)   # JS (for datatable)
#library( iheatmapr)     # iheatmap
library( knitr)
 
#library( pander)        # pander
library( patchwork)     # +/ (ggplot layout)
library( pheatmap)      # pheatmap
library( plotly)
#library( plyr)

library( rmarkdown)
library( scales)        # hue_pal
#library(magick)


#library(genefilter)
library(SeuratWrappers)
#library(scatterplot3d) 

#library(UCell)

#library(sinaplot)


# Single-cell technology
 
#library(rgl)


# Functional Enrichment analysis
library( biomaRt)
#library( clusterProfiler)
#library( org.Mm.eg.db)
#library( org.Hs.eg.db)
#library( rrvgo)

#library(gghalves)
#library(ggdist)
library(ComplexHeatmap)

library(dplyr)
library(plyr)
library(readxl)
library(RColorBrewer)


#For batch correction with Seurat
library(BiocManager)
#library(multtest)
#library(metap)

#library(destiny)


#Harmony
#library(harmony)
library(batchelor)
library(limma)
#library(unix)
#library(DESeq2)
#library(rlist)
library(reshape2)

library(sctransform)
#library(SeuratDisk)
library(viridis)
library(viridisLite)
library(SeuratWrappers)

# Perform the transcriptionnal trajectories analysis (Monocle 3) on NK Cells
library(monocle3)

# Version of seurat Wrapper that actually works
# remotes::install_github('satijalab/seurat-wrappers', ref = "b8feea013e7e19a46e935684b510399ffe0b6740" #

# Perform the transcriptionnal trajectories analysis (Monocle 3) on NK Cells
library(monocle3)
library(Seurat)
library(SeuratData)
 

# Version of seurat Wrapper that actually works
# remotes::install_github('satijalab/seurat-wrappers', ref = "b8feea013e7e19a46e935684b510399ffe0b6740" #


library(magick)
library(magrittr)
#library( magick)
library(cerebroApp)






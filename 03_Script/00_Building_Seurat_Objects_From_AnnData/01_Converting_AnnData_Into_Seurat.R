# Load the data taken from the online portal
# Handle FL (Fetal Liver)
Convert(
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_FL_Science_YS_Atlas, ".h5ad")), 
  "h5seurat", 
  overwrite = TRUE
) 
# Converts the .h5ad file (an AnnData format) to Seurat's .h5seurat format, allowing compatibility with Seurat.

FL_YSpaper = LoadH5Seurat(
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_FL_Science_YS_Atlas, ".h5seurat")), 
  meta.data = FALSE, 
  misc = FALSE
)
# Loads the converted Seurat object, excluding metadata and miscellaneous information.

table_FL = read.csv(
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_FL_Science_YS_Atlas, ".csv"))
)
# Reads a CSV file containing metadata for the Fetal Liver dataset.

rownames(table_FL) = table_FL$X
# Sets the row names of the metadata table to the values in the column `X`.

table_FL$sex = NULL
table_FL$orig.dataset = NULL
# Removes the `sex` and `orig.dataset` columns from the metadata table as they are not needed.

FL_YSpaper@meta.data = table_FL 
# Adds the processed metadata table to the Seurat object's metadata slot.

FL_YSpaper = NormalizeData(FL_YSpaper)
# Normalizes the expression data in the Seurat object (e.g., log-transformation).

FL_YSpaper = SetIdent(FL_YSpaper, value = "cell.labels")
# Sets the default identity class of the Seurat object to the "cell.labels" column from the metadata.

saveRDS(
  FL_YSpaper, 
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_FL_Science_YS_Atlas, ".rds"))
)
# Saves the processed Seurat object as an .rds file for later use.

# Handle YS (Yolk Sac)
Convert(
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_YS_Science_YS_Atlas, ".h5ad")), 
  "h5seurat", 
  overwrite = TRUE
) 
# Converts the .h5ad file for the Yolk Sac dataset into Seurat's .h5seurat format.

YS_YSpaper = LoadH5Seurat(
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_YS_Science_YS_Atlas, ".h5seurat")), 
  meta.data = FALSE, 
  misc = FALSE
)
# Loads the converted Seurat object for the Yolk Sac dataset, excluding metadata and miscellaneous information.

table_YS = read.csv(
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_YS_Science_YS_Atlas, ".csv"))
)
# Reads the metadata CSV file for the Yolk Sac dataset.

rownames(table_YS) = table_YS$X
# Sets the row names of the metadata table to the values in the column `X`.

YS_YSpaper@meta.data = table_YS 
# Adds the processed metadata table to the Seurat object's metadata slot.

YS_YSpaper = NormalizeData(YS_YSpaper)
# Normalizes the expression data in the Seurat object.

YS_YSpaper = SetIdent(YS_YSpaper, value = "LVL3")
# Sets the default identity class of the Seurat object to the "LVL3" column from the metadata.

saveRDS(
  YS_YSpaper, 
  file.path(PATH_EXPERIMENT_REFERENCE_Data, paste0(SAMPLE_YS_Science_YS_Atlas, ".rds"))
)
# Saves the processed Seurat object as an .rds file for later use.

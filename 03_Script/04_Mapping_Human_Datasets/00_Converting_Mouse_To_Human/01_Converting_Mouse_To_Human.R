#Apply scores of Dim, Bright and Adapt to the initial UMAP


#Import signatures
MONITORED_Markers = readRDS(FILE_SIGNATURES_PATH)
  

###Markers for scoring of subsets
# getting data from sheets
List_Monitored= list()
List_Df_Monitored= list()

for( clust in(levels(MONITORED_Markers$cluster))){
  MONITORED_Markers %>%
  dplyr::filter(avg_log2FC > 0) %>%
  dplyr::filter(cluster == clust) %>%
  dplyr::arrange(ARRANGE_BY) %>%
  filter(p_val_adj<FILTER_pvalue_THR)  -> top40_To_Add
  
  #print(clust)
  List_Monitored[[clust]] = top40_To_Add$gene[1:LENGTH_EXTRACT_BEST]
  List_Df_Monitored[[clust]] = top40_To_Add[1:LENGTH_EXTRACT_BEST,]
  #print(top40_To_Add)
  
}




#Convert into human genes
#Sabrina's function to convert to human genes
convertHumanGeneList <- function(x){
  
  require("biomaRt")
  # ensembl <- useMart(biomart = "ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  mouse_attributes <- c("external_gene_name","external_gene_source","hsapiens_homolog_associated_gene_name","hsapiens_homolog_perc_id",
                        "hsapiens_homolog_goc_score", "hsapiens_homolog_orthology_confidence")
  genesV2<-getLDS(attributes = "hgnc_symbol",
                  filters = "hgnc_symbol", values = x, mart = ensembl.human,
                  attributesL = mouse_attributes, martL = ensembl.mouse,uniqueRows = TRUE)
  mousex <- genesV2 %>%
    filter(Gene.name != "") %>%
    group_by(HGNC.symbol) %>%
    slice_max(Human.Gene.order.conservation.score) %>%
    filter(Human.orthology.confidence..0.low..1.high. == 1) %>%
    dplyr::select(HGNC.symbol,Gene.name) %>%
    dplyr::rename(human.gene.name = HGNC.symbol, mouse.gene.name = Gene.name)
  
  return(mousex)
}


convertMouseGeneList <- function(x){
  
  require("biomaRt")
  # ensembl <- useMart(biomart = "ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  human_attributes <- c("external_gene_name","external_gene_source","mmusculus_homolog_associated_gene_name","mmusculus_homolog_perc_id",
                        "mmusculus_homolog_goc_score", "mmusculus_homolog_orthology_confidence")
  genesV2<-getLDS(attributes = "mgi_symbol",
                  filters = "mgi_symbol", values = x, mart = ensembl.mouse,
                  attributesL = human_attributes, martL = ensembl.human,uniqueRows = TRUE)
  
  
  humanx <- genesV2 %>%
    filter(Gene.name != "") %>%     #Remove the filter to prevent NA
    group_by(MGI.symbol) %>%
    slice_max(Mouse.Gene.order.conservation.score) %>%
    #filter(Mouse.orthology.confidence..0.low..1.high. == 1) %>%  #Remove the filter to prevent NA
    dplyr::select(MGI.symbol,Gene.name) %>%
    distinct(MGI.symbol, .keep_all = TRUE) %>% # ensuring that each mouse gene appears only once.
    dplyr::rename(mouse.gene.name = MGI.symbol, human.gene.name = Gene.name)
  
  return(humanx)
}



  #Convert into human

List_Monitored #Have a look at the best genes used for conversion

Human_Monitored = lapply(List_Monitored , convertMouseGeneList) #Convert Mouse to Human, we removed the filter for reliability here


#Convert Mouse to Human


List_name_DF= names(List_Df_Monitored) #Names of the clusters


#Put the human names in the df
for (name_DF in List_name_DF ){
  
  List_Df_Monitored[[name_DF]] <- List_Df_Monitored[[name_DF]] %>%
    left_join(Human_Monitored[[name_DF]] , by = c("gene" = "mouse.gene.name"))
    
}


#Have a look at the human gene NA:
# Assuming your list of data frames is called `list_of_dfs`
filtered_list_of_dfs <- lapply(List_Df_Monitored, function(df) {
  df[is.na(df$human.gene.name), ]
})

# If you want to check the results for a specific dataframe in the list

print(filtered_list_of_dfs)
print(filtered_list_of_dfs[[1]])


#Now keep only the dataframe without NA

List_of_dfs_clean <- lapply(List_Df_Monitored, function(df) {
  df %>%
    filter(!is.na(human.gene.name))
})





#Save a xlsxone df per page
# Create a new Excel workbook
wb <- createWorkbook()

# Loop through the list and add each dataframe as a new sheet
for(i in seq_along(List_of_dfs_clean)) {
  addWorksheet(wb, sheetName = names(List_of_dfs_clean)[i])  # Or use names(List_of_dfs_clean)[i] if your list has named elements
  writeData(wb, sheet = i, List_of_dfs_clean[[i]])
}


# Save the workbook
saveWorkbook(wb, PATH_SAVE_MOUSE_TO_HUMAN_SIGNATURES_ANDpvalue, overwrite = TRUE)






#More results but less reliable (this is not the function used here) Kept as a resource, just in case

convertMouseGeneList_NoFilter <- function(x){
  
  require("biomaRt")
  # ensembl <- useMart(biomart = "ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  ensembl.human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host="https://dec2021.archive.ensembl.org")
  human_attributes <- c("external_gene_name","external_gene_source","mmusculus_homolog_associated_gene_name","mmusculus_homolog_perc_id",
                        "mmusculus_homolog_goc_score", "mmusculus_homolog_orthology_confidence")
  genesV2<-getLDS(attributes = "mgi_symbol",
                  filters = "mgi_symbol", values = x, mart = ensembl.mouse,
                  attributesL = human_attributes, martL = ensembl.human,uniqueRows = TRUE)
  
  
  humanx <- unique(genesV2[, 2])
  
  
  return(humanx)
}








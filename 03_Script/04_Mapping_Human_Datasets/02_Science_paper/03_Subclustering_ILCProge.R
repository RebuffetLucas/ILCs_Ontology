## @knitr SubClustand_Scoring_With_Mouse_FL_Signatures

PBMC3 = subset(PBMC2, idents= SUBSET_TO_SUBCLUSTER)

cat(" \n \n")
cat("## Subclustering Human FL " , SUBSET_TO_SUBCLUSTER , " and scoring with with mouse FL signatures")
cat(" \n \n")



#Prepare the data

if (PREPARE_SUBSET_DATA == TRUE){
  all.genes = rownames(PBMC3)
PBMC3 = ScaleData(PBMC3, features = all.genes)
PBMC3= FindVariableFeatures(PBMC3)
PBMC3 = RunPCA(PBMC3)
}

PBMC3 <- FindNeighbors(PBMC3, reduction= "pca", dims = 1:20)
PBMC3 <- FindClusters(PBMC3, resolution = RESOLUTION_SUBCLUST)

PBMC3 <- RunUMAP(PBMC3, dims = 1:20)


All_Markers = FindAllMarkers(PBMC3 , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)

All_Markers %>%
  group_by(cluster) %>%
  filter(!grepl("RPL|RPS|MT-", gene)) %>% #We use RPS and RPL for scoring but nor for plotind DotPlots
  slice_max(n = FINDMARKERS_SHOWTOP_SUBCLUST, order_by = avg_log2FC) -> top10

p1 = DotPlot(PBMC3, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) & ggtitle( paste0( SUBSET_TO_SUBCLUSTER , " Subclustering spontaneous signatures No RBS / RPS"))

cat(" \n \n")
print(p1)
cat(" \n \n")

All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP_SUBCLUST, order_by = avg_log2FC) -> top10

p2 = DotPlot(PBMC3, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90)) & ggtitle( paste0( SUBSET_TO_SUBCLUSTER , " Subclustering spontaneous signatures"))


cat(" \n \n")
print(p2)
cat(" \n \n")




p_subClust = DimPlot(PBMC3, pt.size = 1.2) & ggtitle( paste0( SUBSET_TO_SUBCLUSTER , "Subclustering"))

cat(" \n \n")
print(p_subClust)
cat(" \n \n")




#Scoring the subclusters

for (names_Sign in names(human_gene_names_best)){
  PBMC3@meta.data[[paste0(names_Sign, "1")]] = NULL #Remove Previous scoring
  PBMC3 = AddModuleScore( PBMC3 , features = human_gene_names_best[names_Sign] , pool= NULL ,name= names_Sign, seed=19) #Add the new scoring
}





cat(" \n \n")
cat("### Applying the different mouse FL scores {.tabset .tabset-fade} \n\n")
cat(" \n \n")

palette_3 = hue_pal()(length(levels(PBMC3@active.ident)))
names(palette_3) = levels(PBMC3@active.ident)




for (names_signatures in paste0(names(human_gene_names_best)  , "1")){
  
  cat(" \n \n")
  cat("#### ",  names_signatures, "\n")
  cat(" \n \n")
  
  p1 = FeaturePlot(PBMC3, features = names_signatures)  &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  p2 = VlnPlot(PBMC3, features = names_signatures ,sort = "decreasing", pt.size = 0, cols = palette_3) & stat_summary(fun.data=data_summary,color="black")
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
}



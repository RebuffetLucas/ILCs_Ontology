## @knitr SubClustand_Scoring_With_Mouse_FL_Signatures

cat(" \n \n")
cat("## Have a look at features of interest {.tabset .tabset-fade} \n\n")
cat(" \n \n")


PBMC3 = PBMC2

DefaultAssay(PBMC3) = "ADT"
PBMC3 = NormalizeData(PBMC3, normalization.method = "CLR", margin = 2)

features_adt = rownames(PBMC3@assays[["ADT"]]@counts)

p1 = DimPlot(PBMC2, group.by = "Cluster", reduction = REDUCTION_OF_INTEREST, cols = palette_2) & ggtitle(" Global Visualization") #First Vizualization

cat(" \n \n")
print(p1)
cat(" \n \n")



for (features_to_check in features_adt){
  
  cat(" \n \n")
  cat("### ",  features_to_check)
  cat(" \n \n")
  
  
  
  
  p1 = FeaturePlot(PBMC3, features = features_to_check) &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds" ) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  
  p2 = VlnPlot(PBMC3, features = features_to_check, cols = palette_2, sort = "decreasing") &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds" ) )
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
}


DefaultAssay(PBMC3) = "RNA"
features_RNA  = c("Tcf7", "Tox", "Id2", "Zbtb16", "Cxcr6", "Il18r1", "Il2rb", "Cxcr5", "Ly6a", "Cd160", "Klrc1", "Klrd1", "Flt3" )


VlnPlot(PBMC3, features= features_to_check)

for (features_to_check in features_RNA){
  
  cat(" \n \n")
  cat("### ",  features_to_check)
  cat(" \n \n")
  
  
  p1 = FeaturePlot(PBMC3, features = features_to_check) &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()
  
  cat(" \n \n")
  print(p1)
  cat(" \n \n")
  
  
  p2 = VlnPlot(PBMC3, features = features_to_check, cols = palette_2, sort = "decreasing") &     scale_colour_gradientn(colours = brewer.pal(n = 9, name = "Reds" ) )
  
  cat(" \n \n")
  print(p2)
  cat(" \n \n")
  
}


cat(" \n \n")
cat("## Have a look at Prot score and Transcripto score \n\n")
cat(" \n \n")

#Normalize and scale RNA and ADT

PBMC3 = PBMC2

DefaultAssay(PBMC3)= "RNA"
PBMC3 = NormalizeData(PBMC3)

PBMC3 = ScaleData(PBMC3, do.center =  DATA_CENTER, do.scale = DATA_SCALE, features = rownames(PBMC3) )


#PBMC3$transcri_score


features_RNA_pos  = c( "Tcf7","Tox", "Id2", "Cxcr6", "Klrd1", "Cd160" , "Ifngr1", "Ifng" , "Cd7", "Klrc1", "Klrk1", "Il18r1")
features_RNA_neg = c( "Zbtb16", "Il18r1", "Il2rb", "Cxcr5", "Il2ra","Flt3" )

transcript_score = colSums( PBMC3@assays[["RNA"]]@scale.data[features_RNA_pos,] ) - colSums (PBMC3@assays[["RNA"]]@scale.data[features_RNA_neg,])
PBMC3$transcript_score = transcript_score


p0= FeaturePlot(PBMC3, features= "transcript_score") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()


#Normalize and scale RNA and ADT
DefaultAssay(PBMC3)= "ADT"
PBMC3 =NormalizeData(PBMC3, normalization.method = "CLR", margin = 2)
PBMC3 = ScaleData(PBMC3, do.center =  DATA_CENTER, do.scale = DATA_SCALE, features = rownames(PBMC3) )



#Scoring based on ADT expression in CITEseq DATA
features_ADT_pos  = c( "CD127-ADT" ,  "a4b7-ADT"  )
features_ADT_neg = c( "CD4-ADT", "CD49b-ADT", "NK11-ADT",  "CD3-ADT" , "CCR6-ADT", "CXCR5-ADT" , "CD25-ADT" , "TCRgd-ADT" ,  "TCRb-ADT" , "NKp46-ADT" , "PD1-ADT",  "ST2-ADT"  ,  "NKG2D-ADT" , "KLRG1-ADT"  )


ADT_score = colSums( PBMC3@assays[["ADT"]]@scale.data[features_ADT_pos,] ) - colSums (PBMC3@assays[["ADT"]]@scale.data[features_ADT_neg,])
PBMC3$ADT_score = ADT_score

p1 = FeaturePlot(PBMC3, features= "ADT_score") &     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu" )) ) & NoAxes()



p2 = FeatureScatter(PBMC3, feature1 = "ADT_score", feature2 = "transcript_score", group.by= "Cluster")

cat(" \n \n")
p0
cat(" \n \n")


cat(" \n \n")
p1
cat(" \n \n")


cat(" \n \n")
p2
cat(" \n \n")


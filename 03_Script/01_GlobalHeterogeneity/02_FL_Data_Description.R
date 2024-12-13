#This scripts gives a quick overview of the FL data
#Annotates the populations and save a new object
#Load the data sent by our colaborators
#Load the data sent by our colaborators
load(paste0(PATH_EXPERIMENT_REFERENCE_Data,"/",SAMPLE_FL))

#Rename and save

levels(immune.combined$seurat_clusters) =  c("C1_ILCP", "C2_Klrd1", "C3_ILC3P", "C4_EILP", "C5_Apoe")
immune.combined = SetIdent(immune.combined, value = "seurat_clusters")


#Saving RDS
SAVE_RDS_FILE = paste0(PATH_EXPERIMENT_OUTPUT_RDSFiles, "/FL.rds")
saveRDS(immune.combined, file =  SAVE_RDS_FILE)


#Plots
p0 = DimPlot(immune.combined, label=  FALSE, cols = palette)
p0


png(file=paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures,"/FL_UMAP.png"), width = 25, height = 20,  units = "cm", res=600 )
print(p0)
dev.off()

pdf(file = paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures,"/FL_UMAP.pdf"), width = 25/2.54, height = 20/2.54)
print(p0)
dev.off()

       
         
table(immune.combined@active.ident)

#Diagnostic with their granularity of clustering

All_Markers = FindAllMarkers(immune.combined , only.pos = FINDMARKERS_ONLYPOS, method= FINDMARKERS_METHOD , min.pct =  FINDMARKERS_MINPCT, logfc.threshold = FINDMARKERS_LOGFC_THR , verbose = TRUE)


All_Markers %>%
  group_by(cluster) %>%
  slice_max(n = FINDMARKERS_SHOWTOP, order_by = avg_log2FC) -> top10


p1 = DotPlot(immune.combined, features = unique(top10$gene) , cols = "RdBu") + NoLegend() + theme(axis.text.x = element_text(angle = 90))

p1



png(file=paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures,"/DotPlot_Top",FINDMARKERS_SHOWTOP,"genesHorizontal.png"), width = 55, height = 10,  units = "cm", res=600 )
print(p1)
dev.off()

pdf(file = paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures,"/DotPlot_Top",FINDMARKERS_SHOWTOP,"genesHorizontal.pdf"), width = 55/2.54, height = 10/2.54)
print(p1)
dev.off()




# Generate a vertical DotPlot for the top markers and save it as a PNG and PDF
p2 = DotPlot(immune.combined, features = unique(top10$gene) , cols = "RdBu") + NoLegend()  +  coord_flip()

p2

png(file=paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures,"/DotPlot_Top",FINDMARKERS_SHOWTOP,"genesVertical.png"), width = 10, height = 50,  units = "cm", res=600 ) 
print(p2)
dev.off()


pdf(file = paste0(PATH_EXPERIMENT_OUTPUT_GlobalHeteroFigures,"/DotPlot_Top",FINDMARKERS_SHOWTOP,"Vertical.pdf"), width = 15/2.54, height = 50/2.54)
print(p2)
dev.off()


#Make an additionnal DotPlot of some well known  genes :
Genes_of_Interest_1 = c(Ncr1, Itga2, Klrg1, Il1rl1, Il17a, Il17f ) #Mature ILCs
Genes_of_Interest_1 = c(Mpo, Apoe, Cd63,  Gata1) #C5
p12 = DotPlot(immune.combined, features =  Genes_of_Interest_ILCP , cols = "RdBu") + NoLegend()  +  coord_flip()
p12




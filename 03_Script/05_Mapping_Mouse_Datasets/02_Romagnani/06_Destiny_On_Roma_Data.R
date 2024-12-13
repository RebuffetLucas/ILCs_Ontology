#Install Libraries if necessary
  #Works
#install.packages("sinaplot")
#BiocManager::install("destiny")
#install.packages("rgl")
#install.packages("magick")
#BiocManager::install("genefilter")
#library(Matrix.utils)




#Loading and preparing the data ###

PBMC = PBMC2
DimPlot(PBMC)


#Remove APOE+ subset
#object= PBMC
object= subset(PBMC, idents= SUBSETS_TO_REMOVE , invert = TRUE)
object$Cluster = droplevels(object$Cluster)

object = SetIdent(object, value =  "Cluster")


#Define palette
col_clust_traj=hue_pal()(length(levels(object$Cluster)))
names(col_clust_traj) = levels(object$Cluster)

DimPlot(object, cols = col_clust_traj)


object = FindVariableFeatures(object)


var.features = VariableFeatures(object)

set.seed(SEED)



#Extract data
log.norm.data=object@assays$RNA@data
log.norm.data.var<-log.norm.data[var.features,]
data<-t(as.matrix(log.norm.data.var))

#Extract Metadata
metadata=object@meta.data


#Batch Correction with Diffusion map (FastMNN)

#object_d2 <- RunFastMNN(object.list = SplitObject(object_d, 
#                                                  split.by = "orig.ident"),
#                        features=var.features,
#                        assay = "RNA")

#matrix=t(object_d2@assays$mnn.reconstructed@data) %>% as.matrix()

matrix= data %>% as.matrix()


dm<-DiffusionMap(data=matrix,
                 censor_val = 30,
                 censor_range = c(30, 40),
                 verbose=TRUE
)


plot(dm)

#Save and reload if necessary
#saveRDS(dm , "/mnt/DOSI/EVLAB/BIOINFO/BIOINFO_PROJECT/ILCs_Ontology/05_Output/03_Trancriptionnal_Trajectories/Destiny/Analysis_Objects/dm.rds")
#dm= readRDS("/mnt/DOSI/SUEVLAB/BIOINFO/BIOINFO_PROJECT/Meta_NK5/05_Output/13_Destiny/DiffusionMap_FastMNN_Subset4_NK1and2.rds")

#plot(dm)
dm_coor=dm@eigenvectors


#Root with automatic choice:
dpt_d <- DPT(dm)

#Root with NK2A as initial
group=factor(metadata[rownames(dm_coor),"Cluster"], levels=levels(object$Cluster))
y=aggregate(dm_coor[,1],by=list(group),median)
y2=y[,2] %>% setNames(y[,1])
y2=sort(y2)
if (names(y2)[1]=="EILP") {
  root=which(dm_coor[,1]==max(dm_coor[,1]))
} else {
  root=which(dm_coor[,1]==max(dm_coor[,1]))   
}

dpt_d <- DPT(dm, tips=root)


pt<- dpt_d[[paste0("DPT",root)]]
pt[pt>quantile(pt,0.99,na.rm=TRUE)]=NA
pt[pt<quantile(pt,0.01,na.rm=TRUE)]=NA

df=data.frame("pseudotime"=pt,cluster=group)
rownames(df)=names(dm$DC1)  

#for colors 
pt_color<- dpt_d[[paste0("DPT",root)]]
pt_color[pt_color>quantile(pt_color,0.99,na.rm=TRUE)]=quantile(pt_color,0.99,na.rm=TRUE)
pt_color[pt_color<quantile(pt_color,0.01,na.rm=TRUE)]=quantile(pt_color,0.01,na.rm=TRUE)

clustering_d=metadata[rownames(dm_coor),"Cluster"] %>%
  as.character() %>% setNames(rownames(dm_coor))



# 2D plots with GGPLOT
pt<- dpt_d[[paste0("DPT",root)]]

  #Reverse DC1 if the UMAP is flipped
#dm_coor[,"DC1"] = -dm_coor[,"DC1"]

#Reverse DC2
#dm_coor[,"DC2"] = -dm_coor[,"DC2"]



object_d=subset(object,cells=rownames(dm_coor))
object_d[["dm"]]=CreateDimReducObject(dm_coor, key="DC")
object_d$pt=pt

df2=data.frame(dm_coor[,1:3],pseudotime=pt,pseudotime2=pt_color, Cluster = object_d$Cluster )


#df2$DC1 = -df2$DC1

# DC 1 &2
p1=DimPlot(object_d,reduction="dm", group.by="Cluster", pt.size = 1.4,
           cols = col_clust_traj) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())  + ggtitle("")

p2=ggplot(df2, aes( DC1, DC2,col = pseudotime2)) +
  geom_point() +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_color_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))


p1 + p2


#Avec les contours noirs

p3 = ggplot(df2, aes( DC1, DC2,fill = Cluster)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) + scale_fill_manual(values= col_clust_traj)

p4=ggplot(df2, aes( DC1, DC2,fill = pseudotime2)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))

p3 + p4


#Avec les contours noirs

p6 = ggplot(df2, aes( DC2, DC3,fill = Cluster)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) + scale_fill_manual(values= col_clust_traj)

p7=ggplot(df2, aes( DC2, DC3,fill = pseudotime2)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))

p6 + p7



#Avec les contours noirs

p8 = ggplot(df2, aes( DC1, DC3,fill = Cluster)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) + scale_fill_manual(values= col_clust_traj)

p9=ggplot(df2, aes( DC1, DC3,fill = pseudotime2)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))

p8 + p9


# Pseudotime Plot

p5 = ggplot(df, aes( pseudotime, cluster,col = cluster)) +  
  ggbeeswarm::geom_quasirandom( alpha=0.7,
                                cex=1,
                                show.legend = FALSE,
                                
                                groupOnX=FALSE)+ylab("")  + 
  ggtitle("Pseudotime in Dataset4") + theme_light()+
  stat_summary(
    aes(group = cluster), fun = median, fun.min = median, fun.max = median,
    geom = "crossbar", color = "black", width = 0.7, lwd = 0.2,
    
    # add this bit here to your stat_summary function
    position=position_dodge(width=0.75)
  ) +
  scale_color_manual(values=col_clust_traj)

p5



if (DO_SAVE_FIGURE == TRUE){
  
#Saving Ouputs in HD
png(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC2.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p3 + p4)
dev.off()

png(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC2_DC3.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p6 + p7)
dev.off()

png(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC3.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p8 + p9)
dev.off()


png(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/Plot_PseudoTime.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p5)
dev.off()


#Saving Ouputs in HD
pdf(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC2.pdf"), width = 30/2.54, height = 15/2.54)
print(p3 + p4)
dev.off()


pdf(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC2_DC3.pdf"), width = 30/2.54, height = 15/2.54)
print(p6 + p7)
dev.off()


pdf(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC3.pdf"), width = 30/2.54, height = 15/2.54)
print(p8 + p9)
dev.off()

pdf(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/Plot_PseudoTime.pdf"), width = 30/2.54, height = 15/2.54 )
print(p5)
dev.off()

}


#Try with the scoring
df3=data.frame(dm_coor[,1:3],pseudotime=pt,pseudotime2=pt_color, Cluster = object_d$Cluster , Score_C1_ILCP = object_d$Cluster_C1_ILCP , Score_C2_Klrd1 = object_d$Cluster_C2_Klrd1, Score_C3_ILC3P = object_d$Cluster_C3_ILC3P , Score_C4_EILP = object_d$Cluster_C4_EILP)


p3 = ggplot(df3, aes( DC1, DC3,fill = Cluster)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) + scale_fill_manual(values= col_clust_traj)

p4=ggplot(df3, aes( DC1, DC3,fill = Score_C1_ILCP)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::viridis(11, alpha=1, begin = 0, end=1, direction = 1))

p5 = ggplot(df3, aes( DC1, DC3,fill = Score_C2_Klrd1)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::viridis(11, alpha=1, begin = 0, end=1, direction = 1))

p6 = ggplot(df3, aes( DC1, DC3,fill = Score_C3_ILC3P)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::viridis(11, alpha=1, begin = 0, end=1, direction = 1))


p3 + p4
p5 + p6


if (DO_SAVE_FIGURE == TRUE){
  
  #Saving Ouputs in HD
  png(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC3_ILCP.png"), width = 30, height = 15,  units = "cm", res=600 )
  print(p3 + p4)
  dev.off()
  
  png(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC3_KLRD1_ILC3P.png"), width = 30, height = 15,  units = "cm", res=600 )
  print(p5 + p6)
  dev.off()
  
  
  
  #Saving Ouputs in HD
  pdf(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC3_ILCP.pdf"), width = 30/2.54, height = 15/2.54)
  print(p3 + p4)
  dev.off()
  
  
  pdf(file=paste0(PATH_ROMAGNANI_DESTINY_ANALYSIS_OUTPUT,"/DC1_DC3_KLRD1_ILC3P.pdf"), width = 30/2.54, height = 15/2.54)
  print(p5 + p6)
  dev.off()
  
  
}


## @knitr DiffusionMap_On_ILCProge


cat(" \n \n")
cat("## Diffusion Map analysis on all Human ILC_Proge across tissues with mouse FL signatures")
cat(" \n \n")



#Prepare the data
object = PBMC3
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

object_d2 <- RunFastMNN(object.list = SplitObject(object, 
                                                  split.by = "batch"),
                        features=var.features,
                        assay = "RNA")

matrix=t(object_d2@assays$mnn.reconstructed@data) %>% as.matrix()


#matrix= data %>% as.matrix()


dm<-DiffusionMap(data=matrix,
                 censor_val = 30,
                 censor_range = c(30, 40),
                 verbose=TRUE
)


#plot(dm)



#plot(dm)
dm_coor=dm@eigenvectors


#Root with automatic choice:
dpt_d <- DPT(dm)

#Root with ILCP as initial
group=metadata[rownames(dm_coor),"cluster"]

y=aggregate(dm_coor[,1],by=list(group),median)
y2=y[,2] %>% setNames(y[,1])
y2=sort(y2)
if (names(y2)[1]=="ILCP") {
  root=which(dm_coor[,1]==max(dm_coor[,1]))
} else {
  #root=which(dm_coor[,1]==max(dm_coor[,1]))   
  # Find rows in dm_coor that belong to the 'ILCP' group
  ilcp_indices = which(group == "ILCP")
  # Calculate the root for the 'ILCP' group only
  root = ilcp_indices[which.min(dm_coor[ilcp_indices, 2])]
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

clustering_d=metadata[rownames(dm_coor),"cluster"] %>%
  as.character() %>% setNames(rownames(dm_coor))




# 2D plots with GGPLOT
pt<- dpt_d[[paste0("DPT",root)]]


#Reverse DC1
dm_coor[,"DC1"] = -dm_coor[,"DC1"]

object_d=subset(object,cells=rownames(dm_coor))
object_d[["dm"]]=CreateDimReducObject(dm_coor, key="DC")
object_d$pt=pt

df2=data.frame(dm_coor[,1:3],pseudotime=pt,pseudotime2=pt_color, cluster = object_d$cluster )


#df2$DC1 = -df2$DC1

# DC 1 &2
p1=DimPlot(object_d,reduction="dm", group.by="cluster", pt.size = 1.4,
           cols = palette_cluster) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank())  + ggtitle("")

p2=ggplot(df2, aes( DC1, DC2,col = pseudotime2)) +
  geom_point() +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_color_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))


cat(" \n \n")
cat("### Dim 1 and Dim 2")
cat(" \n \n")



#Avec les contours noirs

p3 = ggplot(df2, aes( DC1, DC2,fill = cluster)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) + scale_fill_manual(values= palette_cluster)

p4=ggplot(df2, aes( DC1, DC2,fill = pseudotime2)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))


cat(" \n \n")
p3 + p4
cat(" \n \n")


#Avec les contours noirs

p5 = ggplot(df2, aes( DC1, DC3,fill = cluster)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) + scale_fill_manual(values= palette_cluster)

p6=ggplot(df2, aes( DC1, DC3,fill = pseudotime2)) +
  geom_point(shape=21) +  theme_bw() + theme(panel.grid.major = element_blank(),   panel.grid.minor = element_blank()) +
  scale_fill_gradientn(colours = viridis::plasma(11, alpha=1, begin = 0, end=1, direction = 1))

cat(" \n \n")
p5 + p6
cat(" \n \n")



cat(" \n \n")
cat("### Pseudotime")
cat(" \n \n")





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
  scale_color_manual(values=palette_cluster)

cat(" \n \n")
p5
cat(" \n \n")





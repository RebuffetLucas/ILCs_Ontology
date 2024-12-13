## @knitr PCA_Across_Tissues_viz
#Only vizualisation of the PCA
## This script performs Principal Component Analysis (PCA) across tissues for visualization purposes. 
## It calculates and visualizes PCA results, including information on variance explained by principal components,
## the contributions of variables to the PCs, and the distribution of samples in PCA space.
## Additional visualizations include bar plots for variable contributions, PCA plots colored by organ and cluster, 
## and a covariance heatmap. The script also includes options to save plots in high resolution for further analysis or presentation.



res.pca <- dudi.pca(mat,
                    scannf = FALSE,   # Cacher le scree plot
                    nf = 5 ,          # Nombre d'axes gardés
                    center= TRUE ,
                    scale= TRUE
)


#Information carried by the PCAs

p_infoPCs = fviz_eig(res.pca, addlabels = TRUE, ylim = c(0,30 ))


cat(" \n \n")
cat("## PCA Across tissues analysis")
cat(" \n \n")


cat(" \n \n")
print(p_infoPCs)
cat(" \n \n")


if (DO_SAVE_FIGURE == TRUE){
png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCA_Plot_Organ.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p_infoPCs)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCA_Plot_Organ.pdf"), width = 15/2.54, height = 15/2.54 )
print(p_infoPCs)
dev.off()



}


#Composition of the PCs
  #Composition of the variables soleil all
p_soleil = fviz_pca_var(res.pca,
                        col.var = "cos2", # Color by contributions to the PC
                        gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                        repel = TRUE,     # Avoid text overlapping,
                        #select.var = list(cos2 = 0.4 )
                        
)

sorted_cos2 <- sort(p_soleil$data$cos2, decreasing = TRUE)
cos2_Treshold =  sorted_cos2[NUMBER_PCA_COMPONENT_psoleil]


if (DO_SAVE_FIGURE == TRUE){
png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/Soleil_AllVar.png"), width = 15, height = 15,  units = "cm", res=600 )
print(p_soleil)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/Soleil_AllVar.pdf"), width = 15/2.54, height = 15/2.54 )
print(p_soleil)
dev.off()

}

  #Composition of the variables soleil top

p_soleil2 = fviz_pca_var(res.pca,
                         col.var = "cos2", # Color by contributions to the PC
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE,     # Avoid text overlapping,
                         select.var = list(cos2 = cos2_Treshold)
                         
)

cat(" \n \n")
print(p_soleil2)
cat(" \n \n")


if (DO_SAVE_FIGURE == TRUE){
png(file= paste0(SAVE_FIG_PATH , "/PCA_Analysis/Soleil_", NUMBER_PCA_COMPONENT_psoleil ,"AllVar.png") , width = 15, height = 15,  units = "cm", res=600 )
print(p_soleil2)
dev.off()

pdf(file= paste0(SAVE_FIG_PATH , "/PCA_Analysis/Soleil_", NUMBER_PCA_COMPONENT_psoleil ,"AllVar.pdf") , width = 15/2.54, height = 15/2.54 )
print(p_soleil2)
dev.off()

}

#Composition of the variables Barplots

for (compteur in 1:5){
  
  plot = fviz_cos2(res.pca, choice = "var", axes = compteur, top = NUMBER_VAR_BARPLOTS_PCA_COMPO )
  
  if (DO_SAVE_FIGURE == TRUE){
  
  png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/VAR_CONTRIB_PCA",compteur, "_TOP", NUMBER_VAR_BARPLOTS_PCA_COMPO ,".png"), width = 30, height = 15,  units = "cm", res=600 )
  print(plot)
  dev.off()
  
  pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/VAR_CONTRIB_PCA",compteur, "_TOP", NUMBER_VAR_BARPLOTS_PCA_COMPO ,".pdf"), width = 30/2.54 , height = 15/2.54 )
  print(plot)
  dev.off()
  
  
  }
}


#PCA visualization

plot_fviz = fviz_pca_ind(res.pca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE # Avoid text overlapping (slow if many points)
)

data_frame_p3 = plot_fviz$data



split_list=str_split(rownames(data_frame_p3), "_", 2)
split_df = data.frame( do.call("rbind", split_list))

data_frame_p3$Organ = split_df$X1
data_frame_p3$group = split_df$X2


#Make PCA plot colored by Organ
p3 = ggplot(data_frame_p3, aes(x = x, y= y, col = Organ,  label= name))+ #can add "shape=group"
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  geom_text_repel(size = PCA_PLOT_TEXT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values = c("#DF6589FF", "#3C1053FF", "#0077b6"  ))+ 
  labs(x = "PC 1", y = "PC 2")+
  ggtitle("Organ")


cat(" \n \n")  
print(p3)
cat(" \n \n")



p3bis = ggplot(data_frame_p3, aes(x = x, y= y, col = Organ, shape = group, label= NULL))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values = c("#DF6589FF", "#3C1053FF", "#0077b6"  ))+ 
  labs(x = "PC 1", y = "PC 2")+
  ggtitle("Organ")




if (DO_SAVE_FIGURE == TRUE){
png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCAbyTissue.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p3)
dev.off()

png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCAbis.png"), width = 30, height = 15,  units = "cm", res=600 )
print(p3bis)
dev.off()


pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCAbyTissue.pdf"), width = 30/2.54, height = 15/2.54  )
print(p3)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCAbis.pdf"), width = 30/2.54, height = 15/2.54  )
print(p3bis)
dev.off()


}

#Covariance Heatmap

covMat = cor(t(mat), method = "spearman") #Spearman = nonparametric, pearson = parametric


data_frame_heatmap = as.data.frame(covMat)
split_list_heatmap=str_split(rownames(data_frame_heatmap ), "_", 2)

split_df_heatmap = data.frame( do.call("rbind", split_list_heatmap))
colnames(split_df_heatmap) = c("Organ", "Cluster")


Pheatmap_plot  = pheatmap(covMat, border_color = "white", fontsize = 10, annotation_row = split_df_heatmap, annotation_colors = my_colourPheatmap)

Pheatmap_plot  = pheatmap(covMat, border_color = "white", fontsize = 10, annotation_row = split_df_heatmap)




cat(" \n \n")
print(Pheatmap_plot)
cat(" \n \n")


if (DO_SAVE_FIGURE == TRUE){
  
png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/Covariance_Heatmap.png"), width = 20, height = 15,  units = "cm", res=600 )
print(Pheatmap_plot)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/Covariance_Heatmap.pdf"), width = 20/2.54, height = 15/2.54 )
print(Pheatmap_plot)
dev.off()


}


#Visualize PC3

#PCA visualization

plot_fviz_PC3 = fviz_pca_ind(res.pca, axes = c(2,3) ,col.ind = "cos2", 
                         gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                         repel = TRUE # Avoid text overlapping (slow if many points)
)

data_frame_p3_PC3 = plot_fviz_PC3$data
split_list=str_split(rownames(data_frame_p3_PC3), "_", 2)
split_df = data.frame( do.call("rbind", split_list))


data_frame_p3_PC3$Organ = split_df$X1
data_frame_p3_PC3$group = split_df$X2



p4_PC3 = ggplot(data_frame_p3_PC3, aes(x = x, y= y, col = Organ, label= name))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  geom_text_repel(size = PCA_PLOT_TEXT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values =  c("#F8766D" , "#8494FF", "#0CB702" ))+ 
  labs(x = "PC 2", y = "PC 3")+
  ggtitle("Organ")

cat(" \n \n")
print(p4_PC3)
cat(" \n \n")


p4bis_PC3 = ggplot(data_frame_p3_PC3, aes(x = x, y= y, col = Organ, label= NULL))+
  geom_point(size = PCA_PLOT_POINT_SIZE)+
  theme_bw()+
  geom_hline(yintercept = 0, linetype= "dashed")+
  geom_vline(xintercept = 0, linetype= "dashed")+
  scale_color_manual(values = c("#F8766D" , "#8494FF", "#0CB702" ))+ 
  labs(x = "PC 2", y = "PC 3")+
  ggtitle("Organ")

cat(" \n \n")
print(p4bis_PC3)
cat(" \n \n")

if (DO_SAVE_FIGURE == TRUE){
png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCA_ColbyCluster_PC3.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p4_PC3)
dev.off()

png(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCA_ColbyCluster_NoNAME_PC3.png"), width = 20, height = 15,  units = "cm", res=600 )
print(p4bis_PC3)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCA_ColbyCluster_PC3.pdf"), width = 20/2.54, height = 15/2.54 )
print(p4_PC3)
dev.off()

pdf(file=paste0(SAVE_FIG_PATH , "/PCA_Analysis/PCA_ColbyCluster_NoNAME_PC3.pdf"), width = 20/2.54, height = 15/2.54  )
print(p4bis_PC3)
dev.off()
}


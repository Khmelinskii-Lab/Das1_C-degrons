##### for project : Analysis of degron at C terminal
##### sub text:clutsering SHAP data to infer degrons
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB



# set paths
setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")

path_plot<-paste0(getwd(),"/plot/")
path_dataset<-paste0(getwd(),"/output/")
path_shap_data<-paste0(path_dataset,"/shap_instable/")


# read the dataset 

#read the instable shap data
instable_shap<-read.csv(paste0(path_shap_data,"shap_sim_phi_test_class_2_instable_12x_for_instable_all.csv"))
row.names(instable_shap)<-instable_shap$X
instable_shap$X<-NULL

instable_shap<-as.data.frame((t(instable_shap)))

#read the normalized dataset 
testing_dataset<-read.csv(paste0(path_dataset,"testing_data_after_normalization_for_model.csv"))
training_dataset<-read.csv(paste0(path_dataset,"training_data_after_normalization_for_model.csv"))
row.names(testing_dataset)<-testing_dataset$X
row.names(training_dataset)<-training_dataset$X
testing_dataset$X<-NULL
training_dataset$X<-NULL
dataset<-as.data.frame(rbind(testing_dataset,training_dataset))
# get the instable normaliozed data 
normalized_data_instable<-subset(dataset, row.names(dataset) %in% row.names(instable_shap))


# order both the dataset by their row names
normalized_data_instable<-normalized_data_instable[order(row.names(normalized_data_instable)),]
instable_shap<-instable_shap[order(row.names(instable_shap)),]

names(instable_shap)[277:288]<-c("V1_potentialProteinInteractionIndex","V2_potentialProteinInteractionIndex","V3_potentialProteinInteractionIndex","V4_potentialProteinInteractionIndex",      
                                 "V5_potentialProteinInteractionIndex","V6_potentialProteinInteractionIndex","V7_potentialProteinInteractionIndex","V8_potentialProteinInteractionIndex" ,     
                                 "V9_potentialProteinInteractionIndex","V10_potentialProteinInteractionIndex","V11_potentialProteinInteractionIndex","V12_potentialProteinInteractionIndex" )

names(normalized_data_instable)[277:288]<-c("V1_potentialProteinInteractionIndex","V2_potentialProteinInteractionIndex","V3_potentialProteinInteractionIndex","V4_potentialProteinInteractionIndex",      
                                            "V5_potentialProteinInteractionIndex","V6_potentialProteinInteractionIndex","V7_potentialProteinInteractionIndex","V8_potentialProteinInteractionIndex" ,     
                                            "V9_potentialProteinInteractionIndex","V10_potentialProteinInteractionIndex","V11_potentialProteinInteractionIndex","V12_potentialProteinInteractionIndex" )


names(instable_shap)<-paste0(names(instable_shap),"_Contribution_Score")
dataset<-instable_shap

# start analysis

fviz_nbclust(instable_shap, kmeans, method = "wss", k.max = 150)
k2 <- kmeans(instable_shap, centers = 56, nstart = 150)
k2_cluster<-as.data.frame(k2$cluster)

write.csv(k2_cluster,
          paste0(path_dataset,"k2cluster.csv"))

# for finding representative peptides

k2_cluster<-read.csv(paste0(path_dataset,"k2cluster.csv"))




datasetcomplete<-read.csv(paste0(path_dataset,"complete_dataset_without_normalization.csv"))
dataset<-read.csv(paste0(path_dataset,"merged_file_post_remoal.csv"))

dataset<-subset(dataset,dataset$Codon == "Non Stop Codon")
names(k2_cluster)<-c("translation","k2.cluster")

d<-datasetcomplete[,c("hydrophobicity_overall","PSI")]
d$translation<-datasetcomplete$X
d<-subset(d,d$translation %in% dataset$raw_counts_translation)
d$Label<-0
d$Label<-ifelse(d$PSI < 0.5, 1,d$Label)


d1<-merge(d,k2_cluster, by = "translation",  all.x = TRUE)
d1$k2.cluster[is.na(d1$k2.cluster)]<-0

d$hydrophobicity_overall<-as.numeric(d$hydrophobicity_overall)
d$PSI<-as.numeric(d$PSI)
d$Label<-as.factor(d$Label)


#k2_cluster$translation<-row.names(k2_cluster)
d_unstable<-subset(d,d$translation %in% k2_cluster$translation)
d_unstable<-merge(k2_cluster,d_unstable, by = "translation")
d_unstable$hydrophobicity_overall<-as.numeric(d_unstable$hydrophobicity_overall)
d_unstable$PSI<-as.numeric(d_unstable$PSI)
names(d_unstable)<-c("translation","cluster","hydrophobicity","PSI","Label")

median_hydrophobicity<-d_unstable  %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(median=median(hydrophobicity),
                   mean = mean(hydrophobicity))

median_PSI<-d_unstable  %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarize(median=median(PSI),
                   mean = mean(PSI))

total_peptides<-as.data.frame(table(d_unstable$cluster))

names(median_hydrophobicity)<-c("Cluster","Median_Hydrophobicity","Mean_Hydrophobicity")
names(median_PSI)<-c("Cluster","Median_PSI","Mean_PSI")
names(total_peptides)<-c("Cluster","total_peptides")
names(d_unstable)[2]<-"Cluster"

d_unstable<-merge(d_unstable,median_hydrophobicity, by = "Cluster")
d_unstable<-merge(d_unstable,median_PSI, by = "Cluster")
d_unstable<-merge(d_unstable,total_peptides, by = "Cluster")
d_unstable_sub<-d_unstable[,c("Cluster","Mean_Hydrophobicity","Mean_PSI","total_peptides")]
d_unstable_sub<-unique(d_unstable_sub)


similarity<-as.data.frame(matrix(ncol = 2))
for (i in 1:max(k2_cluster$k2.cluster)) {
  
  a<-as.data.frame(parSeqSim(subset(d_unstable,d_unstable$Cluster == i)[,c("translation")]))
  a$sum<-rowSums(a)
  row.names(a)<-subset(d_unstable,d_unstable$Cluster == i)[,c("translation")]
  s<-c(i,row.names(a[a$sum == max(a$sum),]))
  similarity<-as.data.frame(rbind(similarity,s))
  print(i)
}

s<-similarity[2:nrow(similarity),]
write.csv(s,paste0(path_dataset,"similarity.csv"))
s<-read.csv(paste0(path_dataset,"similarity.csv"))
s$X<-NULL
d<-datasetcomplete[,c("hydrophobicity_overall","PSI")]
d$translation<-datasetcomplete$X
d<-subset(d,d$translation %in% dataset$raw_counts_translation)
d$Label<-0
d$Label<-ifelse(d$PSI < 0.5, 1,d$Label)

d$hydrophobicity_overall<-as.numeric(d$hydrophobicity_overall)
d$PSI<-as.numeric(d$PSI)
d$Label<-as.factor(d$Label)
d_unstable<-subset(d,d$translation %in% k2_cluster$translation)

#k2_cluster$translation<-row.names(k2_cluster)
d_unstable<-merge(k2_cluster,d_unstable, by = "translation")
d_unstable$hydrophobicity_overall<-as.numeric(d_unstable$hydrophobicity_overall)
d_unstable$PSI<-as.numeric(d_unstable$PSI)

names(d_unstable)<-c("translation","cluster","hydrophobicity","PSI","Label")
total_peptides<-as.data.frame(table(d_unstable$cluster))
names(total_peptides)<-c("cluster","TotalPeptides")
d_unstable_sub<-subset(d_unstable,d_unstable$translation %in% s$V2)
d_unstable_sub<-merge(d_unstable_sub,total_peptides, id = "cluster")

d_unstable_sub$peptide_group<-300
d_unstable_sub$peptide_group<-ifelse(d_unstable_sub$TotalPeptides <= 200, 200, d_unstable_sub$peptide_group)
d_unstable_sub$peptide_group<-ifelse(d_unstable_sub$TotalPeptides <= 100, 100, d_unstable_sub$peptide_group)
d_unstable_sub$peptide_group<-ifelse(d_unstable_sub$TotalPeptides <= 50, 50, d_unstable_sub$peptide_group)
d_unstable_sub$peptide_group<-ifelse(d_unstable_sub$TotalPeptides <= 10, 10, d_unstable_sub$peptide_group)
d_unstable_sub$peptide_group<-as.factor(d_unstable_sub$peptide_group)

pdf(paste0(path_plot,"FigureS2/highlighted_cluster.pdf"))

ggplot(d_unstable)+
  geom_point(data = d,aes(x = PSI,y = hydrophobicity_overall), color = "grey", alpha = 0.2)+
  geom_point(data = d_unstable_sub,
             aes(x = PSI,
                 y = hydrophobicity, 
                 size = peptide_group 
             ), color = "red", alpha = 0.4)+
  geom_point(data = subset(d_unstable_sub,d_unstable_sub$cluster %in% c("11","31","45")), 
             aes(x = PSI,
                 y = hydrophobicity,
                 size = peptide_group
             ), 
             color = "blue", 
             alpha = 0.4)+
  geom_text_repel(data = d_unstable_sub,
                  aes(x = PSI,
                      y= hydrophobicity,
                      label = cluster),max.overlaps = 15,
                  segment.color = 'grey50'  )+
  ylim(-4,4)+
  xlim(0,1)+
  theme_bw()+
  ylab("Hydrophobicity")+
  xlab("PSI")+
  labs(title = "Hydrophobicity versus PSI with for Unstable clusters \n represented by the most representative peptide in cluster",
       caption = paste0 ("# Unstable peptides : 4,110
  # Total peptides : 46,152,
                         Red: mean Hydrophobicity versus median PSI for each cluster
                         blue: mean Hydrophobicity versus median PSI for Clusters ending with N (#Peptides: 74)"
       ))



dev.off()

# formation of seqlogo based on contribution score
instable_shap<-instable_shap[order(row.names(instable_shap)),]
normalized_data_instable<-normalized_data_instable[order(row.names(normalized_data_instable)),]
seq<-instable_shap[,1:240] * normalized_data_instable[,1:240]



path_cluster<-paste0(path_plot,"FigureS2/cluster_instable_12x/")

cs1 = make_col_scheme(chars=c("H","M","A","P","T","B"), groups=c('Hydrophobicity','MOC', "AliphaticIndex",'PPIndex', 'AlphaTurn','Bulkiness'), 
                      cols=c('black', 'brown', 'orange','red', 'blue','darkgreen'))








for(i in 1:max(k2_cluster$k2.cluster)){
  
  #i<-1
  #seq_required<-subset(instable_shap,row.names(instable_shap) %in% row.names(subset(k2_cluster,k2_cluster$`k2$cluster` == i)))
  seq_required<-subset(seq,row.names(seq) %in% (subset(k2_cluster,k2_cluster$k2.cluster == i))[,"X"])
  column_means<-colMeans(seq_required)
  
  data_seqlogo<-data.frame(matrix(t(column_means[1:240]), nrow=20, byrow=FALSE))
  names(data_seqlogo)<-as.character(-12:-1)
  row.names(data_seqlogo)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
  #names(data_seqlogo)<-as.numeric(names(data_seqlogo))-13
  
  a<- print(ggseqlogo(as.matrix(data_seqlogo), method='custom', seq_type='aa') +
              theme(axis.text.x = element_blank(),
                    axis.ticks.y = element_line(size =1))+
              ylab('CS Mean')+
              xlab("Position")+
              annotate(geom = "text",
                       x=2.5,
                       y = max(data_seqlogo)+0.015, 
                       label = paste0("Motif : ", i, " , \n # peptides : ",nrow(seq_required)))+
              ylim(min(data_seqlogo)-0.05, max(data_seqlogo)+0.1)+
              annotate(geom = 'segment', y = 0, yend = 0, color = 'black', x = 1, xend = Inf, size = 0.25) )
  
  seq_required<-subset(instable_shap,row.names(instable_shap) %in% (subset(k2_cluster,k2_cluster$k2.cluster == i))[,"X"])
  column_means<-colMeans(seq_required)
  data_seqlogo1<-data.frame(matrix(t(column_means[241:312]), nrow=12, byrow=FALSE))
  names(data_seqlogo1)<-c("H","M","A","P","T","B")
  
  row.names(data_seqlogo1)<-as.numeric(row.names(data_seqlogo1))-13
  data_seqlogo1<-as.data.frame(t(data_seqlogo1))
  #pdf(paste0("Y:/lab data/susmitha/edwin/for_paper/plot/instable/with_all_property/cluster_motif_positional_biophysics/",i,".pdf"))
  b <-(ggseqlogo(as.matrix(data_seqlogo1), method='custom', seq_type='aa',col_scheme=cs1)) +
    theme(axis.text.x = element_blank(),
          axis.ticks.y = element_line(size =1))+
    #theme_logo(base_size = 10)+
    ylab('CS Mean')+
    xlab("Position")+
    ggtitle("Positional biophysical property")+
    ylim(min(data_seqlogo1)-0.1, max(data_seqlogo1)+0.1)
  
  seq_required<-subset(instable_shap,row.names(instable_shap) %in% (subset(k2_cluster,k2_cluster$k2.cluster == i))[,"X"])
  overall_biophysical<-seq_required[313:318]
  names(overall_biophysical)<-c('Hydrophobicity','MOC', "AliphaticIndex",'PPIndex', 'AlphaTurn','Bulkiness')
  overall_biophysical<-melt(overall_biophysical)
  #overall_biophysical$variable<-str_split_fixed(overall_biophysical$variable,"_",4)[,1]
  names(overall_biophysical)<-c("Property","Value")
  
  c<- overall_biophysical %>% ggplot(aes(x = Property, y = Value))+
    geom_boxplot()+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90))+
    ylim(-0.07,0.06)+
    ylab('CS')+
    ggtitle("Overall biophysical property")
  
  d<- ggplot(d1)+
    geom_point(data = d1, aes(
      x = PSI,
      y = hydrophobicity_overall
    ), color = "grey")+
    geom_point(data = subset(d1,d1$k2.cluster == i), aes(
      x = PSI,
      y = hydrophobicity_overall
    ), color = "red")+
    ylim(-4,4)+
    xlim(0,1)+
    ylab("Hydrophobicity")+
    xlab("Protein Stability Index")+
    theme_bw()
  
  pdf(paste0(path_cluster,i,".pdf"))
  print(
    plot_grid(a,c,b,d)
  )
  dev.off()
  print(i)
}


# 
k2_cluster<-read.csv(paste0(path_dataset,"k2cluster.csv"))
k2_cluster$translation<-k2_cluster$X
k2_cluster$X<-NULL




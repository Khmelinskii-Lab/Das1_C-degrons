##### for project : Analysis of degron at C terminal
##### sub text:yeast Cterminal
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB
setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")


yeast_CTer<-read.csv(
          paste0(path_output,"/yeast_Cter_all.csv"))
yeast_CTer = yeast_CTer[!duplicated(yeast_CTer$Translation),]
yeast_CTer<-yeast_CTer[complete.cases(yeast_CTer),]

M<-cor(as.matrix(yeast_CTer[,c("WT_replicate2","WT_replicate3","Das1_replicate2","Das1_replicate3" )]))

M<-as.data.frame(M)
M$dataset<-row.names(M)
M_melt<-melt(M, id = "dataset")
M_melt$variable<-factor(M_melt$variable,
                        levels = (c("WT_replicate2",
                                    "WT_replicate3",
                                    "Das1_replicate3"
                        )))
M_melt$dataset<-factor(M_melt$dataset,
                       levels = (c("WT_replicate2",
                                   "WT_replicate3",
                                   "Das1_replicate2",
                                   "Das1_replicate3"
                       )))
#path_plot<-"Y:/lab data/susmitha/edwin/for_paper/new_dataset/plots/new_dataset_plots/CTer_WT/"
pdf(paste0(path_plot,"FigureS6/correlation_among_CTer_dataset.pdf"))

ggplot(M_melt)+
  geom_tile(data = M_melt, aes(x = dataset, y = variable , fill = value))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_fill_gradientn("Correlation",colours=c("#F7931E", 
                                               "#F5AB53", 
                                               "#F4BB75", 
                                               "#F3CB98", 
                                               "#F1DBBB", 
                                               "#F0EBDE",
                                               "#FFFFFF",
                                               "#C4C0E4", 
                                               "#9996EA", 
                                               "#6D6BF0", 
                                               "#4240F6", 
                                               "#0000FF"),
                       values=scales::rescale(c(-1.00000000 ,-0.81, -0.63, -0.45, -0.27, -0.09, 0 ,0.09,  0.27,  0.45,  0.63,  0.81,1.00000000)),
                       limits=c(-1,1))+
  xlab("")+
  ylab("")
dev.off()


row.names(yeast_CTer)<-yeast_CTer$Translation
q<-as.data.frame(yeast_CTer[,c("WT_replicate2","WT_replicate3")])
WT<-lmFit(q[,c(1,2)])
mean_result_WT<-q
mean_result_WT$effect<-WT$sigma * WT$stdev.unscaled
mean_result_WT$ci_95<-WT$coefficients * qt(0.975 , WT$df.residual)
mean_result_WT$value<-WT$coefficients

p<-as.data.frame(yeast_CTer[,c("Das1_replicate2","Das1_replicate3")])
DT<-lmFit(p[,c(1,2)])
mean_result_DT<-p
mean_result_DT$effect<-DT$sigma * DT$stdev.unscaled
mean_result_DT$ci_95<-DT$coefficients * qt(0.975 , DT$df.residual)
mean_result_DT$value<-DT$coefficients


mean_result_DT$translation<-row.names(mean_result_DT)
mean_result_WT$translation<-row.names(mean_result_WT)
WT_Das1<-merge(mean_result_WT, mean_result_DT, by= "translation")

#pdf(paste0(path_plot,"WT_DAs1_Cter.pdf"))
ggplot(WT_Das1)+
  geom_point(aes(x = value.x, y = value.y))+
  theme_bw()+
  xlab("Wild type ")+
  ylab("Das1")+
  xlim(0,1)+
  ylim(0,1)+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(WT_Das1$value.x,WT_Das1$value.y),2)),
           colour = "red")+
  ggtitle("PSI WT versus Das1 in yeast C terminal after lmFit")
#dev.off()
q_data<-WT_Das1[,c("translation" ,"WT_replicate2","WT_replicate3","Das1_replicate2","Das1_replicate3")]
q_data<-melt(q_data,id = "translation")
q_data$group<-str_split_fixed(q_data$variable,"_",2)[,1]
stat.test <- q_data %>%
  group_by(translation) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

WT_Das1<-merge(stat.test, WT_Das1, id  = "translation")


WT_Das1$difference<-WT_Das1$value.y - WT_Das1$value.x
WT_Das1$class<-ifelse(WT_Das1$value.x < 0.5,"Unstable","Stable")

WT_Das1$color<-ifelse(((WT_Das1$difference) > 0.1 & (WT_Das1$class == "Unstable")), "Potential degron","No effect")
WT_Das1$category<-"No Effect"
WT_Das1$category<-ifelse((WT_Das1$color == "Potential degron" & 
                            WT_Das1$value.y < 0.58),"Intermediate", WT_Das1$category)
WT_Das1$category<-ifelse((WT_Das1$color == "Potential degron" & WT_Das1$value.y > 0.58 ),"Stabilized", WT_Das1$category)
WT_Das1$category<-ifelse((WT_Das1$color == "Potential degron" & WT_Das1$value.x > 0.5 ),"No Effect", WT_Das1$category)
WT_Das1$category<-ifelse(((WT_Das1$p > 0.05) & (WT_Das1$category %in% c("Intermediate","Stabilized"))),
                         "Unstable but no effect" ,WT_Das1$category)
WT_Das1$category<-ifelse(((WT_Das1$color == "No effect") & WT_Das1$value.x < 0.5 ),"Unstable but no effect", WT_Das1$category)


#stabilized: difference>0.1 and psi_das1>0.58,Intermediate: difference>0.1 and psi_das1<0.58(if p value > 0.05, stabilized and intermediate are put as no effect)
#No effect: difference>0.1 and psi_wt>0.58,unstable but no effect: difference>0.1 and psi_wt<0.58

#data_with_stat_test<-merge(stat.test, WT_Das1, id  = "translation")
#data_with_stat_test$category<-ifelse(data_with_stat_test$p > 0.05,"No Effect" ,data_with_stat_test$category)

pdf(paste0(path_plot,"Figure6/WT_Das1_4Category.pdf"))

WT_Das1 %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = category))+
  theme_bw()+
  scale_color_manual(values = c( "Stabilized" = "#ff5400",
                                 "Intermediate" = "#8d99ae",
                                 "No Effect" = "#2b2d42",
                                 "Unstable but no effect" = "#5C2751"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (DAS1)")+
  ggtitle("Categories in yeast C Terminal without multiple testing but significance testing")+
  annotate("text", 
           x = 0.15, y = 0.92, 
           label = paste0("Threshold : \n Difference :0.1, \nPSI WT: 0.5; \nPSI Das1 : 0.58"),
           colour = "red", fontface =2)
dev.off()

WT_Das1$hydrophobicity<-hydrophobicity(WT_Das1$translation)
my_comparisons <- list( c("No Effect","Unstable but no effect"), 
                        c("Stabilized","Intermediate"), 
                        c("No Effect", "Intermediate"),
                        c("No Effect","Stabilized"),
                        c("Stabilized","Unstable but no effect"),
                        c("Unstable but no effect","Intermediate")
)
pdf(paste0(path_plot,"Figure6/WT_Das1_hydrophobicity_4Category.pdf"))
ggplot(WT_Das1)+
  geom_boxplot(
    aes(x = category, 
        y = hydrophobicity),
    coef = 1e30,lwd=0.25
  )+theme_bw()+
  scale_color_manual(values = c( "Stabilized" = "#ff5400",
                                 "Intermediate" = "#8d99ae",
                                 "No Effect" = "#2b2d42",
                                 "Unstable but no effect" = "#5C2751"))+
  ylim(-4,4)+
  guides(col = FALSE)+
  ylab("Hydrophobicity")+
  xlab("")+
  ggtitle("Variation of hydrophobicity")+
  theme_bw()+ 
  stat_compare_means(data = WT_Das1, 
                     aes(x = category, 
                         y = hydrophobicity),
                     method = "wilcox.test",ref.group = "Unstable but no effect", label.y = 4)+
  stat_compare_means(data = WT_Das1, 
                     aes(x = category, 
                         y = hydrophobicity),
                     label = "p.signif", method = "wilcox.test",
                     ref.group = "Unstable but no effect",
                     label.y = 3.5)


dev.off()


# frequency 
dataset<-dataset<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))
dataset<-subset(dataset,dataset$Codon == "Non Stop Codon")
names(dataset)
dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]
dataset_for_frequency<-unique(dataset_for_frequency)
dataset_for_frequency<-dataset_for_frequency[!duplicated(dataset_for_frequency$raw_counts_translation),]
dataset_for_frequency_oneHot<-create_one_hot_db(dataset_for_frequency$raw_counts_translation,12)


CTer_data<-create_one_hot_db(WT_Das1$translation,12)
amino_acid<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
t<-strrep(amino_acid, 12)  

stabilized<-subset(WT_Das1,WT_Das1$category == "Stabilized")
count<-nrow(stabilized)
translation<-c(t,stabilized$translation)
stabilized<-create_one_hot_db(translation,12)
stabilized<-stabilized[21:nrow(stabilized),]
stabilized<-as.data.frame(colSums(stabilized))
stabilized$frequency<-stabilized$`colSums(stabilized)`/count
stabilized$allFreq<-as.data.frame(colSums(CTer_data))[,1]/6403
stabilized$log2NormFreq<-log2(stabilized$frequency/stabilized$allFreq)
stabilized$identifier<-row.names(stabilized)
stabilized$position<-substr(stabilized$identifier,2,(nchar(stabilized$identifier)-1))
stabilized$position<-as.numeric(stabilized$position)-13
stabilized$amino_acid<-substr(stabilized$identifier,(nchar(stabilized$identifier)),(nchar(stabilized$identifier)))
stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
stabilized$amino_acid<-factor(stabilized$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))

stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"Figure6/Stabilized_Frequency_Cter.pdf"))
print(
  plot_grid(
    
    ggplot(stabilized, aes (x = position,
                            y= amino_acid,
                            fill = log2NormFreq))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreq < 2 & 
                                stabilized$log2NormFreq > -2), 
                aes(
                  x = position,
                  y= amino_acid,
                  fill = log2NormFreq
                ))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreq < -2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreq > 2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#0000FF")+
      geom_tile(data = subset(stabilized,is.na(stabilized$log2NormFreq)), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_text(aes(label = ifelse(is.na(stabilized$log2NormFreq), "*", ""),), 
                color = "Black", 
                size = 2) +
      coord_equal()+ theme_bw()+
      theme(text = element_text(size = 8))+
      
      scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_6403)",
                           na.value = "grey",
                           colours=c("#F7931E",
                                     "#F5AB53",
                                     "#F4BB75", 
                                     "#F3CB98", 
                                     "#F1DBBB", 
                                     "#F0EBDE",
                                     "#FFFFFF",
                                     "#C4C0E4", 
                                     "#9996EA", 
                                     "#6D6BF0", 
                                     "#4240F6", 
                                     "#0000FF"),
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    
    
    
    ggplot(stabilized, aes (x = position,
                            y= amino_acid,
                            fill = log2NormFreqAll))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreqAll < 2 & 
                                stabilized$log2NormFreqAll > -2), 
                aes(
                  x = position,
                  y= amino_acid,
                  fill = log2NormFreqAll
                ))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreqAll < -2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreqAll > 2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#0000FF")+
      geom_tile(data = subset(stabilized,is.na(stabilized$log2NormFreqAll)), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_text(aes(label = ifelse(is.na(stabilized$log2NormFreqAll), "*", ""),), 
                color = "Black", 
                size = 2) +
      coord_equal()+
      
      theme_bw()+
      theme(text = element_text(size = 8))+
      scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_12x)",
                           na.value = "grey",
                           colours=c("#F7931E",
                                     "#F5AB53",
                                     "#F4BB75", 
                                     "#F3CB98", 
                                     "#F1DBBB", 
                                     "#F0EBDE",
                                     "#FFFFFF",
                                     "#C4C0E4", 
                                     "#9996EA", 
                                     "#6D6BF0", 
                                     "#4240F6", 
                                     "#0000FF"),
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    
    ncol = 1
  )
)
dev.off()



Intermediate<-subset(WT_Das1,WT_Das1$category == "Intermediate")
count<-nrow(Intermediate)
translation<-c(t,Intermediate$translation)
Intermediate<-create_one_hot_db(translation,12)
Intermediate<-Intermediate[21:nrow(Intermediate),]
Intermediate<-as.data.frame(colSums(Intermediate))

Intermediate$frequency<-Intermediate$`colSums(Intermediate)`/count
Intermediate$allFreq<-as.data.frame(colSums(CTer_data))[,1]/6403
Intermediate$log2NormFreq<-log2(Intermediate$frequency/Intermediate$allFreq)
Intermediate$NormFreq<-(Intermediate$frequency/Intermediate$allFreq)
Intermediate$identifier<-row.names(Intermediate)
Intermediate$position<-substr(Intermediate$identifier,2,(nchar(Intermediate$identifier)-1))
Intermediate$position<-as.numeric(Intermediate$position)-13
Intermediate$amino_acid<-substr(Intermediate$identifier,(nchar(Intermediate$identifier)),(nchar(Intermediate$identifier)))
Intermediate$position<-factor(Intermediate$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
Intermediate$amino_acid<-factor(Intermediate$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
Intermediate[Intermediate$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


Intermediate$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
Intermediate$log2NormFreqAll<-log2(Intermediate$frequency/Intermediate$freqAllData)
Intermediate$NormFreqAll<-(Intermediate$frequency/Intermediate$freqAllData)

Intermediate[Intermediate$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"Figure6/Intermediate_Frequency_Cter.pdf"))
print(
  plot_grid(
    ggplot(Intermediate, aes (x = position,
                              y= amino_acid,
                              fill = log2NormFreq))+
      geom_tile(data = subset(Intermediate,Intermediate$log2NormFreq < 2 & 
                                Intermediate$log2NormFreq > -2), 
                aes(
                  x = position,
                  y= amino_acid,
                  fill = log2NormFreq
                ))+
      geom_tile(data = subset(Intermediate,Intermediate$log2NormFreq < -2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_tile(data = subset(Intermediate,Intermediate$log2NormFreq > 2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#0000FF")+
      geom_tile(data = subset(Intermediate,is.na(Intermediate$log2NormFreq)), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_text(aes(label = ifelse(is.na(Intermediate$log2NormFreq), "*", ""),), 
                color = "Black", 
                size = 2) +
      coord_equal()+
      theme_bw()+
      theme(text = element_text(size = 8))+
      scale_fill_gradientn("log2(Freq_Intermediate)/\n (Freq_6403)",
                           na.value = "grey",
                           colours=c("#F7931E",
                                     "#F5AB53",
                                     "#F4BB75", 
                                     "#F3CB98", 
                                     "#F1DBBB", 
                                     "#F0EBDE",
                                     "#FFFFFF",
                                     "#C4C0E4", 
                                     "#9996EA", 
                                     "#6D6BF0", 
                                     "#4240F6", 
                                     "#0000FF"),
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    
    
    ggplot(Intermediate, aes (x = position,
                              y= amino_acid,
                              fill = log2NormFreqAll))+
      geom_tile(data = subset(Intermediate,Intermediate$log2NormFreqAll < 2 & 
                                Intermediate$log2NormFreqAll > -2), 
                aes(
                  x = position,
                  y= amino_acid,
                  fill = log2NormFreqAll
                ))+
      geom_tile(data = subset(Intermediate,Intermediate$log2NormFreqAll < -2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_tile(data = subset(Intermediate,Intermediate$log2NormFreqAll > 2 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#0000FF")+
      geom_tile(data = subset(Intermediate,is.na(Intermediate$log2NormFreqAll)), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_text(aes(label = ifelse(is.na(Intermediate$log2NormFreqAll), "*", ""),), 
                color = "Black", 
                size = 2) +
      coord_equal()+
      theme_bw()+
      theme(text = element_text(size = 8))+
      scale_fill_gradientn("log2(Freq_Intermediate)/\n (Freq_12x)",
                           na.value = "grey",
                           colours=c("#F7931E",
                                     "#F5AB53",
                                     "#F4BB75", 
                                     "#F3CB98", 
                                     "#F1DBBB", 
                                     "#F0EBDE",
                                     "#FFFFFF",
                                     "#C4C0E4", 
                                     "#9996EA", 
                                     "#6D6BF0", 
                                     "#4240F6", 
                                     "#0000FF"),
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    ncol=1
  )
)
dev.off()

write.csv(WT_Das1,paste0(path_output,"WT_Das1_Cter.csv"))

WT_Das1<-read.csv(paste0(path_output,"WT_Das1_Cter.csv"))
#hydrophobicity versus PSI

# plot for distibution of hydrophobicity versus PSI
PSI_hydro<-WT_Das1 %>%
  ggplot(aes(x = value.x ,y =hydrophobicity))+ 
  geom_point(aes(color = class),alpha = 0.8, shape = 19)+
  theme_bw()+
  scale_color_manual(values=c("#4290C0","#999999"))+
  #geom_smooth(method=lm, se=FALSE, linetype = "dashed", color = "black")+
  xlab("Protein Stability Index")+
  ylab(" Overall Hydrophobicity")+
  xlim(0,1)+
  ylim(-4,4)+
  annotate("text", x = 0.1, y = 3.9, label = paste0("R : ", round(cor(
    WT_Das1$hydrophobicity,
    WT_Das1$value.x, method = c("spearman")
  ),3)))

densi<-WT_Das1 %>%
  ggplot(aes(y = hydrophobicity, color = class))+ 
  geom_density(kernel = "rectangular")+
  theme_bw()+
  scale_color_manual(values=c( "#4290C0","#999999"))+
  xlab("Density")+
  ylab(" Overall Hydrophobicity")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ylim(-4,4)

PSI_hydro <- PSI_hydro + rremove("legend")
#densi <- densi + 
pdf(paste0(path_plot,"FigureS6/Hydrophobicity_PSI_Cter.pdf"))
plot_grid(PSI_hydro,densi)
dev.off()


yeastmine2<-read.csv(paste0(path_input,"yeastmine_2.csv"), sep = ",", header = FALSE)
names(yeastmine2)<-c("GeneID","Systematic_Name","Organism","Standard_Name","Gene_Name","Protein_Standard_Name","Protein_Systematic_Name","Gene_Qualifier","Feature_Type","Residue")
yeastmine2$Residue<-substr(yeastmine2$Residue,1,(nchar(yeastmine2$Residue)-1))
yeastmine2$translation<-substr(yeastmine2$Residue,(nchar(yeastmine2$Residue)-11),(nchar(yeastmine2$Residue)))

yeastmine2<-yeastmine2[,c("Protein_Systematic_Name","Gene_Qualifier","translation")]
ORF<-merge(WT_Das1,yeastmine2, by = "translation")

stabilized_ORF<-subset(ORF,ORF$category == "Stabilized")
Intermediate_ORF<-subset(ORF,ORF$category == "Intermediate")
NoEffect_ORF<-subset(ORF,ORF$category == "No Effect")
UnstableButNoEffect_ORF<-subset(ORF,ORF$category == "Unstable but no effect")

a<-as.data.frame(table(stabilized_ORF$Gene_Qualifier))
b<-as.data.frame(table(Intermediate_ORF$Gene_Qualifier))
c<-as.data.frame(table(NoEffect_ORF$Gene_Qualifier))
d<-as.data.frame(table(UnstableButNoEffect_ORF$Gene_Qualifier))
q<-as.data.frame(cbind(a,b,c,d))

q<-q[,c(1,2,4,6,8)]
names(q)<-c("Var1","Stabilized","Intermediate","NoEffect","UnstableButNoEffect")
row.names(q)<-q$Var1
q$Var1<-NULL

q$freq_stabil<-as.numeric(q$Stabilized)/199
q$freq_inter<-as.numeric(q$Intermediate)/49
q$freq_UnNoeffect<-as.numeric(q$UnstableButNoEffect)/219
q$fre_noeffecct<-as.numeric(q$NoEffect)/6119
q[3,9:12]<-q[3,5:8]/(5244/6586)
q[2,9:12]<-q[2,5:8]/(699/6586)
q[1,9:12]<-q[1,5:8]/(673/6586)

names(q)<-c("Stabilized","Intemediate","NoEffect","UnstableButNoEffect","freq_Stabilized","freq_Intermediate","freq_unstableNoEffect","freq_NoEffect","Norm_Stabilized","Norm_Intermediate","Norm_Noeffect","Norm_Unstable_NoEffect")
q$Stabilized<-as.numeric(q$Stabilized)
q$Intemediate<-as.numeric(q$Intemediate)
q$NoEffect<-as.numeric(q$NoEffect)
q$UnstableButNoEffect<-as.numeric(q$UnstableButNoEffect)
f<-chisq.test(as.data.frame(t(q[,1:4])))
q$ORF<-row.names(q)
q_melt<-melt(q[,9:13], id = "ORF")

ggplot(q_melt, aes(fill=variable , y=value, x=ORF)) + 
  geom_bar(position="dodge", stat="identity")+
  theme_bw()+
  ylab("ORF Type")+
  xlab("Normalized Value")+
  scale_fill_manual(values=c("#2E2D4D","#4290C0","#677DB7","#999999"))+
  ggtitle("Enrichment of ORF type")



# go 


yeastmine2<-read.csv(paste0(path_input,"yeastmine_3.csv"), sep = ",", header = FALSE)

names(yeastmine2)<-c("GeneID","Gene_Systematic_Name","Organism","Standard_Name","Gene_Name","Gene_Qualifier","Ontology_Term_ID","Ontology_Term_Name","Ontology_Term_Namespace","Protein_Standard_Name","Residue","Protein_Systematic_Name")
yeastmine2$Residue<-substr(yeastmine2$Residue,1,(nchar(yeastmine2$Residue)-1))
yeastmine2$translation<-substr(yeastmine2$Residue,(nchar(yeastmine2$Residue)-11),(nchar(yeastmine2$Residue)))

yeastmine2<-yeastmine2[,c("GeneID","Protein_Systematic_Name","Gene_Qualifier","translation","Ontology_Term_ID","Ontology_Term_Name","Ontology_Term_Namespace")]
yeastmine2<-unique(yeastmine2)

GO<-merge(WT_Das1,yeastmine2, id = "translation")
GO$gotermid<-paste0(GO$Ontology_Term_ID, "_", GO$Ontology_Term_Name,"_",GO$Ontology_Term_Namespace)
stabilized_for_go<-subset(GO,GO$category == "Stabilized")
stabilized_for_go<-stabilized_for_go[,c("translation","gotermid","Protein_Systematic_Name")]
stabilized_for_go<-unique(stabilized_for_go)

#stabilized_for_go$gotermid<-paste0(stabilized_for_go$GOId, "_", stabilized_for_go$GOtTerm)
a<-as.data.frame(table(stabilized_for_go$gotermid))
b<-as.data.frame(table(GO$gotermid))
names(a)<-c("GO","Frequency_stabilized")
names(b)<-c("GO","Frequency_all")
stable_merged_GO<-merge(a,b, id = "Var1", all.x = TRUE)

stable_merged_GO$norm_freq<-(stable_merged_GO$Frequency_stabilized/199)/(stable_merged_GO$Frequency_all/6580)


stable_merged_GO<-subset(stable_merged_GO,stable_merged_GO$Frequency_stabilized > 5)
stable_merged_GO$GOID<-str_split_fixed(stable_merged_GO$GO,"_",3)[,1]
stable_merged_GO$GOName<-str_split_fixed(stable_merged_GO$GO,"_",3)[,2]
stable_merged_GO$GONapspace<-str_split_fixed(stable_merged_GO$GO,"_",3)[,3]
stable_merged_GO$GONapspace<-ifelse(stable_merged_GO$GONapspace == "process_biological_process","biological_process",stable_merged_GO$GONapspace)
stable_merged_GO$GONapspace<-ifelse(stable_merged_GO$GONapspace == "function_molecular_function","molecular_function",stable_merged_GO$GONapspace)
stable_merged_GO$GONapspace<-ifelse(stable_merged_GO$GONapspace == "component_cellular_component","cellular_component",stable_merged_GO$GONapspace)

stable_merged_GO$status<-ifelse(stable_merged_GO$norm_freq< 1, "Not Enriched","Enriched")
stable_merged_GO<-stable_merged_GO[order(stable_merged_GO$GONapspace),]


pdf("GO_term_with_class.pdf")

ggplot(data=stable_merged_GO, aes(x=norm_freq , y=GOName,fill=status)) +
  geom_bar(stat="identity" ,position="dodge")+
  theme_bw()+
  geom_text(aes(label=paste0(Frequency_stabilized," / ",Frequency_all)), vjust=0.3, hjust = 0,color="black", size=2)+
  theme(text = element_text(size = 8))+
  #theme(axis.text.y=element_text(aes(color=GONapspace)))+
  facet_wrap(~GONapspace)+
  xlab("Normalized Frequency")+
  ylab("Gene Ontology Name")+
  xlim(0,15)+
  ggtitle("GO enrichment for stabilized group")+
  scale_fill_manual(values = c( "Enriched" = "#ff5400",
                                "Not Enriched" = "#8d99ae"))

dev.off()



#### marking

#marking peptides:
d<-WT_Das1[,c("translation","value.x","value.y")]
names(d)<-c("translation","PSI_WT","PSI_Das1")
peptides<-c("VRGAVGGWRLVG","ARWRVGLWRGAS","PSLAASFWRVVG",
            "CTSCGYKFRTNN","LKILRQKMNHQN"
)

d$Label1<-ifelse(d$translation %in% peptides, "tested degron","Others")
d$Label2<-"Others"
d$Label2<-ifelse(d$translation == "CTSCGYKFRTNN", "Rpa12c",d$Label2)
d$Label2<-ifelse(d$translation == "LKILRQKMNHQN", "Atg1c",d$Label2)



d_sub<-subset(d,!(d$Label2 == "Others"))

pdf(paste0(path_plot,"Figure6/highlighted_ATgc_RPL.pdf"))
print(
  ggplot(d)+
    geom_point(data = subset(d,d$Label1 == "Others"), 
               aes(x = PSI_WT, y = PSI_Das1), 
               color = "grey", 
               alpha = 0.2)+
    geom_point(data = d_sub, 
               aes(x = PSI_WT, y = PSI_Das1), 
               color = "red", 
               alpha = 0.6)+
    geom_text_repel(data=d_sub,aes(x=PSI_WT,y=PSI_Das1,label=Label2), color='red', size = 3.5)+
    theme_bw()+
    #scale_color_manual(values=c("#999999", "#4290C0"))+
    #geom_smooth(method=lm, se=FALSE, linetype = "dashed", color = "black")+
    xlab("PSI(WT)")+
    ylab(" PSI(Das1)")+
    xlim(0,1)+
    ylim(0,1)+
    
    ggtitle("PSI WT versus Das1 - Cter"))
dev.off()

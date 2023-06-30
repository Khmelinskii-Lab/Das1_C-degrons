# set the directory


##### for project : Analysis of degron at C terminal
##### sub text:merging with other datasets - do10,das1,scf
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")

#### -------------------------------------- for correlation among different dataset --------------------------------------------------------------###
merged_file<-read.csv(paste0(path_output,"/merged_file_for_rep2and3.csv"))
deg<-read.csv("Y:/lab data/susmitha/edwin/March23/output_files/degron_effect_d1_WT.csv")
merged_file_non_stop<-subset(merged_file,merged_file$Codon == "Non Stop Codon")

merged_file_subset<-subset(merged_file_non_stop, merged_file_non_stop$raw_counts_translation %in% deg$raw_counts_translation)
merged_file_subset$sub_experiment<-paste0(merged_file_subset$replicate,"_",merged_file_subset$raw_counts_sub_experiment)


a<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_WT_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
b<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_d1_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
c<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_d10_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
d<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_WT_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
e<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_d1_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
f<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_c34_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
g<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_c53_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]

h<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_WT_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
i<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_d1_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
j<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_d10_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
k<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_WT_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
l<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_d1_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
m<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_c34_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]
n<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_c53_37C_un")[,c("raw_counts_translation","byaa_pooled_PSI")]


a<-unique(a)
b<-unique(b)
c<-unique(c)
d<-unique(d)
e<-unique(e)
f<-unique(f)
g<-unique(g)
h<-unique(h)
i<-unique(i)
j<-unique(j)
k<-unique(k)
l<-unique(l)
m<-unique(m)
n<-unique(n)

PSI_merged<-merge(a,b, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,c, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,d, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,e, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,f, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,g, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,h, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,i, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,j, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,k, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,l, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,m, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged,n, by = "raw_counts_translation")
PSI_merged<-merge(PSI_merged, deg, by = "raw_counts_translation")

names(PSI_merged)<-c("Translations",
                     "Unstable_WT_rep3",
                     "Unstable_Das1_rep3",
                     "Unstable_Doa10_rep3",
                     "Unstable_WT_SCF_rep3",
                     "Unstable_Das1_SCF_rep3",
                     "Unstable_cdc34_SCF_rep3",
                     "Unstable_cdc53_SCF_rep3",
                     "Unstable_WT_rep2",
                     "Unstable_Das1_rep2",
                     "Unstable_Doa10_rep2",
                     "Unstable_WT_SCF_rep2",
                     "Unstable_Das1_SCF_rep2",
                     "Unstable_cdc34_SCF_rep2",
                     "Unstable_cdc53_SCF_rep2",
                     "X",
                     "Unstable_Das1_rep1",
                     "Unstable_WT_rep1",
                     "Effect_rep1")

PSI_merged$X<-NULL
PSI_merged$Effect_rep1<-NULL
PSI_merged<-PSI_merged[!duplicated(PSI_merged$Translations), ]
row.names(PSI_merged)<-PSI_merged$Translations
PSI_merged$Translations<-NULL

dir.create(paste0(path_plot,"Figure6/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

dir.create(paste0(path_plot,"FigureS5/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)
dir.create(paste0(path_plot,"FigureS6/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

PSI_merged<-read.csv(
          paste0(path_output,
            "/PSI_file_for_other_dataset.csv"))
row.names(PSI_merged)<-PSI_merged$X
PSI_merged$X<-NULL
PSI_merged<-PSI_merged[complete.cases(PSI_merged),]
PSI_merged$type<-NULL

M = cor(as.matrix(PSI_merged))

M<-as.data.frame(M)
M$dataset<-row.names(M)
M_melt<-melt(M, id = "dataset")
M_melt$variable<-factor(M_melt$variable,
                        levels = (c("Unstable_WT_rep1",
                                    "Unstable_WT_rep2",
                                    "Unstable_WT_rep3",
                                    "Unstable_WT_SCF_rep3",
                                    "Unstable_WT_SCF_rep2",
                                    "Unstable_Doa10_rep2",
                                    "Unstable_Doa10_rep3",
                                    "Unstable_Das1_rep1",
                                    "Unstable_Das1_rep2",
                                    "Unstable_Das1_rep3",
                                    "Unstable_Das1_SCF_rep2",
                                    "Unstable_Das1_SCF_rep3" ,
                                    "Unstable_cdc34_SCF_rep2",
                                    "Unstable_cdc34_SCF_rep3",
                                    "Unstable_cdc53_SCF_rep2",
                                    "Unstable_cdc53_SCF_rep3"
                        )))
M_melt$dataset<-factor(M_melt$dataset,
                       levels = (c("Unstable_WT_rep1",
                                   "Unstable_WT_rep2",
                                   "Unstable_WT_rep3",
                                   "Unstable_WT_SCF_rep3",
                                   "Unstable_WT_SCF_rep2",
                                   "Unstable_Doa10_rep2",
                                   "Unstable_Doa10_rep3",
                                   "Unstable_Das1_rep1",
                                   "Unstable_Das1_rep2",
                                   "Unstable_Das1_rep3",
                                   "Unstable_Das1_SCF_rep2",
                                   "Unstable_Das1_SCF_rep3" ,
                                   "Unstable_cdc34_SCF_rep2",
                                   "Unstable_cdc34_SCF_rep3",
                                   "Unstable_cdc53_SCF_rep2",
                                   "Unstable_cdc53_SCF_rep3"
                       )))

pdf(paste0(path_plot,"FigureS5/correlation_among_all_dataset.pdf"))

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





#### -------------------------------------- for correlation among different dataset --------------------------------------------------------------###

q<-as.data.frame(PSI_merged[,c("Unstable_WT_rep1","Unstable_WT_rep2","Unstable_WT_rep3")])
WT<-lmFit(q[,c(1,2,3)])
mean_result_WT<-q
mean_result_WT$effect<-WT$sigma * WT$stdev.unscaled
mean_result_WT$ci_95<-WT$coefficients * qt(0.975 , WT$df.residual)
mean_result_WT$value<-WT$coefficients

p<-as.data.frame(PSI_merged[,c("Unstable_Das1_rep1","Unstable_Das1_rep2","Unstable_Das1_rep3")])
DT<-lmFit(p[,c(1,2,3)])
mean_result_DT<-p
mean_result_DT$effect<-DT$sigma * DT$stdev.unscaled
mean_result_DT$ci_95<-DT$coefficients * qt(0.975 , DT$df.residual)
mean_result_DT$value<-DT$coefficients


mean_result_DT$translation<-row.names(mean_result_DT)
mean_result_WT$translation<-row.names(mean_result_WT)
WT_DT<-merge(mean_result_WT, mean_result_DT, by= "translation")

# variation of Das1 with WT
ggplot(WT_DT)+
  geom_point(aes(x = value.x, y = value.y))+
  theme_bw()+
  xlab("Wild type ")+
  ylab("Das1")+
  xlim(0,1)+
  ylim(0,1)+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(WT_DT$value.x,WT_DT$value.y),2)),
           colour = "red")+
  ggtitle("PSI Das1 versus WT after lmFit")




### highlighting the control peptides
only_instable<-readxl::read_xlsx(paste0(path_input,"/other_dataset/4726.xlsx"), sheet = 1)
WT_DT$group<-"Stable"
WT_DT$group<-ifelse(WT_DT$translation %in% only_instable$translation, "Unstable", WT_DT$group)


pdf(paste0(path_plot,"FigureS5/Stable_highlights_das1.pdf"))
WT_DT %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = group))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stable" = "#ff5400",
    #"Intermediate" = "#8d99ae",
    "Unstable" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (DAS1)")+
  ggtitle("Highlighting stable peptides in lmFit data")
dev.off()

### remove dataset with stable controls in the plots
WT_DT_unstable<-subset(WT_DT,WT_DT$group == "Unstable")


#### -------------------------------------- different category --------------------------------------------------------------###

q_data<-WT_DT_unstable[,c("translation" ,"Unstable_WT_rep1","Unstable_WT_rep2","Unstable_WT_rep3","Unstable_Das1_rep1","Unstable_Das1_rep2","Unstable_Das1_rep3")]
q_data<-melt(q_data,id = "translation")
q_data$group<-str_split_fixed(q_data$variable,"_",3)[,2]
stat.test <- q_data %>%
  group_by(translation) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

WT_DT_unstable$differnce<-WT_DT_unstable$value.y - WT_DT_unstable$value.x

WT_DT_unstable$color<-ifelse(WT_DT_unstable$differnce > 0.1, "Potential degron","No effect")
WT_DT_unstable$category<-"No Effect"
WT_DT_unstable$category<-ifelse((WT_DT_unstable$color == "Potential degron" & 
                          WT_DT_unstable$value.y < 0.58),"Intermediate", WT_DT_unstable$category)
WT_DT_unstable$category<-ifelse((WT_DT_unstable$color == "Potential degron" & WT_DT_unstable$value.y > 0.58 ),"Stabilized", WT_DT_unstable$category)
#WT_DT_unstable$category<-ifelse((WT_DT_unstable$color == "Potential degron" & WT_DT_unstable$value.y > 0.58 ),"Strongly Stabilized", WT_DT_unstable$category)
WT_DT_unstable$category<-ifelse((WT_DT_unstable$color == "Potential degron" & WT_DT_unstable$value.x > 0.5 ),"No Effect", WT_DT_unstable$category)



data_with_stat_test<-merge(stat.test, WT_DT_unstable, id  = "translation")
data_with_stat_test$category<-ifelse(data_with_stat_test$p.adj > 0.05,"No Effect" ,data_with_stat_test$category)

pdf(paste0(path_plot,"Figure5/3Category_WT_versus_das1.pdf"))
data_with_stat_test %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = category))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
                                 "Stabilized" = "#ff5400",
                                 "Intermediate" = "#8d99ae",
                                 "No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (DAS1)")+
  ggtitle("Putting statistically insignificant points as No Effect  with multiple testing")+
  annotate("text", 
           x = 0.15, y = 0.92, 
           label = paste0("Threshold : \n Difference :0.1, \nPSI WT: 0.5; \nPSI Das1 : 0.58"),
           colour = "red", fontface =2)
dev.off()
# hydrophobicity

data_with_stat_test$hydrophobicity<-hydrophobicity(data_with_stat_test$translation)
#data_with_stat_test1<-data_with_stat_test[,c("category","hydrophobicity")]


my_comparisons <- list( c("No Effect", "Stabilized"), c("Stabilized", "Intermediate"), c("No Effect", "Intermediate") )
pdf(paste0(path_plot,"Figure5/3Category_hydrophobicity.pdf"))
ggplot(data_with_stat_test)+
  geom_boxplot(
    aes(x = category, 
        y = hydrophobicity, color = category),
    coef = 1e30,lwd=0.25
  )+
  ylim(-4,4)+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stabilized" = "#ff5400",
    "Intermediate" = "#8d99ae",
    "No Effect" = "#2b2d42"))+
  guides(col = FALSE)+
  ylab("Hydrophobicity")+
  xlab("")+
  ggtitle("Variation of hydrophobicity")+
  theme_bw()+
  stat_compare_means(data = data_with_stat_test, 
                     aes(x = category, 
                         y = hydrophobicity),
                     method = "wilcox.test", 
                     ref.group = "No Effect",
                     label.y = 4)+
  stat_compare_means(data = data_with_stat_test, 
                     aes(x = category, 
                         y = hydrophobicity),
                     label = "p.signif", method = "wilcox.test",
                     ref.group = "No Effect",
                     label.y = 3.5)
  



dev.off()

# frequency distribution
dataset<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))
dataset<-subset(dataset,dataset$Codon == "Non Stop Codon")
names(dataset)
dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]
dataset_for_frequency<-unique(dataset_for_frequency)
dataset_for_frequency<-dataset_for_frequency[!duplicated(dataset_for_frequency$raw_counts_translation),]
dataset_for_frequency_oneHot<-create_one_hot_db(dataset_for_frequency$raw_counts_translation,12)

WT_DT_unstable_data<-create_one_hot_db(WT_DT_unstable$translation,12)
amino_acid<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
t<-strrep(amino_acid, 12)  


stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
count<-nrow(stabilized)
translation<-c(t,stabilized$translation)
stabilized<-create_one_hot_db(translation,12)
stabilized<-stabilized[21:nrow(stabilized),]
stabilized<-as.data.frame(colSums(stabilized))

stabilized$frequency<-stabilized$`colSums(stabilized)`/count
stabilized$allFreq<-as.data.frame(colSums(WT_DT_unstable_data))[,1]/4709
stabilized$log2NormFreq<-log2(stabilized$frequency/stabilized$allFreq)
stabilized$NormFreq<-(stabilized$frequency/stabilized$allFreq)
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

pdf(paste0(path_plot,"Figure5/Stabilized_frequency.pdf"))
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
      
      scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_Unstable)",
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


Intermediate<-subset(data_with_stat_test,data_with_stat_test$category == "Intermediate")
count<-nrow(Intermediate)
translation<-c(t,Intermediate$translation)
Intermediate<-create_one_hot_db(translation,12)
Intermediate<-Intermediate[21:nrow(Intermediate),]
Intermediate<-as.data.frame(colSums(Intermediate))

Intermediate$frequency<-Intermediate$`colSums(Intermediate)`/count
Intermediate$allFreq<-as.data.frame(colSums(WT_DT_unstable_data))[,1]/4709
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
pdf(paste0(path_plot,"Figure5/Intermediate_frequency.pdf"))
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
      scale_fill_gradientn("log2(Freq_Intermediate)/\n (Freq_Unstable)",
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




write.csv(data_with_stat_test,
          paste0(path_output,"data_with_stat_test_wt_das1.csv")
          )

############ doa10 -  wt

q<-as.data.frame(PSI_merged[,c("Unstable_WT_rep1","Unstable_WT_rep2","Unstable_WT_rep3")])
WT<-lmFit(q[,c(1,2,3)])
mean_result_WT<-q
mean_result_WT$effect<-WT$sigma * WT$stdev.unscaled
mean_result_WT$ci_95<-WT$coefficients * qt(0.975 , WT$df.residual)
mean_result_WT$value<-WT$coefficients

p<-as.data.frame(PSI_merged[,c("Unstable_Doa10_rep2","Unstable_Doa10_rep3")])
Doa10<-lmFit(p[,c(1,2)])
mean_result_Doa10<-p
mean_result_Doa10$effect<-Doa10$sigma * Doa10$stdev.unscaled
mean_result_Doa10$ci_95<-Doa10$coefficients * qt(0.975 , Doa10$df.residual)
mean_result_Doa10$value<-Doa10$coefficients


mean_result_Doa10$translation<-row.names(mean_result_Doa10)
mean_result_WT$translation<-row.names(mean_result_WT)
WT_Doa10<-merge(mean_result_WT, mean_result_Doa10, by= "translation")

ggplot(WT_Doa10)+
  geom_point(aes(x = value.x, y = value.y))+
  theme_bw()+
  xlab("Wild type ")+
  ylab("Doa10")+
  xlim(0,1)+
  ylim(0,1)+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(WT_Doa10$value.x,WT_Doa10$value.y),2)),
           colour = "red")+
  ggtitle("PSI Doa10 versus WT after lmFit")

#only_instable<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/4726.xlsx")
WT_Doa10$group<-"Stable"
WT_Doa10$group<-ifelse(WT_Doa10$translation %in% only_instable$translation, "Unstable", WT_Doa10$group)

pdf(paste0(path_plot,"FigureS5/Stable_highlights_doa10.pdf"))
WT_Doa10 %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = group))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stable" = "#ff5400",
    #"Intermediate" = "#8d99ae",
    "Unstable" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (Doa10)")+
  ggtitle("Highlighting stable peptides in lmFit data")
dev.off()

### remove dataset with stable controls in the plots
WT_Doa10_unstable<-subset(WT_Doa10,WT_Doa10$group == "Unstable")


#### -------------------------------------- different category --------------------------------------------------------------###

q_data<-WT_Doa10_unstable[,c("translation" ,"Unstable_WT_rep1","Unstable_WT_rep2","Unstable_WT_rep3","Unstable_Doa10_rep2","Unstable_Doa10_rep3")]
q_data<-melt(q_data,id = "translation")
q_data$group<-str_split_fixed(q_data$variable,"_",3)[,2]
stat.test <- q_data %>%
  group_by(translation) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

WT_Doa10_unstable$differnce<-WT_Doa10_unstable$value.y - WT_Doa10_unstable$value.x

WT_Doa10_unstable$color<-ifelse(WT_Doa10_unstable$differnce > 0.1, "Potential degron","No effect")
WT_Doa10_unstable$category<-"No Effect"
WT_Doa10_unstable$category<-ifelse((WT_Doa10_unstable$color == "Potential degron" & 
                                      WT_Doa10_unstable$value.y < 0.58),"Intermediate", WT_Doa10_unstable$category)
WT_Doa10_unstable$category<-ifelse((WT_Doa10_unstable$color == "Potential degron" & WT_Doa10_unstable$value.y > 0.58 ),"Stabilized", WT_Doa10_unstable$category)
#WT_Doa10_unstable$category<-ifelse((WT_Doa10_unstable$color == "Potential degron" & WT_Doa10_unstable$value.y > 0.58 ),"Strongly Stabilized", WT_Doa10_unstable$category)
WT_Doa10_unstable$category<-ifelse((WT_Doa10_unstable$color == "Potential degron" & WT_Doa10_unstable$value.x > 0.5 ),"No Effect", WT_Doa10_unstable$category)



data_with_stat_test<-merge(stat.test, WT_Doa10_unstable, id  = "translation")
data_with_stat_test$category<-ifelse(data_with_stat_test$p > 0.05,"No Effect" ,data_with_stat_test$category)
pdf(paste0(path_plot,"Figure5/WT_Doa10_3Category.pdf"))
data_with_stat_test %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = category))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stabilized" = "#ff5400",
    "Intermediate" = "#8d99ae",
    "No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (Doa10)")+
  ggtitle("Putting statistically insignificant points as No Effect  without multiple testing")+
  annotate("text", 
           x = 0.15, y = 0.92, 
           label = paste0("Threshold : \n Difference :0.1, \nPSI WT: 0.5; \nPSI Doa10 : 0.58"),
           colour = "red", fontface =2)
dev.off()
# hydrophobicity


data_with_stat_test$hydrophobicity<-hydrophobicity(data_with_stat_test$translation)

my_comparisons <- list( c("No Effect", "Stabilized"), c("Stabilized", "Intermediate"), c("No Effect", "Intermediate") )
pdf(paste0(path_plot,"Figure5/3Category_hydrophobicity_DOA10Wt.pdf"))
ggplot(data_with_stat_test)+
  geom_boxplot(
    aes(x = category, 
        y = hydrophobicity, 
        color = category),
    coef = 1e30,
    lwd=0.25
  )+
  ylim(-4,4)+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stabilized" = "#ff5400",
    "Intermediate" = "#8d99ae",
    "No Effect" = "#2b2d42"))+
  guides(col = FALSE)+
  ylab("Hydrophobicity")+
  xlab("")+
  ggtitle("Variation of hydrophobicity")+
  theme_bw()+
  stat_compare_means(data = data_with_stat_test, 
                     aes(x = category, 
                         y = hydrophobicity),
                     method = "wilcox.test", 
                     ref.group = "No Effect",
                     label.y = 4)+
  stat_compare_means(data = data_with_stat_test, 
                     aes(x = category, 
                         y = hydrophobicity),
                     label = "p.signif", method = "wilcox.test",
                     ref.group = "No Effect",
                     label.y = 3.5)
dev.off()

data_with_stat_test_Das1<-read.csv(paste0(path_output,"data_with_stat_test_wt_das1.csv"))

Das1_group<-data_with_stat_test_Das1[,c("translation","category")]
names(Das1_group)<-c("translation","category_Das1")
WT_Doa10_unstable<-merge(WT_Doa10_unstable,Das1_group, by = "translation")

# seeing if there are any doa10 substrates that are also das1 substrates
ggplot(WT_Doa10_unstable)+
  
  geom_point(data = WT_Doa10_unstable,aes(x = value.x, y = value.y), color = "black")+
  geom_point(data = subset(WT_Doa10_unstable,WT_Doa10_unstable$category_Das1 == "Intermediate"),
             aes(x = value.x, y = value.y), color = "grey")+
  geom_point(data = subset(WT_Doa10_unstable,WT_Doa10_unstable$category_Das1 == "Stabilized"),
             aes(x = value.x, y = value.y), color = "red", alpha = 0.4)+
  
  theme_bw()+
  #scale_color_manual(values = c( 
  # #"Strongly Stabilized" = "#98473E",
  #"Stabilized" = "#ff5400",
  #"Intermediate" = "#8d99ae",
  #"No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (Doa10)")+
  ggtitle("Marking Das1 substrates in Doa10-WT variation")

# frequency 


WT_Doa10_unstable_data<-create_one_hot_db(WT_Doa10_unstable$translation,12)
amino_acid<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
t<-strrep(amino_acid, 12)  

WT_Doa10_unstable$differnce<-WT_Doa10_unstable$value.y - WT_Doa10_unstable$value.x
stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
count<-nrow(stabilized)
translation<-c(t,stabilized$translation)
stabilized<-create_one_hot_db(translation,12)
stabilized<-stabilized[21:nrow(stabilized),]
stabilized<-as.data.frame(colSums(stabilized))

stabilized$frequency<-stabilized$`colSums(stabilized)`/count
stabilized$allFreq<-as.data.frame(colSums(WT_Doa10_unstable_data))[,1]/4709
stabilized$log2NormFreq<-log2(stabilized$frequency/stabilized$allFreq)
stabilized$NormFreq<-(stabilized$frequency/stabilized$allFreq)
stabilized$identifier<-row.names(stabilized)
stabilized$position<-substr(stabilized$identifier,2,(nchar(stabilized$identifier)-1))
stabilized$position<-as.numeric(stabilized$position)-13
stabilized$amino_acid<-substr(stabilized$identifier,(nchar(stabilized$identifier)),(nchar(stabilized$identifier)))
stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
stabilized$amino_acid<-as.character(stabilized$amino_acid)
stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"FigureS5/Stabilized_Frequency_Doa10.pdf"))
print(
  plot_grid(
    
    ggplot(stabilized, aes (x = position,
                            y= amino_acid,
                            fill = log2NormFreq))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreq < 2.5 & 
                                stabilized$log2NormFreq > -2.5), 
                aes(
                  x = position,
                  y= amino_acid,
                  fill = log2NormFreq
                ))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreq < -2.5 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreq > 2.5 ), 
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
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
      scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_Unstable)",
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
                           values=scales::rescale(c(-2.5 ,-2.25 ,-2.0, -1.5 ,-1.0 ,-0.5 , 0.0 , 0.5 , 1.0  ,1.5 , 2.0 , 2.25 , 2.5)),
                           limits=c(-2.5,2.5)),
    
    
    
    
    ggplot(stabilized, aes (x = position,
                            y= amino_acid,
                            fill = log2NormFreqAll))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreqAll < 2.5 & 
                                stabilized$log2NormFreqAll > -2.5), 
                aes(
                  x = position,
                  y= amino_acid,
                  fill = log2NormFreqAll
                ))+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreqAll < -2.5 ), 
                aes(
                  x = position,
                  y= amino_acid
                ),
                fill = "#F7931E")+
      geom_tile(data = subset(stabilized,stabilized$log2NormFreqAll > 2.5 ), 
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
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
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
                           values=scales::rescale(c(-2.50 ,-2.25, -2.0, -1.5 ,-1.0 ,-0.5 , 0.0,  0.5,  1.0  ,1.5 , 2.0 , 2.25 , 2.5)),
                           limits=c(-2.5,2.5)),
    
    
    ncol = 1
  )
)
dev.off()


Intermediate<-subset(data_with_stat_test,data_with_stat_test$category == "Intermediate")
count<-nrow(Intermediate)
translation<-c(t,Intermediate$translation)
Intermediate<-create_one_hot_db(translation,12)
Intermediate<-Intermediate[21:nrow(Intermediate),]
Intermediate<-as.data.frame(colSums(Intermediate))

Intermediate$frequency<-Intermediate$`colSums(Intermediate)`/count
Intermediate$allFreq<-as.data.frame(colSums(WT_Doa10_unstable_data))[,1]/4709
Intermediate$log2NormFreq<-log2(Intermediate$frequency/Intermediate$allFreq)
Intermediate$NormFreq<-(Intermediate$frequency/Intermediate$allFreq)
Intermediate$identifier<-row.names(Intermediate)
Intermediate$position<-substr(Intermediate$identifier,2,(nchar(Intermediate$identifier)-1))
Intermediate$position<-as.numeric(Intermediate$position)-13
Intermediate$amino_acid<-substr(Intermediate$identifier,(nchar(Intermediate$identifier)),(nchar(Intermediate$identifier)))
Intermediate$position<-factor(Intermediate$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
Intermediate$amino_acid<-as.character(Intermediate$amino_acid)
Intermediate[Intermediate$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


Intermediate$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
Intermediate$log2NormFreqAll<-log2(Intermediate$frequency/Intermediate$freqAllData)
Intermediate$NormFreqAll<-(Intermediate$frequency/Intermediate$freqAllData)

Intermediate[Intermediate$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"FigureS5/Intermediate_Frequency_Doa10.pdf"))
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
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
      scale_fill_gradientn("log2(Freq_Intermediate)/\n (Freq_Unstable)",
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
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
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

pdf(paste0(path_plot,"FigureS5/Intermediate_sequence_logo_doa10.pdf"))
print(
  ggseqlogo(subset(data_with_stat_test,data_with_stat_test$category == "Intermediate")[,c("translation")])
)
dev.off()
pdf(paste0(path_plot,"FigureS5/Stabilized_sequence_logo_doa10.pdf"))
print(
  ggseqlogo(subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")[,c("translation")])
)
dev.off()

write.csv(data_with_stat_test,
          paste0(path_output,"data_with_stat_test_wt_doa10.csv")
)

write.csv(WT_Doa10_unstable,
          paste0(path_output,"WT_Doa10_unstable.csv")
)

### wt -cdc 34
q<-as.data.frame(PSI_merged[,c("Unstable_WT_SCF_rep2","Unstable_WT_SCF_rep3")])
WT<-lmFit(q[,c(1,2)])
mean_result_WT<-q
mean_result_WT$effect<-WT$sigma * WT$stdev.unscaled
mean_result_WT$ci_95<-WT$coefficients * qt(0.975 , WT$df.residual)
mean_result_WT$value<-WT$coefficients

p<-as.data.frame(PSI_merged[,c("Unstable_cdc34_SCF_rep2","Unstable_cdc34_SCF_rep3")])
cdc34<-lmFit(p[,c(1,2)])
mean_result_cdc34<-p
mean_result_cdc34$effect<-cdc34$sigma * cdc34$stdev.unscaled
mean_result_cdc34$ci_95<-cdc34$coefficients * qt(0.975 , cdc34$df.residual)
mean_result_cdc34$value<-cdc34$coefficients


mean_result_cdc34$translation<-row.names(mean_result_cdc34)
mean_result_WT$translation<-row.names(mean_result_WT)
WT_cdc34<-merge(mean_result_WT, mean_result_cdc34, by= "translation")

#pdf(paste0(path_plot,"PSI_cdc34_WT_afterLMFit.pdf"))
ggplot(WT_cdc34)+
  geom_point(aes(x = value.x, y = value.y))+
  theme_bw()+
  xlab("Wild type ")+
  ylab("cdc34")+
  xlim(0,1)+
  ylim(0,1)+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(WT_cdc34$value.x,WT_cdc34$value.y),2)),
           colour = "red")+
  ggtitle("PSI cdc34 versus WT after lmFit")
#dev.off()

#only_instable<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/4726.xlsx")
WT_cdc34$group<-"Stable"
WT_cdc34$group<-ifelse(WT_cdc34$translation %in% only_instable$translation, "Unstable", WT_cdc34$group)
#pdf(paste0(path_plot,"WT_cdc34_stable_highlighted.pdf"))
WT_cdc34 %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = group))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stable" = "#8d99ae",
    #"Intermediate" = "#8d99ae",
    "Unstable" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (cdc34)")+
  ggtitle("Highlighting stable peptides in lmFit data")
#dev.off()

### remove dataset with stable controls in the plots
WT_cdc34_unstable<-subset(WT_cdc34,WT_cdc34$group == "Unstable")


#### -------------------------------------- different category --------------------------------------------------------------###

q_data<-WT_cdc34_unstable[,c("translation" ,"Unstable_WT_SCF_rep2","Unstable_WT_SCF_rep3","Unstable_cdc34_SCF_rep2","Unstable_cdc34_SCF_rep3")]
q_data<-melt(q_data,id = "translation")
q_data$group<-str_split_fixed(q_data$variable,"_",3)[,2]
stat.test <- q_data %>%
  group_by(translation) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

WT_cdc34_unstable$differnce<-WT_cdc34_unstable$value.y - WT_cdc34_unstable$value.x

WT_cdc34_unstable$color<-ifelse(WT_cdc34_unstable$differnce > 0.1, "Potential degron","No effect")
WT_cdc34_unstable$category<-"No Effect"
WT_cdc34_unstable$category<-ifelse((WT_cdc34_unstable$color == "Potential degron" & 
                                      WT_cdc34_unstable$value.y < 0.67),"Intermediate", WT_cdc34_unstable$category)
WT_cdc34_unstable$category<-ifelse((WT_cdc34_unstable$color == "Potential degron" & WT_cdc34_unstable$value.y > 0.67 ),"Stabilized", WT_cdc34_unstable$category)
#WT_cdc34_unstable$category<-ifelse((WT_cdc34_unstable$color == "Potential degron" & WT_cdc34_unstable$value.y > 0.58 ),"Strongly Stabilized", WT_cdc34_unstable$category)
WT_cdc34_unstable$category<-ifelse((WT_cdc34_unstable$color == "Potential degron" & WT_cdc34_unstable$value.x > 0.5 ),"No Effect", WT_cdc34_unstable$category)



data_with_stat_test<-merge(stat.test, WT_cdc34_unstable, id  = "translation")
data_with_stat_test$category<-ifelse(data_with_stat_test$p.adj < 0.05,"No Effect" ,data_with_stat_test$category)
# hydrophobicity

data_with_stat_test$hydrophobicity<-hydrophobicity(data_with_stat_test$translation)



data_with_stat_test_Das1<-read.csv(paste0(path_output,"data_with_stat_test_wt_das1.csv"))

Das1_group<-data_with_stat_test_Das1[,c("translation","category")]
names(Das1_group)<-c("translation","category_Das1")
WT_cdc34_unstable<-merge(WT_cdc34_unstable,Das1_group, by = "translation")

pdf(paste0(path_plot,"FigureS5/Das1Peptides_highlighted_in_cdc34WT.pdf"))
ggplot(WT_cdc34_unstable)+
  
  geom_point(data = subset(WT_cdc34_unstable,WT_cdc34_unstable$category_Das1 == "No Effect"),
             aes(x = value.x, y = value.y), color = "black")+
  geom_point(data = subset(WT_cdc34_unstable,WT_cdc34_unstable$category_Das1 == "Intermediate"),
             aes(x = value.x, y = value.y), color = "grey")+
  geom_point(data = subset(WT_cdc34_unstable,WT_cdc34_unstable$category_Das1 == "Stabilized"),
             aes(x = value.x, y = value.y), color = "red")+
  
  theme_bw()+
  #scale_color_manual(values = c( 
  # #"Strongly Stabilized" = "#98473E",
  #"Stabilized" = "#ff5400",
  #"Intermediate" = "#8d99ae",
  #"No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (cdc34)")+
  ggtitle("Marking Das1 substrates peptides in cdc34-WT variation")
dev.off()

write.csv(data_with_stat_test,
          paste0(path_output,"data_with_stat_test_wt_cdc34.csv")
)
write.csv(WT_cdc34_unstable,
          paste0(path_output,"WT_cdc34_unstable.csv")
)

# cdc 34 das1
q<-as.data.frame(PSI_merged[,c("Unstable_Das1_SCF_rep2","Unstable_Das1_SCF_rep3")])
Das1<-lmFit(q[,c(1,2)])
mean_result_Das1<-q
mean_result_Das1$effect<-Das1$sigma * Das1$stdev.unscaled
mean_result_Das1$ci_95<-Das1$coefficients * qt(0.975 , Das1$df.residual)
mean_result_Das1$value<-Das1$coefficients

p<-as.data.frame(PSI_merged[,c("Unstable_cdc34_SCF_rep2","Unstable_cdc34_SCF_rep3")])
cdc34<-lmFit(p[,c(1,2)])
mean_result_cdc34<-p
mean_result_cdc34$effect<-cdc34$sigma * cdc34$stdev.unscaled
mean_result_cdc34$ci_95<-cdc34$coefficients * qt(0.975 , cdc34$df.residual)
mean_result_cdc34$value<-cdc34$coefficients


mean_result_cdc34$translation<-row.names(mean_result_cdc34)
mean_result_Das1$translation<-row.names(mean_result_Das1)
Das1_cdc34<-merge(mean_result_Das1, mean_result_cdc34, by= "translation")

pdf(paste0(path_plot,"FigureS5/PSI_cdc34_Das1_afterLMFit.pdf"))
ggplot(Das1_cdc34)+
  geom_point(aes(x = value.x, y = value.y))+
  theme_bw()+
  xlab("Das1 ")+
  ylab("cdc34")+
  xlim(0,1)+
  ylim(0,1)+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(Das1_cdc34$value.x,Das1_cdc34$value.y),2)),
           colour = "red")+
  ggtitle("PSI cdc34 versus Das1 after lmFit")
dev.off()

#only_instable<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/4726.xlsx")
Das1_cdc34$group<-"Stable"
Das1_cdc34$group<-ifelse(Das1_cdc34$translation %in% only_instable$translation, "Unstable", Das1_cdc34$group)
#pdf(paste0(path_plot,"Das1_cdc34_stable_highlighted.pdf"))
Das1_cdc34 %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = group))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stable" = "#8d99ae",
    #"Intermediate" = "#8d99ae",
    "Unstable" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (Das1)")+
  ylab ("PSI (cdc34)")+
  ggtitle("Highlighting stable peptides in lmFit data")
#dev.off()

### remove dataset with stable controls in the plots
Das1_cdc34_unstable<-subset(Das1_cdc34,Das1_cdc34$group == "Unstable")


#### -------------------------------------- different category --------------------------------------------------------------###

q_data<-Das1_cdc34_unstable[,c("translation" ,"Unstable_Das1_SCF_rep2","Unstable_Das1_SCF_rep3","Unstable_cdc34_SCF_rep2","Unstable_cdc34_SCF_rep3")]
q_data<-melt(q_data,id = "translation")
q_data$group<-str_split_fixed(q_data$variable,"_",3)[,2]
stat.test <- q_data %>%
  group_by(translation) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

Das1_cdc34_unstable$differnce<-Das1_cdc34_unstable$value.y - Das1_cdc34_unstable$value.x

Das1_cdc34_unstable$color<-ifelse(Das1_cdc34_unstable$differnce > 0.1, "Potential degron","No effect")
Das1_cdc34_unstable$category<-"No Effect"
Das1_cdc34_unstable$category<-ifelse((Das1_cdc34_unstable$color == "Potential degron" & 
                                        Das1_cdc34_unstable$value.y < 0.67),"Intermediate", Das1_cdc34_unstable$category)
Das1_cdc34_unstable$category<-ifelse((Das1_cdc34_unstable$color == "Potential degron" & Das1_cdc34_unstable$value.y > 0.67 ),"Stabilized", Das1_cdc34_unstable$category)
#Das1_cdc34_unstable$category<-ifelse((Das1_cdc34_unstable$color == "Potential degron" & Das1_cdc34_unstable$value.y > 0.58 ),"Strongly Stabilized", Das1_cdc34_unstable$category)
Das1_cdc34_unstable$category<-ifelse((Das1_cdc34_unstable$color == "Potential degron" & Das1_cdc34_unstable$value.x > 0.5 ),"No Effect", Das1_cdc34_unstable$category)



data_with_stat_test<-merge(stat.test, Das1_cdc34_unstable, id  = "translation")
data_with_stat_test$category<-ifelse(data_with_stat_test$p.adj < 0.05,"No Effect" ,data_with_stat_test$category)

write.csv(data_with_stat_test,
          paste0(path_output,"data_with_stat_test_Das1_cdc34.csv")
)
write.csv(Das1_cdc34_unstable,
          paste0(path_output,"Das1_cdc34_unstable.csv")
)

# cdc53 -wt
q<-as.data.frame(PSI_merged[,c("Unstable_WT_SCF_rep2","Unstable_WT_SCF_rep3")])
WT<-lmFit(q[,c(1,2)])
mean_result_WT<-q
mean_result_WT$effect<-WT$sigma * WT$stdev.unscaled
mean_result_WT$ci_95<-WT$coefficients * qt(0.975 , WT$df.residual)
mean_result_WT$value<-WT$coefficients

p<-as.data.frame(PSI_merged[,c("Unstable_cdc53_SCF_rep2","Unstable_cdc53_SCF_rep3")])
cdc53<-lmFit(p[,c(1,2)])
mean_result_cdc53<-p
mean_result_cdc53$effect<-cdc53$sigma * cdc53$stdev.unscaled
mean_result_cdc53$ci_95<-cdc53$coefficients * qt(0.975 , cdc53$df.residual)
mean_result_cdc53$value<-cdc53$coefficients


mean_result_cdc53$translation<-row.names(mean_result_cdc53)
mean_result_WT$translation<-row.names(mean_result_WT)
WT_cdc53<-merge(mean_result_WT, mean_result_cdc53, by= "translation")

#pdf(paste0(path_plot,"PSI_cdc53_WT_afterLMFit.pdf"))
ggplot(WT_cdc53)+
  geom_point(aes(x = value.x, y = value.y))+
  theme_bw()+
  xlab("Wild type ")+
  ylab("cdc53")+
  xlim(0,1)+
  ylim(0,1)+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(WT_cdc53$value.x,WT_cdc53$value.y),2)),
           colour = "red")+
  ggtitle("PSI cdc53 versus WT after lmFit")
#dev.off()

#only_instable<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/4726.xlsx")
WT_cdc53$group<-"Stable"
WT_cdc53$group<-ifelse(WT_cdc53$translation %in% only_instable$translation, "Unstable", WT_cdc53$group)
#pdf(paste0(path_plot,"WT_cdc53_stable_highlighted.pdf"))
WT_cdc53 %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = group))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stable" = "#8d99ae",
    #"Intermediate" = "#8d99ae",
    "Unstable" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (cdc53)")+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(WT_cdc53$value.x,WT_cdc53$value.y),2)),
           colour = "red")+
  ggtitle("Highlighting stable peptides in lmFit data")
#dev.off()

### remove dataset with stable controls in the plots
WT_cdc53_unstable<-subset(WT_cdc53,WT_cdc53$group == "Unstable")
q_data<-WT_cdc53_unstable[,c("translation" ,"Unstable_WT_SCF_rep2","Unstable_WT_SCF_rep3","Unstable_cdc53_SCF_rep2","Unstable_cdc53_SCF_rep3")]
q_data<-melt(q_data,id = "translation")
q_data$group<-str_split_fixed(q_data$variable,"_",3)[,2]
stat.test <- q_data %>%
  group_by(translation) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

WT_cdc53_unstable$differnce<-WT_cdc53_unstable$value.y - WT_cdc53_unstable$value.x

WT_cdc53_unstable$color<-ifelse(WT_cdc53_unstable$differnce > 0.1, "Potential degron","No effect")
WT_cdc53_unstable$category<-"No Effect"
WT_cdc53_unstable$category<-ifelse((WT_cdc53_unstable$color == "Potential degron" & 
                                      WT_cdc53_unstable$value.y < 0.67),"Intermediate", WT_cdc53_unstable$category)
WT_cdc53_unstable$category<-ifelse((WT_cdc53_unstable$color == "Potential degron" & WT_cdc53_unstable$value.y > 0.67 ),"Stabilized", WT_cdc53_unstable$category)
#WT_cdc53_unstable$category<-ifelse((WT_cdc53_unstable$color == "Potential degron" & WT_cdc53_unstable$value.y > 0.58 ),"Strongly Stabilized", WT_cdc53_unstable$category)
WT_cdc53_unstable$category<-ifelse((WT_cdc53_unstable$color == "Potential degron" & WT_cdc53_unstable$value.x > 0.5 ),"No Effect", WT_cdc53_unstable$category)



data_with_stat_test_Das1<-read.csv(paste0(path_output,"data_with_stat_test_wt_das1.csv"))

Das1_group<-data_with_stat_test_Das1[,c("translation","category")]
names(Das1_group)<-c("translation","category_Das1")
WT_cdc53_unstable<-merge(WT_cdc53_unstable,Das1_group, by = "translation")

pdf(paste0(path_plot,"FigureS5/Das1Peptides_highlighted_in_cdc53WT.pdf"))
ggplot(WT_cdc53_unstable)+
  
  geom_point(data = subset(WT_cdc53_unstable,WT_cdc53_unstable$category_Das1 == "No Effect"),
             aes(x = value.x, y = value.y), color = "black")+
  geom_point(data = subset(WT_cdc53_unstable,WT_cdc53_unstable$category_Das1 == "Intermediate"),
             aes(x = value.x, y = value.y), color = "grey")+
  geom_point(data = subset(WT_cdc53_unstable,WT_cdc53_unstable$category_Das1 == "Stabilized"),
             aes(x = value.x, y = value.y), color = "red")+
  
  theme_bw()+
  #scale_color_manual(values = c( 
  # #"Strongly Stabilized" = "#98473E",
  #"Stabilized" = "#ff5400",
  #"Intermediate" = "#8d99ae",
  #"No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (WT)")+
  ylab ("PSI (cdc53)")+
  ggtitle("Marking Das1 substrates peptides in cdc53-WT variation")
dev.off()

#cdc53-das1

q<-as.data.frame(PSI_merged[,c("Unstable_Das1_SCF_rep2","Unstable_Das1_SCF_rep3")])
Das1<-lmFit(q[,c(1,2)])
mean_result_Das1<-q
mean_result_Das1$effect<-Das1$sigma * Das1$stdev.unscaled
mean_result_Das1$ci_95<-Das1$coefficients * qt(0.975 , Das1$df.residual)
mean_result_Das1$value<-Das1$coefficients

p<-as.data.frame(PSI_merged[,c("Unstable_cdc53_SCF_rep2","Unstable_cdc53_SCF_rep3")])
cdc53<-lmFit(p[,c(1,2)])
mean_result_cdc53<-p
mean_result_cdc53$effect<-cdc53$sigma * cdc53$stdev.unscaled
mean_result_cdc53$ci_95<-cdc53$coefficients * qt(0.975 , cdc53$df.residual)
mean_result_cdc53$value<-cdc53$coefficients


mean_result_cdc53$translation<-row.names(mean_result_cdc53)
mean_result_Das1$translation<-row.names(mean_result_Das1)
Das1_cdc53<-merge(mean_result_Das1, mean_result_cdc53, by= "translation")

#pdf(paste0(path_plot,"PSI_cdc53_Das1_afterLMFit.pdf"))
ggplot(Das1_cdc53)+
  geom_point(aes(x = value.x, y = value.y))+
  theme_bw()+
  xlab("Wild type ")+
  ylab("cdc53")+
  xlim(0,1)+
  ylim(0,1)+
  annotate("text",
           x = 0.05, y = 0.9, 
           label = paste0("R  : ",round(cor(Das1_cdc53$value.x,Das1_cdc53$value.y),2)),
           colour = "red")+
  ggtitle("PSI cdc53 versus Das1 after lmFit")
#dev.off()

#only_instable<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/4726.xlsx")
Das1_cdc53$group<-"Stable"
Das1_cdc53$group<-ifelse(Das1_cdc53$translation %in% only_instable$translation, "Unstable", Das1_cdc53$group)
#pdf(paste0(path_plot,"Das1_cdc53_stable_highlighted.pdf"))
Das1_cdc53 %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = group))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stable" = "#8d99ae",
    #"Intermediate" = "#8d99ae",
    "Unstable" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (Das1)")+
  ylab ("PSI (cdc53)")+
  ggtitle("Highlighting stable peptides in lmFit data")
#dev.off()

### remove dataset with stable controls in the plots
Das1_cdc53_unstable<-subset(Das1_cdc53,Das1_cdc53$group == "Unstable")


#### -------------------------------------- different category --------------------------------------------------------------###

q_data<-Das1_cdc53_unstable[,c("translation" ,"Unstable_Das1_SCF_rep2","Unstable_Das1_SCF_rep3","Unstable_cdc53_SCF_rep2","Unstable_cdc53_SCF_rep3")]
q_data<-melt(q_data,id = "translation")
q_data$group<-str_split_fixed(q_data$variable,"_",3)[,2]
stat.test <- q_data %>%
  group_by(translation) %>%
  t_test(value ~ group) %>%
  adjust_pvalue(method = "BH") %>%
  add_significance()

Das1_cdc53_unstable$differnce<-Das1_cdc53_unstable$value.y - Das1_cdc53_unstable$value.x

Das1_cdc53_unstable$color<-ifelse(Das1_cdc53_unstable$differnce > 0.1, "Potential degron","No effect")
Das1_cdc53_unstable$category<-"No Effect"
Das1_cdc53_unstable$category<-ifelse((Das1_cdc53_unstable$color == "Potential degron" & 
                                        Das1_cdc53_unstable$value.y < 0.67),"Intermediate", Das1_cdc53_unstable$category)
Das1_cdc53_unstable$category<-ifelse((Das1_cdc53_unstable$color == "Potential degron" & Das1_cdc53_unstable$value.y > 0.67 ),"Stabilized", Das1_cdc53_unstable$category)
#Das1_cdc53_unstable$category<-ifelse((Das1_cdc53_unstable$color == "Potential degron" & Das1_cdc53_unstable$value.y > 0.58 ),"Strongly Stabilized", Das1_cdc53_unstable$category)
Das1_cdc53_unstable$category<-ifelse((Das1_cdc53_unstable$color == "Potential degron" & Das1_cdc53_unstable$value.x > 0.5 ),"No Effect", Das1_cdc53_unstable$category)



data_with_stat_test<-merge(stat.test, Das1_cdc53_unstable, id  = "translation")
data_with_stat_test$category<-ifelse(data_with_stat_test$p.adj < 0.05,"No Effect" ,data_with_stat_test$category)
#pdf(paste0(path_plot,"Das1_cdc53_3Category.pdf"))
data_with_stat_test %>% ggplot(aes(x = value.x, y = value.y))+
  geom_point(aes(color = category))+
  theme_bw()+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stabilized" = "#ff5400",
    "Intermediate" = "#8d99ae",
    "No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (Das1)")+
  ylab ("PSI (cdc53)")+
  ggtitle("Putting statistically insignificant points as No Effect  with multiple testing")+
  annotate("text", 
           x = 0.15, y = 0.92, 
           label = paste0("Threshold : \n Difference :0.1, \nPSI Das1: 0.5; \nPSI cdc53 : 0.67"),
           colour = "red", fontface =2)
#dev.off()
# hydrophobicity


data_with_stat_test$hydrophobicity<-hydrophobicity(data_with_stat_test$translation)

my_comparisons <- list( c("No Effect", "Stabilized"), c("Stabilized", "Intermediate"), c("No Effect", "Intermediate") )
#pdf(paste0(path_plot,"3Category_hydrophobicity.pdf"))
ggplot(data_with_stat_test)+
  geom_boxplot(
    aes(x = category, 
        y = hydrophobicity, color = category),
    coef = 1e30,lwd=0.25
  )+
  ylim(-4,4)+
  scale_color_manual(values = c( 
    #"Strongly Stabilized" = "#98473E",
    "Stabilized" = "#ff5400",
    "Intermediate" = "#8d99ae",
    "No Effect" = "#2b2d42"))+
  guides(col = FALSE)+
  ylab("Hydrophobicity")+
  xlab("")+
  ggtitle("Variation of hydrophobicity")+
  theme_bw()
#dev.off()

data_with_stat_test_Das1<-read.csv(paste0(path_output,"data_with_stat_test_WT_das1.csv"))

Das1_group<-data_with_stat_test_Das1[,c("translation","category")]
names(Das1_group)<-c("translation","category_Das1")
Das1_cdc53_unstable<-merge(Das1_cdc53_unstable,Das1_group, by = "translation")

#pdf(paste0(path_plot,"Das1Peptides_highlighted_in_cdc53Das1.pdf"))
ggplot(Das1_cdc53_unstable)+
  
  geom_point(data = Das1_cdc53_unstable,aes(x = value.x, y = value.y), color = "black")+
  geom_point(data = subset(Das1_cdc53_unstable,Das1_cdc53_unstable$category_Das1 == "Intermediate"),
             aes(x = value.x, y = value.y), color = "grey")+
  geom_point(data = subset(Das1_cdc53_unstable,Das1_cdc53_unstable$category_Das1 == "Stabilized"),
             aes(x = value.x, y = value.y), color = "red")+
  #geom_point(data = subset(Das1_cdc53_unstable,Das1_cdc53_unstable$category_Das1 == "Intermediate"),
  #          aes(x = value.x, y = value.y), color = "grey")+
  theme_bw()+
  #scale_color_manual(values = c( 
  # #"Strongly Stabilized" = "#98473E",
  #"Stabilized" = "#ff5400",
  #"Intermediate" = "#8d99ae",
  #"No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (Das1)")+
  ylab ("PSI (cdc53)")+
  ggtitle("Marking Das1 substrates peptides in cdc53-Das1 variation")
#dev.off()

#pdf(paste0(path_plot,"Das1Peptides_highlighted_in_cdc53Das1.pdf"))
ggplot(Das1_cdc53_unstable)+
  
  geom_point(data = subset(Das1_cdc53_unstable,Das1_cdc53_unstable$category_Das1 == "No Effect"),
             aes(x = value.x, y = value.y), color = "black")+
  geom_point(data = subset(Das1_cdc53_unstable,Das1_cdc53_unstable$category_Das1 == "Stabilized"),
             aes(x = value.x, y = value.y), color = "red")+
  geom_point(data = subset(Das1_cdc53_unstable,Das1_cdc53_unstable$category_Das1 == "Intermediate"),
             aes(x = value.x, y = value.y), color = "grey")+
  theme_bw()+
  #scale_color_manual(values = c( 
  # #"Strongly Stabilized" = "#98473E",
  #"Stabilized" = "#ff5400",
  #"Intermediate" = "#8d99ae",
  #"No Effect" = "#2b2d42"))+
  xlim(0,1)+
  ylim(0,1)+
  xlab("PSI (Das1)")+
  ylab ("PSI (cdc53)")+
  ggtitle("Marking Das1 Intermediate peptides in cdc53-Das1 variation")
#dev.off()

write.csv(data_with_stat_test,
          paste0(path_output,"data_with_stat_test_Das1_cdc53.csv")
)
write.csv(Das1_cdc53_unstable,
          paste0(path_output,"Das1_cdc53_unstable.csv")
)

##########


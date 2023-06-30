##### for project : Analysis of degron at C terminal
##### sub text:plotting various distributions unstable peptides
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB


# set the directory
setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")

dataset<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))

dataset<-subset(dataset,dataset$Codon == "Non Stop Codon")
dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]
dataset_for_frequency<-unique(dataset_for_frequency)
dataset_for_frequency$CTerminal<-substr(dataset_for_frequency$raw_counts_translation,12,12)

dir.create(paste0(path_plot,"Figure1/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

# boxplot with C terminal PSI
pdf(paste0(path_plot,"Figure1/boxplot_PSI_Cterminal.pdf"))
ggplot(dataset_for_frequency,aes(x = CTerminal, y = byaa_pooled_PSI))+
  geom_boxplot(coef=1e30)+
  theme_bw()+
  labs(title = "Variation of PSI with respect to C terminal - 12x",
       
       caption = paste0 ("# Peptides with following C Terminal - A : 2205, C : 1364, D : 998, E : 1329, F :2277, G :4371, H :1004, I :1764, K :1775,       
                         L :5038, M :560, V :900, P :6170, Q :824, R :3903, S :3875, T :1898, V :3669, W :1013, Y :855        
                         "))+
  ylab("protein Stability Index")+
  xlab("Amino Acid at C Terminal")
dev.off()

# for frequency of unstable
dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]
dataset_for_frequency<-unique(dataset_for_frequency)
amino_acid_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  amino_acid_frequency[,i]<-
    table(
      substr(dataset_for_frequency$raw_counts_translation,i,i)
    )
}
amino_acid_frequency$AminoAcid<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
amino_acid_frequency<-melt(amino_acid_frequency, id = "AminoAcid")
amino_acid_frequency$value<-as.numeric(amino_acid_frequency$value)
amino_acid_frequency$relativeFrequencyAll<-amino_acid_frequency$value/nrow(dataset_for_frequency)
amino_acid_frequency$variable<-paste0(amino_acid_frequency$variable,amino_acid_frequency$AminoAcid)
amino_acid_frequency$AminoAcid<-NULL
amino_acid_frequency$value<-NULL

dataset_for_frequency$Label<-ifelse(dataset_for_frequency$byaa_pooled_PSI <0.5, 
                                    "Unstable",
                                    "Stable")
Unstable<-subset(dataset_for_frequency,
                 dataset_for_frequency$Label == "Unstable")
Unstable_aa_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  Unstable_aa_frequency[,i]<-
    table(
      substr(Unstable$raw_counts_translation,i,i)
    )
}
Unstable_aa_frequency$AminoAcid<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
Unstable_aa_frequency<-melt(Unstable_aa_frequency, id = "AminoAcid")
Unstable_aa_frequency$value<-as.numeric(Unstable_aa_frequency$value)
Unstable_aa_frequency$relativeFrequencyUnstable<-Unstable_aa_frequency$value/nrow(Unstable)
Unstable_aa_frequency$variable<-paste0(Unstable_aa_frequency$variable,Unstable_aa_frequency$AminoAcid)
Unstable_aa_frequency$AminoAcid<-NULL
Unstable_aa_frequency$value<-NULL

heatmapUnstable<-merge(amino_acid_frequency,Unstable_aa_frequency, id = "AminoAcid")
heatmapUnstable$normalizedFrequency<-heatmapUnstable$relativeFrequencyUnstable/heatmapUnstable$relativeFrequencyAll
heatmapUnstable$position<-substr(heatmapUnstable$variable,2,(nchar(heatmapUnstable$variable)-1))
heatmapUnstable$position<-as.numeric(heatmapUnstable$position)-13
heatmapUnstable$AminoAcid<-substr(heatmapUnstable$variable,(nchar(heatmapUnstable$variable)),nchar(heatmapUnstable$variable))

heatmapUnstable$AminoAcid<-as.character(heatmapUnstable$AminoAcid)
heatmapUnstable$position<-factor(heatmapUnstable$position, 
                                 levels = c("-12","-11","-10","-9","-8","-7","-6","-5","-4","-3","-2","-1"))
heatmapUnstable$logEF<-log2(heatmapUnstable$normalizedFrequency)

# Heatmap with Normalized Frequency (normalized to 12x) for unstable peptides

pdf(paste0(path_plot,"Figure1/NormalizedFrequency.pdf"))
ggplot(heatmapUnstable, aes (x = position,
                          y= AminoAcid,
                          fill = logEF))+
  geom_tile(data = subset(heatmapUnstable,heatmapUnstable$logEF < 1.5 & 
                            heatmapUnstable$logEF > -1.5), 
            aes(
              x = position,
              y= AminoAcid,
              fill = logEF
            ))+

  geom_tile(data = subset(heatmapUnstable,heatmapUnstable$logEF < -1.5 ), 
            aes(
              x = position,
              y= AminoAcid
            ),
            fill = "#F7931E")+
  geom_tile(data = subset(heatmapUnstable,heatmapUnstable$logEF > 1.5 ), 
            aes(
              x = position,
              y= AminoAcid
            ),
            fill = "#0000FF")+
  geom_tile(data = subset(heatmapUnstable,is.na(heatmapUnstable$logEF)), 
            aes(
              x = position,
              y= AminoAcid
            ),
            fill = "#F7931E")+
  geom_text(aes(label = ifelse(is.na(heatmapUnstable$logEF), "*", ""),), 
            color = "Black", 
            size = 2) +
  coord_equal()+
  theme_bw()+
  theme(text = element_text(size = 8))+
  scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
  scale_fill_gradientn("log2(Freq_heatmapUnstable)/\n (Freq_12x)",
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
                                 "#0000FF",
                                 "#0000FF"),
                       values=scales::rescale(c(-1.50 ,-1.25, -1.00, -0.75, -0.50, -0.25,  0.00  ,0.25  ,0.50  ,0.75  ,1.00  ,1.25  ,1.50)),
                       limits=c(-1.5,1.5))


dev.off()
#### variation of hydrophobicity with PSI


datasetcomplete<-read.csv(paste0(path_output,"complete_dataset_without_normalization.csv"))
d<-datasetcomplete[,c("hydrophobicity_overall","PSI")]
d$translation<-datasetcomplete$X
d<-subset(d,d$translation %in% dataset$raw_counts_translation)
d$Label<-0
d$Label<-ifelse(d$PSI < 0.5, 1,d$Label)

d$hydrophobicity_overall<-as.numeric(d$hydrophobicity_overall)
d$PSI<-as.numeric(d$PSI)
d$Label<-as.factor(d$Label)

# plot for distribution of hydrophobicity versus PSI
PSI_hydro<-d %>%
  ggplot(aes(x = PSI ,y =hydrophobicity_overall))+ 
  geom_point(aes(color = Label),alpha = 0.8, shape = 19)+
  theme_bw()+
  scale_color_manual(values=c("#999999", "#4290C0"))+
  xlab("Protein Stability Index")+
  ylab(" Overall Hydrophobicity")+
  xlim(0,1)+
  ylim(-4,4)+
  annotate("text", x = 0.1, y = 3.9, label = paste0("R : ", round(cor(
    d$hydrophobicity_overall,
    d$PSI, method = c("spearman")
  ),3)))

densi<-d %>%
  ggplot(aes(y = hydrophobicity_overall, color = Label))+ 
  geom_density(kernel = "rectangular")+
  theme_bw()+
  scale_color_manual(values=c("#999999", "#4290C0"))+
  xlab("Density")+
  ylab(" Overall Hydrophobicity")+
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  ylim(-4,4)

PSI_hydro <- PSI_hydro + rremove("legend")

pdf(paste0(path_plot,"Figure1/Hydrophobicity_versus_PSI.pdf"))
print(
  plot_grid(PSI_hydro,densi)
)
dev.off()

###################### plotted Figure 1




# for C terminal di peptide frequency unstable

dir.create(paste0(path_plot,"Figure2/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]
dataset_for_frequency<-unique(dataset_for_frequency)
dataset_for_frequency$Label<-ifelse(dataset_for_frequency$byaa_pooled_PSI <0.5, "Unstable","Stable")
Unstable<-subset(dataset_for_frequency,dataset_for_frequency$Label == "Unstable")

dipeptide_12x<-as.data.frame(table(
  substr(dataset_for_frequency$raw_counts_translation,11,12)
)
)
dipeptide_4k<-as.data.frame(table(
  substr(Unstable$raw_counts_translation,11,12)
)
)
di<-merge(dipeptide_4k,dipeptide_12x, by = "Var1", all.y = TRUE)
names(di)<-c("dipeptide","Freq_Unstable","Freq_12x")
di$normalizedFreq_unstable<-di$Freq_Unstable/nrow(Unstable)
di$normalizedFreq_all<-di$Freq_12x/nrow(dataset_for_frequency)
di$log2NF<-log2(di$normalizedFreq_unstable/di$normalizedFreq_all)
di$NF<-(di$normalizedFreq_unstable/di$normalizedFreq_all)
di$Amino_acid_1<-substr(di$dipeptide,1,1)
di$Amino_acid_2<-substr(di$dipeptide,2,2)
di$Amino_acid_1<-factor(di$Amino_acid_1, levels = (c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
di$Amino_acid_2<-factor(di$Amino_acid_2, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
di$logNF<-log2(di$NF)
di2<-di

# normalized frequency for dipeptide for unstable
pdf(paste0(path_plot,"Figure2/dipeptide_NormalizedFrequency.pdf"))
ggplot(di, aes (x = Amino_acid_1,
                             y= Amino_acid_2,
                             fill = logNF))+
  geom_tile(data = subset(di,di$logNF < 2 & 
                            di$logNF > -2), 
            aes(
              x = Amino_acid_1,
              y= Amino_acid_2,
              fill = logNF
            ))+
  
  geom_tile(data = subset(di,di$logNF < -2 ), 
            aes(
              x = Amino_acid_1,
              y= Amino_acid_2
            ),
            fill = "#F7931E")+
  geom_tile(data = subset(di,di$logNF > 2 ), 
            aes(
              x = Amino_acid_1,
              y= Amino_acid_2
            ),
            fill = "#0000FF")+
  geom_tile(data = subset(di,is.na(di$logNF)), 
            aes(
              x = Amino_acid_1,
              y= Amino_acid_2
            ),
            fill = "#F7931E")+
  geom_text(aes(label = ifelse(is.na(di$logNF), "*", ""),), 
            color = "Black", 
            size = 2) +
  xlab("Amino Acid (-2)")+
  ylab("Amino Acid (-1)")+
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
                                 "#0000FF",
                                 "#0000FF"),
                       values=scales::rescale(c(-2.00, -1.67, -1.33 ,-1.00 ,-0.67 ,-0.33 , 0.00 , 0.33  ,0.67  ,1.00 , 1.33,  1.67 , 2.00)),
                       limits=c(-2,2))
dev.off()

# other plots are from Figure S2


#############



# plot for the selected unstable sequence  


peptides<-c("LVVYLSCVACIN","KKHKKGQKQKIN","PDKRNHLIKLLN",
            "QGRRKEKAVTLN","RLIQPLPVQFLN","KCKIWFREVSMN",
            "DSSRRKRLICMN","CLFHQFITLLVN","FSMRWTFRYWVN",
            "IIGCIAICLMVN","ELSRLLPASFSL"
)

d$Label1<-ifelse(d$translation %in% peptides, "tested degron","Others")
d$Label2<-"Others"
d$Label2<-ifelse(d$translation == "LVVYLSCVACIN", "IN1",d$Label2)
d$Label2<-ifelse(d$translation == "KKHKKGQKQKIN", "IN2",d$Label2)
d$Label2<-ifelse(d$translation == "PDKRNHLIKLLN", "LN1",d$Label2)
d$Label2<-ifelse(d$translation == "QGRRKEKAVTLN", "LN2",d$Label2)
d$Label2<-ifelse(d$translation == "RLIQPLPVQFLN", "LN3",d$Label2)
d$Label2<-ifelse(d$translation == "KCKIWFREVSMN", "MN1",d$Label2)
d$Label2<-ifelse(d$translation == "DSSRRKRLICMN", "MN2",d$Label2)
d$Label2<-ifelse(d$translation == "CLFHQFITLLVN", "VN1",d$Label2)
d$Label2<-ifelse(d$translation == "FSMRWTFRYWVN", "VN2",d$Label2)
d$Label2<-ifelse(d$translation == "IIGCIAICLMVN", "VN3",d$Label2)
d$Label2<-ifelse(d$translation == "ELSRLLPASFSL", "XS",d$Label2)

d_sub<-subset(d,!(d$Label2 == "Others"))
dir.create(paste0(path_plot,"Figure3/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

pdf(paste0(path_plot,"Figure3/highlighted_N.pdf"))
print(
ggplot(d)+
  geom_point(data = subset(d,d$Label1 == "Others"), 
             aes(x = PSI, y = hydrophobicity_overall), 
             color = "grey", 
             alpha = 0.2)+
  geom_point(data = d_sub, 
             aes(x = PSI, y = hydrophobicity_overall), 
             color = "red", 
             alpha = 0.6)+
  geom_text_repel(data=d_sub,aes(x=PSI,y=hydrophobicity_overall,label=Label2), color='red', size = 3.5)+
  theme_bw()+
  #scale_color_manual(values=c("#999999", "#4290C0"))+
  #geom_smooth(method=lm, se=FALSE, linetype = "dashed", color = "black")+
  xlab("Protein Stability Index")+
  ylab(" Overall Hydrophobicity")+
  xlim(0,1)+
  ylim(-4,4)+
  annotate("text", x = 0.1, y = 3.9, label = paste0("R : ", round(cor(
    d$hydrophobicity_overall,
    d$PSI, method = c("spearman")
  ),3)))+
  ggtitle("PSI versus Hydrophobicity - 12x"))
dev.off()



# --------------------------------------------------------------------------------------------------------------------------------------------------------
# --------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------------------------------------------
### plot for checking the distribution across different positions - need not run
dataset_post_removal<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))
dataset_post_removal<-subset(dataset_post_removal,dataset_post_removal$Codon == "Non Stop Codon")
# --------------------------------------------------------------------------------------------------------------------------------------------------------
#pdf("Y:/lab data/susmitha/edwin/for_paper/plot/CTer_Unstable_NormFreq.pdf")


dipeptide_12x_1<-as.data.frame(table(
  substr(dataset_for_frequency$raw_counts_translation,10,11)
)
)
dipeptide_4k_1<-as.data.frame(table(
  substr(Unstable$raw_counts_translation,10,11)
)
)
di<-merge(dipeptide_4k_1,dipeptide_12x_1, by = "Var1", all.y = TRUE)
names(di)<-c("dipeptide","Freq_Unstable","Freq_12x")
di$normalizedFreq_unstable<-di$Freq_Unstable/nrow(Unstable)
di$normalizedFreq_all<-di$Freq_12x/nrow(dataset_for_frequency)
di$log2NF<-log2(di$normalizedFreq_unstable/di$normalizedFreq_all)
di$NF<-(di$normalizedFreq_unstable/di$normalizedFreq_all)
di$Amino_acid_1<-substr(di$dipeptide,1,1)
di$Amino_acid_2<-substr(di$dipeptide,2,2)
di$Amino_acid_1<-factor(di$Amino_acid_1, levels = (c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
di$Amino_acid_2<-factor(di$Amino_acid_2, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))


dataset_post_removal$Cterdi<-substr(dataset_post_removal$byaa_translation,11,12)
dataset_post_removal$a1<-substr(dataset_post_removal$byaa_translation,10,11)
dataset_post_removal$a2<-substr(dataset_post_removal$byaa_translation,9,10)
dataset_post_removal$a3<-substr(dataset_post_removal$byaa_translation,8,9)
dataset_post_removal$a4<-substr(dataset_post_removal$byaa_translation,7,8)
dataset_post_removal$a5<-substr(dataset_post_removal$byaa_translation,6,7)
dataset_post_removal$a6<-substr(dataset_post_removal$byaa_translation,5,6)
dataset_post_removal$a7<-substr(dataset_post_removal$byaa_translation,4,5)
dataset_post_removal$a8<-substr(dataset_post_removal$byaa_translation,3,4)
dataset_post_removal$a9<-substr(dataset_post_removal$byaa_translation,2,3)
dataset_post_removal$a10<-substr(dataset_post_removal$byaa_translation,1,2)

CTer_psi<-dataset_post_removal %>% group_by(Cterdi) %>% summarise(mean = mean(byaa_pooled_PSI))

a1_psi<-dataset_post_removal %>% group_by(a1) %>% summarise(mean_a1 = mean(byaa_pooled_PSI))
a2_psi<-dataset_post_removal %>% group_by(a2) %>% summarise(mean_a2 = mean(byaa_pooled_PSI))
a3_psi<-dataset_post_removal %>% group_by(a3) %>% summarise(mean_a3 = mean(byaa_pooled_PSI))
a4_psi<-dataset_post_removal %>% group_by(a4) %>% summarise(mean_a4 = mean(byaa_pooled_PSI))
a5_psi<-dataset_post_removal %>% group_by(a5) %>% summarise(mean_a5 = mean(byaa_pooled_PSI))
a6_psi<-dataset_post_removal %>% group_by(a6) %>% summarise(mean_a6 = mean(byaa_pooled_PSI))
a7_psi<-dataset_post_removal %>% group_by(a7) %>% summarise(mean_a7 = mean(byaa_pooled_PSI))
a8_psi<-dataset_post_removal %>% group_by(a8) %>% summarise(mean_a8 = mean(byaa_pooled_PSI))
a9_psi<-dataset_post_removal %>% group_by(a9) %>% summarise(mean_a9 = mean(byaa_pooled_PSI))
a10_psi<-dataset_post_removal %>% group_by(a10) %>% summarise(mean_a10 = mean(byaa_pooled_PSI))


names(CTer_psi)<-c("Di","Mean_PSI_Cter")
names(a1_psi)<-c("Di","Mean_PSI_a1")
names(a2_psi)<-c("Di","Mean_PSI_a2")
names(a3_psi)<-c("Di","Mean_PSI_a3")
names(a4_psi)<-c("Di","Mean_PSI_a4")
names(a5_psi)<-c("Di","Mean_PSI_a5")
names(a6_psi)<-c("Di","Mean_PSI_a6")
names(a7_psi)<-c("Di","Mean_PSI_a7")
names(a8_psi)<-c("Di","Mean_PSI_a8")
names(a9_psi)<-c("Di","Mean_PSI_a9")
names(a10_psi)<-c("Di","Mean_PSI_a10")


psi<-merge(CTer_psi,a1_psi, id = "Di")
psi<-merge(psi,a2_psi, id = "Di")
psi<-merge(psi,a3_psi, id = "Di")
psi<-merge(psi,a4_psi, id = "Di")
psi<-merge(psi,a5_psi, id = "Di")
psi<-merge(psi,a6_psi, id = "Di")
psi<-merge(psi,a7_psi, id = "Di")
psi<-merge(psi,a8_psi, id = "Di")
psi<-merge(psi,a9_psi, id = "Di")
psi<-merge(psi,a10_psi, id = "Di")
psi$means<-rowMeans(psi[,3:ncol(psi)])

psi$difference<-psi$Mean_PSI_Cter - psi$means
psi$Amino_acid_1<-substr(psi$Di,1,1)
psi$Amino_acid_2<-substr(psi$Di,2,2)
di$Amino_acid_1<-factor(psi$Amino_acid_1, levels = (c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
psi$Amino_acid_2<-factor(psi$Amino_acid_2, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
#psi$diff_req<-psi$Mean_PSI_Cter - psi$Mean_PSI_a5
#psi$diff_req1<-psi$Mean_PSI_Cter - psi$Mean_PSI_a6

psi %>% 
  ggplot(aes(x = Amino_acid_1, 
             y = Amino_acid_2))+
  geom_tile(aes(fill = difference), alpha = 3)+
  scale_fill_gradientn("PSI-Cter",colours=c(
    "#F7931E", 
    "#F5AB53", 
    "#F4BB75", 
    "#F3CB98", 
    "#F1DBBB", 
    "#FFFFFF", 
    "#C4C0E4", 
    "#9996EA", 
    "#6D6BF0", 
    "#4240F6", 
    "#0000FF"
  ), na.value = "grey98",
  values = scales::rescale(round(seq(-0.2,0.2,0.4/10),2)),
  limits = c(-0.2,0.2))+
  #scale_fill_continuous(, )+
  coord_equal()+
  theme_bw()+
  xlab("Amino Acid 1 (X)")+
  ylab("Amino Acid 2 (Z)")

write.csv(psi, 
          paste0("Y:/lab data/susmitha/edwin/for_paper/new_dataset/plots/new_dataset_plots/","psi_diff.csv") ) 

psi<-read.csv( paste0("Y:/lab data/susmitha/edwin/for_paper/new_dataset/plots/new_dataset_plots/","psi_diff.csv"))

# --------------------------------------------------------------------------------------------------------------------------------------------------------

# --------------------------------------------------------------------------------------------------------------------------------------------------------
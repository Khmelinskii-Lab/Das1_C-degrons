
##### for project : Analysis of degron at C terminal
##### sub text:inital analysis
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")

dir.create(paste0(path_plot,"FigureS1/"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

replicate_1<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))
replicate_1<-subset(replicate_1,replicate_1$Codon == "Non Stop Codon")
replicate_2<-read.csv(paste0(path_output,"merged_file_post_remoal_replicate2.csv"))
replicate_2<-subset(replicate_2,replicate_2$Codon == "Non Stop Codon")

common<-subset(replicate_2,replicate_2$raw_counts_translation %in% replicate_1$raw_counts_translation)

common<-common[,c("byaa_translation","byaa_pooled_PSI")]
common2<-replicate_1[,c("byaa_translation","byaa_pooled_PSI")]
common<-unique(common)
common2<-unique(common2)

for_reproducibility<-merge(x = common, y = common2, by = c("byaa_translation"))
names(for_reproducibility)<-c("byaa_translation","PSI_replicate2","PSI_replicate1")
for_reproducibility<-for_reproducibility[complete.cases(for_reproducibility),]
for_reproducibility<-for_reproducibility[!duplicated(for_reproducibility$byaa_translation),]

correlation<-cor(for_reproducibility$PSI_replicate1,for_reproducibility$PSI_replicate2, method = c("pearson"))
pdf(paste0(path_plot,"FigureS1/PSI_replicate_12x.pdf"))
ggplot(for_reproducibility)+
  aes(x = PSI_replicate1,
      y = PSI_replicate2)+
  geom_point()+
  theme_bw()+
  xlab("PSI (Replicate 1)")+
  ylab("PSI (Replicate 2)")+
  annotate(geom = "text",
           x = 0.1,
           y = 0.9,
           label = paste0("#Peptides : ",
                          nrow(for_reproducibility), 
                          " \n ","R : "  , 
                          round(correlation,3)
           ))
dev.off()

dataset<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))

dataset<-subset(dataset,dataset$Codon == "Non Stop Codon")
names(dataset)
dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]
dataset_for_frequency<-unique(dataset_for_frequency)

amino_acid_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  amino_acid_frequency[,i]<-
    table(
      substr(dataset_for_frequency$raw_counts_translation,i,i)
    )
}

amino_acid_frequency<-as.data.frame(amino_acid_frequency[,12])
amino_acid_frequency$AminoAcid<-c( "A" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,"N", "P" ,"Q", "R", "S" ,"T" ,"V" ,"W" ,"Y")
amino_acid_frequency$Var1<-NULL
names(amino_acid_frequency)<-c("FreqAll","AminoAcid")
amino_acid_frequency$FreqAll<-as.numeric(amino_acid_frequency$FreqAll)
amino_acid_frequency$relativeFreqAll<-amino_acid_frequency$FreqAll/46152


Codon<-data.frame(
  Full_Name = c("Isoleucine", "Leucine","Valine","Phenylalanine","Methionine","Cysteine","Alanine","Glycine","Proline","Threonine","Serine","Tyrosine","Tryptophan","Glutamine","Asparagine","Histidine","Glutamic acid","Aspartic acid","Lysine","Arginine"),
  AminoAcid = c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R"),
  Codon_Number = c(3,6,4,2,1,2,4,4,4,4,6,2,1,2,2,2,2,2,2,6)
)
Codon$expectedFrequency<-(Codon$Codon_Number * 46152)/61

Codon<-Codon[order(Codon$AminoAcid),]
amino_acid_frequency<-merge(amino_acid_frequency, Codon, by = "AminoAcid")
amino_acid_frequency<-amino_acid_frequency[,c("AminoAcid","FreqAll","expectedFrequency")]
amino_acid_frequency<-melt(amino_acid_frequency, id = "AminoAcid")

amino_acid_frequency %>% ggplot(aes(x = AminoAcid, y = value, color = variable,group = variable))+
  geom_point()+
  geom_line()+
  theme_bw()+
  scale_color_manual(values = c("#F7931E","#0000FF"))+
  xlab("Amino Acid")+
  ylab("Value")+
  #labs(fill = "Frequency")+
  theme(legend.text = element_text(colour="black", size=10, 
                                   face="bold"),
        legend.position = c(0.15, 0.87))+
  ggtitle("Comparision of C-Terminal amino acid \n frequency versus expected frequency for 12x")



amino_acid_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  amino_acid_frequency[,i]<-
    table(
      substr(dataset_for_frequency$raw_counts_translation,i,i)
    )
}

amino_acid_frequency$AminoAcid<-c( "A" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,"N", "P" ,"Q", "R", "S" ,"T" ,"V" ,"W" ,"Y")
amino_acid_frequency$Full_Name<-NULL
amino_acid_frequency$Codon_Number<-NULL
for (i in 2:12) {
  amino_acid_frequency[,i]<-as.numeric(amino_acid_frequency[,i])
}
amino_acid_frequency<-melt(amino_acid_frequency, id = "AminoAcid")
amino_acid_frequency$variable<-as.character(amino_acid_frequency$variable)
amino_acid_frequency$position<-substr(amino_acid_frequency$variable,2,(nchar(amino_acid_frequency$variable)) )

amino_acid_frequency<-amino_acid_frequency[,c("AminoAcid","value","position")]

Codon_data<-Codon[,c("AminoAcid","expectedFrequency")]
Codon_data$Position<-"EF"
names(Codon_data)<-c("AminoAcid","Frequency","Position")
names(amino_acid_frequency)<-c("AminoAcid","Frequency","Position")

amino_acid_frequency$Position<-as.numeric(amino_acid_frequency$Position)-13
amino_acid_frequency<-as.data.frame(rbind(amino_acid_frequency,Codon_data))
amino_acid_frequency$Position<-factor(amino_acid_frequency$Position, 
                                      levels = c("EF","-12","-11","-10","-9","-8","-7","-6","-5","-4","-3","-2","-1"))

pdf(paste0(path_output,"FigureS1/Frequency_ExpectedFrequency12x.pdf"))
ggplot(amino_acid_frequency)+
  geom_bar(
    aes(x=Position , 
        y=Frequency  , 
        fill=AminoAcid, 
        label = AminoAcid),
    stat="identity")+
  scale_fill_manual(values=c("#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2",
                             "#808080",
                             "#e2e2e2"
                             
                             
  ))+
  ggtitle("Frequencies and Expected Frequency - All Positions - 12x")+
  theme_bw()
dev.off()


###

dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]

dataset_for_frequency<-unique(dataset_for_frequency)

all_properties<-as.data.frame(cbind(
  hydrophobicity(dataset_for_frequency$raw_counts_translation,scale = "KyteDoolittle"),
  charge(dataset_for_frequency$raw_counts_translation, pH = 7, pKscale = "Lehninger"),
  mw(dataset_for_frequency$raw_counts_translation),
  mz(dataset_for_frequency$raw_counts_translation),
  aIndex(dataset_for_frequency$raw_counts_translation),
  pI(dataset_for_frequency$raw_counts_translation),
  boman(dataset_for_frequency$raw_counts_translation)
))
a<-crucianiProperties(dataset_for_frequency$raw_counts_translation)
for (i in 1:nrow(dataset_for_frequency)) {
  all_properties[i,8]<-a[[i]][1]
  all_properties[i,9]<-a[[i]][3]
}
a<-fasgaiVectors(dataset_for_frequency$raw_counts_translation)
for (j in 1:nrow(dataset_for_frequency)) {
  all_properties[j,10]<-a[[j]][2]
  all_properties[j,11]<-a[[j]][3]
  all_properties[j,12]<-a[[j]][4]
  all_properties[j,13]<-a[[j]][5]
  all_properties[j,14]<-a[[j]][6]
  
}
a<-kideraFactors(dataset_for_frequency$raw_counts_translation)
for (j in 1:nrow(dataset_for_frequency)) {
  all_properties[j,15]<-a[[j]][1]
  all_properties[j,16]<-a[[j]][2]
  all_properties[j,17]<-a[[j]][3]
  all_properties[j,18]<-a[[j]][5]
  all_properties[j,19]<-a[[j]][6]
  all_properties[j,20]<-a[[j]][7]
  all_properties[j,21]<-a[[j]][8]
  all_properties[j,22]<-a[[j]][9]
  all_properties[j,23]<-a[[j]][10]
}
names(all_properties)<-c("hydrophobicity_overall","charge_overall","molecularWeight_overall","massOverCharge_overall",
                         "aliphaticIndex_overall","isoelectricPoint_overall",
                         "potentialProteinInteractionIndex_overall","polarity_overall",
                         "hydrogenBonds_overall","alphaAndTurnPropensity_overall",
                         "bulkyProperties_overall","compositionalCharactersticIndex_overall",
                         "localFlexibility_overall","electronicProperties_overall",
                         "helixBendPreference_overall","sideChainSize_overall",
                         "extendedStructurePreference_overall","doubleBendPreference_overall",
                         "partialSpecificVolume_overall","flatExtendedPreference_overall",
                         "occuranceAlphaRegion_overall","pkc_overall","surroundingHydrophobicity_overall"
                         
)

all_properties$PSI<-dataset_for_frequency$byaa_pooled_PSI

correlation<-as.data.frame(matrix(nrow = 23))
for (i in 1:23) {
  correlation[i,]<-round(cor(all_properties[,i], all_properties[,24]),3)
  
}
correlation$property<-names(all_properties)[1:23]
correlation$absoluteCorrelation<-abs(correlation$V1)

correlation$Choosen<-ifelse(correlation$property %in%  c("hydrophobicity_overall","massOverCharge_overall","aliphaticIndex_overall",
                                                         "potentialProteinInteractionIndex_overall","alphaAndTurnPropensity_overall","bulkyProperties_overall"),
                            "for_model",
                            "others")
correlation$Property<-str_split_fixed(correlation$property,"_",2)[,1]
correlation$absoluteCorrelation<-ifelse((correlation$Property %in% c("extendedStructurePreference","helixBendPreference")),
                                        NA,
                                        correlation$absoluteCorrelation)
correlation<-correlation[complete.cases(correlation),]
correlation$Choosen<-factor(correlation$Choosen, levels = c("others","for_model"))
correlation<-correlation[order(correlation$Choosen),]
correlation$Property<-factor(correlation$Property, levels = c( "bulkyProperties",                  "alphaAndTurnPropensity",           "potentialProteinInteractionIndex" ,"aliphaticIndex" ,                 
                                                               "massOverCharge",                   "hydrophobicity"        ,           "surroundingHydrophobicity"        ,"pkc"             ,                
                                                               "occuranceAlphaRegion",             "flatExtendedPreference" ,          "partialSpecificVolume"            ,"doubleBendPreference",            
                                                               "sideChainSize",                    "electronicProperties"    ,         'localFlexibility'                 ,"compositionalCharactersticIndex", 
                                                               "hydrogenBonds" ,                   "polarity"                 ,        "isoelectricPoint"                 ,"molecularWeight"               ,  
                                                               "charge"))




pdf(paste0(path_plot,"FigureS1/BiophysicalProperty_comparision.pdf"))
correlation %>%
  ggplot(aes(
    x = (Property),
    y = V1,
    fill = Choosen
  ))+
  geom_bar(
    stat=  "identity",
    position = "dodge"
  )+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  xlab("Property")+
  ylab("Absolute Correlation")+
  ggtitle("Biophysical property correlation with PSI")+
  scale_fill_manual(values = c( "for_model" = "#ff5400",
                                "others" = "#8d99ae"))

dev.off()


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


instable_shap_position<-instable_shap[,1:240]
instable_shap_positionBiophysical<-instable_shap[,241:312]
instable_shap_Biophysical<-instable_shap[,313:318]

instable_shap_position<-melt(instable_shap_position)
instable_shap_position$Group<-"Position"
instable_shap_position$Group1<-"Position"

instable_shap_positionBiophysical<-melt(instable_shap_positionBiophysical)
instable_shap_positionBiophysical$Group<-str_split_fixed(instable_shap_positionBiophysical$variable,"_",4)[,2]
instable_shap_positionBiophysical$Group1<-"PositionalBiophysical"

instable_shap_Biophysical<-melt(instable_shap_Biophysical)
instable_shap_Biophysical$Group<-str_split_fixed(instable_shap_Biophysical$variable,"_",4)[,1]
instable_shap_Biophysical$Group1<-"OverallBiophysical"

data_instable<-as.data.frame(rbind(instable_shap_position, instable_shap_positionBiophysical, instable_shap_Biophysical))


positional<-ggplot(data_instable)+
  geom_boxplot(data = subset(data_instable,data_instable$Group1 == "Position"),
               aes(x = Group, 
                   y = value),
               coef = 1e30, color = "red"
  )+
  ylim(-0.4,0.4)+
  guides(col = FALSE)+
  ylab("Contribution Score")+
  xlab("")+
  ggtitle("Sequence")+
  #scale_x_discrete(guide = guide_axis(angle = 90))+
  theme_bw()

positional
positionalBiophysical<-ggplot(data_instable)+
  geom_boxplot(data = subset(data_instable,data_instable$Group1 == "PositionalBiophysical"),
               aes(x = Group, 
                   y = value),color = "blue",
               coef = 1e30
  )+ 
  theme_bw()+
  theme(#remove x axis labels
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank()  #remove y axis ticks
  )+
  ggtitle("Positional Biophysical")+
  guides(col = FALSE)+
  #scale_x_discrete(guide = guide_axis(angle = 90))+
  ylim(-0.4,0.4)

overallBiophysical<-ggplot(data_instable)+
  geom_boxplot(data = subset(data_instable,data_instable$Group1 == "OverallBiophysical"),
               aes(x = Group, 
                   y = value,),
               coef = 1e30
  )+
  xlab("")+
  theme_bw()+
  theme(#remove x axis labels
    axis.title.y=element_blank(),
    axis.text.y=element_blank(),  #remove y axis labels
    axis.ticks.y=element_blank()  #remove y axis ticks
  )+
  labs(title = "Overall Biophysical")+
  
  
  ylim(-0.4,0.4)+
  guides(col = FALSE)

title <- ggdraw() + 
  draw_label(
    " A : Aliphatic Index, B: Bulkiness, H: Hydrophobicity, M: Mass Over Charge,  \n P: Potential Protein Interaction Index, T:  Alpha And Turn Propensity    
    ",fontface = "plain",size = 10,
    x = 0,
    hjust = 0
  ) +
  theme(
    # add margin on the left of the drawing canvas,
    # so title is aligned with left edge of first plot
    plot.margin = margin(0, 0, 0, 7)
  )

plot<-plot_grid(
  positional,positionalBiophysical,overallBiophysical, ncol = 3)

pdf(paste0(path_plot,"FigureS1/SHAP_data_variation.pdf"))
plot_grid(
  plot,title,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.1)
)
dev.off()


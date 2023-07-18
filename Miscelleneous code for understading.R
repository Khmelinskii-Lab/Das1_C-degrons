
##### for project : Analysis of degron at C terminal
##### sub text:Miscelleneous
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

yeastmine2<-read.csv(paste0(path_input,"yeastmine_3.csv"), sep = ",", header = FALSE)
names(yeastmine2)<-c("GeneID","Gene_Systematic_Name","Organism","Standard_Name","Gene_Name","Gene_Qualifier","Ontology_Term_ID","Ontology_Term_Name","Ontology_Term_Namespace","Protein_Standard_Name","Residue","Protein_Systematic_Name")
yeastmine2$Residue<-substr(yeastmine2$Residue,1,(nchar(yeastmine2$Residue)-1))
yeastmine2$translation<-substr(yeastmine2$Residue,(nchar(yeastmine2$Residue)-11),(nchar(yeastmine2$Residue)))

yeastmine2<-yeastmine2[,c("Protein_Systematic_Name","Gene_Qualifier","translation")]
yeastmine2<-unique(yeastmine2)

WT_Das1<-read.csv(paste0(path_output,"WT_Das1_Cter.csv"))




yeastmine_translation<-as.data.frame(unique(yeastmine2$translation))
names(yeastmine_translation)<-c("translation")

amino_acid_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  amino_acid_frequency[,i]<-
    table(
      substr(yeastmine_translation$translation,i,i)
    )
}

amino_acid_frequency<-as.data.frame(amino_acid_frequency[,12])
amino_acid_frequency$AminoAcid<-c( "A" ,"C" ,"D" ,"E" ,"F" ,"G" ,"H" ,"I" ,"K" ,"L" ,"M" ,"N", "P" ,"Q", "R", "S" ,"T" ,"V" ,"W" ,"Y")
amino_acid_frequency$Var1<-NULL
names(amino_acid_frequency)<-c("FreqAll","AminoAcid")
amino_acid_frequency$FreqAll<-as.numeric(amino_acid_frequency$FreqAll)

amino_acid_frequency$relativeFreqAll<-amino_acid_frequency$FreqAll/6410


Codon<-data.frame(
  Full_Name = c("Isoleucine", "Leucine","Valine","Phenylalanine","Methionine","Cysteine","Alanine","Glycine","Proline","Threonine","Serine","Tyrosine","Tryptophan","Glutamine","Asparagine","Histidine","Glutamic acid","Aspartic acid","Lysine","Arginine"),
  AminoAcid = c("I","L","V","F","M","C","A","G","P","T","S","Y","W","Q","N","H","E","D","K","R"),
  Codon_Number = c(3,6,4,2,1,2,4,4,4,4,6,2,1,2,2,2,2,2,2,6)
)
Codon$expectedFrequency<-(Codon$Codon_Number * 6410)/61

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
  ggtitle("Comparision of C-Terminal amino acid \n frequency versus expected frequency for yeast C termini from yeastmine")



amino_acid_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  amino_acid_frequency[,i]<-
    table(
      substr(yeastmine_translation$raw_counts_translation,i,i)
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

#pdf(paste0(path_plot,"FigureS1/Frequency_ExpectedFrequency12x.pdf"))
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
#dev.off()

# for frequency of unstable
dataset_for_frequency<-WT_Das1[,c("translation","value.x")]
dataset_for_frequency<-unique(dataset_for_frequency)
amino_acid_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  amino_acid_frequency[,i]<-
    table(
      substr(dataset_for_frequency$translation,i,i)
    )
}
amino_acid_frequency$AminoAcid<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
amino_acid_frequency<-melt(amino_acid_frequency, id = "AminoAcid")
amino_acid_frequency$value<-as.numeric(amino_acid_frequency$value)
amino_acid_frequency$relativeFrequencyAll<-amino_acid_frequency$value/nrow(dataset_for_frequency)
amino_acid_frequency$variable<-paste0(amino_acid_frequency$variable,amino_acid_frequency$AminoAcid)
amino_acid_frequency$AminoAcid<-NULL
amino_acid_frequency$value<-NULL




Unstable<-subset(WT_Das1,
                 WT_Das1$class == "Unstable")
Unstable_aa_frequency<-as.data.frame(matrix(nrow = 20))
for (i in 1:12) {
  Unstable_aa_frequency[,i]<-
    table(
      substr(Unstable$translation,i,i)
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

#pdf(paste0(path_plot,"Figure1/NormalizedFrequency.pdf"))
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
  scale_fill_gradientn("log2(Freq_Degron)/\n (Freq_YeastAll)",
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


#dev.off()


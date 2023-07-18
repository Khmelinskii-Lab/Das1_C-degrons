# finding the das1 motifs

##### for project : Analysis of degron at C terminal
##### sub text:cluster of amino acid
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")

data_with_stat_test<-read.csv(paste0(path_output,"data_with_stat_test_wt_das1.csv"))


dataset<-dataset<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))
dataset<-subset(dataset,dataset$Codon == "Non Stop Codon")
names(dataset)
dataset_for_frequency<-dataset[,c("raw_counts_translation","byaa_pooled_PSI")]
dataset_for_frequency<-unique(dataset_for_frequency)
dataset_for_frequency$hydrophobicity<-hydrophobicity(dataset_for_frequency$raw_counts_translation)
dataset_for_frequency_oneHot<-create_one_hot_db(dataset_for_frequency$raw_counts_translation,12)

WT_DT_unstable_data<-create_one_hot_db(data_with_stat_test$translation,12)
amino_acid<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
t<-strrep(amino_acid, 12)  

data_with_stat_test<-read.csv(paste0(path_output,"data_with_stat_test_wt_das1.csv"))


nrow(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized")& (substr(data_with_stat_test$translation,8,8) =="W"))))
nrow(subset(data_with_stat_test,((data_with_stat_test$category == "No Effect")& (substr(data_with_stat_test$translation,8,8) =="W"))))
nrow(subset(dataset,((dataset$byaa_pooled_PSI > 0.5)& (substr(dataset$raw_counts_translation,8,8) =="W"))))


stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
stabilized<-subset(stabilized,substr(stabilized$translation,8,8) =="W")
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
#stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
#stabilized$amino_acid<-factor(stabilized$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"Figure5/ClusterWMinus5.pdf"))
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
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
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    
    ggseqlogo(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                            (substr(data_with_stat_test$translation,8,8) =="W")))[,c("translation")])+
      theme(axis.ticks.x = element_line(size =1), axis.ticks.y = element_line(size =1))+
      annotate(geom = "text", 
               x= 2,
               y=3,
               label = paste0("#peptides",
                              nrow((subset(dataset_for_frequency,
                                           dataset_for_frequency$raw_counts_translation %in% 
                                             subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                           (substr(data_with_stat_test$translation,8,8) =="W")))[,c("translation")]))))),
    
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
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
    
    ggplot(dataset_for_frequency)+
      geom_point(data = dataset_for_frequency, aes(
        x = byaa_pooled_PSI,
        y = hydrophobicity
      ), color = "grey")+
      geom_point(data = subset(dataset_for_frequency,
                               dataset_for_frequency$raw_counts_translation %in% 
                                 subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                               (substr(data_with_stat_test$translation,8,8) =="W")))[,c("translation")]), 
                 aes(
                   x = byaa_pooled_PSI,
                   y = hydrophobicity
                 ), color = "red")+
      ylim(-4,4)+
      xlim(0,1)+
      ylab("Hydrophobicity")+
      xlab("Protein Stability Index")+
      theme_bw()+
      annotate(geom = "text", 
               x= 0.1,
               y=3.8,
               label = paste0("#peptides",
                 nrow((subset(dataset_for_frequency,
                                                 dataset_for_frequency$raw_counts_translation %in% 
                                                   subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                                 (substr(data_with_stat_test$translation,8,8) =="W")))[,c("translation")])))))
    
   
    
    
    
 
  )
)

dev.off()
#stabilize
d1<-stabilized
#stabilized1$log2NormFreqAll<-ifelse((stabilized1$log2NormFreqAll< 1 & stabilized1$log2NormFreqAll> -1),0,stabilized1$log2NormFreqAll)

#data_seqlogo<-data.frame(matrix(t(stabilized1[1:240,10]), nrow=20, byrow=FALSE))
#names(data_seqlogo)<-as.character(-12:-1)
#row.names(data_seqlogo)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#names(data_seqlogo)<-as.numeric(names(data_seqlogo))-13
#data_seqlogo[data_seqlogo<2.5,]<-0

#a<- print(ggseqlogo(as.matrix(data_seqlogo), method='custom', seq_type='aa') +
 #           theme(axis.text.x = element_blank(),
  #                axis.ticks.y = element_line(size =1))+
   #         ylab('log2 normalized to 12x')+
    #        xlab("Position")+
     #       annotate(geom = "text",
      #               x=.5,
       #              y = max(data_seqlogo)+0.05, 
        #             label = paste0("Motif : ", 1, " , \n # peptides : ",count))+
         #   ylim(min(data_seqlogo)-0.05, max(data_seqlogo)+0.1)+
          #  annotate(geom = 'segment', y = 0, yend = 0, color = 'black', x = 1, xend = Inf, size = 0.25) )


###############

nrow(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized")& (substr(data_with_stat_test$translation,10,10) =="K"))))
nrow(subset(data_with_stat_test,((data_with_stat_test$category == "No Effect")& (substr(data_with_stat_test$translation,10,10) =="K"))))
nrow(subset(dataset,((dataset$byaa_pooled_PSI > 0.5)& (substr(dataset$raw_counts_translation,10,10) =="K"))))

stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
stabilized<-subset(stabilized,substr(stabilized$translation,10,10) =="K")
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
#stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
#stabilized$amino_acid<-factor(stabilized$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"Figure5/cluster_KMinus3.pdf"))
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
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
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    
    ggseqlogo(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                            (substr(data_with_stat_test$translation,10,10) =="K")))[,c("translation")])+
      theme(axis.ticks.x = element_line(size =1), axis.ticks.y = element_line(size =1))+ 
      annotate(geom = "text", 
                                                                                                  x= 2,
                                                                                                  y=3,
                                                                                                  label = paste0("#peptides",
                                                                                                                 nrow((subset(dataset_for_frequency,
                                                                                                                              dataset_for_frequency$raw_counts_translation %in% 
                                                                                                                                subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                                                                                                              (substr(data_with_stat_test$translation,10,10) =="K")))[,c("translation")]))))),
    
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
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
    
    ggplot(dataset_for_frequency)+
      geom_point(data = dataset_for_frequency, aes(
        x = byaa_pooled_PSI,
        y = hydrophobicity
      ), color = "grey")+
      geom_point(data = subset(dataset_for_frequency,
                               dataset_for_frequency$raw_counts_translation %in% subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                                                               (substr(data_with_stat_test$translation,10,10) =="K")))[,c("translation")]), 
                 aes(
                   x = byaa_pooled_PSI,
                   y = hydrophobicity
                 ), color = "red")+
      ylim(-4,4)+
      xlim(0,1)+
      ylab("Hydrophobicity")+
      xlab("Protein Stability Index")+
      theme_bw()+
      
      annotate(geom = "text", 
               x= 0.1,
               y=3.8,
               label = paste0("#peptides",
                              nrow((subset(dataset_for_frequency,
                                           dataset_for_frequency$raw_counts_translation %in% 
                                             subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                           (substr(data_with_stat_test$translation,10,10) =="K")))[,c("translation")])))))
    
    
    
    
    
    
    
  )
)

dev.off()
#stabilized1<-stabilized
#stabilized1$log2NormFreqAll<-ifelse((stabilized1$log2NormFreqAll< 1.5 & stabilized1$log2NormFreqAll> -1.5),0,stabilized1$log2NormFreqAll)

#data_seqlogo<-data.frame(matrix(t(stabilized1[1:240,10]), nrow=20, byrow=FALSE))
#names(data_seqlogo)<-as.character(-12:-1)
#row.names(data_seqlogo)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#names(data_seqlogo)<-as.numeric(names(data_seqlogo))-13
#data_seqlogo[data_seqlogo<2.5,]<-0

#a<- print(ggseqlogo(as.matrix(data_seqlogo), method='custom', seq_type='aa') +
 #           theme(axis.text.x = element_blank(),
  #                axis.ticks.y = element_line(size =1))+
   #         ylab('log2 normalized to 12x')+
    #        xlab("Position")+
     #       annotate(geom = "text",
      #               x=.5,
       #              y = max(data_seqlogo)+0.05, 
        #             label = paste0("Motif : ", 1, " , \n # peptides : ",count))+
         #   ylim(min(data_seqlogo)-0.05, max(data_seqlogo)+0.1)+
          #  annotate(geom = 'segment', y = 0, yend = 0, color = 'black', x = 1, xend = Inf, size = 0.25) )






###############

nrow(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized")& (substr(data_with_stat_test$translation,11,11) =="I"))))
nrow(subset(data_with_stat_test,((data_with_stat_test$category == "No Effect")& (substr(data_with_stat_test$translation,11,11) =="I"))))
nrow(subset(dataset,((dataset$byaa_pooled_PSI > 0.5)& (substr(dataset$raw_counts_translation,11,11) =="I"))))

stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
stabilized<-subset(stabilized,substr(stabilized$translation,11,11) =="I")
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
#stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
#stabilized$amino_acid<-factor(stabilized$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"Figure5/Cluster_IMinus2.pdf"))
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
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
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    
    ggseqlogo(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                            (substr(data_with_stat_test$translation,11,11) =="I")))[,c("translation")])+
      theme(axis.ticks.x = element_line(size =1), axis.ticks.y = element_line(size =1))+
      annotate(geom = "text", 
               x= 2,
               y=3,
               label = paste0("#peptides",
                              nrow((subset(dataset_for_frequency,
                                           dataset_for_frequency$raw_counts_translation %in% 
                                             subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                           (substr(data_with_stat_test$translation,11,11) =="I")))[,c("translation")]))))),
    
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
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
    
    ggplot(dataset_for_frequency)+
      geom_point(data = dataset_for_frequency, aes(
        x = byaa_pooled_PSI,
        y = hydrophobicity
      ), color = "grey")+
      geom_point(data = subset(dataset_for_frequency,
                               dataset_for_frequency$raw_counts_translation %in% subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                                                               (substr(data_with_stat_test$translation,11,11) =="I")))[,c("translation")]), 
                 aes(
                   x = byaa_pooled_PSI,
                   y = hydrophobicity
                 ), color = "red")+
      ylim(-4,4)+
      xlim(0,1)+
      ylab("Hydrophobicity")+
      xlab("Protein Stability Index")+
      theme_bw()+
      annotate(geom = "text", 
               x= 0.1,
               y=3.8,
               label = paste0("#peptides",
                              nrow((subset(dataset_for_frequency,
                                           dataset_for_frequency$raw_counts_translation %in% 
                                             subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                           (substr(data_with_stat_test$translation,11,11) =="I")))[,c("translation")])))))
    
    
    
    
    
    
    
    
  )
)

dev.off()
#stabilized1<-stabilized
#stabilized1$log2NormFreqAll<-ifelse((stabilized1$log2NormFreqAll< 1 & stabilized1$log2NormFreqAll> -1),0,stabilized1$log2NormFreqAll)

#data_seqlogo<-data.frame(matrix(t(stabilized1[1:240,10]), nrow=20, byrow=FALSE))
#names(data_seqlogo)<-as.character(-12:-1)
#row.names(data_seqlogo)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#names(data_seqlogo)<-as.numeric(names(data_seqlogo))-13

#a<- print(ggseqlogo(as.matrix(data_seqlogo), method='custom', seq_type='aa') +
 #           theme(axis.text.x = element_blank(),
  #                axis.ticks.y = element_line(size =1))+
   #         ylab('log2 normalized to 12x')+
    #        xlab("Position")+
     #       annotate(geom = "text",
      #               x=.5,
       #              y = max(data_seqlogo)+0.05, 
        #             label = paste0("Motif : ", 1, " , \n # peptides : ",count))+
         #   ylim(min(data_seqlogo)-0.05, max(data_seqlogo)+0.1)+
          #  annotate(geom = 'segment', y = 0, yend = 0, color = 'black', x = 1, xend = Inf, size = 0.25) )


###############

nrow(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized")& (substr(data_with_stat_test$translation,12,12) =="N"))))
nrow(subset(data_with_stat_test,((data_with_stat_test$category == "No Effect")& (substr(data_with_stat_test$translation,12,12) =="N"))))
nrow(subset(dataset,((dataset$byaa_pooled_PSI > 0.5)& (substr(dataset$raw_counts_translation,12,12) =="N"))))


stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
stabilized<-subset(stabilized,substr(stabilized$translation,12,12) =="N")
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
#stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
#stabilized$amino_acid<-factor(stabilized$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

pdf(paste0(path_plot,"Figure5/Cluster_NMinus1.pdf"))
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
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
                           values=scales::rescale(c(-2,-1.6,-1.3,-1,-0.6,-0.3,0,0.3,0.61,1.3,1.6,2)),
                           limits=c(-2,2)),
    
    
    ggseqlogo(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                            (substr(data_with_stat_test$translation,12,12) =="N")))[,c("translation")])+
      theme(axis.ticks.x = element_line(size =1), axis.ticks.y = element_line(size =1))+
      annotate(geom = "text", 
               x= 2,
               y=3,
               label = paste0("#peptides:",
                              nrow((subset(dataset_for_frequency,
                                           dataset_for_frequency$raw_counts_translation %in% 
                                             subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                           (substr(data_with_stat_test$translation,12,12) =="N")))[,c("translation")]))))),
    
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
      scale_x_continuous(breaks=seq(-12, -1, 1))+
      scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
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
    
    ggplot(dataset_for_frequency)+
      geom_point(data = dataset_for_frequency, aes(
        x = byaa_pooled_PSI,
        y = hydrophobicity
      ), color = "grey")+
      geom_point(data = subset(dataset_for_frequency,
                               dataset_for_frequency$raw_counts_translation %in% subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                                                               (substr(data_with_stat_test$translation,12,12) =="N")))[,c("translation")]), 
                 aes(
                   x = byaa_pooled_PSI,
                   y = hydrophobicity
                 ), color = "red")+
      ylim(-4,4)+
      xlim(0,1)+
      ylab("Hydrophobicity")+
      xlab("Protein Stability Index")+
      theme_bw()+
      annotate(geom = "text", 
               x= 0.1,
               y=3.8,
               label = paste0("#peptides:",
                              nrow((subset(dataset_for_frequency,
                                           dataset_for_frequency$raw_counts_translation %in% 
                                             subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                           (substr(data_with_stat_test$translation,12,12) =="N")))[,c("translation")])))))
    
    
    
    
    
    
    
    
  )
)

dev.off()
#stabilized1<-stabilized
#stabilized1$log2NormFreqAll<-ifelse((stabilized1$log2NormFreqAll< 1.5 & stabilized1$log2NormFreqAll> -1.5),0,stabilized1$log2NormFreqAll)

#data_seqlogo<-data.frame(matrix(t(stabilized1[1:240,10]), nrow=20, byrow=FALSE))
#names(data_seqlogo)<-as.character(-12:-1)
#row.names(data_seqlogo)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#names(data_seqlogo)<-as.numeric(names(data_seqlogo))-13

#a<- print(ggseqlogo(as.matrix(data_seqlogo), method='custom', seq_type='aa') +
 #           theme(axis.text.x = element_blank(),
  #                axis.ticks.y = element_line(size =1))+
   #         ylab('log2 normalized to 12x')+
    #        xlab("Position")+
     #       annotate(geom = "text",
      #               x=.5,
       #              y = max(data_seqlogo)+0.05, 
        #             label = paste0("Motif : ", 1, " , \n # peptides : ",count))+
         #   ylim(min(data_seqlogo)-0.05, max(data_seqlogo)+0.1)+
          #  annotate(geom = 'segment', y = 0, yend = 0, color = 'black', x = 1, xend = Inf, size = 0.25) )







###############
nrow(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized")& (substr(data_with_stat_test$translation,12,12) =="C"))))
nrow(subset(data_with_stat_test,((data_with_stat_test$category == "No Effect")& (substr(data_with_stat_test$translation,12,12) =="C"))))
nrow(subset(dataset,((dataset$byaa_pooled_PSI > 0.5)& (substr(dataset$raw_counts_translation,12,12) =="C"))))


stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
stabilized<-subset(stabilized,substr(stabilized$translation,12,12) =="C")
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
#stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
#stabilized$amino_acid<-factor(stabilized$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

stabilized$log2NormFreq<-ifelse(stabilized$log2NormFreq < -2,-2,stabilized$log2NormFreq)
stabilized$log2NormFreq<-ifelse(stabilized$log2NormFreq > 2,2,stabilized$log2NormFreq)

stabilized$log2NormFreqAll<-ifelse(stabilized$log2NormFreqAll < -2,-2,stabilized$log2NormFreqAll)
stabilized$log2NormFreqAll<-ifelse(stabilized$log2NormFreqAll > 2,2,stabilized$log2NormFreqAll)

pdf(paste0(path_plot,"Figure5/Cluster_CMinus1.pdf"))
print(plot_grid(
  ggplot(stabilized, aes (x = position,
                          y= amino_acid,
                          fill = log2NormFreq))+
    geom_tile(data = stabilized, 
              aes(
                x = position,
                y= amino_acid,
                fill = log2NormFreq
              ))+
    coord_equal()+
    theme_bw()+
    theme(text = element_text(size = 8))+
    scale_x_continuous(breaks=seq(-12, -1, 1))+
    scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
    geom_text(aes(label = ifelse(is.na(stabilized$log2NormFreq), "*", ""),), 
              color = "Black", 
              size = 2)+
    scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_Unstable)",
                         na.value = "#F7931E",
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
  
  
  
  ggseqlogo(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                          (substr(data_with_stat_test$translation,12,12) =="C")))[,c("translation")])+
    theme(axis.ticks.x = element_line(size =1), axis.ticks.y = element_line(size =1))+
    annotate(geom = "text", 
             x= 2,
             y=3,
             label = paste0("#peptides:",
                            nrow((subset(dataset_for_frequency,
                                         dataset_for_frequency$raw_counts_translation %in% 
                                           subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                         (substr(data_with_stat_test$translation,12,12) =="C")))[,c("translation")])))))
  ,
  
  
  
  
  ggplot(stabilized, aes (x = position,
                          y= amino_acid,
                          fill = log2NormFreqAll))+
    geom_tile(data = stabilized, 
              aes(
                x = position,
                y= amino_acid,
                fill = log2NormFreqAll
              ))+
    coord_equal()+
    theme_bw()+
    theme(text = element_text(size = 8))+
    scale_x_continuous(breaks=seq(-12, -1, 1))+
    scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
    geom_text(aes(label = ifelse(is.na(stabilized$log2NormFreqAll), "*", ""),), 
              color = "Black", 
              size = 2)+
    scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_12x)",
                         na.value = "#F7931E",
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
  
  
  
  ggplot(dataset_for_frequency)+
    geom_point(data = dataset_for_frequency, aes(
      x = byaa_pooled_PSI,
      y = hydrophobicity
    ), color = "grey")+
    geom_point(data = subset(dataset_for_frequency,
                             dataset_for_frequency$raw_counts_translation %in% subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                                                             (substr(data_with_stat_test$translation,12,12) =="C")))[,c("translation")]), 
               aes(
                 x = byaa_pooled_PSI,
                 y = hydrophobicity
               ), color = "red")+
    ylim(-4,4)+
    xlim(0,1)+
    ylab("Hydrophobicity")+
    xlab("Protein Stability Index")+
    theme_bw()+
    annotate(geom = "text", 
             x= 0.1,
             y=3.8,
             label = paste0("#peptides:",
                            nrow((subset(dataset_for_frequency,
                                         dataset_for_frequency$raw_counts_translation %in% 
                                           subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                         (substr(data_with_stat_test$translation,12,12) =="C")))[,c("translation")])))))
  
))


dev.off()
#stabilized1<-stabilized
#stabilized1$log2NormFreqAll<-ifelse((stabilized1$log2NormFreqAll< 1 & stabilized1$log2NormFreqAll> -1),0,stabilized1$log2NormFreqAll)

#data_seqlogo<-data.frame(matrix(t(stabilized1[1:240,10]), nrow=20, byrow=FALSE))
#names(data_seqlogo)<-as.character(-12:-1)
#row.names(data_seqlogo)<-c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
#names(data_seqlogo)<-as.numeric(names(data_seqlogo))-13

#a<- print(ggseqlogo(as.matrix(data_seqlogo), method='custom', seq_type='aa') +
#            theme(axis.text.x = element_blank(),
 #                 axis.ticks.y = element_line(size =1))+
  #          ylab('log2 normalized to 12x')+
   #         xlab("Position")+
    #        annotate(geom = "text",
     #                x=.5,
      #               y = max(data_seqlogo)+0.05, 
       #              label = paste0("Motif : ", 1, " , \n # peptides : ",count))+
        #    ylim(min(data_seqlogo)-0.05, max(data_seqlogo)+0.1)+
         #   annotate(geom = 'segment', y = 0, yend = 0, color = 'black', x = 1, xend = Inf, size = 0.25) )
#






stabilized<-subset(data_with_stat_test,data_with_stat_test$category == "Stabilized")
stabilized<-subset(stabilized,substr(stabilized$translation,11,11) %in% c("I","L","M"))
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
#stabilized$position<-factor(stabilized$position, levels = c(-12,-11,-10,-9,-8,-7,-6,-5,-4,-3,-2,-1))
#stabilized$amino_acid<-factor(stabilized$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
stabilized[stabilized$log2NormFreq == "-Inf",c("log2NormFreq")]<-NA


stabilized$freqAllData<-as.data.frame(colSums(dataset_for_frequency_oneHot))[,1]/46152
stabilized$log2NormFreqAll<-log2(stabilized$frequency/stabilized$freqAllData)
stabilized$NormFreqAll<-(stabilized$frequency/stabilized$freqAllData)

stabilized[stabilized$log2NormFreqAll == "-Inf",c("log2NormFreqAll")]<-NA

stabilized$log2NormFreq<-ifelse(stabilized$log2NormFreq < -2,-2,stabilized$log2NormFreq)
stabilized$log2NormFreq<-ifelse(stabilized$log2NormFreq > 2,2,stabilized$log2NormFreq)

stabilized$log2NormFreqAll<-ifelse(stabilized$log2NormFreqAll < -2,-2,stabilized$log2NormFreqAll)
stabilized$log2NormFreqAll<-ifelse(stabilized$log2NormFreqAll > 2,2,stabilized$log2NormFreqAll)

pdf(paste0(path_plot,"Figure5/Cluster_HydroMinus2.pdf"))
print(plot_grid(
  ggplot(stabilized, aes (x = position,
                          y= amino_acid,
                          fill = log2NormFreq))+
    geom_tile(data = stabilized, 
              aes(
                x = position,
                y= amino_acid,
                fill = log2NormFreq
              ))+
    coord_equal()+
    theme_bw()+
    theme(text = element_text(size = 8))+
    scale_x_continuous(breaks=seq(-12, -1, 1))+
    scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
    geom_text(aes(label = ifelse(is.na(stabilized$log2NormFreq), "*", ""),), 
              color = "Black", 
              size = 2)+
    scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_Unstable)",
                         na.value = "#F7931E",
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
  
  
  
  ggseqlogo(subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                          (substr(data_with_stat_test$translation,11,11) %in% c("I","L","M"))))[,c("translation")])+
    
    theme(axis.ticks.x = element_line(size =1), axis.ticks.y = element_line(size =1))+
    annotate(geom = "text", 
             x= 2,
             y=3,
             label = paste0("#peptides:",
                            nrow((subset(dataset_for_frequency,
                                         dataset_for_frequency$raw_counts_translation %in% 
                                           subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                         (substr(data_with_stat_test$translation,11,11) %in% c("I","L","M"))))[,c("translation")])))))
  ,
  
  
  
  
  ggplot(stabilized, aes (x = position,
                          y= amino_acid,
                          fill = log2NormFreqAll))+
    geom_tile(data = stabilized, 
              aes(
                x = position,
                y= amino_acid,
                fill = log2NormFreqAll
              ))+
    coord_equal()+
    theme_bw()+
    theme(text = element_text(size = 8))+
    scale_x_continuous(breaks=seq(-12, -1, 1))+
    scale_y_discrete(limits = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))+
    geom_text(aes(label = ifelse(is.na(stabilized$log2NormFreqAll), "*", ""),), 
              color = "Black", 
              size = 2)+
    scale_fill_gradientn("log2(Freq_stabilized)/\n (Freq_12x)",
                         na.value = "#F7931E",
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
  
  
  
  ggplot(dataset_for_frequency)+
    geom_point(data = dataset_for_frequency, aes(
      x = byaa_pooled_PSI,
      y = hydrophobicity
    ), color = "grey")+
    geom_point(data = subset(dataset_for_frequency,
                             dataset_for_frequency$raw_counts_translation %in% subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                                                             (substr(data_with_stat_test$translation,11,11) %in% c("I","L","M"))))[,c("translation")]), 
               aes(
                 x = byaa_pooled_PSI,
                 y = hydrophobicity
               ), color = "red")+
    ylim(-4,4)+
    xlim(0,1)+
    ylab("Hydrophobicity")+
    xlab("Protein Stability Index")+
    theme_bw()+

    annotate(geom = "text", 
             x= 0.1,
             y=3.8,
             label = paste0("#peptides:",
                            nrow((subset(dataset_for_frequency,
                                         dataset_for_frequency$raw_counts_translation %in% 
                                           subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                                         (substr(data_with_stat_test$translation,11,11) %in% c("I","L","M"))))[,c("translation")])))))
  
))

dev.off()


X<-list(WMinus5 = subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                (substr(data_with_stat_test$translation,8,8) =="W")))[,c("translation")],
        HydroMinus2 = subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                    (substr(data_with_stat_test$translation,11,11) %in% c("I","L","M"))))[,c("translation")],
        
        KMinus3 = subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                (substr(data_with_stat_test$translation,10,10) =="K")))[,c("translation")],
        NMinus1 = subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                (substr(data_with_stat_test$translation,12,12) =="N")))[,c("translation")],
        CMinus1 = subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                (substr(data_with_stat_test$translation,12,12) =="C")))[,c("translation")],
        IMinus2 = subset(data_with_stat_test,((data_with_stat_test$category == "Stabilized") & 
                                                (substr(data_with_stat_test$translation,11,11) =="I")))[,c("translation")]
)

pdf(paste0(path_plot,"FigureS5/Venn.pdf"))
venn(X[1:5], ilab=TRUE, ilcs = 0.8,opacity = 0.5,zcolor = "style",
     plotsize = 10,ellipse = TRUE,ilabels = TRUE,sncs = 0.7,
     #snames = c("W-1","K-3","I-2","N-1",'C-1'),
     borders = TRUE, box = TRUE, par = TRUE,ggplot = TRUE)

dev.off()
#

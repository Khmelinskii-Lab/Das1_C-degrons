# merging replicate 1

##### for project : Analysis of degron at C terminal
##### sub text:merging all raw data into desired form and creating dataset for deep learning analysis
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB


rm(list = ls())

# set the directory
setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")

# reading all the data


CDeg_12_counts<-read.csv(paste0(path_input,"12_X/raw_counts_Cdeg_12_rep1.txt"), 
                         header = TRUE, 
                         sep = "\t")
CDeg_12_byaa<-read.csv(paste0(path_input,"12_X/raw_PSI_byaa_Cdeg_12_rep1.txt"), 
                       header = TRUE, 
                       sep = "\t")
CDeg_12_bynuc<-read.csv(paste0(path_input,"12_X/raw_PSI_bynuc_Cdeg_12_rep1.txt"), 
                        header = TRUE, 
                        sep = "\t")

#########################################################################################################################################
####################################################integrating data for 12 amino acid long sequence#####################################################################################
#########################################################################################################################################


names(CDeg_12_counts)<-paste0("raw_counts_",
                              names(CDeg_12_counts))
names(CDeg_12_bynuc)<-paste0("bynuc_",
                             names(CDeg_12_bynuc))
names(CDeg_12_byaa)<-paste0("byaa_",
                            names(CDeg_12_byaa))

#########################################################################################################################################
CDeg_12_counts<-CDeg_12_counts[,c("raw_counts_experiment",
                                  "raw_counts_sub_experiment",
                                  "raw_counts_bin",
                                  "raw_counts_fraction",
                                  "raw_counts_f",
                                  "raw_counts_sequence",
                                  "raw_counts_sum",
                                  "raw_counts_translation",
                                  "raw_counts_count",
                                  "raw_counts_normalized_counts")]
CDeg_12_bynuc<-CDeg_12_bynuc[,c("bynuc_experiment",
                                "bynuc_sub_experiment",
                                "bynuc_sequence",
                                "bynuc_translation",
                                "bynuc_PSI")]
CDeg_12_byaa<-CDeg_12_byaa[,c("byaa_experiment",
                              "byaa_sub_experiment",
                              "byaa_translation",
                              "byaa_pooled_PSI")]

CDeg_12_counts$identifier<-paste0(CDeg_12_counts$raw_counts_experiment,"_",
                                  CDeg_12_counts$raw_counts_sub_experiment ,"_",
                                  CDeg_12_counts$raw_counts_sequence,"_",
                                  CDeg_12_counts$raw_counts_translation)

CDeg_12_bynuc$identifier<-paste0(CDeg_12_bynuc$bynuc_experiment,"_",
                                 CDeg_12_bynuc$bynuc_sub_experiment ,"_",
                                 CDeg_12_bynuc$bynuc_sequence,"_",
                                 CDeg_12_bynuc$bynuc_translation)

merged_file<-merge(CDeg_12_counts,
                   CDeg_12_bynuc,
                   by= "identifier")

merged_file$identifier<-paste0(merged_file$raw_counts_experiment,"_",
                               merged_file$raw_counts_sub_experiment,"_",
                               merged_file$bynuc_translation)

CDeg_12_byaa$identifier<-paste0(CDeg_12_byaa$byaa_experiment,"_",
                                CDeg_12_byaa$byaa_sub_experiment ,"_",
                                CDeg_12_byaa$byaa_translation)

merged_file<-merge(merged_file,
                   CDeg_12_byaa,
                   by="identifier")

merged_file$raw_counts_translation<-as.factor(unlist(merged_file$raw_counts_translation))
merged_file$bynuc_translation<-as.factor(unlist(merged_file$bynuc_translation))
merged_file$byaa_translation<-as.factor(unlist(merged_file$byaa_translation))

merged_file$Codon<-ifelse(grepl("\\*",merged_file$raw_counts_translation),"Stop Codon","Non Stop Codon")
#some manual checks
merged_file$group_raw_PSI<-cut(merged_file$bynuc_PSI,
                               c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                               c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))

merged_file$group_raw_pooled_PSI<-cut(merged_file$byaa_pooled_PSI,
                                      c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                                      c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))

q<-merged_file %>% 
  group_by(raw_counts_translation) %>% 
  summarise(sum = sum(raw_counts_count))

merged_file<-merge(merged_file,
                   q, 
                   id="raw_counts_translation")

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
# removing translation with in total smaller reads

pdf(paste0(path_plot,"miscelleneous/raw_histogram_for_threshold.pdf"))
print((histogram_for_raw_count(merged_file,
                               merged_file$raw_counts_count,
                               merged_file$Codon)+
         scale_x_continuous(trans = 'log10')+
         labs(title = "histogram plot for raw count\n", 
              x = "Log of raw read counts", 
              y = " Frequency", 
              color = "Codon class\n")+
         theme_bw()+
         theme(axis.title = element_text(face="bold"))))
dev.off()

pdf(paste0(path_plot,"miscelleneous/raw_histogram_for_threshold_normalized.pdf"))
print((histogram_for_raw_count(merged_file,
                               merged_file$raw_counts_normalized_counts,
                               merged_file$Codon)+
         scale_x_continuous(trans = 'log10')+
         labs(title = "histogram plot for normalized raw count\n", 
              x = "Log of normalized read counts", 
              y = " Frequency", 
              color = "Codon class\n")+
         theme_bw()+
         theme(axis.title = element_text(face="bold"))))
dev.off()

percentage_of_stop_codon<-sum(subset(merged_file,merged_file$Codon == "Stop Codon")[,c("raw_counts_normalized_counts")])/sum(merged_file[,c("raw_counts_normalized_counts")])

merged_file_subset<-subset(merged_file,merged_file$sum >= 10)

pdf(paste0(path_plot,"miscelleneous/raw_histogram_after_threshold.pdf"))
print((histogram_for_raw_count(merged_file_subset,
                               merged_file_subset$raw_counts_count,
                               merged_file_subset$Codon)+
         scale_x_continuous(trans = 'log10')+
         labs(title = "histogram plot for raw count\n", 
              x = "Log of raw read Counts", 
              y = " Frequency", 
              color = "Codon class\n")+
         theme_bw()+
         theme(axis.title = element_text(face="bold"))))
dev.off()

pdf(paste0(path_plot,"miscelleneous/raw_histogram_after_threshold_normalized.pdf"))
print((histogram_for_raw_count(merged_file_subset,
                               merged_file_subset$raw_counts_normalized_counts,
                               merged_file_subset$Codon)+
         scale_x_continuous(trans = 'log10')+
         labs(title = "histogram plot for raw normalized count\n", 
              x = "Log of normalized raw read Counts", 
              y = " Frequency", 
              color = "Codon class\n")+
         theme_bw()+
         theme(axis.title = element_text(face="bold"))))
dev.off()
percentage_of_stop_codon<-sum(subset(merged_file_subset,merged_file_subset$Codon == "Stop Codon")[,c("raw_counts_normalized_counts")])/sum(merged_file_subset[,c("raw_counts_normalized_counts")])

write.csv(merged_file_subset,
          paste0(path_output,"merged_file_post_remoal.csv"))

merged_file_1<-merged_file_subset
merged_file_1$replicate<-"Replicate1"


a<-subset(merged_file_1,merged_file_1$Codon == "Non Stop Codon")
################### merging for replicate 2

# merging replicate 2


CDeg_12_counts<-read.csv(paste0(path_input,"12_X/raw_counts_rep2.txt"), 
                         header = TRUE, 
                         sep = "\t")
CDeg_12_byaa<-read.csv(paste0(path_input,"12_X/raw_PSI_byaa_rep2.txt"), 
                       header = TRUE, 
                       sep = "\t")
CDeg_12_bynuc<-read.csv(paste0(path_input,"12_X/raw_PSI_bynuc_rep2.txt"), 
                        header = TRUE, 
                        sep = "\t")

#########################################################################################################################################
####################################################integrating data for 12 amino acid long sequence#####################################################################################
#########################################################################################################################################


names(CDeg_12_counts)<-paste0("raw_counts_",
                              names(CDeg_12_counts))
names(CDeg_12_bynuc)<-paste0("bynuc_",
                             names(CDeg_12_bynuc))
names(CDeg_12_byaa)<-paste0("byaa_",
                            names(CDeg_12_byaa))

#########################################################################################################################################
CDeg_12_counts<-CDeg_12_counts[,c("raw_counts_experiment",
                                  "raw_counts_sub_experiment",
                                  "raw_counts_bin",
                                  "raw_counts_fraction",
                                  "raw_counts_f",
                                  "raw_counts_sequence",
                                  "raw_counts_sum",
                                  "raw_counts_translation",
                                  "raw_counts_count",
                                  "raw_counts_normalized_counts")]
CDeg_12_bynuc<-CDeg_12_bynuc[,c("bynuc_experiment",
                                "bynuc_sub_experiment",
                                "bynuc_sequence",
                                "bynuc_translation",
                                "bynuc_PSI")]
CDeg_12_byaa<-CDeg_12_byaa[,c("byaa_experiment",
                              "byaa_sub_experiment",
                              "byaa_translation",
                              "byaa_pooled_PSI")]

CDeg_12_counts$identifier<-paste0(CDeg_12_counts$raw_counts_experiment,"_",
                                  CDeg_12_counts$raw_counts_sub_experiment ,"_",
                                  CDeg_12_counts$raw_counts_sequence,"_",
                                  CDeg_12_counts$raw_counts_translation)

CDeg_12_bynuc$identifier<-paste0(CDeg_12_bynuc$bynuc_experiment,"_",
                                 CDeg_12_bynuc$bynuc_sub_experiment ,"_",
                                 CDeg_12_bynuc$bynuc_sequence,"_",
                                 CDeg_12_bynuc$bynuc_translation)

merged_file<-merge(CDeg_12_counts,
                   CDeg_12_bynuc,
                   by= "identifier")

merged_file$identifier<-paste0(merged_file$raw_counts_experiment,"_",
                               merged_file$raw_counts_sub_experiment,"_",
                               merged_file$bynuc_translation)

CDeg_12_byaa$identifier<-paste0(CDeg_12_byaa$byaa_experiment,"_",
                                CDeg_12_byaa$byaa_sub_experiment ,"_",
                                CDeg_12_byaa$byaa_translation)

merged_file<-merge(merged_file,
                   CDeg_12_byaa,
                   by="identifier")

merged_file$raw_counts_translation<-as.factor(unlist(merged_file$raw_counts_translation))
merged_file$bynuc_translation<-as.factor(unlist(merged_file$bynuc_translation))
merged_file$byaa_translation<-as.factor(unlist(merged_file$byaa_translation))

merged_file$Codon<-ifelse(grepl("\\*",merged_file$raw_counts_translation),"Stop Codon","Non Stop Codon")

merged_file$group_raw_PSI<-cut(merged_file$bynuc_PSI,
                               c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                               c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))

merged_file$group_raw_pooled_PSI<-cut(merged_file$byaa_pooled_PSI,
                                      c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                                      c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))

q<-merged_file %>% 
  group_by(raw_counts_translation) %>% 
  summarise(sum = sum(raw_counts_count))

merged_file<-merge(merged_file,
                   q, 
                   id="raw_counts_translation")

merged_file_subset<-subset(merged_file,merged_file$sum >= 10)

write.csv(merged_file_subset,
          paste0(path_output,"merged_file_post_remoal_replicate2.csv"))

merged_file_2<-merged_file_subset
merged_file_2$replicate<-"Replicate2"


b<-subset(merged_file_2,merged_file_2$Codon == "Non Stop Codon" & merged_file_2$raw_counts_sub_experiment == "12Xb")
b<-b[,c("raw_counts_translation","byaa_pooled_PSI")]
b<-unique(b)

a<-subset(merged_file_1,merged_file_1$Codon == "Non Stop Codon" )
a<-a[,c("raw_counts_translation","byaa_pooled_PSI")]
a<-unique(a)


names(a)<-c("Translation","PSI_rep1")
names(b)<-c("Translation","PSI_rep2")

d<-merge(a,b,id = "Translation",all = TRUE)

write.csv(d,
          paste0(path_output,"dataset_to_submit/Library_12_x.csv"))

deg_list<-subset(merged_file_1,merged_file_1$byaa_pooled_PSI < 0.5 & merged_file_1$Codon == "Non Stop Codon" & merged_file_1$sum >=10)
###########################################################
# for other dataset


### -------------------------------------- merge for replicate 2 --------------------------------------------------------------###


raw_counts<-read.csv(paste0(path_input,"other_dataset/raw_counts_rep2.txt"), 
                     sep = "\t", 
                     header = TRUE)
byaa<-read.csv(paste0(path_input,"other_dataset/raw_PSI_byaa_rep2.txt"), 
               sep = "\t", 
               header = TRUE)
bynuc<-read.csv(paste0(path_input,"other_dataset/raw_PSI_bynuc_rep2.txt"), 
                sep = "\t", 
                header = TRUE)

# merge_file

names(raw_counts)<-paste0("raw_counts_",
                          names(raw_counts))
names(bynuc)<-paste0("bynuc_",
                     names(bynuc))
names(byaa)<-paste0("byaa_",
                    names(byaa))

raw_counts<-raw_counts[,c("raw_counts_experiment",
                          "raw_counts_sub_experiment",
                          "raw_counts_bin",
                          "raw_counts_fraction",
                          "raw_counts_f",
                          "raw_counts_sequence",
                          "raw_counts_sum",
                          "raw_counts_translation",
                          "raw_counts_count",
                          "raw_counts_normalized_counts"  )]
bynuc<-bynuc[,c("bynuc_experiment",
                "bynuc_sub_experiment",
                "bynuc_sequence",
                "bynuc_translation",
                "bynuc_PSI" )]
byaa<-byaa[,c("byaa_experiment",
              "byaa_sub_experiment",
              "byaa_translation",
              "byaa_pooled_PSI" )]

raw_counts$identifier<-paste0(raw_counts$raw_counts_experiment,"_",
                              raw_counts$raw_counts_sub_experiment ,"_",
                              raw_counts$raw_counts_sequence,"_",
                              raw_counts$raw_counts_translation)

bynuc$identifier<-paste0(bynuc$bynuc_experiment,"_",
                         bynuc$bynuc_sub_experiment ,"_",
                         bynuc$bynuc_sequence,"_",
                         bynuc$bynuc_translation)

merged_file<-merge(raw_counts,
                   bynuc,
                   by= "identifier")

merged_file$identifier<-paste0(merged_file$raw_counts_experiment,"_",
                               merged_file$raw_counts_sub_experiment,"_",
                               merged_file$bynuc_translation)

byaa$identifier<-paste0(byaa$byaa_experiment,"_",
                        byaa$byaa_sub_experiment ,"_",
                        byaa$byaa_translation)

merged_file<-merge(merged_file,
                   byaa,
                   by="identifier")

merged_file$raw_counts_translation<-as.factor(unlist(merged_file$raw_counts_translation))
merged_file$bynuc_translation<-as.factor(unlist(merged_file$bynuc_translation))
merged_file$byaa_translation<-as.factor(unlist(merged_file$byaa_translation))

merged_file$Codon<-ifelse(grepl("\\*",merged_file$raw_counts_translation),"Stop Codon","Non Stop Codon")

#some manual checks
merged_file$group_raw_PSI<-cut(merged_file$bynuc_PSI,
                               c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                               c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))

merged_file$group_raw_pooled_PSI<-cut(merged_file$byaa_pooled_PSI,
                                      c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                                      c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))


merged_file$replicate<-"Replicate2"
merged_file_1<-merged_file

### -------------------------------------- merge for replicate 3 --------------------------------------------------------------###


raw_counts<-read.csv(paste0(path_input,"other_dataset/raw_counts_rep3.txt"), 
                     sep = "\t", 
                     header = TRUE)
byaa<-read.csv(paste0(path_input,"other_dataset/raw_PSI_byaa_rep3.txt"), 
               sep = "\t", 
               header = TRUE)
bynuc<-read.csv(paste0(path_input,"other_dataset/raw_PSI_bynuc_rep3.txt"), 
                sep = "\t", 
                header = TRUE)

# merge_file

names(raw_counts)<-paste0("raw_counts_",
                          names(raw_counts))
names(bynuc)<-paste0("bynuc_",
                     names(bynuc))
names(byaa)<-paste0("byaa_",
                    names(byaa))

raw_counts<-raw_counts[,c("raw_counts_experiment",
                          "raw_counts_sub_experiment",
                          "raw_counts_bin",
                          "raw_counts_fraction",
                          "raw_counts_f",
                          "raw_counts_sequence",
                          "raw_counts_sum",
                          "raw_counts_translation",
                          "raw_counts_count",
                          "raw_counts_normalized_counts"  )]
bynuc<-bynuc[,c("bynuc_experiment",
                "bynuc_sub_experiment",
                "bynuc_sequence",
                "bynuc_translation",
                "bynuc_PSI" )]
byaa<-byaa[,c("byaa_experiment",
              "byaa_sub_experiment",
              "byaa_translation",
              "byaa_pooled_PSI" )]

raw_counts$identifier<-paste0(raw_counts$raw_counts_experiment,"_",
                              raw_counts$raw_counts_sub_experiment ,"_",
                              raw_counts$raw_counts_sequence,"_",
                              raw_counts$raw_counts_translation)

bynuc$identifier<-paste0(bynuc$bynuc_experiment,"_",
                         bynuc$bynuc_sub_experiment ,"_",
                         bynuc$bynuc_sequence,"_",
                         bynuc$bynuc_translation)

merged_file<-merge(raw_counts,
                   bynuc,
                   by= "identifier")

merged_file$identifier<-paste0(merged_file$raw_counts_experiment,"_",
                               merged_file$raw_counts_sub_experiment,"_",
                               merged_file$bynuc_translation)

byaa$identifier<-paste0(byaa$byaa_experiment,"_",
                        byaa$byaa_sub_experiment ,"_",
                        byaa$byaa_translation)

merged_file<-merge(merged_file,
                   byaa,
                   by="identifier")

merged_file$raw_counts_translation<-as.factor(unlist(merged_file$raw_counts_translation))
merged_file$bynuc_translation<-as.factor(unlist(merged_file$bynuc_translation))
merged_file$byaa_translation<-as.factor(unlist(merged_file$byaa_translation))

merged_file$Codon<-ifelse(grepl("\\*",merged_file$raw_counts_translation),"Stop Codon","Non Stop Codon")

merged_file$group_raw_PSI<-cut(merged_file$bynuc_PSI,
                               c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                               c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))

merged_file$group_raw_pooled_PSI<-cut(merged_file$byaa_pooled_PSI,
                                      c(-Inf,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
                                      c("Between_0_and_0.1","Between_0.1_and_0.2","Between_0.2_and_0.3","Between_0.3_and_0.4","Between_0.4_and_0.5","Between_0.5_and_0.6","Between_0.6_and_0.7","Between_0.7_and_0.8","Between_0.8_and_0.9","Between_0.9_and_1"))
merged_file$replicate<-"Replicate3"
merged_file_2<-merged_file

### -------------------------------------- merge both replicates and write the file in drive --------------------------------------------------------------###

merged_file<-as.data.frame(rbind(merged_file_1, merged_file_2))
write.csv(merged_file,
          paste0(path_output,"/merged_file_for_rep2and3.csv"))


####


#### -------------------------------------- for correlation among different dataset --------------------------------------------------------------###
merged_file<-read.csv(paste0(path_output,"/merged_file_for_rep2and3.csv"))
deg<-read.csv( paste0(path_input,"/other_dataset/degron_effect_d1_WT.csv"))
merged_file_non_stop<-subset(merged_file,merged_file$Codon == "Non Stop Codon")

merged_file_subset<-subset(merged_file_non_stop, merged_file_non_stop$raw_counts_translation %in% deg$raw_counts_translation)
merged_file_subset$sub_experiment<-paste0(merged_file_subset$replicate,"_",merged_file_subset$raw_counts_sub_experiment)

#taking all dataset to put them in column format after check

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

PSI_merged<-merge(a,b, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,c, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,d, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,e, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,f, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,g, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,h, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,i, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,j, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,k, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,l, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,m, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged,n, by = "raw_counts_translation",all = TRUE)
PSI_merged<-merge(PSI_merged, deg, by = "raw_counts_translation",all = TRUE)

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


deg_list<-readxl::read_xlsx(paste0(path_input,"/other_dataset/4726.xlsx"), sheet = 1)
PSI_merged$type<-ifelse(row.names(PSI_merged) %in% deg_list$translation, "degron","Stable_control")

write.csv(PSI_merged,
          paste0(path_output,"/PSI_file_for_other_dataset.csv"))


# for yeast Cterminal
merged_file_subset<-merged_file_non_stop
merged_file_subset$sub_experiment<-paste0(merged_file_subset$replicate,"_",merged_file_subset$raw_counts_sub_experiment)
merged_file_subset<-merged_file_subset[,c("sub_experiment","raw_counts_translation","byaa_pooled_PSI")]

merged_file_subset_non_stop<-subset(merged_file_subset,merged_file_subset$byaa_pooled_PSI < 0.5)
unique(merged_file_subset$sub_experiment)

o<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_WT_Cter")[,c("raw_counts_translation","byaa_pooled_PSI")]

p<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate2_d1_Cter")[,c("raw_counts_translation","byaa_pooled_PSI")]
q<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_WT_Cter")[,c("raw_counts_translation","byaa_pooled_PSI")]
r<-subset(merged_file_subset,merged_file_subset$sub_experiment == "Replicate3_d1_Cter")[,c("raw_counts_translation","byaa_pooled_PSI")]

o <- o[complete.cases(o), ] 
o<-unique(o)
o<-o[!duplicated(o$raw_counts_translation),]

p <- p[complete.cases(p), ] 
p<-unique(p)
p<-p[!duplicated(p$raw_counts_translation), ]

q <- q[complete.cases(q), ] 
q<-unique(q)
q<-q[!duplicated(q$raw_counts_translation), ]

r <- r[complete.cases(r), ] 
r<-unique(r)
r<-r[!duplicated(r$raw_counts_translation), ]


CTer_WT<-merge(o,q, by = "raw_counts_translation", all = TRUE)
CTer_das1<-merge(p,r, by = "raw_counts_translation", all = TRUE)
CTer<-merge(CTer_WT,CTer_das1, by = "raw_counts_translation", all = TRUE)

names(CTer)<-c("Translation","WT_replicate2","WT_replicate3","Das1_replicate2","Das1_replicate3")

yeast<-read.csv(paste0(path_input,"yeast_prot.tsv"), sep = "\t")
yeast$raw_counts_translation<-substr(yeast$Sequence,(nchar(yeast$Sequence)-11),(nchar(yeast$Sequence)))

yeast_CTer<-subset(CTer,CTer$Translation %in% yeast$raw_counts_translation)
names(yeast)[11]<-"Translation"
yeast_CTer<-merge(yeast, yeast_CTer, by = "Translation", all.x = TRUE)
yeast_CTer$type<-ifelse(row.names(yeast_CTer) %in% deg_list$translation, "degron","Stable_control")

write.csv(yeast_CTer,
          paste0(path_output,"/yeast_Cter_all.csv"))

yeast_CTer<-yeast_CTer[,c("Translation",
                          "Entry",
                          "Entry.Name",
                          "Gene.Names",
                          "WT_replicate2",
                          "WT_replicate3",                        
                          "Das1_replicate2",
                          "Das1_replicate3",
                          "type")]
names(yeast_CTer)<-c("Translation","UniprotID","Gene_Name","Synonymous_Gene_Name",
                     "WT_replicate2",
                     "WT_replicate3",                        
                     "Das1_replicate2",
                     "Das1_replicate3",
                     "type")


write.csv(yeast_CTer,
          paste0(path_output,"/dataset_to_submit/yeast_Cter.csv"))





##################################################################################################


#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
## code for creating dataset for deep learning
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
merged_file<-read.csv(paste0(path_output,"merged_file_post_remoal.csv"))
integrated_file_non_stop_codon<-subset(merged_file,
                                       merged_file$Codon == "Non Stop Codon")
integrated_file_non_stop_codon$label<-ifelse(integrated_file_non_stop_codon$byaa_pooled_PSI > 0.5,0,1)
dataset<-integrated_file_non_stop_codon[,c("byaa_translation","label","byaa_pooled_PSI")]
dataset<-unique(dataset)

amino_acid<-dataset$byaa_translation
dataset_creation<-create_one_hot_db(amino_acid,12)
row.names(dataset_creation)<-amino_acid



## creating property dataset

hydrophobicity<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
charge<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
molecular_weight<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
mass_over_charge<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
aliphatic_index<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
isoelectric_point<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
protein_interaction_index<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))

for (i in 1:12) {
  hydrophobicity[,i]<-hydrophobicity(substr(row.names(dataset_creation),i,i),scale = "KyteDoolittle")
  charge[,i]<-charge(substr(row.names(dataset_creation),i,i), pH = 7, pKscale = "Lehninger")
  molecular_weight[,i]<-mw(substr(row.names(dataset_creation),i,i))
  mass_over_charge[,i]<-mz(substr(row.names(dataset_creation),i,i))
  aliphatic_index[,i]<-aIndex(substr(row.names(dataset_creation),i,i))
  isoelectric_point[,i]<-pI(substr(row.names(dataset_creation),i,i))
  protein_interaction_index[,i]<-boman(substr(row.names(dataset_creation),i,i))
  
}


row.names(hydrophobicity)<-row.names(dataset_creation)
row.names(charge)<-row.names(dataset_creation)
row.names(molecular_weight)<-row.names(dataset_creation)
row.names(mass_over_charge)<-row.names(dataset_creation)
row.names(aliphatic_index)<-row.names(dataset_creation)
row.names(isoelectric_point)<-row.names(dataset_creation)
row.names(protein_interaction_index)<-row.names(dataset_creation)

names(hydrophobicity)<-paste0(names(hydrophobicity),"_hydrophobicity")
names(charge)<-paste0(names(charge),"_charge")
names(molecular_weight)<-paste0(names(molecular_weight),"_molecularWeight")
names(mass_over_charge)<-paste0(names(mass_over_charge),"_massOverCharge")
names(aliphatic_index)<-paste0(names(aliphatic_index),"_aliphaticIndex")
names(isoelectric_point)<-paste0(names(isoelectric_point),"_isoelectricPoint")
names(protein_interaction_index)<-paste0(names(protein_interaction_index),"potentialProteinInteractionIndex")

polarity<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
hydrogen_bonding<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
for (i in 1:12) {
  a<-crucianiProperties(substr(row.names(dataset_creation),i,i))
  for (j in 1:nrow(dataset_creation)) {
    polarity[j,i]<-a[[j]][1]
    hydrogen_bonding[j,i]<-a[[j]][3]
  }
  print(i)
}
row.names(polarity)<-row.names(dataset_creation)
row.names(hydrogen_bonding)<-row.names(dataset_creation)

names(polarity)<-paste0(names(polarity),"_polarity")
names(hydrogen_bonding)<-paste0(names(hydrogen_bonding),"_hydrogenBonding")

alpha_and_propensity<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
bulky_properties<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
Compositional_characterstic_index<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
local_flexibility<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
electronic_properties<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
for (i in 1:12) {
  a<-fasgaiVectors(substr(row.names(dataset_creation),i,i))
  for (j in 1:nrow(dataset_creation)) {
    alpha_and_propensity[j,i]<-a[[j]][2]
    bulky_properties[j,i]<-a[[j]][3]
    Compositional_characterstic_index[j,i]<-a[[j]][4]
    local_flexibility[j,i]<-a[[j]][5]
    electronic_properties[j,i]<-a[[j]][6]
    
  }
  print(i)
}

row.names(alpha_and_propensity)<-row.names(dataset_creation)
row.names(bulky_properties)<-row.names(dataset_creation)
row.names(Compositional_characterstic_index)<-row.names(dataset_creation)
row.names(local_flexibility)<-row.names(dataset_creation)
row.names(electronic_properties)<-row.names(dataset_creation)

names(alpha_and_propensity)<-paste0(names(alpha_and_propensity),"_alphaAndTurnPropensity")
names(bulky_properties)<-paste0(names(bulky_properties),"_bulkyProperties")
names(Compositional_characterstic_index)<-paste0(names(Compositional_characterstic_index),"_compositionalCharactersticIndex")
names(local_flexibility)<-paste0(names(local_flexibility),"_localFlexibility")
names(electronic_properties)<-paste0(names(electronic_properties),"_electronicProperties")

helix_bend_preference<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
side_chain_size<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
extended_structure_preference<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
doule_bend_preference<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
partial_specific_volume<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
flat_extended_preference<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
occurance_alpha_region<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
pkc<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))
surrounding_hydrophobicity<-as.data.frame(matrix(nrow = nrow(dataset_creation),ncol = 12))

for (i in 1:12) {
  a<-kideraFactors(substr(row.names(dataset_creation),i,i))
  for (j in 1:nrow(dataset_creation)) {
    helix_bend_preference[j,i]<-a[[j]][1]
    side_chain_size[j,i]<-a[[j]][1]
    extended_structure_preference[j,i]<-a[[j]][3]
    doule_bend_preference[j,i]<-a[[j]][5]
    partial_specific_volume[j,i]<-a[[j]][6]
    flat_extended_preference[j,i]<-a[[j]][7]
    occurance_alpha_region[j,i]<-a[[j]][8]
    pkc[j,i]<-a[[j]][9]
    surrounding_hydrophobicity[j,i]<-a[[j]][10]
    
    
  }
  print(i)
}

row.names(helix_bend_preference)<-row.names(dataset_creation)
row.names(side_chain_size)<-row.names(dataset_creation)
row.names(extended_structure_preference)<-row.names(dataset_creation)
row.names(doule_bend_preference)<-row.names(dataset_creation)
row.names(partial_specific_volume)<-row.names(dataset_creation)
row.names(flat_extended_preference)<-row.names(dataset_creation)
row.names(occurance_alpha_region)<-row.names(dataset_creation)
row.names(pkc)<-row.names(dataset_creation)
row.names(surrounding_hydrophobicity)<-row.names(dataset_creation)

names(helix_bend_preference)<-paste0(names(helix_bend_preference),"_helixBendPreference")
names(side_chain_size)<-paste0(names(side_chain_size),"_sideChainSize")
names(extended_structure_preference)<-paste0(names(extended_structure_preference),"_extendedStructurePreference")
names(doule_bend_preference)<-paste0(names(doule_bend_preference),"_doubleBendPreference")
names(partial_specific_volume)<-paste0(names(partial_specific_volume),"_partialSpecificVolume")
names(flat_extended_preference)<-paste0(names(flat_extended_preference),"_flatExtendedPreference")
names(occurance_alpha_region)<-paste0(names(occurance_alpha_region),"_occuranceAlphaRegion")
names(pkc)<-paste0(names(pkc),"_pkc")
names(surrounding_hydrophobicity)<-paste0(names(surrounding_hydrophobicity),"_surroundingHydrophobicity")

all_properties<-as.data.frame(cbind(
  hydrophobicity(row.names(dataset_creation),scale = "KyteDoolittle"),
  charge(row.names(dataset_creation), pH = 7, pKscale = "Lehninger"),
  mw(row.names(dataset_creation)),
  mz(row.names(dataset_creation)),
  aIndex(row.names(dataset_creation)),
  pI(row.names(dataset_creation)),
  boman(row.names(dataset_creation))
))
a<-crucianiProperties(row.names(dataset_creation))
for (i in 1:nrow(dataset_creation)) {
  all_properties[i,8]<-a[[i]][1]
  all_properties[i,9]<-a[[i]][3]
}
a<-fasgaiVectors(row.names(dataset_creation))
for (j in 1:nrow(dataset_creation)) {
  all_properties[j,10]<-a[[j]][2]
  all_properties[j,11]<-a[[j]][3]
  all_properties[j,12]<-a[[j]][4]
  all_properties[j,13]<-a[[j]][5]
  all_properties[j,14]<-a[[j]][6]
  
}
a<-kideraFactors(row.names(dataset_creation))
for (j in 1:nrow(dataset_creation)) {
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
names(all_properties)<-c("hydrophobicity_overall",
                         "charge_overall",
                         "molecularWeight_overall",
                         "massOverCharge_overall",
                         "aliphaticIndex_overall",
                         "isoelectricPoint_overall",
                         "potentialProteinInteractionIndex_overall",
                         "polarity_overall",
                         "hydrogenBonds_overall",
                         "alphaAndTurnPropensity_overall",
                         "bulkyProperties_overall",
                         "compositionalCharactersticIndex_overall",
                         "localFlexibility_overall",
                         "electronicProperties_overall",
                         "helixBendPreference_overall",
                         "sideChainSize_overall",
                         "extendedStructurePreference_overall",
                         "doubleBendPreference_overall",
                         "partialSpecificVolume_overall",
                         "flatExtendedPreference_overall",
                         "occuranceAlphaRegion_overall",
                         "pkc_overall",
                         "surroundingHydrophobicity_overall"
                         
)

row.names(all_properties)<-row.names(dataset_creation)

dataset_all<-as.data.frame(cbind(dataset_creation[,1:240],
                                 hydrophobicity,
                                 charge,
                                 molecular_weight,
                                 mass_over_charge,
                                 aliphatic_index,
                                 isoelectric_point,
                                 protein_interaction_index,
                                 polarity,hydrogen_bonding,
                                 alpha_and_propensity,
                                 bulky_properties,
                                 Compositional_characterstic_index,
                                 local_flexibility,
                                 electronic_properties,
                                 helix_bend_preference,
                                 side_chain_size,
                                 extended_structure_preference,
                                 doule_bend_preference,
                                 partial_specific_volume,
                                 flat_extended_preference,
                                 occurance_alpha_region,
                                 pkc,
                                 surrounding_hydrophobicity,
                                 all_properties
))


dataset_all$Label<-dataset$label

dataset_all$PSI<-dataset$byaa_pooled_PSI
write.csv(dataset_all,paste0(path_output,"/complete_dataset_without_normalization.csv"))

###
#preparing dataset for deep learning

integrated_file_non_stop_codon$Class<-ifelse(integrated_file_non_stop_codon$byaa_pooled_PSI < 0.5 , 
                                             1, 
                                             0)
instable<-unique(subset(integrated_file_non_stop_codon,integrated_file_non_stop_codon$Class == 1)[,c("byaa_translation")])
stable<-unique(subset(integrated_file_non_stop_codon,integrated_file_non_stop_codon$Class == 0)[,c("byaa_translation")])

sample_training_instable<-instable[sample(1:length(instable), 3750, replace=FALSE)]
sample_training_stable<-stable[sample(1:length(stable), 3750, replace=FALSE)]

sample_testing_instable<-instable[!(instable %in% (sample_training_instable))]
sample_testing_stable<-stable[!(stable %in% (sample_training_stable))]

sample_testing_stable<-sample_testing_stable[sample(1:length(sample_testing_stable), 1000, replace=FALSE)]
all<-(c(as.character(sample_training_instable),as.character(sample_training_stable),as.character(sample_testing_instable),as.character(sample_testing_stable)))

datasetcomplete<-dataset_all
all<-subset(datasetcomplete,datasetcomplete$X %in% all)
row.names(all)<-all$X
all$X<-NULL
normalized_complete_dataset<-BBmisc::normalize(all[,1:539], method = "range", range = c(0,1))
normalized_complete_dataset$Label<-all$Label
normalized_complete_dataset$PSI<-all$PSI

instabletrainingset<-subset(normalized_complete_dataset,
                            row.names(normalized_complete_dataset) %in% sample_training_instable)
instabletrainingset$Label<-1
stabletrainingset<-subset(normalized_complete_dataset,
                          row.names(normalized_complete_dataset) %in% sample_training_stable)
stabletrainingset$Label<-0
training_dataset<-as.data.frame(rbind(instabletrainingset,stabletrainingset))
instabletestingset<-subset(normalized_complete_dataset,
                           row.names(normalized_complete_dataset) %in% sample_testing_instable)
instabletestingset$Label<-1
stabletestingset<-subset(normalized_complete_dataset,
                         row.names(normalized_complete_dataset) %in% sample_testing_stable)
stabletestingset$Label<-0
testing_dataset<-as.data.frame(rbind(instabletestingset,
                                     stabletestingset))

training_dataset<-training_dataset[sample(1:nrow(training_dataset), 
                                          nrow(training_dataset), 
                                          replace=FALSE),]
testing_dataset<-testing_dataset[sample(1:nrow(testing_dataset), 
                                        nrow(testing_dataset), 
                                        replace=FALSE),]

write.csv(training_dataset,
          paste0(path_output,
                 "training_data_after_normalization.csv"))
write.csv(testing_dataset,
          paste0(path_output,"
                 testing_data_after_normalization.csv"))

#training_dataset<-training_data
#testing_dataset<-testing_data

training_dataset<-training_dataset[,c(1:252,277:300,313:324,349:372,517,520,521,523,526,527,540,541)]
testing_dataset<-testing_dataset[,c(1:252,277:300,313:324,349:372,517,520,521,523,526,527,540,541)]

write.csv(training_dataset,
          paste0(path_output,
                 "training_data_after_normalization_for_model.csv"))
write.csv(testing_dataset,
          paste0(path_output,
                 "testing_data_after_normalization_for_model.csv"))

#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

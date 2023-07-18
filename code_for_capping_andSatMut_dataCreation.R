# code for dataset creation
# saturation mutagenesis 

library(reshape2)

IN1_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/IN1.xlsx", sheet =1)
IN2_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/IN2.xlsx", sheet =1)
LN1_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/LN1.xlsx", sheet =1)
LN3_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/LN3.xlsx", sheet =1)
VN1_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/VN1.xlsx", sheet =1)
VN2_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/VN2.xlsx", sheet =1)
VN3_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/VN3.xlsx", sheet =1)


IN1_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/IN1.xlsx", sheet =2)
IN2_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/IN2.xlsx", sheet =2)
LN1_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/LN1.xlsx", sheet =2)
LN3_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/LN3.xlsx", sheet =2)
VN1_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/VN1.xlsx", sheet =2)
VN2_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/VN2.xlsx", sheet =2)
VN3_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/VN3.xlsx", sheet =2)

names(IN1_WT)<-c("Saturation","L1","V1","V2","Y","L2","S","C1","V3","A","C2","I","N")
names(IN1_DT)<-c("Saturation","L1","V1","V2","Y","L2","S","C1","V3","A","C2","I","N")
names(IN2_WT)<-c("Saturation","K1","K2","H","K3","K4","G","Q1","K5","Q2","K6","I","N")
names(IN2_DT)<-c("Saturation","K1","K2","H","K3","K4","G","Q1","K5","Q2","K6","I","N")
names(LN1_WT)<-c("Saturation","P","D","K1","R","N1","H","L1","I","K2","L2","L3","N2")
names(LN1_DT)<-c("Saturation","P","D","K1","R","N1","H","L1","I","K2","L2","L3","N2")
names(LN3_WT)<-c("Saturation","R","L1",	"I",	"Q1","P1",	"L2",	"P2",	"V",	"Q2",	"F",	"L3",	"N")
names(LN3_DT)<-c("Saturation","R","L1",	"I",	"Q1","P1",	"L2",	"P2",	"V",	"Q2",	"F",	"L3",	"N")
names(VN1_WT)<-c("Saturation","C",	"L1",	"F1",	"H",	"Q",	"F2",	"I",	"T",	"L2",	"L3",	"V",	"N")
names(VN1_DT)<-c("Saturation","C",	"L1",	"F1",	"H",	"Q",	"F2",	"I",	"T",	"L2",	"L3",	"V",	"N")
names(VN3_WT)<-c("Saturation","I1"	,"I2"	,"G",	"C1",	"I3",	"A",	"I4",	"C2",	"L",	"M",	"V",	"N")
names(VN3_DT)<-c("Saturation","I1"	,"I2"	,"G",	"C1",	"I3",	"A",	"I4",	"C2",	"L",	"M",	"V",	"N")
names(VN2_WT)<-c("Saturation","F1",  "S" ,     "M"   ,   "R1" , "W1",  "T" ,     "F2" , "R2" , "Y"    ,  "W2" ,"V" ,     "N")
names(VN2_DT)<-c("Saturation","F1",  "S" ,     "M"   ,   "R1" , "W1",  "T" ,     "F2" , "R2" , "Y"    ,  "W2" ,"V" ,     "N")

IN1<-"LVVYLSCVACIN"

IN1_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(IN1_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(IN1,1,i-1),
    IN1_DT$Saturation,
    substr(IN1,i+1,nchar(IN1))
  )
  
  
  b<-IN1_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  IN1_DT_dataset<-as.data.frame(rbind(IN1_DT_dataset,c))
}


IN1_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(IN1_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(IN1,1,i-1),
    IN1_WT$Saturation,
    substr(IN1,i+1,nchar(IN1))
  )
  
  
  b<-IN1_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  IN1_WT_dataset<-as.data.frame(rbind(IN1_WT_dataset,c))
}


IN1_WT_dataset<-IN1_WT_dataset[2:nrow(IN1_WT_dataset),]
IN1_DT_dataset<-IN1_DT_dataset[2:nrow(IN1_DT_dataset),]
IN1_dataset<-as.data.frame(cbind(IN1_WT_dataset, IN1_DT_dataset))
IN1_dataset$Peptides<-NULL
IN1_dataset<-IN1_dataset[,c(2,1,3)]

write.csv(IN1_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/IN1_dataset.csv")



IN2<-"KKHKKGQKQKIN"

IN2_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(IN2_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(IN2,1,i-1),
    IN2_DT$Saturation,
    substr(IN2,i+1,nchar(IN2))
  )
  
  
  b<-IN2_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  IN2_DT_dataset<-as.data.frame(rbind(IN2_DT_dataset,c))
}


IN2_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(IN2_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(IN2,1,i-1),
    IN2_WT$Saturation,
    substr(IN2,i+1,nchar(IN2))
  )
  
  
  b<-IN2_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  IN2_WT_dataset<-as.data.frame(rbind(IN2_WT_dataset,c))
}


IN2_WT_dataset<-IN2_WT_dataset[2:nrow(IN2_WT_dataset),]
IN2_DT_dataset<-IN2_DT_dataset[2:nrow(IN2_DT_dataset),]
IN2_dataset<-as.data.frame(cbind(IN2_WT_dataset, IN2_DT_dataset))
IN2_dataset$Peptides<-NULL
IN2_dataset<-IN2_dataset[,c(2,1,3)]

write.csv(IN2_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/IN2_dataset.csv")




LN3<-"RLIQPLPVQFLN"

LN3_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(LN3_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(LN3,1,i-1),
    LN3_DT$Saturation,
    substr(LN3,i+1,nchar(LN3))
  )
  
  
  b<-LN3_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  LN3_DT_dataset<-as.data.frame(rbind(LN3_DT_dataset,c))
}


LN3_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(LN3_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(LN3,1,i-1),
    LN3_WT$Saturation,
    substr(LN3,i+1,nchar(LN3))
  )
  
  
  b<-LN3_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  LN3_WT_dataset<-as.data.frame(rbind(LN3_WT_dataset,c))
}


LN3_WT_dataset<-LN3_WT_dataset[2:nrow(LN3_WT_dataset),]
LN3_DT_dataset<-LN3_DT_dataset[2:nrow(LN3_DT_dataset),]
LN3_dataset<-as.data.frame(cbind(LN3_WT_dataset, LN3_DT_dataset))
LN3_dataset$Peptides<-NULL
LN3_dataset<-LN3_dataset[,c(2,1,3)]

write.csv(LN3_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/LN3_dataset.csv")



VN1<-"CLFHQFITLLVN"

VN1_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(VN1_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(VN1,1,i-1),
    VN1_DT$Saturation,
    substr(VN1,i+1,nchar(VN1))
  )
  
  
  b<-VN1_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  VN1_DT_dataset<-as.data.frame(rbind(VN1_DT_dataset,c))
}


VN1_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(VN1_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(VN1,1,i-1),
    VN1_WT$Saturation,
    substr(VN1,i+1,nchar(VN1))
  )
  
  
  b<-VN1_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  VN1_WT_dataset<-as.data.frame(rbind(VN1_WT_dataset,c))
}


VN1_WT_dataset<-VN1_WT_dataset[2:nrow(VN1_WT_dataset),]
VN1_DT_dataset<-VN1_DT_dataset[2:nrow(VN1_DT_dataset),]
VN1_dataset<-as.data.frame(cbind(VN1_WT_dataset, VN1_DT_dataset))
VN1_dataset$Peptides<-NULL
VN1_dataset<-VN1_dataset[,c(2,1,3)]

write.csv(VN1_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/VN1_dataset.csv")



VN3<-"IIGCIAICLMVN"

VN3_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(VN3_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(VN3,1,i-1),
    VN3_DT$Saturation,
    substr(VN3,i+1,nchar(VN3))
  )
  
  
  b<-VN3_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  VN3_DT_dataset<-as.data.frame(rbind(VN3_DT_dataset,c))
}


VN3_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(VN3_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(VN3,1,i-1),
    VN3_WT$Saturation,
    substr(VN3,i+1,nchar(VN3))
  )
  
  
  b<-VN3_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  VN3_WT_dataset<-as.data.frame(rbind(VN3_WT_dataset,c))
}


VN3_WT_dataset<-VN3_WT_dataset[2:nrow(VN3_WT_dataset),]
VN3_DT_dataset<-VN3_DT_dataset[2:nrow(VN3_DT_dataset),]
VN3_dataset<-as.data.frame(cbind(VN3_WT_dataset, VN3_DT_dataset))
VN3_dataset$Peptides<-NULL
VN3_dataset<-VN3_dataset[,c(2,1,3)]

write.csv(VN3_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/VN3_dataset.csv")


R1_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R1.xlsx", sheet =1)
R2_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R2.xlsx", sheet =1)
R3_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R3.xlsx", sheet =1)
R4_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R4.xlsx", sheet =1)
R5_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R5.xlsx", sheet =1)
Atg1C_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/Atg1C.xlsx", sheet =1)
Rpa12C_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/Rpa12C.xlsx", sheet =1)
W1_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/W1.xlsx", sheet =1)
W2_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/W2.xlsx", sheet =1)
W3_WT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/W3.xlsx", sheet =1)


R1_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R1.xlsx", sheet =2)
R2_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R2.xlsx", sheet =2)
R3_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R3.xlsx", sheet =2)
R4_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R4.xlsx", sheet =2)
R5_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/R5.xlsx", sheet =2)
Atg1C_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/Atg1C.xlsx", sheet =2)
Rpa12C_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/Rpa12C.xlsx", sheet =2)
W1_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/W1.xlsx", sheet =2)
W2_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/W2.xlsx", sheet =2)
W3_DT<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/Saturation_mutagenesis/W3.xlsx", sheet =2)


names(R1_WT)<-c("Saturation","A","I","T","L1","P1","M","C","P2","L2","L3","L4","R")
names(R1_DT)<-c("Saturation","A","I","T","L1","P1","M","C","P2","L2","L3","L4","R")
names(R2_WT)<-c("Saturation","R1","S","L1","R2","V","L2","L3","P","L4","Y","F","R3")
names(R2_DT)<-c("Saturation","R1","S","L1","R2","V","L2","L3","P","L4","Y","F","R3")
names(R3_WT)<-c("Saturation","C","L1","V","S","L2","P","I","L3","L4","F1","F2","R")
names(R3_DT)<-c("Saturation","C","L1","V","S","L2","P","I","L3","L4","F1","F2","R")
names(R4_WT)<-c("Saturation","R1","V",	"T",	"R2","L",	"D",	"A",	"I1",	"Y",	"F",	"R3",	"I2")
names(R4_DT)<-c("Saturation","R1","V",	"T",	"R2","L",	"D",	"A",	"I1",	"Y",	"F",	"R3",	"I2")
names(R5_WT)<-c("Saturation","G1",	"L",	"Y",	"M",	"C1",	"S",	"V",	"G2",	"C2",	"W",	"F",	"R")
names(R5_DT)<-c("Saturation","G1",	"L",	"Y",	"M",	"C1",	"S",	"V",	"G2",	"C2",	"W",	"F",	"R")
names(Atg1C_WT)<-c("Saturation","L1",	"K1",	"I",	"L2",	"R",	"Q1",	"K2",	"M",	"N1",	"H",	"Q2",	"N2")
names(Atg1C_DT)<-c("Saturation","L1",	"K1",	"I",	"L2",	"R",	"Q1",	"K2",	"M",	"N1",	"H",	"Q2",	"N2")
names(Rpa12C_WT)<-c("Saturation","C1",	"T1",	"S",	"C2",	"G",	"Y",	"K",	"F",	"R",	"T2",	"N1",	"N2")
names(Rpa12C_DT)<-c("Saturation","C1",	"T1",	"S",	"C2",	"G",	"Y",	"K",	"F",	"R",	"T2",	"N1",	"N2")
names(W1_WT)<-c("Saturation","V1"	,"R1"	,"G1",	"A",	"V2",	"G2",	"G3",	"W",	"R2",	"L",	"V3",	"G4")
names(W1_DT)<-c("Saturation","V1"	,"R1"	,"G1",	"A",	"V2",	"G2",	"G3",	"W",	"R2",	"L",	"V3",	"G4")
names(W2_WT)<-c("Saturation","A1"	,"R1"	,"W1",	"R2",	"V",	"G1",	"L",	"W2",	"R3",	"G2",	"A2",	"S")
names(W2_DT)<-c("Saturation","A1"	,"R1"	,"W1",	"R2",	"V",	"G1",	"L",	"W2",	"R3",	"G2",	"A2",	"S")
names(W3_WT)<-c("Saturation","P",  "S1" ,     "L"   ,   "A1" , "A2",  "S2" ,     "F" , "W" , "R"    ,  "V1" ,"V2" ,     "G")
names(W3_DT)<-c("Saturation","P",  "S1" ,     "L"   ,   "A1" , "A2",  "S2" ,     "F" , "W" , "R"    ,  "V1" ,"V2" ,     "G")



R2<-"RSLRVLLPLYFR"

R2_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(R2_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(R2,1,i-1),
    R2_DT$Saturation,
    substr(R2,i+1,nchar(R2))
  )
  
  
  b<-R2_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  R2_DT_dataset<-as.data.frame(rbind(R2_DT_dataset,c))
}


R2_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(R2_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(R2,1,i-1),
    R2_WT$Saturation,
    substr(R2,i+1,nchar(R2))
  )
  
  
  b<-R2_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  R2_WT_dataset<-as.data.frame(rbind(R2_WT_dataset,c))
}


R2_WT_dataset<-R2_WT_dataset[2:nrow(R2_WT_dataset),]
R2_DT_dataset<-R2_DT_dataset[2:nrow(R2_DT_dataset),]
R2_dataset<-as.data.frame(cbind(R2_WT_dataset, R2_DT_dataset))
R2_dataset$Peptides<-NULL
R2_dataset<-R2_dataset[,c(2,1,3)]

write.csv(R2_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R2_dataset.csv")



R3<-"CLVSLPILLFFR"

R3_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(R3_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(R3,1,i-1),
    R3_DT$Saturation,
    substr(R3,i+1,nchar(R3))
  )
  
  
  b<-R3_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  R3_DT_dataset<-as.data.frame(rbind(R3_DT_dataset,c))
}


R3_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(R3_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(R3,1,i-1),
    R3_WT$Saturation,
    substr(R3,i+1,nchar(R3))
  )
  
  
  b<-R3_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  R3_WT_dataset<-as.data.frame(rbind(R3_WT_dataset,c))
}


R3_WT_dataset<-R3_WT_dataset[2:nrow(R3_WT_dataset),]
R3_DT_dataset<-R3_DT_dataset[2:nrow(R3_DT_dataset),]
R3_dataset<-as.data.frame(cbind(R3_WT_dataset, R3_DT_dataset))
R3_dataset$Peptides<-NULL
R3_dataset<-R3_dataset[,c(2,1,3)]

write.csv(R3_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R3_dataset.csv")


R4<-"RVTRLDAIYFRI"

R4_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(R4_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(R4,1,i-1),
    R4_DT$Saturation,
    substr(R4,i+1,nchar(R4))
  )
  
  
  b<-R4_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  R4_DT_dataset<-as.data.frame(rbind(R4_DT_dataset,c))
}


R4_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(R4_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(R4,1,i-1),
    R4_WT$Saturation,
    substr(R4,i+1,nchar(R4))
  )
  
  
  b<-R4_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  R4_WT_dataset<-as.data.frame(rbind(R4_WT_dataset,c))
}


R4_WT_dataset<-R4_WT_dataset[2:nrow(R4_WT_dataset),]
R4_DT_dataset<-R4_DT_dataset[2:nrow(R4_DT_dataset),]
R4_dataset<-as.data.frame(cbind(R4_WT_dataset, R4_DT_dataset))
R4_dataset$Peptides<-NULL
R4_dataset<-R4_dataset[,c(2,1,3)]

write.csv(R4_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R5_dataset.csv")


R5<-"GLYMCSVGCWFR"

R5_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(R5_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(R5,1,i-1),
    R5_DT$Saturation,
    substr(R5,i+1,nchar(R5))
  )
  
  
  b<-R5_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  R5_DT_dataset<-as.data.frame(rbind(R5_DT_dataset,c))
}


R5_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(R5_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(R5,1,i-1),
    R5_WT$Saturation,
    substr(R5,i+1,nchar(R5))
  )
  
  
  b<-R5_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  R5_WT_dataset<-as.data.frame(rbind(R5_WT_dataset,c))
}


R5_WT_dataset<-R5_WT_dataset[2:nrow(R5_WT_dataset),]
R5_DT_dataset<-R5_DT_dataset[2:nrow(R5_DT_dataset),]
R5_dataset<-as.data.frame(cbind(R5_WT_dataset, R5_DT_dataset))
R5_dataset$Peptides<-NULL
R5_dataset<-R5_dataset[,c(2,1,3)]

write.csv(R5_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R4_dataset.csv")


Atg1C<-"LKILRQKMNHQN"

Atg1C_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(Atg1C_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(Atg1C,1,i-1),
    Atg1C_DT$Saturation,
    substr(Atg1C,i+1,nchar(Atg1C))
  )
  
  
  b<-Atg1C_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  Atg1C_DT_dataset<-as.data.frame(rbind(Atg1C_DT_dataset,c))
}


Atg1C_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(Atg1C_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(Atg1C,1,i-1),
    Atg1C_WT$Saturation,
    substr(Atg1C,i+1,nchar(Atg1C))
  )
  
  
  b<-Atg1C_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  Atg1C_WT_dataset<-as.data.frame(rbind(Atg1C_WT_dataset,c))
}


Atg1C_WT_dataset<-Atg1C_WT_dataset[2:nrow(Atg1C_WT_dataset),]
Atg1C_DT_dataset<-Atg1C_DT_dataset[2:nrow(Atg1C_DT_dataset),]
Atg1C_dataset<-as.data.frame(cbind(Atg1C_WT_dataset, Atg1C_DT_dataset))
Atg1C_dataset$Peptides<-NULL
Atg1C_dataset<-Atg1C_dataset[,c(2,1,3)]

write.csv(Atg1C_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/Atg1C_dataset.csv")



Rpa12C<-"CTSCGYKFRTNN"

Rpa12C_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(Rpa12C_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(Rpa12C,1,i-1),
    Rpa12C_DT$Saturation,
    substr(Rpa12C,i+1,nchar(Rpa12C))
  )
  
  
  b<-Rpa12C_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  Rpa12C_DT_dataset<-as.data.frame(rbind(Rpa12C_DT_dataset,c))
}


Rpa12C_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(Rpa12C_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(Rpa12C,1,i-1),
    Rpa12C_WT$Saturation,
    substr(Rpa12C,i+1,nchar(Rpa12C))
  )
  
  
  b<-Rpa12C_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  Rpa12C_WT_dataset<-as.data.frame(rbind(Rpa12C_WT_dataset,c))
}


Rpa12C_WT_dataset<-Rpa12C_WT_dataset[2:nrow(Rpa12C_WT_dataset),]
Rpa12C_DT_dataset<-Rpa12C_DT_dataset[2:nrow(Rpa12C_DT_dataset),]
Rpa12C_dataset<-as.data.frame(cbind(Rpa12C_WT_dataset, Rpa12C_DT_dataset))
Rpa12C_dataset$Peptides<-NULL
Rpa12C_dataset<-Rpa12C_dataset[,c(2,1,3)]

write.csv(Rpa12C_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/Rpa12C_dataset.csv")

W1<-"VRGAVGGWRLVG"

W1_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(W1_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(W1,1,i-1),
    W1_DT$Saturation,
    substr(W1,i+1,nchar(W1))
  )
  
  
  b<-W1_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  W1_DT_dataset<-as.data.frame(rbind(W1_DT_dataset,c))
}


W1_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(W1_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(W1,1,i-1),
    W1_WT$Saturation,
    substr(W1,i+1,nchar(W1))
  )
  
  
  b<-W1_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  W1_WT_dataset<-as.data.frame(rbind(W1_WT_dataset,c))
}


W1_WT_dataset<-W1_WT_dataset[2:nrow(W1_WT_dataset),]
W1_DT_dataset<-W1_DT_dataset[2:nrow(W1_DT_dataset),]
W1_dataset<-as.data.frame(cbind(W1_WT_dataset, W1_DT_dataset))
W1_dataset$Peptides<-NULL
W1_dataset<-W1_dataset[,c(2,1,3)]

write.csv(W1_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/W1_dataset.csv")



W3<-"PSLAASFWRVVG"

W3_DT_dataset<-as.data.frame(matrix(ncol = 2))
names(W3_DT_dataset)<-c("Peptides","PSI_Das1")
for (i in 1:12) {
  a<-paste0(
    substr(W3,1,i-1),
    W3_DT$Saturation,
    substr(W3,i+1,nchar(W3))
  )
  
  
  b<-W3_DT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_Das1")
  
  W3_DT_dataset<-as.data.frame(rbind(W3_DT_dataset,c))
}


W3_WT_dataset<-as.data.frame(matrix(ncol = 2))
names(W3_WT_dataset)<-c("Peptides","PSI_WT")
for (i in 1:12) {
  a<-paste0(
    substr(W3,1,i-1),
    W3_WT$Saturation,
    substr(W3,i+1,nchar(W3))
  )
  
  
  b<-W3_WT[,i+1]
  c<-as.data.frame(cbind(a,b))
  names(c)<-c("Peptides","PSI_WT")
  
  W3_WT_dataset<-as.data.frame(rbind(W3_WT_dataset,c))
}


W3_WT_dataset<-W3_WT_dataset[2:nrow(W3_WT_dataset),]
W3_DT_dataset<-W3_DT_dataset[2:nrow(W3_DT_dataset),]
W3_dataset<-as.data.frame(cbind(W3_WT_dataset, W3_DT_dataset))
W3_dataset$Peptides<-NULL
W3_dataset<-W3_dataset[,c(2,1,3)]

write.csv(W3_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/W3_dataset.csv")


write.csv(W3_dataset,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/W3_dataset.csv")


#R4 and R5 are swapped

R1<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =1)
names(R1)<-c("Amino_Acid","WT","Doa10")
R1$Peptides<-paste0("AITLPMCPLLLR", R1$Amino_Acid)
write.csv(R1,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R1_capping.csv")
R2<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =2)
names(R2)<-c("Amino_Acid","WT","Doa10")
R2$Peptides<-paste0("RSLRVLLPLYFR", R2$Amino_Acid)
write.csv(R2,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R2_capping.csv")
R3<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =3)
names(R3)<-c("Amino_Acid","WT","Doa10")
R3$Peptides<-paste0("CLVSLPILLFFR", R3$Amino_Acid)
write.csv(R3,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R3_capping.csv")
R4<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =4)
names(R4)<-c("Amino_Acid","WT","Doa10")
R4$Peptides<-paste0("RVTRLDAIYFRI", R4$Amino_Acid)
write.csv(R4,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R5_capping.csv")
R5<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =5)
names(R5)<-c("Amino_Acid","WT","Doa10")
R5$Peptides<-paste0("GLYMCSVGCWFR", R5$Amino_Acid)
write.csv(R5,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/R4_capping.csv")
W1<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =8)
names(W1)<-c("Amino_Acid","WT","Das1")
W1$Peptides<-paste0("VRGAVGGWRLVG", W1$Amino_Acid)
write.csv(W1,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/W1_capping.csv")
W2<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =9)
names(W2)<-c("Amino_Acid","WT","Das1")
W2$Peptides<-paste0("ARWRVGLWRGAS", W2$Amino_Acid)
write.csv(W2,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/W2_capping.csv")
W3<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =10)
names(W3)<-c("Amino_Acid","WT","Das1")
W3$Peptides<-paste0("PSLAASFWRVVG", W3$Amino_Acid)
write.csv(W3,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/W3_capping.csv")

Atg1C<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =6)
names(Atg1C)<-c("Amino_Acid","WT","Das1")
Atg1C$Peptides<-paste0("LKILRQKMNHQN", Atg1C$Amino_Acid)
write.csv(Atg1C,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/Atg1C_capping.csv")

Rpa12C<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/new_data2/capping.xlsx", sheet =7)
names(Rpa12C)<-c("Amino_Acid","WT","Das1")
Rpa12C$Peptides<-paste0("CTSCGYKFRTNN", Rpa12C$Amino_Acid)
write.csv(Rpa12C,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/Rpa12C_capping.csv")


IN1<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/capping/Cdeg_cap_all.xlsx", sheet =1)
names(IN1)<-c("Amino_Acid","WT","Das1")
IN1$Peptides<-paste0("LVVYLSCVACIN", IN1$Amino_Acid)
write.csv(IN1,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/IN1_capping.csv")
IN2<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/capping/Cdeg_cap_all.xlsx", sheet =2)
names(IN2)<-c("Amino_Acid","WT","Das1")
IN2$Peptides<-paste0("KKHKKGQKQKIN", IN2$Amino_Acid)
write.csv(IN2,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/IN2_capping.csv")
LN1<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/capping/Cdeg_cap_all.xlsx", sheet =3)
names(LN1)<-c("Amino_Acid","WT","Das1")
LN1$Peptides<-paste0("PDKRNHLIKLLN", LN1$Amino_Acid)
write.csv(LN1,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/LN1_capping.csv")
LN3<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/capping/Cdeg_cap_all.xlsx", sheet =4)
names(LN3)<-c("Amino_Acid","WT","Das1")
LN3$Peptides<-paste0("RLIQPLPVQFLN", LN3$Amino_Acid)
write.csv(LN3,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/LN3_capping.csv")
VN1<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/capping/Cdeg_cap_all.xlsx", sheet =5)
names(VN1)<-c("Amino_Acid","WT","Das1")
VN1$Peptides<-paste0("CLFHQFITLLVN", VN1$Amino_Acid)
write.csv(VN1,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/VN1_capping.csv")
VN2<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/capping/Cdeg_cap_all.xlsx", sheet =6)
names(VN2)<-c("Amino_Acid","WT","Das1")
VN2$Peptides<-paste0("FSMRWTFRYWVN", VN2$Amino_Acid)
write.csv(VN2,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/VN2_capping.csv")
VN3<-readxl::read_xlsx("Y:/lab data/susmitha/edwin/for_paper/input_data/capping/Cdeg_cap_all.xlsx", sheet =7)
names(VN3)<-c("Amino_Acid","WT","Das1")
VN3$Peptides<-paste0("IIGCIAICLMVN", VN3$Amino_Acid)
write.csv(VN3,"Y:/lab data/susmitha/edwin/for_paper/output/dataset_for_submission/VN3_capping.csv")

### 
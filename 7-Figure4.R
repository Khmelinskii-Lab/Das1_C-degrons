##### for project : Analysis of degron at C terminal
##### sub text:capping and saturation mutagenesis for IN,vn,mn,ln
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")
path_model<-paste0(getwd(),"/model")

IN1<-readxl::read_xlsx(paste0(path_input,"/capping/Cdeg_cap_all.xlsx"), sheet =1)
names(IN1)<-c("Amino_Acid","WT","Das1")
IN2<-readxl::read_xlsx(paste0(path_input,"/capping/Cdeg_cap_all.xlsx"), sheet =2)
names(IN2)<-c("Amino_Acid","WT","Das1")
LN1<-readxl::read_xlsx(paste0(path_input,"/capping/Cdeg_cap_all.xlsx"), sheet =3)
names(LN1)<-c("Amino_Acid","WT","Das1")
LN3<-readxl::read_xlsx(paste0(path_input,"/capping/Cdeg_cap_all.xlsx"), sheet =4)
names(LN3)<-c("Amino_Acid","WT","Das1")
VN1<-readxl::read_xlsx(paste0(path_input,"/capping/Cdeg_cap_all.xlsx"), sheet =5)
names(VN1)<-c("Amino_Acid","WT","Das1")
VN2<-readxl::read_xlsx(paste0(path_input,"/capping/Cdeg_cap_all.xlsx"), sheet =6)
names(VN2)<-c("Amino_Acid","WT","Das1")
VN3<-readxl::read_xlsx(paste0(path_input,"/capping/Cdeg_cap_all.xlsx"), sheet =7)
names(VN3)<-c("Amino_Acid","WT","Das1")


IN1_melt<-melt(IN1, id = "Amino_Acid")
IN1_melt$value_log2<-log2(IN1_melt$value)
IN1_melt$group<-"IN1"

IN2_melt<-melt(IN2, id = "Amino_Acid")
IN2_melt$value_log2<-log2(IN2_melt$value)
IN2_melt$group<-"IN2"

LN1_melt<-melt(LN1, id = "Amino_Acid")
LN1_melt$value_log2<-log2(LN1_melt$value)
LN1_melt$group<-"LN1"

LN3_melt<-melt(LN3, id = "Amino_Acid")
LN3_melt$value_log2<-log2(LN3_melt$value)
LN3_melt$group<-"LN3"

VN1_melt<-melt(VN1, id = "Amino_Acid")
VN1_melt$value_log2<-log2(VN1_melt$value)
VN1_melt$group<-"VN1"

VN2_melt<-melt(VN2, id = "Amino_Acid")
VN2_melt$value_log2<-log2(VN2_melt$value)
VN2_melt$group<-"VN2"

VN3_melt<-melt(VN3, id = "Amino_Acid")
VN3_melt$value_log2<-log2(VN3_melt$value)
VN3_melt$group<-"VN3"

capping<-as.data.frame(rbind(IN1_melt,IN2_melt,LN1_melt,LN3_melt,VN1_melt,VN2_melt,VN3_melt))


IN1$delta<-log2(IN1$Das1/IN1$WT)
IN1$group<-"IN1"
IN2$delta<-log2(IN2$Das1/IN2$WT)
IN2$group<-"IN2"
LN1$delta<-log2(LN1$Das1/LN1$WT)
LN1$group<-"LN1"
LN3$delta<-log2(LN3$Das1/LN3$WT)
LN3$group<-"LN3"
VN1$delta<-log2(VN1$Das1/VN1$WT)
VN1$group<-"VN1"
VN2$delta<-log2(VN2$Das1/VN2$WT)
VN2$group<-"VN2"
VN3$delta<-log2(VN3$Das1/VN3$WT)
VN3$group<-"VN3"

delta<-(rbind(IN1,IN2,LN1,LN3,VN1,VN2,VN3))


capping$Amino_Acid<-factor(capping$Amino_Acid,
                          levels = rev(c("Stop","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))

delta$Amino_Acid<-factor(delta$Amino_Acid,
                            levels = rev(c("Stop","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
capping$value<-ifelse(capping$value > 0.75,0.75,capping$value)

dir.create(
  paste0(path_plot,"Figure4/"),
  showWarnings = TRUE,
  recursive = FALSE,
  mode = "0777"
)

pdf(paste0(path_plot,"Figure4/IN1_capping.pdf"))
print(
plot_grid(
ggplot(subset(capping,capping$group== "IN1"))+
  geom_tile(aes(
    x = variable  ,
    y=Amino_Acid,
    fill = value 
  ))+coord_equal()+
  theme_bw()+
  theme(text   = element_text(size = 6)) +

  scale_fill_gradientn("PSI",colours=c(
    "#00A650", 
    "#42B977", 
    "#6DC590", 
    "#98D2AA", 
    "#C4DEC4", 
    "#F0EBDE",
    "#EAC5CF", 
    "#E4A0BF", 
    "#DD7AB0", 
    "#D755A0", 
    "#CE1C89"
  ), na.value = "grey98",
  values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
  limits = c(0,0.75))+
  
  ylab("Mutations")+
  xlab("Saturated Sequence ")+
  labs(title =  "IN1 PSI "
  ),

ggplot(delta[delta$group %in% c('IN1'),])+
  geom_tile(aes(
    x = group  ,
    y=Amino_Acid,
    fill = delta 
  ))+theme_bw()+
  coord_equal()+
  theme(text   = element_text(size = 6)) +
  scale_fill_gradientn("Log2(Delta)",colours=c(
    "#F7931E", 
    "#F5AB53", 
    "#F4BB75", 
    "#F3CB98", 
    "#F1DBBB", 
    "#F0EBDE", 
    "#C4C0E4", 
    "#9996EA", 
    "#6D6BF0", 
    "#4240F6", 
    "#0000FF"
  ), na.value = "grey98",
  values = scales::rescale(c(-4,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,4)),
  limits = c(-4,4))+
  
  ylab("Mutations")+
  xlab("Saturated Sequence ")+
  labs(title =  "IN1 delta "
  )
)
)
dev.off()


######## In2

pdf(paste0(path_plot,"Figure4/IN2_capping.pdf"))
print(
  plot_grid(
    ggplot(subset(capping,capping$group== "IN2"))+
      geom_tile(aes(
        x = variable  ,
        y=Amino_Acid,
        fill = value 
      ))+coord_equal()+
      theme_bw()+
      theme(text   = element_text(size = 6)) +
      
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "IN2 fixed scale "
      ),
    
    ggplot(delta[delta$group %in% c('IN2'),])+
      geom_tile(aes(
        x = group  ,
        y=Amino_Acid,
        fill = delta 
      ))+theme_bw()+
      coord_equal()+
      theme(text   = element_text(size = 6)) +
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F5AB53", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#4240F6", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(c(-4,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,4)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "IN2 delta "
      )
  )
)
dev.off()

######## LN1


pdf(paste0(path_plot,"Figure4/LN1_capping.pdf"))
print(
  plot_grid(
    ggplot(subset(capping,capping$group== "LN1"))+
      geom_tile(aes(
        x = variable  ,
        y=Amino_Acid,
        fill = value 
      ))+coord_equal()+
      theme_bw()+
      theme(text   = element_text(size = 6)) +
      
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "LN1 fixed scale "
      ),
    
    ggplot(delta[delta$group %in% c('LN1'),])+
      geom_tile(aes(
        x = group  ,
        y=Amino_Acid,
        fill = delta 
      ))+coord_equal()+
      theme_bw()+theme(text   = element_text(size = 6)) +
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F5AB53", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#4240F6", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(c(-4,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,4)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "LN1  delta"
      )
  )
)
dev.off()

######## LN3


pdf(paste0(path_plot,"Figure4/LN3_capping.pdf"))
print(
  plot_grid(
    ggplot(subset(capping,capping$group== "LN3"))+
      geom_tile(aes(
        x = variable  ,
        y=Amino_Acid,
        fill = value 
      ))+
      theme_bw()+coord_equal()+
      theme(text   = element_text(size = 6)) +
      
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "LN3 fixed scale "
      ),
    
    ggplot(delta[delta$group %in% c('LN3'),])+
      geom_tile(aes(
        x = group  ,
        y=Amino_Acid,
        fill = delta 
      ))+coord_equal()+
      theme_bw()+theme(text   = element_text(size = 6)) +
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F5AB53", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#4240F6", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(c(-4,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,4)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "LN3  delta"
      )
  )
)
dev.off()

######## VN1


pdf(paste0(path_plot,"Figure4/VN1_capping.pdf"))
print(
  plot_grid(
    ggplot(subset(capping,capping$group== "VN1"))+
      geom_tile(aes(
        x = variable  ,
        y=Amino_Acid,
        fill = value 
      ))+coord_equal()+
      theme_bw()+
      theme(text   = element_text(size = 6)) +
      
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "VN1 fixed scale "
      ),
    
    ggplot(delta[delta$group %in% c('VN1'),])+
      geom_tile(aes(
        x = group  ,
        y=Amino_Acid,
        fill = delta 
      ))+coord_equal()+
      theme_bw()+theme(text   = element_text(size = 6)) +
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F5AB53", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#4240F6", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(c(-4,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,4)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "VN1 delta "
      )
  )
)
dev.off()

######## VN2


pdf(paste0(path_plot,"Figure4/VN2_capping.pdf"))
print(
  plot_grid(
    ggplot(subset(capping,capping$group== "VN2"))+
      geom_tile(aes(
        x = variable  ,
        y=Amino_Acid,
        fill = value 
      ))+coord_equal()+
      theme_bw()+
      theme(text   = element_text(size = 6)) +
      
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "VN2 fixed scale "
      ),
    
    ggplot(delta[delta$group %in% c('VN2'),])+
      geom_tile(aes(
        x = group  ,
        y=Amino_Acid,
        fill = delta 
      ))+coord_equal()+
      theme_bw()+theme(text   = element_text(size = 6)) +
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F5AB53", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#4240F6", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(c(-4,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,4)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "VN2 delta "
      )
  )
)
dev.off()

######## VN3


pdf(paste0(path_plot,"Figure4/VN3_capping.pdf"))
print(
  plot_grid(
    ggplot(subset(capping,capping$group== "VN3"))+
      geom_tile(aes(
        x = variable  ,
        y=Amino_Acid,
        fill = value 
      ))+coord_equal()+
      theme_bw()+
      theme(text   = element_text(size = 6)) +
      
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "VN3 fixed scale "
      ),
    
    ggplot(delta[delta$group %in% c('VN3'),])+
      geom_tile(aes(
        x = group  ,
        y=Amino_Acid,
        fill = delta 
      ))+coord_equal()+
      theme_bw()+theme(text   = element_text(size = 6)) +
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F5AB53", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#4240F6", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(c(-4,-3.2,-2.4,-1.6,-0.8,0,0.8,1.6,2.4,3.2,4)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "VN3 delta "
      )
  )
)
dev.off()

#### Saturation Mutagenesis

# saturation mutagenesis 



IN1_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/IN1.xlsx"), sheet =1)
IN2_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/IN2.xlsx"), sheet =1)
LN1_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/LN1.xlsx"), sheet =1)
LN3_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/LN3.xlsx"), sheet =1)
VN1_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/VN1.xlsx"), sheet =1)
VN2_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/VN2.xlsx"), sheet =1)
VN3_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/VN3.xlsx"), sheet =1)


IN1_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/IN1.xlsx"), sheet =2)
IN2_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/IN2.xlsx"), sheet =2)
LN1_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/LN1.xlsx"), sheet =2)
LN3_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/LN3.xlsx"), sheet =2)
VN1_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/VN1.xlsx"), sheet =2)
VN2_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/VN2.xlsx"), sheet =2)
VN3_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/VN3.xlsx"), sheet =2)

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


### IN1
IN1_WT<-melt(IN1_WT, id = "Saturation")
IN1_WT$value<-as.numeric(IN1_WT$value)
IN1_WT$value<-round(IN1_WT$value,4)
IN1_WT$value_norm<-IN1_WT$value/(subset(IN1_WT,IN1_WT$Saturation == "N" & IN1_WT$variable == "N")[,c("value")])
IN1_WT$value_log2<-log2(IN1_WT$value_norm)
IN1_WT$group<-"IN1"


IN1_DT<-melt(IN1_DT, id = "Saturation")
IN1_DT$value<-as.numeric(IN1_DT$value)
IN1_DT$value<-round(IN1_DT$value,4)
IN1_DT$value_norm<-IN1_DT$value/(subset(IN1_DT,IN1_DT$Saturation == "N" & IN1_DT$variable == "N")[,c("value")])
IN1_DT$value_log2<-log2(IN1_DT$value_norm)
IN1_DT$group<-"IN1"


IN1_WT$identifier<-paste0(IN1_WT$Saturation,IN1_WT$variable)
IN1_DT$identifier<-paste0(IN1_DT$Saturation,IN1_DT$variable)
IN1<-merge(IN1_WT,IN1_DT, by = "identifier")
IN1<-IN1[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(IN1)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
              "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
IN1$Log2Delta<-log2(IN1$PSI_Das1/IN1$PSI_WT)

### IN2
IN2_WT<-melt(IN2_WT, id = "Saturation")
IN2_WT$value<-as.numeric(IN2_WT$value)
IN2_WT$value<-round(IN2_WT$value,4)
IN2_WT$value_norm<-IN2_WT$value/(subset(IN2_WT,IN2_WT$Saturation == "N" & IN2_WT$variable == "N")[,c("value")])
IN2_WT$value_log2<-log2(IN2_WT$value_norm)
IN2_WT$group<-"IN2"


IN2_DT<-melt(IN2_DT, id = "Saturation")
IN2_DT$value<-as.numeric(IN2_DT$value)
IN2_DT$value<-round(IN2_DT$value,4)
IN2_DT$value_norm<-IN2_DT$value/(subset(IN2_DT,IN2_DT$Saturation == "N" & IN2_DT$variable == "N")[,c("value")])
IN2_DT$value_log2<-log2(IN2_DT$value_norm)
IN2_DT$group<-"IN2"


IN2_WT$identifier<-paste0(IN2_WT$Saturation,IN2_WT$variable)
IN2_DT$identifier<-paste0(IN2_DT$Saturation,IN2_DT$variable)
IN2<-merge(IN2_WT,IN2_DT, by = "identifier")
IN2<-IN2[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(IN2)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
              "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
IN2$Log2Delta<-log2(IN2$PSI_Das1/IN2$PSI_WT)

### LN1

LN1_WT<-melt(LN1_WT, id = "Saturation")
LN1_WT$value<-as.numeric(LN1_WT$value)
LN1_WT$value<-round(LN1_WT$value,4)

LN1_WT$value_norm<-LN1_WT$value/(subset(LN1_WT,LN1_WT$Saturation == "N" & LN1_WT$variable == "N2")[,c("value")])
LN1_WT$value_log2<-log2(LN1_WT$value_norm)
LN1_WT$group<-"LN1"


LN1_DT<-melt(LN1_DT, id = "Saturation")
LN1_DT$value<-as.numeric(LN1_DT$value)
LN1_DT$value<-round(LN1_DT$value,4)
LN1_DT$value_norm<-LN1_DT$value/(subset(LN1_DT,LN1_DT$Saturation == "N" & LN1_DT$variable == "N2")[,c("value")])
LN1_DT$value_log2<-log2(LN1_DT$value_norm)
LN1_DT$group<-"LN1"


LN1_WT$identifier<-paste0(LN1_WT$Saturation,LN1_WT$variable)
LN1_DT$identifier<-paste0(LN1_DT$Saturation,LN1_DT$variable)
LN1<-merge(LN1_WT,LN1_DT, by = "identifier")
LN1<-LN1[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(LN1)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
              "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
LN1$Log2Delta<-log2(LN1$PSI_Das1/LN1$PSI_WT)

### LN3

LN3_WT<-melt(LN3_WT, id = "Saturation")
LN3_WT$value<-as.numeric(LN3_WT$value)
LN3_WT$value<-round(LN3_WT$value,4)
LN3_WT$value_norm<-LN3_WT$value/(subset(LN3_WT,LN3_WT$Saturation == "N" & LN3_WT$variable == "N")[,c("value")])
LN3_WT$value_log2<-log2(LN3_WT$value_norm)
LN3_WT$group<-"LN3"


LN3_DT<-melt(LN3_DT, id = "Saturation")
LN3_DT$value<-as.numeric(LN3_DT$value)
LN3_DT$value<-round(LN3_DT$value,4)
LN3_DT$value_norm<-LN3_DT$value/(subset(LN3_DT,LN3_DT$Saturation == "N" & LN3_DT$variable == "N")[,c("value")])
LN3_DT$value_log2<-log2(LN3_DT$value_norm)
LN3_DT$group<-"LN3"


LN3_WT$identifier<-paste0(LN3_WT$Saturation,LN3_WT$variable)
LN3_DT$identifier<-paste0(LN3_DT$Saturation,LN3_DT$variable)
LN3<-merge(LN3_WT,LN3_DT, by = "identifier")
LN3<-LN3[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(LN3)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
              "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
LN3$Log2Delta<-log2(LN3$PSI_Das1/LN3$PSI_WT)

### VN1

VN1_WT<-melt(VN1_WT, id = "Saturation")
VN1_WT$value<-as.numeric(VN1_WT$value)
VN1_WT$value<-round(VN1_WT$value,4)
VN1_WT$value_norm<-VN1_WT$value/(subset(VN1_WT,VN1_WT$Saturation == "N" & VN1_WT$variable == "N")[,c("value")])
VN1_WT$value_log2<-log2(VN1_WT$value_norm)
VN1_WT$group<-"VN1"


VN1_DT<-melt(VN1_DT, id = "Saturation")
VN1_DT$value<-as.numeric(VN1_DT$value)
VN1_DT$value<-round(VN1_DT$value,4)
VN1_DT$value_norm<-VN1_DT$value/(subset(VN1_DT,VN1_DT$Saturation == "N" & VN1_DT$variable == "N")[,c("value")])
VN1_DT$value_log2<-log2(VN1_DT$value_norm)
VN1_DT$group<-"VN1"


VN1_WT$identifier<-paste0(VN1_WT$Saturation,VN1_WT$variable)
VN1_DT$identifier<-paste0(VN1_DT$Saturation,VN1_DT$variable)
VN1<-merge(VN1_WT,VN1_DT, by = "identifier")
VN1<-VN1[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(VN1)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
              "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
VN1$Log2Delta<-log2(VN1$PSI_Das1/VN1$PSI_WT)

### VN2

VN2_WT<-melt(VN2_WT, id = "Saturation")
VN2_WT$value<-as.numeric(VN2_WT$value)
VN2_WT$value<-round(VN2_WT$value,4)
VN2_WT$value_norm<-VN2_WT$value/(subset(VN2_WT,VN2_WT$Saturation == "N" & VN2_WT$variable == "N")[,c("value")])
VN2_WT$value_log2<-log2(VN2_WT$value_norm)
VN2_WT$group<-"VN2"


VN2_DT<-melt(VN2_DT, id = "Saturation")
VN2_DT$value<-as.numeric(VN2_DT$value)
VN2_DT$value<-round(VN2_DT$value,4)
VN2_DT$value_norm<-VN2_DT$value/(subset(VN2_DT,VN2_DT$Saturation == "N" & VN2_DT$variable == "N")[,c("value")])
VN2_DT$value_log2<-log2(VN2_DT$value_norm)
VN2_DT$group<-"VN2"


VN2_WT$identifier<-paste0(VN2_WT$Saturation,VN2_WT$variable)
VN2_DT$identifier<-paste0(VN2_DT$Saturation,VN2_DT$variable)
VN2<-merge(VN2_WT,VN2_DT, by = "identifier")
VN2<-VN2[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(VN2)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
              "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
VN2$Log2Delta<-log2(VN2$PSI_Das1/VN2$PSI_WT)

### VN3

VN3_WT<-melt(VN3_WT, id = "Saturation")
VN3_WT$value<-as.numeric(VN3_WT$value)
VN3_WT$value<-round(VN3_WT$value,4)
VN3_WT$value_norm<-VN3_WT$value/(subset(VN3_WT,VN3_WT$Saturation == "N" & VN3_WT$variable == "N")[,c("value")])
VN3_WT$value_log2<-log2(VN3_WT$value_norm)
VN3_WT$group<-"VN3"


VN3_DT<-melt(VN3_DT, id = "Saturation")
VN3_DT$value<-as.numeric(VN3_DT$value)
VN3_DT$value<-round(VN3_DT$value,4)
VN3_DT$value_norm<-VN3_DT$value/(subset(VN3_DT,VN3_DT$Saturation == "N" & VN3_DT$variable == "N")[,c("value")])
VN3_DT$value_log2<-log2(VN3_DT$value_norm)
VN3_DT$group<-"VN3"


VN3_WT$identifier<-paste0(VN3_WT$Saturation,VN3_WT$variable)
VN3_DT$identifier<-paste0(VN3_DT$Saturation,VN3_DT$variable)
VN3<-merge(VN3_WT,VN3_DT, by = "identifier")
VN3<-VN3[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(VN3)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
              "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
VN3$Log2Delta<-log2(VN3$PSI_Das1/VN3$PSI_WT)


####

data<-as.data.frame(rbind(IN1,IN2,LN1,LN3,VN1,VN2,VN3))
data$Saturation<-factor(data$Saturation,
                        levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
data$PSI_WT<-ifelse(data$PSI_WT > 0.75,0.75,data$PSI_WT)
data$PSI_Das1<-ifelse(data$PSI_Das1 > 0.75,0.75,data$PSI_Das1)
### in1
IN1<-subset(data,data$group== "IN1")
pdf(paste0(path_plot,"Figure4/IN1_SaturationMutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(IN1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("L1" = "L",
                                "V1" = "V",
                                "V2" = "V",
                                "Y" = "Y",
                                "L2" = "L",
                                "S" = "S",
                                "C1" = "C",
                                "V3" = "V",
                                "A" = "A",
                                "C2" = "C",
                                "I" = "I",
                                "N" = "N"))+
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "IN1 WT "
      ),
    
    ggplot(IN1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("L1" = "L",
                                "V1" = "V",
                                "V2" = "V",
                                "Y" = "Y",
                                "L2" = "L",
                                "S" = "S",
                                "C1" = "C",
                                "V3" = "V",
                                "A" = "A",
                                "C2" = "C",
                                "I" = "I",
                                "N" = "N"))+
      scale_fill_gradientn("PSI",colours=c(
        "#00A650", 
        "#42B977", 
        "#6DC590", 
        "#98D2AA", 
        "#C4DEC4", 
        "#F0EBDE",
        "#EAC5CF", 
        "#E4A0BF", 
        "#DD7AB0", 
        "#D755A0", 
        "#CE1C89"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
      limits = c(0,0.75))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "IN1 Das1 "
      ),
    
    
    
    
    ggplot(IN1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("L1" = "L",
                                "V1" = "V",
                                "V2" = "V",
                                "Y" = "Y",
                                "L2" = "L",
                                "S" = "S",
                                "C1" = "C",
                                "V3" = "V",
                                "A" = "A",
                                "C2" = "C",
                                "I" = "I",
                                "N" = "N"))+
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F5AB53", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#4240F6", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(c(-6,-3.00, -2.33, -1.66, -0.99 ,0  ,0.99 , 1.66 , 2.33 , 3,6)),
      limits = c(-6,6))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "IN1 Delta "
      )
    
  ))

dev.off()

## in2
IN2<-subset(data,data$group== "IN2")
IN2$variable<-factor(IN2$variable,
                     levels = c("K1","K2","H","K3","K4","G","Q1","K5","Q2","K6","I","N"))
pdf(paste0(path_plot,"Figure4/IN2_SaturationMutagenesis.pdf"))
print(plot_grid(
  
  
  
  ggplot(IN2)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_WT 
    ))+coord_equal()+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("K1" = "K",
                              "K2" = "K",
                              "H" = "H",
                              "K3" = "K",
                              "K4" = "K",
                              "G" = "G",
                              "Q1" = "Q",
                              "K5" = "K",
                              "Q2" = "Q",
                              "K6" = "K",
                              "I" = "I",
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "IN2 WT "
    ),
  
  ggplot(IN2)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_Das1 
    ))+coord_equal()+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("K1" = "K",
                              "K2" = "K",
                              "H" = "H",
                              "K3" = "K",
                              "K4" = "K",
                              "G" = "G",
                              "Q1" = "Q",
                              "K5" = "K",
                              "Q2" = "Q",
                              "K6" = "K",
                              "I" = "I",
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "IN2 Das1 "
    ),
  
  
  
  ggplot(IN2)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = Log2Delta 
    ))+coord_equal()+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("K1" = "K",
                              "K2" = "K",
                              "H" = "H",
                              "K3" = "K",
                              "K4" = "K",
                              "G" = "G",
                              "Q1" = "Q",
                              "K5" = "K",
                              "Q2" = "Q",
                              "K6" = "K",
                              "I" = "I",
                              "N" = "N"))+
    scale_fill_gradientn("Log2(Delta)",colours=c(
      "#F7931E", 
      "#F5AB53", 
      "#F4BB75", 
      "#F3CB98", 
      "#F1DBBB", 
      "#F0EBDE", 
      "#C4C0E4", 
      "#9996EA", 
      "#6D6BF0", 
      "#4240F6", 
      "#0000FF"
    ), na.value = "grey98",
    values = scales::rescale(c(-6,-3.00, -2.33, -1.66, -0.99 ,0  ,0.99 , 1.66 , 2.33 , 3,6)),
    limits = c(-6,6))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "IN2 Delta "
    )
  
))
dev.off()

###ln1
LN1<-subset(data,data$group== "LN1")
LN1$variable<-factor(LN1$variable,
                     levels = c("P","D","K1","R","N1","H","L1","I","K2","L2","L3","N2"))
pdf(paste0(path_plot,"Figure4/LN1_SaturationMutagenesis.pdf"))

print(plot_grid(
  
  
  
  ggplot(LN1)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_WT 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("P" = "P",
                              "D" = "D",
                              "K1" = "K",
                              "R" = "R",
                              "N1" = "N",
                              "H" = "H",
                              "L1" = "L",
                              "I" = "I",
                              "K2" = "K",
                              "L2" = "L",
                              "L3" = "L",
                              "N2" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "LN1 WT "
    ),
  
  ggplot(LN1)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_Das1 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("P" = "P",
                              "D" = "D",
                              "K1" = "K",
                              "R" = "R",
                              "N1" = "N",
                              "H" = "H",
                              "L1" = "L",
                              "I" = "I",
                              "K2" = "K",
                              "L2" = "L",
                              "L3" = "L",
                              "N2" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "LN1 Das1 "
    ),
  
  
  
  ggplot(LN1)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = Log2Delta 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("P" = "P",
                              "D" = "D",
                              "K1" = "K",
                              "R" = "R",
                              "N1" = "N",
                              "H" = "H",
                              "L1" = "L",
                              "I" = "I",
                              "K2" = "K",
                              "L2" = "L",
                              "L3" = "L",
                              "N2" = "N"))+
    scale_fill_gradientn("Log2(Delta)",colours=c(
      "#F7931E", 
      "#F5AB53", 
      "#F4BB75", 
      "#F3CB98", 
      "#F1DBBB", 
      "#F0EBDE", 
      "#C4C0E4", 
      "#9996EA", 
      "#6D6BF0", 
      "#4240F6", 
      "#0000FF"
    ), na.value = "grey98",
    values = scales::rescale(c(-6,-3.00, -2.33, -1.66, -0.99 ,0  ,0.99 , 1.66 , 2.33 , 3,6)),
    limits = c(-6,6))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "LN1 Delta "
    )
  
))
dev.off()
###LN3
LN3<-subset(data,data$group== "LN3")
LN3$variable<-factor(LN3$variable,
                     levels = c("R","L1",	"I",	"Q1","P1",	"L2",	"P2",	"V",	"Q2",	"F",	"L3",	"N"))
pdf(paste0(path_plot,"Figure4/LN3_SaturationMutagenesis.pdf"))
print(plot_grid(
  
  
  
  ggplot(LN3)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_WT 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("R" = "R",
                              "L1" = "L",	
                              "I" = "I",	
                              "Q1" = "Q",
                              "P1" = "P",	
                              "L2" = "L",	
                              "P2"= "P",	
                              "V" = "V",	
                              "Q2" = "Q",	
                              "F" = "F",	
                              "L3" = "F",	
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "LN3 WT "
    ),
  
  ggplot(LN3)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_Das1 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("R" = "R",
                              "L1" = "L",	
                              "I" = "I",	
                              "Q1" = "Q",
                              "P1" = "P",	
                              "L2" = "L",	
                              "P2"= "P",	
                              "V" = "V",	
                              "Q2" = "Q",	
                              "F" = "F",	
                              "L3" = "F",	
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "LN3 Das1 "
    ),
  
  
  
  ggplot(LN3)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = Log2Delta 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("R" = "R",
                              "L1" = "L",	
                              "I" = "I",	
                              "Q1" = "Q",
                              "P1" = "P",	
                              "L2" = "L",	
                              "P2"= "P",	
                              "V" = "V",	
                              "Q2" = "Q",	
                              "F" = "F",	
                              "L3" = "F",	
                              "N" = "N"))+
    scale_fill_gradientn("Log2(Delta)",colours=c(
      "#F7931E", 
      "#F5AB53", 
      "#F4BB75", 
      "#F3CB98", 
      "#F1DBBB", 
      "#F0EBDE", 
      "#C4C0E4", 
      "#9996EA", 
      "#6D6BF0", 
      "#4240F6", 
      "#0000FF"
    ), na.value = "grey98",
    values = scales::rescale(c(-6,-3.00, -2.33, -1.66, -0.99 ,0  ,0.99 , 1.66 , 2.33 , 3,6)),
    limits = c(-6,6))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "LN3 Delta "
    )
  
))
dev.off()

###VN1

VN1<-subset(data,data$group== "VN1")
VN1$variable<-factor(VN1$variable,
                     levels = c("C",	"L1",	"F1",	"H",	"Q",	"F2",	"I",	"T",	"L2",	"L3",	"V",	"N"))
pdf(paste0(path_plot,"Figure4/VN1_SaturationMutagenesis.pdf"))
print(plot_grid(
  
  
  
  ggplot(VN1)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_WT 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("C" = "C",
                              "L1" = "L",
                              "F1" = "F",
                              "H" = "H",	
                              "Q" = "Q",	
                              "F2" = "F",	
                              "I" = "I",	
                              "T" = "T",	
                              "L2" = "L",	
                              "L3" = "L",	
                              "V" = "V",	
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN1 WT "
    ),
  
  ggplot(VN1)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_Das1 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("C" = "C",
                              "L1" = "L",
                              "F1" = "F",
                              "H" = "H",	
                              "Q" = "Q",	
                              "F2" = "F",	
                              "I" = "I",	
                              "T" = "T",	
                              "L2" = "L",	
                              "L3" = "L",	
                              "V" = "V",	
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN1 Das1 "
    ),
  
  
  
  ggplot(VN1)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = Log2Delta 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("C" = "C",
                              "L1" = "L",
                              "F1" = "F",
                              "H" = "H",	
                              "Q" = "Q",	
                              "F2" = "F",	
                              "I" = "I",	
                              "T" = "T",	
                              "L2" = "L",	
                              "L3" = "L",	
                              "V" = "V",	
                              "N" = "N"))+
    scale_fill_gradientn("Log2(Delta)",colours=c(
      "#F7931E", 
      "#F5AB53", 
      "#F4BB75", 
      "#F3CB98", 
      "#F1DBBB", 
      "#F0EBDE", 
      "#C4C0E4", 
      "#9996EA", 
      "#6D6BF0", 
      "#4240F6", 
      "#0000FF"
    ), na.value = "grey98",
    values = scales::rescale(c(-6,-3.00, -2.33, -1.66, -0.99 ,0  ,0.99 , 1.66 , 2.33 , 3,6)),
    limits = c(-6,6))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN1 Delta "
    )
  
  
  
))
dev.off()

### vn2
VN2<-subset(data,data$group== "VN2")
VN2$variable<-factor(VN2$variable,
                     levels = c("F1",  "S" ,     "M"   ,   "R1" , "W1",  "T" ,     "F2" , "R2" , "Y"    ,  "W2" ,"V" ,     "N"))
pdf(paste0(path_plot,"Figure4/VN2_SaturationMutagenesis.pdf"))
print(plot_grid(
  
  
  
  ggplot(VN2)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_WT 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("F1" = "F",  
                              "S"  = "S",     
                              "M" ="M",   
                              "R1"  = "R", 
                              "W1" = "W",  
                              "T" = "T",     
                              "F2"  = "F", 
                              "R2" = "R", 
                              "Y" ="Y",  
                              "W2" = "W",
                              "V" ="V",     
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN2 WT "
    ),
  
  ggplot(VN2)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_Das1 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("F1" = "F",  
                              "S"  = "S",     
                              "M" ="M",   
                              "R1"  = "R", 
                              "W1" = "W",  
                              "T" = "T",     
                              "F2"  = "F", 
                              "R2" = "R", 
                              "Y" ="Y",  
                              "W2" = "W",
                              "V" ="V",     
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN2 Das1 "
    ),
  
  
  
  ggplot(VN2)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = Log2Delta 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("F1" = "F",  
                              "S"  = "S",     
                              "M" ="M",   
                              "R1"  = "R", 
                              "W1" = "W",  
                              "T" = "T",     
                              "F2"  = "F", 
                              "R2" = "R", 
                              "Y" ="Y",  
                              "W2" = "W",
                              "V" ="V",     
                              "N" = "N"))+
    scale_fill_gradientn("Log2(Delta)",colours=c(
      "#F7931E", 
      "#F5AB53", 
      "#F4BB75", 
      "#F3CB98", 
      "#F1DBBB", 
      "#F0EBDE", 
      "#C4C0E4", 
      "#9996EA", 
      "#6D6BF0", 
      "#4240F6", 
      "#0000FF"
    ), na.value = "grey98",
    values = scales::rescale(c(-6,-3.00, -2.33, -1.66, -0.99 ,0  ,0.99 , 1.66 , 2.33 , 3,6)),
    limits = c(-6,6))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN2 Delta "
    )
  
))
dev.off()

###VN3
VN3<-subset(data,data$group== "VN3")
VN3$variable<-factor(VN3$variable,
                     levels = c("I1"	,"I2"	,"G",	"C1",	"I3",	"A",	"I4",	"C2",	"L",	"M",	"V",	"N"))
pdf(paste0(path_plot,"Figure4/VN3_SaturationMutagenesis.pdf"))
print(plot_grid(
  
  
  
  ggplot(VN3)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_WT 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("I1" = "I"	,
                              "I2" = "I"	,
                              "G" = "G",	
                              "C1" = "C",	
                              "I3" = "I",	
                              "A" = "A",	
                              "I4" = "I",	
                              "C2" = "C",	
                              "L" = "L",	
                              "M" = "M",	
                              "V" = "V",	
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN3 WT "
    ),
  
  ggplot(VN3)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = PSI_Das1 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("I1" = "I"	,
                              "I2" = "I"	,
                              "G" = "G",	
                              "C1" = "C",	
                              "I3" = "I",	
                              "A" = "A",	
                              "I4" = "I",	
                              "C2" = "C",	
                              "L" = "L",	
                              "M" = "M",	
                              "V" = "V",	
                              "N" = "N"))+
    scale_fill_gradientn("PSI",colours=c(
      "#00A650", 
      "#42B977", 
      "#6DC590", 
      "#98D2AA", 
      "#C4DEC4", 
      "#F0EBDE",
      "#EAC5CF", 
      "#E4A0BF", 
      "#DD7AB0", 
      "#D755A0", 
      "#CE1C89"
    ), na.value = "grey98",
    values = scales::rescale(round(seq(0,0.75, 0.75/12),2)),
    limits = c(0,0.75))+
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN3 Das1 "
    ),
  
  
  
  ggplot(VN3)+
    geom_tile(aes(
      x = variable  ,
      y=Saturation,
      fill = Log2Delta 
    ))+
    theme_bw()+
    theme(text   = element_text(size = 8)) +
    scale_x_discrete(labels=c("I1" = "I"	,
                              "I2" = "I"	,
                              "G" = "G",	
                              "C1" = "C",	
                              "I3" = "I",	
                              "A" = "A",	
                              "I4" = "I",	
                              "C2" = "C",	
                              "L" = "L",	
                              "M" = "M",	
                              "V" = "V",	
                              "N" = "N"))+
    scale_fill_gradientn("Log2(Delta)",colours=c(
      "#F7931E", 
      "#F5AB53", 
      "#F4BB75", 
      "#F3CB98", 
      "#F1DBBB", 
      "#F0EBDE", 
      "#C4C0E4", 
      "#9996EA", 
      "#6D6BF0", 
      "#4240F6", 
      "#0000FF"
    ), na.value = "grey98",
    values = scales::rescale(c(-6,-3.00, -2.33, -1.66, -0.99 ,0  ,0.99 , 1.66 , 2.33 , 3,6)),
    limits = c(-6,6))+
    
    
    ylab("Mutations")+
    xlab("Saturated Sequence ")+
    labs(title =  "VN3 Delta "
    )
  
))
dev.off()


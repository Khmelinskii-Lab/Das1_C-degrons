##### for project : Analysis of degron at C terminal
##### sub text:capping and saturation mutagenesis
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")
path_model<-paste0(getwd(),"/model")




R1<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =1)
names(R1)<-c("Amino_Acid","WT","Doa10")
R2<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =2)
names(R2)<-c("Amino_Acid","WT","Doa10")
R3<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =3)
names(R3)<-c("Amino_Acid","WT","Doa10")
R4<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =4)
names(R4)<-c("Amino_Acid","WT","Doa10")
R5<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =5)
names(R5)<-c("Amino_Acid","WT","Doa10")
Atg1C<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =6)
names(Atg1C)<-c("Amino_Acid","WT","Das1")
Rpa12C<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =7)
names(Rpa12C)<-c("Amino_Acid","WT","Das1")
W1<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =8)
names(W1)<-c("Amino_Acid","WT","Das1")
W2<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =9)
names(W2)<-c("Amino_Acid","WT","Das1")
W3<-readxl::read_xlsx(paste0(path_input,"/capping/capping2.xlsx"), sheet =10)
names(W3)<-c("Amino_Acid","WT","Das1")

R1_melt<-melt(R1, id = "Amino_Acid")
R1_melt$value_log2<-log2(R1_melt$value)
R1_melt$group<-"R1"

R2_melt<-melt(R2, id = "Amino_Acid")
R2_melt$value_log2<-log2(R2_melt$value)
R2_melt$group<-"R2"

R3_melt<-melt(R3, id = "Amino_Acid")
R3_melt$value_log2<-log2(R3_melt$value)
R3_melt$group<-"R3"

R4_melt<-melt(R4, id = "Amino_Acid")
R4_melt$value_log2<-log2(R4_melt$value)
R4_melt$group<-"R4"

R5_melt<-melt(R5, id = "Amino_Acid")
R5_melt$value_log2<-log2(R5_melt$value)
R5_melt$group<-"R5"

Atg1C_melt<-melt(Atg1C, id = "Amino_Acid")
Atg1C_melt$value_log2<-log2(Atg1C_melt$value)
Atg1C_melt$group<-"Atg1C"

Rpa12C_melt<-melt(Rpa12C, id = "Amino_Acid")
Rpa12C_melt$value_log2<-log2(Rpa12C_melt$value)
Rpa12C_melt$group<-"Rpa12C"

W1_melt<-melt(W1, id = "Amino_Acid")
W1_melt$value_log2<-log2(W1_melt$value)
W1_melt$group<-"W1"

W2_melt<-melt(W2, id = "Amino_Acid")
W2_melt$value_log2<-log2(W2_melt$value)
W2_melt$group<-"W2"

W3_melt<-melt(W3, id = "Amino_Acid")
W3_melt$value_log2<-log2(W3_melt$value)
W3_melt$group<-"W3"

capping<-as.data.frame(rbind(R1_melt,R2_melt,R3_melt,R4_melt,R5_melt,Atg1C_melt,Rpa12C_melt,W1_melt,W2_melt,W3_melt))
R1$delta<-log2(R1$Doa10/R1$WT)
R1$group<-"R1"
R2$delta<-log2(R2$Doa10/R2$WT)
R2$group<-"R2"
R3$delta<-log2(R3$Doa10/R3$WT)
R3$group<-"R3"
R4$delta<-log2(R4$Doa10/R4$WT)
R4$group<-"R4"
R5$delta<-log2(R5$Doa10/R5$WT)
R5$group<-"R5"
Atg1C$delta<-log2(Atg1C$Das1/Atg1C$WT)
Atg1C$group<-"Atg1C"
Rpa12C$delta<-log2(Rpa12C$Das1/Rpa12C$WT)
Rpa12C$group<-"Rpa12C"
W1$delta<-log2(W1$Das1/W1$WT)
W1$group<-"W1"
W2$delta<-log2(W2$Das1/W2$WT)
W2$group<-"W2"
W3$delta<-log2(W3$Das1/W3$WT)
W3$group<-"W3"

delta<-(rbind(R1,R2,R3,R4,R5))
capping$Amino_Acid<-factor(capping$Amino_Acid,
                           levels = rev(c("Stop","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))

delta$Amino_Acid<-factor(delta$Amino_Acid,
                         levels = rev(c("Stop","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))



pdf(paste0(path_plot,"FigureS5/Doa10_R1.pdf"))
print(
  plot_grid(
    ggplot(subset(capping,capping$group== "R1"))+
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
      labs(title =  "R1 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('R1'),])+
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
      values = scales::rescale(round(seq(-4,4,8/10),2)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R1 delta "
      )
  )
)
dev.off()

pdf(paste0(path_plot,"FigureS5/Doa10_R2.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "R2"))+
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
      labs(title =  "R2 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('R2'),])+
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
      values = scales::rescale(round(seq(-4,4,8/10),2)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R2 delta "
      )
  )
)
dev.off()


pdf(paste0(path_plot,"FigureS5/Doa10_R3.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "R3"))+
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
      labs(title =  "R3 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('R3'),])+
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
      values = scales::rescale(round(seq(-4,4,8/10),2)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R3 delta "
      )
  )
)
dev.off()

# R4 and R5 are swapped for easier writing purposes
pdf(paste0(path_plot,"FigureS5/Doa10_R5.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "R4"))+
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
      labs(title =  "R4 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('R4'),])+
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
      values = scales::rescale(round(seq(-4,4,8/10),2)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R4 delta "
      )
  )
)
dev.off()

pdf(paste0(path_plot,"FigureS5/Doa10_R4.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "R5"))+
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
      labs(title =  "R5 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('R5'),])+
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
      values = scales::rescale(round(seq(-4,4,8/10),2)),
      limits = c(-4,4))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R5 delta "
      )
  )
)
dev.off()



delta<-(rbind(Atg1C,Rpa12C,W1,W2,W3))
capping$Amino_Acid<-factor(capping$Amino_Acid,
                           levels = rev(c("Stop","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))

delta$Amino_Acid<-factor(delta$Amino_Acid,
                         levels = rev(c("Stop","A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))
# Rpa12C
pdf(paste0(path_plot,"Figure6/Das_Rpa12C.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "Rpa12C"))+
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
      labs(title =  "Rpa12C PSI "
      ),
    
    ggplot(delta[delta$group %in% c('Rpa12C'),])+
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
      values = scales::rescale(round(seq(-3,3,6/10),2)),
      limits = c(-3,3))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "Rpa12C delta "
      )
  )
)
dev.off()

# Atg1C
pdf(paste0(path_plot,"Figure6/Das_Atg1C.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "Atg1C"))+
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
      labs(title =  "Atg1C PSI "
      ),
    
    ggplot(delta[delta$group %in% c('Atg1C'),])+
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
      values = scales::rescale(round(seq(-3,3,6/10),2)),
      limits = c(-3,3))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "Atg1C delta "
      )
  )
)
dev.off()

#RPL12c


# W1
pdf(paste0(path_plot,"Figure5/Das_W1.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "W1"))+
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
      labs(title =  "W1 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('W1'),])+
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
      values = scales::rescale(round(seq(-3,3,6/10),2)),
      limits = c(-3,3))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "W1 delta "
      )
  )
)
dev.off()

pdf(paste0(path_plot,"Figure5/Das_W2.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "W2"))+
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
      labs(title =  "W2 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('W2'),])+
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
      values = scales::rescale(round(seq(-3,3,6/10),2)),
      limits = c(-3,3))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "W2 delta "
      )
  )
)
dev.off()




pdf(paste0(path_plot,"Figure5/Das_W3.pdf"))

print(
  plot_grid(
    ggplot(subset(capping,capping$group== "W3"))+
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
      labs(title =  "W3 PSI "
      ),
    
    ggplot(delta[delta$group %in% c('W3'),])+
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
      values = scales::rescale(round(seq(-3,3,6/10),2)),
      limits = c(-3,3))+
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "W3 delta "
      )
  )
)
dev.off()

################# saturation mutagenesis

# saturation mutagenesis 



R1_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R1.xlsx"), sheet =1)
R2_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R2.xlsx"), sheet =1)
R3_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R3.xlsx"), sheet =1)
R4_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R4.xlsx"), sheet =1)
R5_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R5.xlsx"), sheet =1)
Atg1C_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/Atg1C.xlsx"), sheet =1)
Rpa12C_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/Rpa12C.xlsx"), sheet =1)
W1_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/W1.xlsx"), sheet =1)
W2_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/W2.xlsx"), sheet =1)
W3_WT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/W3.xlsx"), sheet =1)


R1_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R1.xlsx"), sheet =2)
R2_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R2.xlsx"), sheet =2)
R3_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R3.xlsx"), sheet =2)
R4_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R4.xlsx"), sheet =2)
R5_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/R5.xlsx"), sheet =2)
Atg1C_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/Atg1C.xlsx"), sheet =2)
Rpa12C_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/Rpa12C.xlsx"), sheet =2)
W1_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/W1.xlsx"), sheet =2)
W2_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/W2.xlsx"), sheet =2)
W3_DT<-readxl::read_xlsx(paste0(path_input,"/Saturation_mutagenesis/W3.xlsx"), sheet =2)


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


### R1
R1_WT<-melt(R1_WT, id = "Saturation")
R1_WT$value<-as.numeric(R1_WT$value)
R1_WT$value<-round(R1_WT$value,4)
R1_WT$value_norm<-R1_WT$value/(subset(R1_WT,R1_WT$Saturation == "R" & R1_WT$variable == "R")[,c("value")])
R1_WT$value_log2<-log2(R1_WT$value_norm)
R1_WT$group<-"R1"


R1_DT<-melt(R1_DT, id = "Saturation")
R1_DT$value<-as.numeric(R1_DT$value)
R1_DT$value<-round(R1_DT$value,4)
R1_DT$value_norm<-R1_DT$value/(subset(R1_DT,R1_DT$Saturation == "R" & R1_DT$variable == "R")[,c("value")])
R1_DT$value_log2<-log2(R1_DT$value_norm)
R1_DT$group<-"R1"


R1_WT$identifier<-paste0(R1_WT$Saturation,R1_WT$variable)
R1_DT$identifier<-paste0(R1_DT$Saturation,R1_DT$variable)
R1<-merge(R1_WT,R1_DT, by = "identifier")
R1<-R1[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(R1)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
R1$Log2Delta<-log2(R1$PSI_Das1/R1$PSI_WT)

### R2
R2_WT<-melt(R2_WT, id = "Saturation")
R2_WT$value<-as.numeric(R2_WT$value)
R2_WT$value<-round(R2_WT$value,4)
R2_WT$value_norm<-R2_WT$value/(subset(R2_WT,R2_WT$Saturation == "R" & R2_WT$variable == "R3")[,c("value")])
R2_WT$value_log2<-log2(R2_WT$value_norm)
R2_WT$group<-"R2"


R2_DT<-melt(R2_DT, id = "Saturation")
R2_DT$value<-as.numeric(R2_DT$value)
R2_DT$value<-round(R2_DT$value,4)
R2_DT$value_norm<-R2_DT$value/(subset(R2_DT,R2_DT$Saturation == "R" & R2_DT$variable == "R3")[,c("value")])
R2_DT$value_log2<-log2(R2_DT$value_norm)
R2_DT$group<-"R2"


R2_WT$identifier<-paste0(R2_WT$Saturation,R2_WT$variable)
R2_DT$identifier<-paste0(R2_DT$Saturation,R2_DT$variable)
R2<-merge(R2_WT,R2_DT, by = "identifier")
R2<-R2[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(R2)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
R2$Log2Delta<-log2(R2$PSI_Das1/R2$PSI_WT)

### R3

R3_WT<-melt(R3_WT, id = "Saturation")
R3_WT$value<-as.numeric(R3_WT$value)
R3_WT$value<-round(R3_WT$value,4)

R3_WT$value_norm<-R3_WT$value/(subset(R3_WT,R3_WT$Saturation == "R" & R3_WT$variable == "R")[,c("value")])
R3_WT$value_log2<-log2(R3_WT$value_norm)
R3_WT$group<-"R3"


R3_DT<-melt(R3_DT, id = "Saturation")
R3_DT$value<-as.numeric(R3_DT$value)
R3_DT$value<-round(R3_DT$value,4)
R3_DT$value_norm<-R3_DT$value/(subset(R3_DT,R3_DT$Saturation == "R" & R3_DT$variable == "R")[,c("value")])
R3_DT$value_log2<-log2(R3_DT$value_norm)
R3_DT$group<-"R3"


R3_WT$identifier<-paste0(R3_WT$Saturation,R3_WT$variable)
R3_DT$identifier<-paste0(R3_DT$Saturation,R3_DT$variable)
R3<-merge(R3_WT,R3_DT, by = "identifier")
R3<-R3[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(R3)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
R3$Log2Delta<-log2(R3$PSI_Das1/R3$PSI_WT)

### R4

R4_WT<-melt(R4_WT, id = "Saturation")
R4_WT$value<-as.numeric(R4_WT$value)
R4_WT$value<-round(R4_WT$value,4)
R4_WT$value_norm<-R4_WT$value/(subset(R4_WT,R4_WT$Saturation == "I" & R4_WT$variable == "I2")[,c("value")])
R4_WT$value_log2<-log2(R4_WT$value_norm)
R4_WT$group<-"R4"


R4_DT<-melt(R4_DT, id = "Saturation")
R4_DT$value<-as.numeric(R4_DT$value)
R4_DT$value<-round(R4_DT$value,4)
R4_DT$value_norm<-R4_DT$value/(subset(R4_DT,R4_DT$Saturation == "I" & R4_DT$variable == "I2")[,c("value")])
R4_DT$value_log2<-log2(R4_DT$value_norm)
R4_DT$group<-"R4"


R4_WT$identifier<-paste0(R4_WT$Saturation,R4_WT$variable)
R4_DT$identifier<-paste0(R4_DT$Saturation,R4_DT$variable)
R4<-merge(R4_WT,R4_DT, by = "identifier")
R4<-R4[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(R4)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
R4$Log2Delta<-log2(R4$PSI_Das1/R4$PSI_WT)

### R5

R5_WT<-melt(R5_WT, id = "Saturation")
R5_WT$value<-as.numeric(R5_WT$value)
R5_WT$value<-round(R5_WT$value,4)
R5_WT$value_norm<-R5_WT$value/(subset(R5_WT,R5_WT$Saturation == "R" & R5_WT$variable == "R")[,c("value")])
R5_WT$value_log2<-log2(R5_WT$value_norm)
R5_WT$group<-"R5"


R5_DT<-melt(R5_DT, id = "Saturation")
R5_DT$value<-as.numeric(R5_DT$value)
R5_DT$value<-round(R5_DT$value,4)
R5_DT$value_norm<-R5_DT$value/(subset(R5_DT,R5_DT$Saturation == "R" & R5_DT$variable == "R")[,c("value")])
R5_DT$value_log2<-log2(R5_DT$value_norm)
R5_DT$group<-"R5"


R5_WT$identifier<-paste0(R5_WT$Saturation,R5_WT$variable)
R5_DT$identifier<-paste0(R5_DT$Saturation,R5_DT$variable)
R5<-merge(R5_WT,R5_DT, by = "identifier")
R5<-R5[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y","group.y" )]
names(R5)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)

R5$Log2Delta<-log2(R5$PSI_Das1/R5$PSI_WT)












####

data<-as.data.frame(rbind(R1,R2,R3,R4,R5))
data$Saturation<-factor(data$Saturation,
                        levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))

R1<-subset(data,data$group== "R1")

pdf(paste0(path_plot,"FigureS5/Doa10_R1_Saturation_Mutagenesis.pdf"))
#pdf("Y:/lab data/susmitha/edwin/for_paper/new_data2/SaturationMutagenesis/Doa10_R1.pdf")
print(
  plot_grid(
    
    ggplot(R1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("A" = "A",
                                "I" = "I",
                                "T" = "T",
                                "L1" = "L",
                                "P1" = "P",
                                "M" = "M",
                                "C" = "C",
                                "P2" = "P",
                                "L2" = "L",
                                "L3" = "L",
                                "L4" = "L",
                                "R" = "R"))+
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
      labs(title =  "R1 WT "
      ),
    
    ggplot(R1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("A" = "A",
                                "I" = "I",
                                "T" = "T",
                                "L1" = "L",
                                "P1" = "P",
                                "M" = "M",
                                "C" = "C",
                                "P2" = "P",
                                "L2" = "L",
                                "L3" = "L",
                                "L4" = "L",
                                "R" = "R"))+
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
      labs(title =  "R1 Doa10 "
      ),
    
    
    
    
    ggplot(R1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("A" = "A",
                                "I" = "I",
                                "T" = "T",
                                "L1" = "L",
                                "P1" = "P",
                                "M" = "M",
                                "C" = "C",
                                "P2" = "P",
                                "L2" = "L",
                                "L3" = "L",
                                "L4" = "L",
                                "R" = "R"))+
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
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R1 Delta "
      )
    
  ))



dev.off()

pdf(paste0(path_plot,"FigureS5/Doa10_R2_Saturation_Mutagenesis.pdf"))

print(
  plot_grid(
    
    ggplot(R2)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("R1" = "R",
                                "S" = "S",
                                "L1" = "L",
                                "R2" = "R",
                                "V" = "V",
                                "L2" = "L",
                                "L3" = "L",
                                "P" = "P",
                                "L4" = "L",
                                "Y" = "Y",
                                "F" = "F",
                                "R3" = "R"))+
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
      labs(title =  "R2 WT " ),
    
    ggplot(R2)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("R1" = "R",
                                "S" = "S",
                                "L1" = "L",
                                "R2" = "R",
                                "V" = "V",
                                "L2" = "L",
                                "L3" = "L",
                                "P" = "P",
                                "L4" = "L",
                                "Y" = "Y",
                                "F" = "F",
                                "R3" = "R"))+
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
      labs(title =  "R2 Doa10 "),
    
    
    
    
    ggplot(R2)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("R1" = "R",
                                "S" = "S",
                                "L1" = "L",
                                "R2" = "R",
                                "V" = "V",
                                "L2" = "L",
                                "L3" = "L",
                                "P" = "P",
                                "L4" = "L",
                                "Y" = "Y",
                                "F" = "F",
                                "R3" = "R"))+
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R2 Delta ")
    
  ))

dev.off()

R3<-subset(data,data$group== "R3")
R3$variable<-factor(R3$variable,
                    levels = c("C","L1","V","S","L2","P","I","L3","L4","F1","F2","R"))
R3$Log2Delta<-ifelse(R3$Log2Delta > 3.5,3.5,R3$Log2Delta)
R3$Log2Delta<-ifelse(R3$Log2Delta < -3.5,-3.5,R3$Log2Delta)
#"C","L1","V","S","L2","P","I","L3","L4","F1","F2","R"

pdf(paste0(path_plot,"FigureS5/Doa10_R3_Saturation_Mutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(R3)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("C" = "C",
                                "L1" = "L",
                                "V" = "V",
                                "S" = "S",
                                "L2" = "L",
                                "P" = "P",
                                "I" = "I",
                                "L3" = "L",
                                "L4" = "L",
                                "F1" = "F",
                                "F2" = "F",
                                "R" = "R"))+
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
      labs(title =  "R3 WT " ),
    
    ggplot(R3)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("C" = "C",
                                "L1" = "L",
                                "V" = "V",
                                "S" = "S",
                                "L2" = "L",
                                "P" = "P",
                                "I" = "I",
                                "L3" = "L",
                                "L4" = "L",
                                "F1" = "F",
                                "F2" = "F",
                                "R" = "R"))+
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
      labs(title =  "R3 Doa10 "),
    
    
    
    
    ggplot(R3)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("C" = "C",
                                "L1" = "L",
                                "V" = "V",
                                "S" = "S",
                                "L2" = "L",
                                "P" = "P",
                                "I" = "I",
                                "L3" = "L",
                                "L4" = "L",
                                "F1" = "F",
                                "F2" = "F",
                                "R" = "R"))+
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R3 Delta ")
    
  ))
dev.off()

# swapping R4 and R5
R4<-subset(data,data$group== "R4")
R4$variable<-factor(R4$variable,
                    levels = c("R1","V",	"T",	"R2","L",	"D",	"A",	"I1",	"Y",	"F",	"R3",	"I2"))

pdf(paste0(path_plot,"FigureS5/Doa10_R4_Saturation_Mutagenesis.pdf"))
#"C","L1","V","S","L2","P","I","L3","L4","F1","F2","R"
print(
  plot_grid(
    
    ggplot(R4)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("R1" = "R",
                                "V" = "V",	
                                "T" = "T",	
                                "R2" = "R",
                                "L" = "L",	
                                "D" = "D",	
                                "A" = "A",	
                                "I1" = "I",	
                                "Y" = "Y",	
                                "F" = "F",	
                                "R3" = "R",	
                                "I2" = "I"))+
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
      labs(title =  "R4 WT " ),
    
    ggplot(R4)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("R1" = "R",
                                "V" = "V",	
                                "T" = "T",	
                                "R2" = "R",
                                "L" = "L",	
                                "D" = "D",	
                                "A" = "A",	
                                "I1" = "I",	
                                "Y" = "Y",	
                                "F" = "F",	
                                "R3" = "R",	
                                "I2" = "I"))+
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
      labs(title =  "R4 Doa10 "),
    
    
    
    
    ggplot(R4)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("R1" = "R",
                                "V" = "V",	
                                "T" = "T",	
                                "R2" = "R",
                                "L" = "L",	
                                "D" = "D",	
                                "A" = "A",	
                                "I1" = "I",	
                                "Y" = "Y",	
                                "F" = "F",	
                                "R3" = "R",	
                                "I2" = "I"))+
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R4 Delta ")
    
  ))

dev.off()

R5<-subset(data,data$group== "R5")
R5$variable<-factor(R5$variable,
                    levels = c("G1",	"L",	"Y",	"M",	"C1",	"S",	"V",	"G2",	"C2",	"W",	"F",	"R"))
#"C","L1","V","S","L2","P","I","L3","L4","F1","F2","R"
pdf(paste0(path_plot,"FigureS5/Doa10_R4_Saturation_Mutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(R5)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("G1" = "G",	
                                "L" = "L",
                                "Y" = "Y",
                                "M" = "M",
                                "C1" = "C",
                                "S" = "S",
                                "V"="V",
                                "G2" = "G",
                                "C2" = "C",
                                "W" = "W",
                                "F" = "F",
                                "R" = "R"))+
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
      labs(title =  "R5 WT " ),
    
    ggplot(R5)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("G1" = "G",	
                                "L" = "L",
                                "Y" = "Y",
                                "M" = "M",
                                "C1" = "C",
                                "S" = "S",
                                "V"="V",
                                "G2" = "G",
                                "C2" = "C",
                                "W" = "W",
                                "F" = "F",
                                "R" = "R"))+
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
      labs(title =  "R5 Doa10 "),
    
    
    
    
    ggplot(R5)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("G1" = "G",	
                                "L" = "L",
                                "Y" = "Y",
                                "M" = "M",
                                "C1" = "C",
                                "S" = "S",
                                "V"="V",
                                "G2" = "G",
                                "C2" = "C",
                                "W" = "W",
                                "F" = "F",
                                "R" = "R"))+
      scale_fill_gradientn("Log2(Delta)",colours=c(
        "#F7931E", 
        "#F4BB75", 
        "#F3CB98", 
        "#F1DBBB", 
        "#F0EBDE", 
        "#C4C0E4", 
        "#9996EA", 
        "#6D6BF0", 
        "#0000FF"
      ), na.value = "grey98",
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "R5 Delta ")
    
  ))

dev.off()


Atg1C_WT<-melt(Atg1C_WT, id = "Saturation")
Atg1C_WT$value<-as.numeric(Atg1C_WT$value)
Atg1C_WT$value<-round(Atg1C_WT$value,4)
Atg1C_WT$value_norm<-Atg1C_WT$value/(subset(Atg1C_WT,Atg1C_WT$Saturation == "N" & Atg1C_WT$variable == "N2")[,c("value")])
Atg1C_WT$value_log2<-log2(Atg1C_WT$value_norm)
Atg1C_WT$group<-"Atg1C"


Atg1C_DT<-melt(Atg1C_DT, id = "Saturation")
Atg1C_DT$value<-as.numeric(Atg1C_DT$value)
Atg1C_DT$value<-round(Atg1C_DT$value,4)
Atg1C_DT$value_norm<-Atg1C_DT$value/(subset(Atg1C_DT,Atg1C_DT$Saturation == "N" & Atg1C_DT$variable == "N2")[,c("value")])
Atg1C_DT$value_log2<-log2(Atg1C_DT$value_norm)
Atg1C_DT$group<-"Atg1C"


Atg1C_WT$identifier<-paste0(Atg1C_WT$Saturation,Atg1C_WT$variable)
Atg1C_DT$identifier<-paste0(Atg1C_DT$Saturation,Atg1C_DT$variable)
Atg1C<-merge(Atg1C_WT,Atg1C_DT, by = "identifier")
Atg1C<-Atg1C[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(Atg1C)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
                "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
Atg1C$Log2Delta<-log2(Atg1C$PSI_Das1/Atg1C$PSI_WT)

Rpa12C_WT<-melt(Rpa12C_WT, id = "Saturation")
Rpa12C_WT$value<-as.numeric(Rpa12C_WT$value)
Rpa12C_WT$value<-round(Rpa12C_WT$value,4)
Rpa12C_WT$value_norm<-Rpa12C_WT$value/(subset(Rpa12C_WT,Rpa12C_WT$Saturation == "N" & Rpa12C_WT$variable == "N2")[,c("value")])
Rpa12C_WT$value_log2<-log2(Rpa12C_WT$value_norm)
Rpa12C_WT$group<-"Rpa12C"


Rpa12C_DT<-melt(Rpa12C_DT, id = "Saturation")
Rpa12C_DT$value<-as.numeric(Rpa12C_DT$value)
Rpa12C_DT$value<-round(Rpa12C_DT$value,4)
Rpa12C_DT$value_norm<-Rpa12C_DT$value/(subset(Rpa12C_DT,Rpa12C_DT$Saturation == "N" & Rpa12C_DT$variable == "N2")[,c("value")])
Rpa12C_DT$value_log2<-log2(Rpa12C_DT$value_norm)
Rpa12C_DT$group<-"Rpa12C"


Rpa12C_WT$identifier<-paste0(Rpa12C_WT$Saturation,Rpa12C_WT$variable)
Rpa12C_DT$identifier<-paste0(Rpa12C_DT$Saturation,Rpa12C_DT$variable)
Rpa12C<-merge(Rpa12C_WT,Rpa12C_DT, by = "identifier")
Rpa12C<-Rpa12C[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(Rpa12C)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
                 "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
Rpa12C$Log2Delta<-log2(Rpa12C$PSI_Das1/Rpa12C$PSI_WT)
######
W1_WT<-melt(W1_WT, id = "Saturation")
W1_WT$value<-as.numeric(W1_WT$value)
W1_WT$value<-round(W1_WT$value,4)
W1_WT$value_norm<-W1_WT$value/(subset(W1_WT,W1_WT$Saturation == "G" & W1_WT$variable == "G4")[,c("value")])
W1_WT$value_log2<-log2(W1_WT$value_norm)
W1_WT$group<-"W1"


W1_DT<-melt(W1_DT, id = "Saturation")
W1_DT$value<-as.numeric(W1_DT$value)
W1_DT$value<-round(W1_DT$value,4)
W1_DT$value_norm<-W1_DT$value/(subset(W1_DT,W1_DT$Saturation == "G" & W1_DT$variable == "G4")[,c("value")])
W1_DT$value_log2<-log2(W1_DT$value_norm)
W1_DT$group<-"W1"


W1_WT$identifier<-paste0(W1_WT$Saturation,W1_WT$variable)
W1_DT$identifier<-paste0(W1_DT$Saturation,W1_DT$variable)
W1<-merge(W1_WT,W1_DT, by = "identifier")
W1<-W1[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(W1)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
W1$Log2Delta<-log2(W1$PSI_Das1/W1$PSI_WT)

W2_WT<-melt(W2_WT, id = "Saturation")
W2_WT$value<-as.numeric(W2_WT$value)
W2_WT$value<-round(W2_WT$value,4)
W2_WT$value_norm<-W2_WT$value/(subset(W2_WT,W2_WT$Saturation == "S" & W2_WT$variable == "S")[,c("value")])
W2_WT$value_log2<-log2(W2_WT$value_norm)
W2_WT$group<-"W2"


W2_DT<-melt(W2_DT, id = "Saturation")
W2_DT$value<-as.numeric(W2_DT$value)
W2_DT$value<-round(W2_DT$value,4)
W2_DT$value_norm<-W2_DT$value/(subset(W2_DT,W2_DT$Saturation == "S" & W2_DT$variable == "S")[,c("value")])
W2_DT$value_log2<-log2(W2_DT$value_norm)
W2_DT$group<-"W2"


W2_WT$identifier<-paste0(W2_WT$Saturation,W2_WT$variable)
W2_DT$identifier<-paste0(W2_DT$Saturation,W2_DT$variable)
W2<-merge(W2_WT,W2_DT, by = "identifier")
W2<-W2[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(W2)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
W2$Log2Delta<-log2(W2$PSI_Das1/W2$PSI_WT)

W3_WT<-melt(W3_WT, id = "Saturation")
W3_WT$value<-as.numeric(W3_WT$value)
W3_WT$value<-round(W3_WT$value,4)
W3_WT$value_norm<-W3_WT$value/(subset(W3_WT,W3_WT$Saturation == "G" & W3_WT$variable == "G")[,c("value")])
W3_WT$value_log2<-log2(W3_WT$value_norm)
W3_WT$group<-"W3"


W3_DT<-melt(W3_DT, id = "Saturation")
W3_DT$value<-as.numeric(W3_DT$value)
W3_DT$value<-round(W3_DT$value,4)
W3_DT$value_norm<-W3_DT$value/(subset(W3_DT,W3_DT$Saturation == "G" & W3_DT$variable == "G")[,c("value")])
W3_DT$value_log2<-log2(W3_DT$value_norm)
W3_DT$group<-"W3"


W3_WT$identifier<-paste0(W3_WT$Saturation,W3_WT$variable)
W3_DT$identifier<-paste0(W3_DT$Saturation,W3_DT$variable)
W3<-merge(W3_WT,W3_DT, by = "identifier")
W3<-W3[,c("Saturation.x" ,"variable.x" ,  "value.x"  ,    "value_norm.x" ,"value_log2.x","value.y"    ,  "value_norm.y", "value_log2.y" ,"group.y")]
names(W3)<-c("Saturation","variable","PSI_WT","value_norm_PSI","value_log2_PSI",
             "PSI_Das1","value_norm_Das1","value_log2_Das1","group"
)
W3$Log2Delta<-log2(W3$PSI_Das1/W3$PSI_WT)


data1<-as.data.frame(rbind(Atg1C,Rpa12C,W1,W2,W3))
data1$Saturation<-factor(data1$Saturation,
                         levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))

#data1$amino_acid<-factor(data1$amino_acid, levels = rev(c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")))

data1$PSI_WT<-ifelse(data1$PSI_WT>0.75,0.75,data1$PSI_WT)
data1$PSI_Das1<-ifelse(data1$PSI_Das1>0.75,0.75,data1$PSI_Das1)
## atg1c
#pdf("Y:/lab data1/susmitha/edwin/for_paper/plot/SaturationMutagenesis/SatMutv2/IN1.pdf")
Atg1C<-subset(data1,data1$group== "Atg1C")
Atg1C$variable<-factor(Atg1C$variable,
                       levels = c("V1"	,"R1"	,"G1",	"A",	"V2",	"G2",	"G3",	"W",	"R2",	"L",	"V3",	"G4"))
pdf(paste0(path_plot,"Figure6/Das1_Atg1C_Saturation_Mutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(Atg1C)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      labs(title =  "Atg1C WT "
      ),
    
    ggplot(Atg1C)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      labs(title =  "Atg1C Das1 "
      ),
    
    
    
    
    ggplot(Atg1C)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "Atg1C Delta "
      )
    
  ))

dev.off()


## rpa12c

Rpa12C<-subset(data1,data1$group== "Rpa12C")
Rpa12C$variable<-factor(Rpa12C$variable,
                        levels = c("V1"	,"R1"	,"G1",	"A",	"V2",	"G2",	"G3",	"W",	"R2",	"L",	"V3",	"G4"))
pdf(paste0(path_plot,"Figure6/Das1_Rpa12C_Saturation_Mutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(Rpa12C)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      labs(title =  "Rpa12C WT "
      ),
    
    ggplot(Rpa12C)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      labs(title =  "Rpa12C Das1 "
      ),
    
    
    
    
    ggplot(Rpa12C)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "Rpa12C Delta "
      )
    
  ))

dev.off()
##

## w1

W1<-subset(data1,data1$group== "W1")
W1$variable<-factor(W1$variable,
                    levels = c("V1"	,"R1"	,"G1",	"A",	"V2",	"G2",	"G3",	"W",	"R2",	"L",	"V3",	"G4"))
pdf(paste0(path_plot,"Figure5/Das1_W1_Saturation_Mutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(W1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      labs(title =  "W1 WT "
      ),
    
    ggplot(W1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      labs(title =  "W1 Das1 "
      ),
    
    
    
    
    ggplot(W1)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("V1" = "V"	,
                                "R1" = "R"	,
                                "G1" = "G",
                                "A" = "A",
                                "V2" = "V",
                                "G2" = "G",
                                "G3" = "G",
                                "W" = "W",
                                "R2" = "R",
                                "L" = "L",
                                "V3" = "V",
                                "G4" = "G"))+
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
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "W1 Delta "
      )
    
  ))

dev.off()

W2<-subset(data1,data1$group== "W2")
W2$variable<-factor(W2$variable,
                    levels = c("A1"	,
                               "R1"	,
                               "W1",
                               "R2",	
                               "V" ,	
                               "G1",
                               "L" ,
                               "W2",	
                               "R3",
                               "G2",
                               "A2",
                               "S"))
pdf(paste0(path_plot,"Figure5/Das1_W2_Saturation_Mutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(W2)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("A1" = "A"	,
                                "R1" = "R"	,
                                "W1" = "W",
                                "R2" = "R",	
                                "V" = "V",	
                                "G1" = "G",
                                "L" = "L",
                                "W2" = "W",	
                                "R3" = "R",
                                "G2" = "G",
                                "A2" = "A",
                                "S" = "S"))+
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
      labs(title =  "W2 WT "
      ),
    
    ggplot(W2)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("A1" = "A"	,
                                "R1" = "R"	,
                                "W1" = "W",
                                "R2" = "R",	
                                "V" = "V",	
                                "G1" = "G",
                                "L" = "L",
                                "W2" = "W",	
                                "R3" = "R",
                                "G2" = "G",
                                "A2" = "A",
                                "S" = "S"))+
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
      labs(title =  "W2 Das1 "
      ),
    
    
    
    
    ggplot(W2)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("A1" = "A"	,
                                "R1" = "R"	,
                                "W1" = "W",
                                "R2" = "R",	
                                "V" = "V",	
                                "G1" = "G",
                                "L" = "L",
                                "W2" = "W",	
                                "R3" = "R",
                                "G2" = "G",
                                "A2" = "A",
                                "S" = "S"))+
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
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "W2 Delta "
      )
    
  ))

dev.off()

W3<-subset(data1,data1$group== "W3")
W3$variable<-factor(W3$variable,
                    levels = c("P",  "S1" ,     "L"   ,   "A1" , "A2",  "S2" ,     "F" , "W" , "R"    ,  "V1" ,"V2" ,     "G"))
pdf(paste0(path_plot,"Figure5/Das1_W3_Saturation_Mutagenesis.pdf"))
print(
  plot_grid(
    
    ggplot(W3)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_WT 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("P" = "P",
                                "S1" = "S",
                                "L" = "L",
                                "A1" = "A", 
                                "A2" = "A",
                                "S2" = "S" ,
                                "F" = "F" ,
                                "W" = "W",
                                "R" = "R" ,
                                "V1" = "V" ,
                                "V2" = "V" ,
                                "G" = "G"))+
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
      labs(title =  "W3 WT "
      ),
    
    ggplot(W3)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = PSI_Das1 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("P" = "P",
                                "S1" = "S",
                                "L" = "L",
                                "A1" = "A", 
                                "A2" = "A",
                                "S2" = "S" ,
                                "F" = "F" ,
                                "W" = "W",
                                "R" = "R" ,
                                "V1" = "V" ,
                                "V2" = "V" ,
                                "G" = "G"))+
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
      labs(title =  "W3 Das1 "
      ),
    
    
    
    
    ggplot(W3)+
      geom_tile(aes(
        x = variable  ,
        y=Saturation,
        fill = Log2Delta 
      ))+
      theme_bw()+
      theme(text   = element_text(size = 8)) +
      scale_x_discrete(labels=c("P" = "P",
                                "S1" = "S",
                                "L" = "L",
                                "A1" = "A", 
                                "A2" = "A",
                                "S2" = "S" ,
                                "F" = "F" ,
                                "W" = "W",
                                "R" = "R" ,
                                "V1" = "V" ,
                                "V2" = "V" ,
                                "G" = "G"))+
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
      values = scales::rescale(round(seq(-3.5,3.5,7/10),2)),
      limits = c(-3.5,3.5))+
      
      
      ylab("Mutations")+
      xlab("Saturated Sequence ")+
      labs(title =  "W3 Delta "
      )
    
  ))

dev.off()



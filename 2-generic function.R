##### for project : Analysis of degron at C terminal
##### sub text:generic functiona
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

histogram_for_raw_count<-function(dataset,raw,class){
  
  Quality_Control_histogram<-dataset %>% ggplot(aes(x = raw)) +
    geom_histogram(aes(color = class, fill = class), 
                   position = "identity")+
    labs(title = ("Histogram depicting frequency of samples for a given raw count"), x = "Raw count", y = "Frequency")+
    theme(text = element_text(size = 8,angle = 0)) +
    theme(axis.text.y = element_text(angle = 0))+
    theme(axis.text.x = element_text(angle = 90))
  #Quality_Control_histogram<-Quality_Control_histogram + colScale 
  return(Quality_Control_histogram)
  
}


bin_wise_cluster_translation<-function(dataset,nraw,rev,translation){
  bin_cluster_translation<-dataset %>%
    ggplot(aes(y = nraw, x = rev, color = translation, group = translation)) + 
    geom_line()+
    geom_point()+
    theme_bw()+
    theme(axis.title = element_text(face="bold"))+
    labs(title = "Cluster based analysis for translations", x = "Bins", y = "Normalized Count")
  return(bin_cluster_translation)
}


graph_for_distribution_of_CS_property<-function(dataset,xaxis,yaxis,titlePlot,minimum_y,maximum_y,xaxixlabel,yaxislabel){
  (
    ggplot(dataset, aes(factor(xaxis), yaxis)) + 
      geom_boxplot() +
      theme_bw()+
      theme(axis.text.x = element_text(angle = 90)) +
      #scale_x_discrete(labels=c("1"="-12","2"="-11", "3"="-10","4"="-9","5"="-8", "6"="-7","7"="-6","8"="-5", "9"="-4","10"="-3","11"="-2", "12"="-1"))+
      labs(title = titlePlot, x = xaxixlabel, y = yaxislabel)+
      ylim(minimum_y,maximum_y)
  )
}

CS_per_property_position_wise<-function(dataset){
  
  dataset<-melt(dataset)
  property_name<-unique(str_split_fixed(dataset$variable,"_",3)[,2])
  dataset$variable<-str_split_fixed(dataset$variable,"_",2)[,1]
  dataset$variable<-factor(dataset$variable, levels = as.character(-12:-1))
  plot_propertyCS_per_position<-graph_for_distribution_of_CS_property(dataset,
                                                                      dataset$variable,
                                                                      dataset$value,
                                                                      paste0("Contribution Score for ",property_name," per position"),
                                                                      min(dataset$value),
                                                                      max(dataset$value),
                                                                      "Position","Contribution Score"
                                                                      )
  
  return(plot_propertyCS_per_position)  
}
create_one_hot_db<-function(amino_acid,threshold){
  converted_string<-as.data.frame(lapply(amino_acid,function(x) substr(as.character(x),1,threshold)))
  a<-as.data.frame((as.data.frame(lapply(converted_string[1,],function(x) strsplit(x, "")))))
  a<-as.data.frame(t(as.data.frame(lapply(converted_string[1,],function(x) strsplit(x, "")))))
  row.names(a)<-amino_acid
  df <- dummy.data.frame(a, sep="")
  return(df)
}

set.seed(123)
##### for project : Analysis of degron at C terminal
##### sub text:SHAP dataset creation for deep learning model
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB



setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/input_data/")
path_plot<-paste0(getwd(),"/plot/")
path_output<-paste0(getwd(),"/output/")
path_model<-paste0(getwd(),"/model")

training_data<-read.csv(paste0(path_output,"/training_data_after_normalization_for_model.csv"))
testing_data<-read.csv(paste0(path_output,"/testing_data_after_normalization_for_model.csv"))
row.names(training_data)<-training_data$X
training_data$X<-NULL
row.names(testing_data)<-testing_data$X
testing_data$X<-NULL

all<-as.data.frame(rbind(training_data,testing_data))
allData<-all[,1:318]
all_label<-all[,319]
allData<-as.matrix(allData)
allLabel<-to_categorical(all_label)

model_without_iMet<-load_model_tf(paste0(path_model,"/model_12_x"))
predictor.without.iMet<-Predictor$new(
  model = model_without_iMet,
  data = as.data.frame(allData),
  y = allLabel
)



All<-model_without_iMet %>% predict((allData))
testing_class_iMet<-as.data.frame(All)
row.names(testing_class_iMet)<-row.names(allData)
testing_class_iMet$max<-colnames(testing_class_iMet)[apply(testing_class_iMet,1,which.max)]
testing_class_iMet$class<-ifelse(testing_class_iMet$max == "V1",0,1)
testing_class_iMet$actual_label<-all_label
testing_class_iMet<-testing_class_iMet[,4:5]
names(testing_class_iMet)<-c("Predicted_Class","Actual_Class")
correctly_predicted_testing<-subset(testing_class_iMet,testing_class_iMet$Predicted_Class == testing_class_iMet$Actual_Class)
correctly_predicted_testing_instable<-subset(correctly_predicted_testing,correctly_predicted_testing$Predicted_Class == 1)
correctly_predicted_testing_stable<-subset(correctly_predicted_testing,correctly_predicted_testing$Predicted_Class == 0)

# for correctly predicted instable

testing<-subset(allData, row.names(allData) %in% row.names(correctly_predicted_testing_instable))
testing<-as.data.frame(testing)
shap_sim_phi_test<-as.data.frame(matrix(nrow=636))
shap_sim_fea_val_test<-as.data.frame(matrix(nrow=636))
for(i in 1:4110){
  shap.sequence.test<-Shapley$new(predictor.without.iMet,x.interest = as.data.frame(testing[i,1:318]))
  name<-row.names(testing)[(i)]
  shap_sim_phi_test[,(i)]<-shap.sequence.test$results[,3]
  shap_sim_fea_val_test[,(i)]<-shap.sequence.test$results[,5]
  names(shap_sim_phi_test)[(i)]<-name
  names(shap_sim_fea_val_test)[(i)]<-name
  print((i))
}
path_dataset
shap_sim_fea_val_test_class_1<-shap_sim_fea_val_test[1:318,]
row.names(shap_sim_fea_val_test_class_1)<-shap.sequence.test$results$feature[1:318]
shap_sim_fea_val_test_class_2<-shap_sim_fea_val_test[319:636,]
row.names(shap_sim_fea_val_test_class_2)<-shap.sequence.test$results$feature[319:636]
shap_sim_phi_test_class_1<-shap_sim_phi_test[1:318,]
row.names(shap_sim_phi_test_class_1)<-shap.sequence.test$results$feature[1:318]
shap_sim_phi_test_class_2<-shap_sim_phi_test[319:636,]
row.names(shap_sim_phi_test_class_2)<-shap.sequence.test$results$feature[319:636]

write.csv(shap_sim_fea_val_test_class_1,paste0(path_output,"/shap_instable/shap_sim_fea_val_test_class_1_stable_12x_for_instable_all.csv"))
write.csv(shap_sim_fea_val_test_class_2,paste0(path_output,"/shap_instable/shap_sim_fea_val_test_class_2_stable_12x_for_instable_all.csv"))
write.csv(shap_sim_phi_test_class_1,paste0(path_output,"/shap_instable/shap_sim_phi_test_class_1_stable_12x_for_instable_all.csv"))
write.csv(shap_sim_phi_test_class_2,paste0(path_output,"/shap_instable/shap_sim_phi_test_class_1_stable_12x_for_instable_all.csv"))



### for correctly predicted stable

testing<-subset(allData, row.names(allData) %in% row.names(correctly_predicted_testing_stable))
testing<-as.data.frame(testing)
shap_sim_phi_test<-as.data.frame(matrix(nrow=636))
shap_sim_fea_val_test<-as.data.frame(matrix(nrow=636))
for(i in 1:4485){
  shap.sequence.test<-Shapley$new(predictor.without.iMet,x.interest = as.data.frame(testing[i,1:318]))
  name<-row.names(testing)[(i)]
  shap_sim_phi_test[,(i)]<-shap.sequence.test$results[,3]
  shap_sim_fea_val_test[,(i)]<-shap.sequence.test$results[,5]
  names(shap_sim_phi_test)[(i)]<-name
  names(shap_sim_fea_val_test)[(i)]<-name
  print((i))
}
path_dataset
shap_sim_fea_val_test_class_1<-shap_sim_fea_val_test[1:318,]
row.names(shap_sim_fea_val_test_class_1)<-shap.sequence.test$results$feature[1:318]
shap_sim_fea_val_test_class_2<-shap_sim_fea_val_test[319:636,]
row.names(shap_sim_fea_val_test_class_2)<-shap.sequence.test$results$feature[319:636]
shap_sim_phi_test_class_1<-shap_sim_phi_test[1:318,]
row.names(shap_sim_phi_test_class_1)<-shap.sequence.test$results$feature[1:318]
shap_sim_phi_test_class_2<-shap_sim_phi_test[319:636,]
row.names(shap_sim_phi_test_class_2)<-shap.sequence.test$results$feature[319:636]

write.csv(shap_sim_fea_val_test_class_1,paste0(path_output,"/shap_stable/shap_sim_fea_val_test_class_1_stable_12x_for_stable_all.csv"))
write.csv(shap_sim_fea_val_test_class_2,paste0(path_output,"/shap_stable/shap_sim_fea_val_test_class_2_stable_12x_for_stable_all.csv"))
write.csv(shap_sim_phi_test_class_1,paste0(path_output,"/shap_stable/shap_sim_phi_test_class_1_stable_12x_for_stable_all.csv"))
write.csv(shap_sim_phi_test_class_2,paste0(path_output,"/shap_stable/shap_sim_phi_test_class_1_stable_12x_for_stable_all.csv"))



##### for project : Analysis of degron at C terminal
##### sub text:merging deep learning model
##### Code written by: Susmitha Shankar
##### Lab : Khmelinskii Lab, IMB

setwd("Y:/lab data/susmitha/edwin/DAS1_SCF_Paper/")
path_input<-paste0(getwd(),"/output/")
set.seed(123)
training_data<-read.csv(paste0(path_input,"training_data_after_normalization_for_model.csv"))
testing_data<-read.csv(paste0(path_input,"testing_data_after_normalization_for_model.csv"))

row.names(training_data)<-training_data$X
row.names(testing_data)<-testing_data$X

training_data$X<-NULL
testing_data$X<-NULL

instable_test<-subset(testing_data,testing_data$Label == 1)
stable_test<-subset(testing_data,testing_data$Label == 0)

test1<-as.data.frame(rbind(instable_test[1:300,], stable_test[1:300,]))
test2<-as.data.frame(rbind(instable_test[301:600,], stable_test[301:600,]))
test3<-as.data.frame(rbind(instable_test[601:976,], stable_test[601:1000,]))

test1<-test1[sample(1:nrow(test1), nrow(test1), replace=FALSE),]
test2<-test2[sample(1:nrow(test2), nrow(test2), replace=FALSE),]
test3<-test3[sample(1:nrow(test3), nrow(test3), replace=FALSE),]

write.csv(test1, paste0(path_input, "test1_post_correction.csv"))
write.csv(test2, paste0(path_input, "test2_post_correction.csv"))
write.csv(test3, paste0(path_input, "test3_post_correction.csv"))

all<-as.data.frame(rbind(training_data,testing_data))
allData<-all[,1:318]
all_label<-all[,319]
allData<-as.matrix(allData)
allLabel<-to_categorical(all_label)

trainingData<-training_data[,1:318]
training_label<-training_data[,319]
trainingData<-as.matrix(trainingData)
trainLabel<-to_categorical(training_label)
 50
testingData<-testing_data[,1:318]
testing_label<-testing_data[,319]
testingData<-as.matrix(testingData)
testLabel<-to_categorical(testing_label)

testingData1<-test1[,1:318]
testing_label1<-test1[,319]
testingData1<-as.matrix(testingData1)
testLabel1<-to_categorical(testing_label1)

testingData2<-test2[,1:318]
testing_label2<-test2[,319]
testingData2<-as.matrix(testingData2)
testLabel2<-to_categorical(testing_label2)

testingData3<-test3[,1:318]
testing_label3<-test3[,319]
testingData3<-as.matrix(testingData3)
testLabel3<-to_categorical(testing_label3)

FCNN_model <-keras_model_sequential() 
FCNN_model %>%
  layer_dense(units = 50, 
              activation = "relu", 
              input_shape = c(318)) %>%
  layer_dropout(rate = 0.3)%>%
  layer_dense(units = 12, 
              activation = "relu",
              regularizer_l1_l2(l1 = 0.001, 
                                l2 = 0.001)) %>%
  layer_dropout(rate = 0.3)%>%
  layer_dense(units = 5, 
              activation = "relu", 
              regularizer_l1_l2(l1 = 0.001, 
                                l2 = 0.001)) %>%
  layer_dropout(rate = 0.3)%>%
  layer_dense(units = 2, 
              activation = "softmax")
summary(FCNN_model)
FCNN_model %>% compile(
  loss = 'binary_crossentropy',
  optimizer = optimizer_adam() ,
  metrics = c('accuracy',
              tf$keras$metrics$AUC(),
              tf$keras$metrics$Precision(),
              tf$keras$metrics$Recall(),
              tf$keras$metrics$MeanSquaredError())
)

FCNN_model_history <- FCNN_model %>% fit(
  trainingData, trainLabel,
  batch_size = 500,
  epochs = 30,
  validation_split = 0.2)

testing_1<-FCNN_model %>% 
  evaluate(testingData1,
           testLabel1)
testing_2<-FCNN_model %>% 
  evaluate(testingData2,
           testLabel2)
testing_3<-FCNN_model %>% 
  evaluate(testingData3,
           testLabel3)
testing_<-FCNN_model %>% 
  evaluate(testingData,
           testLabel)

dir.create(paste0(getwd(),"/model"), 
           showWarnings = TRUE, 
           recursive = FALSE, 
           mode = "0777"
)

path_model<-paste0(getwd(),"/model")
FCNN_model %>% save_model_tf(paste0(path_model,"/model_12_x"))


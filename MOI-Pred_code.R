#Import libraries.
library(caret) #v. 6.0.84
library(MLmetrics) #v. 1.1.1
library(data.table) #v. 1.14.0
library(tidyverse) #v. 1.2.1
library(cvAUC) #v. 1.1.0
library(pROC) #v. 1.14.0

#Import datasets.
path = "your_path_to_working_directory"
Discovery = fread(file = file.path(path, "file_with_annotated_ExoVar_and_gnomAD_variants.txt"), data.table = F)
Validation = fread(file = file.path(path, "file_with_annotated_GEM_and_expert_reviewed_variants.txt"), data.table = F)

#Files should have 23 variant-level features, 7 gene-level features and mode-of-inheritance label ("Pheno").
#Variant-level features: SIFT_score, MutationTaster_score, MutationAssessor_score, PROVEAN_score, M-CAP_score, MutPred_score, GenoCanyon_score, integrated_fitCons_score, GERP++_RS, phyloP100way_vertebrate, phyloP20way_mammalian, phastCons100way_vertebrate, phastCons20way_mammalian, SiPhy_29way_logOdds, Eigen, FATHMM_noncoding, GWAVA_tss_score, Kaviar_AF, Polyphen2_HDIV_score, LRT_score, FATHMM_score, VEST3_score, CADD_phred
#Gene-level features: Shet, oe, episcore, AD_rank, HI, StringDB_AD, StringDB_AR

metrics_recessive = list()
metrics_dominant = list()
metrics_benign = list()

for( i in 1:100 ) {

#Define variants to use for training, internal testing and external validation.
I = sample(which(Discovery$Pheno == "Dominant"), size = round(length(which(Discovery$Pheno == "Dominant")) * 0.90)) 
I = c(I, sample(which(Discovery$Pheno == "Recessive"), size = length(I)))
I = c(I, sample(which(Discovery$Pheno == "Benign"), size = (length(I) / 2)))
Discovery_train = Discovery[I,]
Discovery_test = Discovery[-I,]
I = which(Discovery_test$Pheno == "Dominant")
I = c(I, sample(which(Discovery_test$Pheno == "Recessive"), size = length(I)))
I = c(I, sample(which(Discovery_test$Pheno == "Benign"), size = (length(I) / 2)))
Discovery_test = Discovery_test[I,]
J = c(length(which(Discovery_test$Pheno[which(Discovery_test$gene %in% Discovery_train$gene)] == "Benign")), 
length(which(Discovery_test$Pheno[which(Discovery_test$gene %in% Discovery_train$gene)] == "Dominant")), 
length(which(Discovery_test$Pheno[which(Discovery_test$gene %in% Discovery_train$gene)] == "Recessive")))
I = which(Discovery_test$gene %in% Discovery_train$gene)
Discovery_train = rbind(Discovery_train, Discovery_test[I,]) 
Discovery_test = Discovery_test[-I,]
A = as.data.frame.table(table(Discovery_train$gene))
A = A[sort(A$Freq, decreasing = FALSE, index.return=TRUE)$ix,]
A[,1] = as.character(A[,1])
A2 = A[which(A$Freq == 2),]
P = Discovery_train[which(Discovery_train$gene %in% A[which(A[,2] == 1),1]),]
P2 = Discovery_train[which(Discovery_train$gene %in% A[which(A[,2] == 2),1]),]
I = c()
if( length(which(P$Pheno == "Benign")) >= J[1] ) { I = rbind(I, P[sample(which(P$Pheno == "Benign"), size = J[1]),])
J[1] = 0 } else { I = rbind(I, P[which(P$Pheno == "Benign"),])
J[1] = J[1] - length(which(P$Pheno == "Benign")) }
if( length(which(P$Pheno == "Dominant")) >= J[2] ) { I = rbind(I, P[sample(which(P$Pheno == "Dominant"), size = J[2]),])
J[2] = 0 } else { I = rbind(I, P[which(P$Pheno == "Dominant"),])
J[2] = J[2] - length(which(P$Pheno == "Dominant")) }
if( length(which(P$Pheno == "Recessive")) >= J[3] ) { I = rbind(I, P[sample(which(P$Pheno == "Recessive"), size = J[3]),])
J[3] = 0 } else { I = rbind(I, P[which(P$Pheno == "Recessive"),])
J[3] = J[3] - length(which(P$Pheno == "Recessive")) }
P = c()
for( g in unique(P2$gene) ) { if( length(names(table(P2$Pheno[which(P2$gene == g)]))) == 1 ) { P = c(P, g) } }
P = P2[which(P2$gene %in% P),]
if( J[1] != 0 ) { A = A2[which(A2[,1] %in% P$gene[which(P$Pheno == "Benign")]),1]
I = rbind(I, P[which(P$gene %in% A[1:(J[1] / 2)]),]) }
if( J[2] != 0 ) { A = A2[which(A2[,1] %in% P$gene[which(P$Pheno == "Dominant")]),1]
I = rbind(I, P[which(P$gene %in% A[1:(J[2] / 2)]),]) }
if( J[3] != 0 ) { A = A2[which(A2[,1] %in% P$gene[which(P$Pheno == "Recessive")]),1]
I = rbind(I, P[which(P$gene %in% A[1:(J[3] / 2)]),]) }
Discovery_test = rbind(Discovery_test, I)
Discovery_train = Discovery_train[-which(Discovery_train$coord %in% I$coord),]
I = which(Validation$Pheno == "Dominant")
I = c(I, sample(which(Validation$Pheno == "Recessive"), size = length(I)))
I = c(I, sample(which(Validation$Pheno == "Benign"), size = (length(I) / 2)))
Validation_test = Validation[I,]

#Scaling.
fit_scale = preProcess(Discovery_train[,-which(colnames(Discovery_train) %in% c("coord", "gene", "Pheno"))], method = c("center", "scale"))
Discovery_train[,-which(colnames(Discovery_train) %in% c("coord", "gene", "Pheno"))] = 
predict(fit_scale, Discovery_train[,-which(colnames(Discovery_train) %in% c("coord", "gene", "Pheno"))])
Discovery_test[,-which(colnames(Discovery_test) %in% c("coord", "gene", "Pheno"))] = 
predict(fit_scale, Discovery_test[,-which(colnames(Discovery_test) %in% c("coord", "gene", "Pheno"))])
Validation_test[,-which(colnames(Validation_test) %in% c("coord", "gene", "Pheno"))] = 
predict(fit_scale, Validation_test[,-which(colnames(Validation_test) %in% c("coord", "gene", "Pheno"))])

#Feature selection.
BB = rfe(x = Discovery_train[,-which(colnames(Discovery_train) %in% c("coord", "gene", "Pheno"))], y = as.factor(Discovery_train$Pheno), 
sizes = c(1:30), rfeControl = rfeControl(functions=rfFuncs, method="cv", number=10))
W = which(colnames(Discovery_train) %in% BB$optVariables)
Discovery_train = Discovery_train[,c(W, (dim(Discovery_train)[2] - 2):dim(Discovery_train)[2])]
Discovery_test = Discovery_test[,which(colnames(Discovery_test) %in% colnames(Discovery_train))]
Validation_test = Validation_test[,which(colnames(Validation_test) %in% colnames(Discovery_train))]


#Train random forest model.
fitControl_10CV = trainControl(method = "cv", number = 10, savePredictions="final", classProbs=T, summaryFunction=defaultSummary)
fit_3way = train(Pheno ~ ., data = Discovery_train[,-which(colnames(Discovery_train) %in% c("coord", "gene"))], method = "rf", trControl=fitControl_10CV)

#Predict MOI for variants in internal test set.
YtestPredRaw = predict(fit_3way, Discovery_test[,-c((dim(Discovery_test)[2] - 2):dim(Discovery_test)[2])], type = "raw")
YtestPredProb = predict(fit_3way, Discovery_test[,-c((dim(Discovery_test)[2] - 2):dim(Discovery_test)[2])], type = "prob")
YtestTrue = Discovery_test$Pheno

#Predict MOI for variants in external validation set.
YvalPredRaw = predict(fit_3way, Validation_test[,-c((dim(Validation_test)[2] - 2):dim(Validation_test)[2])], type = "raw")
YvalPredProb = predict(fit_3way, Validation_test[,-c((dim(Validation_test)[2] - 2):dim(Validation_test)[2])], type = "prob") 
YvalTrue = Validation_test$Pheno

#Compute performance metrics.
Y_recessive_raw = as.character(fit_3way$pred$pred)
Y_recessive_raw[which(Y_recessive_raw != "Recessive")] = "0"
Y_recessive_raw[which(Y_recessive_raw == "Recessive")] = "1"
Y_recessive_true = as.character(fit_3way$pred$obs)
Y_recessive_true[which(Y_recessive_true != "Recessive")] = "0"
Y_recessive_true[which(Y_recessive_true == "Recessive")] = "1"
ConfM_recessive=confusionMatrix(data=as.factor(Y_recessive_raw),reference=as.factor(Y_recessive_true),positive="1")
Y_recessive_prob = as.character(fit_3way$pred$Recessive)
AUC = auc(as.numeric(Y_recessive_true), as.numeric(Y_recessive_prob))
metrics = data.frame(F1 = ConfM_recessive$byClass[[7]], AUC = AUC, Accuracy = ConfM_recessive$overall[[1]], 
Sensitivity = ConfM_recessive$byClass[[1]], Specificity = ConfM_recessive$byClass[[2]], NPV = ConfM_recessive$byClass[[4]], 
PPV = ConfM_recessive$byClass[[3]])
Y_recessive_raw = as.character(YtestPredRaw)
Y_recessive_raw[which(Y_recessive_raw != "Recessive")] = "0"
Y_recessive_raw[which(Y_recessive_raw == "Recessive")] = "1"
Y_recessive_true = as.character(YtestTrue)
Y_recessive_true[which(Y_recessive_true != "Recessive")] = "0"
Y_recessive_true[which(Y_recessive_true == "Recessive")] = "1"
ConfM_recessive=confusionMatrix(data=as.factor(Y_recessive_raw),reference=as.factor(Y_recessive_true),positive="1")
Y_recessive_prob = YtestPredProb$Recessive
AUC = auc(as.numeric(Y_recessive_true), Y_recessive_prob)
I = data.frame(F1 = ConfM_recessive$byClass[[7]], AUC = AUC, Accuracy = ConfM_recessive$overall[[1]], 
Sensitivity = ConfM_recessive$byClass[[1]], Specificity = ConfM_recessive$byClass[[2]], NPV = ConfM_recessive$byClass[[4]], 
PPV = ConfM_recessive$byClass[[3]])
metrics = rbind(metrics, I)
Y_recessive_raw = as.character(YvalPredRaw)
Y_recessive_raw[which(Y_recessive_raw != "Recessive")] = "0"
Y_recessive_raw[which(Y_recessive_raw == "Recessive")] = "1"
Y_recessive_true = as.character(YvalTrue)
Y_recessive_true[which(Y_recessive_true != "Recessive")] = "0"
Y_recessive_true[which(Y_recessive_true == "Recessive")] = "1"
ConfM_recessive=confusionMatrix(data=as.factor(Y_recessive_raw),reference=as.factor(Y_recessive_true),positive="1")
Y_recessive_prob = YvalPredProb$Recessive
AUC = auc(as.numeric(Y_recessive_true), Y_recessive_prob)
I = data.frame(F1 = ConfM_recessive$byClass[[7]], AUC = AUC, Accuracy = ConfM_recessive$overall[[1]], 
Sensitivity = ConfM_recessive$byClass[[1]], Specificity = ConfM_recessive$byClass[[2]], NPV = ConfM_recessive$byClass[[4]], 
PPV = ConfM_recessive$byClass[[3]])
metrics = rbind(metrics, I)
metrics_recessive = metrics

Y_dominant_raw = as.character(fit_3way$pred$pred)
Y_dominant_raw[which(Y_dominant_raw != "Dominant")] = "0"
Y_dominant_raw[which(Y_dominant_raw == "Dominant")] = "1"
Y_dominant_true = as.character(fit_3way$pred$obs)
Y_dominant_true[which(Y_dominant_true != "Dominant")] = "0"
Y_dominant_true[which(Y_dominant_true == "Dominant")] = "1"
ConfM_dominant=confusionMatrix(data=as.factor(Y_dominant_raw),reference=as.factor(Y_dominant_true),positive="1")
Y_dominant_prob = as.character(fit_3way$pred$Dominant)
AUC = auc(as.numeric(Y_dominant_true), as.numeric(Y_dominant_prob))
metrics = data.frame(F1 = ConfM_dominant$byClass[[7]], AUC = AUC, Accuracy = ConfM_dominant$overall[[1]], 
Sensitivity = ConfM_dominant$byClass[[1]], Specificity = ConfM_dominant$byClass[[2]], NPV = ConfM_dominant$byClass[[4]], 
PPV = ConfM_dominant$byClass[[3]])
Y_dominant_raw = as.character(YtestPredRaw)
Y_dominant_raw[which(Y_dominant_raw != "Dominant")] = "0"
Y_dominant_raw[which(Y_dominant_raw == "Dominant")] = "1"
Y_dominant_true = as.character(YtestTrue)
Y_dominant_true[which(Y_dominant_true != "Dominant")] = "0"
Y_dominant_true[which(Y_dominant_true == "Dominant")] = "1"
ConfM_dominant=confusionMatrix(data=as.factor(Y_dominant_raw),reference=as.factor(Y_dominant_true),positive="1")
Y_dominant_prob = YtestPredProb$Dominant
AUC = auc(as.numeric(Y_dominant_true), Y_dominant_prob)
I = data.frame(F1 = ConfM_dominant$byClass[[7]], AUC = AUC, Accuracy = ConfM_dominant$overall[[1]], 
Sensitivity = ConfM_dominant$byClass[[1]], Specificity = ConfM_dominant$byClass[[2]], NPV = ConfM_dominant$byClass[[4]], 
PPV = ConfM_dominant$byClass[[3]])
metrics = rbind(metrics, I)
Y_dominant_raw = as.character(YvalPredRaw)
Y_dominant_raw[which(Y_dominant_raw != "Dominant")] = "0"
Y_dominant_raw[which(Y_dominant_raw == "Dominant")] = "1"
Y_dominant_true = as.character(YvalTrue)
Y_dominant_true[which(Y_dominant_true != "Dominant")] = "0"
Y_dominant_true[which(Y_dominant_true == "Dominant")] = "1"
ConfM_dominant=confusionMatrix(data=as.factor(Y_dominant_raw),reference=as.factor(Y_dominant_true),positive="1")
Y_dominant_prob = YvalPredProb$Dominant
AUC = auc(as.numeric(Y_dominant_true), Y_dominant_prob)
I = data.frame(F1 = ConfM_dominant$byClass[[7]], AUC = AUC, Accuracy = ConfM_dominant$overall[[1]], 
Sensitivity = ConfM_dominant$byClass[[1]], Specificity = ConfM_dominant$byClass[[2]], NPV = ConfM_dominant$byClass[[4]], 
PPV = ConfM_dominant$byClass[[3]])
metrics = rbind(metrics, I)
metrics_dominant = metrics

Y_benign_raw = as.character(fit_3way$pred$pred)
Y_benign_raw[which(Y_benign_raw != "Benign")] = "0"
Y_benign_raw[which(Y_benign_raw == "Benign")] = "1"
Y_benign_true = as.character(fit_3way$pred$obs)
Y_benign_true[which(Y_benign_true != "Benign")] = "0"
Y_benign_true[which(Y_benign_true == "Benign")] = "1"
ConfM_benign=confusionMatrix(data = as.factor(Y_benign_raw), reference= as.factor(Y_benign_true), positive = "1")
Y_benign_prob = as.character(fit_3way$pred$Benign)
AUC = auc(as.numeric(Y_benign_true), as.numeric(Y_benign_prob))
metrics = data.frame(F1 = ConfM_benign$byClass[[7]], AUC = AUC, Accuracy = ConfM_benign$overall[[1]], 
Sensitivity = ConfM_benign$byClass[[1]], Specificity = ConfM_benign$byClass[[2]], NPV = ConfM_benign$byClass[[4]], 
PPV = ConfM_benign$byClass[[3]])
Y_benign_raw = as.character(YtestPredRaw)
Y_benign_raw[which(Y_benign_raw != "Benign")] = "0"
Y_benign_raw[which(Y_benign_raw == "Benign")] = "1"
Y_benign_true = as.character(YtestTrue)
Y_benign_true[which(Y_benign_true != "Benign")] = "0"
Y_benign_true[which(Y_benign_true == "Benign")] = "1"
ConfM_benign=confusionMatrix(data = as.factor(Y_benign_raw), reference= as.factor(Y_benign_true), positive = "1")
Y_benign_prob = YtestPredProb$Benign
AUC = auc(as.numeric(Y_benign_true), Y_benign_prob)
I = data.frame(F1 = ConfM_benign$byClass[[7]], AUC = AUC, Accuracy = ConfM_benign$overall[[1]], 
Sensitivity = ConfM_benign$byClass[[1]], Specificity = ConfM_benign$byClass[[2]], NPV = ConfM_benign$byClass[[4]], 
PPV = ConfM_benign$byClass[[3]])
metrics = rbind(metrics, I)
Y_benign_raw = as.character(YvalPredRaw)
Y_benign_raw[which(Y_benign_raw != "Benign")] = "0"
Y_benign_raw[which(Y_benign_raw == "Benign")] = "1"
Y_benign_true = as.character(YvalTrue)
Y_benign_true[which(Y_benign_true != "Benign")] = "0"
Y_benign_true[which(Y_benign_true == "Benign")] = "1"
ConfM_benign=confusionMatrix(data = as.factor(Y_benign_raw), reference= as.factor(Y_benign_true), positive = "1")
Y_benign_prob = YvalPredProb$Benign
AUC = auc(as.numeric(Y_benign_true), Y_benign_prob)
I = data.frame(F1 = ConfM_benign$byClass[[7]], AUC = AUC, Accuracy = ConfM_benign$overall[[1]], 
Sensitivity = ConfM_benign$byClass[[1]], Specificity = ConfM_benign$byClass[[2]], NPV = ConfM_benign$byClass[[4]], 
PPV = ConfM_benign$byClass[[3]])
metrics = rbind(metrics, I)
metrics_benign = metrics

#Save models and vectors.
write.table(x = fit_3way$pred, file = file.path(path, "/Train/",i,"Ytrain.txt"), quote = F, sep="\t", row.names=T)
write.table(x = as.data.frame(cbind(YtestPredRaw, YtestTrue)), file = file.path(path, "/Test/",i,"YtestRaw.txt"), quote = F, sep = "\t", row.names = T)
write.table(x = as.data.frame(cbind(YtestPredProb, YtestTrue)), file = file.path(path, "/Test/",i,"YtestProb.txt"), quote = F, sep = "\t", row.names = T)
write.table(x = as.data.frame(cbind(YvalPredRaw, YvalTrue)), file = file.path(path, "/Validation/",i,"YvalRaw.txt"), quote = F, sep = "\t", row.names = T)
write.table(x = as.data.frame(cbind(YvalPredProb, YvalTrue)), file = file.path(path, "/Validation/",i,"YvalProb.txt"), quote = F, sep = "\t", row.names = T)
saveRDS(object = Discovery_train$coord, file = file.path(path, "/Train/",i,"IDs_train.RDS"))
saveRDS(object = Discovery_test$coord, file = file.path(path, "/Test/",i,"IDs_test.RDS"))
saveRDS(object = Validation_test$coord, file = file.path(path, "/Validation/",i,"IDs_val.RDS"))

write.table(metrics_benign, file = file.path(path, "/metrics/",i,"metrics_benign.txt"), sep="\t", row.names=F, col.names = T, quote = F)
write.table(metrics_dominant, file = file.path(path, "/metrics/",i,"metrics_dominant.txt"), sep="\t", row.names=F, col.names = T, quote = F)
write.table(metrics_recessive, file = file.path(path, "/metrics/",i,"metrics_recessive.txt"), sep="\t", row.names=F, col.names = T, quote = F)

saveRDS(object = fit_3way, file = file.path(path, "/models/",i,"fit_3way.RDS"))
saveRDS(object = fit_scale, file = file.path(path, "/models/",i,"fit_scale.RDS")) }

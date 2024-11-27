library(caret)
library(ranger)
library(pROC)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(knitr)


setwd("tarabrennan/masters_proj/final_data/")

meta_file <- "meta_fil.csv"
asv_file <- "asv_relativeabundance.csv"

#add in desired metavariable variables here:
meta_vars <- c("age_dental_visit",
               "bmi", "gender.1.female.2.male", "daily.alcohol.consumption..1.2..3.from.drkst12r.0.less.than.daily..4.1.daily.",
               "avepa..av.hours.week.exercise.", "smokstat", "brush", "floss", "misspr..missing.teeth.", "dentist")

#specify desired response variables here:
PDresponse_vars <- paste0("PD", 1:8)
response_vars <- c(PDresponse_vars, "PD5_BoP20", "PD5_BoP40", "PD6_BoP20", "PD6_BoP40", 
                   "PD7_BoP20", "PD7_BoP40", "PD8_BoP20","PD8_BoP40")

prepare_data <- function(meta_file, asv_file, meta_vars, response_vars) {
  meta_data <- read.csv(meta_file, row.names = 1)
  asv_data <- read.csv(asv_file, row.names = 1)
  varsfrommeta = c(meta_vars,response_vars)
  
  meta_data <- meta_data[, varsfrommeta, drop = FALSE]
  combined_data <- merge(asv_data, meta_data, by = 0)
  
  combined_data <- na.omit(combined_data)
  rownames(combined_data) <- combined_data$Row.names
  combined_data <- combined_data[, !colnames(combined_data) %in% c("Row.names")]
  
  response = combined_data[, colnames(combined_data) %in% response_vars]
  data = combined_data[,!colnames(combined_data) %in% response_vars]
  
  if (!all(rownames(response) == rownames(data))) {
    stop("Row names don't match")
  }
    return(list(predictor_data = data, response_data = response))
}

modeldata <- prepare_data(meta_file, asv_file, meta_vars, response_vars)
predictor_data <- modeldata$predictor_data
response_data = modeldata$response_data

############## run model ################
roc_list <- list()

for (response_var in colnames(response_data)) {
  response <- factor(ifelse(response_data[[response_var]] == 0, 'no', 'yes'))
  data <- data.frame(response = response, predictor_data)
  
  train_index <- createDataPartition(data$response, p = 0.75, list = FALSE)
  training_set <- data[train_index, ]
  testing_set <- data[-train_index, ]
  
  train_control <- trainControl(
    method = "cv",
    number = 5,
    sampling = "up",
    classProbs = TRUE,
    summaryFunction = twoClassSummary
  )
  
  tune_grid <- expand.grid(
    mtry = c(10, 20, 30),
    min.node.size = c(5, 10, 15),
    splitrule = c("gini", "extratrees")
  )
  
  cv_model <- train(
    response ~ .,
    data = training_set,
    method = "ranger",
    trControl = train_control,
    #tuneLength = 1,
    tuneGrid = tune_grid,
    importance = 'impurity',
    metric = "ROC"
  )
  
  print(response_var)
  print(cv_model)
  print(cv_model$results)
  
  predictions <- predict(cv_model, newdata = testing_set, type = "prob")
  
  roc_obj <- roc(
    response = testing_set$response,
    predictor = predictions$yes,
    levels = rev(levels(testing_set$response))
  )
  
  roc_list[[response_var]] <- roc_obj
}

####### roc plots and tables #########
roc_data <- bind_rows(
  lapply(names(roc_list), function(var_name) {
    roc_curve <- roc_list[[var_name]]
    data.frame(
      tpr = roc_curve$sensitivities,
      fpr = 1 - roc_curve$specificities,
      variable = var_name
    )
  })
)

roc_data_sub = roc_data[roc_data$variable %in% colnames(response_data),] 
                                                 
desired_order <- c(unique(roc_data_sub$variable))
roc_data_sub$variable <- factor(roc_data_sub$variable, levels = desired_order)
colors <- colorRampPalette(brewer.pal(12, "Paired"))(length(desired_order))

rocplot <- ggplot(roc_data_sub, aes(x = fpr, y = tpr, color = variable)) +
  geom_line(linewidth = 0.5) +
  scale_color_manual(values = colors) +
  labs(
    title = "ROC Curves for Clinical Targets: All Meta Variables",
    x = "False Positive Rate (1 - Specificity)",
    y = "True Positive Rate (Sensitivity)",
    color = "Variable"
  ) +
  theme_minimal() +
  theme(legend.position = "right") 
ggsave(rocplot, file = "/Users/tarabrennan/masters_proj/final_figures/allmeta_RF_ROCcurve.pdf", width = 8, height = 6)


auc_values <- sapply(roc_list, function(roc_obj) {
  auc(roc_obj)
})
auc_table <- data.frame(t(auc_values))
colnames(auc_table) <- colnames(response_data)
rownames(auc_table) <- "AUC"
kable(auc_table, format = "markdown", digits = 3)



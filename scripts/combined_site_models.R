library(rsample)      # For data splitting
library(caret)        # For training control and modeling
library(Cubist)       # For Cubist models
library(iml)          # For feature importance
library(Metrics)      # For RMSE calculation
library(dplyr)        # For data manipulation
library(doParallel)   # For parallel processing
library(caretEnsemble)  # Caret ensemble

# Load the configuration file
source("../config.R")

# Set up parallel processing
cl <- makePSOCKcluster(parallelly::availableCores() - 1)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths())

# Data loader function
data_loader <- function(return_max_months = FALSE, return_with_site = TRUE) {
  p_val <- 0.75
  hyytiala <- within(read.csv(file = paste0(base_data_path, "smear2008-18_allmonths.csv")), rm(X, SoilTempB, DiffuseFract))
  varrio  <- within(read.csv(file = paste0(base_data_path, "varrio_2015-2019_allmonths.csv")), rm(X, DiffuseFract))
  hyytiala_2019 <- within(read.csv(file = paste0(base_data_path, "smear2019-21_allmonths.csv")), rm(X, SoilTempB, DiffuseFract))
  
  if (return_with_site) {
    hyytiala$site      <- 'hyytiala'
    varrio$site        <- 'varrio'
    hyytiala_2019$site <- 'hyytiala_thinning'
  }
  
  split_varrio <- initial_split(varrio, prop = p_val)
  train_var <- training(split_varrio)
  test_var  <- testing(split_varrio)
  
  split_hyytiala <- initial_split(hyytiala, prop = p_val)
  train_hyy <- training(split_hyytiala)
  test_hyy  <- testing(split_hyytiala)
  
  split_hyytiala_2019 <- initial_split(hyytiala_2019, prop = p_val)
  train_hyy_n <- training(split_hyytiala_2019)
  test_hyy_n  <- testing(split_hyytiala_2019)
  
  train <- rbind(train_hyy, train_var, train_hyy_n)
  test  <- rbind(test_hyy, test_var, test_hyy_n)
  
  if (return_with_site) {
    dummy  <- dummyVars(" ~ .", data=within(train, rm(Time)))
    train  <- data.frame(predict(dummy, newdata = within(train, rm(Time))), Time = train$Time)
    test   <- data.frame(predict(dummy, newdata = within(test,  rm(Time))), Time = test$Time)
  }
  
  if (return_max_months) {
    train <- train[as.numeric(strftime(train$Time, "%m")) %in% 7:8,]
    test  <- test[as.numeric(strftime(test$Time, "%m")) %in% 7:8,]
  }
  
  if (!return_with_site) {
    train <- train %>% select(-starts_with("site"))
    test  <- test  %>% select(-starts_with("site"))
  }
  
  return(list(train = train, test = test))
}

# Hyperparameters
hyperparams_with_site <- list(
  cubist = list(committees = 100, neighbors = 9),
  rf =     list(mtry = 11, min.node.size = 5),
  avNNet = list(size = 13, decay = 0.1, bag = FALSE)
)

hyperparams_without_site <- list(
  cubist = list(committees = 100, neighbors = 6),
  rf =     list(mtry = 8, min.node.size = 5),
  avNNet = list(size = 13, decay = 0.1, bag = FALSE)
)

# Feature importance function
feature_importance <- function(model, y, X) {
  mod <- Predictor$new(model, data = X, y = y)
  imp <- FeatureImp$new(mod, loss = "rmse", n.repetitions = 4)
  return(imp)
}

ale_calculate <- function(model, X) {
  mod <- Predictor$new(model, data = X)
  ale <- FeatureEffects$new(mod)
  return(ale)
}

fi_df <- data.frame()
ale_df <- data.frame()

# Train models function
train_models <- function(train, test, model_names, seed, return_with_site, name) {
  set.seed(seed)
  
  trainControl <- trainControl(method="repeatedcv", 
                               number=10,
                               repeats=5,
                               index = createFolds(train$NEE, 10),
                               savePredictions="final",
                               allowParallel = TRUE)
  
  params <- if (return_with_site) hyperparams_with_site else hyperparams_without_site
  
  modelTypes <- list(
    cubist = caretModelSpec(method="cubist", tuneGrid=expand.grid(committees=params$cubist$committees, neighbors=params$cubist$neighbors)),
    rf = caretModelSpec(method="ranger", tuneGrid=expand.grid(mtry=params$rf$mtry, min.node.size=params$rf$min.node.size, splitrule="extratrees")),
    avNNet = caretModelSpec(method="avNNet", MaxNWts=13 * (ncol(train)) + 13 + 1, maxit=10000, linout=TRUE, 
                            tuneGrid=expand.grid(size=params$avNNet$size, decay=params$avNNet$decay, bag=params$avNNet$bag)),
    lm = caretModelSpec(method="lm")
  )
  
  models <- caretList(
    NEE~., data=within(train, rm(Time)),
    trControl=trainControl,
    metric = "RMSE",
    tuneList = modelTypes
  )
  
  results_df <- data.frame(value = numeric(), score_f = character(), seed = integer(), t = character(), model = character())
  
  for (j in model_names) {
    predicted <- predict(models[j], test)
    rmse_test <- rmse(test$NEE, predicted)
    r2_test <- R2(test$NEE, predicted)
    
    predicted <- predict(models[j], train)
    rmse_train <- rmse(train$NEE, predicted)
    r2_train <- R2(train$NEE, predicted)
    
    results_df <- rbind(results_df, data.frame(value = rmse_test, score_f = "RMSE", seed = seed, t = "Test", model = j))
    results_df <- rbind(results_df, data.frame(value = r2_test, score_f = "R2", seed = seed, t = "Test", model = j))
    results_df <- rbind(results_df, data.frame(value = rmse_train, score_f = "RMSE", seed = seed, t = "Train", model = j))
    results_df <- rbind(results_df, data.frame(value = r2_train, score_f = "R2", seed = seed, t = "Train", model = j))
    
    fi <- feature_importance(models[j], test$NEE, within(test, rm("NEE", "Time")))$results
    fi$seed <- seed
    fi$model <- j
    fi$site <- name
    fi_df <<- rbind(fi_df, fi)
    
    ale_calculate <- ale(models[j], within(test, rm("NEE", "Time")))$results
    ale$seed <- seed
    ale$model <- j
    ale$site <- names[i]
    ale_df <<- rbind(ale_df, ale)  
  }
  
  return(list(models_list = models, results_df = results_df))
}

# Main execution
data_configs <- list(
  list(return_max_months = FALSE, return_with_site = TRUE),
  list(return_max_months = TRUE, return_with_site = TRUE),
  list(return_max_months = FALSE, return_with_site = FALSE),
  list(return_max_months = TRUE, return_with_site = FALSE)
)

names_param <- c("all_with_site", "max_with_site", "all_without_site", "max_without_site")
model_names_param <- c("cubist", "rf", "avNNet", "lm")
seeds <- c(1234, 5678, 91011, 34569, 9898)

results <- list()
all_results_df <- data.frame(value = numeric(), score_f = character(), seed = integer(), model = character(), config = character())

for (i in 1:length(data_configs)) {
  config <- data_configs[[i]]
  cat("Running for config:", names_param[i], "\n")
  
  for (seed in seeds) {
    cat("Running for seed:", seed, "\n")
    start.time <- Sys.time()
    
    data <- data_loader(return_max_months = config$return_max_months, return_with_site = config$return_with_site)
    train <- data$train
    test <- data$test
    
    trained_models <- train_models(train, test, model_names_param, seed, config$return_with_site, names_param[i])
    end.time <- Sys.time()
    time.taken <- round(end.time - start.time, 2)
    
    cat("Time taken:", time.taken, "\n")
    
    results[[paste(names_param[i], seed, sep = "_")]] <- list(trained_models$models_list, time.taken)
    all_results_df <- rbind(all_results_df, cbind(trained_models$results_df, config = names_param[i]))
  }
}

# Save results
save(results, file = paste0(base_output_path, "trained_models_combined.RData"))
write.csv(all_results_df, file = paste0(base_output_path, "all_results_df_combined.csv"), row.names = FALSE)
write.csv(fi_df, file = paste0(base_output_path, "fi_combined.csv"), row.names = FALSE)
write.csv(ale_df, file = paste0(base_output_path, "ale_combined.csv"), row.names = FALSE)


# Output time taken for each configuration and seed
sink(file = paste0(base_output_path, "test_output_combined.txt"))
for (result_name in names(results)) {
  cat("Time taken for config", result_name, ":", results[[result_name]][[2]], "\n")
}
sink(file = NULL)

# Stop the parallel cluster
stopCluster(cl)
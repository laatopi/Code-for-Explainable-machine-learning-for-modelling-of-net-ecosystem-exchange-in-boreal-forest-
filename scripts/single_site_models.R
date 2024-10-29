library(rsample)      # For data splitting
library(caret)        # For training control and modeling
library(Cubist)       # For Cubist models
library(iml)          # For feature importance
library(Metrics)      # For RMSE calculation
library(dplyr)        # For data manipulation
library(doParallel)   # For parallel processing
library(caretEnsemble)  # Caret ensemble

# Load the configuration file
source("./config.R")

# Set up parallel processing
cl <- makePSOCKcluster(parallelly::availableCores() - 1)
registerDoParallel(cl)
clusterEvalQ(cl, .libPaths())

# Load datasets
hyytiala_all <- within(read.csv(file = paste0(base_data_path, "smear2008-18_allmonths.csv")), rm(X))
varrio_all  <- within(read.csv(file = paste0(base_data_path, "varrio_2015-2019_allmonths.csv")), rm(X))
hyytiala_2019_all <- within(read.csv(file = paste0(base_data_path, "smear2019-21_allmonths.csv")), rm(X))

hyytiala_max <- within(read.csv(file = paste0(base_data_path, "smear2008-18_maxgrow.csv")), rm(X))
varrio_max  <- within(read.csv(file = paste0(base_data_path, "varrio_2015-2019_maxgrow.csv")), rm(X))
hyytiala_2019_max <- within(read.csv(file = paste0(base_data_path, "smear2019-21_maxgrow.csv")), rm(X))

# Hyperparameters for each dataset
hyperparams <- list(
  hyytiala_all = list(
    cubist = list(committees = 100, neighbors = 9),
    rf =     list(mtry = 3, min.node.size = 5),
    avNNet = list(size = 13, decay = 0.1, bag = FALSE)),
  
  varrio_all = list(
    cubist = list(committees = 100, neighbors = 6),
    rf =     list(mtry = 6, min.node.size = 5),
    avNNet = list(size = 13, decay = 0.1, bag = FALSE)),
  
  hyytiala_max = list(
    cubist = list(committees = 90, neighbors = 9),
    rf =     list(mtry = 3, min.node.size = 5),
    avNNet = list(size = 13, decay = 0.1, bag = FALSE)),
  
  varrio_max = list(
    cubist = list(committees = 100, neighbors = 3),
    rf =     list(mtry = 2, min.node.size = 5),
    avNNet = list(size = 13, decay = 0.1, bag = FALSE))
)

# Feature importance function
feature_importance <- function(model, y, X) {
  mod <- Predictor$new(model, data = X, y = y)
  imp <- FeatureImp$new(mod, loss = "rmse", n.repetitions = 3)
  return(imp)
}

fi_df <- data.frame()

# Function to perform data splitting and model training
train_models <- function(data, names, model_names, seed, p_val, hyperparams) {
  
  set.seed(seed)
  
  split <- initial_split(data[[1]], prop = p_val, strata = NEE)
  train_hyy_all <- training(split)
  test_hyy_all  <- testing(split)
  hyy_all <- list(train_hyy_all, test_hyy_all)
  
  split <- initial_split(data[[2]], prop = p_val, strata = NEE)
  train_var_all <- training(split)
  test_var_all  <- testing(split)
  var_all <- list(train_var_all, test_var_all)
  
  split <- initial_split(data[[3]], prop = p_val, strata = NEE)
  train_hyy_max <- training(split)
  test_hyy_max  <- testing(split)
  hyy_max <- list(train_hyy_max, test_hyy_max)
  
  split <- initial_split(data[[4]], prop = p_val, strata = NEE)
  train_var_max <- training(split)
  test_var_max  <- testing(split)
  var_max <- list(train_var_max, test_var_max)
  
  data_param <- list(hyy_all, var_all, hyy_max, var_max)
  
  models_list <- list()
  results_df <- data.frame(value = numeric(), score_f = character(), seed = integer(), t = character())
  
  for (i in 1:4) {
    
    train <- data_param[[i]][[1]]
    test  <- data_param[[i]][[2]]
    
    trainControl <- trainControl(method="repeatedcv", 
                                 number=10,
                                 repeats=5,
                                 index = createFolds(train$NEE, 10),
                                 savePredictions="final",
                                 allowParallel = TRUE)
    
    params <- hyperparams[[names[i]]]
    
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
    
    for (j in model_names) {
      predicted <- predict(models[j], test)
      rmse_test <- rmse(test$NEE, predicted)
      r2_test <- R2(test$NEE, predicted)
      
      predicted <- predict(models[j], train)
      rmse_train <- rmse(train$NEE, predicted)
      r2_train <- R2(train$NEE, predicted)
      
      results_df <- rbind(results_df, data.frame(value = rmse_test, score_f = "RMSE", seed = seed, t = "Test", model = j, site = names[i]))
      results_df <- rbind(results_df, data.frame(value = r2_test, score_f = "R2", seed = seed, t = "Test", model = j, site = names[i]))
      results_df <- rbind(results_df, data.frame(value = rmse_train, score_f = "RMSE", seed = seed, t = "Train", model = j, site = names[i]))
      results_df <- rbind(results_df, data.frame(value = r2_train, score_f = "R2", seed = seed, t = "Train", model = j, site = names[i]))
      
      fi <- feature_importance(models[j], test$NEE, within(test, rm("NEE", "Time")))$results
      fi$seed <- seed
      fi$model <- j
      fi$site <- names[i]
      
      fi_df <<- rbind(fi_df, fi)
      
      print(fi_df)
    }
    
    models_list[names[[i]]] <- list(models)
    
  }
  return(list(models_list = models_list, results_df = results_df))
}

# Main execution
data_param <- list(hyytiala_all, varrio_all, hyytiala_max, varrio_max)
names_param <- c("hyytiala_all", "varrio_all", "hyytiala_max", "varrio_max")
model_names_param <- c("cubist", "rf", "avNNet", "lm")
p_val = .8

# Run for different seeds
results <- list()
all_results_df <- data.frame(value = numeric(), score_f = character(), seed = integer(), model = character())
seeds <- c(1234, 5678, 91011, 404, 742)

for (seed in seeds) {
  cat("Running for seed:", seed, "\n")
  start.time <- Sys.time()
  trained_models <- train_models(data_param, names_param, model_names_param, seed, p_val, hyperparams)
  end.time <- Sys.time()
  time.taken <- round(end.time - start.time, 2)
  
  results[[as.character(seed)]] <- list(trained_models$models_list, time.taken)
  all_results_df <- rbind(all_results_df, trained_models$results_df)
}

# Save results to output directory
save(results, file = paste0(base_output_path, "trained_models_seeds.RData"))
write.csv(all_results_df, file = paste0(base_output_path, "all_results_df.csv"), row.names = FALSE)
write.csv(fi_df, file = paste0(base_output_path, "fi_seed.csv"), row.names = FALSE)

# Output time taken for each seed
sink(file = paste0(base_output_path, "test_output_seeds.txt"))
for (seed in seeds) {
  cat("Time taken for seed", seed, ":", results[[as.character(seed)]][[2]], "\n")
}
sink(file = NULL)

# Stop the parallel cluster
stopCluster(cl)
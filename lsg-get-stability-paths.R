# Solution - Stability

# Load already implemented functions
source("get-stability-paths.R")

# Implementation of missing functions:
############## resample ########################################################
sample_without_replacement <- function(nrows, strata = NULL, fraction = 0.5) {
  # @nrows (count): sample size
  # @strata (vector or NULL): Division into strata
  # @fraction (floating point between [0,1]): Fraction used for stratification
  
  # Calculates the random sample without replayement 
  # (implementation follows: sample_with_replacement)
  
  # Return: vector of drawn sample
  
  if (is.null(strata)) {
    return(sample(nrows, size = ceiling(fraction * nrows), replace = FALSE)) 
    # --> early exit!
  }
  rows <- tapply(
    X = seq_len(nrows), 
    INDEX = strata, FUN = sample_strata_without_replacement,
    fraction = fraction
  )
  as.vector(rows)
}

sample_strata_without_replacement <- function(nrows, fraction) {
  # @nrows (count): entire sample size
  # @fraction (floating point between [0,1]): Fraction used for stratification
  
  # Draws ceiling(fraction * length(nrows)) many samples within one strata
  
  # Return (vector): Drawn samples
  return(sample(nrows, size = ceiling(fraction * length(nrows)),
                replace = FALSE))
}
############## get_selected ####################################################
get_selected <- function(new_model) {
  # @new_model (regsubsets): Fitted model for each penalization
  
  # Calculates the selected covariate for each penalization
  
  # Return (matrix):  Columns corresponds to the covariates
  #                   Rows corresponds to the used penalization
  #                   (0) not selected \\ (1) selected
  
  summary_outmat <- as.matrix(summary(new_model)[["outmat"]])
  
  # Change values to numeric by changing "*" to 1 and "" to 0
  selected_mat <- structure(vapply(summary_outmat, FUN = change_to_numeric,
                                   numeric(1)), dim = dim(summary_outmat))
  
  # Adding a zero row
  selected_mat <- rbind(rep(0, dim(selected_mat)[[2]]), selected_mat)
  
  # Adding penalization for row names, obtained by the new_model attributes
  rownames(selected_mat) <- 
    c(0, attributes(summary(new_model)[["which"]])[["dimnames"]][[1]])
  # Adding covariates for column names
  colnames(selected_mat) <- colnames(summary_outmat)
  
  return(selected_mat)
}

change_to_numeric <- function(char) {
  # @char (character): "*" or " "
  # switch "*" -> 1 and " " -> 0
  # Returns 0 (covariate was  not selected) or 1 (covariate was selected)
  
  if (identical(char, "*")) {
    return(1) # --> early exit!
  }
  return(0)
}


############## make_paths ######################################################
make_paths <- function(selected) {
  # @selected (list of matrices): Each matrice:
  #                   Columns corresponds to the covariates
  #                   Rows corresponds to the used penalization
  #                   (0) not selected \\ (1) selected
  
  # Calculates for each penalization and covariate the percentage how often it 
  # was selected
  
  # Return (matrix):  Columns corresponds to the covariates
  #                   Rows corresponds to the used penalization
  #                   each entry gives the percentage how ofthen it was selected
  
  # The sum of all matrix entries for each matrix element, correspondes to how 
  # often the covariate was used for fitting with this method
  # https://stackoverflow.com/questions/11641701/sum-a-list-of-matrices
  selected_sum <- Reduce('+', selected)
  # calculating the proportion to number of repetitions
  path <- (1 / length(selected)) * selected_sum
  return(path)
}


############## plot_stability_paths ############################################
library(ggplot2)
plot_stability_paths <- function(stability_paths){
  # @stability_paths (matrix): Columns corresponds to the covariates
  #                   Rows corresponds to the used penalization
  #                   each entry gives the percentage how ofthen it was selected
  
  # Plots the stability path
  covariate <- names(sort(stability_paths[(nrow(stability_paths) - 1), ],
                          decreasing = TRUE))
  regularization <- rownames(stability_paths)
  plot_data <- data.frame(regularization = rep(regularization, 
                                          times = length(covariate)),
                     covariate = rep(covariate, each = length(regularization)),
                     path_values = as.vector(stability_paths))
  
  # Sort data
  plot_data[["covariate"]] <- factor(plot_data[["covariate"]],
                                     levels = covariate)
  
  

  ggplot(data = plot_data, aes(x = regularization, 
                               y = as.numeric(path_values),
                               group = covariate, colour = covariate)) +
    geom_line() + 
    geom_point() + 
    labs(x = "# covariate", y = expression(Pi)) + 
    theme(legend.title = element_blank()) + 
    scale_fill_discrete(breaks = levels(plot_data[["covariate"]]))
}

############## Test ############################################################
library(MASS)

data(prostate, package = "ElemStatLearn")
data <- prostate

max_formula <- lpsa ~ (. - train)
model <- leaps::regsubsets(max_formula,
                           data = data, nbest = 1, nvmax = 8,
                           really.big = TRUE
)


set.seed(20141020)
stability_paths <- get_stability_paths(model, data, reps = 1000)
stability_paths

plot_stability_paths(stability_paths)

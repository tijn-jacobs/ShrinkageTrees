#' ShrinkageTrees: Bayesian Tree Ensembles for Survival Analysis and Causal Inference
#'
#' Bayesian regression tree ensembles for survival analysis and causal
#' inference. Implements BART, DART, Bayesian Causal Forests (BCF), and
#' Horseshoe Forests models. Supports right-censored survival outcomes via
#' accelerated failure time (AFT) formulations. Designed for high-dimensional
#' prediction and heterogeneous treatment effect estimation in causal inference.
#'
#' @docType "_PACKAGE"
#' @name ShrinkageTrees-package
#' @aliases ShrinkageTrees
#'
#' @importFrom stats quantile
#' @importFrom utils head
NULL

# Suppress R CMD check NOTE for .data used inside ggplot2 aes() calls.
utils::globalVariables(".data")

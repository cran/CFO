#' Generate operating characteristics of phase I drug-combination trials in multiple simulations
#' 
#' Based on the toxicity outcomes, this function is used to conduct multiple simulations of phase I drug-combination trials and obtain relevant the operating characteristics.
#'
#' @usage CFO2d.oc(nsimu = 1000, target, p.true, init.level = c(1,1), ncohort, cohortsize,
#'                prior.para = list(alp.prior = target, bet.prior = 1 - target), 
#'                cutoff.eli = 0.95, early.stop = 0.95, seeds = NULL)
#'
#' @param nsimu the total number of trials to be simulated. The default value is 1000.
#' @param target the target DLT rate.
#' @param p.true a matrix representing the true DIL rates under the different dose levels.
#' @param init.level a numeric vector of length 2 representing the initial dose level (default is \code{c(1,1)}).
#' @param ncohort the total number of cohorts.
#' @param cohortsize the number of patients of each cohort. 
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of (\code{cutoff.eli = 0.95}) for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' @param seeds A vector of random seeds for each simulation, for example, \code{seeds = 1:nsimu} (default is \code{NULL}).
#'
#' @note In the example, we set \code{nsimu = 10} for testing time considerations. In reality, \code{nsimu} 
#'          is typically set to 1000 or 5000 to ensure the accuracy of the results.
#'
#' @return The \code{CFO.oc()} function returns basic setup of ($simu.setup) and the operating 
#'         characteristics of the design: \cr
#' \itemize{
#'   \item p.true: the matrix of the true DLT rates under the different dose levels.
#'   \item selpercent: the matrix of the selection percentage of each dose level.
#'   \item npatients: a matrix of the averaged number of patients allocated to different doses in one simulation.
#'   \item ntox: a matrix of the averaged number of DLT observed for different doses in one simulation.
#'   \item MTDsel: the percentage of the correct selection of the MTD.
#'   \item MTDallo: the averaged percentage of patients assigned to the target DLT rate.
#'   \item oversel: the percentage of selecting a dose above the MTD.
#'   \item overallo: the averaged percentage of patients assigned to dose levels with a DLT rate greater than the target.
#'   \item averDLT: the averaged total number of DLTs observed.
#'   \item percentstop: the percentage of early stopping without selecting the MTD.
#'   \item simu.setup: the parameters for the simulation set-up.
#' }
#' 
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Wang W, Jin H, Zhang Y, Yin G (2023). Two-dimensional calibration-free odds (2dCFO)
#'             design for phase I drug-combination trials. \emph{Frontiers in Oncology}, 13, 1294258.
#' @import pbapply
#' @export
#' @examples
#' ## Simulate a two-dimensional dose-finding trial with 20 cohorts of size 3 for 10 replications.
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#' 0.10, 0.15, 0.30, 0.45, 0.55,
#' 0.15, 0.30, 0.45, 0.50, 0.60), 
#' nrow = 3, ncol = 5, byrow = TRUE)
#' target <- 0.3; ncohort <- 12; cohortsize <- 3
#' CFO2doc <- CFO2d.oc(nsimu = 5, target, p.true, init.level = c(1,1), ncohort, cohortsize, 
#'                     seeds = 1:5)
#' summary(CFO2doc)
#' plot(CFO2doc)

CFO2d.oc <- function(nsimu = 1000, target, p.true, init.level = c(1,1), ncohort, cohortsize,
                     prior.para = list(alp.prior = target, bet.prior = 1 - target), 
                     cutoff.eli = 0.95, early.stop = 0.95, seeds = NULL){
  
  # Run the CFO2d.simu function nsimu times using lapply
  results <- pblapply(1:nsimu, function(i) {
    CFO2d.simu(target, p.true, init.level, ncohort, cohortsize, prior.para, cutoff.eli=cutoff.eli, early.stop=early.stop, seed = seeds[i])
  })
  
  selpercent <- matrix(0, dim(p.true)[1], dim(p.true)[2])
  overallo <- 0; sumPatients <- 0; sumTox <- 0
  for (i in 1:nsimu) {
    MTD <- results[[i]]$MTD
    if(MTD[1]<= dim(p.true)[1] & MTD[2] <= dim(p.true)[2]){
      selpercent[MTD] <- selpercent[MTD] + 1
    }
    
    if (MTD[1]== dim(p.true)[1] & MTD[2] == dim(p.true)[2]){
      overallo <- overallo
    }else{
      overallo <- overallo + sum(results[[i]]$npatients)*results[[i]]$ptoxic
    }
    sumTox <- sumTox + sum(results[[i]]$ntox)
    sumPatients <- sumPatients + sum(results[[i]]$npatients)
  }
  
  # Compute the average of the results
  nearlystop <- sum(sapply(results, `[[`, "earlystop"))
  avg_results <- list()
  avg_results$p.true <- p.true
  avg_results$selpercent <- selpercent / (nsimu)
  avg_results$npatients <- Reduce('+', lapply(results, `[[`, "npatients")) / (nsimu)
  avg_results$ntox <- Reduce('+', lapply(results, `[[`, "ntox")) / (nsimu)
  avg_results$MTDsel <- mean(sapply(results, `[[`, "correct"))
  avg_results$MTDallo <- mean(sapply(results, `[[`, "npercent"))
  avg_results$oversel <- sum(avg_results$selpercent[p.true > target])
  avg_results$overallo <- overallo/sumPatients
  avg_results$averDLT <- sumTox/sumPatients
  avg_results$percentstop <- mean(sapply(results, `[[`, "earlystop"))
  avg_results$simu.setup <- data.frame(target = target, ncohort = ncohort, cohortsize = cohortsize, design = "2dCFO", nsimu = nsimu)
  
  class(avg_results) <- c("cfo_oc","cfo")
  
  return(avg_results)
}


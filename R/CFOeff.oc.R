#' Generate operating characteristics of phase I/II trials single-drug trials in multiple simulations.
#' 
#' Based on the toxicity outcomes and efficacy outcomes, this function is used to perform multiple simulations for phase I/II single-drug trials and obtain relevant operating characteristics.
#' @usage CFOeff.oc(target, p.true=p.true, pE.true=pE.true, prior.para = 
#'                  list(alp.prior = target, bet.prior = 1 - target, 
#'                  alp.prior.eff = 0.5, bet.prior.eff = 0.5),  
#'                  init.level = 1, cohortsize=cohortsize, ncohort=ncohort, 
#'                  nsimu, cutoff.eli=0.95, 
#'                  early.stop=0.95, effearly.stop = 0.9, mineff,
#'                  seeds = NULL)
#' @param target the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param pE.true the true efficacy rates under the different dose levels.
#' @param prior.para the prior parameters for two beta distributions, where set as \code{list(alp.prior = target, 
#'                  bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5)} by default. \code{alp.prior} and \code{bet.prior} 
#'                  represent the parameters of the prior distribution for the true DLT rate at any dose level. This prior distribution 
#'                  is specified as Beta(\code{alpha.prior}, \code{beta.prior}). \code{alp.eff.prior} and \code{bet.eff.prior}
#'                  represent the parameters of the Jeffreys' prior distribution for the efficacy probability at any dose level.
#'                  This prior distribution is specified as Beta(\code{alpha.eff.prior}, \code{beta.eff.prior}).
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is 1.
#' @param cohortsize the number of patients in each cohort.
#' @param ncohort the total number of cohorts.
#' @param nsimu the total number of trials to be simulated.
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping due to overly toxic. The default value \code{early.stop = 0.95}
#'                   generally works well.
#' @param effearly.stop the threshold value for early stopping due to low efficacy. The trial would be terminated
#'                      early if \eqn{Pr(q_k<\psi |y_k,m_k \ge 3)} is smaller than the value of \code{effearly.stop} where \eqn{q_k, y_k} and \eqn{m_k}
#'                      are the efficacy probability, the number of efficacy outcomes and the number of patients at dose level \eqn{k}. 
#'                      \eqn{\psi} is the the lowest acceptable efficacy rate which is set by \code{mineff} here. 
#'                      By default, \code{effearly.stop} is set as \code{0.9}.
#' @param mineff the lowest acceptable efficacy rate.
#' @param seeds a vector of random seeds for each simulation, for example, \code{seeds = 1:nsimu} (default is \code{NULL}).
#' @note In the example, we set \code{nsimu = 3} for testing time considerations. 
#' In reality, \code{nsimu} is typically set as 5000 to ensure the accuracy of the results.
#' @return The \code{CFOeff.oc()} function returns a list object, which includes the basic setup (\code{simu.setup}), comprising the following components:
#'         \itemize{
#' \item p.true: the true DLT rates under the different dose levels.
#' \item pE.true: the true efficacy rates under the different dose levels.
#' \item selpercent: the selection percentage at each dose level.
#' \item npatients: the averaged number of patients treated at each dose level in one simulation.
#' \item ntox: the averaged number of toxicity observed at each dose level in one simulation.
#' \item neff: the averaged number of efficacy outcome at each dose level in one simulation.
#' \item OBDsel: the percentage of correct selection of the OBD.
#' \item OBDallo: the percentage of patients allocated to the OBD.
#' \item averDLT: the percentage of the patients suffering DLT.
#' \item avereff: the percentage of the patients with efficacy outcomes.
#' \item percentstop: the percentage of early stopping without selecting the OBD.
#' \item simu.setup: the parameters for the simulation set-up.
#' \item class: the phase of the trial.
#' }
#' 
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#' @importFrom dplyr transmute 
#' @import pbapply
#' @export
#' @examples 
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 12; nsimu = 3; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.05, 0.07, 0.1, 0.12, 0.16)
#' pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75)
#' result <- CFOeff.oc (target, p.true, pE.true, prior.para, 
#'           init.level,cohortsize, ncohort, nsimu, mineff = mineff, seeds = 1:nsimu)
#' summary(result)
#' plot(result)
#' \donttest{#earlystop for overly tox
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 12; nsimu = 3; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.75, 0.77, 0.81, 0.82, 0.86)
#' pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75)
#' result <- CFOeff.oc (target, p.true, pE.true, prior.para, 
#'           init.level,cohortsize, ncohort, nsimu, mineff = mineff, seeds = 1:nsimu)
#' summary(result)
#' plot(result)
#' }
#' \donttest{#earlystop for lower efficacy
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 20; nsimu = 3; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.05, 0.07, 0.1, 0.12, 0.16)
#' pE.true=c(0.001, 0.001, 0.001, 0.002, 0.003)
#' result <- CFOeff.oc (target, p.true, pE.true, prior.para, 
#'           init.level,cohortsize, ncohort, nsimu, mineff = mineff, seeds = 1:nsimu)
#' summary(result)
#' plot(result)
#' }
CFOeff.oc <- function(target, p.true=p.true, pE.true=pE.true, prior.para = 
                        list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5),  
                      init.level = 1, cohortsize=cohortsize, ncohort=ncohort, nsimu, cutoff.eli=0.95, 
                      early.stop=0.95, effearly.stop = 0.9, mineff,
                      seeds = NULL) {
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  OBD.level <- function(phi, mineff, p.true, pE.true) {
    if (p.true[1] > phi + 0.1) {
      OBD <- 99
      return(OBD)
    }
    OBD <- which.min(abs(phi - p.true))
    eff.idxs <- mineff > pE.true[1:OBD]
    if (sum(eff.idxs) == OBD) {
      OBD <- 99
      return(OBD)
    }
    
    OBD <- which.max(pE.true[1:OBD])
    return(OBD)
  }
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################  
  run.fn <- function(i) {
    res <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
                       prior.para, cutoff.eli, early.stop, effearly.stop,
                       mineff, seed = seeds[i])
    ress <- list(
      res = res,
      paras = list(p.true = p.true,
                   pE.true = pE.true,
                   obd = tobd, 
                   prior.para = prior.para,
                   target = target,
                   min_eff = mineff,
                   ncohort = ncohort,
                   cohortsize = cohortsize)
    )
    return(ress)
  }
  
  tobd <- OBD.level(target, mineff, p.true, pE.true)
  results <- pblapply(1:nsimu, run.fn)
  results_nopara <- lapply(1:nsimu, function(i)results[[i]]$res)
  ndose <- length(results_nopara[[1]]$npatients)
  
  Perc <- rep(0, ndose)
  nPatients <- rep(0, ndose); nTox <- rep(0, ndose); nEff <- rep(0, ndose)
  sumPatients <- 0; sumTox <- 0; sumEff <- 0
  nonErrStops <- 0
  OBDsel <- 0; OBDallo <- 0; oversel <- 0; overallo <- 0
  totaltime <- 0
  for (res in results_nopara){
    if (res$OBD != 99){
      nonErrStops <- nonErrStops + 1
      Perc[res$OBD] <- Perc[res$OBD] + 1
    }
    
    if (!is.null(res$totaltime)) {
      totaltime <- totaltime + res$totaltime
    }
    
    OBDsel <- OBDsel + sum(res$OBD == tobd)
    OBDallo <- OBDallo + res$npatients[tobd]
    sumTox <- sumTox + sum(res$ntox)
    sumEff <- sumEff + sum(res$neff)
    sumPatients <- sumPatients + sum(res$npatients)
    nPatients <- nPatients + res$npatients
    nTox <- res$ntox + nTox
    nEff <- res$neff + nEff
  }
  
  
  selpercent <- Perc/nsimu
  OBDsel <- OBDsel/nsimu
  OBDallo <- OBDallo/sumPatients
  averDLT <- sumTox/sumPatients
  averEff <- sumEff/sumPatients
  errStop <- nsimu - nonErrStops
  npatients = nPatients/nsimu
  ntox = nTox/nsimu
  neff = nEff/nsimu
  
  
  out <- list(p.true = p.true, pE.true = pE.true, selpercent = selpercent, npatients = npatients, ntox = ntox, neff = neff,
              OBDsel = OBDsel, OBDallo = OBDallo, averDLT = averDLT, avereff = averEff, percentstop = errStop/nsimu,
              simu.setup = data.frame(target = target,ncohort = ncohort, 
                                      cohortsize = cohortsize, nsimu = nsimu), class = "phaseI/II")
  class(out) <- c("cfo_eff_oc","cfo")
  return(out)
}
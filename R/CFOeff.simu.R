#' Conduct one simulation using the calibration-free odds (CFO) design for phase I/II trials
#'
#' In the CFO design for phase I/II trials, the function is used to conduct one single simulation and find the optimal biological dose (OBD).
#' @usage CFOeff.simu(target, p.true, pE.true, ncohort=10, init.level=1,  cohortsize=3,
#'                    prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                    alp.prior.eff = 0.5, bet.prior.eff = 0.5), 
#'                    cutoff.eli = 0.95, early.stop = 0.95, 
#'                    effearly.stop = 0.9, mineff, seed = NULL)
#' @param target the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param pE.true the true efficacy rates under the different dose levels.
#' @param ncohort the total number of cohorts.
#' @param init.level the dose level assigned to the first cohort. The default value of \code{init.level} is 1.
#' @param cohortsize the number of patients of each cohort.
#' @param prior.para the prior parameters for two beta distributions, where set as \code{list(alp.prior = target, 
#'                  bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5)} by default. \code{alp.prior} and \code{bet.prior} 
#'                  represent the parameters of the prior distribution for the true DLT rate at any dose level. This prior distribution 
#'                  is specified as Beta(\code{alpha.prior}, \code{beta.prior}). \code{alp.eff.prior} and \code{bet.eff.prior}
#'                  represent the parameters of the Jeffreys' prior distribution for the efficacy probability at any dose level.
#'                  This prior distribution is specified as Beta(\code{alpha.eff.prior}, \code{beta.eff.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                   the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping due to overly toxic. The default value \code{early.stop = 0.95}
#'                   generally works well.
#' @param effearly.stop the threshold value for early stopping due to low efficacy. The trial would be terminated
#'                      early if \eqn{Pr(q_k<\psi |y_k,m_k \ge 3)} is smaller than the value of \code{effearly.stop} where \eqn{q_k, y_k} and \eqn{m_k}
#'                      are the efficacy probability, the number of efficacy outcomes and the number of patients at dose level \eqn{k}. 
#'                      \eqn{\psi} is the the lowest acceptable efficacy rate which is set by \code{mineff} here. 
#'                      By default, \code{effearly.stop} is set as \code{0.9}.
#' @param mineff the lowest acceptable efficacy rate.
#' @param seed an integer to be set as the seed of the random number generator for reproducible results. The default value is set to \code{NULL}.
#' @note The \code{CFOeff.simu} function is designed to conduct a single CFO simulation for phase I/II trials. The dose elimination rule is the
#'          same as the case in phase I (refer to the function \code{CFO.simu}). As for early stopping rule, compared to the case of phase I, the rule
#'          in this case further considers the efficacy data to terminate the trial early if none of the admissible dose levels show adequate 
#'          efficacious effect.
#' @return  The \code{CFOeff.simu} function returns a list object comprising the following components:
#'         \itemize{
#'         \item OBD: the selected OBD. \code{OBD = 99} indicates that the simulation is terminated due to early stopping.
#'         \item target: the target DLT rate.
#'         \item npatients: the total number of patients allocated to all dose levels.
#'         \item neff: the total number of efficacy outcomes for all dose levels.
#'         \item ntox: the total number of DLTs observed for all dose levels.
#'         \item pE.true: the true efficacy rates under the different dose levels.
#'         \item p.true: the true DLT rates under the different dose levels.
#'         \item cohortdose: a vector including the dose level assigned to each cohort.
#'         \item ptoxic: the percentage of subjects assigned to dose levels with a DLT rate greater than the target.
#'         \item patientDLT: a vector including the DLT outcome observed for each patient.
#'         \item patienteff: a vector including the efficacy outcome observed for each patient.
#'         \item over.doses: a vector indicating whether each dose is overdosed or not (1 for yes).
#'         \item under.eff: a vector indicating whether the efficacy of each dose is lower than
#'               acceptable efficacy rate (1 for yes).
#'         \item correct: a binary indicator of whether the recommended dose level matches the correct OBD (1 for yes).
#'               The correct OBD is the dose level in the admissible set with the upper bound being the correct MTD, 
#'               which has the highest true efficacy probability.
#'         \item sumDLT: the total number of DLT observed.
#'         \item sumeff: the total number of efficacy outcome observed.
#'         \item earlystop: a binary indicator of whether the trial is early stopped (1 for yes).
#'         \item stopreason: the reason for earlystop. \code{overly_toxic} represents the trial was terminated 
#'         beacuse all tested doses were overly toxic. \code{low_efficacy} represents the trial was terminated 
#'         because all tested doses show low efficacy.
#'         \item class: the phase of the trial.
#'         }
#'
#' @author Jialu Fang, Wenliang Wang, Ninghao Zhang, and Guosheng Yin
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#' @export
#' 
#'
#' @examples 
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 20; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.05, 0.07, 0.1, 0.12, 0.16)
#' pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75)
#' result <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
#'                        prior.para, mineff = mineff, seed = 1)
#' summary(result)
#' plot(result)
#' \donttest{### overly toxic
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 20; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.55, 0.57, 0.61, 0.62, 0.66)
#' pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75)
#' result <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
#'                        prior.para, mineff = mineff, seed = 1)
#' summary(result)
#' plot(result)
#' }
#' \donttest{### low efficacy
#' target <- 0.30; mineff <- 0.30; cohortsize = 3; ncohort = 20; init.level = 1
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.05, 0.07, 0.1, 0.12, 0.16)
#' pE.true=c(0.001, 0.003, 0.004, 0.005, 0.006)
#' result <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
#'                        prior.para, mineff = mineff, seed = 1)
#' summary(result)
#' plot(result)
#' }
CFOeff.simu <- function(target, p.true, pE.true, ncohort=10, init.level=1, cohortsize=3,
                        prior.para = list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, 
                                          bet.prior.eff = 0.5), cutoff.eli=0.95, early.stop=0.95, 
                        effearly.stop = 0.9, mineff, seed = NULL) {
  
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  post.prob.fn <- function(target, y, n, alp.prior=0.1, bet.prior=0.1){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(target, alp, bet)
  }
  
  
  under.eff.fn <- function(mineff, effearly.stop,prior.para=list())
  {
    args <- c(list(target = mineff), prior.para)
    x <- prior.para$x
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior.eff
    bet.prior <- prior.para$bet.prior.eff
    ppE <- 1 - post.prob.fn(mineff, x, n, alp.prior, bet.prior)
    if ((ppE >= effearly.stop) & (n >= 3)) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
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
  
  # compute the marginal prob when lower < phiL/phiC/phiR < upper
  # i.e., Pr(Y=y|lower<target<upper)
  overdose.fn <- function(target, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(target, y, n, alp.prior, bet.prior)
    if ((pp >= threshold) & (prior.para$n >= 3)) {
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  
  moveprobs <- function(ad.xs, ad.ns, alp.prior, bet.prior){
    alps <- ad.xs + alp.prior
    bets <- ad.ns - ad.xs + bet.prior
    nd <- length(ad.xs)
    
    Nsps <- 10000
    sps.list <- list() 
    for (i in 1:nd){
      sps.list[[i]] <- rbeta(Nsps, alps[i], bets[i])
    }
    
    spss <- do.call(rbind, sps.list)
    argMaxs <- apply(spss, 2, which.max)
    probs <- as.vector(table(argMaxs))/Nsps
    
    probs
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################  
  # target: Target DIL rate
  # mineff: The minimal efficacy rate, only used for early stop
  # p.true: True DIL rates under the different dose levels
  # pE.true: True efficacy probs under the different dose levels
  # ncohort: The number of cohorts
  # cohortsize: The sample size in each cohort
  # alp.prior, bet.prior: prior parameters
  
  
  set.seed(seed)
  tadd.args <- prior.para
  earlystop <- 0
  stopreason <- NULL
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  cidx <- init.level
  
  tys <- rep(0, ndose) # number of DLT responses for different doses.
  txs <- rep(0, ndose) # number of efficacy responses for different doses.
  tns <- rep(0, ndose) # number of subject for different doses.
  tover.doses <- rep(0, ndose) # Whether each dose is too toxic or not, 1 yes.
  tunder.effs <- rep(0, ndose) # Whether the dose is not efficacious or not, 1 yes
  DLTlist <- c()
  EFFlist <- c()
  # if a dose is not efficacious enough or it is too toxic, it is would be eliminated from the admissible set.
  
  for (i in 1:ncohort) {
    pc <- p.true[cidx] 
    pEc <- pE.true[cidx] 
    doselist[i] <- cidx
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    cEres <- rbinom(cohortsize, 1, pEc)
    DLTlist <- c(DLTlist, cres)
    EFFlist <- c(EFFlist, cEres)
    
    # update results
    tys[cidx] <- tys[cidx] + sum(cres)#The number of observed DLT
    txs[cidx] <- txs[cidx] + sum(cEres)#The number of efficacy outcomes
    tns[cidx] <- tns[cidx] + cohortsize#The number of patient at dose level k
    
    cy <- tys[cidx]
    cx <- txs[cidx]
    cn <- tns[cidx]
    
    prior.para <- c(list(y = cy, n = cn, x = cx, tys = tys, txs = txs, tns = tns, cidx = cidx), tadd.args)
    
    if (overdose.fn(target, cutoff.eli, prior.para)) {
      tover.doses[cidx:ndose] <- 1
    }
    
    if (under.eff.fn(mineff, effearly.stop, prior.para)) {
      tunder.effs[cidx] <- 1
    }else{
      tunder.effs[cidx] <- 0
    }
    
    
    if (tover.doses[1] == 1) {
      stopreason <- "overly_toxic"
      earlystop <- 1
      break()
    }
    if (sum(tunder.effs[tover.doses == 0]) == sum(tover.doses == 0)){
      stopreason <- "low_efficacy"
      earlystop <- 1
      break()
      
    }
    
    nextinfo <- CFOeff.next(target, txs, tys, tns, cidx, prior.para, cutoff.eli, early.stop, effearly.stop, mineff)
    cidx <- nextinfo$nextdose
    if (cidx == 99){
      if (nextinfo$decision == "stop_for_tox"){
        stopreason <- "overly_toxic"
      }else{
        stopreason <- "low_efficacy"
      }
      earlystop <- 1
      break()
    }
    
    
  }
  
  
  if (earlystop == 0) {
    MTD <- CFO.selectmtd(target, tns, tys)$MTD
    if ( (MTD == 99) | (sum(tunder.effs[1:MTD]) == MTD)) {
      OBD <- 99
    }else{
      OBD.probs <- moveprobs(txs[1:MTD], tns[1:MTD], prior.para$alp.prior.eff, prior.para$bet.prior.eff)
      OBD <- which.max(OBD.probs)
    }
    
  }else{
    OBD <- 99
  }
  
  tobd = OBD.level(target, mineff, p.true, pE.true)
  correct <- 0
  if (OBD == tobd){
    correct <- 1
  }
  
  
  
  ptoxic <- sum(tns[which(p.true > target)])/(ncohort*cohortsize)
  
  out <- list(OBD = OBD, target = target, npatients = tns, neff = txs, ntox = tys, pE.true = pE.true, p.true = p.true, 
              cohortdose = doselist, ptoxic = ptoxic, patientDLT = DLTlist, patienteff = EFFlist, 
              over.doses = tover.doses, under.eff = tunder.effs, correct = correct, 
              sumDLT = sum(DLTlist), sumeff = sum(EFFlist), earlystop = earlystop, stopreason = stopreason, 
              class = "phaseI/II")
  class(out) <- c("cfo_eff_trial", "cfo")
  return(out)
  
}
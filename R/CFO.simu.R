#' Conduct one simulation using the Calibration-free odds (CFO) or accumulative CFO (aCFO) design.
#' 
#' In the CFO and aCFO designs, the function is used to conduct one single simulation and find the maximum tolerated dose (MTD).
#'
#' @usage CFO.simu(design, target, p.true, init.level = 1, ncohort, cohortsize,
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target), 
#'        cutoff.eli = 0.95, early.stop = 0.95, seed = NULL)
#'
#' @param design option for selecting different designs, which can be set as \code{'CFO'} and \code{'aCFO'}.
#' @param target the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is 1.
#' @param ncohort the total number of cohorts.
#' @param cohortsize the number of patients of each cohort.
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' @param seed an integer to be set as the seed of the random number generator for reproducible results. The default value is set to \code{NULL}.
#'                            
#' @note The \code{CFO.simu()} function is designed to conduct a single CFO or aCFO simulation. If \code{design = 'CFO'}, it corresponds 
#'          to the CFO design. If \code{design = 'aCFO'}, it corresponds to the aCFO design. \cr
#'          The early stopping and dose elimination rules are incorporated into the CFO or aCFO design 
#'          to ensure patient safety and benefit. If there is substantial evidence indicating that the current dose level 
#'          exhibits excessive toxicity, we exclude the current dose level as well as higher dose levels from the trial. If the lowest dose level is overly toxic, the trial will be terminated 
#'          according to the early stopping rule. Upon the predefined maximum sample size is reached or the lowest dose 
#'          level is over-toxicity, the experiment is concluded, and the MTD is determined using isotonic regression.
#'
#' @return The \code{CFO.simu} function returns a list object comprising the following components:
#'         \itemize{
#'         \item target: the target DLT rate.
#'         \item MTD: the selected MTD. \code{MTD = 99} indicates that the simulation is terminated due to early stopping.
#'         \item correct: a binary indicator of whether the recommended dose level matches the target DLT rate (1 for yes).
#'         \item npatients: the total number of patients allocated to all dose levels.
#'         \item ntox: the total number of DLTs observed for all dose levels.
#'         \item npercent: the percentage of subjects assigned to the target DLT rate.
#'         \item over.doses: a vector indicating whether each dose is overdosed or not (1 for yes).
#'         \item cohortdose: a vector including the dose level assigned to each cohort.
#'         \item ptoxic: the percentage of subjects assigned to dose levels with a DLT rate greater than the target.
#'         \item patientDLT: a vector including the DLT outcome observed for each patient.
#'         \item sumDLT: the total number of DLT observed.
#'         \item earlystop: a binary indicator of whether the trial is early stopped (1 for yes).
#'         }
#' 
#' @author Jialu Fang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#'
#' @examples
#' target <- 0.2; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <- c(0.01, 0.07, 0.20, 0.35, 0.50, 0.65, 0.80)
#' ### find the MTD for a single CFO simulation
#' CFOtrial <- CFO.simu(design = 'CFO', target, p.true, init.level, ncohort, cohortsize, seed = 1)
#' summary(CFOtrial)
#' plot(CFOtrial)
#' ### find the MTD for a single aCFO simulation
#' aCFOtrial <- CFO.simu(design = 'aCFO', target, p.true, init.level, ncohort, cohortsize, seed = 1)
#' summary(aCFOtrial)
#' plot(aCFOtrial)
#' @export
CFO.simu <- function(design, target, p.true, init.level=1, ncohort, cohortsize,
                     prior.para=list(alp.prior=target, bet.prior=1-target),
                     cutoff.eli=0.95, early.stop=0.95, seed=NULL){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.9){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  overdose.fn <- function(phi, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    # print(data.frame("prob of overdose" = pp))
    if ((pp >= threshold) & (prior.para$n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
  }
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  set.seed(seed)
  earlystop <- 0
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  currdose <- init.level
  
  ays <- rep(0, ndose) # number of responses for different doses.
  ans <- rep(0, ndose) # number of subject for different doses.
  tover.doses <- rep(0, ndose) # Whether each dose is overdosed or not, 1 yes
  DLTlist <- c()
  
  for (i in 1:ncohort){
    pc <- p.true[currdose]
    doselist[i] <- currdose
    
    # sample from current dose
    cres <- rbinom(cohortsize, 1, pc)
    DLTlist <- c(DLTlist, cres)
    
    # update results
    ays[currdose] <- ays[currdose] + sum(cres)
    ans[currdose] <- ans[currdose] + cohortsize
    
    cy <- ays[currdose]
    cn <- ans[currdose]
    
    prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
    
    if (overdose.fn(target, cutoff.eli, prior.para)){
      tover.doses[currdose:ndose] <- 1
    }
    
    if (currdose == 1){
      if (cutoff.eli != early.stop) {
        cy <- ays[1]
        cn <- ans[1]
        prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
        if (overdose.fn(target, early.stop, prior.para)){
          tover.doses[1:ndose] <- 1
        }
      }
    }
    
    
    if (tover.doses[1] == 1){
      earlystop <- 1
      break()
    }
    
    prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (design == 'CFO'){
      # the results for current 3 dose levels
      if (currdose!=1){
        cys <- ays[(currdose-1):(currdose+1)]
        cns <- ans[(currdose-1):(currdose+1)]
      }else{
        cys <- c(NA, ays[1:(currdose+1)])
        cns <- c(NA, ans[1:(currdose+1)])
      }
      currdose <- CFO.next(target, cys, cns, currdose, prior.para, cutoff.eli, early.stop)$nextdose
    }else if (design == 'aCFO'){
      currdose <- aCFO.next(target, ays, ans, currdose, prior.para, cutoff.eli, early.stop)$nextdose
    }else{
      stop("The input design is invalid; it can only be set as 'CFO' or 'aCFO'.")
    }
  }
  
  if (earlystop==0){
    MTD <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose=FALSE)$MTD
  }else{
    MTD <- 99
  }
  
  correct <- 0
  if (MTD == target){
    correct <- 1
  }
  
  npercent <- ans[which(p.true == target)]/(ncohort*cohortsize)
  ptoxic <- sum(ans[which(p.true > target)])/(ncohort*cohortsize)
  
  out<-list(target=target, MTD=MTD, correct=correct, npatients=ans, ntox=ays, npercent=npercent, 
            over.doses=tover.doses, cohortdose=doselist, ptoxic=ptoxic, patientDLT=DLTlist,
            sumDLT=sum(DLTlist), earlystop=earlystop)
  class(out) <- "cfo"
  return(out)
}

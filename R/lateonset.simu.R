#' Conduct one simulation using the calibration-free odds type (CFO-type) design with late-onset toxicity for phase I trials.
#'
#' Based on the toxicity outcomes of the enrolled cohorts, the function is used to conduct one single simulation and find the 
#' maximum tolerated dose (MTD) for the CFO-type designs with late-onset toxicities for phase I trials, specifically, 
#' including time-to-event CFO (TITE-CFO) design, fractional CFO (fCFO) design, benchmark CFO design, 
#' time-to-event accumulative CFO (TITE-aCFO) design, fractional accumulative CFO (f-aCFO) design and benchmark aCFO design.
#'
#' @usage lateonset.simu(design, target, p.true, init.level = 1, ncohort, cohortsize,
#'        assess.window, tte.para, accrual.rate, accrual.dist,  
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target), 
#'        cutoff.eli = 0.95, early.stop = 0.95, seed = NULL)
#'
#' @param design option for selecting different designs, which can be set as \code{'TITE-CFO'}, \code{'TITE-aCFO'}, 
#'               \code{'fCFO'}, \code{'f-aCFO'}, \code{'bCFO'}, and \code{'b-aCFO'}. Specifically, \code{'bCFO'} refers 
#'               to the benchmark CFO design, and \code{'b-aCFO'} denotes the benchmark aCFO design.
#' @param target the target DLT rate.
#' @param p.true the true DLT rates under the different dose levels.
#' @param ncohort the total number of cohorts.
#' @param cohortsize the number of patients of each cohort. 
#' @param assess.window the maximal assessment window size.
#' @param tte.para the parameter related with the distribution of the time to DLT events. The time to DLT is sampled from a Weibull 
#'                 distribution, with \code{tte.para} representing the proportion of DLTs occurring within the first half of the 
#'                 assessment window.
#' @param accrual.rate the accrual.rate rate, i.e., the number of patients accrued per unit time.
#' @param accrual.dist the distribution of the arrival times of patients. When \code{accrual.dist = 'fix'}, it corresponds to all 
#'                     patients in each cohort arriving simultaneously at a given accrual rate. When \code{accrual.dist = 'unif'}, 
#'                     it corresponds to a uniform distribution, and when \code{accrual.dist = 'exp'}, it corresponds to an 
#'                     exponential distribution.
#' @param init.level the dose level assigned to the first cohort. The default value \code{init.level} is 1.
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' @param seed an integer to set as the seed of the random number generator for reproducible results. The default value is set to NULL.
#' 
#' @note The early stopping and dose elimination rules are incorporated into the design 
#'        to ensure patient safety and benefit.
#'
#' @return The \code{lateonset.simu()} function returns a list object comprising the following components: 
#' \itemize{
#' \item target: the target DLT rate.
#' \item MTD: the selected MTD. \code{MTD = 99} indicates that this trial is terminated due to early stopping.
#' \item correct: a binary indicator of whether the recommended dose level matches the correct MTD (1 for yes).
#'       The correct MTD is the dose level at which the true DLT rate is closest to the target DLT rate.
#' \item npatients: the total number of patients allocated to all dose levels
#' \item ntox: the total number of DLTs observed for all dose levels.
#' \item over.doses: a vector indicating whether each dose is overdosed or not (1 for yes).
#' \item cohortdose: a vector including the dose level assigned to each cohort.
#' \item ptoxic: the percentage of subjects assigned to dose levels with a DLT rate greater than the target.
#' \item patientDLT: a vector including the DLT outcome observed for each patient.
#' \item sumDLT: the total number of DLT observed.
#' \item earlystop: a binary indicator of whether the trial is early stopped (1 for yes).
#' \item p_est: the isotonic estimate of the DLT probablity at each dose and associated \eqn{95\%} credible interval.
#'      \code{p_est = NA} if all tested doses are overly toxic.
#' \item p_overdose: p_overdose: the probability of overdosing defined as \eqn{Pr(toxicity > \code{target}|data)}.
#'      \code{p_overdose = NA} if all tested doses are overly toxic.
#' \item totaltime: the duration of the trial.
#' \item entertimes: the time that each participant enters the trial.
#' \item DLT.times: the time to DLT for each subject in the trial. If no DLT occurs for a certain subject, 
#'                  \code{DLT.times} is 0.
#' }
#' 
#'         
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials. 
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Jin H, Yin G (2023). Time‐to‐event calibration‐free odds design: A robust efficient design for 
#'             phase I trials with late‐onset outcomes. \emph{Pharmaceutical Statistics}. 22(5), 773–783.\cr
#'             Yin G, Zheng S, Xu J (2013). Fractional dose-finding methods with late-onset toxicity in 
#'             phase I clinical trials. \emph{Journal of Biopharmaceutical Statistics}, 23(4), 856-870. \cr
#'             Fang J, Yin G (2024). Fractional accumulative calibration‐free odds (f‐aCFO) design for delayed toxicity 
#'             in phase I clinical trials. \emph{Statistics in Medicine}.
#' @export
#'
#' @examples
#' target <- 0.2; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <- c(0.01, 0.07, 0.20, 0.35, 0.50, 0.65, 0.80)
#' assess.window <- 3; accrual.rate <- 2; tte.para <- 0.5; accrual.dist <- 'unif'
#' ## find the MTD for a single TITE-CFO simulation
#' TITECFOtrial <- lateonset.simu (design = 'TITE-CFO', target, p.true, init.level,  
#'                 ncohort, cohortsize, assess.window, tte.para, accrual.rate, accrual.dist, seed = 1)
#' summary(TITECFOtrial)
#' plot(TITECFOtrial)
#' ## find the MTD for a single TITE-aCFO simulation
#' TITEaCFOtrial <- lateonset.simu (design = 'TITE-aCFO', target, p.true, init.level,  
#'                 ncohort, cohortsize, assess.window, tte.para, accrual.rate, accrual.dist, seed = 1)
#' summary(TITEaCFOtrial)
#' plot(TITEaCFOtrial)
#' ## find the MTD for a single fCFO simulation
#' fCFOtrial <- lateonset.simu (design = 'fCFO', target, p.true, init.level,  
#'                 ncohort, cohortsize, assess.window, tte.para, accrual.rate, accrual.dist, seed = 1)
#' summary(fCFOtrial)
#' plot(fCFOtrial)
#' ## find the MTD for a single f-aCFO simulation
#' faCFOtrial <- lateonset.simu (design = 'f-aCFO', target, p.true, init.level,  
#'                 ncohort, cohortsize, assess.window, tte.para, accrual.rate, accrual.dist, seed = 1)
#' summary(faCFOtrial)
#' plot(faCFOtrial)
lateonset.simu <- function(design, target, p.true, init.level=1, ncohort, cohortsize,
                           assess.window, tte.para, accrual.rate, accrual.dist, 
                           prior.para=list(alp.prior=target, bet.prior=1-target), 
                           cutoff.eli=0.95, early.stop=0.95, seed=NULL){
  
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  
  # The function is to obtain the DLT results (with TITE) for each subject
  gen.tite<-function(n, pi, assess.window=1, alpha=0.5){
    #args:
    #   n: Num of subjects to generate
    #   pi: Target DLT rate, pi=Pr(T<=assess.window)
    #   assess.window: Maximal window size
    #   alpha: Parameter for generate time
    #Return:
    #   if no DLT, tox.t=0
    ############ subroutines ############
    weib<-function(n, pi, pihalft)
    {
      ## solve parameters for Weibull given pi=1-S(T) and phalft=1-S(T/2)
      alpha = log(log(1-pi)/log(1-pihalft))/log(2);
      lambda = -log(1-pi)/(assess.window^alpha);
      t = (-log(runif(n))/lambda)^(1/alpha);
      return(t);
    }
    ############ end of subroutines ############
    tox = rep(0, n);
    t.tox = rep(0, n);
    #### Weibull
    pihalft = alpha*pi;  # alpha*100% event in (0, 1/2T)
    t.tox = weib(n, pi, pihalft);
    tox[t.tox<=assess.window]=1;
    ntox.st = sum(tox);
    t.tox[tox==0]=0;
    return(list(tox=tox, t.tox=t.tox, ntox.st=ntox.st));
  }
  
  MTD.level <- function(phi, p.true){
    if (p.true[1]>phi+0.1){
      MTD <- 99
      return(MTD)
    }
    MTD <- which.min(abs(phi - p.true))
    return(MTD)
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
  ndose <- length(p.true)
  doselist <- rep(0, ncohort)
  
  earlystop <- 0
  enter.times <- NULL # enter time of each subject
  dlt.times <- NULL # dlt time of each subject
  dlts <- NULL # dlt event for each subject
  doses <- NULL # dose level for each subject
  current.t<- 0
  currdose <- init.level  #current dose level
  
  over.doses <- rep(0, ndose)
  
  for (i in 1:ncohort){
    curP <- p.true[currdose]
    doselist[i] <- currdose
    
    if (accrual.dist=='fix'){
      delta.times <- rep(0, cohortsize)
    }else if (accrual.dist == 'unif'){
      delta.times <- cumsum(c(0, runif(cohortsize-1, 0, 2/accrual.rate)))
    }else if (accrual.dist == 'exp'){
      delta.times <- cumsum(c(0, rexp(cohortsize-1, rate=accrual.rate)))
    }

    enter.times <- c(enter.times, current.t+delta.times)
    
    # obtain the results of the patients
    obscohort <- gen.tite(cohortsize, curP, alpha=tte.para, assess.window=assess.window);
    dlt.times <- c(dlt.times, obscohort$t.tox);
    dlts <- c(dlts, obscohort$tox);
    doses <- c(doses, rep(currdose, cohortsize));
    
    # Move to next cohort 
    if (i != ncohort){
      if (accrual.dist=='fix'){
        delta.time <- cohortsize/accrual.rate
      }else if (accrual.dist == 'unif'){
        delta.time <- runif(1, 0, 2/accrual.rate)
      }else if (accrual.dist == 'exp'){
        delta.time <- rexp(1, rate=accrual.rate)
      }
    }else{
      delta.time <- assess.window
    }
    current.t<- enter.times[length(enter.times)] + delta.time;
    
    if (design == 'bCFO' || design == 'b-aCFO'){
      current.t <- enter.times[length(enter.times)] + assess.window
      res <- lateonset.next(design, target, ndose, currdose, assess.window, enter.times, dlt.times, current.t, doses, 
                            prior.para, cutoff.eli, early.stop)
      over.doses <- res$over.doses
      overTox <- res$overTox
      current.t <- current.t + delta.time
    }else{
      res <- lateonset.next(design, target, ndose, currdose, assess.window, enter.times, dlt.times, current.t, doses, 
                            prior.para, cutoff.eli, early.stop)
      over.doses <- res$over.doses
      overTox <- res$overTox
    }
    
    if (over.doses[1] == 1){
      earlystop <- 1
      break()
    } else{
      currdose <- res$nextdose
    }
  }
  
  ans <- NULL
  ays <- NULL
  assess.t <- enter.times + assess.window
  y.raw <- (dlt.times!=0)*1
  for (j in 1:ndose){
    ans <- c(ans, sum(doses==j))
    ays <- c(ays, sum(y.raw[doses==j]))
  }
  
  result <- CFO.selectmtd(target, ans, ays, prior.para, cutoff.eli, early.stop, verbose = TRUE)
  p_est <- result$p_est
  p_overdose <- result$p_overdose
  
  if (earlystop==0){
    MTD <- result$MTD
  }else{
    MTD <- 99
  }
  
  tmtd <- MTD.level(target, p.true)
  correct <- 0
  if (MTD == tmtd){
    correct <- 1
  }
  
  ptoxic <- sum(ans[which(p.true > target)])/(ncohort*cohortsize)
  
  out <- list(target=target, MTD=MTD, correct=correct, npatients=ans, ntox=ays, 
              over.doses=over.doses, cohortdose=doselist, ptoxic=ptoxic, patientDLT = dlts,
              sumDLT=sum(dlts), earlystop=earlystop, p_est = p_est, p_overdose = p_overdose,
              totaltime=assess.t[length(assess.t)], entertimes=enter.times, DLTtimes=dlt.times)
  class(out) <- c("cfo_trial", "cfo")
  return(out)
}

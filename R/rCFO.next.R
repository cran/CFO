#' 
#' Determination of the dose level for next cohort in the randomized calibration-free odds (rCFO) design for phase I trials
#' 
#' In the rCFO design for phase I trials, the function is used to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.
#'
#' @usage rCFO.next(target, cys, cns, currdose, 
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target),
#'        cutoff.eli = 0.95, early.stop = 0.95, seed)
#'
#' @param target the target DLT rate.
#' @param cys the cumulative numbers of DLTs observed at the left, current, and right dose levels.
#' @param cns the cumulative numbers of patients treated at the left, current, and right dose levels.
#' @param currdose the current dose level.
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' @param seed an integer to be set as the seed of the random number generator for reproducible results. The default value is set to \code{NULL}.
#'
#' @details The original CFO design makes deterministic dose movement by constructing two odds ratios, \eqn{\pi_L =O_C/ \overline{O}_{L}}
#'          and \eqn{\pi_R =\overline{O}_{C}/ O_R}, and comparing them against thresholds \eqn{\gamma_L} and \eqn{\gamma_R}, respectively.
#'          The rCFO design introduces a randomization scheme, normalizes odds ratios, \eqn{\pi_L}, and \eqn{\pi_R} into probabilities, and constructs
#'          probabilities for dose escalation, de-escalation, and staying at the same dose.
#'          
#' @note    When the current dose level is the lowest or highest (i.e., at the boundary), the parts in \code{cys} and 
#'          \code{cns} where there is no data are filled with \code{NA}. \cr
#'          The dose level indicated by \code{overtox} and all the dose levels above experience over-toxicity, and these dose levels will be eliminated.
#'          
#' @return The \code{rCFO.next()} function returns a list object comprising the following elements:
#' \itemize{
#'   \item target: the target DLT rate.
#'   \item cys: the cumulative counts of DLTs observed at the left, current, and right dose levels.
#'   \item cns: the cumulative counts of patients treated at the left, current, and right dose levels.
#'   \item decision: the decision in the CFO design, where \code{left}, \code{stay}, and \code{right} represent the 
#'   movement directions, and \code{stop} indicates stopping the experiment.
#'   \item currdose: the current dose level.
#'   \item nextdose: the recommended dose level for the next cohort. \code{nextdose = 99} indicates that the trial is 
#'   terminated due to early stopping.
#'   \item overtox: the situation regarding which positions experience over-toxicity. The dose level indicated 
#'   by \code{overtox} and all the dose levels above experience over-toxicity. \code{overtox = NA} signifies that 
#'   the occurrence of over-toxicity did not happen.
#'   \item toxprob: the expected toxicity probability, \eqn{Pr(p_k > \phi | x_k, m_k)}, at the left, current, and 
#'   right dose levels, where \eqn{p_k}, \eqn{x_k}, and \eqn{m_k} is the dose-limiting toxicity (DLT) rate, the 
#'   numbers of observed DLTs, and the numbers of patients at dose level \eqn{k}. \code{NA} indicates that there 
#'   are no patients at the corresponding dose level.
#' }
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#' 
#' @examples
#' ## determine the dose level for the next cohort of new patients
#' cys <- c(0, 1, 0); cns <- c(3, 6, 0)
#' decision <- rCFO.next(target=0.2, cys=cys, cns=cns, currdose=3)
#' summary(decision)
#' 
#' cys <- c(NA, 3, 0); cns <- c(NA, 3, 0)
#' decision <- rCFO.next(target=0.2, cys=cys, cns=cns, currdose=1)
#' summary(decision)
#' 
#' cys <- c(0, 3, NA); cns <- c(3, 3, NA)
#' decision <- rCFO.next(target=0.2, cys=cys, cns=cns, currdose=7)
#' summary(decision)
#' 
#' @import stats
#' @export
rCFO.next <- function(target, cys, cns, currdose, prior.para=list(alp.prior=target, bet.prior=1-target),
                          cutoff.eli=0.95, early.stop=0.95, seed=NULL){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
    if(n != 0){
      alp <- alp.prior + y 
      bet <- bet.prior + n - y
      res <- 1 - pbeta(phi, alp, bet)
    }else{
      res <- NA
    }
    return(res)
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
  
  prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- bet.prior + n1 - y1
    bet2 <- bet.prior + n2 - y2
    fn.min <- function(x){
      dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2)) 
    }
    fn.max <- function(x){
      pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    const.min <- integrate(fn.min, lower=0, upper=0.999, subdivisions=1000, rel.tol = 1e-10)$value
    const.max <- integrate(fn.max, lower=0, upper=0.999, rel.tol = 1e-10)$value
    p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
    
    list(p1=p1, p2=p2)
  }
  
  OR.values <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior, type){
    ps <- prob.int(phi, y1, n1, y2, n2, alp.prior, bet.prior)
    if (type=="L"){
      pC <- 1 - ps$p2
      pL <- 1 - ps$p1
      oddsC <- pC/(1-pC)
      oddsL <- pL/(1-pL)
      OR <- oddsC*oddsL
      
    }else if (type=="R"){
      pC <- 1 - ps$p1
      pR <- 1 - ps$p2
      oddsC <- pC/(1-pC)
      oddsR <- pR/(1-pR)
      OR <- (1/oddsC)/oddsR
    }
    return(OR)
  }
  
  All.OR.table <- function(phi, n1, n2, type, alp.prior, bet.prior){
    ret.mat <- matrix(rep(0, (n1+1)*(n2+1)), nrow=n1+1)
    for (y1cur in 0:n1){
      for (y2cur in 0:n2){
        ret.mat[y1cur+1, y2cur+1] <- OR.values(phi, y1cur, n1, y2cur, n2, alp.prior, bet.prior, type)
      }
    }
    ret.mat
  }
  
  # compute the marginal prob when lower < phiL/phiC/phiR < upper
  # i.e., Pr(Y=y|lower<phi<upper)
  margin.phi <- function(y, n, lower, upper){
    if (upper > 1){upper <- 1}
    C <- 1/(upper-lower)
    fn <- function(phi) {
      dbinom(y, n, phi)*C
    }
    integrate(fn, lower=lower, upper=upper)$value
  }
  
  # Obtain the table of marginal distribution of (y1, y2) 
  # after intergrate out (phi1, phi2)
  # under H0 and H1
  # H0: phi1=phi, phi < phi2 < 2phi
  # H1: phi2=phi, 0   < phi1 < phi
  margin.ys.table <- function(n1, n2, phi, hyperthesis){
    if (hyperthesis=="H0"){
      p.y1s <- dbinom(0:n1, n1, phi)
      p.y2s <- sapply(0:n2, margin.phi, n=n2, lower=phi, upper=2*phi)
    }else if (hyperthesis=="H1"){
      p.y1s <- sapply(0:n1, margin.phi, n=n1, lower=0, upper=phi)
      p.y2s <- dbinom(0:n2, n2, phi)
    }
    p.y1s.mat <- matrix(rep(p.y1s, n2+1), nrow=n1+1)
    p.y2s.mat <- matrix(rep(p.y2s, n1+1), nrow=n1+1, byrow=TRUE)
    margin.ys <- p.y1s.mat * p.y2s.mat
    margin.ys
  }
  
  # Obtain the optimal gamma for the hypothesis test
  optim.gamma.fn <- function(n1, n2, phi, type, alp.prior, bet.prior){
    OR.table <- All.OR.table(phi, n1, n2, type, alp.prior, bet.prior) 
    ys.table.H0 <- margin.ys.table(n1, n2, phi, "H0")
    ys.table.H1 <- margin.ys.table(n1, n2, phi, "H1")
    
    argidx <- order(OR.table)
    sort.OR.table <- OR.table[argidx]
    sort.ys.table.H0 <- ys.table.H0[argidx]
    sort.ys.table.H1 <- ys.table.H1[argidx]
    n.tol <- length(sort.OR.table)
    
    if (type=="L"){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H0[1:i])
        err2 <- sum(sort.ys.table.H1[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
    }else if (type=='R'){
      errs <- rep(0, n.tol-1)
      for (i in 1:(n.tol-1)){
        err1 <- sum(sort.ys.table.H1[1:i])
        err2 <- sum(sort.ys.table.H0[(i+1):n.tol])
        err <- err1 + err2
        errs[i] <- err
      }
      min.err <- min(errs)
      if (min.err > 1){
        gam <- 0
        min.err <- 1
      }else {
        minidx <- which.min(errs)
        gam <- sort.OR.table[minidx]
      }
      
    }
    list(gamma=gam, min.err=min.err)
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
  
  cover.doses <- c(0,0,0)
  
  for (i in 1:3){
    cy <- cys[i]
    cn <- cns[i]
    if (is.na(cn)){
      cover.doses[i] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
      if (overdose.fn(target, cutoff.eli, prior.para)){
        cover.doses[i:3] <- 1
        break()
      }
    }
  }
  
  cover.prob <- c(0,0,0)
  for (i in 1:3){
    cy <- cys[i]
    cn <- cns[i]
    if (is.na(cn)){
      cover.prob[i] <- NA
    }else{
      cover.prob[i] <- post.prob.fn(target, cy, cn, alp.prior, bet.prior)
    }
  }
  
  idx <- if (currdose == 2) 1 else if (currdose == 1) 2 else NA
  if (!is.na(idx) & (cutoff.eli != early.stop)) {
    cy <- cys[idx]
    cn <- cns[idx]
    if (is.na(cn)){
      cover.doses[idx] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
      if (overdose.fn(target, early.stop, prior.para)){
        cover.doses[idx:3] <- 1
      }
    }
  }
  
  cover.doses <- ifelse(is.na(cys), NA, cover.doses)
  
  position <- which(cover.doses == 1)[1]
  overtox <- c(-1, 0, 1)[position] + currdose
  prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
  if ((cover.doses[2] == 1)&(currdose == 1)){
    index <- NA
    decision <- "stop"
  } else {
    if (cover.doses[2] == 1){
      index <- -1
      decision <- "de-escalation"
    }
    else{
      if (is.na(cys[1]) & (cover.doses[3]==1)){
        index <- 0
        decision <- "stay"
      }
      else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
        gam2 <- optim.gamma.fn(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma 
        OR.v2 <- OR.values(target, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        if (OR.v2>gam2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma 
        OR.v1 <- OR.values(target, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        if (OR.v1>gam1){
          index <- -1
          decision <- "de-escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], target, "L", alp.prior, bet.prior)$gamma 
        gam2 <- optim.gamma.fn(cns[2], cns[3], target, "R", alp.prior, bet.prior)$gamma 
        OR.v1 <- OR.values(target, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        OR.v2 <- OR.values(target, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        pL <- OR.v1/(OR.v1+OR.v2)
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          if (runif(1) < pL){
            index <- -1
            decision <- "de-escalation"
          } else {
            index <- 0
            decision <- "stay"
          }
        }else if (!v1 & v2){
          if (runif(1) < pL){
            index <- 0
            decision <- "stay"
          } else {
            index <- 1
            decision <- "escalation"
          }
        }else if (v1 & v2){
          if (OR.v1 == OR.v2){
            print('equal occur')
            index <- 0
            decision <- "stay"
          } else {
            if (runif(1) < pL){
              index <- -1
              decision <- "de-escalation"
            } else {
              index <- 1
              decision <- "escalation"
            }
          }
        }else{
          index <- 0
          decision <- "stay"
        }
      }
    }
  }
  
  if (decision=='stop'){
    nextdose <- 99
  }else{
    nextdose <- currdose+index
  }
  
  out <- list(target=target, cys=cys, cns=cns, decision=decision, currdose = currdose, 
              nextdose=nextdose, overtox=overtox, toxprob=cover.prob)
  class(out) <- c("cfo_decision", "cfo")
  return(out)
}
#' Determination of the dose level for next cohort in the calibration-free odds (CFO) design for phase I/II trials
#' 
#' In the CFO design for phase I/II trials, the function is used to determine the dose movement 
#' based on the toxicity outcomes and efficacy outcomes of the enrolled cohorts.
#' @usage CFOeff.next(target, axs, ays, ans, currdose, 
#'                    prior.para=list(alp.prior = target, bet.prior = 1 - target, 
#'                    alp.prior.eff = 0.5, bet.prior.eff = 0.5),  
#'                    cutoff.eli=0.95, early.stop=0.95, effearly.stop = 0.9, mineff)
#' @param target the target DLT rate.
#' @param axs the cumulative counts of efficacy outcomes at all dose levels.
#' @param ays the cumulative counts of DLTs observed at all dose levels.
#' @param ans the cumulative counts of patients treated at all dose levels.
#' @param currdose the current dose level.
#' @param prior.para the prior parameters for two beta distributions, where set as \code{list(alp.prior = target, 
#'                  bet.prior = 1 - target, alp.prior.eff = 0.5, bet.prior.eff = 0.5)} by default. \code{alp.prior} and \code{bet.prior} 
#'                  represent the parameters of the prior distribution for the true DLT rate at any dose level. This prior distribution 
#'                  is specified as Beta(\code{alpha.prior}, \code{beta.prior}). \code{alp.eff.prior} and \code{bet.eff.prior}
#'                  represent the parameters of the Jeffreys' prior distribution for the efficacy probability at any dose level.
#'                  This prior distribution is specified as Beta(\code{alpha.eff.prior}, \code{beta.eff.prior}).
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
#' 
#' @details
#' The CFO design for phase I/II trials will determine admissible set \eqn{A_n} through the dose escalation rules for the MTD. The current dose is set as 
#' \eqn{d_n}. If the decision is to de-escalate the dose, the set \eqn{A_n} will be \eqn{\{1,\dots,d_n-1\}}. If the decision is to stay at the 
#' current dose, then the admissible set \eqn{A_n} will be \eqn{\{1,\dots,d_n\}}. If the decision is to escalate the dose, then \eqn{A_n} will be
#' \eqn{\{1,\dots,d_n+1\}}. The dose level \eqn{d_{n+1}} for the next cohort will be selected from \eqn{A_n} by using the rule:
#' \eqn{d_{n+1} = argmax_{k\in A_n}Pr(q_k = max_{j\in A_n}\{q_j\}| D_n)} where \eqn{D_n} and \eqn{q_k} are the current data and the 
#' efficacy probability for dose level \eqn{k}.
#' 
#' @return The \code{CFOeff.next()} function returns a list object comprising the following elements:
#' \itemize{
#'   \item target: the target DLT rate.
#'   \item axs: the cumulative counts of efficacy outcomes at all dose levels.
#'   \item ays: the cumulative counts of DLTs observed at all dose levels.
#'   \item ans: the cumulative counts of patients treated at all dose levels.
#'   \item decision: the decision in the CFO design, where \code{de-escalation}, \code{stay}, and \code{escalation} represent the 
#'   movement directions of the dose level, \code{stop_for_tox} indicates stopping the experiment because the lowest dose level 
#'   is overly toxic and \code{stop_for_low_eff} indicates that all dose level in the admissible set shows low efficacy.
#'   \item currdose: the current dose level.
#'   \item nextdose: the recommended dose level for the next cohort. \code{nextdose = 99} indicates that the trial is 
#'   terminated due to early stopping.
#'   \item overtox: the situation regarding which positions experience over-toxicity. The dose level indicated 
#'   by \code{overtox} and all the dose levels above experience over-toxicity. \code{overtox = NA} signifies that 
#'   the occurrence of over-toxicity did not happen.
#'   \item toxprob: the expected toxicity probability, \eqn{Pr(p_k > \phi | x_k, m_k)}, for doses in admissible set,
#'   where \eqn{p_k}, \eqn{x_k}, and \eqn{m_k} are the dose-limiting toxicity (DLT) rate, the 
#'   numbers of observed DLTs, and the numbers of patients at dose level \eqn{k}.
#'   \item effprob: the empirical probability of \eqn{Pr(q_k=max_{j\in A_n}\{q_j\}|D_n)} for doses in admissible set, 
#'   where \eqn{q_k} is efficacy probability at dose level \eqn{k}. \eqn{A_n} is the admissible set determined through 
#'   the dose escalation rules for the MTD and \eqn{D_n} is the current cumulative dataset.
#'   \item admset: the admissible set \eqn{A_n}. The dose level for the next cohort will be selected from \eqn{A_n}.
#'   \item class: the phase of the trial.
#' }
#' 
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#' 
#' @export
#'
#' @examples 
#' axs = c(3, 1, 7, 11, 26); ays = c(0, 0, 0, 0, 6); ans = c(6, 3, 12, 17, 36)
#' target <- 0.4
#' decision <- CFOeff.next(target,axs,ays,ans,currdose = 3, mineff = 0.3)
#' summary(decision)
#' \donttest{#early stop for overly toxic
#' axs = c(13, 11, 7, 11, 26); ays = c(25, 18, 12, 17, 26); ans = c(36, 23, 22, 27, 36)
#' target <- 0.4
#' decision <- CFOeff.next(target,axs,ays,ans,currdose = 1, mineff = 0.3)
#' summary(decision)
#' }
#' \donttest{#early stop for low efficacy
#' axs = c(0, 0, 0, 0, 0); ays = c(2, 1, 1, 1, 6); ans = c(36, 23, 22, 27, 36)
#' target <- 0.4
#' decision <- CFOeff.next(target,axs,ays,ans,currdose = 1, mineff = 0.3)
#' summary(decision)
#' }

CFOeff.next <- function(target, axs, ays, ans, currdose, prior.para=list(alp.prior = target, bet.prior = 1 - target, alp.prior.eff = 0.5, 
                                                                         bet.prior.eff = 0.5), 
                        cutoff.eli=0.95, early.stop=0.95, effearly.stop=0.9, mineff){
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
  
  overdose.fn <- function(phi, threshold, prior.para=list()){
    y <- prior.para$y
    n <- prior.para$n
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
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
  # i.e., Pr(Y=y|lower<phi<upper) upper = 1 if upper > 1
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
  #the results for current 3 dose levels
  if (currdose != 1) {
    cys <- ays[(currdose - 1):(currdose + 1)]
    cns <- ans[(currdose - 1):(currdose + 1)]
  }else{
    cys <- c(NA, ays[1:(currdose + 1)])
    cns <- c(NA, ans[1:(currdose + 1)])
  }
  
  
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  alp.prior.eff <- prior.para$alp.prior.eff
  bet.prior.eff <- prior.para$bet.prior.eff
  
  cover.doses <- c(0,0,0)
  cunder.effs <- c(0,0,0)
  
  for (i in 1:3){
    cy <- cys[i]
    cn <- cns[i]
    if (is.na(cn)){
      cover.doses[i] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior, 
                                            alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff))
      if (overdose.fn(target, cutoff.eli, prior.para)){
        cover.doses[i:3] <- 1
        break()
      }
    }
  }
  cover.prob <- rep(0, length(ays))
  for (i in 1:length(ays)){
    ty <- ays[i]
    tn <- ans[i]
    if (is.na(tn)){
      cover.prob[i] <- NA
    }else{
      cover.prob[i] <- post.prob.fn(target, ty, tn, alp.prior, bet.prior)
    }
  }
  
  idx <- if (currdose == 2) 1 else if (currdose == 1) 2 else NA
  if (!is.na(idx) & (cutoff.eli != early.stop)) {
    cy <- cys[idx]
    cn <- cns[idx]
    if (is.na(cn)){
      cover.doses[idx] <- NA
    }else{
      prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior, 
                                            alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff))
      if (overdose.fn(target, early.stop, prior.para)){
        cover.doses[idx:3] <- 1
      }
    }
  }
  
  cover.doses <- ifelse(is.na(cys), NA, cover.doses)
  
  
  position <- which(cover.doses == 1)[1]
  overtox <- c(-1, 0, 1)[position] + currdose
  prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior, 
                       alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff))
  if ((cover.doses[2] == 1)&(currdose == 1)){
    index <- NA
    decision <- "stop_for_tox"
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
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          index <- -1
          decision <- "de-escalation"
        }else if (!v1 & v2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
    }
  }
  

  if (decision == 'stop_for_tox'){
    set <- NULL
    probs <- NULL
    nextdose <- 99
  }else{
    up.idx <- currdose + index
  if (up.idx == 1) {
    nextdose <- 1
    probs <- 1
    set <- 1
    cover.prob <- post.prob.fn(target, ays[1], ans[1], alp.prior, bet.prior)
  }else{
    low.idx <- 1
    ad.xs <- axs[low.idx:up.idx]
    ad.ys <- ays[low.idx:up.idx]
    ad.ns <- ans[low.idx:up.idx]
    set <- c(low.idx:up.idx)
    
    for (dose in set){
      ax <- ad.xs[dose]
      an <- ad.ns[dose]
      if (is.na(cn)){
        cover.doses[dose] <- NA
      }else{
        prior.para <- c(list(x=ax, n=an),list(alp.prior=alp.prior, bet.prior=bet.prior, 
                                              alp.prior.eff = alp.prior.eff, bet.prior.eff = bet.prior.eff))
        if (under.eff.fn(mineff, effearly.stop, prior.para)){
          cunder.effs[dose] <- 1
        }
      }
    }
    probs <- moveprobs(ad.xs, ad.ns, prior.para$alp.prior.eff, prior.para$bet.prior.eff)
    if (sum(cunder.effs) == length(set)){
      nextdose <- 99
      decision = 'stop_for_low_eff'
    }
    else {
      if (length(ad.xs) == 1) {
        nextdose <- low.idx
      }else{
        nextdose <- which.max(probs)
        
      }
    }
    cover.prob = cover.prob[low.idx:up.idx]

  }
  }
  if(decision == "stop_for_tox"){#This only happen when current dose level is 1
    cover.prob = cover.prob[1]
  }
  
  out <- list(target = target, axs = axs, ays = ays, ans = ans, decision = decision, currdose = currdose, 
              nextdose = nextdose, overtox = overtox, toxprob = cover.prob, effprob = probs, 
              admset = set, class = "phaseI/II")
  class(out) <- c("cfo_eff_decision", "cfo")
  return(out)
}


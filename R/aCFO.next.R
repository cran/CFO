#' Determination of the dose level for next cohort in the accumulative calibration-free odds (aCFO) design for phase I trials
#' 
#' In the aCFO design for phase I trials, the function is used to determine the dose movement based on the toxicity outcomes of the enrolled cohorts.
#'
#' @usage aCFO.next(target, ays, ans, currdose, 
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target),
#'        cutoff.eli = 0.95, early.stop = 0.95)
#'
#' @param target the target DLT rate.
#' @param ays the cumulative numbers of DLTs observed in patients for all dose levels.
#' @param ans the cumulative numbers of patients for all dose levels.
#' @param currdose the current dose level.
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' 
#'
#' @details The aCFO design is an extension of the CFO design. It integrates dose information from all positions (ranging 
#' from the lowest to the highest dose levels) into the decision-making process of the trial. Before assigning the dose level 
#' for a new cohort, aCFO compares the evidence from the current dose level with all doses to its left and right. In contrast, 
#' the original CFO design makes dose allocation by examining one dose level above and one below the current dose level. 
#' Consequently, the aCFO design enhances the utilization of information while maintaining the characteristics of the CFO 
#' design (model-free and calibration-free). Additionally, the aCFO design preserves the same early stopping and dose 
#' elimination criteria as the CFO design.
#'          
#' @note The dose level indicated by \code{overtox} and all the dose levels above experience over-toxicity, and these dose levels will be eliminated.
#'          
#' @return The \code{aCFO.next()} function returns a list object comprising the following elements: 
#' \itemize{
#'   \item target: the target DLT rate.
#'   \item ays: the cumulative counts of DLTs observed at all dose levels.
#'   \item ans: the cumulative counts of patients treated at all dose levels.
#'   \item decision: the decision in the aCFO design, where \code{left}, \code{stay}, and \code{right} represent the 
#'   movement directions, and \code{stop} indicates stopping the experiment.
#'   \item currdose: the current dose level.
#'   \item nextdose: the recommended dose level for the next cohort. \code{nextdose = 99} indicates that the trial is 
#'   terminated due to early stopping.
#'   \item overtox: the situation regarding which position experiences over-toxicity. The dose level indicated by 
#'   \code{overtox} and all the dose levels above experience over-toxicity. \code{overtox = NA} signifies that the 
#'   occurrence of over-toxicity did not happen.
#'   \item toxprob: the expected toxicity probability, \eqn{Pr(p_k > \phi | x_k, m_k)}, at all dose
#'   levels, where \eqn{p_k}, \eqn{x_k}, and \eqn{m_k} is the dose-limiting toxicity (DLT) rate, the 
#'   numbers of observed DLTs, and the numbers of patients at dose level \eqn{k}. \code{NA} indicates that there 
#'   are no patients at the corresponding dose level.
#' }
#'  
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin 
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Fang J, Yin G (2024). Fractional accumulative calibration‐free odds (f‐aCFO) design for delayed toxicity 
#'             in phase I clinical trials. \emph{Statistics in Medicine}, 43(17), 3210-3226.
#'
#' @examples
#' ## determine the dose level for the next cohort of new patients
#' ays <- c(0, 0, 1, 0, 0, 0, 0); ans <- c(3, 3, 6, 0, 0, 0, 0)
#' decision <- aCFO.next(target = 0.2, ays = ays, ans = ans, currdose = 3, 
#'             prior.para = list(alp.prior = 0.2, bet.prior = 0.8))
#' summary(decision)
#' 
#' ays <- c(3, 0, 0, 0, 0, 0, 0); ans <- c(3, 0, 0, 0, 0, 0, 0)
#' decision <- aCFO.next(target = 0.2, ays = ays, ans = ans, currdose = 1,
#'             prior.para = list(alp.prior = 0.2, bet.prior = 0.8))
#' summary(decision)
#' 
#' ays <- c(0, 0, 0, 0, 0, 0, 3); ans <- c(3, 3, 3, 3, 3, 3, 3)
#' decision <- aCFO.next(target = 0.2, ays = ays, ans = ans, currdose = 7,
#'             prior.para = list(alp.prior = 0.2, bet.prior = 0.8))
#' summary(decision)
#' 
#' @import stats
#' @export
aCFO.next <- function(target, ays, ans, currdose, prior.para=list(alp.prior=target, bet.prior=1-target),
                      cutoff.eli=0.95, early.stop=0.95){
  
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
  
  OR.union.values <- function(phi, cans, cays, alp.prior, bet.prior, type){
    ndose <- length(cays)
    if (type=="L"){
      OR.list <- rep(0, ndose-1)
      for (i in 1:(ndose-1)){
        OR.list[i] <- OR.values(phi, cays[i], cans[i], cays[ndose], cans[ndose], alp.prior, bet.prior, type)
      }
    }else if (type=="R"){
      OR.list <- rep(0, ndose-1)
      for (i in 2:ndose){
        OR.list[i-1] <- OR.values(phi, cays[1], cans[1], cays[i], cans[i], alp.prior, bet.prior, type)
      }
    }
    return(sum(OR.list))
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
  # i.e., Pr(Y=y|lower<phi<upper); upper = 1 if upper > 1
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
  
  optim.gamma.union.fn <- function(cans, phi, type, alp.prior, bet.prior){
    ndose <- length(cans)
    if (type == "L"){
      gamma.list <- rep(0, ndose-1)
      for (i in 1:(ndose-1)){
        gamma.list[i] <- optim.gamma.fn(cans[i], cans[ndose], phi, type, alp.prior, bet.prior)$gamma
      }
    }else if (type == "R"){
      gamma.list <- rep(0, ndose-1)
      for (i in 2:ndose){
        gamma.list[i-1] <- optim.gamma.fn(cans[1], cans[i], phi, type, alp.prior, bet.prior)$gamma
      }
    }
    return(sum(gamma.list))
  }
  
  ###############################################################################
  ############################MAIN DUNCTION######################################
  ###############################################################################
  ndose <- length(ays)
  
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  tover.doses <- rep(0, ndose)
  for (i in 1:ndose){
    cy <- ays[i]
    cn <- ans[i]
    prior.para <- c(list(y=cy, n=cn), list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (overdose.fn(target, cutoff.eli, prior.para)){
      tover.doses[i:ndose] <- 1
      break()
    }
  }
  
  tover.prob <- rep(0, ndose)
  for (i in 1:ndose){
    cy <- ays[i]
    cn <- ans[i]
    tover.prob[i] <- post.prob.fn(target, cy, cn, alp.prior, bet.prior)
  }
  
  if (cutoff.eli != early.stop) {
    cy <- ays[1]
    cn <- ans[1]
    prior.para <- c(list(y=cy, n=cn),list(alp.prior=alp.prior, bet.prior=bet.prior))
    if (overdose.fn(target, early.stop, prior.para)){
      tover.doses[1:ndose] <- 1
    }
  }
  
  position <- which(tover.doses == 1)[1]
  prior.para <- c(list(alp.prior=alp.prior, bet.prior=bet.prior))
  if ((tover.doses[1] == 1) & (position == 1)){
    index <- NA
    decision <- "stop"
  } else {
    if (currdose!=1){
      cys <- ays[(currdose-1):(currdose+1)]
      cns <- ans[(currdose-1):(currdose+1)]
      cover.doses <- tover.doses[(currdose-1):(currdose+1)]
      #cover.doses <- c(0, 0, 0) # No elimination rule
    }else{
      cys <- c(NA, ays[1:(currdose+1)])
      cns <- c(NA, ans[1:(currdose+1)])
      cover.doses <- c(NA, tover.doses[1:(currdose+1)])
      #cover.doses <- c(NA, 0, 0) # No elimination rule
    }
    
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
        OR.v2 <- OR.union.values(target, ans[currdose:ndose], ays[currdose:ndose], alp.prior, bet.prior, type="R")
        gam2 <- optim.gamma.union.fn(ans[currdose:ndose], target, "R", alp.prior, bet.prior)
        if (OR.v2>gam2){
          index <- 1
          decision <- "escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- optim.gamma.union.fn(ans[1:currdose], target, "L", alp.prior, bet.prior)
        OR.v1 <- OR.union.values(target, ans[1:currdose], ays[1:currdose], alp.prior, bet.prior, type="L")
        if (OR.v1>gam1){
          index <- -1
          decision <- "de-escalation"
        }else{
          index <- 0
          decision <- "stay"
        }
      }
      else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- optim.gamma.union.fn(ans[1:currdose], target, "L", alp.prior, bet.prior)
        gam2 <- optim.gamma.union.fn(ans[currdose:ndose], target, "R", alp.prior, bet.prior)
        OR.v1 <- OR.union.values(target, ans[1:currdose], ays[1:currdose], alp.prior, bet.prior, type="L")
        OR.v2 <- OR.union.values(target, ans[currdose:ndose], ays[currdose:ndose], alp.prior, bet.prior, type="R")
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
  
  if (decision=='stop'){
    nextdose <- 99
  }else{
    nextdose <- currdose+index
  }
  
  out <- list(target=target, ays=ays, ans=ans, decision=decision, currdose = currdose, 
              nextdose=nextdose, overtox=position, toxprob=tover.prob)
  class(out) <- c("cfo_decision", "cfo")
  return(out)
}


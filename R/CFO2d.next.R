#' Determinate the dose level for the next cohort in the two-dimensional calibration-free odds (2dCFO) design.
#'
#' This function is used to determine the next dose level for the next cohort in the 2dCFO design.
#'
#' @usage CFO2d.next(target, cys, cns, currdose, 
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target), 
#'        cutoff.eli = 0.95, early.stop = 0.95, seed = NULL)
#'
#' @param target the target DLT rate.
#' @param cys a matrix of the number of DLTs observed for each dose combination.
#' @param cns a matrix of the number of patients allocated to each dose combination.
#' @param currdose a vector of the current dose indices in the horizontal and vertical direction.
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' @param seed an integer to be set as the seed of the random number generator for reproducible results. The default value is set to \code{NULL}.
#' 
#' @details In the 2dCFO design, decision-making within the two-dimensional toxicity probability space is conducted by performing two independent one-dimensional 
#'          CFO analyses along both the horizontal and vertical axes (Wang et al. 2023).
#'  
#' @note    When the current dose level is the lowest or highest (i.e., at the boundary), the parts in \code{cys} and 
#'          \code{cns} where there is no data are filled with \code{NA}. \cr
#'          The dose level indicated by \code{overtox} and all the dose levels above experience over-toxicity, and these dose levels will be eliminated.
#' 
#' @return The \code{CFO2d.next()} function returns a list with the following components:
#' \itemize{
#'   \item target: the target DLT rate.
#'   \item cys: a 3 by 3 matrix of the number of DLT observed for each dose combination at and around the current dose.
#'   \item cns: a 3 by 3 matrix of the number of patients allocated to each dose combination at and around the current dose.
#'   \item decision: a vector of length 2 representing the recommended decisions for vertical and horizontal 
#'   directions, and \code{stop} indicates stopping the experiment.
#'   \item currdose: the current dose combination.
#'   \item nextdose: the recommended dose combination for the next cohort. \code{nextdose = (99, 99)} indicates that the trial is 
#'   terminated due to early stopping.
#'   \item overtox: the situation regarding which positions experience over-toxicity. The dose level indicated 
#'   by \code{overtox} and all the dose levels above experience over-toxicity. \code{overtox = NA} signifies that the 
#'   occurrence of over-toxicity did not happen.
#' }
#' 
#' @author Jialu Fang, Wenliang Wang, and Guosheng Yin
#' 
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Wang W, Jin H, Zhang Y, Yin G (2023). Two-dimensional calibration-free odds (2dCFO)
#'             design for phase I drug-combination trials. \emph{Frontiers in Oncology}, 13, 1294258.
#' 
#' @export
#' @examples
#' cns <- matrix(c(3, 3, 0,
#'                 0, 6, 0,
#'                 0, 0, 0), 
#'               nrow = 3, ncol = 3, byrow = TRUE)
#' 
#' cys <- matrix(c(0, 1, 0,
#'                 0, 2, 0,
#'                 0, 0, 0), 
#'               nrow = 3, ncol = 3, byrow = TRUE)
#' currdose <- c(2,3)
#' decision <- CFO2d.next(target = 0.3, cys, cns, currdose = currdose, seed = 1)
#' summary(decision)


CFO2d.next <- function(target, cys, cns, currdose, prior.para=list(alp.prior=target, bet.prior=1-target), cutoff.eli=0.95, early.stop=0.95, seed=NULL){
  cidx.A <- 0
  cidx.B <- 0
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  set.seed(seed)
  cover.doses=matrix(0,3,3)
  overtox <- NA
  
  # posterior probability of pj >= phi given data
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.9){
    alp <- alp.prior + y 
    bet <- bet.prior + n - y
    1 - pbeta(phi, alp, bet)
  }
  
  prob.int <- function(phi, y1, n1, y2, n2, alp.prior, bet.prior){
    alp1 <- alp.prior + y1
    alp2 <- alp.prior + y2
    bet1 <- alp.prior + n1 - y1
    bet2 <- alp.prior + n2 - y2
    fn.min <- function(x){
      dbeta(x, alp1, bet1)*(1-pbeta(x, alp2, bet2))
    }
    fn.max <- function(x){
      pbeta(x, alp1, bet1)*dbeta(x, alp2, bet2)
    }
    const.min <- integrate(fn.min, lower=0, upper=1)$value
    const.max <- integrate(fn.max, lower=0, upper=1)$value
    p1 <- integrate(fn.min, lower=0, upper=phi)$value/const.min
    p2 <- integrate(fn.max, lower=0, upper=phi)$value/const.max
    
    list(p1=p1, p2=p2)
  }
  
  overdose.fn <- function(phi, threshold, y, n, prior.para=list(alp.prior=phi, bet.prior=1-phi)){
    alp.prior <- prior.para$alp.prior
    bet.prior <- prior.para$bet.prior
    pp <- post.prob.fn(phi, y, n, alp.prior, bet.prior)
    if ((pp >= 0.95) & (n>=3)){
      return(TRUE)
    }else{
      return(FALSE)
    }
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
    }else if (type=="D"){
      pC <- 1 - ps$p2
      pD <- 1 - ps$p1
      oddsC <- pC/(1-pC)
      oddsD <- pD/(1-pD)
      OR <- oddsC*oddsD
    }else if (type=="U"){
      pC <- 1 - ps$p1
      pU <- 1 - ps$p2
      oddsC <- pC/(1-pC)
      oddsU <- pU/(1-pU)
      OR <- (1/oddsC)/oddsU
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
  
  make.decision.1dCFO.fn <- function(phi, cys, cns, alp.prior, bet.prior, cover.doses, diag=FALSE){
    if (cover.doses[2] == 1){
      return(1)
    }else{
      if (is.na(cys[1]) & (cover.doses[3]==1)){
        return(2)
      }else  if (is.na(cys[1]) & (!(cover.doses[3]==1))){
        gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma
        OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        if (OR.v2>gam2){
          return(3)
        }else{
          return(2)
        }
      }else  if (is.na(cys[3]) | (cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        if (OR.v1>gam1){
          return(1)
        }else{
          return(2)
        }
        
      }else  if (!(is.na(cys[1]) | is.na(cys[3]) | cover.doses[3]==1)){
        gam1 <- optim.gamma.fn(cns[1], cns[2], phi, "L", alp.prior, bet.prior)$gamma
        gam2 <- optim.gamma.fn(cns[2], cns[3], phi, "R", alp.prior, bet.prior)$gamma
        OR.v1 <- OR.values(phi, cys[1], cns[1], cys[2], cns[2], alp.prior, bet.prior, type="L")
        OR.v2 <- OR.values(phi, cys[2], cns[2], cys[3], cns[3], alp.prior, bet.prior, type="R")
        v1 <- OR.v1 > gam1
        v2 <- OR.v2 > gam2
        if (v1 & !v2){
          return(1)
        }else if (!v1 & v2){
          return(3)
        }else{
          return(2)
        }
      }
    }
  }
  
  if (overdose.fn(target, cutoff.eli, cys[2,2], cns[2,2], prior.para)){
    cover.doses[2,2] <- 1
    overtox <- currdose
  }
  if (!is.na(cns[2,3])){
    if (overdose.fn(target, cutoff.eli, cys[2,3], cns[2,3], prior.para)){
      cover.doses[2,3] <- 1
      cover.doses[3,3] <- 1
      overtox <- currdose + c(0,1)
    }
  } else {
    cover.doses[2,3] <- NA
    cover.doses[3,3] <- NA
  }
  if (!is.na(cns[3,2])){
    if (overdose.fn(target, cutoff.eli, cys[3,2], cns[3,2], prior.para)){
      cover.doses[3,2] <- 1
      cover.doses[3,3] <- 1
      overtox <- currdose + c(1,0)
    }
  } else {
    cover.doses[3,2] <- NA
    cover.doses[3,3] <- NA
  }
  if (!is.na(cns[2,3])&!is.na(cns[3,2])){
    if(overdose.fn(target, cutoff.eli, cys[2,3], cns[2,3], prior.para)&overdose.fn(target, cutoff.eli, cys[3,2], cns[3,2], prior.para)){
      overtox <- currdose
    }
  }
  
  if (cutoff.eli != early.stop) {
    if (currdose==c(1,1) & overdose.fn(target, early.stop, cys[1,1], cns[1,1], prior.para)){
      cover.doses[1,1] <- 1
      out <- list(target=target, cys=cys, cns=cns, decision="stop", currdose = currdose, nextdose = c(99,99), overtox = c(1,1))
      class(out) <- "cfo"
      return(out)
    }
  }
  
  
  # horizontal direction
  idx.chg.A <- make.decision.1dCFO.fn(target, cys[2,], cns[2,], alp.prior, bet.prior, cover.doses[2,]) - 2
  # vertical direction
  idx.chg.B <- make.decision.1dCFO.fn(target, cys[,2], cns[,2], alp.prior, bet.prior, cover.doses[,2]) - 2
  
  if (idx.chg.A == 1 & idx.chg.B == 1){
    ### horizontal and vertical only
    
    OR.R <- OR.values(target, cys[2,2], cns[2,2], cys[2,3], cns[2,3], alp.prior, bet.prior, type="R")
    OR.U <- OR.values(target, cys[2,2], cns[2,2], cys[3,2], cns[3,2], alp.prior, bet.prior, type="R")
    
    if (OR.R == OR.U){
      rand <- rbinom(1,1,0.5)
      if(rand == 0){
        cidx.A <- 1
      } else {
        cidx.B <- 1
      }
    } else if (OR.R > OR.U){
      cidx.B <- 1
    } else {
      cidx.A <- 1
    }
    
  } else if (idx.chg.A == -1 & idx.chg.B == -1){
    if (is.na(cys[2,1]) & is.na(cys[1,2])){
      cidx.A <- 0
      cidx.B <- 0
    } else if (is.na(cys[2,1])){
      cidx.A <- -1
    } else if (is.na(cys[1,2])){
      cidx.B <- -1
    } else {
      OR.L <- OR.values(target, cys[2,2], cns[2,2], cys[2,1], cns[2,1], alp.prior, bet.prior, type="L")
      OR.D <- OR.values(target, cys[2,2], cns[2,2], cys[1,2], cns[1,2], alp.prior, bet.prior, type="L")
      
      if (OR.L == OR.D){
        rand <- rbinom(1,1,0.5)
        if(rand == 0){
          cidx.A <- -1
        } else {
          cidx.B <- -1
        }
      } else if (OR.L > OR.D){
        cidx.B <- -1
      } else {
        cidx.A <- -1
      }
    }
  } else if (idx.chg.A == 1 & idx.chg.B == -1){
    DCR <- make.decision.1dCFO.fn(target, c(cys[1,2],cys[2,2],cys[2,3]), c(cns[1,2],cns[2,2],cns[2,3]), alp.prior, 
                                  bet.prior, c(cover.doses[1,2],cover.doses[2,2],cover.doses[2,3])) - 2
    if (DCR == 1){
      cidx.B <- 1
    } else if (DCR == -1){
      cidx.A <- -1
    }
  } else if (idx.chg.A == -1 & idx.chg.B == 1){
    LCU <- make.decision.1dCFO.fn(target, c(cys[2,1],cys[2,2],cys[3,2]), c(cns[2,1],cns[2,2],cns[3,2]), alp.prior, 
                                  bet.prior, c(cover.doses[2,1],cover.doses[2,2],cover.doses[3,2])) - 2
    if (LCU == 1){
      cidx.A <- 1
    } else if (LCU == -1){
      cidx.B <- -1
    }
  } else if (idx.chg.A == 1 & idx.chg.B == 0){
    cidx.B <- 1
  } else if (idx.chg.A == 0 & idx.chg.B == 1){
    cidx.A <- 1
  } else if (idx.chg.A == -1 & idx.chg.B == 0){
    cidx.B <- -1
  } else if (idx.chg.A == 0 & idx.chg.B == -1){
    cidx.A <- -1
  }
  
  nextdose <- currdose+c(cidx.A, cidx.B)
  decision_values <- c("de-escalation", "stay", "escalation")
  decision <- decision_values[match(c(cidx.A, cidx.B), c(-1, 0, 1))]
  out <- list(target=target, cys=cys, cns=cns, decision=decision, currdose = currdose, nextdose = nextdose, overtox = overtox)
  class(out) <- "cfo"
  return(out)
}

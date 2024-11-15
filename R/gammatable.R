#'Generating table of threshold \eqn{\gamma_L} and \eqn{\gamma_R} in the calibration-free odds (CFO) design
#'
#'Generate all the possible thresholds under different \eqn{m_C}, \eqn{m_L} and \eqn{m_R}
#' @usage gammatable(npatient, target, 
#'             para.prior = list(alp.prior = target, bet.prior = 1 - target))
#' @param npatient the numbers of patients involved in the trial.
#' @param target the target DLT rate.
#' @param para.prior the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}). 
#'
#' @return The \code{gammatable()} function returns a list object comprising the following elements:
#' \itemize{
#'    \item gammatb.left: the table of threshold \eqn{\gamma_L} under different \eqn{m_L} 
#'    and \eqn{m_C} where \eqn{m_C} and \eqn{m_L} represent the number of patients at current dose level and left dose level.
#'    \item gammatb.right: the table of threshold \eqn{\gamma_R} under different \eqn{m_R} 
#'    and \eqn{m_C} where \eqn{m_C} and \eqn{m_R} represent the number of patients at current dose level and right dose level.
#' }
#' @note This function generate two matrices. \code{gammatb.left} contains the threshold \eqn{\gamma_L}, 
#' and \code{gammatb.right} contains the threhold \eqn{\gamma_R}. For matrix \code{gammatb.left}, the row index represent the number of patients 
#' at left dose level, and the column index represent the number of patients at current dose level. For matrix \code{gammatb.right}, the row index represent the number of patients 
#' at right dose level, and the column index represent the number of patients at current dose level. 
#' For example, if you want to get the threshold \eqn{\gamma_L} in the case of \eqn{m_C = 12, m_L = 13}, you can reach it by \code{result$gammatb.left[13,12]} 
#' 
#' @import utils 
#' @author Jialu Fang, Ninghao Zhang, Wenliang Wang, and Guosheng Yin
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066.
#' @examples
#' npatient <- 3; target <- 0.3
#' para.prior = list(alp.prior = target, bet.prior = 1 - target)
#' result <- gammatable(npatient, target, para.prior)
#' plot(result)
#' \donttest{#This example may cost you a long time to run
#' npatient <- 30; target <- 0.3
#' para.prior = list(alp.prior = target, bet.prior = 1 - target)
#' result <- gammatable(npatient, target, para.prior)
#' plot(result)
#' }
#' @export
gammatable <- function(npatient, target, para.prior = list(alp.prior = target, bet.prior = 1 - target)){
  ###############################################################################
  ###############define the functions used for main function#####################
  ###############################################################################
  post.prob.fn <- function(phi, y, n, alp.prior=0.1, bet.prior=0.1){
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
    const.min <- integrate(fn.min, lower=0, upper=0.99, subdivisions=1000, rel.tol = 1e-10)$value
    const.max <- integrate(fn.max, lower=0, upper=1, rel.tol = 1e-10)$value
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
  gamtableL <- matrix(nrow = npatient, ncol = npatient)
  minerrtableL <- matrix(nrow = npatient, ncol = npatient)
  gamtableR <- matrix(nrow = npatient, ncol = npatient)
  minerrtableR <- matrix(nrow = npatient, ncol = npatient)
  pb <- txtProgressBar(style = 3)
  nrun <- 0 
  for (i in 1:npatient){
    for (j in 1:npatient){
      resL <- optim.gamma.fn(i, j, target, "L",  para.prior$alp.prior, para.prior$bet.prior)
      resR <- optim.gamma.fn(j, i, target, "R",  para.prior$alp.prior, para.prior$bet.prior)
      gamtableL[i,j] <- resL$gamma
      gamtableR[i,j] <- resR$gamma
      setTxtProgressBar(pb,((i - 1)*npatient +j)/(npatient*npatient))
    }
  }
  close(pb)
  out <- list(gammatb.left = gamtableL, gammatb.right = gamtableR)
  class(out) <- c("cfo_decision", "cfo")
  return(out)
}  

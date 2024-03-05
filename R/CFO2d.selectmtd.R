#' Select the maximum tolerated dose (MTD) for the real drug combination trial
#'
#' Select the maximum tolerated dose (MTD) when the real drug combination trial is completed
#'
#' @usage CFO2d.selectmtd(target, npts, ntox, 
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target), 
#'        cutoff.eli = 0.95, early.stop = 0.95, verbose = TRUE)
#'
#' @param target the target DLT rate.
#' @param npts a matrix containing the number of patients treated at each dose level.
#' @param ntox a matrix containing the number of patients who experienced DLT at each dose level.
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' @param verbose set \code{verbose = TRUE} to return more details of the results.
#'
#' @details \code{CFO2d.selectmtd()} selects the MTD based on isotonic estimates of toxicity
#'          probabilities. \code{CFO2d.selectmtd()} selects as the MTD dose \eqn{j^*}, for which the
#'          isotonic estimate of the DLT rate is closest to the target. If there
#'          are ties, we select from the ties the highest dose level when the estimate
#'          of the DLT rate is smaller than the target, or the lowest dose level
#'          when the estimate of the DLT rate is greater than the target. The
#'          isotonic estimates are obtained by the pooled-adjacent-violators algorithm
#'          (PAVA).
#'          
#' @note  The MTD selection and dose escalation/deescalation rule are two independent
#'        components of the trial design. Isotonic regression is employed to select the MTD after the completion of the trial.
#'        When appropriate, another dose selection procedure (e.g., based on a fitted logistic model) can be used to select
#'        the MTD after the completion of the trial using the 2dCFO design.          
#'
#' @return  \code{CFO2d.selectmtd()} returns 
#' \itemize{
#'   \item target: the target DLT rate.
#'   \item MTD: the selected MTD. \code{MTD = (99, 99)} indicates that all tested doses are overly toxic.
#'   \item p_est: the isotonic estimate of the DLT probablity at each dose and associated \eqn{95\%} credible interval.
#'    \code{p_est = NA} if all tested doses are overly toxic.
#'   \item p_est_CI: the credible interval for the isotonic estimate.
#'    \code{p_est_CI = NA} if all tested doses are overly toxic.
#' }
#'
#'
#' @author Jialu Fang, Wenliang Wang, and Guosheng Yin
#'
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Wang W, Jin H, Zhang Y, Yin G (2023). Two-dimensional calibration-free odds (2dCFO)
#'             design for phase I drug-combination trials. \emph{Frontiers in Oncology}, 13, 1294258. \cr
#'             Bril G, Dykstra R, Pillers C, Robertson T (1984). Algorithm AS 206: Isotonic regression in two independent variables. 
#'             \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, 33(3), 352–357.
#'
#' @examples
#' ntox <- matrix(c(0, 0, 2, 0, 0,
#'                 0, 2, 7, 0, 0,
#'                 0, 2, 0, 0, 0), 
#'               nrow = 3, ncol = 5, byrow = TRUE)
#' 
#' npts <- matrix(c(3,  0, 12, 0, 0,
#'                 3, 12, 24, 0, 0,
#'                 3,  3,  0, 0, 0), 
#'               nrow = 3, ncol = 5, byrow = TRUE)
#' selmtd <- CFO2d.selectmtd(target=0.3, npts=npts, ntox=ntox)
#' summary(selmtd)
#' plot(selmtd)
#'
#' @import Iso 
#' @export

CFO2d.selectmtd <- function (target, npts, ntox, prior.para=list(alp.prior=target, bet.prior=1-target), cutoff.eli = 0.95,
                             early.stop = 0.95, verbose = TRUE)
{
  y = ntox
  n = npts
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  if (nrow(n) > ncol(n) | nrow(y) > ncol(y)) {
    stop("npts and ntox should be arranged in a way (i.e., rotated) such that for each of them, the number of rows is less than or equal to the number of columns.")
    
  }
  elimi = matrix(0, dim(n)[1], dim(n)[2])
  if (cutoff.eli != early.stop) {
    if (n[1, 1] >= 3) {
      if (1 - pbeta(target, y[1,1] + alp.prior, n[1,1] - y[1,1] + bet.prior) > early.stop) {
        elimi[, ] = 1
      }
    }
  }
  for (i in 1:dim(n)[1]) {
    for (j in 1:dim(n)[2]) {
      if (n[i, j] >= 3) {
        if (1 - pbeta(target, y[i,j] + alp.prior, n[i,j] - y[i,j] + bet.prior) > cutoff.eli) {
          elimi[i:dim(n)[1], j] = 1
          elimi[i, j:dim(n)[2]] = 1
          break
        }
      }
    }
  }
  
  selectdose=NULL
  
  if (elimi[1] == 1) {
    selectdose = c(99, 99)
    selectdoses = matrix(selectdose, nrow = 1)
  }else {
    phat = (y + alp.prior)/(n + alp.prior + bet.prior)
    phat = round(Iso::biviso(phat, n + alp.prior + bet.prior, warn = TRUE)[, ],2)
    # phat.out = phat
    lower.mat=qbeta(0.025,y+alp.prior,n-y+bet.prior)
    lower.mat=round(Iso::biviso(lower.mat),2)
    
    upper.mat=qbeta(0.975,y+alp.prior,n-y+bet.prior)
    upper.mat=round(Iso::biviso(upper.mat),2)
    phat.out<-matrix(paste0(format(phat,digits=1),"(",lower.mat,", ",upper.mat,")"),byrow=FALSE,nrow=dim(phat)[1])
    colnames(phat.out)=paste0("B",1:dim(n)[2])
    rownames(phat.out)=paste0("A",1:dim(n)[1])
    phat.out.noCI=round(phat,2)
    phat.out[n == 0] = "NA"
    phat[elimi == 1] = 1.1
    phat = phat * (n != 0) + (1e-05) * (matrix(rep(1:dim(n)[1],
                                                   each = dim(n)[2], len = length(n)), dim(n)[1], byrow = T) +
                                          matrix(rep(1:dim(n)[2], each = dim(n)[1], len = length(n)),
                                                 dim(n)[1]))
    
    if(is.null(selectdose)){
      phat[n == 0] = 10
      selectdose = which(abs(phat - target) == min(abs(phat -
                                                         target)), arr.ind = TRUE)
      
      
      if (length(selectdose) > 2)
        selectdose = selectdose[1, ]
      aa = function(x) as.numeric(as.character(x))

      selectdoses = matrix(99, nrow = 1, ncol = 2)
      selectdoses[1, ] = matrix(selectdose, nrow = 1)

      selectdoses = matrix(selectdoses[selectdoses[, 2] !=
                                         99, ], ncol = 2)
    }
    
    colnames(selectdoses) = c("DoseA", "DoseB")
    
  }

  if (verbose == TRUE) {
    if (selectdoses[1, 1] == 99 && selectdoses[1, 2] == 99) {
      message("All tested doses are overly toxic. No MTD is selected! \n")
      out=list(target = target, MTD = selectdoses, p_est = NA, p_est_CI = NA)
    }else{
      out=list(target = target, MTD = selectdoses, p_est=phat.out.noCI, p_est_CI = phat.out)
    }
  }
  else {
    if (selectdoses[1, 1] == 99 && selectdoses[1, 2] == 99) {
      message("All tested doses are overly toxic. No MTD is selected! \n")
    }
    out = list(target = target, MTD = selectdoses)
  }
  class(out)<-c("cfo_sel","cfo")
  return(out)
}


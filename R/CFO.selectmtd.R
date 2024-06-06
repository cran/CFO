#'
#' Select the maximum tolerated dose (MTD) for the real single-drug trial
#'
#' Select the maximum tolerated dose (MTD) when the real single-drug trial is completed
#'
#' @usage CFO.selectmtd(target, npts, ntox, 
#'        prior.para = list(alp.prior = target, bet.prior = 1 - target), 
#'        cutoff.eli = 0.95, early.stop = 0.95, verbose = TRUE)
#'
#' @param target the target DLT rate.
#' @param npts a vector containing the number of patients treated at each dose level.
#' @param ntox a vector containing the number of patients who experienced DLT at each dose level.
#' @param prior.para the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)} 
#'                  by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for 
#'                  the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).
#' @param cutoff.eli the cutoff to eliminate overly toxic doses for safety. We recommend
#'                    the default value of \code{cutoff.eli = 0.95} for general use.
#' @param early.stop the threshold value for early stopping. The default value \code{early.stop = 0.95}
#'                generally works well.
#' @param verbose set \code{verbose=TRUE} to return more details of the results.
#'
#' @details \code{CFO.selectmtd()} selects the MTD based on isotonic estimates of toxicity
#'          probabilities. \code{CFO.selectmtd()} selects as the MTD dose \eqn{j^*}, for which the
#'          isotonic estimate of the DLT rate is closest to the target. If there
#'          are ties, we select from the ties the highest dose level when the estimate
#'          of the DLT rate is smaller than the target, or the lowest dose level
#'          when the estimate of the DLT rate is greater than the target. The
#'          isotonic estimates are obtained by the pooled-adjacent-violators algorithm
#'          (PAVA).
#'          
#'
#' @return  \code{CFO.selectmtd()} returns 
#' \itemize{
#'   \item target: the target DLT rate.
#'   \item MTD: the selected MTD. \code{MTD = 99} indicates that all tested doses are overly toxic.
#'   \item p_est: the isotonic estimate of the DLT probablity at each dose and associated \eqn{95\%} credible interval.
#'    \code{p_est = NA} if all tested doses are overly toxic.
#'   \item p_overdose: the probability of overdosing defined as \eqn{Pr(toxicity > \code{target}|data)}.
#'    \code{p_overdose = NA} if all tested doses are overly toxic.
#' }
#'
#'
#' @note  The MTD selection and dose escalation/de-escalation rule are two independent
#'        components of the trial design. Isotonic regression is employed to select the MTD after the completion of the trial.
#'        When appropriate, another dose selection procedure (e.g., based on a fitted logistic model) can be used to select
#'        the MTD after the completion of the trial using the CFO-type design.
#'
#'
#' @author Jialu Fang, Wenliang Wang, and Guosheng Yin
#'
#' @references Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
#'             \emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
#'             Bril G, Dykstra R, Pillers C, Robertson T (1984). Algorithm AS 206: Isotonic regression in two independent variables. 
#'             \emph{Journal of the Royal Statistical Society. Series C (Applied Statistics)}, 33(3), 352–357. \cr
#'             Fang J, Yin G (2024). Fractional accumulative calibration‐free odds (f‐aCFO) design for delayed toxicity 
#'             in phase I clinical trials. \emph{Statistics in Medicine}.
#'
#'
#' @examples
#'
#' ### select the MTD for the CFO-type single-drug trial
#' n <- c(3,3,27,3,0,0,0)
#' y <- c(0,0,4,2,0,0,0)
#' selmtd <- CFO.selectmtd(target=0.2, npts=n, ntox=y)
#' summary(selmtd)
#' plot(selmtd)
#'
#' @export
CFO.selectmtd <- function (target, npts, ntox, prior.para=list(alp.prior=target, bet.prior=1-target), 
                        cutoff.eli = 0.95, early.stop = 0.95, verbose = TRUE)
{
  pava <- function(x, wt = rep(1, length(x))) {
    n <- length(x)
    if (n <= 1)
      return(x)
    if (any(is.na(x)) || any(is.na(wt))) {
      stop("Missing values in 'x' or 'wt' not allowed")
    }
    lvlsets <- (1:n)
    repeat {
      viol <- (as.vector(diff(x)) < 0)
      if (!(any(viol)))
        break
      i <- min((1:(n - 1))[viol])
      lvl1 <- lvlsets[i]
      lvl2 <- lvlsets[i + 1]
      ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)
      x[ilvl] <- sum(x[ilvl] * wt[ilvl])/sum(wt[ilvl])
      lvlsets[ilvl] <- lvl1
    }
    x
  }
  
  if (is.null(prior.para$alp.prior)){
    prior.para <- c(prior.para, list(alp.prior=target, bet.prior=1-target))
  }
  alp.prior <- prior.para$alp.prior
  bet.prior <- prior.para$bet.prior
  
  y = ntox
  n = npts
  ndose = length(n)
  elimi = rep(0, ndose)
  for (i in 1:ndose) {
    if (n[i] >= 3) {
      if (1 - pbeta(target, y[i] + alp.prior, n[i] - y[i] + bet.prior) >
          cutoff.eli) {
        elimi[i:ndose] = 1
        break
      }
    }
  }
  if (cutoff.eli != early.stop) {
    if (n[1] >= 3) {
      if (1 - pbeta(target, y[1] + alp.prior, n[1] - y[1] + bet.prior) >
          early.stop) {
        elimi[1:ndose] = 1
      }
    }
  }
  if (elimi[1] == 1 || sum(n[elimi == 0]) == 0){
    selectdose = 99
  }else {
    adm.set = (n != 0) & (elimi == 0)
    adm.index = which(adm.set == T)
    y.adm = y[adm.set]
    n.adm = n[adm.set]
    phat = (y.adm + alp.prior)/(n.adm + alp.prior + bet.prior)
    phat.var = (y.adm + alp.prior) * (n.adm - y.adm + bet.prior)/((n.adm +
                          alp.prior + bet.prior)^2 * (n.adm + alp.prior + bet.prior + 1))
    phat = pava(phat, wt = 1/phat.var)
    phat = phat + (1:length(phat)) * 1e-10
    selectd = sort(abs(phat - target), index.return = T)$ix[1]
    selectdose = adm.index[selectd]
  }
  
  if (verbose == TRUE) {
    trtd = (n != 0)
    poverdose = pava(1 - pbeta(target, y[trtd] + alp.prior, n[trtd] -
                                 y[trtd] + bet.prior))
    phat.all = pava((y[trtd] + alp.prior)/(n[trtd] + alp.prior + bet.prior), wt = 1/((y[trtd] +
                          alp.prior) * (n[trtd] - y[trtd] + bet.prior)/((n[trtd] + alp.prior + bet.prior)^2 *
                                 (n[trtd] + alp.prior + bet.prior + 1))))
    A1 = A2 = A3 = A4 = NULL
    k = 1
    for (i in 1:ndose) {
      if (n[i] > 0) {
        A1 = append(A1, formatC(phat.all[k], digits = 2,
                                format = "f"))
        A2 = append(A2, formatC(qbeta(0.025, y[i] + alp.prior,
                                      n[i] - y[i] + bet.prior), digits = 2, format = "f"))
        A3 = append(A3, formatC(qbeta(0.975, y[i] + alp.prior,
                                      n[i] - y[i] + bet.prior), digits = 2, format = "f"))
        A4 = append(A4, formatC(poverdose[k], digits = 2,
                                format = "f"))
        k = k + 1
      }
      else {
        A1 = append(A1, "----")
        A2 = append(A2, "----")
        A3 = append(A3, "----")
        A4 = append(A4, "----")
      }
    }
    p_est = data.frame(cbind(dose = 1:length(npts), phat = A1,
                             CI = paste("(", A2, ",", A3, ")", sep = "")))
    if (selectdose == 99) {
      message("All tested doses are overly toxic. No MTD is selected! \n")
      out = list(target = target, MTD = selectdose, p_est = NA, p_overdose = NA)
    } else {
      out = list(target = target, MTD = selectdose, p_est = p_est, p_overdose = A4)
    }
  }
  else {
    if (selectdose == 99) {
      message("All tested doses are overly toxic. No MTD is selected! \n")
    }
    out = list(target = target, MTD = selectdose)
  }
  class(out)<-c("cfo_sel","cfo")
  return(out)
}


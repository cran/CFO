% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CFO.simu.R
\name{CFO.simu}
\alias{CFO.simu}
\title{Conduct one simulation using the calibration-free odds (CFO) or accumulative CFO (aCFO) design for phase I trials.}
\usage{
CFO.simu(design, target, p.true, init.level = 1, ncohort, cohortsize,
       prior.para = list(alp.prior = target, bet.prior = 1 - target), 
       cutoff.eli = 0.95, early.stop = 0.95, seed = NULL)
}
\arguments{
\item{design}{option for selecting different designs, which can be set as \code{'CFO'} and \code{'aCFO'}.}

\item{target}{the target DLT rate.}

\item{p.true}{the true DLT rates under the different dose levels.}

\item{init.level}{the dose level assigned to the first cohort. The default value \code{init.level} is 1.}

\item{ncohort}{the total number of cohorts.}

\item{cohortsize}{the number of patients of each cohort.}

\item{prior.para}{the prior parameters for a beta distribution, where set as \code{list(alp.prior = target, bet.prior = 1 - target)}
by default, \code{alp.prior} and \code{bet.prior} represent the parameters of the prior distribution for
the true DLT rate at any dose level. This prior distribution is specified as Beta(\code{alpha.prior}, \code{beta.prior}).}

\item{cutoff.eli}{the cutoff to eliminate overly toxic doses for safety. We recommend
the default value of \code{cutoff.eli = 0.95} for general use.}

\item{early.stop}{the threshold value for early stopping. The default value \code{early.stop = 0.95}
generally works well.}

\item{seed}{an integer to be set as the seed of the random number generator for reproducible results. The default value is set to \code{NULL}.}
}
\value{
The \code{CFO.simu} function returns a list object comprising the following components:
\itemize{
\item target: the target DLT rate.
\item MTD: the selected MTD. \code{MTD = 99} indicates that the simulation is terminated due to early stopping.
\item correct: a binary indicator of whether the recommended dose level matches the correct MTD (1 for yes).
The correct MTD is the dose level at which the true DLT rate is closest to the target DLT rate.
\item npatients: the total number of patients allocated to all dose levels.
\item ntox: the total number of DLTs observed for all dose levels.
\item over.doses: a vector indicating whether each dose is overdosed or not (1 for yes).
\item cohortdose: a vector including the dose level assigned to each cohort.
\item ptoxic: the percentage of subjects assigned to dose levels with a DLT rate greater than the target.
\item patientDLT: a vector including the DLT outcome observed for each patient.
\item sumDLT: the total number of DLT observed.
\item earlystop: a binary indicator of whether the trial is early stopped (1 for yes).
}
}
\description{
In the CFO and aCFO designs for phase I trials, the function is used to conduct one single simulation and find the maximum tolerated dose (MTD).
}
\note{
The \code{CFO.simu()} function is designed to conduct a single CFO or aCFO simulation. If \code{design = 'CFO'}, it corresponds
to the CFO design. If \code{design = 'aCFO'}, it corresponds to the aCFO design. \cr
The early stopping and dose elimination rules are incorporated into the CFO or aCFO design
to ensure patient safety and benefit. If there is substantial evidence indicating that the current dose level
exhibits excessive toxicity, we exclude the current dose level as well as higher dose levels from the trial. If the lowest dose level is overly toxic, the trial will be terminated
according to the early stopping rule. Upon the predefined maximum sample size is reached or the lowest dose
level is over-toxicity, the experiment is concluded, and the MTD is determined using isotonic regression.
}
\examples{
target <- 0.2; ncohort <- 12; cohortsize <- 3; init.level <- 1
p.true <- c(0.01, 0.07, 0.20, 0.35, 0.50, 0.65, 0.80)
### find the MTD for a single CFO simulation
CFOtrial <- CFO.simu(design = 'CFO', target, p.true, init.level, ncohort, cohortsize, seed = 1)
summary(CFOtrial)
plot(CFOtrial)
\donttest{
# This test may take longer than 5 seconds to run
# It is provided for illustration purposes only
# Users can run this code directly
### find the MTD for a single aCFO simulation
aCFOtrial <- CFO.simu(design = 'aCFO', target, p.true, init.level, ncohort, cohortsize, seed = 1)
summary(aCFOtrial)
plot(aCFOtrial)
}
}
\references{
Jin H, Yin G (2022). CFO: Calibration-free odds design for phase I/II clinical trials.
\emph{Statistical Methods in Medical Research}, 31(6), 1051-1066. \cr
Fang J, Yin G (2024). Fractional accumulative calibration‐free odds (f‐aCFO) design for delayed toxicity
in phase I clinical trials. \emph{Statistics in Medicine}.
}
\author{
Jialu Fang, Wenliang Wang, Ninghao Zhang, and Guosheng Yin
}

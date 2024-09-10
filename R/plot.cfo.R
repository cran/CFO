#' Plot the results by other functions
#'
#' Plot the objects returned by other functions, including (1) dose allocation of a single trial;
#' (2) the estimate of toxicity probability for each dose and corresponding 95% credible interval;
#' (3) operating characteristics of multiple simulations, including MTD selection percentage, 
#' the averaged number of patients allocated to different doses in one simulation and the averaged 
#' number of DLT observed for different doses in one simulation.
#' 
#' @param x the object returned by other functions
#' @param ... ignored arguments
#' @param name the name of the object to be plotted.
#'             User does not need to input this parameter.
#'
#' @return \code{plot()} returns a figure or a series of figures depending on the object entered.
#' 
#' @note In the example, we set \code{nsimu = 5} for testing time considerations. In reality, \code{nsimu} 
#'       is typically set to 5000 to ensure the accuracy of the results.
#' 
#' 
#' @author Jialu Fang, Wenliang Wang, Ninghao Zhang, and Guosheng Yin
#' 
#' @importFrom grDevices dev.flush dev.hold devAskNewPage
#' @importFrom graphics axis barplot mtext par plot
#' @importFrom graphics abline arrows legend points
#' @import ggplot2
#' @export
#'
#' @examples
#' 
#' 
#' ## settings for 1dCFO
#' nsimu <- 5; ncohort <- 12; cohortsize <- 3; init.level <- 1
#' p.true <- c(0.02, 0.05, 0.20, 0.28, 0.34, 0.40, 0.44); target <- 0.2
#' assess.window <- 3; accrual.rate <- 2; tte.para <- 0.5; accrual.dist <- 'unif'
#' 
#' ## plot the object returned by CFO.simu()
#' CFOtrial <- CFO.simu(design = 'CFO', target, p.true, init.level, ncohort, cohortsize, seed = 1)
#' plot(CFOtrial)
#' 
#' ## plot the object returned by CFO.selectmtd()
#' selmtd <- CFO.selectmtd(target=0.2, npts=c(3,3,27,3,0,0,0), ntox=c(0,0,4,2,0,0,0))
#' plot(selmtd)
#' 
#' \donttest{
#' # This test may take longer than 5 seconds to run
#' # It is provided for illustration purposes only
#' # Users can run this code directly
#' 
#' ## plot the object returned by lateonset.simu()
#' ## f-aCFO design
#' faCFOtrial <- lateonset.simu (design = 'f-aCFO', target, p.true, init.level,  
#'                 ncohort, cohortsize, assess.window, tte.para, accrual.rate, accrual.dist, seed = 1)
#' plot(faCFOtrial)
#' 
#' ## summarize the object returned by CFO.oc()
#' faCFOoc <- CFO.oc (nsimu, design = 'f-aCFO', target, p.true, init.level, ncohort, cohortsize,
#'         assess.window, tte.para, accrual.rate, accrual.dist, seeds = 1:nsimu)
#' plot(faCFOoc)
#' 
#' ## settings for 2dCFO
#' p.true <- matrix(c(0.05, 0.10, 0.15, 0.30, 0.45,
#'                    0.10, 0.15, 0.30, 0.45, 0.55,
#'                    0.15, 0.30, 0.45, 0.50, 0.60), 
#'                  nrow = 3, ncol = 5, byrow = TRUE)
#' target <- 0.3; ncohort <- 12; cohortsize <- 3
#' 
#' ## plot the single simulation returned by CFO2d.simu()
#' CFO2dtrial <- CFO2d.simu(target, p.true, init.level = c(1,1), ncohort, cohortsize, seed = 1)
#' plot(CFO2dtrial)
#' 
#' ## plot the multiple simulation returned by CFO2d.oc()
#' CFO2doc <- CFO2d.oc(nsimu = 5, target, p.true, init.level = c(1,1), ncohort, cohortsize, 
#'                     seeds = 1:5)
#' plot(CFO2doc)
#' 
#' ## select a MTD based on the trial data
#' ntox <- matrix(c(0, 0, 2, 0, 0, 0, 2, 7, 0, 0, 0, 2, 0, 0, 0), nrow = 3, ncol = 5, byrow = TRUE)
#' npts <- matrix(c(3, 0, 12, 0, 0, 3, 12, 24, 0, 0, 3, 3, 0, 0, 0), nrow = 3, ncol = 5, byrow = TRUE)
#' selmtd <- CFO2d.selectmtd(target=0.3, npts=npts, ntox=ntox)
#' plot(selmtd)
#' 
#' ## summarize the object returned by CFOeff.next()
#' decision <- CFOeff.next(target=0.4,txs=c(3,1,7,11,26),tys=c(0,0,0,0,6),
#'               tns= c(6, 3, 12, 17, 36), currdose = 3, mineff = 0.3)
#' plot(decision)
#' 
#' ## summarize the object returned by CFOeff.simu()
#' target <- 0.30; mineff <- 0.30
#' prior.para = list(alp.prior = target, bet.prior = 1 - target, 
#'                   alp.prior.eff = 0.5, bet.prior.eff = 0.5)
#' p.true=c(0.05, 0.07, 0.1, 0.12, 0.16)
#' pE.true=c(0.35, 0.45, 0.5, 0.55, 0.75)
#' result <- CFOeff.simu(target, p.true, pE.true, ncohort, init.level, cohortsize,
#'                        prior.para, mineff = mineff, seed = 1)
#' plot(result)
#' 
#' ## summarize the object returned by CFOeff.oc()
#' nsimu = 10
#' result <- CFOeff.oc(target, p.true, pE.true, prior.para, 
#'           init.level,cohortsize, ncohort, nsimu, mineff = mineff, seeds = 1:nsimu)
#' plot(result)
#' 
#'}
#' 

plot.cfo<- function (x,..., name = deparse(substitute(x)))
{
  new.obj = unlist(strsplit(name, split = "\\$"))
  strpattern = "none"
  if (length(new.obj) >= 2){
    strpattern = new.obj[2]
  }
  assign("objectPlot", get(new.obj[1]))
  if (!is.element(strpattern, c("none", names(objectPlot)))) {
    warning("Please double check and specify the variable to be plotted...\n")
  }
  else {
    ###############################################################################
    ############################plot for CFO.oc()###############################
    ###############################################################################
    if (!is.null(objectPlot$simu.setup)) { #plot for CFOeff.oc
      oldpar <- par(no.readonly = TRUE) 
      on.exit(par(oldpar))
      if(!is.null(objectPlot$class)) {
        attributesToPlot <- c("selpercent", "npatients", "ntox", "neff")
        titles <- c("OBD selection", "Average patients allocation", "Average DLT observed", "Average efficacy outcome observed")
        ylabels <- c("Percentage (%)", "Number of patients", "Number of DLTs", "Number of efficacy")
        
        par(mfrow = c(2, 2))
        
        # Loop through each attribute and create a plot
        for (i in seq_along(attributesToPlot)) {
          attr <- attributesToPlot[i]
          # Check if the attribute exists in the objectPlot
          if (!is.null(objectPlot[[attr]])) {
            # Extract the vector
            vectorToPlot <- objectPlot[[attr]]
            
            # Convert to percentages only for selPercent
            if (attr == "selpercent") {
              vectorToPlot <- vectorToPlot * 100
            }
            
            # Create the bar plot with horizontal x-axis labels
            bplot <- barplot(vectorToPlot, ylab = ylabels[i], main = titles[i], xlab = "Dose level",
                             cex.names = 1, xaxt = "n", ylim = c(0, max(vectorToPlot))*1.3,
                             cex.lab = 1.3)
            axis(1, at = bplot, labels = seq(1, length(objectPlot[[attr]])))
          }
        }
      }
      else if(is.null(dim(objectPlot$selpercent))){
        attributesToPlot <- c("selpercent", "npatients", "ntox")
        titles <- c("MTD selection", "Average patients allocation", "Average DLT observed")
        ylabels <- c("Percentage (%)", "Number of patients", "Number of DLTs")
        
        par(mfrow = c(3, 1))
        
        # Loop through each attribute and create a plot
        for (i in seq_along(attributesToPlot)) {
          attr <- attributesToPlot[i]
          # Check if the attribute exists in the objectPlot
          if (!is.null(objectPlot[[attr]])) {
            # Extract the vector
            vectorToPlot <- objectPlot[[attr]]
            
            # Convert to percentages only for selPercent
            if (attr == "selpercent") {
              vectorToPlot <- vectorToPlot * 100
            }
            
            # Create the bar plot with horizontal x-axis labels
            bplot <- barplot(vectorToPlot, ylab = ylabels[i], main = titles[i], xlab = "Dose level",
                             cex.names = 1, xaxt = "n", ylim = c(0, max(vectorToPlot))*1.3,
                             cex.lab = 1.3)
            axis(1, at = bplot, labels = seq(1, length(objectPlot[[attr]])))
          }
        }
      }
      else if(length(dim(objectPlot$selpercent))==2) {
        attributesToPlot <- c("selpercent", "npatients", "ntox")
        titles <- c("MTD selection", "Average patients allocation", "Average DLT observed")
        ylabels <- c("Percentage (%)", "Number of patients", "Number of DLTs")
        
        par(mfrow = c(3, 1))
        
        # Loop through each attribute and create a plot
        for (i in seq_along(attributesToPlot)) {
          attr <- attributesToPlot[i]
          # Check if the attribute exists in the objectPlot
          if (!is.null(objectPlot[[attr]])) {
            # Extract the matrix
            matrixToPlot <- objectPlot[[attr]]
            
            # Convert the matrix to a vector by column
            matrixVector <- as.vector(matrixToPlot)
            
            # Convert to percentages only for selpercent
            if (attr == "selpercent") {
              matrixVector <- matrixVector * 100
            }
            
            # Create x-axis labels
            dimMatrix <- dim(matrixToPlot)
            xLabels <- expand.grid(row = 1:dimMatrix[1], col = 1:dimMatrix[2])
            xLabels <- apply(xLabels, 1, function(x) paste("(", x[1], ",", x[2], ")", sep = ""))
            
            # Create the bar plot with horizontal x-axis labels
            barplot(matrixVector, names.arg = xLabels, las = 2,
                    xlab = "Combined dose level", ylab = ylabels[i], main = titles[i])
          }
        }
      }
    }
    
    
    ###############################################################################
    #########################plot for XXX.simu()###################################
    ###############################################################################
    else if (!is.null(objectPlot$correct)) { 
      if (length(objectPlot$OBD) == 1){#plot for CFOeff
        if (objectPlot$OBD == 99){
          if (objectPlot$stopreason == "overly_toxic"){
            warning("All tested doses are overly toxic. No OBD should be selected! \n\n")
          }else if (objectPlot$stopreason == "low_efficacy"){
            warning("All tested doses show low efficacy. No OBD should be selected! \n\n")
          }
        }else{
          dose <- objectPlot$cohortdose
          DLT <- objectPlot$patientDLT
          EFF <- objectPlot$patienteff
          ncohort <- length(objectPlot$cohortdose)
          cohortsize <- sum(objectPlot$npatients)/ncohort
          
          # Generate y_labels
          y_labels <- seq(1, max(dose))
          
          # Generate sequences for each patient
          sequences <- 1:(ncohort * cohortsize)
          
          # Generate dose_levels for each patient
          dose_levels <- rep(dose, each = cohortsize)
          
          # Generate DLT_observed for each patient
          DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
          
          # Generate EFF_observed for each patient
          EFF_observed <- matrix(EFF, nrow = cohortsize, ncol = ncohort)
          
          dfDLT <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
          dfEFF <- data.frame(sequence = sequences, dose_levels = dose_levels, EFF_observed = EFF_observed)
          
          p <- ggplot() +
            # Plot toxicity data
            geom_point(data = dfDLT, aes(x = sequence, y = dose_levels + 0.05, fill = as.factor(DLT_observed)), 
                       size = 2.5, color = "black", shape = 22) +
            # Plot efficacy data
            geom_point(data = dfEFF, aes(x = sequence, y = dose_levels - 0.05, shape = as.factor(EFF_observed)), 
                       size = 2.5, color = "black") +
            scale_fill_manual(values = c('white', 'black'), labels = c('DLT not observed', 'DLT observed')) +
            scale_shape_manual(values = c(13, 21), labels = c('No efficacy', 'Efficacy')) +
            scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
            labs(x = "Sequence of patients treated", 
                 y = "Dose level",
                 fill = 'DLT observed') +
            theme_minimal() +
            theme(legend.position = "bottom",
                  legend.title = element_blank(),
                  legend.key = element_rect(fill = NA, color = "black"))
          # Display the plot
          print(p)
        }
      }
      else if (length(objectPlot$MTD) == 1) {
        if (objectPlot$MTD == 99) {
          warning("All tested doses are overly toxic. No MTD should be selected! \n\n")
        }
        else {
          if (!is.null(objectPlot$totaltime)){ #plot for lateonset.simu()
            dose <- objectPlot$cohortdose
            DLT <- objectPlot$patientDLT
            ncohort <- length(objectPlot$cohortdose)
            cohortsize <- sum(objectPlot$npatients)/ncohort
            
            # Generate y_labels
            y_labels <- seq(1, max(dose))
            
            # Generate sequences for each patient
            sequences <- objectPlot$entertimes
            
            # Generate dose_levels for each patient
            dose_levels <- rep(dose, each = cohortsize)
            
            # Generate DLT_observed for each patient
            DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
            
            new_seq <- ifelse(objectPlot$DLTtimes!=0, sequences+objectPlot$DLTtimes, NA)
            new_y <- ifelse(objectPlot$DLTtimes!=0, dose_levels, NA)
            
            add_noise <- function(vec) {
              counts <- table(vec)
              counts <- table(names(counts))
              result <- numeric(length(vec))
              for (i in seq_along(vec)) {
                if (!is.na(vec[i])) {
                  result[i] <- vec[i] + 0.1 * counts[as.character(vec[i])]  # add 0.05 for unique value
                  counts[as.character(vec[i])] <- counts[as.character(vec[i])] + 1
                }
              }
              return(result)
            }
            
            new_y <- add_noise(new_y)
            
            null_data <- rep(NA, cohortsize)
            sequences <- c(sequences, null_data)
            dose_levels <- c(dose_levels, null_data)
            DLT_observed <- cbind(DLT_observed,rep(2, cohortsize))
            
            df <- data.frame(sequences = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
            dfnew <- data.frame(sequences = as.vector(na.omit(sequences)), dose_levels = as.vector(na.omit(dose_levels)), new_seq = new_seq, new_y = new_y)
            dfnew <- na.omit(dfnew)
            
            suppressWarnings({
              # Create the plot
              p <- ggplot(df, aes(x = sequences, y = dose_levels)) +
                geom_point(aes(shape = factor(DLT_observed,levels=c(0,1,2))), color = 'black', size = 2.5) +
                geom_step(direction = 'hv', color = 'black') +
                scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
                labs(x = "Time (in months)", 
                     y = "Dose level",
                     fill = 'DLT observed') +
                theme_minimal() +
                theme(text = element_text(size = 12), 
                      legend.title=element_blank(), 
                      legend.position = 'top', legend.margin = margin(0, 0, 0, 0)) +
                scale_shape_manual(values = c(1, 16, 4), 
                                   labels = c('DLT not observed', 'DLT observed', 'DLT time'), 
                                   drop = FALSE) 
              
              for (row in 1:(nrow(dfnew))){
                xuse=c(dfnew[row,"sequences"],dfnew[row,"new_seq"])
                yuse=c(dfnew[row,"dose_levels"],dfnew[row,"new_y"])
                dfuse <-data.frame(xuse=xuse, yuse=yuse)
                p <- p + 
                  annotate("point", x = xuse[2], y = yuse[2], shape = 4,size = 2.5) +
                  geom_step(aes(x = xuse, y = yuse), data = dfuse,direction = 'vh',
                            linetype = 2)
              }
              print(p)})
          }
          
          
          else{ #plot for CFO.simu()
            if (objectPlot$MTD == 99) {
              warning("All tested doses are overly toxic. No MTD should be selected! \n\n")
            }
            else {
              dose <- objectPlot$cohortdose
              DLT <- objectPlot$patientDLT
              ncohort <- length(objectPlot$cohortdose)
              cohortsize <- sum(objectPlot$npatients)/ncohort
              
              # Generate y_labels
              y_labels <- seq(1, max(dose))
              
              # Generate sequences for each patient
              sequences <- 1:(ncohort * cohortsize)
              
              # Generate dose_levels for each patient
              dose_levels <- rep(dose, each = cohortsize)
              
              # Generate DLT_observed for each patient
              DLT_observed <- matrix(DLT, nrow = cohortsize, ncol = ncohort)
              
              df <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
              
              # Create the plot
              p <- ggplot(df, aes(x = sequence, y = dose_levels)) +
                geom_point(aes(fill = as.factor(DLT_observed)), color = 'black', shape = 21, size = 2.5) +
                geom_step(direction = 'hv', color = 'black') +
                scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
                labs(x = "Sequence of patients treated", 
                     y = "Dose level",
                     fill = 'DLT observed') +
                theme_minimal() +
                theme(text = element_text(size = 12), legend.title=element_blank(), 
                      legend.position= 'top', legend.margin = margin(0, 0, 0, 0)) +
                scale_fill_manual(values = c('white', 'black'), labels = c('DLT not observed', 'DLT observed'))
              
              # Display the plot
              print(p)
            }
          }
        }
      }
      
      else{
        if (objectPlot$MTD[1] == 99 | objectPlot$MTD[2] == 99) {
          warning("All tested doses are overly toxic. No MTD should be selected! \n\n")
        }
        else {
          dose <- objectPlot$cohortdose
          DLT <- objectPlot$patientDLT   ###need to change!!!!
          ncohort <- dim(objectPlot$cohortdose)[1]
          cohortsize <- sum(objectPlot$npatients)/ncohort
          dim <- dim(objectPlot$ntox)
          
          # Generate y_labels
          y_labels <- expand.grid(1:dim[1], 1:dim[2])
          y_labels <- apply(y_labels, 1, function(x) paste('(', x[1], ',', x[2], ')'))
          
          # Generate sequences for each patient
          sequences <- 1:(ncohort * cohortsize)
          
          # Generate dose_levels for each patient
          dose_levels <- rep(match(apply(dose, 1, function(x) paste('(', x[1], ',', x[2], ')')), y_labels), each = cohortsize)
          
          # Generate DLT_observed for each patient
          DLT_observed <- t(DLT)
          
          df <- data.frame(sequence = sequences, dose_levels = dose_levels, DLT_observed = DLT_observed)
          
          # Create the plot
          p <- ggplot(df, aes(x = sequence, y = dose_levels)) +
            geom_point(aes(fill = as.factor(DLT_observed)), color = 'black', shape = 21, size = 2.5) +
            geom_step(direction = 'hv', color = 'black') +
            scale_y_continuous(breaks = 1:length(y_labels), labels = y_labels) +
            labs(x = "Sequence of patients treated", 
                 y = "Combined dose level",
                 fill = 'DLT observed') +
            theme_minimal() +
            theme(text = element_text(size = 12), legend.title=element_blank(), 
                  legend.position = 'top', legend.margin = margin(0, 0, 0, 0)) +
            scale_fill_manual(values = c('white', 'black'), labels = c('DLT not observed', 'DLT observed'))
          # Display the plot
          print(p)
        }
      }
    }
    
    
    ###############################################################################
    #########################plot for CFO.selectmtd()###################################
    ###############################################################################
    else if (!is.null(objectPlot$p_est)){
      if (objectPlot$MTD[1] == 99) {
        warning("All tested doses are overly toxic. No MTD is selected!\n")
      }
      else {
        if (!is.null(objectPlot$p_est)) {
          
          if (length(objectPlot$MTD) >= 2) {
            p_est.comb=objectPlot$p_est
            rownames(p_est.comb)=1:dim(p_est.comb)[1]
            colnames(p_est.comb)=1:dim(p_est.comb)[2]
            barplot(p_est.comb,beside=TRUE,ylab="DLT rate",
                    ylim=c(0,round(max(p_est.comb,na.rm=TRUE)*1.5,1)),xlab="Drug B",legend.text=rownames(p_est.comb),
                    args.legend=list(title="Drug A",horiz=TRUE,x="top"))
          }
          else {
            p_est = objectPlot$p_est
            p_hat = p_est[, 2]
            ci = p_est[, 3]
            ci = gsub("[\\(\\)]", "", ci)
            conf.intv = matrix(unlist(strsplit(ci, ",")),
                               byrow = TRUE, ncol = 2)
            if (p_est[1, 2] == "----") {
              warning("The trial is stopped since the lowest dose is too toxic.\n")
            }
            else {
              numbs = ifelse(sum(p_hat == "----") ==
                               0, length(p_hat), min(which(p_hat ==
                                                             "----")) - 1)
              numbs2 = length(p_hat)
              phatx = as.numeric(as.character(p_hat[1:numbs]))
              lwr = as.numeric(as.character(conf.intv[1:numbs,
                                                      1]))
              upr = as.numeric(as.character(conf.intv[1:numbs,
                                                      2]))
              plot(1:numbs2, ylim = c(0, 1), xlab = "Dose level",
                   ylab = "DLT rate", pch = "", xaxt = "n",
                   cex.lab = 1.3)
              axis(1, at = 1:numbs2, labels = 1:numbs2)
              abline(h = objectPlot$target, lty = 2,
                     col = 2)
              points(1:numbs, phatx, pch = 19)
              arrows(x0 = 1:numbs, x1 = 1:numbs, y0 = lwr,
                     y1 = upr, code = 3, angle = 90, length = 0.1)
              if (numbs < numbs2) {
                points((numbs + 1):numbs2, seq(min(1,
                                                   max(phatx, na.rm = T) + 0.05), min(max(phatx,
                                                                                          na.rm = T) + 0.2, 1), length = numbs2 -
                                                 numbs), pch = "*", cex = 1.5)
                legend("topleft", "*   no patient treated")
              }
            }
          }
        }
      }
      
      
    }
    
    ###############################################################################
    #########################plot for others###################################
    ###############################################################################
    else{
      warning("The variable cannot be plotted. \n 
              Please double check and specify the variable to be plotted...\n")
    }
  }
}


################################################################################
################################################################################
################################################################################
##
## Name: FMradio
## Author: Carel F.W. Peeters
##
## Maintainer: Carel F.W. Peeters
##             Statistics for Omics Research Unit
##             Dept. of Epidemiology & Biostatistics
##             Amsterdam Public Health research institute
##             Amsterdam University medical centers,
##             Location VU University medical center
##             Amsterdam, the Netherlands
## Email:	     cf.peeters@vumc.nl
##
## Version: 1.1
## Last Update:	29/04/2019
## Description:	Pipeline (support) for prediction with radiomic data compression
##
################################################################################
################################################################################
################################################################################



##------------------------------------------------------------------------------
##
## Hidden Support Functions
##
##------------------------------------------------------------------------------

.corLW <- function(R, lambda){
  ##############################################################################
  # Function that calculates the Ledoit-Wolf type regularized correlation
  # matrix for a given penalty value
  # - R      > correlation matrix
  # - lambda > value penalty parameter
  ##############################################################################
  
  # Dependencies:
  # require("base")
  
  # Determine and return
  return((1-lambda) * R + lambda * diag(dim(R)[1]))
}



.LL <- function(S, R){
  ##############################################################################
  # Function that computes the value of the (negative) log-likelihood
  # - S > sample correlation matrix
  # - R > (possibly regularized) correlation matrix
  ##############################################################################
  
  # Dependencies:
  # require("base")
  
  # Evaluate and return
  LL <- log(det(R)) + sum(S*solve(R))
  return(LL)
}



.kcvl <- function(lambda, X, folds){
  ##############################################################################
  # Function that calculates a cross-validated negative log-likelihood score
  # for a single penalty value
  # - lambda > value penalty parameter
  # - X      > (standardized) data matrix, observations in rows
  # - folds  > cross-validation sample splits
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("stats")
  
  # Evaluate
  cvLL <- 0
  for (f in 1:length(folds)){
    R    <- cor(X[-folds[[f]], , drop = FALSE])
    S    <- cor(X[folds[[f]], , drop = FALSE])
    nf   <- dim(X[folds[[f]], , drop = FALSE])[1]
    cvLL <- cvLL + nf*.LL(S, .corLW(R, lambda))
  }
  
  # Return
  return(cvLL/length(folds))
}



.LH <- function(R, LM, UM){
  ##############################################################################
  # Evaluates term proportional to -2 times the log-likelihood of the FA-model
  # R  > (regularized) covariance or correlation matrix
  # LM > loadings matrix
  # UM > diagonal uniquenesses matrix
  #
  # NOTES:
  # - The UM matrix is assumed to be positive definite
  # - The LM matrix is assumed to be of full column rank
  # - The LM matrix is assumed to have the observed features in the rows
  ##############################################################################
  
  # Dependencies:
  # require("base")
  
  # Evaluate
  Rfit <- LM %*% t(LM) + UM
  LH   <- log(det(Rfit)) + sum(diag(solve(Rfit) %*% R))
  
  # Return
  return(LH)
}



.DF <- function(R, LM, UM){
  ##############################################################################
  # Evaluates the regular FA discrepancy function
  # R  > (regularized) covariance or correlation matrix
  # LM > loadings matrix
  # UM > diagonal uniquenesses matrix
  #
  # NOTES:
  # - The UM matrix is assumed to be positive definite
  # - The LM matrix is assumed to be of full column rank
  # - The LM matrix is assumed to have the observed features in the rows
  ##############################################################################

  # Dependencies:
  # require("base")
    
  # Evaluate
  DF <- .LH(R, LM, UM) - log(det(R)) - dim(R)[1]
  
  # Return
  return(DF)
}



.IC <- function(R, LM, UM, m, n, type = "BIC"){
  ##############################################################################
  # Evaluates the regular FA discrepancy function
  # R    > (regularized) covariance or correlation matrix
  # LM   > loadings matrix
  # UM   > diagonal uniquenesses matrix
  # m    > dimension of latent vector
  # n    > sample size
  # type > character specifying the penalty type: either BIC or AIC
  #
  # NOTES:
  # - The UM matrix is assumed to be positive definite
  # - The LM matrix is assumed to be of full column rank
  # - The LM matrix is assumed to have the observed features in the rows
  ##############################################################################

  # Dependencies:
  # require("base")
  
  # Fit determination
  p   <- dim(R)[1]
  fit <- n * (p * log(2 * pi) + .LH(R, LM, UM))
  
  # Penalty determination
  pars <- p*(m + 1) - (m*(m - 1))/2
  if (type == "BIC"){
    pen <- log(n) * pars
  }
  if (type == "AIC"){
    pen <- 2 * pars
  }
  
  # Evaluate
  IC <- fit + pen
  
  # Return
  return(IC)
}



.rank <- function(R){
  ##############################################################################
  # Evaluates the rank of a symmetric (correlation) matrix
  # R > symmetric (correlation) matrix
  #
  # NOTES:
  # - Determines the rank as the number of eigenvalues that equal or exceed the
  #   threshold set at the tolerance times the largest eigenvalue
  # - The tolerance is determined as the number of features times the machine 
  #   epsilon
  # - Corresponds to the default method of the rankMatrix function from the
  #   Matrix package, which itself corresponds to the MATLAB default
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("stats")
  
  # Evaluate rank
  eval <- Re(eigen(R)$values)
  tol  <- dim(R)[1] * .Machine$double.eps
  rank <- sum(eval >= tol * eval[1])
  
  # Return
  return(rank)
}



.tr <- function(M){
  ##############################################################################
  # - Internal function to compute the trace of a matrix
  # - M > matrix input
  ##############################################################################
  
  return(sum(diag(M)))
}




##------------------------------------------------------------------------------
##
## Function for heatmap visualization
##
##------------------------------------------------------------------------------

if (getRversion() >= "2.15.1") utils::globalVariables(c("X1", "X2", "value"))

radioHeat <- function(R, lowColor = "blue", highColor = "red", labelsize = 10,
                      diag = TRUE, threshold = FALSE, threshvalue = .95,
                      values = FALSE, textsize = 10, legend = TRUE, main = ""){
  ##############################################################################
  # Dedicated heatmapping of radiomics-based correlation matrix
  # - R           > (regularized) correlation matrix
  # - lowColor    > determines color scale in the negative range, default = "blue"
  # - highColor   > determines color scale in the positive range, default = "red"
  # - labelsize   > set textsize row and column labels, default = 10
  # - diag        > logical determining treatment diagonal elements R. If FALSE,
  #                 then the diagonal elements are given the midscale color of
  #                 white; only when R is a square matrix
  # - threshold   > logical determining if only values above a certain 
  #                 (absolute) threshold should be visualized
  # - threshvalue > value used when threshold = TRUE
  # - values      > optional inclusion of cell-values, default = FALSE
  # - textsize    > set textsize cell values
  # - legend      > optional inclusion of color legend, default = TRUE
  # - main        > character specifying the main title, default = ""
  # 
  # NOTES:
  # - The thresholding arguments can be used to (i) visually assess redundancy
  #   and (ii) check if the redundancy filter has done its work
  # - The cell-values returned when values = TRUE are the values after possible
  #   thresholding
  ##############################################################################
  
  # Dependencies
  # require("base")
  # require("ggplot2")
  # require("reshape")
  
  # Checks
  if (!is.matrix(R)){
    stop("Supply 'R' as matrix")
  }
  if (class(lowColor) != "character"){
    stop("Input (lowColor) is of wrong class")
  }
  if (length(lowColor) != 1){
    stop("Length lowColor must be one")
  }
  if (class(highColor) != "character"){
    stop("Input (highColor) is of wrong class")
  }
  if (length(highColor) != 1){
    stop("Length highColor must be one")
  }
  if (class(labelsize) != "numeric"){
    stop("Input (labelsize) is of wrong class")
  }
  else if (length(labelsize) != 1){
    stop("Length labelsize must be one")
  }
  if (labelsize <= 0){
    stop("labelsize must be positive")
  }
  if (class(diag) != "logical"){
    stop("Input (diag) is of wrong class")
  }
  if (class(threshold) != "logical"){
    stop("Input (threshold) is of wrong class")
  }
  if (class(legend) != "logical"){
    stop("Input (legend) is of wrong class")
  }
  if (class(main) != "character"){
    stop("Input (main) is of wrong class")
  }
  if (class(values) != "logical"){
    stop("Input (values) is of wrong class")
  }
  
  # Conditional checks
  if (threshold){
    if(class(threshvalue) != "numeric"){
      stop("Input (threshvalue) is of wrong class")
    }
    if(length(threshvalue) != 1){
      stop("Length input (threshvalue) must be one")
    }
    if(threshvalue < 0){
      stop("Input (threshvalue) cannot be negative")
    }
  }
  if (values){
    if(class(textsize) != "numeric"){
      stop("Input (textsize) is of wrong class")
    }
    if(length(textsize) != 1){
      stop("Length input (textsize) must be one")
    }
    if(textsize <= 0){
      stop("Input (textsize) must be positive")
    }
  }
  
  # Thresholding
  if (threshold){
    R[abs(R) < threshvalue] <- 0
  }
  
  # Put matrix in data format
  if (nrow(R) == ncol(R) & !diag) {diag(R) <- 0}
  Mmelt    <- melt(R)
  Mmelt$X1 <- factor(as.character(Mmelt$X1), 
                     levels = unique(Mmelt$X1), ordered = TRUE)
  Mmelt$X2 <- factor(as.character(Mmelt$X2), 
                     levels = unique(Mmelt$X2), ordered = TRUE)
  
  # Visualize
  ggplot(Mmelt, aes(X2, X1, fill = value)) + 
    geom_tile() +
    xlab(" ") + 
    ylab(" ") +    
    ylim(rev(levels(Mmelt$X1))) +
    ggtitle(main) +
    theme(axis.ticks = element_blank()) +
    theme(axis.text.y = element_text(size = labelsize)) +
    theme(axis.text.x = element_text(angle = -90, 
                                     vjust = .5,
                                     hjust = 0, 
                                     size = labelsize)) +
    scale_fill_gradient2("", low = lowColor,  
                         mid = "white",
                         high = highColor, 
                         midpoint = 0) +
    {if (values) geom_text(aes(label = round(value, 3)), size = textsize)} +
    {if (!legend) theme(legend.position = "none")}
}




##------------------------------------------------------------------------------
## 
## Function for Redundancy Filtering
##
##------------------------------------------------------------------------------

RF <- function(R, t = .95){
  ##############################################################################
  # Performs redundancy filtering (RF) of square (correlation) matrix
  # R > (correlation) matrix
  # t > absolute value for thresholding
  #
  # NOTES:
  # - When the input matrix R is a correlation matrix, then t should satisfy
  #   -1 < t < 1, for the return matrix to be sensical for further analysis
  ##############################################################################
  
  # Dependencies:
  # require("base")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (class(t) != "numeric"){
    stop("Input (t) is of wrong class")
  }
  if (length(t) != 1){
    stop("Length input (t) must be one")
  }
  
  # Needed
  GO = TRUE
  HM <- numeric()
  
  # Loop
  while(GO){
    for(j in 1:dim(R)[1]){
      HM[j] <- sum(abs(R[j,, drop = FALSE]) >= t)
    }
    Mx <- which(HM == max(HM))
    if(max(HM) < 2){
      GO = FALSE
    } else {
      R  <- R[-Mx[1],-Mx[1]]
      HM <- numeric()
    }
  }
  
  # Return
  return(R)
}




##------------------------------------------------------------------------------
## 
## Function for Subsetting
##
##------------------------------------------------------------------------------

subSet <- function(X, Rf){
  ##############################################################################
  # Convenience function subsetting data matrix or expression set
  # Subsets to retain only the features given in a filtered object
  # X  > data matrix or expression set
  # Rf > filtered correlation matrix
  #
  # NOTES:
  # - If X is a matrix, the observations are expected to be in the rows
  # - Subsetted data are needed for penalty parameter selection for 
  #   the regularized correlation estimator
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("Biobase")
  
  # Checks
  if (!inherits(X, "matrix") & !inherits(X, "ExpressionSet")){
    stop("Input (X) should be either of class 'matrix' or 'ExpressionSet'")
  }
  if (!is.matrix(Rf)){
    stop("Input (Rf) should be a matrix")
  }
  if (nrow(Rf) != ncol(Rf)){
    stop("Input (Rf) should be square matrix")
  }
  
  # Filter
  if (class(X) == "matrix"){
    These <- which(colnames(X) %in% colnames(Rf))
    return(X[,These])
  }
  if (class(X) == "ExpressionSet"){
    These <- which(featureNames(X) %in% colnames(Rf))
    return(X[These,])
  }
}




##------------------------------------------------------------------------------
## 
## Function for Regularized Correlation Matrix Estimation
##
##------------------------------------------------------------------------------

regcor <- function(X, fold = 5, verbose = TRUE){
  ##############################################################################
  # Function determining the optimal penalty value and, subsequently, the 
  # optimal Ledoit-Wolf type regularized correlation matrix using 
  # K-fold cross validation of the negative log-likelihood
  # X       > centered and scaled (subsetted) data matrix, observations in rows
  # fold    > cross-validation fold, default gives 5-fold CV
  # verbose > logical indicating if function should run silently
  #
  # NOTES:
  # - The estimator is a Ledoit-Wolf type estimator
  # - The K-fold cv procedure makes use of the Brent algorithm
  ##############################################################################
  
  # Dependencies:
  # require("base")
  # require("stats")
  
  # Checks
  if (!is.matrix(X)){
    stop("Input (X) should be a matrix")
  }
  if (class(fold) != "numeric" & class(fold) != "integer"){ 
    stop("Input (fold) is of wrong class") 
  }
  if (fold <=  1){ 
    stop("Input (fold) should be at least 2") 
    }
  if (fold > nrow(X)){ 
    stop("Input (fold) cannot exceed the sample size") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  
  # List with K-folds
  if (verbose){cat("Determining folds...", "\n")}
  fold    <- max(min(ceiling(fold), nrow(X)), 2)
  fold    <- rep(1:fold, ceiling(nrow(X)/fold))[1:nrow(X)]
  shuffle <- sample(1:nrow(X), nrow(X))
  folds   <- split(shuffle, as.factor(fold))
  
  # Determine optimal penalty value
  if (verbose){cat("Determining optimal penalty value...", "\n")}
  optLambda <- optim(.5, .kcvl, method = "Brent", lower = 0,
                     upper = 1, X = X, folds = folds)$par
  
  # Return
  return(list(optPen = optLambda, optCor = .corLW(cor(X), optLambda)))
}




##------------------------------------------------------------------------------
## 
## Function for Assessing Factorability
##
##------------------------------------------------------------------------------

SA <- function(R){
  ##############################################################################
  # Calculate Kaiser-Meyer-Olkin measure of feature-sampling adequacy
  # R > (regularized) covariance or correlation matrix
  #
  # NOTES:
  # - The KMO index provides a practical measure for the assessment of 
  #   factorability
  # - Factorability refers to the ability to identify coherent latent features
  # - A KMO index equalling or exceeding .9 would be considered to indicate
  #   great factorability
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  
  # Preliminaries
  Ri      <- solve(R)
  Sc      <- diag(sqrt((1/diag(Ri))))
  P       <- Sc %*% Ri %*% Sc
  diag(R) <- 0
  diag(P) <- 0
  
  # Calculate overall KMO index
  KMO <- sum(R^2)/(sum(R^2) + sum(P^2))
  
  # Calculate KMO index per feature
  KMOfeature <- colSums(R^2)/(colSums(R^2) + colSums(P^2))
  
  # Return
  return(list(KMO = KMO, KMOfeature = KMOfeature))
}




##------------------------------------------------------------------------------
## 
## Functions for Factor Analytic Dimensionality Assessment
##
##------------------------------------------------------------------------------

dimGB <- function(R, graph = TRUE, verbose = TRUE){
  ##############################################################################
  # Calculates the first, second, and third Guttman lower-bounds
  # R        > (regularized) correlation matrix
  # graph    > logical indicating if the results should be visualized
  # verbose  > logical indicating if function should run silently
  #
  # NOTES:
  # - Output corresponds to the Guttman bounds for the minimum rank of a 
  #   reduced correlation matrix
  # - This minimum rank can be used to make a choice on the dimensionality
  #   of the latent vector
  # - The first lower-bound corresponds to kaiser's rule
  # - Note that there is an ordering in the bounds in the sense that:
  #   bound 1 <= bound 2 <= bound 3
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("graphics")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (class(graph) != "logical"){ 
    stop("Input (graph) is of wrong class") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  
  # Preliminaries
  p  <- dim(R)[1]
  Ip <- diag(p)
  
  # Lower-bound 1
  if (verbose){cat("Calculating Guttman bounds...", "\n")}
  LB1e <- eigen(R - Ip)$values
  LB1  <- sum(LB1e > 0)
  
  # Lower-bound 2
  rmax <- apply((abs(R) - Ip), 2, max)
  rmax <- 1 - rmax^2
  LB2e <- eigen(R - diag(rmax))$values
  LB2  <- sum(LB2e > 0)
  
  # Lower-bound 3
  Ri   <- solve(R)
  Ris  <- diag(1/diag(Ri))
  LB3e <- eigen(R - Ris)$values
  LB3  <- sum(LB3e > 0)
  
  # Summarize
  LBT        <- as.table(c(LB1, LB2, LB3))
  names(LBT) <- c("First.lower-bound", "Second.lower-bound", "Third.lower-bound")
  
  # Visualize
  if (graph){
    if (verbose){cat("Visualizing...", "\n")}
    dims <- c(1:p)
    LBe <- c(LB1e, LB2e, LB3e)
    plot(dims, LB1e, axes = FALSE, type = "l", 
         col = "red", xlab = "eigenvalue no.", 
         ylab = "Eigenvalue", main = "",
         ylim = c(min(LBe) - 1, max(LBe)))
    lines(dims, LB2e, col = "blue")
    lines(dims, LB3e, col = "green")
    axis(2, col = "black", lwd = 1)
    axis(1, xlim = c(1,p), col = "black", lwd = 1, tick = TRUE)
    abline(h = 0)
    legend("topright",  
           legend = c("First reduced correlation matrix", 
                      "Second reduced correlation matrix",
                      "Third reduced correlation matrix"), 
           col = c("red", "blue", "green"), 
           lty = 1)
  }
  
  # Return
  return(LBT)
}



dimIC <- function(R, n, maxdim, Type = "BIC", graph = TRUE, verbose = TRUE){
  ##############################################################################
  # Performs dimensionality assessment by way of the AIC or BIC
  # R       > (regularized) covariance or correlation matrix
  # n       > sample size
  # maxdim  > maximum number of latent factors to be assessed
  # Type    > character specifying the penalty type: either BIC or AIC
  # graph   > logical indicating if output should also be visualized
  # verbose > logical indicating if function should run silently
  #
  # NOTES:
  # - maxdim cannot exceed the Ledermann-bound
  # - usually, one wants to set maxdim (much) lower than the Ledermann-bound
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  # require("graphics")
  
  # Preliminaries for checks
  p    <- ncol(R)
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (class(n) != "numeric" & class(n) != "integer"){
    stop("Input (n) is of wrong class")
  }
  if (length(n) != 1){
    stop("Length input (n) must be one")
  }
  if (n <= 1){
    stop("Input (n) must be a strictly positive (numeric) integer")
  }
  if (class(maxdim) != "numeric" & class(maxdim) != "integer"){
    stop("Input (maxdim) is of wrong class")
  }
  if (length(maxdim) != 1){
    stop("Length input (maxdim) must be one")
  }
  if (maxdim <= 1){
    stop("Input (maxdim) cannot be lower than 1")
  }
  if (maxdim > mmax){
    stop("Input (maxdim) is too high")
  }
  if (!(Type %in% c("BIC", "AIC"))){
    stop("Input (Type) should be one of {'BIC', 'AIC'}")
  }
  if (class(graph) != "logical"){ 
    stop("Input (graph) is of wrong class") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  
  # Preliminaries
  IC  <- numeric()
  dim <- seq(from = 1, to = maxdim, by = 1)
  
  # Calculate ICs
  for (j in 1:length(dim)){
    m          <- dim[j]
    fit        <- factanal(factors = m, covmat = R, rotation = "none")
    loadings   <- fit$loadings[1:p,]
    Uniqueness <- diag(fit$uniquenesses)
    IC[j]      <- .IC(R, loadings, Uniqueness, m = m, n = n, type = Type)
    if (verbose){cat(paste("Calculating IC for dimension m = ", 
                           m, " done\n", sep = ""))}
  }
  
  # Visualization
  if (Type == "BIC"){ylabel = "BIC score"}
  if (Type == "AIC"){ylabel = "AIC score"}
  if (graph){
    if (verbose){cat("Visualizing...", "\n")}
    MT <- seq(from = 1, to = maxdim, by = 2)
    plot(dim, IC, axes = FALSE, type = "l", 
         col = "red", xlab = "dimension of factor solution", 
         ylab = ylabel)
    axis(2, ylim = c(min(IC),max(IC)), col = "black", lwd = 1)
    axis(1, xlim = c(1,maxdim), col = "black", lwd = 1, tick = TRUE, at = MT)
  }
  
  # Return object
  ICt <- as.data.frame(cbind(dim,IC))
  colnames(ICt) <- c("latent dimension", ylabel)
  
  ## Return
  return(ICt)
}



dimLRT <- function(R, X, maxdim, rankDOF = TRUE, graph = TRUE, 
                   alpha = .05, Bartlett = FALSE, verbose = TRUE){
  ##############################################################################
  # Performs dimensionality assessment using likelihood ratio testing
  # R        > (regularized) covariance or correlation matrix
  # X        > centered and scaled (subsetted) data matrix, observations in rows
  # maxdim   > maximum number of latent factors to be assessed
  # rankDOF  > logical indicating if the degrees of freedom should be based on
  #            the rank of the raw correlation matrix
  # graph    > logical indicating if output should also be visualized
  # alpha    > numeric giving alpha level, only used when graph = TRUE
  # Bartlett > logical indicating if Bartlett correction should be applied
  # verbose  > logical indicating if function should run silently
  #
  # NOTES:
  # - maxdim cannot exceed the Ledermann-bound
  # - usually, one wants to set maxdim (much) lower than the Ledermann-bound
  # - note that, if p > n, the maximum rank of the raw correlation matrix
  #   is n - 1. In this case there is an alternative Ledermann-bound when 
  #   rankDOF = TRUE. The number of information points in the correlation matrix
  #   is then given as (n-1)*n/2 and this number must exceed 
  #   p*maxdim + p - (maxdim*(maxdim - 1))/2, putting more restrictions on 
  #   maxdim
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  # require("graphics")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (!is.matrix(X)){
    stop("Input (X) should be a matrix")
  }
  
  # Preliminaries for further checks
  p    <- ncol(R)
  n    <- dim(X)[1]
  pr   <- .rank(cor(X))
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  
  # (Conditional) Checks
  if (class(maxdim) != "numeric" & class(maxdim) != "integer"){
    stop("Input (maxdim) is of wrong class")
  }
  if (length(maxdim) != 1){
    stop("Length input (maxdim) must be one")
  }
  if (maxdim <= 1){
    stop("Input (maxdim) cannot be lower than 1")
  }
  if (maxdim > mmax){
    stop("Input (maxdim) is too high")
  }
  if (class(rankDOF) != "logical"){ 
    stop("Input (rankDOF) is of wrong class") 
  }
  if (class(graph) != "logical"){ 
    stop("Input (graph) is of wrong class") 
  }
  if (graph){
    if (class(alpha) != "numeric"){
      stop("Input (alpha) is of wrong class") 
    }
    if (length(alpha) != 1){
      stop("Length input (alpha) must be one")
    }
    if (alpha < 0){
      stop("Input (alpha) should be strictly positive") 
    }
    if (alpha > 1){
      stop("Input (alpha) cannot exceed unity") 
    }
  }
  if (class(Bartlett) != "logical"){ 
    stop("Input (Bartlett) is of wrong class") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  if (p >= n & rankDOF == TRUE){
    if ((pr*(pr + 1))/2 - (p*maxdim + p - (maxdim*(maxdim - 1))/2) < 1){
      stop("Input (maxdim) implies negative degrees of freedom for LRT")
    }
  }
  
  # Preliminaries
  LRT  <- numeric()
  pval <- numeric()
  dim  <- seq(from = 1, to = maxdim, by = 1)
  
  # Calculate LRTs and p-values
  for (j in 1:length(dim)){
    m          <- dim[j]
    fit        <- factanal(factors = m, covmat = R, rotation = "none")
    loadings   <- fit$loadings[1:p,]
    Uniqueness <- diag(fit$uniquenesses)
    if (Bartlett){
      LRT[j] <- (n - 1 - (2 * p + 5)/6 - (2 * m)/3) * .DF(R, loadings, Uniqueness)
    } else {
      LRT[j] <- (n - 1) * .DF(R, loadings, Uniqueness)
    }
    if (rankDOF){
      dof <- (pr*(pr + 1))/2 - (p*m + p - (m*(m - 1))/2)
    } else {
      dof <- ((p - m)^2 - (p + m))/2
    }
    pval[j] <- pchisq(LRT[j], df = dof, lower.tail = FALSE)
    if (verbose){cat(paste("Calculating likelihood ratio test for dimension m = ", 
                           m, " done\n", sep = ""))}
  }
  
  # Visualization
  if (graph){
    if (verbose){cat("Visualizing...", "\n")}
    MT <- seq(from = 1, to = maxdim, by = 2)
    plot(dim, pval, axes = FALSE, type = "l", 
         col = "red", xlab = "dimension of factor solution", 
         ylab = "p-value LR test")
    axis(2, ylim = c(0, 1), col = "black", lwd = 1)
    axis(1, xlim = c(1,maxdim), col = "black", lwd = 1, tick = TRUE, at = MT)
    abline(h = alpha, col = "blue")
  }
  
  # Return object
  LRTt <- as.data.frame(cbind(dim, LRT, pval))
  colnames(LRTt) <- c("latent dimension", "Statistic", "p.value")
  
  ## Return
  return(LRTt)
}




##------------------------------------------------------------------------------
##
## Functions for Assessing the Choice of Factor Dimensionality
##
##------------------------------------------------------------------------------

SMC <- function(R, LM){
  ##############################################################################
  # Juxtaposing SMCs and estimates of communalities
  # R  > (regularized) covariance or correlation matrix
  # LM > (rotated) loadings matrix
  #
  # NOTES:
  # - Squared multiple correlations (SMCs) refer to the proportion of variance 
  #   in feature j that is explained by the remaining p - 1 features
  # - Note that the choice of orthogonal rotation does not affect the 
  #   communality estimates
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (!is.matrix(LM)){
    stop("Input (LM) should be a matrix")    
  } 
  if (ncol(R) != nrow(LM)){
    stop("Column (or row) dimension input (R) incompatible with row-dimension input (LM)")    
  }
  
  # Calculate SMCs and communalities
  SMC           <- 1 - 1/diag(solve(R))
  Communalities <- diag((LM) %*% t(LM))
  
  # Return
  return(cbind(SMC, Communalities))
}



dimVAR <- function(R, maxdim, graph = TRUE, verbose = TRUE){
  ##############################################################################
  # Support function assessing proportion and cumulative variances for a 
  # range of factor solutions
  # R        > (regularized) covariance or correlation matrix
  # maxdim   > maximum number of latent factors to be assessed
  # graph    > logical indicating if output should also be visualized
  # verbose  > logical indicating if function should run silently
  #
  # NOTES:
  # - maxdim cannot exceed the Ledermann-bound
  # - usually, one wants to set maxdim (much) lower than the Ledermann-bound
  # - factor solutions are assessed ranging from 1 to maxdim
  # - SS, proportion variance and cumulative variance are tabulated for each
  #   factor solution
  # - these are based on the unrotated solution
  # - note that the cumulative variance does not depend on the choice of 
  #   (orthogonal) rotation
  # - when graph = TRUE, a graph is returned visualizing the total cumulative
  #   variance against the dimension of the factor solution
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  # require("graphics")
  
  # Preliminaries for checks
  p    <- ncol(R)
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (class(maxdim) != "numeric" & class(maxdim) != "integer"){
    stop("Input (maxdim) is of wrong class")
  }
  if (length(maxdim) != 1){
    stop("Length input (maxdim) must be one")
  }
  if (maxdim <= 1){
    stop("Input (maxdim) cannot be lower than 1")
  }
  if (maxdim > mmax){
    stop("Input (maxdim) is too high")
  }
  if (class(graph) != "logical"){ 
    stop("Input (graph) is of wrong class") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  
  # Preliminaries
  dim <- seq(from = 1, to = maxdim, by = 1)
  L   <- list()
  mCV <- numeric()
  
  # Calculate SS, proportion variance, cumulative variance
  for (j in 1:length(dim)){
    m      <- dim[j]
    fit    <- factanal(factors = m, covmat = R, rotation = "none")
    SS     <- diag(t(fit$loadings)%*%(fit$loadings))
    PV     <- SS/p
    CV     <- cumsum(PV)
    mCV[j] <- max(CV)
    L[[j]] <- rbind(SS,PV,CV)
    names(L)[j] <- paste("dimension =", m)
    if (verbose){cat(paste("Calculating statistics for dimension m = ", 
                           m, " done\n", sep = ""))}
  }
  
  # Visualization
  if (graph){
    if (verbose){cat("Visualizing...", "\n")}
    MT <- seq(from = 1, to = maxdim, by = 2)
    plot(dim, mCV, axes = FALSE, type = "l", 
         col = "red", xlab = "dimension of factor solution", 
         ylab = "total cumulative variance")
    axis(2, ylim = c(0, max(mCV)), col = "black", lwd = 1)
    axis(1, xlim = c(1,maxdim), col = "black", lwd = 1, tick = TRUE, at = MT)
  }
  
  # Return
  return(list(CumVar = mCV, varianceTables = L))
}




##------------------------------------------------------------------------------
##
## Function for Producing ML Factor Solution under m factors
##
##------------------------------------------------------------------------------

mlFA <- function(R, m){
  ##############################################################################
  # Fits an m-factor model through maximum likelihood estimation
  # Produces a Varimax-rotated solution
  # R > (regularized) covariance or correlation matrix
  # m > dimension of factor solution
  # 
  # NOTES:
  # - this function is basically a wrapper around the factanal (stats) function
  # - its purpose is to produce a factor solution of the chosen dimension by
  #   maximum likelihood
  # - the wrapper ensures the model is fitted under the same circumstances
  #   under which latent dimensionality is assessed
  # - produces a Varimax-rotated solution
  # - its output can be used to produce factor scores by the facScore function
  # - note that these cannot be obtained using the factanal function as the
  #   fit is based on the (regularized) correlation matrix
  # - the function assumes that the chosen dimension m is reasonable
  # - the order of the ouput is the order of the features as given in the R 
  #   object
  # - as the 'loadings-part' of the output is of class "loadings" (factanal)
  #   it can be subjected to the print function, sorting the output to 
  #   emphasize the simple structure: e.g.,
  #   print(fit$Loadings, digits = 2, cutoff = .3, sort = TRUE)
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  
  # Preliminaries for checks
  p    <- ncol(R)
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (class(m) != "numeric" & class(m) != "integer"){
    stop("Input (m) is of wrong class")
  }
  if (length(m) != 1){
    stop("Length input (m) must be one")
  }
  if (m <= 1){
    stop("Input (m) cannot be lower than 1")
  }
  if (m > mmax){
    stop("Input (m) is too high")
  }
  
  # Wrapper
  fit <- factanal(factors = m, covmat = R, rotation = "varimax")
  
  # Return
  rotmatrix <- fit$rotmat
  rownames(rotmatrix) <- colnames(fit$loadings)
  colnames(rotmatrix) <- rownames(rotmatrix)
  return(list(Loadings = fit$loadings, 
              Uniqueness = diag(fit$uniquenesses),
              rotmatrix = rotmatrix))
}




##------------------------------------------------------------------------------
## 
## Functions for Obtaining and Assessing Factor Scores
##
##------------------------------------------------------------------------------

facScore <- function(X, LM, UM, type = "thomson"){
  ##############################################################################
  # Finds factor scores from data and given factor-solution
  # Default is Thomson-type factor scores
  # X    > data matrix, observations in rows
  # LM   > (rotated) loadings matrix
  # UM   > diagonal uniquenesses matrix
  # type > character indicating the type of factor score to calculate
  #
  # NOTES:
  # - The input data are assumed to be scaled (or at least centered)
  # - The UM matrix is assumed to be positive definite
  # - The LM matrix is assumed to be of full column rank
  # - The factor-scores obtained the default way are Thomson-type scores
  # - These will be near-orthogonal
  # - These are close to Bartlett-type scores, which are unbiased
  # - Anderson-type scores are orthogonal
  # - The return-object is a dataframe
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("expm")
  
  # Checks
  if (!is.matrix(X)){
    stop("Input (X) should be a matrix")    
  }
  if (!is.matrix(LM)){
    stop("Input (LM) should be a matrix")    
  } 
  if (ncol(X) != nrow(LM)){
    stop("Column-dimension input (X) incompatible with row-dimension input (LM)")    
  }
  if (!is.matrix(UM)){
    stop("Input (UM) should be a matrix")    
  } 
  if (nrow(UM) != ncol(UM)){
    stop("Input (UM) should be square matrix")
  }
  if (ncol(X) != nrow(UM)){
    stop("Column-dimensions input (X) incompatible with dimension input (UM)")    
  }
  if (!(type %in% c("thomson", "bartlett", "anderson"))){
    stop("Input (type) should be one of {'thomson', 'bartlett', 'anderson'}")
  }

  # Needed
  UMi <- diag(1/diag(UM))
  m   <- diag(dim(LM)[2])
  
  # Obtain factor scores
  if (type == "thomson"){
    scores <- as.data.frame(X %*% UMi %*% LM %*% solve(m + t(LM) %*% UMi %*% LM))
  }
  if (type == "bartlett"){
    scores <- as.data.frame(X %*% UMi %*% LM %*% solve(t(LM) %*% UMi %*% LM))
  }
  if (type == "anderson"){
    G <- t(LM) %*% UMi %*% LM
    scores <- as.data.frame(X %*% UMi %*% LM %*% sqrtm(solve((G %*% (m + G)))))
  }
  
  # Return
  return(scores)
}



facSMC <- function(R, LM){
  ##############################################################################
  # Calculate squared multiple correlations for predicting the common factors
  # R  > (regularized) correlation matrix
  # LM > (rotated) loadings matrix
  #
  # NOTES:
  # - Calculates squared multiple correlations between the observed features
  #   and the common factors
  # - The closer to unity, the better one is able to determine the factors scores
  #   from the observed features
  # - Hence, the closer to unity, the lesser the problem of factor-score 
  #   indeterminacy and the better one is able to uniquely determine the 
  #   factor scores
  # - Note that the computations assume an orthogonal factor model
  # - Hence, only orthogonal rotations of the loadings matrix should be used 
  #   (or no rotation at all)
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  
  # Checks
  if (!is.matrix(R)){
    stop("Input (R) should be a matrix")
  }
  if (nrow(R) != ncol(R)){
    stop("Input (R) should be square matrix")
  }
  if (!is.matrix(LM)){
    stop("Input (LM) should be a matrix")    
  } 
  if (ncol(R) != nrow(LM)){
    stop("Column-dimension input (R) incompatible with row-dimension input (LM)")    
  }
  
  # Determine squared multiple correlations for predicting factors
  return(diag(t(LM) %*% solve(R) %*% LM))
}




##------------------------------------------------------------------------------
##
## Function for Automatic Full Monty
##
##------------------------------------------------------------------------------

autoFMradio <- function(X, t = .95, fold = 5, GB = 1, type = "thomson",
                        verbose = TRUE, printInfo = TRUE, seed = NULL){
  ##############################################################################
  # Function that automatically performs the 3 main steps of the FMradio
  # workflow with minimal user input
  # X         > data matrix or expression set
  # t         > absolute value for thresholding
  # fold      > cross-validation fold, default gives 5-fold CV
  # GB        > numerical indication of which Guttman bound to use
  # type      > character indicating the type of factor score to calculate
  # verbose   > logical indicating if function should run silently
  # printInfo > logical indicating if (run) information should be printed
  #             on-screen
  # seed      > set seed for random number generator
  #
  # NOTES:
  # - Function performs the full monty, working on the correlation scale
  # - If X is a matrix, the observations are expected to be in the rows
  # - t should satisfy -1 < t < 1
  # - If seed = NULL then the starting seed is determined by drawing a single
  #   integer from the integers 1:9e5
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("Biobase")
  # require("stats")
  # require("expm")
  
  # Checks
  if (!inherits(X, "matrix") & !inherits(X, "ExpressionSet")){
    stop("Input (X) should be either of class 'matrix' or 'ExpressionSet'")
  }
  if (class(X) == "ExpressionSet"){
    X <- t(exprs(X))
  }
  if (class(t) != "numeric"){
    stop("Input (t) is of wrong class")
  }
  if (length(t) != 1){
    stop("Length input (t) must be one")
  }
  if (class(fold) != "numeric" & class(fold) != "integer"){ 
    stop("Input (fold) is of wrong class") 
  }
  if (fold <=  1){ 
    stop("Input (fold) should be at least 2") 
  }
  if (fold > nrow(X)){ 
    stop("Input (fold) cannot exceed the sample size") 
  }
  GB <- as.integer(GB)
  if (length(GB) != 1){
    stop("Length input (GB) must be one")
  }
  if (!(GB %in% as.integer(c(1,2,3)))){
    stop("Input (GB) must be either 1, 2, or 3")
  }
  if (!(type %in% c("thomson", "bartlett", "anderson"))){
    stop("Input (type) should be one of {'thomson', 'bartlett', 'anderson'}")
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  if (class(printInfo) != "logical"){ 
    stop("Input (printInfo) is of wrong class") 
  }
  
  # Preliminaries
  Rraw <- cor(X)
  if (is.null(seed)){
    seed <- sample(1:9e5, 1)
  }
  set.seed(seed)
  
  # Redundancy filtering
  if (verbose){cat("Step 1.1: Performing redundancy filtering...", "\n")}
  Rfilter <- RF(Rraw, t = t)
  Xfilter <- subSet(X, Rfilter)
  
  # Regularized correlation matrix
  if (verbose){cat("Step 1.2: Calculating regularized correlation matrix...", "\n")}
  Rreg <- regcor(Xfilter, fold = fold, verbose = FALSE)
  
  # Assessing latent dimensionality
  if (verbose){cat("Step 2.1: Assessing latent dimensionality...", "\n")}
  m <- dimGB(Rreg$optCor, graph = FALSE, verbose = FALSE)[GB]
  
  # Maximum likelihood factor analysis
  if (verbose){cat("Step 2.2: Performing ML factor analysis...", "\n")}
  FA <- mlFA(Rreg$optCor, m = m)
  
  # Computing factor scores
  if (verbose){cat("Step 3: Computing factor scores...", "\n")}
  Scores <- facScore(Xfilter, FA$Loadings, FA$Uniqueness, type = type)
  
  # Computing additional information
  p      <- dim(Xfilter)[2]
  SMC    <- facSMC(Rreg$optCor, FA$Loadings)
  SS     <- diag(t(FA$Loadings)%*%(FA$Loadings))
  PV     <- SS/p
  CV     <- cumsum(PV)
  
  # Printing information
  if (printInfo){
    KMO <- SA(Rreg$optCor)$KMO
    if (GB == as.integer(1)){GBc <- "first"}
    if (GB == as.integer(2)){GBc <- "second"}
    if (GB == as.integer(3)){GBc <- "third"}
    cat("\n")
    cat("\n")
    cat("Step 1: \n")
    cat(paste("Redundancy filtering at threshold value:", t, "\n"))
    cat(paste("      Features retained after filtering:", p, "\n"))
    cat(paste("    Number of folds in cross-validation:", fold, "\n"))
    cat(paste("        Optimal value penalty parameter:", Rreg$optPen, "\n"))
    cat(paste("                              KMO index:", KMO, "\n"))
    cat("Step 2: \n")
    cat(paste("Number of latent factors determined by:", GBc, "Guttman bound", "\n"))
    cat(paste("              Number of latent factors:", m, "\n"))
    cat(paste("      Proportion of explained variance:", CV[m], "\n"))
    cat("Step 3: \n")
    cat(paste("    Type of factor score returned:", type, "\n"))
    cat(paste("Minimum determinacy factor scores:", min(SMC), "\n"))
  }
  
  # Return
  return(list(Scores = Scores,
              FilteredData = Xfilter,
              FilteredCor = Rfilter,
              optPen = Rreg$optPen,
              optCor = Rreg$optCor,
              m = m,
              Loadings = FA$Loadings, 
              Uniqueness = FA$Uniqueness,
              Exvariance = CV,
              determinacy = SMC,
              used.seed = seed))
}




##------------------------------------------------------------------------------
##
## Function for Generating Data under the Common Factor Model
##
##------------------------------------------------------------------------------

FAsim <- function(p, m, n, simplestructure = TRUE, balanced = TRUE,
                  loadingfix = TRUE, loadingnegative = TRUE,
                  loadingvalue = .8, loadingvaluelow = .2, numloadings,
                  loadinglowerH = .7, loadingupperH = .9, 
                  loadinglowerL = .1, loadingupperL = .3){
  ##############################################################################
  # Simulate data according to the factor analytic model
  # - p               > feature dimension
  # - m               > dimension of latent vector
  # - n               > number of samples
  # - simplestructure > logical indicating if factor structure should be 
  #                     factorially pure
  # - balanced        > logical indicating if the 'significant' loadings
  #                     should be divided evenly over the respective factors
  # - loadingfix      > logical indicating if the loadings should have a 
  #                     fixed value
  # - loadingnegative > logical indicating if, next to positive, also negative
  #                     loadings should be present
  # - loadingvalue    > value for high loadings, used when loadingfix = TRUE
  # - loadingvaluelow > value for low loadings, used when loadingfix = TRUE &
  #                     simplestructure = FALSE
  # - numloadings     > vector with length equalling argument m, indicating the
  #                     number of 'significant' loadings per factor. Used when
  #                     balanced = FALSE
  # - loadinglowerH   > lower-bound of 'significant' (high) loadings, used when 
  #                     loadingfix = FALSE
  # - loadingupperH   > upper-bound of 'significant' (high) loadings, used when 
  #                     loadingfix = FALSE
  # - loadinglowerL   > lower-bound of 'non-significant' (low) loadings, used 
  #                     when loadingfix = FALSE & simplestructure = FALSE
  # - loadingupperL   > upper-bound of 'non-significant' (low) loadings, used 
  #                     when loadingfix = FALSE & simplestructure = FALSE
  #
  # NOTES:
  # - Produces a standardized data matrix of size n x p
  # - Also output the correlation matrix based on the generated data and the
  #   loadings matrix and uniquenesses vector on which the data-generation was 
  #   based
  # - A uniform distribution is assumed when generating draws between 
  #   loadinglowerH and loadingupperH
  # - A uniform distribution is assumed when generating draws between 
  #   loadinglowerL and loadingupperL
  ##############################################################################  
  
  # Dependencies:
  # require("base")
  # require("stats")
  # require("MASS")
  
  # Checks
  if (class(p) != "numeric" & class(p) != "integer"){ 
    stop("Input (p) is of wrong class") 
  }
  if (p <=  1){ 
    stop("Input (p) should be at least 2") 
  }
  if (length(p) != 1){
    stop("Length input (p) must be one")
  }
  if (class(m) != "numeric" & class(m) != "integer"){ 
    stop("Input (m) is of wrong class") 
  }
  if (m <=  0){ 
    stop("Input (m) should be strictly positive") 
  }
  if (length(m) != 1){
    stop("Length input (m) must be one")
  }
  mmax <- floor((2*p + 1 - sqrt(8*p + 1))/2)
  if (m > mmax){
    stop("Input (m) is too high")
  }
  if (class(n) != "numeric" & class(n) != "integer"){ 
    stop("Input (n) is of wrong class") 
  }
  if (n <=  1){ 
    stop("Input (n) should be at least 2") 
  }
  if (length(n) != 1){
    stop("Length input (n) must be one")
  }
  if (class(simplestructure) != "logical"){ 
    stop("Input (simplestructure) is of wrong class") 
  }
  if (class(balanced) != "logical"){ 
    stop("Input (balanced) is of wrong class") 
  }
  if (class(loadingfix) != "logical"){ 
    stop("Input (loadingfix) is of wrong class") 
  }
  if (class(loadingnegative) != "logical"){ 
    stop("Input (loadingnegative) is of wrong class") 
  }
  
  # Conditional checks
  if (loadingfix){
    if (class(loadingvalue) != "numeric"){ 
      stop("Input (loadingvalue) is of wrong class") 
    }
    if (length(loadingvalue) != 1){ 
      stop("Length input (loadingvalue) must be one") 
    }
    if (class(loadingvaluelow) != "numeric"){ 
      stop("Input (loadingvaluelow) is of wrong class") 
    }
    if (length(loadingvaluelow) != 1){ 
      stop("Length input (loadingvaluelow) must be one") 
    }
  }
  if (!loadingfix){
    if (class(loadinglowerH) != "numeric"){ 
      stop("Input (loadinglowerH) is of wrong class") 
    }
    if (length(loadinglowerH) != 1){ 
      stop("Length input (loadinglowerH) must be one") 
    }
    if (class(loadingupperH) != "numeric"){ 
      stop("Input (loadingupperH) is of wrong class") 
    }
    if (length(loadingupperH) != 1){ 
      stop("Length input (loadingupperH) must be one") 
    }
    if (!simplestructure){
      if (class(loadinglowerL) != "numeric"){ 
        stop("Input (loadinglowerL) is of wrong class") 
      }
      if (length(loadinglowerL) != 1){ 
        stop("Length input (loadinglowerL) must be one") 
      }
      if (class(loadingupperL) != "numeric"){ 
        stop("Input (loadingupperL) is of wrong class") 
      }
      if (length(loadingupperL) != 1){ 
        stop("Length input (loadingupperL) must be one") 
      }
    }
  }
  if (!balanced){
    if (class(numloadings) != "numeric" & class(numloadings) != "integer"){
      stop("Input (numloadings) is of wrong class") 
    }
    if (length(numloadings) != m){
      stop("Length of input (numloadings) should equal argument m") 
    }
    if (sum(numloadings) != p){
      stop("Sum of inputs (numloadings) should equal argument p") 
    }
  }
  
  # Loadings matrix
  Lambda <- matrix(0,p,m)
  if (balanced){
    h   <- floor(p/m)
    hi1 <- 1 - h
    hi2 <- 0
    for (i in 1:m){
      hi1 <- hi1 + h
      hi2 <- hi2 + h
      if (loadingfix){
        Lambda[(hi1):(hi2),i] <- loadingvalue
      }
      if (!loadingfix){
        Lambda[(hi1):(hi2),i] <- 
          runif(length((hi1):(hi2)), 
                min = loadinglowerH, 
                max = loadingupperH)
      }
      if (loadingnegative){
        Polarity <- sign(runif(length((hi1):(hi2)), min = -1, max = 1))
        Lambda[(hi1):(hi2),i] <- Lambda[(hi1):(hi2),i] * Polarity
      }
    }
    if (hi2 != p){
      if (loadingfix){
        Lambda[((hi2)+1):p,i] <- loadingvalue
      }
      if (!loadingfix){
        Lambda[((hi2)+1):p,i] <-
          runif(length(((hi2)+1):p), 
                min = loadinglowerH, 
                max = loadingupperH)
      }
    }
    if (!simplestructure){
      if (loadingfix){
        if (loadingnegative){
          Polarity <- sign(runif(length(Lambda[Lambda == 0]), min = -1, max = 1))
          loadingvaluelow <- loadingvaluelow * Polarity
        }
        Lambda[Lambda == 0] <- loadingvaluelow
      }
      if (!loadingfix){
        nbs <- runif(length(Lambda[Lambda == 0]), 
                     min = loadinglowerL, 
                     max = loadingupperL)
        if (loadingnegative){
          Polarity <- sign(runif(length(Lambda[Lambda == 0]), min = -1, max = 1))
          nbs      <- nbs * Polarity
        }
        Lambda[Lambda == 0] <- nbs
      }
    }
  }
  if (!balanced){
    low  <- cumsum(numloadings)
    low  <- c(1,low[-length(low)] + 1)
    high <- cumsum(numloadings)
    for (i in 1:length(numloadings)){
      if (loadingfix){
        Lambda[low[i]:high[i],i] <- loadingvalue
      }
      if (!loadingfix){
        Lambda[low[i]:high[i],i] <- 
          runif(length(low[i]:high[i]), 
                min = loadinglowerH, 
                max = loadingupperH)
      }
      if (loadingnegative){
        Polarity <- sign(runif(length(low[i]:high[i]), min = -1, max = 1))
        Lambda[low[i]:high[i],i] <- Lambda[low[i]:high[i],i] * Polarity
      }
    }
    if (!simplestructure){
      if (loadingfix){
        if (loadingnegative){
          Polarity <- sign(runif(length(Lambda[Lambda == 0]), min = -1, max = 1))
          loadingvaluelow <- loadingvaluelow * Polarity
        }
        Lambda[Lambda == 0] <- loadingvaluelow
      }
      if (!loadingfix){
        nbs <- runif(length(Lambda[Lambda == 0]), 
                     min = loadinglowerL, 
                     max = loadingupperL)
        if (loadingnegative){
          Polarity <- sign(runif(length(Lambda[Lambda == 0]), min = -1, max = 1))
          nbs      <- nbs * Polarity
        }
        Lambda[Lambda == 0] <- nbs
      }
    }
  }
  
  # Uniqueness matrix
  Psi <- diag(p) - diag(diag(Lambda %*% t(Lambda)))
  
  # Generate data
  SigmaModel <- Lambda %*% t(Lambda) + Psi
  Scaling    <- diag(sqrt(1/diag(SigmaModel)))
  SigmaModel <- Scaling %*% SigmaModel %*% Scaling
  Data       <- mvrnorm(n = n, mu = rep(0,p), Sigma = SigmaModel, 
                        tol = 1e-6, empirical = FALSE, EISPACK = FALSE)
  
  # Automatic naming of features
  prefix         <- "feature"
  suffix         <- seq(1:p)
  colnames(Data) <- paste(prefix, suffix, sep = ".")
  
  # Return
  data <- scale(Data)
  return(list(data = data, loadings = Lambda, 
              Uniqueness = diag(Psi), cormatrix = cor(data)))
}




##------------------------------------------------------------------------------
## 
## Hidden Gems
##
##------------------------------------------------------------------------------

.Airwolf <- function() {
  ##############################################################################
  # - Bet you did not know that
  ##############################################################################
  
  cat("
      ###################################################
      ###################################################
      Airwolf's scanner was powered by FMradio!
      - Stringfellow Hawke
      ###################################################
      ################################################### \n")
}



.Airwolf2 <- function() {
  ##############################################################################
  # - Truth
  ##############################################################################
  
  cat("
      ###################################################
      ###################################################
      My tune is better than KITT's!
      - Stringfellow Hawke
      ###################################################
      ################################################### \n")
}

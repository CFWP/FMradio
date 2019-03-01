\name{dimIC}
\alias{dimIC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dimIC(x)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


dimIC <- function(R, n, maxdim, Type = "BIC", graph = TRUE, verbose = TRUE){
  ##############################################################################
  # Performs dimensionality assessment by way of the BIC
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
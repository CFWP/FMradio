\name{scree}
\alias{scree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
scree(x)
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

scree <- function(R, m, type = "ML", vertical = TRUE, verbose = TRUE){
  ##############################################################################
  # Gives scree plot under an m-factor solution
  # R        > (regularized) covariance or correlation matrix
  # m        > dimension of factor solution
  # type     > character indicating which eigenvalues to calculate
  # vertical > logical indicating if value of argument m should be indicated
  #            in graph by a vertical line
  # verbose  > logical indicating if function should run silently
  #
  # NOTES:
  # - Argument 'm' only used under type = ML and type = reduced
  # - If type = ML, the eigenvalues under the m-factor solution are returned
  # - If type = reduced, the eigenvalues of the reduced correlation matrix
  #   under the m-factor solution are returned
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
  if (!(type %in% c("ML", "reduced"))){
    stop("Input (type) should be one of {'ML', 'reduced'}")
  }
  if (class(vertical) != "logical"){ 
    stop("Input (vertical) is of wrong class") 
  }
  if (class(verbose) != "logical"){ 
    stop("Input (verbose) is of wrong class") 
  }
  
  # Fit m-factor model
  fito <- factanal(factors = m, covmat = R, rotation = "none")
  Psi  <- diag(fito$uniquenesses)
  
  # Calculate eigenvalues
  if (verbose){cat("Calculating eigenvalues...", "\n")}
  if (type == "ML"){
    PsiI <- diag(1/sqrt(diag(Psi)))
    Te   <- PsiI %*% R %*% PsiI
    EVs  <- eigen(Te)$values
    yLab <- "Eigenvalue ML factor solution"
  }
  if (type == "reduced"){
    Te  <- R - Psi
    EVs <- eigen(Te)$values
    yLab <- "Eigenvalue reduced correlation matrix"
  }
  
  # Visualize
  if (verbose){cat("Visualizing...", "\n")}
  dims <- c(1:dim(R)[1])
  plot(dims, EVs, axes = FALSE, type = "l", 
       col = "red", xlab = "eigenvalue no.", 
       ylab = yLab, main = paste("Dimension of factor solution:", m))
  axis(2, ylim = c(0, max(EVs)), col = "black", lwd = 1)
  axis(1, xlim = c(1,max(dims)), col = "black", lwd = 1, tick = TRUE)
  abline(v = m)
  
  # Return object
  EVt <- as.data.frame(cbind(dims, EVs))
  colnames(EVt) <- c("latent dimension", "Eigenvalue")
  
  # Return
  return(EVt)
}
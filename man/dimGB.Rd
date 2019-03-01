\name{dimGB}
\alias{dimGB}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dimGB(x)
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


dimGB <- function(R, graph = TRUE, verbose = TRUE){
  ##############################################################################
  # Calculates the first, second, and third Guttman lower-bounds
  # R        > (regularized) covariance or correlation matrix
  # verbose  > logical indicating if function should run silently
  #
  # NOTES:
  # - First lower-bound corresponds to kaiser's rule
  # - What are the lower-bounds
  # - There is ordering in bounds
  # - 
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

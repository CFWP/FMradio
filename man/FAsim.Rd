\name{FAsim}
\alias{FAsim}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FAsim(x)
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
  # - loadinglowerH   > lower-bound of 'significant' (high) loading, used when 
  #                     loadingfix = FALSE
  # - loadingupperH   > upper-bound of 'significant' (high) loading, used when 
  #                     loadingfix = FALSE
  # - loadinglowerL   > lower-bound of 'non-significant' (low) loading, used 
  #                     when loadingfix = FALSE & simplestructure = FALSE
  # - loadingupperL   > upper-bound of 'non-significant' (low) loading, used 
  #                     when loadingfix = FALSE & simplestructure = FALSE
  #
  # NOTES:
  # - Produces a standardized data matrix of size n x p
  # - Also output the correlation matrix based on the generated data and the
  #   loadings and uniquenesses matrices on which the data-generation was based
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

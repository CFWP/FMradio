\name{autoFMradio}
\alias{autoFMradio}
\title{
Wrapper for automated workflow
}
\description{
\code{autoFMradio} is a wrapper function that automates the three main steps of the \code{FMradio} workflow.
}
\usage{
autoFMradio(X, t = .95, fold = 5, GB = 1, type = "thomson",
            verbose = TRUE, printInfo = TRUE, seed = NULL)
}
\arguments{
  \item{X}{
   A data \code{matrix} or an \code{ExpressionSet} object. 
  }
  \item{t}{
  A scalar \code{numeric} indicating the absolute value for thresholding.
  }
  \item{fold}{
  A \code{numeric} integer or \code{integer} indicating the number of folds to use in cross-validation.
  }
  \item{GB}{
  A \code{numeric} integer or \code{integer} indicating which Guttman bound to use for determining the number of latent features to retain.\cr
  Must be either 1, 2, or 3.
  }
  \item{type}{
  A \code{character} indicating the type of factor score to calculate.\cr
  Must be one of: "thomson", "bartlett", "anderson".
  }
  \item{verbose}{
  A \code{logical} indicating if function should run silently.\cr
  Runs silently when \code{verbose = FALSE}.
  }
  \item{printInfo}{
  A \code{logical} indicating if additional information should be printed on-screen.\cr
  Suppresses printing when \code{verbose = FALSE}.
  }
  \item{seed}{
  A \code{numeric} integer or \code{integer} indicating the seed for the random number generator.
  }
}
\details{
The \code{autoFMradio} function automates the three main steps of the workflow by providing a wrapper around all core functions.

Step 1 (regularized correlation matrix estimation) is performed using the \code{X}, \code{t}, and \code{fold} arguments.
The raw correlation matrix based on data \code{X} is redundancy-filtered using the threshold provided in \code{t}.
Subsequently, a regularized estimate of the correlation matrix (on the possibly filtered feature set) is computed with the optimal penalty value determined by cross-validation.
The number of folds is set by the \code{fold} argument.
For more information on Step 1 see \code{\link{RF}}, \code{\link{subSet}}, and \code{\link{regcor}}.

Step 2 (factor analytic data compression) is performed using the \code{GB} argument.
With this argument one can use either the first, second, or third Guttman bound to select the intrinsic dimensionality of the latent vector.
This bound, together with the regularized correlation matrix, is used in a maximum likelihood factor analysis with simple-structure rotation.
For more information on Step 2, see \code{\link{dimGB}} and \code{\link{mlFA}}.

Step 3 (obtaining factor scores) is performed using the \code{type} argument.
It determines factor scores: the score each object/individual would obtain on each of the latent factors.
The \code{type} argument determines the type of factor score that is calculated.
For more information on Step 3, see \code{\link{facScore}}.

When \code{printInfo = TRUE} additional information is printed on-screen after the full procedure has run its course. 
This additional information pertains to each of the steps mentioned above.
For Step 1 it reiterates the thresholding value for redundancy filtering and gives the number of features retained after this filtering.
It also reiterates the number of folds used in determining the optimal penalty value as well as this value itself.
Moreover, it provides the value of the Kaiser-Meyer-Olkin index on the optimal regularized correlation matrix estimate (see \code{\link{SA}}).
For Step 2 it reiterates which Guttman bound was used in determining the number of latent factors as well as the number of latent factors retained.
It also gives the proportion of explained variance under the factor solution of the chosen latent dimension (see \code{\link{dimVAR}}).
For step 3 it reiterates the type of factor score that was calculated.
Also, it prints the lowest `determinacy score' amongst the latent factors (see \code{\link{facSMC}}).

The factor scores in the \code{$Scores} slot of the output (see below) can be directly used as input features in any prediction or classification procedure.
In case of external (rather than internal) validation one can use the parameter matrices in the \code{$Loadings} and \code{$Uniqueness} slots in combination with fresh data to provide a validation factor projection based on the training solution.
See Peeters \emph{et al.} (2019).
}
\value{
The function returns an object of class \code{list}:
\item{$Scores}{An object of class \code{data.frame} containing the factor scores. Observations are represented in the rows. Each column represent a latent factor.}
\item{$FilteredData}{Subsetted data \code{matrix} containing only those features retained after redundancy filtering.}
\item{$FilteredCor}{A correlation \code{matrix} based on the data in the \code{$FilteredData} slot.}
\item{$optPen}{A \code{numeric} scalar representing the optimal value for the penalty parameter.}
\item{$optCor}{A \code{matrix} representing the regularized correlation matrix under the optimal penalty-value.}
\item{$m}{An \code{integer} correspond to number of latent factors retained under the chosen Guttman bound.}
\item{$Loadings}{A matrix of class \code{loadings} representing the loadings matrix in which in which each element \eqn{\lambda_{jk}} is the loading of the \eqn{j}th feature on the \eqn{k}th latent factor.}
\item{$Uniqueness}{A \code{matrix} representing the diagonal matrix carrying the unique variances.}
\item{$Exvariance}{A \code{numeric} vector representing the cumulative variance for each respective latent feature.}
\item{$determinacy}{A \code{numeric} vector indicating, for each factor, the squared multiple correlation between the observed features and the common latent factor.}
\item{$used.seed}{A \code{numeric} or \code{integer} used as the starting seed in random number generation.}
}
\references{
Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.
}
\author{
Carel F.W. Peeters <carel.peeters@wur.nl>
}
\note{
When \code{seed = NULL} the starting seed is determined by drawing a single integer from the integers \code{1:9e5}. This non-user-supplied seed is also found in the \code{$used.seed} slot of the output.
}
\seealso{
\code{\link{RF}}, \code{\link{subSet}}, \code{\link{regcor}}, \code{\link{dimGB}}, \code{\link{mlFA}}, \code{\link{facScore}}
}
\examples{
## Simulate some data according to a factor model with 3 latent factors
simDAT <- FAsim(p = 24, m = 3, n = 40, loadingvalue = .9)
X <- simDAT$data

## Perform the lot
FullMonty <- autoFMradio(X, GB = 1, seed = 303)
}

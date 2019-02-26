\name{radioHeat}
\alias{radioHeat}
\title{
Visualize a (correlation) matrix as a heatmap
}
\description{
Function that provides dedicated heatmapping of a radiomics-based correlation matrix
Can be used to visually assess the elements of a (possibly thresholded) matrix.
It supports the assessment of collinearity.
}
\usage{
radioHeat(R, lowColor = "blue", highColor = "red", labelsize = 10,
          diag = TRUE, threshold = FALSE, threshvalue = .95,
          values = FALSE, textsize = 10, legend = TRUE, main = "")
}
\arguments{
  \item{R}{
  (regularized) correlation \code{matrix}
  }
  \item{lowColor}{
  A \code{character} that determines the color scale in the negative range.
  }
  \item{highColor}{
  A \code{character} that determines the color scale in the positive range.
  }
  \item{labelsize}{
  A \code{numeric} that sets the textsize of row and column labels.
  }
  \item{diag}{
  A \code{logical} determining if the diagonal elements of the matrix should be included in the
  color scaling. 
  This argument is only used when \code{R} is a square \code{matrix}.
  }
  \item{threshold}{
  A \code{logical} determining if only values above a certain (absolute) threshold should be
  visualized. 
  }
  \item{threshvalue}{
  A \code{numeric} indicating the absolute thresholding value when \code{threshold = TRUE}.
  }
  \item{values}{
  A \code{logical} determining the optional inclusion of cell-values.
  }
  \item{textsize}{
  A \code{numeric} indicating the textsize of the cell-values when \code{values = TRUE}.
  }
  \item{legend}{
  A \code{logical} indicating whether a color legend should be included.
  }
  \item{main}{
  A \code{character} giving the main figure title.
  }
}
\details{
This function utilizes \href{https://cran.r-project.org/package=ggplot2}{ggplot2} (Wickham, 2009) to visualize a matrix as a heatmap: a false color plot in which the individual matrix entries are represented by colors. 
\code{lowColor} determines the color scale for matrix entries in the negative range.
\code{highColor} determines the color scale for matrix entries in the positive range. 
For the colors supported by the arguments \code{lowColor} and \code{highColor}, see \url{www.stat.columbia.edu/~tzheng/files/Rcolor.pdf}. 
White entries in the plot represent the midscale value of 0. 
One can opt to set the diagonal entries to the midscale color of white when one is interested in (heatmapping) the off-diagonal elements only. 
To achieve this, set \code{diag = FALSE}. 
Naturally, the \code{diag} argument is only used when the input matrix \code{M} is a square matrix.

The intended usage is to visualize a correlation matrix on radiomic features as a heatmap.
Such a heatmap may be used to support the assessment of strong collinearity or even redundancy amongst the features.
To this end, it is also possible to visualize a thresholded correlation matrix when \code{threshold = TRUE} based on the absolute thresholding value given in the \code{threshvalue} argument (hence the thresholding is done internally).
This enables easier visual access to (blocks of) collinearity in radiomic-feature-based correlation matrices. 
}
\references{
Wickham, H. (2009). ggplot2: elegant graphics for data analysis. New York: Springer.
}
\author{
Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\note{
\itemize{
\item While geared towards the visualization of correlation matrices, the function is quite general, in the sense that it can represent any \code{matrix} as a heatmap.
\item When \code{values = TRUE} and \code{threshold = TRUE} the cell-values are those of the thresholded matrix.
}
}
\seealso{
\code{\link{RF}}, \code{\link{regcor}}
}
\examples{
## Generate some (high-dimensional) data
p = 25
n = 10
set.seed(333)
X = matrix(rnorm(n*p), nrow = n, ncol = p)
colnames(X)[1:25] = letters[1:25]
R <- cor(X)

## Visualize the correlation matrix as a heatmap
radioHeat(R)

## Remove diagonal entries from visualization
radioHeat(R, diag = FALSE)

## Additionally, visualize only those entries whose absolute value exceed .5
radioHeat(R, diag = FALSE, threshold = TRUE, threshvalue = .5)

## Additionally, include cell values
radioHeat(R, diag = FALSE, threshold = TRUE, threshvalue = .5, 
          values = TRUE, textsize = 3)
}
\concept{radiomic-feature heatmap}

\name{FMradio-package}
\alias{FMradio-package}
\alias{FMradio}
\docType{package}
\title{
Factor modeling for radiomic data
}
\description{
The \code{FMradio} package provides a workflow that uses factor modeling to project the high-dimensional and collinear radiomic feature-space onto a lower-dimensional orthogonal meta-feature space that retains most of the information contained in the full data set. 
These projected meta-features can be directly used as robust and stable covariates in any downstream
prediction or classification model.
}
\details{
Radiomics refers to the mining of large numbers of quantitative features from standard-of-care
clinical images.
\code{FMradio} aims to provide support for stable prediction and classification modeling with radiomics data, irrespective of imaging modality (such as MRI, PET, or CT).
The workflow has 3 main steps that ultimately enable stable prediction and classification.

\bold{Step 1: Regularized correlation matrix estimation.}
Radiomic data are often high-dimensional in the sense that there are more features than observations.
Moroever, radiomic data are often highly collinear, in the sense that collections of features may be highly correlated (in the absolute sense).
This results in the correlation matrix on the radiomic features to be ill-conditioned or even singular.
It is also this combination of characteristics that proves difficult to predictive modeling.
As the factor-analytic procedure is based on the modeling of moment structures such as the correlation matrix, the first step is to obtain a regularized, well-conditioned estimate of the correlation matrix.
The following functions are then of use:
\itemize{
  \item \code{\link{radioHeat}}
  \item \code{\link{RF}}
  \item \code{\link{subSet}}
  \item \code{\link{regcor}}
}
The \code{radioHeat} function can be used to visualize (a possibly regularized) correlation matrix as a heatmap.
It can also be used to visually assess feature-redundancy.
The \code{RF} function provides functionality for filtering features that are so collinear that they are deemed redundant.
The \code{suBSet} function provides functionality to subset data objects to those features retained after possible filtering.
The \code{regcor} function subsequently provides a regularized estimate of the correlation matrix (on the possibly filtered feature set).

\bold{Step 2: Factor analytic data compression.}
The next step would be to project the collinear and high-dimensional radiomic feature-space onto a lower-dimensional orthogonal meta-feature space.
Factor analysis can be used for this purpose.
The following functions are then of use:
\itemize{
  \item \code{\link{SA}}
  \item \code{\link{dimGB}}
  \item \code{\link{dimVAR}}
  \item \code{\link{SMC}}
  \item \code{\link{mlFA}}
}
The \code{SA} function assesses if performing a factor analysis on the (possibly regularized) correlation matrix would be appropriate.
The \code{dimGB} function can be used to determine the number of latent factors (i.e., to determine the intrinsic dimensionality of the meta-feature space).
The \code{dimVAR} and \code{dimSMC} functions can be used to provide additional decision support with respect to the output of the \code{dimGB} function.
The \code{mlFA} function then performs a maximum likelihood factor analysis using the (possibly regularized) correlation matrix and the choice of intrinsic dimensionality as inputs.

\bold{Step 3: Obtaining factor scores.}
The third step would be to use the factor analytic solution to obtain factor scores: the score each object/individual would obtain on each of the latent factors.
The following functions are then of use:
\itemize{
  \item \code{\link{facScore}}
  \item \code{\link{facSMC}}
}
The \code{facScore} function provides several options for computing factors scores.
The determinacy of these scores can be assessed with the \code{facSMC} function.

\bold{Step 4: Prediction and classification.}
The factor scores obtained with Step 3 can be directly used as (low-dimensional and orthogonal) covariates in any prediction, classification or learning procedure.
One may use the full flexibility provided by the CRAN repository for this step.
\cr
\cr
\emph{Additional functionality.}
The package also provides additional functionality.
These are contained in the following (convenience) functions:
\itemize{
  \item \code{\link{dimLRT}}
  \item \code{\link{dimIC}}
  \item \code{\link{FAsim}}
}
The \code{dimLRT} and \code{dimIC} functions provide alternative options for assessing the number of latent factors using likelihood ratio testing and information criteria, respectively.
These are only recommended when the sample size is large relative to the number of features.
\code{FAsim} provides a flexible function for generating data according to the orthogonal common factor analytic model. 
All these functions may be of use in comparative exercises.
The package also provides a wrapper function that automates the 3 main steps of the workflow:
\itemize{
  \item \code{\link{autoFMradio}}
}
}
\author{
Carel F.W. Peeters [cre, aut]\cr
Caroline Ubelhor [ctb]
\cr
\cr
\emph{Maintainer}: Carel F.W. Peeters <cf.peeters@vumc.nl>
}
\references{
Peeters, C.F.W. \emph{et al.} (2019). Stable prediction with radiomics data.
\href{https://arxiv.org/abs/1903.11696}{arXiv:1903.11696 [stat.ML]}.
}

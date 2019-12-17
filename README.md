![FMradio](https://github.com/CFWP/FMradio/blob/master/inst/FMradioLOGO.png)


[![Project Status: Active - The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![CRAN version](http://www.r-pkg.org/badges/version/FMradio)](https://cran.r-project.org/package=FMradio)
[![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/FMradio)](https://cran.r-project.org/package=FMradio/index.html)
[![Total CRAN downloads](http://cranlogs.r-pkg.org/badges/grand-total/FMradio)](http://www.r-pkg.org/pkg/FMradio)


**FMradio**
---------------
The R-package **FMradio** supports stable prediction and classification with radiomics data through factor-analytic modeling.
This support can be invoked irrespective of the imaging modality (such as, e.g., MRI, PET, CT) used to produce the radiomics data.


## Installation

If you wish to install the latest version of **FMradio** directly from the master branch here at GitHub, run

```R
#install.packages("devtools")  # Uncomment if devtools is not installed
devtools::install_github("CFWP/FMradio")
```
Be sure that you have the [package development prerequisites](http://www.rstudio.com/ide/docs/packages/prerequisites) if you wish to install the package from the source.


## Manual and other documentation

* A pdf-version of the manual can be found [here](https://cfwp.github.io/PDFs/FMradio.pdf).

* The R-script used to produce the simulations and results contained in Peeters, C.F.W., *et al*. (2019) [see references below] can be found [here](https://cfwp.github.io/PDFs/Simulations&Analysis.R).


## References

Relevant publications to **FMradio** include:

 1. Peeters, C.F.W. (2019). 
    *"FMradio: Factor modeling for radiomic data"*. 
    R package, version 1.1.1
 2. Peeters, C.F.W., *et al*. (2019)
    *"Stable prediction with radiomics data"*.
    [arXiv:1903.11696 [stat.ML]](https://arxiv.org/abs/1903.11696)

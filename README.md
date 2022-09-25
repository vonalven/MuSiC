MuSiC
=============================================
`MuSiC` is an analysis toolkit for single-cell RNA-Seq experiments. To use this package, you will need the R statistical computing environment (version 3.0 or later) and several packages available through Bioconductor and CRAN.

The original release of `MuSiC` is a deconvolution method that utilizes cross-subject scRNA-seq to estimate cell type proportions in bulk RNA-seq data.
![MuSiC\_pipeline](FigureMethod.jpg)

`MuSiC2` is an iterative algorithm aiming to improve cell type deconvolution for bulk RNA-seq data using scRNA-seq data as reference when the bulk data are generated from samples with multiple clinical conditions where at least one condition is different from the scRNA-seq reference.
![MuSiC\_music2](MuSiC2.jpg)


How to cite `MuSiC`
-------------------
Please cite the following publications:

> *Bulk tissue cell type deconvolution with multi-subject single-cell expression reference*<br />
> <small>X. Wang, J. Park, K. Susztak, N.R. Zhang, M. Li<br /></small>
> Nature Communications. 2019 Jan 22 [https://doi.org/10.1038/s41467-018-08023-x](https://doi.org/10.1038/s41467-018-08023-x) 

> *MuSiC2: cell type deconvolution for multi-condition bulk RNA-seq data*<br />
> <small>J. Fan, Y. Lyu, Q. Zhang, X. Wang, R. Xiao, M. Li<br /></small>
> biorXiv. 2022 [https://www.biorxiv.org/content/10.1101/2022.05.08.491077v1](doi: https://doi.org/10.1101/2022.05.08.491077)


Installation
------------
Both `MuSiC` and `MuSiC2` functions are available in one package.
``` r
# install devtools if necessary
install.packages('devtools')

# install the MuSiC package
devtools::install_github('xuranw/MuSiC')

# load
library(MuSiC)
```

More Information
-----------------
Please see Tutorials for [MuSiC](http://xuranw.github.io/MuSiC/articles/MuSiC.html) and [MuSiC2](http://xuranw.github.io/MuSiC/articles/pages/MuSiC2.html).

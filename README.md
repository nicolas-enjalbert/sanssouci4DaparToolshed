# sanssouci4DaparToolshed

<!-- badges: start -->
[![R-CMD-check](https://github.com/nicolas-enjalbert/sanssouci4DaparToolshed/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/nicolas-enjalbert/sanssouci4DaparToolshed/actions)
[![coverage](https://github.com/sanssouci-org/sanssouci4DaparToolshed/actions/workflows/test_covr.yaml/badge.svg)](https://github.com/sanssouci-org/sanssouci4DaparToolshed/actions/workflows/test_covr.yaml)
<!-- badges: end -->

## Overview

**sanssouci4DaparToolshed** provides helper functions to seamlessly 
integrate the multiple testing framework implemented in the 
[**sanssouci**](https://github.com/sanssouci-org/sanssouci) 
package into the [**DaparToolshed**](https://github.com/edyp-lab/DaparToolshed) 
workflow.

In this workflow, data preparation (normalization, imputation, etc.) is performed
on the peptide dataset. Once the data has been aggregated to form the protein dataset,
differential analysis is performed using a *post hoc* method. This type of method
provides statistical guarantees regarding the *false discovery rate* for
simultaneous, data-driven selections on a volcano plot. 

## Installation

You can install the development version from GitHub:

```r
# install.packages("pak")
pak::pak("nicolas-enjalbert/sanssouci4DaparToolshed")
```

or

```r
# install.packages("remotes")
remotes::install_github("nicolas-enjalbert/sanssouci4DaparToolshed")
```

## Dependencies

The package relies on:

- DaparToolshed (>= 0.99.37),
- ggplot2,
- plotly,
- limma,
- S4Vectors,
- sanssouci (>= 0.16.3),
- SummarizedExperiment

Please use R >= 4.5.0.

## Citation

If you use this package in your work, please follow `sanssouci` and 
`DaparToolshed` recommendation. You can use 

```r
citation("sanssouci")
citation("DaparToolshed")
citation("sanssouci4DaparToolshed")
```

## License

GPL (>= 3)

## Bug reports

Please report issues or request new features through the GitHub issue tracker.


## Example 

We illustrate the package with on 
[vignette](vignettes/sanssouci4DaparToolshed.Rmd). The analyse is made on 
single-cell proteomic data. The same workflow can be used on bulk proteomic 
data.
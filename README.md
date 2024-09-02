# About
Possibilistic meta-analysis with finite-sample validity guarantees

This is a repository for the R package to perform possibilistic meta-analysis in the inferential models framework. The R package's main files are:

- src/compute_funcs-20240828.cpp: this file defines the Rcpp functions that compute the possibility contours for the individual, optimal, meta-analysis and equilibrated possibility contours.
- R/meta_analysis_funcs-20240902.R: this file defines the function for the meta-analysis possibility contour and helper functions.

The IMeta man file contains a examples for running the regression models from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The IMeta R package can be installed in one of two ways:

- from the downloaded gzipped tarball as R CMD INSTALL IMeta_1.0-1.tar.gz

- from the downloaded and renamed IMeta folder as R CMD build IMeta and R CMD INSTALL IMeta_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the IMeta R package, please consider citing the relevant manuscript: E.C. Hector, L. Cella and R. Martin (2024+). Efficient possibilistic meta-analysis with finite-sample validity guarantees.

# References

G.V. Glass. (1976). Primary, secondary, and meta-analysis of research. Educational Researcher, 5(10):3–8.

R. Martin and C. Liu. (2015). Inferential models: reasoning with uncertainty. CRC Press.

# About
Possibilistic divide-and-conquer analysis with finite-sample validity guarantees

This is a repository for the R package to perform possibilistic divide-and-conquer analysis in the inferential models framework. The R package's main files are:

- src/compute_funcs.cpp: this file defines the Rcpp functions that compute the possibility contours for the individual, optimal, divide-and-conquer analysis and equilibrated possibility contours.
- R/meta_analysis_funcs.R: this file defines the function for the divide-and-conquer analysis possibility contour and helper functions.

The IMdac man file contains a examples for running the models from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The IMdac R package can be installed in one of two ways:

- from the downloaded gzipped tarball as R CMD INSTALL IMdac_1.0-1.tar.gz

- from the downloaded and renamed IMdac folder as R CMD build IMdac and R CMD INSTALL IMdac_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the IMdac R package, please consider citing the relevant manuscript: E.C. Hector, L. Cella and R. Martin (2025+). 

# References

G.V. Glass. (1976). Primary, secondary, and meta-analysis of research. Educational Researcher, 5(10):3–8.

E.C. Hector, L. Tang, L. Zhou, and P.X.-K. Song (2024). Handbook on Bayesian, Fiducial and Frequentist Inference, Chapter "Data integration and fusion in the Bayesian and Frequentist frameworks." Chapman and Hall/CRC Press.

R. Martin and C. Liu. (2015). Inferential models: reasoning with uncertainty. CRC Press.

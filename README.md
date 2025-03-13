# About
Divide-and-conquer with finite sample sizes: valid and efficient possibilistic inference

This is a repository for the R package to perform possibilistic divide-and-conquer analysis in the inferential models framework. The R package's main file is the R/dac_funcs-20250221.R file. It defines the function for the divide-and-conquer possibility contour and helper functions.

The IMdac man file and the examples folder contain examples for running the models from the paper.

Please email ehector@ncsu.edu with any questions or bug-reports.

# Installation

The IMdac R package can be installed in one of two ways:

- from the downloaded gzipped tarball as R CMD INSTALL IMdac_1.0-1.tar.gz

- from the downloaded and renamed IMdac folder as R CMD build IMdac and R CMD INSTALL IMdac_1.0-1.tar.gz

Please make sure to have all packages listed in the DESCRIPTION file already installed. If you encounter a library not found error for lgfortran, please try installing gfortran from here: https://cran.r-project.org/bin/macosx/tools/.

# Citation

If you use the IMdac R package, please consider citing the relevant manuscript: E.C. Hector, L. Cella and R. Martin (2025+). Divide-and-conquer with finite sample sizes: valid and efficient possibilistic inference.

# References

Hector, E. C., Tang, L., Zhou, L., and Song, P. X.-K. (2024). Handbook on Bayesian, Fiducial and Frequentist Inference, chapter Data integration and fusion in the Bayesian and Frequentist frameworks. Chapman and Hall/CRC Press.

R. Martin and C. Liu. (2015). Inferential models: reasoning with uncertainty. CRC Press.

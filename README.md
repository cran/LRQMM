LRQMM : Fitting Linear Quantile Regression Mixed Models  
September 5, 2019
Title: Fitting Linear Quantile Regression Mixed Models
Version: 1.1.0
Author: Sayyed Reza Alavian
Maintainer: Sayyed Reza Alavian <s.rezaalavian@mail.um.ac.ir>
Description: Fit a quantile regression mixed model using a sparse implementation of the 
             Frisch-Newton interior-point algorithm as described in
             Portnoy and Koenker (1977, Statistical Science) <https://www.jstor.org/stable/2246217>.
License: GPL-2 | GPL-3
Encoding: UTF-8
Imports: GeneticsPed, SparseM, MASS, quantreg, Matrix, MasterBayes, MCMCglmm
__________________________________________________________________________________________________________
We used "GeneticsPed", "Matrix", "MasterBayes", "MCMCglmm", "MASS", "SparseM" and "quantreg"
packages in this function. befor using "lrqmm" function be sure from installation this packages.
"GeneticsPed" available in <https://bioconductor.org/packages/release/bioc/src/contrib/GeneticsPed_1.46.0.tar.gz>.
other packages are available in CRAN.
___________________________________________________________________________________________________________
If you cnnot install GeneticsPed, can see <http://bioconductor.org/packages/release/bioc/html/GeneticsPed.html>
or using this code in R befor installing LRQMM:
>if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
>BiocManager::install("GeneticsPed")
___________________________________________________________________________________________________________
After insatllotion can run this example:
>library(LRQMM)
>data(Cow)
>with(lrqmm(id,sire,dam,HERD,LACTLENGHT,2,0.5,Factor=TRUE),data=Cow)
___________________________________________________________________________________________________________
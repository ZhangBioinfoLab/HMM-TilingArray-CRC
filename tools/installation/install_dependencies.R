#This R program is used to install all dependencies in the QHMM toolset
#Please run this program prior to any other program

#set repo
setRepositories(ind = c(1:6,8))

#Install CRAN package dependencies
pkgTest <- function(x)
{
    if (!require(x,character.only = TRUE))
    {
        install.packages(x,dep=TRUE,repos='http://cran.rstudio.com')
        if(!require(x,character.only=TRUE)) stop("Package not found")
    }
}

pkgTest('optparse')
pkgTest('caret')
pkgTest('tools')
pkgTest('cluster')
#pkgTest('AffyTiling')
pkgTest('cluster')
pkgTest('parallel')
pkgTest('ggplot2')
pkgTest('limma')
pkgTest('Biobase')
pkgTest('doMC')
pkgTest('pamr')
pkgTest('plyr')
#pkgTest('sva')
pkgTest('HiddenMarkov')
pkgTest('markovchain')
pkgTest('MASS')
pkgTest('scales')
pkgTest('ROCR')
pkgTest('pROC')


##Installation Requirement
All tools in QHMM are developed in R 3.1.1, Python 2.7.10 and tested in R 3.1.1, 3.2.2, and Python 2.7.10

##Install dependencies
Please run the following command to install R dependencies before any other tools:
```bash
Rscript install_dependencies.R 
```
The R script will install the following packages to your system:
1. optparse
2. caret
3. tools
4. cluster
5. parallel
6. ggplot2
7. limma
8. Biobase
9. doMC
10. pamr
11. plyr
12. HiddenMarkov
13. markovchain
14. MASS
15. scales
16. ROCR
17. pROC

If any of the above package is failed to be installed, please manually install that package in R environment. 

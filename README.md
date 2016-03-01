#QHMM: a program for cancer biomarker discovery with circulating DNA methylation tiling array

##Synopsis

This repository provides all program tools and source code used in the project "QHMM: a program for cancer biomarker discovery with circulating DNA methylation tiling array". The pipeline for this project includes:
1. Data Cleaning(Preprocessing, this step we provide recommended techniques, but users are free to choose their own preprocessing method)
2. Feature Selection
3. Classification
4. Analysis
Interested users or readers can refer to the original study and supplementary file for further details.

##Motivation
Traditional cancer diagnosis is based on assessing the morphology of cancer cells. However, in some cancers whose cells are not easily accessible, diagnosis by cancer 
cells requires tumour biopsies and this would be a huge pain for patients. Identification of cancer-specific epigenetic DNA alterations in the cell-free circulating 
DNA (cirDNA) is a new opportunity for cancer diagnosis and screening, especially for cancers whose cells are not easily accessi- ble, such as colorectal 
cancer(CRC). Therefore, there is a critical need for new analytical methods and tools for the circulating DNA biomarker discovery research. In this project, we propose a feature selection approach called "QHMM", which is a combination of Hidden Markov Model and a sliding window approach to investigate and select diagnostic marker candidates.

##Tests and Code Examples
Here we show a typical pipeline for analyzing Tiling array data with our tools. The data we used here is a toy dataset from examples folder. We highly encourge users to try qhmm on the toy dataset and follow the walkthrough.txt file, which provides an example of running main pipeline of qhmm

##Installation
System requirement: R 3.2.2 or R 3.1.1, python 2.7.10 or higher verion
Most of the tools used in this project can be running under bash terminal. Some dependencies packages are required to be installed for R and python. No other special installation is needed. (please check walkthrough.txt for installation section)

##Usage

Please refer to the MANUAL.txt for detailed usage information for different tools.

##Contributors
Junjiang Lin, Department of Computer Science, Terrence Donnelly Centre for Cellular and Bimolecular Research, University of Toronto, Toronto, Ontario,Canada
Yue Li, Computer Science and Artificial Intelligence Laboratory, Massachusetts Institute of Technology, Cambridge, MA, United States
Juozas Gordevicius, Institute of Mathematics and Informatics, Vilnius University, Vilnius, Lithuania
Arturas Petronis, Krembil Family Epigenetics Laboratory, Campbell Family Mental Health Research Institute, Centre for Addition and Mental Health, Toronto, Ontario, Canada
Zhaolei Zhang, Department of Computer Science, Terrence Donnelly Centre for Cellular and Bimolecular Research, University of Toronto, Toronto, Ontario,Canada, Banting and Best Department of Medical Research, Department of Molecular Genetics, University of Toronto, Ontario, Canada

##License
GNU license is used for this project. All users have the freedoms to run, study, share, or modify the software.


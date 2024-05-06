# A new method for the reproducible development of aptamers (Neomers) - Visualizations
---
## Overview 
This repository contains the R scripts and necessary data needed to reproduce the visualizations from the PLos One submission titled, "A new method for the reproducible development of aptamers (Neomers)". 

For further information on the neomer method and result visualizations, please refer to our preprint: [A new method for the reproducible development of aptamers (Neomers)](https://biorxiv.org/cgi/content/short/2023.12.19.572437v1)(2023).

## Contents

### 1. Binding assay
This folder contains an R script and the necessary data used to generate the binding assay visualizations for the IL-6 candidate aptamers. The selective binding capabilities of the candidate aptamers were determined and visualized with evaluation against IL-6 and HSA respectively. 

To run, place all files in the Data tables folder in the same directory as the *Binding_assay_neomer_paper.R* script. 

### 2. IL-6 analysis
Data containing the top 10,000 selected IL-6 candidate aptamers, frequencies and fold values generated using the Neomer method are contained within the IL-6 analysis folder. A corresponding R script used to plot the analysis visualizations is also contained in this folder.   
To run, place all files in the Data tables folder in the same directory as the *IL6 analysis plotting.R* script. 

### 3. Secondary structure analysis

#### 1. Predicted structures of the candidate IL6 aptamers:
The secondary structures of the selected candidate aptamer sequences and the Neomer template sequence were generated using [RNAfold](http://rna.tbi.univie.ac.at/cgi-bin/RNAWebSuite/RNAfold.cgi). Sequence information and RNAfold input parameters can be found in the Candidate sequence information and RNAfold and forna parameters folders respectively. The structures were visualized for publication using [forna](http://rna.tbi.univie.ac.at/forna/).

#### 2. Structural diversity visualization
Data obtained from a structural diversity analysis for the Neomer and SELEX libraries were used alongside an R script to visualize the comparisons and differences between the two different libraries. Please refer to the preprint for the link to python script for the structural diversity analysis.
To run, place all files in the Output data folder in the same directory as the *Secondary structure plotting.R* script. 

## Contact
Please contact the corresponding author with any questions.

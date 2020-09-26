# alphadiversitygraphs
## Evenness-Richness alpha diversity scatter plots

Shannon's diversity index is a popular alpha diversity metric because it estimates both richness and evenness in a single equation. 
However, since its value is dependent of both those parameters, there is theoretically an infinite number of richness / evenness value 
combinations translating into the same index score. By decoupling both components measured by Shannon's index, two communities having 
identical indices can be differentiated by mapping richness and evenness coordinates on a scatter plot. In such graphs, confidence ellipses 
would allow testing significant differences between groups of samples. Multivariate statistical tests such as PERMANOVA can be performed on 
distance matrices calculated from richness and evenness coordinates and detect statistically significant differences that would have remained unforeseen otherwise.

*************

## Description
This GitHub repository contains the source code used in the publication "Evenness-Richness scatter plots: a graphical, intuitive approach for alpha diversity analysis".
The article (preprint, not peer-reviewed) is available on BioRxiv: https://www.biorxiv.org/content/10.1101/2020.09.23.310045v1

## Dependencies
 * R v3.4.0+ and RStudio v0.9+
 * ggplot2
 * phyloseq
 
## Usage
Run `source_code.Rmd` in RStudio with all required dependencies installed.

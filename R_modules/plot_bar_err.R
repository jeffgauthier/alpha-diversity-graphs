################################################
### plot_barr_err : Making an abundance barplot with error bars
###                 from a phyloseq object
### By Jeff Gauthier - v0.9 beta 3 - 10 mar 2021
################################################
### This function makes a barplot of the top N OTUs 
### grouped in a given taxonomic rank, with standard 
### error bars for each taxa. Facetting can be used 
### with the "formula" option.
################################################
### Edits:
### 
### 10mar21:  - Replaced SEM error bars by SDM error bars
###
### 30oct18:  - Added an option to choose between top OTUs or top taxa
###           - Sample counts are now normalized
###
################################################
### Usage:
###
### plot_bar_err(physeq, top, by, taxrank, formula, ncol)
###
###    physeq: any formal phyloseq object
###        by: name of the variable by which data is grouped (default "SampleType")
###         n: number of OTUs or taxa to include in plot (default 10)
###  top.type: choose between top OTUs or top taxa (possible values "OTU" or "taxa")
###     facet: a formula for facetting (default ~SampleType)
###      ncol: number of facet columns (default 1)
###
###############################################


################################################
### BEGIN function
################################################

plot_bar_err <- function(physeq, by="SampleType", n=10, top.type="OTU", taxrank=NULL, facet=NULL, ncol=1) {
  
  ## load dependencies
  library(phyloseq)
  library(plyr)
  library(ggplot2)
  library(gmodels)
  
  
  ## Extract top taxa summary statistics (mean and SE) in a data frame
  ## Two possibilities: top N by OTUs or top N by taxonomic rank
  
  if(top.type=="OTU"){
    
    top_names <- names(sort(taxa_sums(physeq), decreasing=TRUE)[1:n])        ## default mode is ascending, the opposite of what we want to do.
    physeq_top <- prune_taxa(top_names, physeq)
    physeq_top_norm <- transform_sample_counts(physeq_top, function (x) {x/sum(x)})
    df <- psmelt(physeq_top_norm)
    avgs <- ddply(df, c(by, taxrank), 
                  function(x) c(mean=mean(x$Abundance),
                                ci_lower=as.numeric(ci(x$Abundance, confidence = 0.95)[2]),
                                ci_upper=as.numeric(ci(x$Abundance, confidence = 0.95)[3])))
    
  } else if(top.type=="taxa") {
    
    physeq_glom <- tax_glom(physeq, taxrank=taxrank)                           ## agglomerate OTUs by a given taxonomic rank prior to summarizing
    top_names <- names(sort(taxa_sums(physeq_glom), decreasing=TRUE)[1:n])   ## default mode is ascending, the opposite of what we want to do.
    physeq_top <- prune_taxa(top_names, physeq_glom)
    physeq_top_norm <- transform_sample_counts(physeq_top, function (x) {x/sum(x)})
    df <- psmelt(physeq_top_norm)
    avgs <- ddply(df, c(by, taxrank), 
                  function(x) c(mean=mean(x$Abundance), 
                                ci_lower=as.numeric(ci(x$Abundance, confidence = 0.95)[2]),
                                ci_upper=as.numeric(ci(x$Abundance, confidence = 0.95)[3])))
  }

  
  ## make the plot
  p <- ggplot(avgs, aes_string(x=taxrank, y="mean", fill=taxrank)) + 
    geom_bar(position='dodge', stat="identity") +
    geom_errorbar(aes(ymin=ci_lower, ymax=ci_upper)) +
    facet_wrap(facet, ncol=ncol) +
    theme(axis.text.x=element_blank()) # to remove the superposed x-labels (there is a legend anyway)
  return(p)
  
}

##############################################
### END function
##############################################

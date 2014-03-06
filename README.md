## Multi-marker and linear mixed model approaches to mapping disease and quantitative trait loci: Module for 2014 Kyoto Course on Bioinformatics

Second Kyoto Course and Symposium on Bioinformatics for
Next-generation Sequencing with Applications in Human Genetics<br>
Center for Genomic Medicine<br>
Kyoto University<br>
March 10-12, 2014

###License

Copyright (C) 2014 [Peter Carbonetto](http://www.cs.ubc.ca/spider/pcarbo)<br>
Dept. of Human Genetics<br>
University of Chicago

This program is free software: you can redistribute it and/or modify
it under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html) as
published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but
**without any warranty**; without even the implied warranty of
**merchantability** or **fitness for a particular purpose**. See
[LICENSE](LICENSE) for more details.

###Introduction

The goal of this module is to explore linear mixed models (LMMs) and
multi-marker models for identifying genetic variants correlated with a
quantitative trait. LMMs and multi-marker mapping methods have been
used to map the genetic loci that contribute to risk of complex human
diseases, but here we will illustrate these methods for mapping
quantitative trait loci (QTLs) in mice instead. Specifically, we will
be working with data from advanced intercrosses of (inbred) lab
strains. In this module, our goal is to develop a better appreciation
for the features of LMMs and multi-marker methods;when and why these
approaches might be useful; and how to interpret the results.

The focus is on cultivating a better understanding of these methods
by example; we will not explore the mathematical foundations of
these methods, at least not in detail.

An advanced intercross population has several features that make it
well-suited for exploring data analysis procedures for genetic
association studies: (1) the patterns of linkage disequilibrium are
more predictable; (2) all observed alleles are common; (3)
environmental conditions are expected to play a smaller role given
that the environment is well-controlled; and (4) perhaps most
significantly for the purposes of this module, in all cases we can
trace an allele to one of the two inbred founders. As a result, the
genetic markers tell us directly whether alleles from a given pair of
individuals are identical-by-descent (IBD), while in general (e.g. in
human studies) we cannot make such an inference.

![Transmission of alleles in F2 cross.](figures/intercross.gif)
**Illustration of the transmission of alleles in an F2 intercross.**
From *Broman and Sen, A Guide to QTL Mapping with R/qtl, 2009.*

In a traditional intercross, two inbred strains are crossed to obtain
the first generation (F1) of mice.  In our case, the inbred founders
are the [LG/J](http://jaxmice.jax.org/strain/000675.html) and
[SM/J](http://jaxmice.jax.org/strain/000687.html) inbred strains
obtained from The Jackson Laboratories, two very widely used inbred
strains in mouse genetics. Since the mother and father have identical
chromosomes, each F1 animal is genetically identical (see the
figure). Therefore, if we denote the LG/J and SM/J inbred strains by A
and B, all the F1 mice have genotype AB on autosomal chromosomes.

The next generation (F2) is obtained by randomly crossing the F1 mice,
so there three possible genotypes appear in the F2 mice, AA, AB and
BB, on autosomal chromosomes. This random mating ("outbreeding")
repeats in subsequent generations to produce an advanced intercross
line. In the data set we will examine, we have phenotype and genotype
data from F2 mice, and in mice from the 34th generation cross. While
understanding the features of advanced intercross lines is not the
central objective of this module, it will be useful to observe them as
we progress.

The traits we will map in this study are those measured as part of
conditioned fear tests, which is intended to assess anxiety-related
behaviour. The amount of "freezing" is measured in three separate days
after exposure to paired tones (the conditioned stimulus) and shocks
(the unconditioned stimulus). The larger goal here to advance our
understanding of the genes underlying anxiety disorders. The fear
conditioning traits are expected to have a highly complex genetic
basis, and any individual loci are expected to explain only a small
proportion of the variance in these traits.

We have also recorded coat colour, which is a major confounder for
these tests, because the albino mice are not tracked as well with the
video equipment. (LG/J mice are albino, and SM/J mice are agouti.)
This also provides an opportunity to attempt to map the loci for these
Mendelian traits.

###Objectives

1. Get some exposure to [github](http://github.com), an excellent
resource for sharing data and source code, and collaborating on
projects.

###Prerequisites

1. You are reasonably comfortable working with using the R interface,
executing R scripts, manipulating variables, working with functions,
and making small changes to scripts.

2. You are familiar, or at least exposed to, the basic notions and
terminology used to describe mapping QTLs, and the analysis of genetic
association studies.

###Getting started

You will need to complete these steps to work on the exercises for
this module.

First, make sure you have a recent version of
[R](http://www.r-project.org) running on your computer.

Next, you will need to download and install R packages
[qtl](http://github.com/kbroman/qtl) and
[QTLRel](http://github.com/pcarbo/QTLRel) using the
**install.packages** function in R, if you don't already have these
packages installed on your computer.

We will also make use of the [varbvs](http://github.com/pcarbo/varbvs)
package. While you can also install this package using
**install.packages**, I recommend getting the most up-to-date version
from github. Follow the "Quick start for R" instructions to install
this package by downloading the source code from github.

Finally, you will need to download the source code for this
module. The simplest way to do this is to [download this repository
as a ZIP archive](http://github.com/pcarbo/kyoto/archive/master.zip).
Alternatively, if you have a github account, you can
[fork](http://help.github.com/articles/fork-a-repo) the repository,
and clone it on your local machine.

Once you have completed these steps, you are ready to move on to the
exercises below.

###Overview of data files

Here is a brief summary of the files in the [data](data) directory:

+ [pheno.csv](data/pheno.csv) Phenotype data from 3-day fear
conditioning study for mice from the F2 and F34 generations of the
LG/J x SM/J advanced intercross line. Includes other information such
as gender, age and coat colour.

+ [geno.csv](data/geno.csv) Genotypes at 4608 SNPs for mice from the
F2 cross and the F34 cross.

+ [map.csv](data/map.csv) Information about the 4608 SNPs, including
genetic distances and chromosomal positions (in bases).

###Overview of R source code files

Here is a brief summary of the files in the [R](R) directory:

+ [read.data.R](R/read.data.R) Defines several functions for
reading experimental cross data from files in comma-delimited ("csv")
format.

+ [data.manip.R](R/data.manip.R) Functions for manipulating the
experimental cross data, and for converting the data into the
various formats used by R/qtl and QTLRel.

+ [mapping.tools.R](R/mapping.tools.R) Functions for analyzing the
QTL experiment data, and for computing marker-based estimates of
pairwise relatedness.

+ [mvnpermute.R](R/mvnpermute.R) Function written by Mark Abney to
execute a permutation-based test with multivariate normally
distributed data.

+ [misc.R](R/misc.R) Various R functions that do not fit in the other
files.

+ [map.qtls.R](R/map.qtls.R) This script maps QTLs across the genome
in a single filial generation of an advanced intercross line (AIL),
using two different "single-marker" methods: one that ignores unequal
relatedness in the mouse population (from the R/qtl library), and
another that attempts to correct for this (from the QTLRel library).

+ [multi.map.qtls.R](R/multi.map.qtls.R) Description goes here.

###Part A

In Part A of this module, we investigate the linear mixed model for
mapping QTLs in the advanced intercross. I've made some of the more
involved questions optional so that we have time to move on to Part B,
and we can return to the more difficult questions later. For all of
Part A, we will work with the [map.qtls.R](R/map.qtls.R) script in R.

**Note:** Because some of the computations can take 10-15 minutes to
complete, I recommend working in teams of 2-4 so that each person can
run the script with different parameters, and then you can compare the
results obtained with your team members.

####Support for association with and without accounting for a polygenic effect

Here will compare genome-wide scans for a polygenic trait in the F2
and F34 samples with and without inclusion of the polygenic effect to
account for population structure. We will start by assessing support
for SNPs that explain variance in freezing after exposure to the tone
on third say of the conditioned fear tests ("freezing to cue"). Set
the script parameters as follows:

    phenotype    <- "freezetocue"
    num.perm.qtl <- 100
    num.perm.rel <- 1
    threshold    <- 0.05
    covariates   <- c("sex","age","albino","agouti")

Set **generation = "F2"**, and reexecute the script with **generation
= "F34"**.

We will contrast these results against a Mendelian trait: whether or
not the mouse is albino. (This is a binary trait, but we can still
treat it as a continuous variable and attempt to map QTLs for these
trait using a linear regression. It would be better to use a logistic
regression.) For this analysis, set **phenotype <- "albino"** and
**covariates <- NULL**.

The script calculates two set of LOD scores for all available SNPs on
chromosomes 1-19, and shows them in a single figure. These LOD scores
are stored in two data frames: **gwscan.qtl**, the output from the qtl
function scanone (the light blue line in the figure); and
**gwscan.rel**, the output from the analogous function in QTLRel,
scanOne (the dark blue line in the figure).

To determine whether or not a LOD score constitutes "significant"
support for an association between genotype and phenotype, we
calculated a threshold for significance by estimating the distribution
of LOD scores under the null hypothesis, and then took the threshold
to be the 100(1 – alpha)th percentile of this distribution, with alpha
= 0.05. It is recommended that this be done with a large number of
replicates (at least 1000) to ensure that the threshold is fairly
stable, but here I use a smaller number (*num.perm.qtl = 100*) so that
the computations can be completed in a reasonable amount of time. This
threshold is shown by the dotted red line in the figure.

**Note:** This procedure fails to account for differences in genetic
sharing among the AIL mice. I also have implemented a test that
accounts for the covariate structure of the polygenic effect. However,
this permutation-based test is much slower, so I set **num.perm.rel =
1**. You are invited to investigate on your own time these two
different methods for permutation-based tests: (1) the method that
accounts for the covariance structure in the AIL when permuting the
data; (2) the standard method that assumes all mice are equally
related.

**Questions**

+ Which QTLs would you report as significant with the basic linear
  regression (qtl), and with the LMM (QTLRel), in the F2 and F34 mice?

+ How would you characterize the support for association using the
basic linear regression, and using the LMM, in the F2 and F34 cohorts?

+ How do the two methods behave in the F2 and F34 populations?

+ It is also useful to compare the genome-wide scans in the F2 and F34
generations, because the patterns of linkage disequilibrium are
different. (The F34 mice have had a much greater opportunity to
accumulate recombinations.) Based on the results in freezetocue and
albino, what can you say about the F2 and F34 mice in terms of: (1)
ability to identify QTLs, (2) ability to pinpoint the gene or genetic
loci underlying these traits?

+ Optional: What locus to you identify for the albino trait, and look
up the associated SNPs in the
[UCSC Genome Browser](http://genome.ucsc.edu) (Mouse Genome Assembly
37) to see whether the region overlaps a known gene for this trait.

+ For the LMM analysis, this script fits the LMM to the data
separately for each chromosome. We can investigate the parameters
corresponding to the variance components of this model. How do these
parameter estimates differ among the chromosomes, and can you observe
any trend? These parameters are stored in the matrix **vcparams**. See
**help(estVC)** for a brief explanation of what variance componentese
these numbers correspond to.

+ Optional (though highly recommended!): Investigate "proximal
contamination" by modifying the script so that **R** is only
calculated once outside the loop that runs over each chromosome. What
happens to the genome-wide scan if instead we compute **R** using
*all* markers? Why does this happen?

####Realized relatedness

Working with the marker-based estimates of genetic sharing, or
relatedness, gives us an opportunity to examine these estimates more
closely. Using the function **rr.matrix** returns an n-by-n matrix,
where n is the number of samples. This matrix is **R** in the script.
Each entry of this matrix is simply the number of alleles that share
the same state, averaged over all SNPs. For a given SNP, this is 0 if
the genotypes of individuals i and j are homozygous and different; 2
if both genotypes are homozygous and the same; and 1 in all other
cases. (Note that, in an AIL, this is equivalent to the number of
alleles that are IBD since all alleles originate from two inbred
founders.) To account for uncertainty in the genotype estimates
whenever the genotypes are missing, we used the following formula for
the expected number of shared alleles.

**Questions**

+ The entries of the realized relatedness matrix are also the kinship
coefficients times 2. The identity coefficients for (i,i) F2 pairs are
d1 = 1/2, d7 = 1/2, and for (i,j) F2 pairs they are d = (1/8, 1/8,
1/4, 0, 1/4, 0, 1/4, 0, 0), and the kinship coefficients are 3/4 and
1/2 for the diagonal and off-diagonal entries of the F2 kinship
matrix, respectively). Looking at the histograms of the diagonal (i,i)
and off-diagonal (i,j) entries of the realized relatedness matrix
**R**, how do the marker-based estimates in the F2 generation compare
to what we would expect, and what does that tell about

+ Compare the relatedness coefficients **R** estimated in the F2 and
F34 mice. What do these relatedness coefficients tell you about
population structure and inbreeding in the these mice?

+ Optional: When we have the probabilities of genotypes AA, AB and BB,
the what is the formula for the expected number of shared alleles
between two individuals?

+ Optional: In human studies, people typically use a different
realized relatedness matrix. Instead of calculating kinship
coefficients based on the genotypes, they calculate the realized
matrix as **R = X % * % t(X)**, where **X** is the n x p genotype
matrix, where n is the number of samples, and p is the number of
markers; that is, the matrix populated with the allele counts (0, 1 or
2), or the mean allele counts if there is some uncertainty in the
genotypes. In my lecture, I claimed that this matrix would yield
identical results. Implement this in the script, and demonstrate
empirically that this relatedness matrix yields the same LOD scores at
all SNPs.

###Part B

In Part B of this module, we will contrast our experiences so far with
the basic linear regression and LMM methods---methods that are based
on a single-marker linear regression---with multi-marker regression
methods that simultaneously consider all markers as potential
predictors of the phenotype. Our goal is not only to understand the
benefits of multi-marker mapping, but critically we need to understand
how to interpret the results. For all of Part B, we will work with the
[map.qtls.R](R/multi.map.qtls.R) script in R.

Here we will focus on mapping QTLs for freezetocue since the albino
trait is not particularly challenging to map, so won't reveal anything
particularly interesting about multi-marker mapping. (However, whis
shouldn't prevent you from trying mapping the albino trait if you
would like to do so.) 

As implemented for this module, the multi-marker analysis consists of
an inner loop and and outer loop. In the outer loop, we try different
combinations of the hyperparameters in order to identify a combination
that fits better than the other combinations. This is sort of like EM,
in the sense that the goal is to arrive at parameter estimates that
maximize the likelihood, except that maximization is done in a not
very intelligent fashion. The combinations of hyperparameters that are
assessed are determined according to **sigma**, **sa** and
**log10odds** which are set near the beginning of the script.

In the inner loop, for each setting of the hyperparameters we
calculate approximate posterior probabilities for the regression
coefficients corresponding to all the markers. The inner loop
calculations are done efficiently using a *variational
approximation*. (Alternatively, we could use MCMC to get better
estimates of the posterior probabilities, but this is much slower.)



**Questions**

###Exit slip

They will return this back to me so that I can get feedback about the
module.

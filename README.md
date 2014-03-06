## Multi-marker and linear mixed model approaches to mapping disease and quantitative trait loci: Module for 2014 Kyoto Course on Bioinformatics

Second Kyoto Course and Symposium on Bioinformatics for
Next-generation Sequencing with Applications in Human Genetics<br>
[Center for Genomic Medicine](http://www.genome.med.kyoto-u.ac.jp)<br>
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

The goal of this module is to explore the use of linear mixed models
(LMMs) and multi-marker models to identify regions of the genome that
are correlated with quantitative traits. While LMMs and multi-marker
mapping methods have been used to map the genetic loci that contribute
to risk of complex diseases, here we will investigate these methods
for mapping quantitative trait loci (QTLs) in mice, and specifically
advanced intercrosses of inbred mouse laboratory strains. We will try
to develop a better appreciation for the features of these approaches,
when and why these approaches might be useful, and how to interpret
the results of the analyses.

What we will **not** do in this module is explore the mathematics of
these methods, at least not in detail. Understanding the mathematical
fundamentals of these methods is certainly important, but here we are
exploring the features of these methods **by example**.

An advanced intercross population has several features that make it
well-suited for exploring data analysis methods for genetic
association studies. These features include: (1) the patterns of
linkage disequilibrium are fairly predictable; (2) all observed
alleles are common; (3) environmental conditions are expected to play
a smaller role given that the environment is well-controlled; and (4)
perhaps most significantly for the purposes of this module, in all
cases we can trace an allele to one of the two inbred founders. As a
result, the genetic markers tell us directly whether alleles from a
given pair of individuals are identical-by-descent (IBD), while in
general (e.g. in human studies) we cannot make such an inference.

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

Outline objectives here.

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
mapping QTLs in the advanced intercross.

###Exit slip

They will return this back to me so that I can get feedback about the
module.

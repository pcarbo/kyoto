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
and multi-marker models to identify regions of the genome that are
correlated with quantitative traits. While these statistical methods
have been used in many studies to map genetic loci harbouring risk
factors for complex diseases, here we will investigate these methods
for mapping quantitative trait loci (QTLs) in advanced intercrosses of
inbred mouse laboratory strains. We will try to get a better
appreciation for the features of these approaches, and when they might
be useful, and how to interpret the results of the analyses.

An advanced intercross population has several features that make it
well-suited for exploring data analysis methods for genetic
association studies. These features include: (1) the patterns of
linkage disequilibrium are fairly predictable; (2) all observed
alleles are common; and (3) in all cases we can trace an allele to one
of the two inbred founders.

![Transmission of alleles in an F2 intercross.](figures/intercross.gif)

In a traditional intercross, two inbred strains are crossed to obtain
the first generation (F1) of mice.  In our case, the inbred founders
are the LG/J and SM/J inbred strains obtained from The Jackson
Laboratories, two very widely used inbred strains in mouse
genetics. Since the mother and father have identical chromosomes, each
F1 animal is genetically identical (see the figure). Therefore, if we
denote the LG/J and SM/J inbred strains by A and B, all the F1 mice
have genotype AB on autosomal chromosomes.

The next generation (F2) is obtained by randomly crossing the F1 mice,
so there three possible genotypes appear in the F2 mice, AA, AB and
BB, on autosomal chromosomes. This random mating ("outbreeding")
repeats in subsequent generations to produce an advanced intercross
line. In the data set we will examine, we have phenotype and genotype
data from F2 mice, and in mice from the 34th generation cross. While
understanding the features of advanced intercross lines is not the
central objective of this module, it will be useful to observe them as
we progress.

###Getting started

Instructions for getting started go here.

###Overview of data files

Here is a brief summary of the files in the [data](data) directory:

+ [pheno.csv](data/pheno.csv) Phenotype data from 3-day fear
conditioning study for mice from the F2 and F34 generations of the
advanced intercross line. Includes other information such as gender,
age and coat colour.

###Overview of R source code files

Here is a brief summary of the files in the [R](R) directory:

+ [read.data.R](code/read.data.R) Defines several functions for
reading experimental cross data from files in comma-delimited ("csv")
format.

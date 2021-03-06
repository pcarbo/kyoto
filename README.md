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
used to map genetic loci that contribute to risk of complex human
diseases. But here I demonstrate the the use of these methods for
mapping quantitative trait loci (QTLs) in mice. Specifically, we will
work with data from an "advanced intercross" of inbred lab strains. In
this module, our goal is to develop a better appreciation for the
features of LMMs and multi-marker methods, when and why these
approaches might be useful, and how to correctly interpret the
results. The focus is on cultivating a better understanding of these
methods by example; in the lecture we will explore the mathematical
foundations of these methods.

An advanced intercross population has several features that make it
well-suited for exploring data analysis procedures for genetic
association studies:

1. The patterns of linkage disequilibrium are
more predictable than in humans;

2. All the alleles are common

3. The environment is better controlled, so environmental conditions
are expected to play a smaller role in complex traits than in humans;
and

4. We can trace each allele to one of the two inbred founders.

As a result of 4, imputing missing genotypes is straightforward. Also,
we can determine whether the alleles from any two individuals at a
given locus are identical-by-descent (IBD); in human studies, this is
usually not possible, or not so straightforward.

![Transmission of alleles in F2 cross.](figures/intercross.gif)
**Illustration of the transmission of alleles in an F2 intercross.**
Source: *Broman and Sen, A Guide to QTL Mapping with R/qtl, 2009.*

In a traditional intercross, two inbred strains are crossed to obtain
the first generation (F1) of mice.  In our case, the inbred founders
are the [LG/J](http://jaxmice.jax.org/strain/000675.html) and
[SM/J](http://jaxmice.jax.org/strain/000687.html) inbred strains
obtained from Jackson Laboratories. Since the mother (A) and father
(B) have identical chromosomes, each F1 mouse is genetically identical
(see the figure). If we denote the LG/J and SM/J alleles by A and B,
all F1 mice have genotype AB at loci on autosomal chromosomes.

The next generation (F2) is obtained by randomly crossing the F1
mice. The F2 mice have three possible genotype outcomes on autosomal
chromosomes, AA, AB and BB. This random mating is repeated in
subsequent generations to produce the advanced intercross. In the data
set used for this module, we have phenotype and genotype data from F2
mice, and in mice from the 34th generation of the advanced
intercross. Understanding the features of an advanced intercross is
not a central objective of this module (I've included some optional
questions about this below).

The quantitative traits we analyze were measured as part of a study of
anxiety-like behaviour in mice. In the "conditioned fear" paradigm to
assess anxiety, the amount of "freezing" is measured over three
separate days after exposure to paired tones (the conditioned
stimulus) and shocks (the unconditioned stimulus). The goal of this
study was to advance our understanding of the genes underlying anxiety
disorders. The fear conditioning traits are expected to have a highly
complex genetic basis; any individual locus is expected to explain
only a small proportion of the total variation in these traits.

We have also recorded coat colour; LG/J mice are albino, SM/J mice are
agouti. Coat colour is a confounder for the fear conditioning
phenotypes because albino mice are not tracked as well with the video
equipment. (If you do not include coat colour as a covariate for the
fear conditioning trait, you will see that the QTL mapping reveals a
strong association near a known pigmentation gene.) This also provides
an opportunity to attempt to map QTLs for a Mendelian trait.

###Objectives

1. Learn how to apply a LMM to a genetic association study, and
discover some of the features of LMMs through practice.

2. Learn how to apply multi-marker methods to a genetic association
study, and discover some the features of multi-marker methods through
practice.

3. Work with genotype and phenotype data in R from a mouse advanced
intercross line.

4. Learn how to use R packages qtl, QTLRel and varbvs.

5. Get some exposure to [github](http://github.com), an online tool
for sharing data and source code, and collaborating on projects (and
creating teaching modules for workshops!).

###Prerequisites

1. You are comfortable using the R interface and executing R
scripts. It would also be helpful if you know how to manipulate
variables in R (e.g. vectors, matrices, lists, data frames), work with
R functions, and make small modifications to R scripts.

2. You are familiar, or at least have exposure to, the basic notions
and standard terminology used to describe statistical analysis of
genetic association studies.

3. You are familiar with core concepts in genetics and genetic
analysis.

###Getting started

To work on the exercises in this module, you will first need to
complete these steps.

First, make sure you have a recent version of
[R](http://www.r-project.org) installed on your computer.

Next, download and install R packages
[qtl](http://github.com/kbroman/qtl) and
[QTLRel](http://github.com/pcarbo/QTLRel) using the
**install.packages** function in R. We will also use the
[varbvs](http://github.com/pcarbo/varbvs) package. You can also
install this package using **install.packages**, but I recommend
getting the most up-to-date version from github, which has a few bug
fixes. Follow the "Quick start for R" instructions to install this
package from the source code on github. In summary, enter the
following commands in R:

    Sys.setenv(http_proxy = "http://proxy.kuins.net:8080")
	install.packages("qtl")
	install.packages("QTLRel")
	install.packages("varbvs")

The first command configures your computer for the network proxy.

Finally, you will need to download the source code for this
module. The simplest way to do this is to [download this repository
as a ZIP archive](http://github.com/pcarbo/kyoto/archive/master.zip).
Alternatively, if you have a github account, you can
[fork](http://help.github.com/articles/fork-a-repo) the repository,
and clone it on your local machine.

Once you have completed these steps, you should be able to run the
scripts in R by starting the program in the **kyoto/R** directory, or
by changing the working directory in R to **kyoto/R** using the
**setwd** function in R. (Here, **kyoto** is the name of the folder
containing all the files you downloaded from this reposistory).

###Overview of data files

Here is a brief summary of the files in the [data](data) directory:

+ [pheno.csv](data/pheno.csv) Phenotype data from 3-day fear
conditioning study for mice from the F2 and F34 generations of the
LG/J x SM/J advanced intercross line. Includes other information such
as gender, age and coat colour.

+ [geno.csv](data/geno.csv) Genotypes at 4608 SNPs for mice from 
F2 and F34 crosses.

+ [map.csv](data/map.csv) Information about the 4608 SNPs, including
genetic distances and chromosomal positions.

###Overview of R source code files

Here is a brief summary of the files in the [R](R) directory:

+ [read.data.R](R/read.data.R) Defines several functions for
reading experimental cross data from files in comma-delimited ("csv")
format.

+ [data.manip.R](R/data.manip.R) Functions for manipulating the
experimental cross data, and for converting the data into the
various formats used by R/qtl and QTLRel.

+ [mapping.tools.R](R/mapping.tools.R) Functions for analyzing the
QTL experiment data, and for computing the marker-based estimates of
pairwise relatedness.

+ [mvnpermute.R](R/mvnpermute.R) Function written by Mark Abney to
execute a permutation-based test with multivariate normally
distributed data.

+ [misc.R](R/misc.R) Various utility functions that do not fit in the
other files.

+ [map.qtls.R](R/map.qtls.R) This script maps QTLs across the genome
in a single filial generation of an advanced intercross line (AIL),
using two different "single-marker" methods: a method that ignores
unequal relatedness in the mouse population ("scanone" from the R/qtl
library); and a method that attempts to correct for this ("scanOne"
from the QTLRel library).

+ [multi.map.qtls.R](R/multi.map.qtls.R) This script maps QTLs across
the genome by simultaneously considering all markers as candidate
predictors for the trait ("multi-marker mapping"). This script uses
the same data as map.qtls.R.

###Part A

In Part A of this module, we investigate the linear mixed model (LMM)
for mapping QTLs in the advanced intercross. For all of Part A, we
will work with the R script [map.qtls.R](R/map.qtls.R). You can run
the script in R with the command **source("map.qtls.R")** provided
that you have set the working directory to the directory containing
this script.

**Important note #1:** Some of the computations take a long time to
complete (as long as 15-20 minutes, depending on your computer). For
this reason, I recommend working in teams (2-4 people) so that each
member of your team can run the script with different settings, and
then you can compare the results you generated with your team
members. To save time, I recommend focusing on the complex trait
("freezetocue"), and save analyses of the Mendelian trait ("albino")
for later.

**Important note #2**: R/qtl and QTLRel may not be the best tools for
working with very large genetic data sets (e.g. with more than 100,000
SNPs). For larger data sets, other programs written in C or Python
(e.g. GEMMA) will be more suitable for your data set. Also, only the
MATLAB version of the varbvs package has been tested on large data
sets.

**Important note #3:** Some of the more involved questions are
optional; please focus on the non-optional questions, so that we have
time for Part B of the module. Later, if there is time, we can return
to the optional questions. 

####QTLs with and without a polygenic effect

Here we compare genome-wide scans for a polygenic trait in the F2 and
F34 cohorts using: (1) a linear regression model that does not include
the polygenic effect (the "basic linear regression"); (2) a linear
regression model that includes the "polygenic effect" (with a
covariance matrix that is intended to capture population
structure). Method #1 ignores unequal relatedness, and Method #2
attempts to correct for this. Method #1 is implemented in "scanone"
from the R/qtl library, and Method #2 is implemented in "scanOne" from
the QTLRel library.

We start by assessing support for SNPs that explain variance in a
behavioural trait. The trait is freezing after exposure to tones
("freezing to cue") on the third day of the conditioned fear test. Set
the parameters at the top of the **map.qtls.R** script as follows:

    phenotype    <- "freezetocue"
	generation   <- "F2"
    num.perm.qtl <- 100
    num.perm.rel <- 1
    threshold    <- 0.05
    covariates   <- c("sex","age","albino","agouti")

Also, try running the script with **generation = "F34"**. (This takes
much longer to complete because many more SNPs are genotyped in the
F34 mice, so I recommend starting this as soon as possible.)

We will contrast these results against results on a Mendelian trait:
whether the mouse's coat is white or not. (This "albino" trait is
binary, but we can still treat it as a continuous phenotype and
attempt to map QTLs using a linear regression.) For this trait, set

    phenotype  <- albino
	covariates <- NULL

In total, you should run **map.qtls.R** 4 times with different
parameter settings.

**LOD scores:** The script calculates two sets of LOD scores for all
available SNPs on chromosomes 1-19. Once all the calculations are
completed, it displays the LOD scores for all the SNPs in a single
figure. In the R workspace, the LOD scores are stored in two data
frames: **gwscan.qtl** is the output from the qtl function **scanone**
which does not include the polygenic effect (this is also the light
blue curve in the figure); and **gwscan.rel**, the output from the
analogous function in **QTLRel**, scanOne, which includes the
polygenic effect (this is also the dark blue line in the figure).

**Sigificance thresholds:** To determine whether or not a LOD score
constitutes "significant" support for an association between genotype
and phenotype, we calculate a threshold for significance by estimating
the distribution of LOD scores under the null hypothesis, then we take
the threshold to be the 100(1 – *a*)th percentile of this
distribution, with *a* = 0.05. It is recommended that the null
distribution be estimated with a large number of replicates (at least
1000) to ensure that the threshold is fairly stable. But in this case
I suggest you use a smaller number of replicates (**num.perm.qtl =
100**) so that the computations can be completed before you have to go
home or return to your hotel room. This significance threshold is
shown by the dotted orange line in the figure.

**Note:** This permutation procedure does not account for differences
in genetic sharing among the AIL mice. I have also implemented a
permutation-based test that observes the covariance structure of the
polygenic effect. However, this permutation-based test is much slower,
so I set **num.perm.rel = 1**. On your own time, you are invited to
investigate these two different methods for estimating the null
distribution of LOD scores: (1) the method that generates permutations
observing the covariance structure of the AIL mice; (2) the standard
method that assumes all mice are equally related.

**Questions**

1. Which QTLs are reported as significant when we use the basic linear
regression (qtl), and when we use the LMM (QTLRel), in the F2 mice,
and in the F34 mice? Do we identify different QTLs, or the same QTLs?

2. What trends do you notice about the association signal using the
basic linear regression compared to the LMM, in the F2 and F34
cohorts? How can you explain these trends?

3. It is also useful to compare the genome-wide scans in the F2 and
F34 generations, because the patterns of linkage disequilibrium are
very different. (The F34 mice have accumulated many more
recombinations.)  Based on the results in freezetocue (and albino),
what do you observe about the F2 and F34 mice in terms of: (1) ability
to identify QTLs? (2) overlap in the QTLs? (3) ability to pinpoint the
location of the genetic polymorphisms underlying these traits? Are the
significance thresholds different in the F2 and F34 cohorts? If so,
why? In what way do the two mapping methods behave differently in the
F2 and F34 populations?

4. **Optional:** What locus do you identify for the albino trait? Does
the QTL region overlap a known gene for this trait? Search for the
associated SNPs in the [UCSC Genome Browser](http://genome.ucsc.edu)
(Mouse Genome Assembly 37) to investigate this.

5. For the LMM data analysis, the script fits the LMM to the data
separately for each chromosome. We can investigate the parameters
corresponding to the variance components of this model. How do these
parameter estimates differ among the chromosomes? Do you observe a
trend in these estimates based on looking at the association signal on
these chromosomes? These parameters are stored in matrix
**vcparams**. See **help(estVC)** for a brief explanation of the
corresponding variance components.

6. **Optional, though highly recommended**: Investigate the idea of
"proximal contamination" by modifying the script so that the matrix
**R** is only calculated once using all SNPs. What happens to the LOD
scores if we compute **R** only once using *all* markers? Can you
explain why these LOD scores are different?

####Realized relatedness (optional)

Working with the marker-based estimates of genetic sharing gives us an
opportunity to examine these estimates more closely in the advanced
intercross. The function **rr.matrix** returns an *n*-by-*n* matrix,
where *n* is the number of samples. In the R script, this matrix is
denoted by **R**. Each entry of this matrix is the estimated number of
alleles that share the same state, averaged over all available
SNPs. For a given SNP, this is 0 if the genotypes of individuals *i*
and *j* are homozygous and different; 2 if both genotypes are
homozygous and the same; and 1 in all other cases. (Note that, in an
AIL, this is equivalent to the number of alleles that are
identity-by-descent, or IBD, since all alleles originate from the two
inbred founders.) To account for uncertainty in the genotype estimates
whenever the genotypes are missing, we use a formula for the expected
number of shared alleles.

**Questions**

1. The entries of the realized relatedness matrix are also the kinship
coefficients, times 2. The condensed identity coefficients for (i,i)
F2 pairs are d1 = 1/2, d7 = 1/2; for (i,j) F2 pairs, they are d =
(1/8, 1/8, 1/4, 0, 1/4, 0, 1/4, 0, 0). The kinship coefficients are
3/4 and 1/2 for the diagonal and off-diagonal entries of the F2
kinship matrix, respectively. Looking at the histograms of the
diagonal (i,i) and off-diagonal (i,j) entries of the realized
relatedness matrix **R**, how do the marker-based estimates of genetic
sharing in the F2 generation compare to the *expected* sharing (*i.e.*
the kinship coefficient), and what does this distribution tell us
about genetic relatedness in these mice?

2. Compare the relatedness coefficients **R** estimated in the F2 and
F34 mice. What do these relatedness coefficients tell us about
differences in these populations? Why do we observe these differences?

3. Given probabilities of genotypes AA, AB and BB, the what is the
formula for the expected number of shared alleles between two
individuals?

4. In human studies, people typically use a different realized
relatedness matrix. Instead of calculating kinship coefficients based
on the genotypes, they calculate the realized matrix as **R <-
matrix.square(X)/p**, where **X** is the *n* x *p* genotype matrix,
*n* is the number of samples, and *p* is the number of markers; that
is, **X** is the matrix populated with allele counts (0, 1 or 2), or
the mean allele counts if there is some uncertainty in the
genotypes. Implement this alternative relatedness matrix in the
script, and demonstrate empirically that this relatedness matrix
yields very similar LOD scores.

###Part B

In Part B of this module, we will contrast our experiences so far with
the approaches based on a single-marker linear regression to
"multi-marker" approaches that simultaneously consider all SNPs as
potential predictors of the phenotype. Our objective is to understand
the features of the multi-marker approach, and to understand how to
correctly interpret the results. For all of Part B, we will work with
the [multi.map.qtls.R](R/multi.map.qtls.R) script in R. You can run
the script in R with the command **source("multi.map.qtls.R")**
provided that you have set the working directory to the directory
containing this script.

In Part B, we map QTLs for freezetocue only since the albino trait is
not particularly challenging, as we have seen, and will not reveal
anything particularly interesting about multi-marker mapping
methods. (However, this should not stop you from trying multi-marker
mapping for the albino trait!) Thus, for this module we set the script
parameters to

    phenotype  <- "freezetocue"
    covariates <- c("sex","age","albino","agouti")

As before, we compare the genome-wide scans in the F2 and F34 mice, so
please run the script with **generation = "F2"** and with **generation
= "F34"**.

In order to correctly interpret the multi-marker mapping results, we
first need cover a few important points:

+ **Posterior probabilities**: In the single-marker mapping, we
typically quantify support for association using LOD scores or
p-values. In the multi-marker mapping, we use posterior probabilities;
specifically, for a given SNP, we report the posterior probability
that the coefficient in the linear regression is not zero, or the
probability that the SNP is *included* in the regression (we call this
the "posterior inclusion probability").  This is what is plotted along
the vertical axis of the figure.

+ **Significance thresholds**: Unlike single-marker mapping, there is
no need to separately calculate a threshold for significance. In the
figure, I have drawn a threshold at a posterior probability of 0.9,
but this is an arbitary threshold.

+ **Posterior approximation**: The interpretation of these
probabilities is somewhat complicated by the fact that I calculate the
posterior probabilities *approximately*.  This is done in order to
make these calculations more efficient, so that we can tackle very
large data sets. In particular, this approximation does not accurately
capture the uncertainty in the included SNP when many nearby SNPs are
strongly correlated with each other. (We could use MCMC instead to get
better estimates of the posterior probabilities, but this is slower,
and does not scale as well to large data sets.) Therefore, you should
interpret a SNP with a high posterior probability as meaning that the
QTL is located *somewhere nearby* the SNP.

+ **Model fitting**: A key component of the multi-marker mapping is
that it *fits the priors to the data*. More precisely, we have an
additional set of parameters (the "hyperparameters") that specify the
prior distribution of the regression coefficients. Fitting the priors
to the data allows us to *calibrate* the posterior probabilities so
that they better reflect the probability that the locus truly contains
a causal polymorphism. In the script, we try different combinations of
hyperparameters (these combinations are chosen by hand), and we settle
on the combination that best fits the data; *i.e.* the hyperparameter
setting that maximizes the probability of the data. After running the
script, the best combination of the hyperparameters is **grid[i,]**.

+ Since we often have very little information about the
hyperparameters in genome-wide association studies, we are often
uncertain about which setting of the hyperparameters is
best. Therefore, it is usually better to account for this uncertainty
by *averaging* over different hyperparameter settings. I did not
implement this averaging for this module because it adds another layer
of complexity.

**Questions**

+ To what extent do the QTLs with the strongest support in the
multi-marker mapping (the highest posterior probabilities) overlap the
QTLs with the strongest support in the single-marker mapping? Does the
genome-wide association signal from the multi-marker mapping
(posterior probabilities) more closely resemble the LOD scores from
the LMM (LMM), or the results using the basic linear regression (qtl)?
How can we explain the similarities or differences?

+ The second question concerns interpretation of the
hyperparameters. Recall, the "best fit" hyperparameters are
**grid[i,]**. Given what we you have learned about the prior
genome-wide log-odds (**grid[i,"log10odds"]**), does the "best fit"
value make sense for the F2 and F34 populations? How can you explain
this estimate?

+ **Optional:** **grid[i,"sa"]** gives the "best fit" for the prior
variance of the regression coefficients. Why should we be cautious in
interpreting this estimate?

###Exit slip

My hope is to get some feedback from all of you to find out which
parts of the module where most successful. I would be grateful you
could answer the three questions below, and send me an email at
pcarbo@uchicago.edu with your answers.

1. Were you and your team members able to successfully execute all the
steps in R? If not, which steps did you have most trouble with?

2. During this module, we explored the features of LMMs and
multi-marker modeling by making observations about the results we
obtained using different analysis approaches. Which observation did
you find most surprising? Why was it so surprising? And what insight
did it give you about the different analysis approaches?

3. Is the topic of this module relevant to your current research? If
yes, how is it relevant? If not, do you expect it to be relevant to
your future research?

I hope you enjoyed this module, and the workshop. Thank you for your
time!

###Acknowledgments

Some of the code and data I used for this module was contributed by
Abraham Palmer, Clarissa Parker, Mark Abney and Riyan Cheng. Also
thank you to Fumihiko Matsuda and Mark Lathrop for organizing the
workshop.

\name{rSEA-package}
\alias{rSEA-package}
\docType{package}
\title{Simultaneous Enrichment Analysis (SEA) of all possible feature-sets using the unified null hypothesis
}
\description{
This package uses raw p-values of genomic features as input and evaluates any given list of feature-sets or pathways. For each set the adjusted p-value and TDP lower-bound are calculated. The type of test can be defined by arguments and can be refined as necessary. The p-values are corrected for every possible set of features, making the method flexible in choice of pathway list and test type.
For more details see: Ebrahimpoor, M (2019) <doi:10.1093/bib/bbz074>
}

\details{
The unified null hypothesis is tested using closed testing procedure and all-resolutions inference. It combines the self-contained and ompetitive approaches in one framework. In short, using p-values of the individual features as input, the package can provide an FWER-adjusted p-value along with a lower bound and a point estimate for the proportion of true discoveries per feature-set. The flexibility in revising the choice of feature-sets without inflating type-I error is the most important property of SEA.
}
\author{
Mitra Ebrahimpoor.

Maintainer: Mitra Ebrahimpoor<m.ebrahimpoor@lumc.nl>
}
\references{
Mitra Ebrahimpoor, Pietro Spitali, Kristina Hettne, Roula Tsonaka, Jelle Goeman,
Simultaneous Enrichment Analysis of all Possible Gene-sets: Unifying Self-Contained
and Competitive Methods, Briefings in Bioinformatics,bbz074
https://doi.org/10.1093/bib/bbz074
}

\keyword{rSEA}

\name{makenetwork}
\alias{makenetwork}

\title{
Make the Network Used for Matching with Two Criteria
}
\description{
This function is of limited interest to most users, and is called by the makematch() function in the package.  Makes the network used in the two-criteria matching method of Zhang et al (2023) with the possible addition of edges that permit some treated individuals to be removed rather than matched (Rosenbaum 2012,2024).
}
\usage{
makenetwork(costL, costR, ncontrols = 1, controlcosts = NULL,
            treatedcosts=NULL, large=100)
}
\arguments{
  \item{costL}{
The distance matrix on the left side of the network, used for pairing.
}
  \item{costR}{
The distance matrix on the right side of the network, used for balancing.
}
  \item{ncontrols}{
One positive integer, 1 for pair matching, 2 for matching two controls to each treated individual, etc.
}
  \item{controlcosts}{
An optional vector of nonnegative costs used to penalize the control-control edges.
}
  \item{treatedcosts}{
An optional vector of nonnegative costs used to penalize the treated-deletion edges.  This option is available only if ncontrols=1.  It is described in Rosenbaum (2024)
and is closely related to the ideas in Rosenbaum (2012).
}
\item{large}{
A large positive number.  Used only if treatedcosts=NULL.  See the parallel discussion in the documentation for makematch.
}
}
\details{
This function creates the network depicted in Figure 1 of Zhang et al. (2023).

A minimum cost flow in this network is found by passing net to callrelax() in the package 'rcbalance'.  If you use callrelax(), I strongly suggest you do this with solver set to 'rrelaxiv'.  The 'rrelaxiv' package has an academic license.  The 'rrelaxiv' package uses Fortran code from RELAX IV developed by Bertsekas and Tseng (1988, 1994) based on Bertsekas' (1990) auction algorithm.
}
\value{
\item{idtreated}{Row identifications for treated individuals}
\item{idcontrol}{Control identifications for control individuals}
\item{net}{A network for use with callrelax in the 'rcbalance' package.}
}
\references{
Bertsekas, D. P., Tseng, P. (1988) <doi:10.1007/BF02288322> The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.

Bertsekas, D. P. (1990) <doi:10.1287/inte.20.4.133> The auction algorithm for assignment and other network flow problems: A tutorial. Interfaces, 20(4), 133-149.

Bertsekas, D. P., Tseng, P. (1994)
<http://web.mit.edu/dimitrib/www/Bertsekas_Tseng_RELAX4_!994.pdf> RELAX-IV: A Faster Version of the RELAX Code for Solving Minimum Cost Flow Problems.

Rosenbaum, P. R. (2012) <doi:10.1198/jcgs.2011.09219> Optimal matching of an optimally chosen subset in observational studies. Journal of Computational and Graphical Statistics, 21(1), 57-71.

Rosenbaum, P. R. (2024) Tightening an observational block design to form an
optimally balanced subdesign.  Manuscript.

Zhang, B., D. S. Small, K. B. Lasater, M. McHugh, J. H. Silber, and P. R. Rosenbaum (2023) <doi:10.1080/01621459.2021.1981337> Matching one sample according to two criteria in observational studies. Journal of the American Statistical Association, 118, 1140-1151.
}
\author{
Paul R. Rosenbaum
}




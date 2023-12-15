\name{addMahal}
\alias{addMahal}
\concept{Matching}
\concept{Mahalanobis distance}

\title{
Rank-Based Mahalanobis Distance Matrix
}
\description{
Adds a rank-based Mahalanobis distance to an exisiting distance matrix.
}
\usage{
addMahal(costmatrix, z, X)
}
\arguments{
  \item{costmatrix}{
An existing cost matrix with sum(z) rows and sum(1-z) columns.  The function checks the compatability of costmatrix, z and X; so, it may stop with an error if these are not of appropriate dimensions.  In particular, costmatrix may come from startcost().
}
  \item{z}{
A vector with z[i]=1 if individual i is treated or z[i]=0 if individual i is control.  The rows of costmatrix refer to treated individuals and the columns refer to controls.
}
  \item{X}{
A matrix with length(z) rows containing covariates or a vector with length(z) containing a single covariate.
}
}
\details{
The rank-based Mahalanobis distance is defined in section 9.3 of Rosenbaum (2020).
Individual covariates are replaced by their ranks before computing the
Mahalanobis distance.  Even when ties are present, the untied variance of
ranks is used.  These adjustments improve the distance when: (i) a covariate
contains an extreme outlier, causing its variance to increase, thereby
causing the distance to ignore large differences in that covariate, and (ii)
rare binary covariates that either match or mismatch by 1, and would have
very small variances if the tied variance were used.
}
\value{
A new distance matrix that is the sum of costmatrix and the rank-based Mahalanobis distances.
}
\references{
Rosenbaum, P. R. (2020) <doi:10.1007/978-3-030-46405-9> Design of Observational Studies (2nd Edition).  New York: Springer.

Rubin, D. B. (1980) <doi:10.2307/2529981> Bias reduction using Mahalanobis-metric matching. Biometrics, 36, 293-298.


}
\author{
Paul R. Rosenbaum
}

\examples{
data(aHDLt)

# names and corresponding rownames help when viewing output
rownames(aHDLt)<-aHDLt$SEQN
z<-aHDLt$z
names(z)<-aHDLt$SEQN

#
# First 12 people
aHDLt[1:12,]
#
# Start with a zero distance matrix.
dist<-startcost(z)
dist[1:3,1:9]
#
# Add Mahalanobis distances to the zero distance matrix.
dist<-addMahal(dist,z,cbind(aHDLt$age,aHDLt$education,aHDLt$female))
round(dist[1:3,1:9],2)

dim(dist)
sum(z)
sum(1-z)
}


\name{addNearExact}
\alias{addNearExact}

\title{
Add a Near-exact Penalty to an Exisiting Distance Matrix.
}
\description{
Add a Near-exact Penalty to an Exisiting Distance Matrix.
}
\usage{
addNearExact(costmatrix, z, exact, penalty = 1000)
}

\arguments{
  \item{costmatrix}{
An existing cost matrix with sum(z) rows and sum(1-z) columns.  The function checks the compatability of costmatrix, z and exact; so, it may stop with an error if these are not of appropriate dimensions.  In particular, costmatrix may come from startcost().
}
  \item{z}{
A vector with z[i]=1 if individual i is treated or z[i]=0 if individual i is control.  The rows of costmatrix refer to treated individuals and the columns refer to controls.
}
  \item{exact}{
A vector with the same length as z.  Typically, exact represent a nominal covariate.  Typically, exact is a vector whose coordinates take a small or moderate number of values.
}
  \item{penalty}{
One positive number.
}
}
\details{
If the ith treated individual and the jth control have different values of exact, then the distance between them in costmatrix is increased by adding penalty.
}
\value{
A penalized distance matrix.
}

\author{
Paul R. Rosenbaum
}
\note{
In two-criteria matching, addNearExact has a different effect when used on the
left side of the network than on the right side.  On the left side, it implements near-exact matching, but on the right side it implements fine or near-fine balance.
Details follow.

On the left, a sufficiently large penalty will maximize the number of individuals exactly matched for exact.  A smaller penalty will tend to increase the number of individuals matched exactly, without prioritizing one covariate over all others.

On the right, a sufficiently large penalty will strive for fine balance, and if that is infeasible, it will achieve near-fine balance.  A smaller penalty will tend to increase balance, without prioritizing one covariate over all others.

In effect, Zubizarreta et al. (2011) seek near-fine, near-exact matching for the same covariate, essentially by placing it on both the left and the right, with a much smaller penalty on the left.  Strive to balance the covariate, and if you can also pair for it, then so much the better.

If the left distance matrix is penalized, it will affect pairing and balance; however, if the right distance matrix is penalized it will affect balance only.

Adding several near-exact penalties for different covariates on the right distance matrix implements a Hamming distance on the joint distribution of those covariates, as discussed in Zhang et al. (2023).  The tighten() function has the Hamming distance
as an option and the example illustrates its use.

Near-exact matching for a nominal covariate is discussed and contrasted with exact matching in Sections 10.3 and 10.4 of Rosenbaum (2020).  Near-exact matching is always feasible, because it implements a constraint using a penalty.  Exact matching may be infeasible, but when feasible it may be used to speed up computations.  For an alternative method of
speeding computations, see Yu et al. (2020) who identify feasible constraints very quickly prior to matching with those constraints.
}
\references{
Rosenbaum, P. R. (2020) <doi:10.1007/978-3-030-46405-9> Design of Observational Studies (2nd Edition).  New York: Springer.

Yang, D., Small, D. S., Silber, J. H. and Rosenbaum, P. R. (2012)
<doi:10.1111/j.1541-0420.2011.01691.x> Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes. Biometrics, 68, 628-636. (Extension of fine balance useful when fine balance is infeasible. Comes as close as possible to fine balance. Implemented in makematch() by placing a large near-exact
penalty on a nominal/integer covariate x1 on the right distance matrix.)

Yu, R., Silber, J. H., Rosenbaum, P. R. (2020) <doi:10.1214/19-STS699> Matching methods for observational studies derived from large administrative databases. Statistical Science, 35, 338-355.

Zhang, B., D. S. Small, K. B. Lasater, M. McHugh, J. H. Silber, and P. R. Rosenbaum (2023) <doi:10.1080/01621459.2021.1981337> Matching one sample according to two criteria in observational studies. Journal of the American Statistical Association, 118, 1140-1151.

Zubizarreta, J. R., Reinke, C. E., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2011) <doi:10.1198/tas.2011.11072> Matching for several sparse nominal variables in a case control study of readmission following surgery. The American Statistician, 65(4), 229-238.
}

\examples{
data(aHDLt)
rownames(aHDLt)<-aHDLt$SEQN
z<-aHDLt$z
names(z)<-aHDLt$SEQN
aHDLt[1:12,]
dist<-startcost(z)
dist<-addNearExact(dist,z,aHDLt$education)
dist[1:3,1:9]
}

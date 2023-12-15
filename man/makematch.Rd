\name{makematch}
\alias{makematch}
\concept{Matching}
\concept{Two-criteria matching}
\concept{Propensity score}
\concept{Network optimization}

\title{
Make a Match Using Two Criteria Matching with Optimal Subset Matching
}
\description{
The function makematch() is called by this package's main function, tighten(), as part of tightening an observational block design.
}
\usage{
makematch(dat, costL, costR, ncontrols = 1, controlcosts = NULL,
          treatedcosts=NULL, large=100, solver="rrelaxiv")
}
\arguments{
  \item{dat}{
A data frame.  Typically, this is the entire data set.  Part of it will be
returned as a matched sample with some added variables.
}
  \item{costL}{
The distance matrix on the left side of the network, used for pairing.  This matrix would most often be made by adding distances to a zero distance matrix created by startcost(), for instance, using addMahal().  In Figure 1 of Zhang et al. (2023), these are the costs on the left treated-control edges.
}
  \item{costR}{
The distance matrix on the right side of the network, used for balancing.  This matrix would most often be made by adding distances to a zero distance matrix created by startcost(), for instance, using addNearExact().  If you do not need a right distance matrix, then initialize it to zero using startcost() and do not add additional distances to its intial form.  In Figure 1 of Zhang et al. (2023), these are the costs on the right control-treated edges.
}
  \item{ncontrols}{
One positive integer, 1 for pair matching, 2 for matching two controls to each treated individual, etc.
}
  \item{controlcosts}{
An optional vector of costs used to penalize the control-control edges.  For instance, one might penalize the use of controls with low propensity scores.
}
  \item{treatedcosts}{
An optional vector of costs that penalize the treated-treated edges used for subset matching in Rosenbaum (2012,2024).  This option is available only if ncontrols=1.  If treatedcosts = NULL, then the cost is set to the largest element of costL, costR and controlcosts multiplied-by-large to prevent subset matching.  Otherwise, if treatedcosts is not NULL, then the ith coordinate of treatedcosts is the cost of not matching treated individual i.
}
  \item{large}{A large positive number.  Used only if treatedcosts = NULL.  See treatedcosts.}
  \item{solver}{Determines the network optimization code that is used.  Options are solver="rrelaxiv" and solver="rlemon".  The Relax IV code of Bertsekas and Tseng (1988) is suggested, with solver="rrelaxiv", but is has an academic license, so various issues may arise.  The makematch() function calls the callrelax() function in
Pimentel's rcbalance package, whose documentation provides additional details, if needed.}
}
\details{
Implements the two-criteria matching method of Zhang et al (2022)
with the possible addition of edges that permit some treated individuals to be removed rather than matched (Rosenbaum 2012).  It is helpful to look at Figure 1 in Zhang et al. (2023) before using this function and Figure 4 in Rosenbaum (2024).
}
\value{
Returns a matched data set.  The matched rows of dat are returned with a new variable mset indicating the matched set.  The returned file is sorted by mset and z.
}
\references{
Bertsekas, D. P., Tseng, P. (1988) <doi:10.1007/BF02288322> The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.

Bertsekas, D. P. (1990) <doi:10.1287/inte.20.4.133> The auction algorithm for assignment and other network flow problems: A tutorial. Interfaces, 20(4), 133-149.

Bertsekas, D. P., Tseng, P. (1994)
<http://web.mit.edu/dimitrib/www/Bertsekas_Tseng_RELAX4_!994.pdf> RELAX-IV: A Faster Version of the RELAX Code for Solving Minimum Cost Flow Problems.

Hansen, B. B. (2007)
<https://www.r-project.org/conferences/useR-2007/program/presentations/hansen.pdf> Flexible, optimal matching for observational studies. R News, 7, 18-24. ('optmatch' package)

Pimentel, S. D. (2016) "Large, Sparse Optimal Matching with R Package rcbalance" <https://obsstudies.org/large-sparse-optimal-matching-with-r-package-rcbalance/> Observational Studies, 2, 4-23. (Discusses and illustrates the use of Pimentel's 'rcbalance' package.)


Rosenbaum, P. R. (1989) <doi:10.1080/01621459.1989.10478868> Optimal matching for observational studies.  Journal of the American Statistical Association, 84(408), 1024-1032. (Discusses and illustrates fine balance using minimum cost flow in a network in section 3.2.  This is implemented using makematch() by placing a large near-exact
penalty on a nominal/integer covariate x1 on the right distance matrix.)

Rosenbaum, P. R., Ross, R. N. and Silber, J. H. (2007)
<doi:10.1198/016214506000001059> Minimum distance matched sampling with fine balance in an observational study of treatment for ovarian cancer. Journal of the American Statistical Association, 102, 75-83.

Rosenbaum, P. R. (2012) <doi:10.1198/jcgs.2011.09219> Optimal matching of an optimally chosen subset in observational studies. Journal of Computational and Graphical Statistics, 21(1), 57-71.

Rosenbaum, P. R. (2020) <doi:10.1007/978-3-030-46405-9> Design of Observational Studies (2nd Edition).  New York: Springer.

Rosenbaum, P. R. (2024) Tightening an observational block design to form an
optimally balanced subdesign.  Manuscript.

Yang, D., Small, D. S., Silber, J. H. and Rosenbaum, P. R. (2012)
<doi:10.1111/j.1541-0420.2011.01691.x> Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes. Biometrics, 68, 628-636. (Extension of fine balance useful when fine balance is infeasible. Comes as close as possible to fine balance. Implemented in makematch() by placing a large near-exact
penalty on a nominal/integer covariate x1 on the right distance matrix.)

Yu, R. (2023) <doi:10.1111/biom.13771> How well can fine balance work for covariate balancing? Biometrics, 79(3), 2346-2356.

Zhang, B., D. S. Small, K. B. Lasater, M. McHugh, J. H. Silber, and P. R. Rosenbaum (2023) <doi:10.1080/01621459.2021.1981337> Matching one sample according to two criteria in observational studies. Journal of the American Statistical Association, 118, 1140-1151.
(This is the basic reference for two-criteria matching, which generalizes matching with fine balance.)

Zubizarreta, J. R., Reinke, C. E., Kelz, R. R., Silber, J. H. and Rosenbaum, P. R. (2011) <doi:10.1198/tas.2011.11072> Matching for several sparse nominal variables in a case control study of readmission following surgery. The American Statistician, 65(4), 229-238.
(This paper combines near-exact matching and fine balance for the same covariate.  It is implemented in makematch() by placing the same covariate on the left and the right.)
}
\author{
Paul R. Rosenbaum
}

\examples{

# This is a simple example to illustrate makematch().
# Please see the better match in the documentation for
# tighten() before continuing with this example.
#

###########
# Tighten match from 1-3 to 1-2 to adjust for BMI.
# BMI might be affected by alcohol consumption, so
# the primary comparison did not adjust for it.
###########

# The example below illustrates mechanics.
# See the same example in tighten() for a
# simpler construction of a better match.

data(aHDLt)
z<-aHDLt$z   #treatment indicator

# If you need to debug, it is helpful to have the same names for z
# and rownames for dat.
# These names are then used by functions that create distance matrices,
# including startcost(), addNearExact() and addMahal().
# If you use names, they will be checked for consistency and an error
# may result from incompatible names.

rownames(aHDLt)<-aHDLt$SEQN
names(z)<-aHDLt$SEQN

#  Create a zero cost (i.e., distance) matrix for the left and right
#  sides of the network.
left<-startcost(z)
right<-startcost(z)

left<-addNearExact(left,z,aHDLt$block) # Forces within-block matching.
#  Prefer to retain within blocks people who are closest for age, education.
left<-addMahal(left,z,cbind(aHDLt$education,aHDLt$age))

# Try to balance the categories of BMI which are out of balance.
right<-addNearExact(right,z,aHDLt$ibmi,penalty=20)

m<-makematch(aHDLt,left,right,ncontrols=2,large=10)

# Tightened blocks
table(aHDLt$z)
table(m$z)
table(table(aHDLt$block))
table(table(m$block))
table(table(m$mset))

boxplot(m$bmi[m$z==1],m$bmi[m$z==0],
        names=c("D","C"),ylab="BMI")

# Cost matrix for left side of network.
# Rows are treated, columns are control.
z[1:10]
round(left[1:5,1:4],1)

# Cost matrix for right side of network.
round(right[1:5,1:4],1)

}


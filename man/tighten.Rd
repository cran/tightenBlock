\name{tighten}
\alias{tighten}
\concept{Matching}
\concept{Tightening a block design}
\concept{Subset matching}
\concept{Network optimization}
\concept{Hamming distance}

\title{
Tightening an Observational Block Design
}
\description{
Implements a simple version of optimal balanced tightening of an observational
block design.  Finer control of the tightening may be obtained using the
makematch() function; however, that requires more attention to detail.
}
\usage{
tighten(dat, z, block, x = NULL, f = NULL,  ncontrols = 1,
    subset = NULL, pspace=10, solver="rrelaxiv")
}
\arguments{
  \item{dat}{
A data frame for all data in the block design prior to tightening.
Part of dat will be
returned as a tightened block design with some added variables.
}
  \item{z}{
The vector z indicates treatment/control, where 1 indicates treated, and 0
indicates control.  The number of rows of dat must equal length(z); otherwise,
the tighten will stop with an error message.
}
  \item{block}{
A vector indicating the block.  Every block must contain exactly one treated individual; that is, one individual whose value of z is 1.  The minimum
block size is two.  The length of z must agree with the length of block.  An
error will result if these conditions are not met.
}
  \item{x}{
If not NULL, x is a matrix, dataframe or vector containing covariates to be used
in a robust, rank based Mahalanobis distance (Rosenbaum 2020, Section 9.3).
This is the within-block distance.  If x is NULL, then the within-block distance
is zero, and only between block distances affect tightening.  In you want more information, see the documentation in this package for the function addMahal().
}
  \item{f}{
If not NULL, f determines the between block distances.  In general, f describes
fairly coarse nominal categories, either for one or several variables.
If f is a vector or a
factor, then it is used in fine or near-fine balance.  If f is a matrix or a dataframe, then each column is a nominal variable, and they define a Hamming
distance that counts the number of columns that differ.
If f is a vector or a factor, then its length must
equal the length of z.  Otherwise, the number of rows of f must equal the length
of z.  If you want more information, see the documentation in this package for the function addNearExact().
}
  \item{ncontrols}{
  If the blocks are initially of size J, with one treated and J-1 controls,
  then the tightened design will have blocks of size 1+ncontrols with one
  treated and ncontrols controls.  If J-1 > 1, then ncontrols must be less
  than J-1.  For instance, if J=4, then ncontrols might be either 1 or 2.
  If J=2 so J-1=1, then ncontrols must be 1 and subset must not be
  NULL.
}
  \item{subset}{If subset is NULL, then no blocks are discarded.  Otherwise,
  subset should be a positive number.  Roughly speaking, if subset=50, then
  tighten will prefer to discard a block than face a cost or distance of 50
  by including the block.  The distance here is the sum of within and between
  block distances.  This option is available only if the blocks are of size
  two, i.e., pairs.  For blocks larger than pairs, we prefer to reduce their
  size rather than eliminate them.  An error will result if the subset is not
  NULL when the blocks are larger than pairs.
}
  \item{pspace}{This parameter is slightly technical, so it may be best
  to use the default until forced to do otherwise.  Priorities in the
  matching are enforced by penalized costs.  Most important, with the
  largest penalty, is to retain the original block structure.  Second,
  if subset=NULL, then no blocks are to be deleted.  Third, is the fine
  balance/Hamming distance constraint.  The lowest priority with no penalty is the
  Mahalanobis distance within blocks.  pspace is used to space the
  penalties, so priorities are respected.  Increasing pspace emphasizes the
  relative importance of priorities, but may slow down the optimization.
  In general, the largest lower priority penalized distance is multiplied by
  pspace to produce the new penalty.  The easist way to check that pspace
  is large enough is to increase it; if it is large enough, the match should
  be the same, but may take longer to compute.
}
\item{solver}{Determines the network optimization code that is used.  Options are solver="rrelaxiv" and solver="rlemon".  The Relax IV code of Bertsekas and Tseng (1988) is suggested, with solver="rrelaxiv", but is has an academic license, so various issues may arise.  The makematch() function calls the callrelax() function in
Pimentel's rcbalance package, whose documentation provides additional details, if needed.}
}
\details{
The tighten function produces a simple version of an optimally
tightened block design, combining a Mahalanobis distance within-blocks with
some version of fine balancing between blocks.  You can achieve finer
control and subler effects using the makematch function, whose structure is
closer to the mathematical structure of the tightening problem.  The examples
in the makematch function document these finer features.
}
\value{
Returns a data set for a tightened block design.  The matched rows of dat are returned with a new variable mset indicating the matched set.  The returned file is sorted by mset and z.
}
\references{

Bertsekas, D. P., Tseng, P. (1988) <doi:10.1007/BF02288322> The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.

Bertsekas, D. P. (1990) <doi:10.1287/inte.20.4.133> The auction algorithm for assignment and other network flow problems: A tutorial. Interfaces, 20(4), 133-149.

Bertsekas, D. P., Tseng, P. (1994)
<http://web.mit.edu/dimitrib/www/Bertsekas_Tseng_RELAX4_!994.pdf> RELAX-IV: A Faster Version of the RELAX Code for Solving Minimum Cost Flow Problems.

Hansen, B. B. (2007)
<https://www.r-project.org/conferences/useR-2007/program/presentations/hansen.pdf> Flexible, optimal matching for observational studies. R News, 7, 18-24. ('optmatch' package)

Pimentel, S. D. (2016) "Large, Sparse Optimal Matching with R Package rcbalance" <https://obsstudies.org/large-sparse-optimal-matching-with-r-package-rcbalance/> Observational Studies, 2, 4-23. (Discusses and illustrates the use of Pimentel's 'rcbalance' package.)

Rosenbaum, P. R. (1984) <doi:10.2307/2981697> The consequences of adjustment for a concomitant variable that has been affected by the treatment. Journal of the Royal Statistical Society Series A: Statistics in Society, 147(5), 656-666. (Related
to the BMI example.)

Rosenbaum, P. R. (1989) <doi:10.1080/01621459.1989.10478868> Optimal matching for observational studies.  Journal of the American Statistical Association, 84(408), 1024-1032. (Discusses and illustrates fine balance using minimum cost flow in a network in Section 3.2.)

Rosenbaum, P. R. (2006) <doi:10.1093/biomet/93.3.573> Differential effects and generic biases in observational studies. Biometrika, 93(3), 573-586.  (Related
to the dentist example.)

Rosenbaum, P. R. (2012) <doi:10.1198/jcgs.2011.09219> Optimal matching of an optimally chosen subset in observational studies. Journal of Computational and Graphical Statistics, 21(1), 57-71.

Rosenbaum, P. R. (2020) <doi:10.1007/978-3-030-46405-9> Design of Observational Studies (2nd Edition).  New York: Springer.  (Discusses fine balance and the
rank-based Mahalanobis distance.)

Rosenbaum, P. R. (2024) Tightening an observational block design to form an
optimally balanced subdesign.  Manuscript.

Yang, D., Small, D. S., Silber, J. H. and Rosenbaum, P. R. (2012)
<doi:10.1111/j.1541-0420.2011.01691.x> Optimal matching with minimal deviation from fine balance in a study of obesity and surgical outcomes. Biometrics, 68, 628-636. (Extension of fine balance useful when fine balance is infeasible. Comes as close as possible to fine balance.  )

Yu, R. (2023) <doi:10.1111/biom.13771> How well can fine balance work for covariate balancing? Biometrics, 79(3), 2346-2356.

Zhang, B., D. S. Small, K. B. Lasater, M. McHugh, J. H. Silber, and P. R. Rosenbaum (2023) <doi:10.1080/01621459.2021.1981337> Matching one sample according to two criteria in observational studies. Journal of the American Statistical Association, 118, 1140-1151.  (This method generalizes the concept
of fine balance.  Here, it is used with the Hamming distance when f is a matrix.)
}
\author{
Paul R. Rosenbaum
}

\examples{

# The two examples below are from Rosenbaum (2024).

###########
# Tighten match from 1-3 to 1-2 to adjust for BMI.
# BMI might be affected by alcohol consumption, so
# the primary comparison did not adjust for it.
# See Rosenbaum (1984).
###########

data(aHDLt)
result<-tighten(aHDLt,aHDLt$z,aHDLt$block,
       x=cbind(aHDLt$age,aHDLt$education),
       f=cbind(aHDLt$ibmi,(aHDLt$bmi>22.5)+(aHDLt$bmi>27.5)+(aHDLt$bmi>32.5)),
                ncontrols=2)


omit<-aHDLt[!is.element(aHDLt$SEQN,result$SEQN),]
boxplot(result$bmi[result$z==1],result$bmi[result$z==0],omit$bmi,
    names=c("D","C","O"),ylab="BMI")
boxplot(result$hdl[result$z==1],result$hdl[result$z==0],omit$hdl,
    names=c("D","C","O"),ylab="HDL Cholesterol")

# Tightened blocks
table(aHDLt$z)
table(result$z)
table(table(aHDLt$block))
table(table(result$block))
table(table(result$mset))


##################
# Dentist in 1 Year Differential Effect
##################

#
#  Example of tightening in support of a differential
#  comparison to address a generic bias; see Rosenbaum (2006).
#
#
# First, create a data frame in which each block contains a
# differential comparison.
#
# Identify blocks in which the treated individual did not go
# to the dentist in the past year, but at lease one control did.
# Vector dife picks out the blocks that have this pattern.
data(aHDLt)
dif1<-tapply(((aHDLt$z==1)&(aHDLt$dentist1Y==0)),aHDLt$block,sum)==1
dif0<-tapply(((aHDLt$z==0)&(aHDLt$dentist1Y==1)),aHDLt$block,sum)>=1
dif<-dif1&dif0
dife<-as.vector(rbind(dif,dif,dif,dif))
elig<-((aHDLt$z==1)&(aHDLt$dentist1Y==0))|((aHDLt$z==0)&(aHDLt$dentist1Y==1))

# Now, form a data.frame for the people eligible for this comparison
# The tightened match will occur in this data.frame.
aHDLd<-cbind(aHDLt,elig)[dife,]
aHDLe<-aHDLd[((aHDLd$z==1)&(aHDLd$dentist1Y==0))|
            ((aHDLd$z==0)&(aHDLd$dentist1Y==1)),]
rm(dif,dif0,dif1,dife,elig,aHDLd)

#  With data frame aHDLe, several tightened block designs
#  are constructed.

x<-cbind(aHDLe$age,aHDLe$education)

# No pairs are deleted.  Education is not quite balanced.
m1<-tighten(aHDLe,aHDLe$z,aHDLe$block,x=x,f=aHDLe$education)
table(m1$z,m1$education)
tapply(m1$age,m1$z,summary)

# 6 pairs are deleted.  Education is alsmost balanced.
m2<-tighten(aHDLe,aHDLe$z,aHDLe$block,x=x,f=aHDLe$education,subset=150)
table(m2$z,m2$education)
tapply(m2$age,m2$z,summary)

# 13 pairs are deleted.  Perhaps too many. Education is balanced.
m3<-tighten(aHDLe,aHDLe$z,aHDLe$block,x=x,f=aHDLe$education,subset=50)
table(m3$z,m3$education)
tapply(m3$age,m3$z,summary)

oldpar<-par(mfrow=c(1,3))
barplot(t(table(m1$education,1-m1$z)),beside=TRUE,ylim=c(0,50),
        ylab="Count",xlab="Education Level",main="All 118 Pairs",
        col=gray.colors(2))
barplot(t(table(m2$education,1-m2$z)),beside=TRUE,ylim=c(0,50),
        ylab="Count",xlab="Education Level",main="Best 112 Pairs")
barplot(t(table(m3$education,1-m3$z)),beside=TRUE,ylim=c(0,50),
        ylab="Count",xlab="Education Level",main="Best 105 Pairs")
par(oldpar)

}


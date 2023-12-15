\name{tightenBlock-package}
\alias{tightenBlock-package}
\alias{tightenBlock}
\docType{package}
\title{
\packageTitle{tightenBlock}
}
\description{
\packageDescription{tightenBlock}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{tightenBlock}
\packageIndices{tightenBlock}
}
\author{
\packageAuthor{tightenBlock}

Maintainer: \packageMaintainer{tightenBlock}
}
\references{

Bertsekas, D. P., Tseng, P. (1988) <doi:10.1007/BF02288322> The relax codes for linear minimum cost network flow problems. Annals of Operations Research, 13, 125-190.

Rosenbaum, P. R., Ross, R. N. and Silber, J. H. (2007) <10.1198/016214506000001059> Minimum distance matched sampling with fine balance in an observational study of treatment for ovarian cancer. Journal of the American Statistical Association, 102(477), 75-83.

Rosenbaum, P. R. (2012) <doi:10.1198/jcgs.2011.09219> Optimal matching of an optimally chosen subset in observational studies. Journal of Computational and Graphical Statistics, 21(1), 57-71.

Rosenbaum, P. R. (2024) Tightening an observational block design to form an
optimally balanced subdesign.  Manuscript.

Zhang, B., D. S. Small, K. B. Lasater, M. McHugh, J. H. Silber, and P. R. Rosenbaum (2023) <doi:10.1080/01621459.2021.1981337> Matching one sample according to two criteria in observational studies. Journal of the American Statistical Association, 118, 1140-1151.
}
\keyword{ package }

\examples{
data(aHDLt)
result<-tighten(aHDLt,aHDLt$z,aHDLt$block,
       x=cbind(aHDLt$age,aHDLt$education),
       f=cbind(aHDLt$ibmi,(aHDLt$bmi>22.5)+(aHDLt$bmi>27.5)+(aHDLt$bmi>32.5)),
                ncontrols=2)
}
\concept{Tightening blocks}
\concept{Causal inference}
\concept{Observational study}
\concept{Fine balance}
\concept{Optimal subset matching}
\concept{Two criteria matching}
\concept{Block design}

\name{startcost}
\alias{startcost}

\title{
Initialize a Distance Matrix.
}
\description{
Creates an distance matrix of zeros of dimensions compatible with the treatment indicator vector z.
}
\usage{
startcost(z)
}
\arguments{
  \item{z}{
A vector with z[i]=1 if individual i is treated or z[i]=0 if individual i is control.  The rows of costmatrix refer to treated individuals and the columns refer to controls.  Although not strictly required, it is best that z has names that are the same as the names of the data frame dat that will be used in matching.
}
}
\value{
A matrix of zeros with sum(z) rows and sum(1-z) columns.  If z has names, then they become the row and column names of this matrix.
}
\author{
Paul R. Rosenbaum
}

\examples{
# Although names are not required, they aid clarity.  The cost matrix
# has sum(z) rows and sum(1-z) columns, where z is a binary vector.
# If needed, it is easier to connect an entry in z to a row or column in
# cost if they have the same name.  If z has names, the data matrix should
# have the same names.
data(aHDLt)
rownames(aHDLt)<-aHDLt$SEQN
z<-aHDLt$z
names(z)<-aHDLt$SEQN
dist<-startcost(z)
dist[1:3,1:9]
}


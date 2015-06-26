\name{mvtboost-package}
\alias{mvtboost-package}
\alias{mvtboost}
\docType{package}
\title{
\packageTitle{mvtboost}
}
\description{
\packageDescription{mvtboost}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{mvtboost}
\packageIndices{mvtboost}
The most important function is \code{mvtb}, which fits the multivariate tree boosting model. See \code{?mvtb} for details.
}
\author{
\packageAuthor{mvtboost}

Maintainer: \packageMaintainer{mvtboost}
}
\references{
Miller P.J., Lubke G.H, McArtor D.B., Bergeman C.S. (2015) Finding structure in data: A data mining alternative to multivariate multiple regression. Psychological Methods.
}
\keyword{ Boosting, multivariate responses }
\seealso{
\code{gbm}
}
\examples{
data(trees)
X <- as.matrix(trees[,1:2])
Y <- as.matrix(trees[,3])
mvtb(X=X,Y=Y,n.trees=30)
}

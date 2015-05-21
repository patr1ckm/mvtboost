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
~~ An overview of how to use the package, including the most ~~
~~ important functions ~~
}
\author{
\packageAuthor{mvtboost}

Maintainer: \packageMaintainer{mvtboost}
}
\references{
Miller P.J., Lubke G.H, McArtor D.B., Bergeman C.S. (2015) Finding structure in data: A data mining alternative to multivariate multiple regression. Psychological Methods.
}
~~ Optionally other standard keywords, one per line, from file ~~
~~ KEYWORDS in the R documentation directory ~~
\keyword{ package }
\seealso{
\code{\link[gbm-package]{<pkg>}}
}
\examples{
data(trees)
X <- as.matrix(trees[,1:2])
Y <- as.matrix(trees[,3])
mvtb(X=X,Y=Y,n.trees=30)
}

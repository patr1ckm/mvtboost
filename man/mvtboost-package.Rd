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
data(wellbeing)
Y <- wellbeing[,18:23]
X <- wellbeing[,1:17]
Ys <- scale(Y)
cont.id <- unlist(lapply(X,is.numeric))
Xs <- scale(X[,cont.id])
res <- mvtb(Y=Ys,X=X)
summary(res)
plot(res,predictor.no = 11)
predict(res,newdata=Xs)
mvtb.cluster(res)
mvtb.heat(t(mvtb.ri(res)))
}

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
The most important function is \code{mvtb}, which fits the multivariate tree boosting model. See \code{?mvtb} for details. The fitted model objects have \code{summary}, \code{print}, \code{plot} and \code{predict} methods. Additionally, \code{mvtb.ri} to computes the relative influence of each predictor, and \code{mvtb.covex} computes an estimate of the covariance explained in pairs of outcomes by predictors. These tables can be displayed as heatmaps using \code{mvtb.heat}. Examples for fitting, tuning and interpreting the models are available in the help pages and two package vignettes:

vignette("mvtboost_vignette")
vignette("mvtboost_wellbeing")

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
Y <- wellbeing[,21:26]
X <- wellbeing[,1:20]
Ys <- scale(Y)
cont.id <- unlist(lapply(X,is.numeric))
Xs <- scale(X[,cont.id])

res <- mvtb(Y=Ys,X=Xs)

summary(res)
plot(res,predictor.no = 8)
predict(res,newdata=Xs)

covex <- mvtb.covex(res, Y=Ys, X=Xs)
mvtb.cluster(covex)
par(mar=c(4,7,1,1))
mvtb.heat(covex,cexRow=.8)
par(mar=c(5,5,1,1))
mvtb.heat(t(mvtb.ri(res)),cexRow=.8,cexCol=1,dec=0)
}

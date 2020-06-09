\name{SubgrpID-package}
\alias{SubgrpID-package}
\docType{package}
\title{
Biomarker signature development - Exploratory Subgroup Identification package
}
\description{
Prognostic and preditive biomarker signature development for Exploratory Subgroup Identification in Randomized Clinical Trials
}
\details{
\tabular{ll}{
Package: \tab SubgrpID\cr
Type: \tab Package\cr
Version: \tab 1.06\cr
Date: \tab 2016-05-05\cr
License: \tab Proprietary by Exploratory Statistics, AbbVie Inc.(internal evaluation only)\cr
}
}
\author{
Xin Huang and Yan Sun
Maintainer: Yan Sun <sun.yan@abbvie.com> and Xin Huang <xin.huang@abbvie.com>
}
\references{
  Huang X. et al. (2016)
  Biomarker Signatures for Patient Subgroup Selection in Clinical Drug Development. \emph{Statistics in Medicine} (Submitted)
}
\keyword{ AIM }
%\seealso{
%}
\examples{
data(Sepsis.train)
data(Sepsis.test)

yvar="survival"
xvars=names(Sepsis.train)[2:12]
trtvar="THERAPY"

set.seed(123)
subgrp <- SubgrpID(data.train=Sepsis.train,
                  yvar=yvar,
                  trtvar=trtvar,
                  trtref="active",
                  xvars=xvars,
                  type="b",
                  des.res = "smaller",
                  method="AIM.Rule")
subgrp$res
subgrp$train.stat
subgrp$train.plot

}

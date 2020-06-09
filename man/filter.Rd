\name{filter}
\alias{filter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Filter function for predictive/prognostic biomarker candidates for signature development
}
\description{
Filter function for Prognostic and preditive biomarker signature development for Exploratory Subgroup Identification in Randomized Clinical Trials
}
\usage{
filter(data,type,yvar,xvars,censorvar,trtvar,trtref=1,n.boot=50,cv.iter=20,pre.filter=length(xvars), filter.method=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
input data frame 
}
  \item{type}{
type of response variable: "c" continuous; "s" survival; "b" binary
}
  \item{yvar}{
variable (column) name for response variable
}
  \item{xvars}{
vector of variable names for predictors (covariates)
}
  \item{censorvar}{
variable name for censoring (1: event; 0: censor), default = NULL
}
  \item{trtvar}{
variable name for treatment variable, default = NULL (prognostic signature)
}
  \item{trtref}{
coding (in the column of trtvar) for treatment arm, default = 1 (no use for prognostic signature)
}
  \item{n.boot}{
number of bootstrap for the BATTing procedure
}
  \item{cv.iter}{
Algotithm terminates after cv.iter successful iterations of cross-validation, or after max.iter total iterations, whichever occurs first
}
  \item{pre.filter}{
NULL (default), no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
}
  \item{filter.method}{
NULL (default), no prefiltering; "univariate", univaraite filtering; "glmnet", glmnet filtering
}

}
\details{ 
The function contains two algorithms for filtering high-dimentional multivariate (prognostic/predictive) biomarker candidates via univariate fitering (used p-values of group difference for prognostic case, p-values of interaction term for predictive case); LASSO/Elastic Net method. (Tian L. et al 2012)
}
\value{

  \item{var}{a vector of filter results of variable names}

}
\references{
Tian L, Alizadeh A, Gentles A, Tibshirani R (2012) A Simple Method for Detecting Interactions between a Treatment and a Large Number of Covariates. arXiv
}
\author{
Xin Huang, Yan Sun, Paul Trow, Arunava Chakravartty, Lu Tian, Viswanath Devanarayan
}
\note{
TBE
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
data(Sepsis.train)

yvar="survival"
xvars=names(Sepsis.train)[2:12]
trtvar="THERAPY"
trtref="active"
set.seed(123)
filter.res <- filter(data=Sepsis.train,
                     type="b",
                     yvar=yvar,
                     xvars=xvars,
                     trtvar=trtvar,
                     trtref=trtref,
                     pre.filter=20, 
                     filter.method="univariate")

filter.res

set.seed(123)
filter.res <- filter(data=Sepsis.train,
                     type="b",
                     yvar=yvar,
                     xvars=xvars,
                     trtvar=trtvar,
                     trtref=trtref,
                     pre.filter="opt", 
                     filter.method="glmnet")

filter.res
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ filter }

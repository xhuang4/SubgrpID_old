\name{SubgrpID}
\alias{SubgrpID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Exploratory Subgroup Identification main function
}
\description{
Prognostic and preditive biomarker signature development for Exploratory Subgroup Identification in Randomized Clinical Trials
}
\usage{
SubgrpID(data.train, data.test=NULL, yvar, censorvar=NULL, trtvar=NULL, trtref=NULL, xvars, type="c",n.boot=ifelse(method=="PRIM",0,25), des.res="larger",min.sigp.prcnt=0.20, pre.filter=NULL, filter.method=NULL, k.fold=5, cv.iter=20, max.iter=500, mc.iter=20, method=c("AIM.Rule"), train.percent.prim=0.5, do.cv=FALSE, out.file=FALSE, file.path="", plots=F, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.train}{
data frame for training dataset
}
  \item{data.test}{
data frame for testing dataset, default = NULL
}
  \item{yvar}{
variable (column) name for response variable
}
  \item{censorvar}{
variable name for censoring (1: event; 0: censor), default = NULL
}
  \item{trtvar}{
variable name for treatment variable, default = NULL (prognostic signature)
}
  \item{trtref}{
coding (in the column of trtvar) for treatment arm
}
  \item{xvars}{
vector of variable names for predictors (covariates)
}
  \item{type}{
type of response variable: "c" continuous; "s" survival; "b" binary
}
  \item{n.boot}{
number of bootstrap for batting procedure, or the variable selection procedure for PRIM; for PRIM, when n.boot=0, bootstrapping for variable selection is not conducted
}
  \item{des.res}{
the desired response. "larger": prefer larger response. "smaller": prefer smaller response
}
  \item{min.sigp.prcnt}{
desired proportion of signature positive group size for a given cutoff
}
  \item{pre.filter}{
NULL (default), no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
}
  \item{filter.method}{
NULL (default), no prefiltering; "univariate", univariate filtering; "glmnet", glmnet filtering; "unicart", univariate rpart filtering for prognostic case
}
  \item{k.fold}{
cross-validation folds
}
  \item{cv.iter}{
Algotithm terminates after cv.iter successful iterations of cross-validation, or after max.iter total iterations, whichever occurs first
}
  \item{max.iter}{
total iterations, whichever occurs first
}
  \item{mc.iter}{
number of iterations for the Monte Carlo procedure to get a stable "best number of predictors"
}
  \item{method}{
algorithms performed for subgroup identification, one of the ("AIM", "AIM.Rule", "Seq.BT", "PRIM", "SIDES", "optAUC")
}
  \item{train.percent.prim}{
percentage of the sub-training set used only by PRIM method; if train.percent.prim=1, all data will be used both for sub-training and sub-testing inside the PRIM
  }
  \item{do.cv}{
whether to perform cross validation for performance evaluation. TRUE or FALSE (Default)
}
  \item{out.file}{
Name of output result files. If NULL no output file would be saved
}
  \item{file.path}{
default: current working directory. When specifying a dir, use "/" at the end. e.g. "TEMP/"
}
  \item{plots}{
default: F. whether to save plots
}
  \item{...}{
Optional SIDES arguments:
    \describe{
      \item{M:}{ numberof splits retained}
      \item{gamma:}{ selection parameter on stepwise "improvement"}
      \item{maxncov:}{upper bound of numbers of covariates in each subgroup description}
      \item{minSplit:}{minimal size of subjects in one trt arm}
      \item{minNode:}{minimal size of either half in one split}
      \item{maxCand:}{if not NULL, number of randomly drawn covariates to be used in one split}
      \item{pm.level:}{cutoff level of permutation test}
      \item{pm.repeat:}{size of permuation test}
      \item{ord.bin:}{the maximium number of bins for ordinal/continuous covariates}
    }
    Optional optAUC arguments:
    \describe{
    \item{method.subgrp:}{method for subgroup identification using optAUC: "BATTing" or "youden".}
    \item{scale:}{logical variable indicating if scaling is required.}
    \item{xselect:}{a logical flag. TRUE: use LASSO for variable selection;
    FALSE: use all predictors to calculate the composite score. Default = FALSE.}
    \item{cv.iter.xselect:}{the number of iterations of cross-validation used
    for xselect.}
    \item{k.fold.xselect:}{the number of folds of cross-validation used for
    xselect.}
    \item{method.xselect:}{method used for variable selection: CV (default)
    cross validated AUC/ABC; "aAUC" approximate AUC (only for prognostic case).}
    }
  }
}
\details{
The function contains four algorithms for developing threshold-based multivariate (prognostic/predictive) biomarker signatures via resampled tree-based algorithms (Sequential BATTing),  Monte-Carlo variations of the Adaptive Indexing method (AIM and AIM-Rule), and Patient Rule Induction Method
. Variable selection is automatically built-in to these algorithms. Final signatures are returned with interaction plots for predictive signatures. Cross-validation performance evaluation and testing dataset results are also output.
}
\value{

  \item{res}{list of all results from the algorithm}
  \item{train.stat}{list of subgroup statistics on training dataset}
  \item{test.stat}{list of subgroup statistics on testing dataset}
  \item{cv.res}{list of all results from cross-validation on training dataset}
  \item{train.plot}{interaction plot for training dataset}
  \item{test.plot}{interaction plot for testing dataset}
}
\references{
  Huang X. et al. (2016)
  Biomarker Signatures for Patient Subgroup Selection in Clinical Drug Development. \emph{Statistics in Medicine} (Submitted)

  Chen G. et al. (2015)
  A PRIM approach to predictive-signature development for patient stratification
  \emph{Applied Statistics}, \bold{34}, doi: 10.1002/sim.6343.
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
data(Sepsis.test)

yvar="survival"
xvars=names(Sepsis.train)[2:12]
trtvar="THERAPY"

set.seed(123)
subgrp <- SubgrpID(data.train=Sepsis.train,
                  data.test=Sepsis.test,
                  yvar=yvar,
                  trtvar=trtvar,
                  trtref="active",
                  xvars=xvars,
                  type="b",
                  des.res = "smaller",
                  method="AIM.Rule")
subgrp$res
subgrp$train.stat
subgrp$test.stat
subgrp$train.plot
subgrp$test.plot
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ BATTing }
\keyword{ Sequential-BATTing }
\keyword{ AIM-Rule }
\keyword{ PRIM }

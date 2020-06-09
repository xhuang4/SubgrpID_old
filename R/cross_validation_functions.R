###################################################################################################
#
# Cross-Validation Functions
# Author: Paul Trow
# Last revised: 10/25/14
#
# This file contains general cross-validation code and performance analysis functions.
#
#   1. Cross-validation code
#   2. Performance analysis functions
#
###################################################################################################

###################################################################################################
#
# Cross-validation code
#
###################################################################################################

#' Perform k-fold cross-validation of a model.
#'
#' @param data the CV data
#' @param model.Rfunc Name of the model function.
#' @param model.Rfunc.args List of input arguments to model.Rfunc.
#' @param predict.Rfunc Name of the prediction function, which takes the prediction rule returned by model.Rfunc along with any input data (not necessarily the input data to kfold.cv) and returns a TRUE-FALSE predictionvector specifying the positive and negative classes for the data.
#' @param predict.Rfunc.args List containing input arguments to predict.Rfunc, except for data and predict.rule.
#' @param k.fold Number of folds of the cross-validation.
#' @param cv.iter Number of iterations of the cross-validation. If model.Rfunc returns an error at any of the k.fold calls, the current iteration is aborted. Iterations are repeated until cv.iter successful iterations have occurred.
#' @param strata Stratification vector of length the number of rows of data, usually corresponding to the vector of events.
#' @param max.iter Function stops after max.iter iterations even if cv.iter successful iterations have not occurred.
#'
#' @return List of length 2 with the following fields:
#' @return cv.data - List of length cv.iter. Entry i contains the output of predict.Rfunc at the ith iteration.
#' @return sig.list - list of length cv.iter * k.fold, whose entries are the prediction.rules (signatures) returned by model.Rfunc at each k.fold iteration.
#'
kfold.cv <- function(data, model.Rfunc, model.Rfunc.args, predict.Rfunc, predict.Rfunc.args, k.fold=5, cv.iter=50, strata, max.iter=500) {
  ###### Main kfold.cv Code #################

    n.data <- nrow(data)
    data.index <-  1:n.data

    if(is.null(strata)){
        strata <- rep(1, n.data)
    }

    cv.data <- list()
    sig.list <- list()
    predict.func.errors <- list()
    model.func.errors <- list()
    iter.count <- 1
    # Number of iterations
    iter.success.count <- 0
    # Number of successful iterations (all folds successful)

    while (iter.success.count < cv.iter){
        # Continue until cv.iter successful iterations have been performed (or max.iter reached).

        if (iter.count > max.iter) {
            cat("Maximum number of iterations reached - ending cross-validation.")
            break
        }

        cat("CV iteration ", iter.count, "(Successes: ", iter.success.count, ")\n\n")
        cv.index <- balanced.folds(strata, k.fold)
        train.index <- lapply(cv.index, function(x) setdiff(data.index, x))
        # Indices for training data for each of the k folds.
        cv.vec <- c()
        sig.list.short <- list()
        fold.success.count <- 0
        # Number of successful fold evaluations.

        for (j in 1:k.fold) {
            cat("Fold ", j)
            predict.rule <- try(eval(call(model.Rfunc, data=data[train.index[[j]], ], args=model.Rfunc.args)), silent=TRUE)
            # Apply model.Rfunc to the training data to get a prediction rule.

            if(class(predict.rule)!="try-error") {
                sig.list.short <- c(sig.list.short, list(predict.rule))
                cv.temp <- try(eval(call(predict.Rfunc, data=data[cv.index[[j]], ], predict.rule=predict.rule, args=predict.Rfunc.args)), silent=TRUE)
                # Apply prediction function with prediction rule to the left-out fold of data.

                if (class(cv.temp)!="try-error") {
                    fold.vec <- rep(j, nrow(cv.temp))
                    # Record fold number for current subset of data.
                    cv.temp <- cbind(cv.temp, data.frame("fold"=fold.vec))
                    cv.vec <- rbind(cv.vec, cv.temp)
                    # Append the output of prediction function to cv.vec.
                    fold.success.count <- fold.success.count + 1
                } else {
                    # Error in predict.Rfunc, so quit current iteration.
                    error.tmp <- list(cv.temp)
                    names(error.tmp) <- paste("iter ", iter.count, " fold ", j, sep="")
                    predict.func.errors <- c(predict.func.errors, error.tmp)
                    cat(" ", cv.temp, "\n")
                    break;
                }

                cat("\n")
            } else {
                # Error in model.Rfunc, so quit current iteration.
                error.tmp <- list(predict.rule)
                names(error.tmp) <- paste("iter ", iter.count, " fold ", j, sep="")
                model.func.errors <- c(model.func.errors, error.tmp)
                cat(" ", model.Rfunc, ": ", predict.rule, "\n")
                break;
            }
        }

        iter.count <- iter.count + 1

        if (fold.success.count == k.fold) {
            # Model and prediction functions were successful at all folds.
            # Put cv.vec back in original order.
            cv.vec$row.numbers <- as.numeric(rownames(cv.vec))
            cv.vec <- cv.vec[order(cv.vec$row.numbers), ]
            cv.vec <- subset(cv.vec, select=setdiff(names(cv.vec), "row.numbers"))

            cv.data <- c(cv.data, list(cv.vec))
            # Append predictions for current iteration to cv.data
            iter.success.count <- iter.success.count + 1
            sig.list.tmp <- list(sig.list.short)
            names(sig.list.tmp) <- paste("Iteration ", iter.success.count, sep="")
            sig.list <- c(sig.list, sig.list.tmp)
        }
    }

    error.log <- list(model.func.errors=model.func.errors, predict.func.errors=predict.func.errors)
    results <- list(sig.list=sig.list, cv.data=cv.data, error.log=error.log)
    results
}

###################################################################################################
#
#  Cross-Validation Support Functions
#
###################################################################################################

#' Create balanced folds for cross-validation.
#'
#' @param y the response vector
#' @param nfolds number of folds
#'
#' @return This function returns balanced folds
#'
balanced.folds <- function(y, nfolds = min(min(table(y)),10)) {
    #
    # Create balanced folds for cross-validation.
    #
    y[is.na(y)] <- resample(y[!is.na(y)], size = sum(is.na(y)),
                            replace = TRUE)
    totals <- table(y)

    if (length(totals) < 2) {
        return(cv.folds(length(y), nfolds))
    } else {
        fmax <- max(totals)
        nfolds <- min(nfolds, fmax)
        nfolds <- max(nfolds, 2)
        yids <- split(seq(y), y)
        bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds,
                         length(totals))

        for (i in seq(totals)) {

            if (length(yids[[i]]) > 1) {
                bigmat
                bigmat[seq(totals[i]), i] <- sample(yids[[i]])
            }

            if (length(yids[[i]]) == 1) {
                bigmat[seq(totals[i]), i] <- yids[[i]]
            }
        }

        smallmat <- matrix(bigmat, nrow = nfolds)
        smallmat <- permute.rows(t(smallmat))
        res <- vector("list", nfolds)

        for (j in 1:nfolds) {
            jj <- !is.na(smallmat[, j])
            res[[j]] <- smallmat[jj, j]
        }
        return(res)
    }
}

#' Creates a permutation of given size.
#'
#' @param x the x vector.
#' @param size resampling size.
#' @param ... optional argument.
#'
#' @return A resample of x is returned.
#'
resample <- function(x, size, ...) {

    if (length(x) <= 1) {
        if (!missing(size) && size == 0)
            x[FALSE]
        else x
    } else {
        sample(x, size, ...)
    }
}

#' Cross-validation folds.
#'
#' @param n number of observations.
#' @param folds number of folds.
#'
#' @return a list containing the observation numbers for each fold.
#'
cv.folds <- function(n, folds = 10) {
    split(resample(1:n), rep(1:folds, length = n))
}


#' Ramdomly permute the rows of a matrix.
#'
#' @param A a matrix for which its rows have to be permuted.
#'
#' @return the matrix with permuted rows.
#'
permute.rows <- function(A) {
    B <- t(A)
    row.list <- split(B, rep(1:ncol(B), each = nrow(B)))
    # Convert A into list of rows.
    t(sapply(row.list, permute.vector))
}

#' Randomly permute the entries of a vector.

#'
#' @param x the vector for which its entries have to be permuted
#'
#' @return the permuted vector
#'
permute.vector <- function(x) {
    x[sample(1:length(x))]
}

#' Create a list of variables corresponding to the arguments of the function func.name and assigns values.
#'
#' @param func.name function name
#'
#' @return list of variables corresponding to the arguments of the function
#'
make.arg.list <- function(func.name) {
    arg.list <- list()
    arg.names <- names(as.list(args(func.name)))

    for (i in 1:length(arg.names)) {
        arg.val <- eval(parse(text=arg.names[i]))

        if (!is.null(arg.val)) {
            arg.tmp <- list(arg.val)
            names(arg.tmp) <- arg.names[i]
            arg.list <- c(arg.list, arg.tmp)
        }
    }

    arg.list
}

###################################################################################################
#
#  Cross-validation Performance Evaluation
#
###################################################################################################

# library(Zelig)

#' Cross-validation Performance Evaluation
#' @description  Take the raw output of kfold.cv and calculate performance statistics for each iteration of the cross-validation.
#'
#' @param cv.data output of prediction function from kfold.cv
#' @param y data frame of the response variable from CV data.
#' @param censor.vec data frame indicating censoring for survival data. For binary or continuous data, set censor.vec <- NULL.
#' @param trt.vec data frame indicating whether or not the patient was treated. For the pronostic case, set trt.vec <- NULL.
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c"
#'
#' @return a list containing raw statistics and fold information
#'
evaluate.cv.results <- function(cv.data, y, censor.vec, trt.vec, type) {

  pred.class<-fold<-NULL

    pvals <- data.frame()
    data <- y
    yvar <- names(y)

    if (type=="s" | type=="b") {
        ratios <- data.frame()
    }

    if ((type=="b")&is.null(trt.vec)) {
        bin.stats <- data.frame()
    }

    if (!is.null(censor.vec)) {
        censorvar <- names(censor.vec)
        data <- cbind(data, censor.vec)
    } else {
        censorvar <- NULL
    }

    if (!is.null(trt.vec)) {
        trtvar <- names(trt.vec)
        data <- cbind(data, trt.vec)
        sigpos.trt.stats <- data.frame()
        sigpos.ctrl.stats <- data.frame()
        signeg.trt.stats <- data.frame()
        signeg.ctrl.stats <- data.frame()
    } else {
        trtvar <- NULL
        sigpos.stats <- data.frame()
        signeg.stats <- data.frame()
    }

    if (!is.null(cv.data[[1]]$pred.class)&!is.null(cv.data[[1]]$fold)) {
        pred.classes <- subset(cv.data[[1]], select=pred.class)
        folds <- subset(cv.data[[1]], select=fold)
    } else {
        pred.classes <- NULL
        folds <- NULL
    }


    for (i in 1:length(cv.data)) {

        if (!is.null(cv.data[[i]]$pred.class)& !is.null(cv.data[[i]]$fold)) {

            pred.class.tmp <- subset(cv.data[[i]], select=pred.class)

            if (i > 1) {
                pred.classes <- cbind(pred.classes, pred.class.tmp)
                fold.tmp <- subset(cv.data[[i]], select=fold)
                folds <- cbind(folds, fold.tmp)
            }

            data.tmp <- cbind(data, pred.class.tmp)
            stats.tmp <- cv.pval(yvar=yvar, censorvar=censorvar, trtvar=trtvar, data=data.tmp, type=type)
            pvals <- rbind(pvals, stats.tmp$pval)

            if (type=="s") {
                ratios <- rbind(ratios, stats.tmp$hazard.ratios)
            } else if (type=="b") {
                ratios <- rbind(ratios, stats.tmp$odds.ratios)
            }

            if ((type=="b")&is.null(trtvar)) {
                bin.stats <- rbind(bin.stats, stats.tmp$bin.stats)
            }


            if (is.null(trtvar)) {
                group.stats <- find.prog.stats(data=data.tmp, yvar=yvar, type=type, censorvar=censorvar)
                sigpos.stats <- rbind(sigpos.stats, group.stats["sig+", ])
                signeg.stats <- rbind(signeg.stats, group.stats["sig-", ])
            } else {
                group.stats <- find.pred.stats(data=data.tmp, yvar=yvar, trtvar=trtvar, type=type, censorvar=censorvar)
                sigpos.trt.stats <- rbind(sigpos.trt.stats, group.stats["sig+.trt", ])
                sigpos.ctrl.stats <- rbind(sigpos.ctrl.stats, group.stats["sig+.ctrl", ])
                signeg.trt.stats <- rbind(signeg.trt.stats, group.stats["sig-.trt", ])
                signeg.ctrl.stats <- rbind(signeg.ctrl.stats, group.stats["sig-.ctrl", ])
            }
        }
    }

    if (type=="s" | type=="b") {
        raw.stats <- list(pvals=pvals, ratios=ratios)
    } else {
        raw.stats <- list(pvals=pvals)
    }

    if (!is.null(trtvar)) {
        # Some minor housekeeping - row names clutter up the output.
        row.names(sigpos.trt.stats) <- NULL
        row.names(sigpos.ctrl.stats) <- NULL
        row.names(signeg.trt.stats) <- NULL
        row.names(signeg.ctrl.stats) <- NULL

        raw.stats <- c(raw.stats, list(sigpos.trt.stats=sigpos.trt.stats,
                                       sigpos.ctrl.stats=sigpos.ctrl.stats,
                                       signeg.trt.stats=signeg.trt.stats,
                                       signeg.ctrl.stats=signeg.ctrl.stats))
    } else {
        # Toss out the row names
        row.names(sigpos.stats) <- NULL
        row.names(signeg.stats) <- NULL
        raw.stats <- c(raw.stats, list(sigpos=sigpos.stats,
                                       signeg=signeg.stats))

        if (type=="b") {
            raw.stats <- c(raw.stats, list(bin.stats=bin.stats))
        }

    }

    stats <- list(raw.stats=raw.stats, pred.classes=pred.classes,folds=folds)

}

#' Calculate summary statistics from raw statistics returned by evaluate.cv.results.
#'
#' @param raw.stats raw statistics from evaluate.cv.results
#' @param trtvar treatment variable name
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c"
#'
#' @return a list containing p-values, summary statistics and group statistics.
#'
summarize.cv.stats <- function(raw.stats, trtvar, type) {
    #
    #
    #

    if (!is.null(trtvar)) {
        # Predictive case
        med.pval.inter = median(raw.stats$pvals[, "interaction"], na.rm=TRUE)
        ind.med.inter = which.min(abs(raw.stats$pvals[, "interaction"]-med.pval.inter))

        if (type=="s") {
            summary.pvals = raw.stats$pvals[ind.med.inter,]
            summary.ratios = raw.stats$ratios[ind.med.inter,]
            med.sigpos.trt = raw.stats$sigpos.trt.stats[ind.med.inter,]
            med.sigpos.ctrl = raw.stats$sigpos.ctrl.stats[ind.med.inter,]
            med.signeg.trt = raw.stats$signeg.trt.stats[ind.med.inter,]
            med.signeg.ctrl = raw.stats$signeg.ctrl.stats[ind.med.inter,]

            mad.names=paste("mad.",names(raw.stats$pvals),sep="")
            summary.pvals[,mad.names]=apply(raw.stats$pvals,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$ratios),sep="")
            summary.ratios[,mad.names]=apply(raw.stats$ratios,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$sigpos.trt.stats),sep="")
            med.sigpos.trt[,mad.names]=apply(raw.stats$sigpos.trt.stats,2,mad,na.rm=TRUE)
            med.sigpos.ctrl[,mad.names]=apply(raw.stats$sigpos.ctrl.stats,2,mad,na.rm=TRUE)
            med.signeg.trt[,mad.names]=apply(raw.stats$signeg.trt.stats,2,mad,na.rm=TRUE)
            med.signeg.ctrl[,mad.names]=apply(raw.stats$signeg.ctrl.stats,2,mad,na.rm=TRUE)

        } else if (type=="b") {
            summary.pvals = raw.stats$pvals[ind.med.inter,]
            summary.ratios = raw.stats$ratios[ind.med.inter,]
            med.sigpos.trt = raw.stats$sigpos.trt.stats[ind.med.inter,]
            med.sigpos.ctrl = raw.stats$sigpos.ctrl.stats[ind.med.inter,]
            med.signeg.trt = raw.stats$signeg.trt.stats[ind.med.inter,]
            med.signeg.ctrl = raw.stats$signeg.ctrl.stats[ind.med.inter,]

            mad.names=paste("mad.",names(raw.stats$pvals),sep="")
            summary.pvals[,mad.names]=apply(raw.stats$pvals,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$ratios),sep="")
            summary.ratios[,mad.names]=apply(raw.stats$ratios,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$sigpos.trt.stats),sep="")
            med.sigpos.trt[,mad.names]=apply(raw.stats$sigpos.trt.stats,2,mad,na.rm=TRUE)
            med.sigpos.ctrl[,mad.names]=apply(raw.stats$sigpos.ctrl.stats,2,mad,na.rm=TRUE)
            med.signeg.trt[,mad.names]=apply(raw.stats$signeg.trt.stats,2,mad,na.rm=TRUE)
            med.signeg.ctrl[,mad.names]=apply(raw.stats$signeg.ctrl.stats,2,mad,na.rm=TRUE)

        } else {
            # type="c"
            summary.pvals = raw.stats$pvals[ind.med.inter,]
            summary.ratios=data.frame("ratios"=NA)
            med.sigpos.trt = raw.stats$sigpos.trt.stats[ind.med.inter,]
            med.sigpos.ctrl = raw.stats$sigpos.ctrl.stats[ind.med.inter,]
            med.signeg.trt = raw.stats$signeg.trt.stats[ind.med.inter,]
            med.signeg.ctrl = raw.stats$signeg.ctrl.stats[ind.med.inter,]

            mad.names=paste("mad.",names(raw.stats$pvals),sep="")
            summary.pvals[,mad.names]=apply(raw.stats$pvals,2,mad,na.rm=TRUE)

            summary.ratios[,"mad.ratios"]=NA

            mad.names=paste("mad.",names(raw.stats$sigpos.trt.stats),sep="")
            med.sigpos.trt[,mad.names]=apply(raw.stats$sigpos.trt.stats,2,mad,na.rm=TRUE)
            med.sigpos.ctrl[,mad.names]=apply(raw.stats$sigpos.ctrl.stats,2,mad,na.rm=TRUE)
            med.signeg.trt[,mad.names]=apply(raw.stats$signeg.trt.stats,2,mad,na.rm=TRUE)
            med.signeg.ctrl[,mad.names]=apply(raw.stats$signeg.ctrl.stats,2,mad,na.rm=TRUE)

        }
        summary.group.stats <- rbind(med.sigpos.trt, med.sigpos.ctrl, med.signeg.trt, med.signeg.ctrl)
        row.names(summary.group.stats)=c("sig+.trt", "sig+.ctrl", "sig-.trt", "sig-.ctrl")
        summary <- list(pvals=summary.pvals, ratios=summary.ratios, group.stats=summary.group.stats)
    } else {
        # Prognostic case
        med.pval = median(raw.stats$pvals[, "pval"], na.rm=TRUE)
        ind.med = which.min(abs(raw.stats$pvals[, "pval"]-med.pval))

        if (type=="s") {
            summary.pvals = subset(raw.stats$pvals, subset=(1:dim(raw.stats$pvals)[1])==ind.med)
            summary.ratios = subset(raw.stats$ratios, subset=(1:dim(raw.stats$ratios)[1])==ind.med)
            med.sigpos = raw.stats$sigpos[ind.med,]
            med.signeg = raw.stats$signeg[ind.med,]

            mad.names=paste("mad.",names(raw.stats$pvals),sep="")
            summary.pvals[,mad.names]=apply(raw.stats$pvals,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$ratios),sep="")
            summary.ratios[,mad.names]=apply(raw.stats$ratios,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$sigpos),sep="")
            med.sigpos[,mad.names]=apply(raw.stats$sigpos,2,mad,na.rm=TRUE)
            med.signeg[,mad.names]=apply(raw.stats$signeg,2,mad,na.rm=TRUE)

            summary.group.stats <- rbind(med.sigpos,med.signeg)
            row.names(summary.group.stats)=c("sig+","sig-")
            summary <- list(pvals=summary.pvals, ratios=summary.ratios, group.stats=summary.group.stats)


        } else if (type=="b") {

            summary.pvals = subset(raw.stats$pvals, subset=(1:dim(raw.stats$pvals)[1])==ind.med)
            summary.ratios = subset(raw.stats$ratios, subset=(1:dim(raw.stats$ratios)[1])==ind.med)
            summary.bin.stats = subset(raw.stats$bin.stats, subset=(1:dim(raw.stats$bin.stats)[1])==ind.med)
            med.sigpos = raw.stats$sigpos[ind.med,]
            med.signeg = raw.stats$signeg[ind.med,]

            mad.names=paste("mad.",names(raw.stats$pvals),sep="")
            summary.pvals[,mad.names]=apply(raw.stats$pvals,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$ratios),sep="")
            summary.ratios[,mad.names]=apply(raw.stats$ratios,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$bin.stats),sep="")
            summary.bin.stats[,mad.names]=apply(raw.stats$bin.stats,2,mad,na.rm=TRUE)

            mad.names=paste("mad.",names(raw.stats$sigpos),sep="")
            med.sigpos[,mad.names]=apply(raw.stats$sigpos,2,mad,na.rm=TRUE)
            med.signeg[,mad.names]=apply(raw.stats$signeg,2,mad,na.rm=TRUE)

            summary.group.stats <- rbind(med.sigpos,med.signeg)
            row.names(summary.group.stats)=c("sig+","sig-")
            summary <- list(pvals=summary.pvals, ratios=summary.ratios, group.stats=summary.group.stats, bin.stats=summary.bin.stats)

        } else {
            summary.pvals = subset(raw.stats$pvals, subset=(1:dim(raw.stats$pvals)[1])==ind.med)
            summary.ratios = data.frame("ratios"=NA)
            med.sigpos = raw.stats$sigpos[ind.med,]
            med.signeg = raw.stats$signeg[ind.med,]

            mad.names=paste("mad.",names(raw.stats$pvals),sep="")
            summary.pvals[,mad.names]=apply(raw.stats$pvals,2,mad,na.rm=TRUE)

            mad.names="mad.ratios"
            summary.ratios[,mad.names]=NA

            mad.names=paste("mad.",names(raw.stats$sigpos),sep="")
            med.sigpos[,mad.names]=apply(raw.stats$sigpos,2,mad,na.rm=TRUE)
            med.signeg[,mad.names]=apply(raw.stats$signeg,2,mad,na.rm=TRUE)

            summary.group.stats <- rbind(med.sigpos,med.signeg)
            row.names(summary.group.stats)=c("sig+","sig-")
            summary <- list(pvals=summary.pvals, ratios=summary.ratios, group.stats=summary.group.stats)


        }
    }

}

#' Get statistics for a single set of predictions.
#'
#' @param y data frame of the response variable.
#' @param predict.data output of prediction function from kfold.cv.
#' @param censor.vec data frame indicating censoring for survival data. For binary or continuous data, set censor.vec <- NULL.
#' @param trt.vec data frame indicating whether or not the patient was treated. For the pronostic case, set trt.vec <- NULL.
#' @param trtref treatment reference.
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c".
#'
#' @return a list containing p-value and group statistics.
#'
evaluate.results <- function(y, predict.data, censor.vec=NULL, trt.vec=NULL, trtref=NULL, type) {
    y <- as.data.frame(y)
    predict.data <- as.data.frame(predict.data)

    if (nrow(predict.data) != nrow(y)) {
        stop("Error: prediction data and y have different numbers of rows.")
    }

    yvar <- names(y)
    data <- cbind(y, predict.data)

    if (!is.null(censor.vec)) {
        censor.vec <- as.data.frame(censor.vec)
        censorvar <- names(censor.vec)

        if (nrow(censor.vec) != nrow(y)) {
            stop("Error: censor.vec and y have different numbers of rows.")
        }

        data <- cbind(data, censor.vec)
    } else {
        censorvar <- NULL
    }

    if (!is.null(trt.vec)) {
        trt.vec <- as.data.frame(trt.vec)
        trtvar <- names(trt.vec)
        Trt <- trt.vec[, trtvar]

        if (!is.null(trtref)) {
            # If trtref is provided, create 0-1 vector for trt1 in which 1 corresponds to trtref.
            trt1<-rep(0,length=length(Trt))
            trt1[Trt==trtref]<-1
            trt.vec <- data.frame(trt1)
            names(trt.vec) <- trtvar
            trtref <- NULL
        }

        data <- cbind(data, trt.vec)
    } else {
        trtvar <- NULL
    }

    stats.tmp <- cv.pval(yvar=yvar, censorvar=censorvar, trtvar=trtvar, data=data, type=type)
    pval <- stats.tmp$pval

    if (type=="s") {
        ratios <- stats.tmp$hazard.ratios
    } else if (type=="b") {
        ratios <- stats.tmp$odds.ratios
    }

    if ((type=="b")&is.null(trtvar)) {
        bin.stats <- stats.tmp$bin.stats
    }


    if (is.null(trtvar)) {
        group.stats <- find.prog.stats(data=data, yvar=yvar, type=type, censorvar=censorvar)
    } else {
        group.stats <- find.pred.stats(data=data, yvar=yvar, trtvar=trtvar, type=type, censorvar=censorvar)
    }

    if (type=="b") {
        stats <- list(pval=pval, ratios=ratios, group.stats=group.stats)
        if(is.null(trtvar)) stats=c(stats, list(bin.stats=bin.stats))

    } else if (type=="s") {
        stats <- list(pval=pval, ratios=ratios, group.stats=group.stats)
    } else {
        stats <- list(pval=pval, group.stats=group.stats)
    }

    stats
}

#' Find predictive stats from response and prediction vector
#'
#' @param data data frame with response and prediction vector
#' @param yvar response variable name
#' @param trtvar treatment variable name
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c".
#' @param censorvar censoring variable name
#'
#' @return a data frame of predictive statistics
#'
find.pred.stats <- function(data, yvar, trtvar, type, censorvar){

    sigpos.trt <- data[data$pred.class==TRUE & data[, trtvar]==1,]
    sigpos.ctrl <- data[data$pred.class==TRUE & data[, trtvar]==0,]
    signeg.trt <- data[data$pred.class==FALSE & data[, trtvar]==1,]
    signeg.ctrl <- data[data$pred.class==FALSE & data[, trtvar]==0,]

    n.sigpos.trt <- nrow(sigpos.trt)
    n.sigpos.ctrl <- nrow(sigpos.ctrl)
    n.signeg.trt <- nrow(signeg.trt)
    n.signeg.ctrl <- nrow(signeg.ctrl)

    n.subj <- c(n.sigpos.trt, n.sigpos.ctrl, n.signeg.trt, n.signeg.ctrl)

    if (type=="c") {
        y.median <- c(median(sigpos.trt[, yvar]), median(sigpos.ctrl[, yvar]), median(signeg.trt[, yvar]), median(signeg.ctrl[, yvar]))
        y.mean <- c(mean(sigpos.trt[, yvar]), mean(sigpos.ctrl[, yvar]), mean(signeg.trt[, yvar]), mean(signeg.ctrl[, yvar]))
        y.sd <- c(sd(sigpos.trt[, yvar]), sd(sigpos.ctrl[, yvar]), sd(signeg.trt[, yvar]), sd(signeg.ctrl[, yvar]))
        group.stats <- data.frame(n=n.subj, median=y.median, mean=y.mean, sd=y.sd, row.names=c("sig+.trt", "sig+.ctrl", "sig-.trt", "sig-.ctrl"))
    }

    if (type=="s") {
        model.text <- paste("Surv(", yvar, ",", censorvar, ") ~ pred.class+", trtvar, sep="")
        model.formula <- as.formula(model.text)

        fit=survfit(model.formula,data=data)

        mean.surv.time=(summary(fit,rmean="common")$table)[,"*rmean"][4:1]
        se.mean.surv=(summary(fit,rmean="common")$table)[,"*se(rmean)"][4:1]
        median.surv.time=(summary(fit,rmean="common")$table)[,"median"][4:1]


        model.text <- paste("Surv(", yvar, ",", censorvar, ") ~ factor(", trtvar, ")*factor(pred.class)", sep="")
        model.formula <- as.formula(model.text)
        model <- coxph(model.formula, data=data)
        coefs <- coef(model)
        data.pos <- rbind(c(1, sum(data[, trtvar]==1 & data$pred.class==TRUE), sum(coefs[1:3])),
                          c(0, sum(data[, trtvar]==0 & data$pred.class==TRUE), coefs[2])
        )
        colnames(data.pos) <- c("trt", "n", "log.hazard.ratio")
        data.neg <- rbind(c(1, sum(data[, trtvar]==1 & data$pred.class==FALSE), coefs[1]),
                          c(0, sum(data[, trtvar]==0 & data$pred.class==FALSE), 0)
        )
        colnames(data.neg) <- c("trt", "n", "log.hazard.ratio")
        data.pos[,"log.hazard.ratio"] <- exp(data.pos[,"log.hazard.ratio"])
        data.neg[,"log.hazard.ratio"] <- exp(data.neg[,"log.hazard.ratio"])
        hazard.ratio <- c(data.pos[1, "log.hazard.ratio"], data.pos[2, "log.hazard.ratio"], data.neg[1, "log.hazard.ratio"], data.neg[2, "log.hazard.ratio"])
        group.stats <- data.frame("n"=n.subj, "mean.surv.time"=mean.surv.time, "se.mean.surv"=se.mean.surv, "median.surv.time"=median.surv.time, "hazard ratio" = hazard.ratio,
                                  row.names=c("sig+.trt", "sig+.ctrl", "sig-.trt", "sig-.ctrl"))
    }

    if (type=="b"){
        pos.trt.events <- sum(sigpos.trt[, yvar]==1)
        pos.ctrl.events <- sum(sigpos.ctrl[, yvar]==1)
        neg.trt.events <- sum(signeg.trt[, yvar]==1)
        neg.ctrl.events <- sum(signeg.ctrl[, yvar]==1)
        resp.rate <- c(pos.trt.events/n.sigpos.trt, pos.ctrl.events/n.sigpos.ctrl, neg.trt.events/n.signeg.trt, neg.ctrl.events/n.signeg.ctrl)

        group.stats <- data.frame(n=n.subj)
        model.text <- paste(yvar, "~factor(", trtvar, ") * factor(pred.class)", sep="")
        model.formula <- as.formula(model.text)
        model <- glm(model.formula, family=binomial(link="logit"), data=data)
        coefs <- coef(model)
        data.pos <- rbind(c(1, sum(data[, trtvar]==1 & data$pred.class==TRUE), sum(coefs[2:4])),
                          c(0, sum(data[, trtvar]==0 & data$pred.class==TRUE), coefs[3])
        )
        colnames(data.pos) <- c("trt", "n", "log.odds.ratio")
        data.neg <- rbind(c(1, sum(data$trt==1 & data$pred.class==FALSE), coefs[2]),
                          c(0, sum(data$trt==0 & data$pred.class==FALSE), 0)
        )
        colnames(data.neg) <- c("trt", "n", "log.odds.ratio")
        data.pos[,"log.odds.ratio"] <- exp(data.pos[,"log.odds.ratio"])
        data.neg[,"log.odds.ratio"] <- exp(data.neg[,"log.odds.ratio"])
        odds.ratio <- c(data.pos[1, "log.odds.ratio"], data.pos[2, "log.odds.ratio"], data.neg[1, "log.odds.ratio"], data.neg[2, "log.odds.ratio"])
        group.stats <- data.frame("n"=n.subj, "resp.rate"=resp.rate, "odds.ratio"=odds.ratio,
                                  row.names=c("sig+.trt", "sig+.ctrl", "sig-.trt", "sig-.ctrl"))
    }

    return(group.stats)
}

#' Find prognostic stats from response and prediction vector
#'
#' @param data data frame with response and prediction vector
#' @param yvar response variable name
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c".
#' @param censorvar censoring variable name
#'
#' @return a data frame of predictive statistics
#'
find.prog.stats <- function(data, yvar, type, censorvar) {
    # Find prognostic stats from prediction vector pred, in which neg. pred = 0, pos. pred. = 1
    # events is a column of outcomes (0 or 1)
    #
    sigpos <- data[data$pred.class==TRUE,]
    signeg <- data[data$pred.class==FALSE,]
    n.sigpos <- nrow(sigpos)
    n.signeg <- nrow(signeg)
    n.subj <- c(n.sigpos, n.signeg)


    if (type=="c") {
        y.median <- c(median(sigpos[, yvar]), median(signeg[, yvar]))
        y.mean <- c(mean(sigpos[, yvar]), mean(signeg[, yvar]))
        y.sd <- c(sd(sigpos[, yvar]), sd(signeg[, yvar]))
        group.stats <- data.frame(n=n.subj, median=y.median, mean=y.mean, sd=y.sd, row.names=c("sig+", "sig-"))
    }

    if (type=="s") {
        model.text <- paste("Surv(", yvar, ",", censorvar, ")~ pred.class", sep="")
        model.formula <- as.formula(model.text)

        fit=survfit(model.formula,data=data)

        mean.surv.time=(summary(fit,rmean="common")$table)[,"*rmean"][2:1]
        se.mean.surv=(summary(fit,rmean="common")$table)[,"*se(rmean)"][2:1]
        median.surv.time=(summary(fit,rmean="common")$table)[,"median"][2:1]


        model.text <- paste("Surv(", yvar, ",", censorvar, ") ~ pred.class", sep="")
        model.formula <- as.formula(model.text)
        model <- coxph(model.formula, data=data)
        coefs <- coef(model)
        data.pos <- c(sum(data$pred.class==TRUE), coefs[[1]])
        names(data.pos) <- c("n", "log.hazard.ratio")
        data.neg <- c(sum(data$pred.class==FALSE), 0)
        names(data.neg) <- c("n", "log.hazard.ratio")
        data.pos["log.hazard.ratio"] <- exp(data.pos["log.hazard.ratio"])
        data.neg["log.hazard.ratio"] <- exp(data.neg["log.hazard.ratio"])
        hazard.ratio <- c(data.pos["log.hazard.ratio"], data.neg["log.hazard.ratio"])
        group.stats <- data.frame("n"=n.subj, "mean.surv.time"=mean.surv.time, "se.mean.surv"=se.mean.surv, "median.surv.time"=median.surv.time, "hazard ratio" = hazard.ratio,
                                  row.names=c("sig+", "sig-"))
    }
    if (type=="b") {

        pos.events <- sum(sigpos[, yvar]==1)
        neg.events <- sum(signeg[, yvar]==1)
        resp.rate <- c(pos.events/n.sigpos, neg.events/n.signeg)

        group.stats <- data.frame(n=n.subj)
        model.text <- paste(yvar, "~pred.class", sep="")
        model.formula <- as.formula(model.text)
        model <- glm(model.formula, family=binomial(link="logit"), data=data)
        coefs <- coef(model)
        data.pos <- c(sum(data$pred.class==TRUE), coefs[[2]])
        names(data.pos) <- c("n", "log.odds.ratio")
        data.neg <- c(sum(data$pred.class==FALSE), 0)
        names(data.neg) <- c("n", "log.odds.ratio")
        data.pos["log.odds.ratio"] <- exp(data.pos["log.odds.ratio"])
        data.neg["log.odds.ratio"] <- exp(data.neg["log.odds.ratio"])
        odds.ratio <- c(data.pos["log.odds.ratio"],data.neg["log.odds.ratio"])

        group.stats <- data.frame("n"=n.subj, "resp.rate"=resp.rate, "odds.ratio" = odds.ratio,
                                  row.names=c("sig+", "sig-"))
    }

    return(group.stats)
}


#' p-value calculation for each iteration of cross validation.
#'
#' @param yvar response variable name.
#' @param censorvar censor-variable name.
#' @param trtvar treatment variable name. For prognostic case trtvar=NULL.
#' @param data dataset containing response and predicted output.
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c".
#'
#' @return p-value based on response and prediction vector for each iteration.
#'
cv.pval <- function(yvar, censorvar=NULL, trtvar=NULL, data, type="s"){
    pred.class<-NULL #check
    if (type=="s") {
        if (!is.null(trtvar)) {
            model.text1 <- paste("Surv(", yvar, ",", censorvar, ")~", trtvar, sep="")
            model.formula1 <- as.formula(model.text1)
            model.text2 <- paste("Surv(", yvar, ",", censorvar, ")~pred.class", sep="")
            model.formula2 <- as.formula(model.text2)
            model.text3 <- paste("Surv(", yvar, ",", censorvar, ")~", trtvar , "*pred.class", sep="")
            model.formula3 <- as.formula(model.text3)

            model.temp1.1 <- try(coxph(model.formula1, data=data, subset=pred.class))
            if (class(model.temp1.1)[1] != "try-error") {

                trt.diff.pos.gp <- try(summary(model.temp1.1)$coefficients[trtvar, "Pr(>|z|)"])
                if (class(trt.diff.pos.gp)[1] == "try-error") {
                    trt.diff.pos.gp <- NA
                }
                HR.pos.gp <- try(summary(model.temp1.1)$coefficients[trtvar, "exp(coef)"])
                if (class(HR.pos.gp)=="try-error") {
                    HR.pos.gp <- NA
                }
            } else {
                trt.diff.pos.gp <- NA
                HR.pos.gp <- NA
            }

            model.temp1.2 <- try(coxph(model.formula1, data=data, subset=!pred.class))
            if (class(model.temp1.2)[1] != "try-error") {
                trt.diff.neg.gp <- try(summary(model.temp1.2)$coefficients[trtvar, "Pr(>|z|)"])
                if (class(trt.diff.neg.gp)[1] == "try-error") {
                    trt.diff.neg.gp <- NA
                }
                HR.neg.gp <- try(summary(model.temp1.2)$coefficients[trtvar, "exp(coef)"])
                if (class(HR.neg.gp)=="try-error") {
                    HR.neg.gp <- NA
                }
            } else {
                trt.diff.neg.gp <- NA
                HR.neg.gp <- NA
            }
            trt.id <- data[, trtvar]==1

            model.temp2.1 <- try(coxph(model.formula2, data=data, subset=trt.id))
            if (class(model.temp2.1)[1] != "try-error") {
                gp.diff.trt.arm <- try(summary(model.temp2.1)$coefficients["pred.classTRUE", "Pr(>|z|)"])
                if (class(gp.diff.trt.arm)[1] == "try-error") {
                    gp.diff.trt.arm <- NA
                }
                HR.trt <- try(summary(model.temp2.1)$coefficients["pred.classTRUE", "exp(coef)"])
                if (class(HR.trt)=="try-error") {
                    HR.trt <- NA
                }
            } else {
                gp.diff.trt.arm <- NA
                HR.trt <- NA
            }

            model.temp2.2 <- try(coxph(model.formula2, data=data, subset=!trt.id))
            if (class(model.temp2.2)[1] != "try-error") {
                gp.diff.ctrl.arm <- try(summary(model.temp2.2)$coefficients["pred.classTRUE", "Pr(>|z|)"])
                if (class(gp.diff.ctrl.arm)[1] == "try-error") {
                    gp.diff.ctrl.arm <- NA
                }
                HR.ctrl <- try(summary(model.temp2.2)$coefficients["pred.classTRUE", "exp(coef)"])
                if (class(HR.trt)=="try-error") {
                    HR.ctrl <- NA
                }
            } else {
                gp.diff.ctrl.arm <- NA
                HR.ctrl <- NA
            }

            model.temp3 <- try(coxph(model.formula3, data=data))
            if (class(model.temp3)[1] != "try-error") {
                interaction <- try(summary(model.temp3)$coefficients[paste(trtvar, ":pred.classTRUE", sep=""), "Pr(>|z|)"])
                if (class(interaction)[1] == "try-error") {
                    interaction <- NA
                }
            } else {
                interaction <- NA
            }

            data.diag <- rbind(data[data[, trtvar]==1 & data[, "pred.class"]==TRUE, ], data[data[, trtvar]==0 & data[, "pred.class"]==FALSE, ])
            # Union along the diagonal of the 2-by-2 table for treatment and predictive class.
            model.diag <- try(coxph(model.formula1, data=data.diag))
            if (class(model.diag)!="try-error") {
                trt.pos.gp.by.ctrl.neg.gp <- try(summary(model.diag)$coefficients[trtvar, "Pr(>|z|)"])
                if (class(trt.pos.gp.by.ctrl.neg.gp)=="try-error") {
                    trt.pos.gp.by.ctrl.neg.gp <- NA
                }
                HR.trt.pos.ctrl.neg <- try(summary(model.diag)$coefficients[trtvar, "exp(coef)"])
                if (class(HR.trt.pos.ctrl.neg)=="try-error") {
                    HR.trt.pos.ctrl.neg <- NA
                }
            } else {
                trt.pos.gp.by.ctrl.neg.gp <- NA
                HR.trt.pos.ctrl.neg <- NA
            }

            pval <- data.frame(trt.diff.pos.gp=trt.diff.pos.gp,
                               trt.diff.neg.gp=trt.diff.neg.gp,
                               gp.diff.trt.arm=gp.diff.trt.arm,
                               gp.diff.ctrl.arm=gp.diff.ctrl.arm,
                               interaction=interaction,
                               trt.pos.gp.by.ctrl.neg.gp=trt.pos.gp.by.ctrl.neg.gp)

            hazard.ratios <- data.frame(HR.pos.gp=HR.pos.gp, HR.neg.gp=HR.neg.gp, HR.trt=HR.trt,
                                        HR.ctrl=HR.ctrl, HR.trt.pos.ctrl.neg=HR.trt.pos.ctrl.neg)

            results <- list(pval=pval, hazard.ratios=hazard.ratios)

        } else {
            model.text <- paste("Surv(", yvar, ",", censorvar, ")~pred.class", sep="")
            model.formula <- as.formula(model.text)
            model.temp <- try(coxph(model.formula, data=data))
            if (class(model.temp)[1] !="try-error") {
                pval <- try(summary(model.temp)$coefficients["pred.classTRUE", "Pr(>|z|)"])
                if (class(pval)=="try-error") {
                    pval <- NA
                }

                HR <- try(summary(model.temp)$coefficients["pred.classTRUE", "exp(coef)"])
                if (class(HR)=="try-error") {
                    HR <- NA
                }

            } else {
                pval <- NA
                HR <- NA
            }

            pval <- data.frame(pval=pval)
            hazard.ratios <- data.frame(HR=HR)
            results <- list(pval=pval, hazard.ratios=hazard.ratios)
        }
    }

    if (type=="b") {

        if (!is.null(trtvar)) {
            model.text1 <- paste(yvar, "~", trtvar, sep="")
            model.formula1 <- as.formula(model.text1)
            model.text2 <- paste(yvar, "~pred.class", sep="")
            model.formula2 <- as.formula(model.text2)
            model.text3 <- paste(yvar, "~", trtvar ,"*pred.class",sep="")
            model.formula3 <- as.formula(model.text3)

            model.temp1.1 <- try(glm(model.formula1, family=binomial(link="logit"), data=data, subset=pred.class))
            if (class(model.temp1.1)[1] != "try-error") {
                trt.diff.pos.gp <- try(summary(model.temp1.1)$coefficients[trtvar, "Pr(>|z|)"])
                # In some cases, glm may return a result, but without these coefficients, in which case set trt.diff.pos.gp <- NA
                if (class(trt.diff.pos.gp)[1] == "try-error") {
                    trt.diff.pos.gp <- NA
                }
                OR.pos.gp <- try(exp(summary(model.temp1.1)$coefficients[trtvar, "Estimate"]))
                if (class(OR.pos.gp)=="try-error") {
                    OR.pos.gp <- NA
                }
            } else {
                trt.diff.pos.gp <- NA
                OR.pos.gp <- NA
            }

            model.temp1.2 <- try(glm(model.formula1, family=binomial(link="logit"), data=data, subset=!pred.class))
            if (class(model.temp1.2)[1] != "try-error") {
                trt.diff.neg.gp <- try(summary(model.temp1.2)$coefficients[trtvar, "Pr(>|z|)"])
                if (class(trt.diff.neg.gp)[1] == "try-error") {
                    trt.diff.neg.gp <- NA
                }
                OR.neg.gp <- try(exp(summary(model.temp1.2)$coefficients[trtvar, "Estimate"]))
                if (class(OR.neg.gp)=="try-error") {
                    OR.neg.gp <- NA
                }
            } else {
                trt.diff.neg.gp <- NA
                OR.neg.gp <- NA
            }

            trt.id <- data[, trtvar]==1

            model.temp2.1 <- try(glm(model.formula2, family=binomial(link="logit"), data=data, subset=trt.id))
            if (class(model.temp2.1)[1] != "try-error") {
                gp.diff.trt.arm <- try(summary(model.temp2.1)$coefficients["pred.classTRUE", "Pr(>|z|)"])
                if (class(gp.diff.trt.arm)[1] == "try-error") {
                    gp.diff.trt.arm <- NA
                }
                OR.trt <- try(exp(summary(model.temp2.1)$coefficients["pred.classTRUE", "Estimate"]))
                if (class(OR.trt)=="try-error") {
                    OR.trt <- NA
                }
            } else {
                gp.diff.trt.arm <- NA
                OR.trt <- NA
            }

            model.temp2.2 <- glm(model.formula2, family=binomial(link="logit"), data=data, subset=!trt.id)
            if (class(model.temp2.2)[1] != "try-error") {
                gp.diff.ctrl.arm <- try(summary(model.temp2.2)$coefficients["pred.classTRUE", "Pr(>|z|)"])
                if (class(gp.diff.ctrl.arm)[1] == "try-error") {
                    gp.diff.ctrl.arm <- NA
                }
                OR.ctrl <- try(exp(summary(model.temp2.2)$coefficients["pred.classTRUE", "Estimate"]))
                if (class(OR.trt)=="try-error") {
                    OR.ctrl <- NA
                }
            } else {
                gp.diff.ctrl.arm <- NA
                OR.ctrl <- NA
            }

            model.temp3 <- try(glm(model.formula3, family=binomial(link="logit"), data=data))
            if (class(model.temp3)[1] != "try-error") {
                interaction <- try(summary(model.temp3)$coefficients[paste(trtvar, ":pred.classTRUE", sep=""), "Pr(>|z|)"])
                if (class(interaction)[1] == "try-error") {
                    interaction <- NA
                }
            } else {
                interaction <- NA
            }

            data.diag <- rbind(data[data[, trtvar]==1 & data[, "pred.class"]==TRUE, ], data[data[, trtvar]==0 & data[, "pred.class"]==FALSE, ])
            # Union along the diagonal of the 2-by-2 table for treatment and predictive class.
            model.diag <- try(glm(model.formula1, family=binomial(link="logit"), data=data.diag))
            if (class(model.diag)[1] !="try-error") {
                trt.pos.gp.by.ctrl.neg.gp <- try(summary(model.diag)$coefficients[trtvar, "Pr(>|z|)"])
                if (class(trt.pos.gp.by.ctrl.neg.gp)=="try-error") {
                    trt.pos.gp.by.ctrl.neg.gp <- NA
                }
                OR.trt.pos.ctrl.neg <- try(exp(summary(model.diag)$coefficients[trtvar, "Estimate"]))
                if (class(OR.trt.pos.ctrl.neg)=="try-error") {
                    OR.trt.pos.ctrl.neg <- NA
                }
            } else {
                trt.pos.gp.by.ctrl.neg.gp <- NA
                OR.trt.pos.ctrl.neg <- NA
            }

            pval <- data.frame(trt.diff.pos.gp=trt.diff.pos.gp,
                               trt.diff.neg.gp=trt.diff.neg.gp,
                               gp.diff.trt.arm=gp.diff.trt.arm,
                               gp.diff.ctrl.arm=gp.diff.ctrl.arm,
                               interaction=interaction,
                               trt.pos.gp.by.ctrl.neg.gp=trt.pos.gp.by.ctrl.neg.gp)

            odds.ratios <- data.frame(OR.pos.gp=OR.pos.gp, OR.neg.gp=OR.neg.gp, OR.trt=OR.trt,
                                      OR.ctrl=OR.ctrl, OR.trt.pos.ctrl.neg=OR.trt.pos.ctrl.neg)

            results <- list(pval=pval, odds.ratios=odds.ratios)
        } else {
            model.text <- paste(yvar, "~pred.class", sep="")
            model.formula <- as.formula(model.text)
            model.temp <- try(glm(model.formula, family=binomial(link="logit"), data=data))

            if (class(model.temp)[1] !="try-error") {
                pval <- try(summary(model.temp)$coefficients["pred.classTRUE", "Pr(>|z|)"])

                if (class(pval)=="try-error") {
                    pval <- NA
                }

                OR <- try(exp(summary(model.temp)$coefficients["pred.classTRUE", "Estimate"]))

                if (class(OR)=="try-error") {
                    OR <- NA
                }
            } else {
                pval <- NA
                OR <- NA
            }

            bin.stats = binary.stats(data[,"pred.class"], data[,yvar])
            pval <- data.frame(pval=pval)
            odds.ratios <- data.frame(OR=OR)
            results <- list(pval=pval, odds.ratios=odds.ratios,bin.stats=bin.stats)
        }
    }

    if (type=="c") {
        if (!is.null(trtvar)) {
            model.text1 <- paste(yvar, "~", trtvar, sep="")
            model.formula1 <- as.formula(model.text1)
            model.text2 <- paste(yvar, "~pred.class", sep="")
            model.formula2 <- as.formula(model.text2)
            model.text3 <- paste(yvar, "~",trtvar ,"*pred.class",sep="")
            model.formula3 <- as.formula(model.text3)

            model.temp1.1 <- try(lm(model.formula1, data=data, subset=pred.class))
            if (class(model.temp1.1)[1] != "try-error") {
                trt.diff.pos.gp <- try(summary(model.temp1.1)$coefficients[trtvar, "Pr(>|t|)"])
                if (class(trt.diff.pos.gp)[1] == "try-error") {
                    trt.diff.pos.gp <- NA
                }
            } else {
                trt.diff.pos.gp <- NA
            }

            model.temp1.2 <- try(lm(model.formula1, data=data, subset=!pred.class))
            if (class(model.temp1.2)[1] != "try-error") {
                trt.diff.neg.gp <- try(summary(model.temp1.2)$coefficients[trtvar, "Pr(>|t|)"])
                if (class(trt.diff.neg.gp)[1] == "try-error") {
                    trt.diff.neg.gp <- NA
                }
            } else {
                trt.diff.neg.gp <- NA
            }

            trt.id <- data[, trtvar]==1

            model.temp2.1 <- try(lm(model.formula2, data=data, subset=trt.id))
            if (class(model.temp2.1)[1] != "try-error") {
                gp.diff.trt.arm <- try(summary(model.temp2.1)$coefficients["pred.classTRUE", "Pr(>|t|)"])
                if (class(gp.diff.trt.arm)[1] == "try-error") {
                    gp.diff.trt.arm <- NA
                }
            } else {
                gp.diff.trt.arm <- NA
            }

            model.temp2.2 <- try(lm(model.formula2, data=data, subset=!trt.id))
            if (class(model.temp2.2)[1] != "try-error") {
                gp.diff.ctrl.arm <- try(summary(model.temp2.2)$coefficients["pred.classTRUE", "Pr(>|t|)"])
                if (class(gp.diff.ctrl.arm)[1] == "try-error") {
                    gp.diff.ctrl.arm <- NA
                }
            } else {
                gp.diff.ctrl.arm <- NA
            }

            model.temp3 <- try(lm(model.formula3, data=data))
            if (class(model.temp3)[1] != "try-error") {
                interaction <- try(summary(model.temp3)$coefficients[paste(trtvar, ":pred.classTRUE", sep=""), "Pr(>|t|)"])
                if (class(interaction)[1] == "try-error") {
                    interaction <- NA
                }
            } else {
                interaction <- NA
            }

            data.diag <- rbind(data[data[, trtvar]==1 & data[, "pred.class"]==TRUE, ], data[data[, trtvar]==0 & data[, "pred.class"]==FALSE, ])
            # Union along the diagonal of the 2-by-2 table for treatment and predictive class.
            model.diag <- try(lm(model.formula1, data=data.diag))
            if (class(model.diag)!="try-error") {
                trt.pos.gp.by.ctrl.neg.gp <- try(summary(model.diag)$coefficients[trtvar, "Pr(>|t|)"])
                if (class(trt.pos.gp.by.ctrl.neg.gp)=="try-error") {
                    trt.pos.gp.by.ctrl.neg.gp <- NA
                }
            } else {
                trt.pos.gp.by.ctrl.neg.gp <- NA
            }

            pval <- data.frame(trt.diff.pos.gp=trt.diff.pos.gp,
                               trt.diff.neg.gp=trt.diff.neg.gp,
                               gp.diff.trt.arm=gp.diff.trt.arm,
                               gp.diff.ctrl.arm=gp.diff.ctrl.arm,
                               interaction=interaction,
                               trt.pos.gp.by.ctrl.neg.gp=trt.pos.gp.by.ctrl.neg.gp)

            results <- list(pval=pval)
        } else {
            model.text <- paste(yvar, "~pred.class", sep="")
            model.formula <- as.formula(model.text)
            model.temp <- try(lm(model.formula, data=data))
            if (class(model.temp)[1] !="try-error") {
                pval <- try(summary(model.temp)$coefficients["pred.classTRUE", "Pr(>|t|)"])
                if (class(pval)=="try-error") {
                    pval <- NA
                }
            } else {
                pval <- NA
            }
            pval <- data.frame(pval=pval)
            results <- list(pval=pval)
        }
    }
    results
}


#' A function for interaction plot
#'
#' @param data.eval output of evaluate.results or summarize.cv.stats
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c".
#' @param main title of the plot
#' @param trt.lab treatment label
#'
#' @return A ggplot object.
#'
interaction.plot <- function(data.eval, type, main="Interaction Plot", trt.lab=c("Trt.", "Ctrl.")) {
  group<-resp<-trt<-se<-NULL

  #data.eval: output of evaluate.results or summarize.cv.stats

  group.stats=data.eval$group.stats

  if (dim(group.stats)[1]<4){
    print("No interaction plot for prognostic case.")
    return(NULL)
  }else{
#    library(ggplot2)

    group.stats=group.stats[c("sig+.trt", "sig+.ctrl", "sig-.trt", "sig-.ctrl"),]
    group.stats[,"trt"]=factor(trt.lab[c(1,2,1,2)],levels=trt.lab)
    group.stats[,"group"]=factor(c("sig+","sig+","sig-","sig-"),levels=c("sig-","sig+"))

    if (type=="s"){
      group.stats[,"resp"]=group.stats[,"mean.surv.time"]
      group.stats[,"se"]=group.stats[,"se.mean.surv"]
      ylabel="Restricted Mean Survival Time"
    }

    if (type=="c"){
      group.stats[,"resp"]=group.stats[,"mean"]
      group.stats[,"se"]=group.stats[,"sd"]/(group.stats[,"n"])^0.5
      ylabel="Mean"
    }

    if (type=="b"){
      group.stats[,"resp"]=group.stats[,"resp.rate"]
      group.stats[,"se"]=(group.stats[,"resp.rate"]*(1-group.stats[,"resp.rate"])/group.stats[,"n"])^0.5
      ylabel="Response/Event Rate"
    }


    fig=ggplot(group.stats, aes(x=group, y=resp, group=trt, color=trt, linetype=trt)) +
      geom_errorbar(aes(ymin=resp-se, ymax=resp+se,color=trt), width=0,size=0.8,position=position_dodge(width=0.2)) +
      xlab("") +
      ylab(ylabel) +
      #            labs(color="") +
      geom_line(size=0.8, position=position_dodge(width=0.2)) +
      geom_point(position=position_dodge(width=0.2),size=0.8,aes(color=trt))+
      theme_bw(base_size=15) +
      ggtitle(main) +
      #scale_y_continuous(limits=c(1000,4000),breaks=c(1000,2000,3000,4000),labels=c("1","2","3","4"))+
      theme(text=element_text(size=15), axis.text=element_text(size=15),legend.text=element_text(size=15),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title = element_text(hjust = 0.5))
    fig
  }

}


#' A function for binary statistics
#'
#' @param pred.class predicted output for each subject
#' @param y.vec response vector
#'
#' @return a data frame with sensitivity, specificity, NPV, PPV and accuracy
#'
binary.stats <- function(pred.class, y.vec) {
    table <- table(pred.class, y.vec)
    N <- length(pred.class)
    # These should be OK even if table has only 1 row.

    if (nrow(table)==2 & ncol(table)==2) {
        # pred.class contains 0s and 1s, so table is 2-by-2.
        TP <- table[2,2]
        FP <- table[2,1]
        TN <- table[1,1]
        FN <- table[1,2]
        PPV <- TP / (TP + FP)
        NPV <- TN / (TN + FN)

    } else {
        # pred.class has only 1 distict value or events has only 1 distinct value
        if (length(unique(pred.class)) == 1&&unique(pred.class) == 1) {
            # pred.class is all 1s
            FP <- table[1,1]
            TP <- table[1,2]
            TN <- 0
            FN <- 0
            PPV <- TP / (TP + FP)
            NPV <- 0
            # NPV undefined

        }
        if (length(unique(pred.class)) == 1&&unique(pred.class) == 0) {
            # pred.class is all 0s
            FP <- 0
            TP <- 0
            TN <- table[1,1]
            FN <- table[1,2]
            PPV <- 0
            NPV <- TN / (TN + FN)
            # NPV undefined

        }
    }
    sens <- TP / (TP + FN)
    spec <- TN / (TN + FP)
    acc <- (TP + TN) / N
    #mcc <- (TP * TN - FP * FN) / sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN)

    diag.stats <- data.frame("sens"=sens, "spec"=spec,
                             "PPV"=PPV, "NPV"=NPV, "acc"=acc)
    return(diag.stats)
}

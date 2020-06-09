#
# aim_batting_v11
# Feb 2017
#
#
# The functions in this file implement versions of the Adaptive Index Model (AIM) method
# for predictive modeling.
#
# Note: To use these functions, you must first do the following:
#
# 1. Load the AIM library - library(AIM)
# 2. source the following files:
#
#   * aim_batting.R (this file) - Note: aim_batting overrides some of the functions in the AIM library, so
#                                 do *not* load the AIM library after sourcing aim.batting.
#   * filter.R
#   * cross_validation_functions.R
#
# This file contains the following functions:
#
#   * aim.rule.batting - Implements a modification of the original method, which does the following:
#
#     - Performs BATTing on AIM rules, adding one rule at a time, starting from the most significant rule.
#     -	Finds the best number of rules to separate the positive and negative groups.
#
#   * cv.aim.batting - Performs k-fold cross-validation on aim.batting.
#
#   * cv.aim.rule.batting - Performs k-fold cross-validation on aim.rule.batting.
#
#   * pred.aim - Takes an AIM rule, returned by aim.batting or aim.rule.batting,
#                  and the predictors, and returns a TRUE-FALSE vector whose length is
#                  the number of rows of the original data, giving the predictive class for each row.
#

###################################################################################################
#
#  Functions Shared by aim.batting and aim.rule.batting
#
###################################################################################################

#'  Find score of cutoff (returned by aim.find.cutoff.pred) for predictive case.
#'
#' @param data data frame containing the response, covariate, treatment variable and censoring variable (only for time to event response).
#' @param yvar response variable name.
#' @param censorvar censoring variable name 1:event; 0: censor.
#' @param trtvar treatment variable name.
#' @param trtref code for treatment arm.
#' @param xvar covariate variable name.
#' @param type "c" continuous; "s" survival; "b" binary.
#' @param cutoff cutpoint of interest.
#' @param nsubj number of subjects.
#' @param min.sigp.prcnt desired proportion of signature positive group size.
#'
#' @return AIM score for a single covariate-cutoff combination.
#' 
aim.score.pred <- function(data, yvar, censorvar, trtvar, trtref, xvar, type, cutoff, nsubj,min.sigp.prcnt) {


    id <- 1*(data[, xvar]>cutoff)
    n.grp<-table(id) #*#
    sigp.prcnt<-sum(id)/nsubj #*#

#     n.trt <- table(data[id, trtvar])
#     n.id.trt <- table(id, data[, trtvar])
    id.trt <- 1*(data[, trtvar]==trtref)
#     idbytrt=id*id.trt

    if (type=="s") {
        model.text <- paste("Surv(", yvar, ",", censorvar, ")~","id*id.trt", sep="")
        model.formula <- as.formula(model.text)

        if (sigp.prcnt>=min.sigp.prcnt & length(n.grp)>=2) {
            cph.fit <- try(coxph(model.formula, data=data), silent=TRUE)

            if ((class(cph.fit)!="try-error") & (is.na(summary(cph.fit)$coefficients[3,5])==FALSE)) {
                hr.results <- summary(cph.fit)$coefficients[3,c(4,5)]
                score <- c(hr.results[2], hr.results[2])
            } else {
                score <- c(1,0)
            }
        } else {
            score <- c(1,0)
        }
    }

    if (type=="b") {
        model.text <- paste(yvar, "~", "id*id.trt", sep="")
        model.formula <- as.formula(model.text)

        if (sigp.prcnt>=min.sigp.prcnt & length(n.grp)>=2) {
            #
            # See Criteria 1-3 in SIDES paper
            #
            glm.fit <- try(glm(model.formula, family=binomial,data=data), silent=TRUE)
            or.results <- try(summary(glm.fit)$coefficients[4, c(3,4)], silent=TRUE)
            if (class(glm.fit)!="try-error" & class(or.results)!="try-error") {
                score <- c(or.results[2], or.results[2])
            } else {
                score <- c(1,0)
            }
        } else {
            score <- c(1,0)
        }
    }

    if (type=="c") {
        model.text <- paste(yvar, "~", "id*id.trt", sep="")
        model.formula <- as.formula(model.text)

        if (sigp.prcnt>=min.sigp.prcnt & length(n.grp)>=2) {
            lm.fit <- try(lm(model.formula, data=data), silent=TRUE)
            beta.results <- try(summary(lm.fit)$coefficients[4, c(3,4)],silent=TRUE)
            if (class(lm.fit)!="try-error" & class(beta.results)!= "try-error") {
                score <- c(beta.results[2], beta.results[2])
            } else {
                score <- c(1,0)
            }
        } else {
            score <- c(1,0)
        }
    }
    score
}

#' Find score of cutoff (returned by aim.find.cutoff.pred) for prognostic case.
#'
#' @param data data frame containing the response, covariate, treatment variable and censoring variable (only for time to event response).
#' @param yvar response variable name.
#' @param censorvar censoring variable name 1:event; 0: censor.
#' @param xvar covariate variable name.
#' @param type "c" continuous; "s" survival; "b" binary.
#' @param cutoff cutpoint of interest.
#' @param nsubj number of subjects.
#' @param min.sigp.prcnt desired proportion of signature positive group size.
#'
#' @return AIM score for a single covariate-cutoff combination.
#'
aim.score.prog <- function(data, yvar, censorvar, xvar, type, cutoff, nsubj,min.sigp.prcnt) {

    id <- 1*(data[, xvar]>cutoff)
    n.grp<-table(id) #*#
    sigp.prcnt<-sum(id)/nsubj #*#

    if (type=="s") {
        model.text <- paste("Surv(", yvar, ",", censorvar, ")~", "id", sep="")
        model.formula <- as.formula(model.text)

        if (sigp.prcnt>=min.sigp.prcnt & length(n.grp)>=2) { # Xin:  min(n.id.trt)>20 & length(n.id.trt)>1
            cph1.fit <- try(coxph(model.formula, data=data), silent=TRUE)

            if ((class(cph1.fit)!="try-error") & (is.na(summary(cph1.fit)$coefficients[5])==FALSE)) {
                hr1.results <- summary(cph1.fit)$coefficients[c(4,5)]
                score <- hr1.results[2]
            } else {
                score <- 1
            }
        } else {
            score <- 1
        }

    }

    if (type=="b") {
        model.text <- paste(yvar, "~", "id", sep="")
        model.formula <- as.formula(model.text)

        if (sigp.prcnt>=min.sigp.prcnt & length(n.grp)>=2) { # Xin:  min(n.id.trt)>20 & length(n.id.trt)>1
            glm1.fit <- try(glm(model.formula, family=binomial,data=data), silent=TRUE)

            if (class(glm1.fit)!="try-error") {
                or1.results <- summary(glm1.fit)$coefficients[2, c(3,4)]
                score <- or1.results[2]
            } else {
                score <- 1
            }
        } else {
            score <- 1
        }
    }

    if (type=="c") {
        model.text <- paste(yvar, "~", "id", sep="")
        model.formula <- as.formula(model.text)

        if (sigp.prcnt>=min.sigp.prcnt & length(n.grp)>=2) { # Xin:  min(n.id.trt)>20 & length(n.id.trt)>1
            lm1.fit <- try(lm(model.formula, data=data), silent=TRUE)
            if (class(lm1.fit)!="try-error") {
                beta1.results <- summary(lm1.fit)$coefficients[2, c(3,4)]
                score <- beta1.results[2]
            } else {
                score <- 1
            }
        } else {
            score <- 1
        }
    }
    score
}

###################################################################################################
#
#  AIM BATTing
#
###################################################################################################


#' The main AIM-BATTing function 
#' @description  This function finds the aim score for each subject in the dataset and using aim score as the predictor, performs BATTing to find the best threshold for each predictor.
#'
#' @param y data frame of the response variable.
#' @param x data frame of predictors, each column of which corresponds to a variable.
#' @param censor.vec data frame indicating censoring for survival data. For binary or continuous data, set censor.vec <- NULL.
#' @param trt.vec data frame indicating whether or not the patient was treated. For the pronostic case, set trt.vec <- NULL.
#' @param trtref code for treatment arm.
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c".
#' @param n.boot number of bootstraps in bootstrapping step.
#' @param des.res the desired response. "larger": prefer larger response; "smaller": prefer smaller response.
#' @param min.sigp.prcnt desired proportion of signature positive group size.
#' @param mc.iter # of iterations for the MC procedure to get a stable "best number of predictors".
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected.
#' @param filter.method NULL, no prefiltering, "univariate", univariate filtering; "glmnet", glmnet filtering; "unicart", CART filtering (only for prognostic case).
#'
#' @importFrom AIM cox.main
#'
#' @return A list containing variables in signature and their thresholds.
#'
aim.batting <- function(y, x, censor.vec=NULL, trt.vec=NULL, trtref=NULL, type, n.boot, des.res="larger", min.sigp.prcnt=0.20, mc.iter=1, mincut=0.1, pre.filter=NULL, filter.method=NULL) {


    BATTing <- function(yvar, censorvar=NULL, trtvar=NULL, trtref=NULL, xvar, data, type="c", n.boot=50, min.sigp.prcnt=0.20) {
        # Returns cutoff for a single x variable.
        niter <- n.boot
        nsubj <- nrow(data)
        cutoff.vec <- NULL

        for (i in 1:niter) {
            train.id <- sample(1:nsubj, nsubj, replace=TRUE)
            data.train <- data[train.id, ]
            if (!is.null(trtvar)) {
                cutoff.temp <- aim.find.cutoff.pred(data=data.train, yvar=yvar, censorvar=censorvar, trtvar=trtvar, xvar=xvar, type=type, trtref=trtref,nsubj,min.sigp.prcnt)
            } else {
                cutoff.temp <- aim.find.cutoff.prog(data=data.train, yvar=yvar, censorvar=censorvar, xvar=xvar, type=type, nsubj,min.sigp.prcnt)
            }
            cutoff.vec <- c(cutoff.vec, cutoff.temp)
        }

        if (all(is.na(cutoff.vec))) {
            cutoff.med <- NA
        } else {
            cutoff.vec <- cutoff.vec[!is.na(cutoff.vec)]
            cutoff.unique <- sort(unique(cutoff.vec))
            cutoff.tab <- table(cutoff.vec)
            cutoff.med <- cutoff.unique[which.max(cutoff.tab)]
        }

        cutoff.med
    }

    ###########################################################################
    #
    #                 Predictive BATTing Functions
    #
    ###########################################################################


    aim.find.cutoff.pred <- function(data, yvar, censorvar, trtvar, xvar, type, trtref, nsubj,min.sigp.prcnt) {
        cut.vec <- sort(unique(data[, xvar]))
        cut.score <- lapply(cut.vec, function(x) aim.score.pred(data=data, yvar=yvar, censorvar=censorvar, xvar=xvar, trtvar=trtvar, trtref=trtref, type=type, cutoff=x, nsubj=nsubj, min.sigp.prcnt=min.sigp.prcnt))
        cut.score.mat <- do.call(rbind, cut.score)
        cut.score.mat <- cbind(cut.vec, cut.score.mat)
        #signeg.status <- cut.score.mat[,3]>0.05  Should this be 0?
        signeg.status <- (cut.score.mat[,2]<1)|(cut.score.mat[,3] > 0)

        if (any(signeg.status)) { # If any of the signeg.status is TRUE
            cut.score.mat2 <- cut.score.mat[signeg.status, ]

            if (sum(signeg.status)<=1) {
                cut.val <- cut.score.mat2[1]
            } else {
                sigpos.index <- which.min(cut.score.mat2[, 2])
                cut.val <- cut.score.mat2[sigpos.index, 1]
            }
        } else {
            cut.val <- NA
        }
        cut.val
    }

    ###########################################################################
    #
    #                 Prognostic BATTing Functions
    #
    ###########################################################################

    aim.find.cutoff.prog <- function(data, yvar, censorvar, xvar, type, nsubj, min.sigp.prcnt) {
        # Find optimal cutoff for prognostic case.
        cut.vec <- sort(unique(data[, xvar]))
        cut.score <- lapply(cut.vec, function(x) aim.score.prog(data=data, yvar=yvar, censorvar=censorvar, xvar=xvar, type=type, cutoff=x, nsubj, min.sigp.prcnt))
        cut.score.mat <- do.call(c, cut.score)
        cut.score.mat <- cbind(cut.vec, cut.score.mat)
        sigpos.index <- which.min(cut.score.mat[, 2])
        cut.val <- cut.score.mat[sigpos.index, 1]
        cut.val
    }

    aim.convert <- function(model, nvar, marker.names) {
        aim.model <- model$res[[nvar]]
        var.names <- marker.names[aim.model[, "jmax"]]
        dir <- rep(">", nrow(aim.model))
        dir[aim.model[, "maxdir"]==-1] <- "<"
        aim.out <- data.frame(variable=var.names, direction=dir, cutoff=aim.model[, "cutp"])
        aim.out
    }

    ################### Main aim.batting Code ##################################################

    # Check whether the number of columns of x is 1.
    if (ncol(x)==1) {
        stop("The number of columns of x must be greater than 1.")
    }

    data <- cbind(y, x)
    xvars <- names(x)
    yvar <- names(y)
    marker.data <- as.matrix(x)
    response <- y[, yvar]

    if (!is.null(trt.vec)) {
        trtvar <- names(trt.vec)
        Trt <- trt.vec[, trtvar]

        if (!is.null(trtref)) {
            # If trtref is provided, create 0-1 vector for trt1 in which 1 corresponds to trtref.
            if (is.element(trtref, unique(Trt))) {
                trt1 <- rep(0,length=length(Trt))
                trt1[Trt==trtref] <- 1
            } else {
                stop("trtref must be one of the entries of trt.vec.")
            }
        } else {
            # trtref not provided, so assume Trt is 0-1 vector.
            if (setequal(union(c(0,1), unique(Trt)), c(0,1))) {
                # Elements of Trt are subset of c(0,1).
                trt1 <- Trt
            } else {
                stop("If trt.vec is not a 0-1 vector, trtref must be supplied.")
            }
        }
        trt.vec[, trtvar] <- trt1
        data <- cbind(data, trt.vec)
        trtref <- 1
    } else {
        trtvar <- NULL
    }

    if (!is.null(censor.vec)) {
        data <- cbind(data, censor.vec)
        censorvar <- names(censor.vec)
        censor <- censor.vec[, censorvar]
    } else {
        censorvar <- NULL
        censor <- NULL
    }


    if (!is.null(pre.filter)){
      xvars_new=filter(data=data,type=type,yvar=yvar,xvars=xvars,censorvar=censorvar,trtvar=trtvar,n.boot=50,cv.iter=15,pre.filter=pre.filter, filter.method=filter.method)
      del_ind=which(!(xvars %in% xvars_new))+1
      data=data[,-1*del_ind]
      x=x[,xvars_new]
      xvars_copy=xvars
      xvars=xvars_new
      marker.data=as.matrix(x)
    }

    nvars <- length(xvars)
    nvar=rep(0,mc.iter)
    i.mc=1
    err.mc=0
    min.ianderr=10
    max.properr=0.9 #y

    if (type=="s" & !is.null(trt.vec)) {

        while (i.mc <= mc.iter) {

            aim.cv.model <- try(cv.cox.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE, status=censor, nsteps=nvars - 1, mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc <- i.mc+1
            } else{
                err.mc <- err.mc+1
            }

            if ((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr){
                err.mc="error"
                break
            }
        }

        if (err.mc != "error"){
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            nvar <- nvar.unique[which.max(nvar.tab)]
            aim.model <- cox.interaction(x=marker.data, trt=trt1, delta=censor, y=response, backfit=TRUE, nsteps=nvar,mincut=mincut)
        }
    }

    if (type=="s" & is.null(trt.vec)) {

        while (i.mc <=mc.iter) {
            aim.cv.model <- try(cv.cox.main(x=marker.data, status=censor, y=response,backfit=TRUE, nsteps=nvars - 1,mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc=i.mc+1
            } else {
                err.mc=err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr) {
                err.mc="error"
                break
            }
        }

        if (err.mc != "error"){
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            nvar <- nvar.unique[which.max(nvar.tab)]
            aim.model <- cox.main(x=marker.data, delta=censor, y=response, backfit=TRUE, nsteps=nvar,mincut=mincut)
        }
    }

    if (type=="b" & !is.null(trt.vec)) {

        while (i.mc <=mc.iter) {
            aim.cv.model <- try(cv.logistic.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE, nsteps=nvars - 1,mincut=mincut), silent=TRUE)
            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc=i.mc+1
            } else {
                err.mc=err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr) {
                err.mc="error"
                break
            }
        }

        if (err.mc != "error"){
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            nvar <- nvar.unique[which.max(nvar.tab)]
            aim.model <- logistic.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE, nsteps=nvar,mincut=mincut)
        }
    }

    if (type=="b" & is.null(trt.vec)) {
        while (i.mc <=mc.iter){
            aim.cv.model <- try(cv.logistic.main(x=marker.data, y=response, backfit=TRUE, nsteps=nvars - 1,mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc=i.mc+1
            } else {
                err.mc=err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr){
                err.mc="error"
                break
            }
        }

        if (err.mc != "error") {
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            nvar <- nvar.unique[which.max(nvar.tab)]
            aim.model <- logistic.main(x=marker.data, y=response, backfit=TRUE, nsteps=nvar,mincut=mincut)
        }
    }

    if (type=="c" & !is.null(trt.vec)) {

        while (i.mc <=mc.iter) {
            aim.cv.model <- try(cv.lm.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE,nsteps=nvars - 1, mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc=i.mc+1
            } else {
                err.mc=err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr) {
                err.mc="error"
                break
            }
        }

        if (err.mc != "error"){
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            nvar <- nvar.unique[which.max(nvar.tab)]
            aim.model <- lm.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE, nsteps=nvar,mincut=mincut)
        }
    }

    if (type=="c" & is.null(trt.vec)) {

        while (i.mc <=mc.iter) {
            aim.cv.model <- try(cv.lm.main(x=marker.data, y=response, backfit=TRUE, nsteps=nvars - 1,mincut=mincut), silent=TRUE)
            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc=i.mc+1
            } else {
                err.mc=err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr) {
                err.mc="error"
                break
            }
        }

        if (err.mc != "error") {
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            nvar <- nvar.unique[which.max(nvar.tab)]
            aim.model <- lm.main(x=marker.data, y=response,backfit=TRUE, nsteps=nvar,mincut=mincut)
        }
    }

    if (err.mc != "error") {
        pred.aim <- index.prediction(aim.model$res[[nvar]], marker.data)

        if (!is.null(trt.vec)) {
            dir.temp <- summary(lm(response~trt1*pred.aim))$coefficient[4,1]
        } else {
            dir.temp=cor(response,pred.aim)
        }

        if ((des.res=="smaller" & dir.temp>0) | (des.res=="larger" & dir.temp<0)) {
            for (nvar.temp in c(1:nvar)){
                model.temp <- aim.model$res[[nvar.temp]]
                xvars.temp <- xvars[model.temp[,"jmax"]]
                xvars.temp.value <- lapply(as.data.frame(cbind(NULL,marker.data[,xvars.temp])), function(x) sort(unique(x)))

                for (i.xvars in 1:length(xvars.temp)){
                    pos.temp <- which(xvars.temp.value[[i.xvars]]==model.temp[i.xvars,"cutp"])
                    model.temp[i.xvars,"cutp"]=xvars.temp.value[[i.xvars]][pos.temp+model.temp[i.xvars,"maxdir"]]
                }

                model.temp[,"maxdir"] <- -1*model.temp[,"maxdir"]
                aim.model$res[[nvar.temp]]=model.temp
            }

            pred.aim <- index.prediction(aim.model$res[[nvar]], marker.data)
        }


        if (!is.null(trt.vec)) {
            aim.data <- data.frame(response=response, censor=data[, censorvar], Trt=trt.vec[, trtvar], index=pred.aim)
            batting.results <- BATTing(yvar="response", censorvar="censor", trtvar="Trt", trtref=trtref, xvar="index", data=aim.data, type=type, n.boot=n.boot, min.sigp.prcnt=min.sigp.prcnt)
        } else {
            aim.data <- data.frame(response=response, censor=data[, censorvar], index=pred.aim)
            batting.results <- BATTing(yvar="response", censorvar="censor", xvar="index", data=aim.data, type=type, n.boot=n.boot, min.sigp.prcnt=min.sigp.prcnt)
        }
        aim.out <- aim.convert(aim.model, nvar, xvars)
        aim.model=aim.model$res[[nvar]]

        if(!is.null(pre.filter)){
          xvars.pos=match(xvars[aim.model[, "jmax"]],xvars_copy)
          aim.model[,"jmax"]=xvars.pos
        }

        output <- list(aim.model=aim.model, aim.rule=aim.out, bat.cutoff=batting.results, Note="Sig+ Grp: score > bat.cutoff")
    } else {
        output <- NULL
        print("Error in call to the AIM library.")
    }
    output
}

###################################################################################################


###################################################################################################
#
#  AIM-Rule BATTing
#
###################################################################################################

#' The main AIM-Rule-BATTing function 
#' @description This function first uses AIM to get the candidate rules and then applies Sequential BATTing to get the best rule(s).
#'
#' @param y data frame of the response variable.
#' @param x data frame of predictors, each column of which corresponds to a variable.
#' @param censor.vec data frame indicating censoring for survival data. For binary or continuous data, set censor.vec <- NULL.
#' @param trt.vec data frame indicating whether or not the patient was treated. For the pronostic case, set trt.vec <- NULL.
#' @param trtref code for treatment arm.
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c".
#' @param n.boot number of bootstraps in bootstrapping step.
#' @param des.res the desired response. "larger": prefer larger response; "smaller": prefer smaller response.
#' @param min.sigp.prcnt desired proportion of signature positive group size.
#' @param mc.iter # of iterations for the MC procedure to get a stable "best number of predictors".
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected.
#' @param filter.method NULL, no prefiltering, "univariate", univariate filtering; "glmnet", glmnet filtering; "unicart", CART filtering (only for prognostic case).
#'
#' @return A list of containing variables in signature and their thresholds.
#'
aim.rule.batting <- function(y, x, censor.vec=NULL, trt.vec=NULL, trtref=NULL, type, n.boot, des.res="larger", min.sigp.prcnt=0.2, mc.iter=1, mincut=0.1, pre.filter=NULL, filter.method=NULL) {

    BATTing <- function(yvar, censorvar=NULL, trtvar=NULL, trtref=NULL, data, aim.model, marker.data, type="c", n.boot=50, min.sigp.prcnt=0.20) {
        # Returns cutoff for a single x variable.
        niter <- n.boot
        nsubj <- nrow(data)
        cutoff.vec <- NULL
        for (i in 1:niter) {
            train.id <- sample(1:nsubj, nsubj, replace=TRUE)
            data.train <- data[train.id, ]
            marker.train <- marker.data[train.id,]
            if (!is.null(trtvar)) {
                cutoff.temp <- aim.find.cutoff.pred(data=data.train, aim.model=aim.model, marker=marker.train, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, trtref=trtref, nsubj, min.sigp.prcnt)
            } else {
                cutoff.temp <- aim.find.cutoff.prog(data=data.train, aim.model=aim.model, marker=marker.train, yvar=yvar, censorvar=censorvar, type=type, nsubj, min.sigp.prcnt)
            }
            cutoff.vec <- c(cutoff.vec, cutoff.temp)
        }

        if (all(is.na(cutoff.vec))) {
            cutoff.med <- NA
        } else {
            cutoff.vec <- cutoff.vec[!is.na(cutoff.vec)]
            cutoff.unique <- sort(unique(cutoff.vec))
            cutoff.tab <- table(cutoff.vec)
            cutoff.med <- cutoff.unique[which.max(cutoff.tab)]
        }

        cutoff.med
    }

    ###########################################################################
    #
    #                 Predictive BATTing Functions
    #
    ###########################################################################

    aim.find.cutoff.pred <- function(data, aim.model, marker, yvar, censorvar, trtvar, type, trtref, nsubj, min.sigp.prcnt) {
        # Find optimal cutoff for predictive case.

        nmax <- length(aim.model$res)
        cut.score.mat <- NULL

        for (nvar.i in 1:nmax){
            index <- index.prediction(aim.model$res[[nvar.i]], marker)
            data$index <- index
            cut.score <- aim.score.pred(data=data, yvar=yvar, censorvar=censorvar, xvar="index", trtvar=trtvar, trtref=trtref, type=type, cutoff=nvar.i-1, nsubj=nsubj, min.sigp.prcnt=min.sigp.prcnt)
            cut.score.mat <- rbind(cut.score.mat,c(nvar.i,cut.score))
        }

        #signeg.status <- cut.score.mat[,3]>0.05  Should this be 0?
        signeg.status <- (cut.score.mat[,2]<1)|(cut.score.mat[,3] > 0)

        if (any(signeg.status)) {
            cut.score.mat2 <- cut.score.mat[signeg.status, ]

            if (sum(signeg.status)<=1) {
                cut.val <- cut.score.mat2[1]
            } else {
                sigpos.index <- which.min(cut.score.mat2[, 2])
                cut.val <- cut.score.mat2[sigpos.index, 1]
            }
        } else {
            cut.val <- NA
        }
        cut.val
    }

    ###########################################################################
    #
    #                 Prognostic BATTing Functions
    #
    ###########################################################################

    aim.find.cutoff.prog <- function(data, aim.model, marker,yvar, censorvar, type, nsubj, min.sigp.prcnt) {
        # Find optimal cutoff for prognostic case.
        nmax <- length(aim.model$res)
        cut.score.mat <- NULL

        for (nvar.i in 1:nmax){
            index <- index.prediction(aim.model$res[[nvar.i]], marker)
            data$index <- index
            cut.score <- aim.score.prog(data=data, yvar=yvar, censorvar=censorvar, xvar="index", type=type, cutoff=nvar.i-1, nsubj=nsubj, min.sigp.prcnt = min.sigp.prcnt)
            cut.score.mat=rbind(cut.score.mat,c(nvar.i,cut.score))
        }

        sigpos.index <- which.min(cut.score.mat[, 2])
        cut.val <- cut.score.mat[sigpos.index, 1]
        cut.val
    }

    aim.convert <- function(model, nvar, marker.names) {
        aim.model <- model$res[[nvar]]
        var.names <- marker.names[aim.model[, "jmax"]]
        dir <- rep(">", nrow(aim.model))
        dir[aim.model[, "maxdir"]==-1] <- "<"
        aim.out <- data.frame(variable=var.names, direction=dir, cutoff=aim.model[, "cutp"])
        aim.out
    }

    ################### Main aim.rule.batting Code ##################################################

    data <- cbind(y, x)
    xvars <- names(x)
    yvar <- names(y)
    marker.data <- as.matrix(x)
    response <- y[, yvar]

    if (!is.null(trt.vec)) {
        trtvar <- names(trt.vec)
        Trt <- trt.vec[, trtvar]

        if (!is.null(trtref)) {
            # If trtref is provided, create 0-1 vector for trt1 in which 1 corresponds to trtref.
            if (is.element(trtref, unique(Trt))) {
                trt1 <- rep(0,length=length(Trt))
                trt1[Trt==trtref] <- 1
            } else {
                stop("trtref must be one of the entries of trt.vec.")
            }
        } else {
            # trtref not provided, so assume Trt is 0-1 vector.
            if (setequal(union(c(0,1), unique(Trt)), c(0,1))) {
                # Elements of Trt are subset of c(0,1).
                trt1 <- Trt
            } else {
                stop("If trt.vec is not a 0-1 vector, trtref must be supplied.")
            }
        }
        trt.vec[, trtvar] <- trt1
        data <- cbind(data, trt.vec)
        trtref <- 1
    } else {
        trtvar <- NULL
    }

    if (!is.null(censor.vec)) {
        data <- cbind(data, censor.vec)
        censorvar <- names(censor.vec)
        censor <- censor.vec[, censorvar]
    } else {
        censorvar <- NULL
        censor <- NULL
    }

    if (!is.null(pre.filter)){
      xvars_new=filter(data=data,type=type,yvar=yvar,xvars=xvars,censorvar=censorvar,trtvar=trtvar,n.boot=50,cv.iter=15,pre.filter=pre.filter, filter.method=filter.method)
      del_ind=which(!(xvars %in% xvars_new))+1
      data=data[,-1*del_ind]
      x=x[,xvars_new]
      xvars_copy=xvars
      xvars=xvars_new
      marker.data=as.matrix(x)
    }


    nvars <- length(xvars)
    nvar <- rep(0, mc.iter)
    i.mc <- 1
    err.mc <- 0
    min.ianderr <- 10
    max.properr <- 0.9

    if (type=="s" & !is.null(trt.vec)) {

        while (i.mc <= mc.iter) {


            aim.cv.model <- try(cv.cox.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE,status=censor, nsteps=nvars - 1, mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc <- i.mc+1
            } else{
                err.mc <- err.mc+1
            }

            if ((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr) {
                err.mc <- "error"
                break
            }
        }

        if (err.mc != "error") {
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            max.tab <- nvar.unique[nvar.tab==max(nvar.tab)]
            nvar <- max.tab[length(max.tab)]
            aim.model <- cox.interaction(x=marker.data, trt=trt1, delta=censor, y=response, backfit=TRUE, nsteps=nvar, mincut=mincut)
        }
    }


    if (type=="s" & is.null(trt.vec)) {

        while (i.mc <= mc.iter){
            aim.cv.model <- try(cv.cox.main(x=marker.data, status=censor, y=response,backfit=TRUE, nsteps=nvars - 1, mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc=i.mc+1
            } else {
                err.mc=err.mc+1
            }

            if ((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr) {
                err.mc <- "error"
                break
            }
        }

        if (err.mc != "error"){
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            max.tab <- nvar.unique[nvar.tab==max(nvar.tab)]
            nvar <- max.tab[length(max.tab)]
            aim.model <- cox.main(x=marker.data, delta=censor, y=response, backfit=TRUE, nsteps=nvar,mincut=mincut)
        }
    }

    if (type=="b" & !is.null(trt.vec)) {

        while (i.mc <=mc.iter) {
            aim.cv.model <- try(cv.logistic.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE, nsteps=nvars - 1, mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc <- i.mc+1
            } else {
                err.mc <- err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr){
                err.mc <- "error"
                break
            }
        }

        if (err.mc != "error"){
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            max.tab <- nvar.unique[nvar.tab==max(nvar.tab)]
            nvar <- max.tab[length(max.tab)]
            aim.model <- logistic.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE, nsteps=nvar, mincut=mincut)
        }
    }

    if (type=="b" & is.null(trt.vec)) {
        while (i.mc <=mc.iter) {
            aim.cv.model <- try(cv.logistic.main(x=marker.data, y=response, backfit=TRUE, nsteps=nvars - 1, mincut=mincut), silent=TRUE)
            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc <- i.mc+1
            } else {
                err.mc <- err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr){
                err.mc <- "error"
                break
            }
        }

        if (err.mc != "error") {
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            max.tab <- nvar.unique[nvar.tab==max(nvar.tab)]
            nvar <- max.tab[length(max.tab)]
            aim.model <- logistic.main(x=marker.data, y=response, backfit=TRUE, nsteps=nvar, mincut=mincut)
        }
    }

    if (type=="c" & !is.null(trt.vec)) {

        while (i.mc <=mc.iter) {
            aim.cv.model <- try(cv.lm.interaction(x=marker.data, trt=trt1, y=response,backfit=TRUE, nsteps=nvars - 1, mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc <- i.mc+1
            } else {
                err.mc <- err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr){
                err.mc <- "error"
                break
            }
        }

        if (err.mc != "error"){
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            max.tab <- nvar.unique[nvar.tab==max(nvar.tab)]
            nvar <- max.tab[length(max.tab)]
            aim.model <- lm.interaction(x=marker.data, trt=trt1, y=response, backfit=TRUE, nsteps=nvar, mincut=mincut)
        }
    }

    if (type=="c" & is.null(trt.vec)) {

        while (i.mc <= mc.iter) {
            aim.cv.model <- try(cv.lm.main(x=marker.data, y=response, backfit=TRUE, nsteps=nvars - 1, mincut=mincut), silent=TRUE)

            if (class(aim.cv.model) != "try-error") {
                nvar[i.mc] <- aim.cv.model$kmax
                i.mc <- i.mc+1
            } else {
                err.mc <- err.mc+1
            }

            if((i.mc+err.mc)>=min.ianderr & (err.mc/(i.mc+err.mc))>=max.properr){
                err.mc <- "error"
                break
            }
        }

        if (err.mc != "error") {
            nvar.unique <- sort(unique(nvar))
            nvar.tab <- table(nvar)
            max.tab <- nvar.unique[nvar.tab==max(nvar.tab)]
            nvar <- max.tab[length(max.tab)]
            aim.model <- lm.main(x=marker.data, y=response, backfit=TRUE, nsteps=nvar, mincut=mincut)
        }
    }

    if (err.mc != "error") {
        pred.aim <- index.prediction(aim.model$res[[nvar]], marker.data)

        if (!is.null(trt.vec)) {
          dir.temp <- summary(lm(response~trt1*pred.aim))$coefficient[4,1]
        } else {
          dir.temp=cor(response,pred.aim)
        }

        if ((des.res=="smaller" & dir.temp>0)|(des.res=="larger" & dir.temp<0)) {
            for (nvar.temp in 1:nvar){
                model.temp=aim.model$res[[nvar.temp]]
                xvars.temp=xvars[model.temp[,"jmax"]]
                xvars.temp.value=lapply(as.data.frame(cbind(NULL,marker.data[,xvars.temp])),function(x) sort(unique(x)))


                for (i.xvars in 1:length(xvars.temp)){
                    pos.temp=which(xvars.temp.value[[i.xvars]]==model.temp[i.xvars,"cutp"])
                    model.temp[i.xvars,"cutp"]=xvars.temp.value[[i.xvars]][pos.temp+model.temp[i.xvars,"maxdir"]]
                }
                model.temp[,"maxdir"]=-1*model.temp[,"maxdir"]

                aim.model$res[[nvar.temp]]=model.temp
            }

        }


        if (!is.null(trt.vec)) {
            aim.data <- data.frame(response=response, censor=data[, censorvar], Trt=trt.vec[, trtvar],index=NA)
            batting.results <- BATTing(yvar="response", censorvar="censor", trtvar="Trt", trtref=trtref, data=aim.data,
                                       aim.model=aim.model, marker.data=marker.data, type=type, n.boot=n.boot,min.sigp.prcnt=min.sigp.prcnt)
        } else {
            aim.data <- data.frame(response=response, censor=data[, censorvar],index=NA)
            batting.results <- BATTing(yvar="response", censorvar="censor", data=aim.data, aim.model=aim.model,
                                       marker.data=marker.data, type=type, n.boot=n.boot, min.sigp.prcnt=min.sigp.prcnt)
        }
        aim.out <- aim.convert(aim.model, batting.results, xvars)
        aim.model=aim.model$res[[batting.results]]

        if(!is.null(pre.filter)){
          xvars.pos=match(xvars[aim.model[, "jmax"]],xvars_copy)
          aim.model[,"jmax"]=xvars.pos
        }

        output <- list(aim.model=aim.model, aim.rule=aim.out, bat.cutoff=batting.results-1, Note="Sig+ Group: score>cutoff (in this case, satisfy all rules)")
    } else {
        output <- NULL
        print("Error in call to the AIM library.")
    }
    output
}

###################################################################################################


###################################################################################################
#
#  Corrected vesions of AIM package functions.
#
###################################################################################################


########################################### Predictive ############################################


#' A function for CV in linear AIM with interaction.
#'
#' @param x the predictor matrix.
#' @param trt the treatment indicator vector.
#' @param y the vector of the continuous response variable.
#' @param K.cv number of folds for CV.
#' @param num.replicate number of CV iterations.
#' @param nsteps the maximum number of binary rules to be included in the index.
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param backfit a logical argument indicating whether the existing cutpoints are adjusted after including new binary rule.
#' @param maxnumcut the maximum number of binary splits per predictor.
#' @param dirp a vector for pre-specified direction of the binary split for each of the predictors. 0 represents "no pre-given direction"; 1 represents "(x>cut)"; -1 represents "(x<cut)". Alternatively, "dirp=0" represents that there is no pre-given direction for any of the predictor.
#'
#' @return returns optimal number of binary rules based on CV along with CV score test statistics, pre-validated score test statistics and prevalidated fits for individual observation.
#' 
cv.lm.interaction=function (x, trt, y, K.cv = 5, num.replicate = 1, nsteps, mincut = 0.1, backfit = F, maxnumcut = 1, dirp = 0) {
  x = as.matrix(x)
  n = length(y)
  sc.tot = matrix(0, nrow = K.cv, ncol = nsteps)
  preval.tot = matrix(0, nrow = nsteps, ncol = n)
  for (b in 1:num.replicate) {
    g = sample(rep(1:K.cv, (n + K.cv - 1)/K.cv))
    g = g[1:n]
    cva = vector("list", K.cv)
    sc = matrix(0, nrow = K.cv, ncol = nsteps)
    preval = matrix(0, nrow = nsteps, ncol = n)
    for (i in 1:K.cv) {
      cva[[i]] = lm.interaction(x[g != i, ], trt[g != i],
                                y[g != i], nsteps = nsteps, mincut = mincut,
                                backfit = backfit, maxnumcut = maxnumcut, dirp = dirp)
      for (ii in 1:nsteps) {
        aa = index.prediction(cva[[i]]$res[[ii]], x[g == i, ])
        fit = lm(y[g == i] ~ (trt[g == i])+aa : (trt[g == i]))

        if (class(try(summary(fit)$coef[3, 3],silent=T))=="try-error")
        {
          sc[i, ii]<-0
        } else {
          sc[i,ii]<-summary(fit)$coef[3, 3]
        }
        preval[ii, g == i] = aa
      }
    }
    sc.tot = abs(sc.tot) + sc
    preval.tot = preval.tot + preval
  }
  #meansc = colMeans(sc.tot,na.rm=TRUE)
  meansc = apply(sc.tot,2,function(x) sum(x)/sum(x!=0))
  kmax = which.max(meansc)
  pvfit.score = rep(0, nsteps)
  for (i in 1:nsteps) {
    fit = lm(y ~ trt+preval.tot[i, ] : trt)
    pvfit.score[i] = summary(fit)$coef[3, 3]
  }
  if(length(kmax)==0) stop("length(kmax) is 0")
  return(list(kmax = kmax, meanscore = meansc/num.replicate,
              pvfit.score = pvfit.score, preval = preval.tot/num.replicate))
}


###################################################################################################

#' A function for CV in Cox AIM with interaction.
#'
#' @param x the predictor matrix.
#' @param trt the treatment indicator vector.
#' @param y the vector of the time to event response variable.
#' @param status status indicator: 1=failure 0=alive
#' @param K.cv number of folds for CV.
#' @param num.replicate number of CV iterations.
#' @param nsteps the maximum number of binary rules to be included in the index.
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param backfit a logical argument indicating whether the existing cutpoints are adjusted after including new binary rule.
#' @param maxnumcut the maximum number of binary splits per predictor.
#' @param dirp a vector for pre-specified direction of the binary split for each of the predictors. 0 represents "no pre-given direction"; 1 represents "(x>cut)"; -1 represents "(x<cut)". Alternatively, "dirp=0" represents that there is no pre-given direction for any of the predictor.
#'
#' @return returns optimal number of binary rules based on CV along with CV partial likelihood score test statistics, pre-validated partial likelihood score test statistics and prevalidated fits for individual observation.
#' 
cv.cox.interaction=function (x, trt, y, status, K.cv = 5, num.replicate = 1, nsteps, mincut = 0.1, backfit = F, maxnumcut = 1, dirp = 0) {
  x = as.matrix(x)
  n = length(y)
  sc.tot = matrix(0, nrow = K.cv, ncol = nsteps)
  preval.tot = matrix(0, nrow = nsteps, ncol = n)
  for (b in 1:num.replicate) {
    g = sample(rep(1:K.cv, (n + K.cv - 1)/K.cv))
    g = g[1:n]
    cva = vector("list", K.cv)
    sc = matrix(0, nrow = K.cv, ncol = nsteps)
    preval = matrix(0, nrow = nsteps, ncol = n)
    pv.fit = rep(0, nsteps)
    for (i in 1:K.cv) {
      cva[[i]] = cox.interaction(x[g != i, ], trt[g != i], y[g != i], status[g != i], nsteps = nsteps,
                                 mincut = mincut, backfit = backfit, maxnumcut = maxnumcut,
                                 dirp = dirp)
      for (ii in 1:nsteps) {
        aa = index.prediction(cva[[i]]$res[[ii]], x[g == i, ])
        fit = coxph(Surv(y[g == i], status[g == i]) ~ trt[g == i]+aa : trt[g == i])
        if (class(try(summary(fit)$coef[2, 4],silent = T)) == "try-error")
        {
          sc[i, ii] = 0
        } else {
          sc[i, ii] = summary(fit)$coef[2, 4]
        }

        preval[ii, g == i] = aa
      }
    }
    sc.tot = abs(sc.tot) + sc
    preval.tot = preval.tot + preval
  }
#  meansc = colMeans(sc.tot,na.rm=TRUE)
  meansc = apply(sc.tot,2,function(x) sum(x)/sum(x!=0))
  kmax = which.max(meansc)
  pvfit.score = rep(0, nsteps)
  for (i in 1:nsteps) {
    fit = coxph(Surv(y, status) ~ trt+preval.tot[i, ] : trt)
    pvfit.score[i] = summary(fit)$coef[2, 4]
  }
  if(length(kmax)==0) stop("length(kmax) is 0")
  return(list(kmax = kmax, meanscore = meansc/num.replicate,
              pvfit.score = pvfit.score, preval = preval.tot/num.replicate))
}


###################################################################################################

#' A function for CV in logistic AIM with interaction.
#'
#' @param x the predictor matrix.
#' @param trt the treatment indicator vector.
#' @param y the vector of the binary response variable.
#' @param K.cv number of folds for CV.
#' @param num.replicate number of CV iterations.
#' @param nsteps the maximum number of binary rules to be included in the index.
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param backfit a logical argument indicating whether the existing cutpoints are adjusted after including new binary rule.
#' @param maxnumcut the maximum number of binary splits per predictor.
#' @param dirp a vector for pre-specified direction of the binary split for each of the predictors. 0 represents "no pre-given direction"; 1 represents "(x>cut)"; -1 represents "(x<cut)". Alternatively, "dirp=0" represents that there is no pre-given direction for any of the predictor.
#' @param weight a positive value for the weight given to outcomes. "weight=0" means that all observations are equally weighted.
#'
#' @return returns optimal number of binary rules based on CV along with CV score test statistics and pre-validated score test statistics for the treatment*index interaction and prevalidated fits for individual observation.
#' 
cv.logistic.interaction = function (x, trt, y, K.cv = 5, num.replicate = 1, nsteps, mincut = 0.1, backfit = F, maxnumcut = 1, dirp = 0, weight = 1)
{
  x = as.matrix(x)
  n = length(y)
  sc.tot = matrix(0, nrow = K.cv, ncol = nsteps)
  preval.tot = matrix(0, nrow = nsteps, ncol = n)
  for (b in 1:num.replicate) {
    g = sample(rep(1:K.cv, (n + K.cv - 1)/K.cv))
    g = g[1:n]
    cva = vector("list", K.cv)
    sc = matrix(0, nrow = K.cv, ncol = nsteps)
    preval = matrix(0, nrow = nsteps, ncol = n)
    for (i in 1:K.cv) {
      cva[[i]] = logistic.interaction(x[g != i, ], trt[g != i], y[g != i], nsteps = nsteps, mincut = mincut,
                                      backfit = backfit, maxnumcut = maxnumcut, dirp = dirp,
                                      weight = weight)
      for (ii in 1:nsteps) {
        aa = index.prediction(cva[[i]]$res[[ii]], x[g == i, ])
        fit = glm(y[g == i] ~ (trt[g == i])+aa : (trt[g == i]), family = "binomial")
        if (class(try(summary(fit)$coef[3, 3],silent=T))=="try-error")
        {
          sc[i, ii]<-0
        } else {
          sc[i,ii]<-summary(fit)$coef[3, 3]
        }
        preval[ii, g == i] = aa
      }
    }
    sc.tot = abs(sc.tot) + sc
    preval.tot = preval.tot + preval
  }
#  meansc = colMeans(sc.tot, na.rm=TRUE)
  meansc = apply(sc.tot,2,function(x) sum(x)/sum(x!=0))
  kmax = which.max(meansc)
  pvfit.score = rep(0, nsteps)
  for (i in 1:nsteps) {
    fit = glm(y ~ trt+preval.tot[i, ] : trt, family = "binomial")
    pvfit.score[i] = summary(fit)$coef[3, 3]
  }
  if(length(kmax)==0) stop("length(kmax) is 0")
  return(list(kmax = kmax, meanscore = meansc/num.replicate,
              pvfit.score = pvfit.score, preval = preval.tot/num.replicate))
}

###################################################################################################


########################################### Prognostic ############################################

#' A function for the number of binary rules in the main effect AIM with continuous outcome 
#'
#' @param x the predictor matrix.
#' @param y the vector of the continuous response variable.
#' @param K.cv number of folds for CV.
#' @param num.replicate number of CV iterations.
#' @param nsteps the maximum number of binary rules to be included in the index.
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param backfit a logical argument indicating whether the existing cutpoints are adjusted after including new binary rule.
#' @param maxnumcut the maximum number of binary splits per predictor.
#' @param dirp a vector for pre-specified direction of the binary split for each of the predictors. 0 represents "no pre-given direction"; 1 represents "(x>cut)"; -1 represents "(x<cut)". Alternatively, "dirp=0" represents that there is no pre-given direction for any of the predictor.
#'
#' @return returns optimal number of binary rules based on CV along with CV score test statistics for the main effect, pre-validated score test statistics and prevalidated fits for individual observation.
#' 
cv.lm.main=function(x, y, K.cv=5, num.replicate=1, nsteps, mincut=0.1, backfit=F, maxnumcut=1, dirp=0){
  x=as.matrix(x)
  n=length(y)

  sc.tot=matrix(0, nrow=K.cv, ncol=nsteps)
  preval.tot=matrix(0, nrow=nsteps, ncol=n)

  for(b in 1:num.replicate)
  {g=sample(rep(1:K.cv,(n+K.cv-1)/K.cv))
  g=g[1:n]
  cva=vector("list",K.cv)

  sc=matrix(0, nrow=K.cv, ncol=nsteps)
  preval=matrix(0, nrow=nsteps, ncol=n)
  for(i in 1:K.cv)
  {cva[[i]]=lm.main(x[g!=i ,],  y[g!=i], mincut=mincut, nsteps=nsteps, backfit=backfit, maxnumcut=maxnumcut, dirp=dirp)
  for(ii in 1:nsteps)
  {aa=index.prediction(cva[[i]]$res[[ii]],x[g==i,])
  fit=lm(y[g==i]~aa)
  if (class(try(summary(fit)$coef[2, 3],silent=T))=="try-error")
  {
    sc[i, ii]<-0
  } else {
    sc[i,ii]<-summary(fit)$coef[2, 3]
  }
  preval[ii,g==i]=aa
  }
  }
  sc.tot=abs(sc.tot)+sc
  preval.tot=preval.tot+preval
  }
  #meansc=colMeans(sc.tot)
  meansc = apply(sc.tot,2,function(x) sum(x)/sum(x!=0))
  kmax=which.max(meansc)
  pvfit.score=rep(0, nsteps)
  for(i in 1:nsteps)
  {fit=lm(y~preval.tot[i,])
  pvfit.score[i]=summary(fit)$coef[2,3]
  }


  return(list(kmax=kmax, meanscore=meansc/num.replicate, pvfit.score=pvfit.score, preval=preval.tot/num.replicate))
}

###################################################################################################

#' A function for the number of binary rules in the main effect AIM with binary outcome 
#'
#' @param x the predictor matrix.
#' @param y the vector of the binary response variable.
#' @param K.cv number of folds for CV.
#' @param num.replicate number of CV iterations.
#' @param nsteps the maximum number of binary rules to be included in the index.
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param backfit a logical argument indicating whether the existing cutpoints are adjusted after including new binary rule.
#' @param maxnumcut the maximum number of binary splits per predictor.
#' @param dirp a vector for pre-specified direction of the binary split for each of the predictors. 0 represents "no pre-given direction"; 1 represents "(x>cut)"; -1 represents "(x<cut)". Alternatively, "dirp=0" represents that there is no pre-given direction for any of the predictor.
#' @param weight a positive value for the weight given to outcomes. "weight=0" means that all observations are equally weighted.
#'
#' @return returns optimal number of binary rules based on CV along with CV score test statistics for the main effect, pre-validated score test statistics and prevalidated fits for individual observation.
#' 
cv.logistic.main=function(x, y,  K.cv=5, num.replicate=1, nsteps, mincut=0.1, backfit=F, maxnumcut=1, dirp=0, weight=1){
  x=as.matrix(x)
  n=length(y)

  sc.tot=matrix(0, nrow=K.cv, ncol=nsteps)
  preval.tot=matrix(0, nrow=nsteps, ncol=n)

  for(b in 1:num.replicate)
  {g=sample(rep(1:K.cv,(n+K.cv-1)/K.cv))
  g=g[1:n]
  cva=vector("list",K.cv)

  sc=matrix(0, nrow=K.cv, ncol=nsteps)
  preval=matrix(0, nrow=nsteps, ncol=n)
  for(i in 1:K.cv)
  {cva[[i]]=logistic.main(x[g!=i ,],  y[g!=i], mincut=mincut, nsteps=nsteps, backfit=backfit, maxnumcut=maxnumcut, dirp=dirp, weight=weight)
  for(ii in 1:nsteps)
  {aa=index.prediction(cva[[i]]$res[[ii]],x[g==i,])
  fit=glm(y[g==i]~aa, family="binomial")
  sc[i,ii]=summary(fit)$coef[2,3]
  preval[ii,g==i]=aa
  }
  }
  sc.tot=abs(sc.tot)+sc
  preval.tot=preval.tot+preval
  }

  #meansc=colMeans(sc.tot)
  meansc = apply(sc.tot,2,function(x) sum(x)/sum(x!=0))
  kmax=which.max(meansc)
  pvfit.score=rep(0, nsteps)
  for(i in 1:nsteps)
  {fit=glm(y~preval.tot[i,], family="binomial")
  pvfit.score[i]=summary(fit)$coef[2,3]
  }


  return(list(kmax=kmax, meanscore=meansc/num.replicate, pvfit.score=pvfit.score, preval=preval.tot/num.replicate))
}

###################################################################################################

#' A function for the number of binary rules in the main effect AIM with time to event outcome 
#'
#' @param x the predictor matrix.
#' @param y the vector of the time to event response variable.
#' @param status a logical argument vector indicating status of a patient: 1=failure, 0=alive. 
#' @param K.cv number of folds for CV.
#' @param num.replicate number of CV iterations.
#' @param nsteps the maximum number of binary rules to be included in the index.
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param backfit a logical argument indicating whether the existing cutpoints are adjusted after including new binary rule.
#' @param maxnumcut the maximum number of binary splits per predictor.
#' @param dirp a vector for pre-specified direction of the binary split for each of the predictors. 0 represents "no pre-given direction"; 1 represents "(x>cut)"; -1 represents "(x<cut)". Alternatively, "dirp=0" represents that there is no pre-given direction for any of the predictor.
#'
#' @return returns optimal number of binary rules based on CV along with CV partial likelhood score test statistics for the main effect, pre-validated partial likelhood score test statistics and prevalidated fits for individual observation.
#' 
cv.cox.main=function(x, y, status, K.cv=5, num.replicate=1, nsteps, mincut=0.1, backfit=F, maxnumcut=1, dirp=0){
  x=as.matrix(x)
  n=length(y)

  sc.tot=matrix(0, nrow=K.cv, ncol=nsteps)
  preval.tot=matrix(0, nrow=nsteps, ncol=n)

  for(b in 1:num.replicate)
  {g=sample(rep(1:K.cv,(n+K.cv-1)/K.cv))
  g=g[1:n]
  cva=vector("list",K.cv)

  sc=matrix(0, nrow=K.cv, ncol=nsteps)
  preval=matrix(0, nrow=nsteps, ncol=n)
  for(i in 1:K.cv)
  {cva[[i]]=cox.main(x[g!=i ,],  y[g!=i], status[g!=i],mincut=mincut, nsteps=nsteps, backfit=backfit, maxnumcut=maxnumcut, dirp=dirp)
  for(ii in 1:nsteps)
  {aa=index.prediction(cva[[i]]$res[[ii]],x[g==i,])
  fit=coxph(Surv(y[g==i],status[g==i])~aa)
  sc[i,ii]=sign(fit$coef)*sqrt(fit$scor)
  preval[ii,g==i]=aa
  }
  }
  sc.tot=abs(sc.tot)+sc
  preval.tot=preval.tot+preval
  }

  #meansc=colMeans(sc.tot)
  meansc = apply(sc.tot,2,function(x) sum(x)/sum(x!=0))
  kmax=which.max(meansc)
  pvfit.score=rep(0, nsteps)
  for(i in 1:nsteps)
  {fit=coxph(Surv(y,status)~preval.tot[i,])
  pvfit.score[i]=sign(fit$coef)*sqrt(fit$scor)
  }


  return(list(kmax=kmax, meanscore=meansc/num.replicate, pvfit.score=pvfit.score, preval=preval.tot/num.replicate))
}

###################################################################################################

#' Interaction Cox AIM
#'
#' @param x the predictor matrix.
#' @param trt the treatment indicator vector.
#' @param y the vector of the time to event response variable.
#' @param delta status indicator: 1=failure 0=alive
#' @param nsteps the maximum number of binary rules to be included in the index.
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param backfit a logical argument indicating whether the existing cutpoints are adjusted after including new binary rule.
#' @param maxnumcut the maximum number of binary splits per predictor.
#' @param dirp a vector for pre-specified direction of the binary split for each of the predictors. 0 represents "no pre-given direction"; 1 represents "(x>cut)"; -1 represents "(x<cut)". Alternatively, "dirp=0" represents that there is no pre-given direction for any of the predictor.
#'
#' @return "cox.interaction" returns "maxsc", which is the observed partial likelihood score test statistics for the index*treatment interaction in the fitted model and "res", which is a list with components
#' 
cox.interaction=function (x, trt, y, delta, nsteps = 8, mincut = 0.1, backfit = F,
                          maxnumcut = 1, dirp = 0) {
  n = length(y)
  id = order(y)
  y = y[id]
  x = x[id, ]
  delta = delta[id]
  trt = trt[id]
  x0 = x
  p.true = length(x0[1, ])
  x = cbind(x0, -x0)
  marker.x = rep(1:p.true, 2)
  direction.x = c(rep(-1, p.true), rep(1, p.true))
  if (sum(abs(dirp)) > 0) {
    index1 = (1:p.true)[dirp == -1 | dirp == 0]
    index2 = (1:p.true)[dirp == 1 | dirp == 0]
    x = cbind(x0[, index1], x0[, index2])
    marker.x = c((1:p.true)[index1], (1:p.true)[index2])
    direction.x = c(rep(-1, length(index1)), rep(1, length(index2)))
  }
  p = length(x[1, ])
  if (mincut > 0) {
    effect.range = ceiling(n * mincut):floor(n * (1 - mincut))
  }
  if (mincut == 0) {
    effect.range = 1:n
  }
  score.current = rep(0, n)
  x.out = x
  ntime = length(unique(y))
  time.index = rep(c(1, cumsum(table(y)) + 1)[-(ntime + 1)],
                   table(y))
  num.risk = (n:1)[time.index]
  time.index.plus = rep(cumsum(table(y)), table(y))
  id.in = NULL
  id.out = 1:p.true
  p.out = p
  zvalue = NULL
  cut.value = NULL
  direction = NULL
  imax = NULL
  res = as.list(NULL)
  num.cut = rep(0, p.true)
  fit = coxph(Surv(y, delta) ~ trt, method = "breslow")
  beta = fit$coef
  sigma0 = fit$var
  risk = exp(beta * trt)
  S0 = (rev(cumsum(rev(risk))))[time.index]
  Sr = (rev(cumsum(rev(risk * trt))))[time.index]
  order.index = risk.pool = delta.pool = trt.pool = x.replicate = matrix(0,
                                                                         n, p)
  for (i in 1:p) {
    idx = order(x[, i])
    order.index[, i] = idx
    risk.pool[, i] = risk[idx]
    delta.pool[, i] = delta[idx]
    trt.pool[, i] = trt[idx]
    x.replicate[, i] = rep((1:n)[cumsum(table(x[, i]))],
                           table(x[, i]))
  }
  ds0 = (cumsum(delta/S0))[time.index.plus]
  drs0 = (cumsum(delta * Sr/S0^2))[time.index.plus]
  dss0 = (cumsum(delta/S0^2))[time.index.plus]
  ds0.pool = drs0.pool = dss0.pool = matrix(0, n, p)
  for (i in 1:p) {
    idx = order(x[, i])
    ds0.pool[, i] = ds0[idx]
    dss0.pool[, i] = dss0[idx]
    drs0.pool[, i] = drs0[idx]
  }
  subject = matrix(0, n, n)
  for (i in 1:n) {
    subject[(n - num.risk[i] + 1):n, i] = 1
  }
  i2.diff.pool = apply(trt.pool^2 * risk.pool^2 * dss0.pool,
                       2, cumsum)
  i3.diff.pool = apply(trt.pool^2 * risk.pool * ds0.pool, 2,
                       cumsum)
  i4.diff.pool = apply(trt.pool * risk.pool * drs0.pool, 2,
                       cumsum)
  score.diff.pool = apply(delta.pool * trt.pool, 2, cumsum) -
    apply(trt.pool * risk.pool * ds0.pool, 2, cumsum)
  w1.pool = matrix(0, n, p.out)
  for (i in 1:p.out) {
    idx = order.index[, i]
    subject.id = subject[idx, ]
    temp1 = apply(risk[idx] * trt[idx] * subject.id, 2, cumsum)
    temp1 = rbind(0, temp1[-n, ])
    w1.pool[, i] = diag((apply(t(temp1) * delta/S0^2, 2,
                               cumsum))[time.index.plus, ][idx, ])
  }
  count = 1
  loop = 1
  while (loop == 1) {
    Swwrr0 = (rev(cumsum(rev(risk * trt^2 * score.current^2))))[time.index]
    Swrr0 = (rev(cumsum(rev(risk * trt^2 * score.current))))[time.index]
    Swr0 = (rev(cumsum(rev(risk * trt * score.current))))[time.index]
    score.pool = w.pool = matrix(0, n, p.out)
    for (i in 1:p.out) {
      idx = order.index[, i]
      score.pool[, i] = score.current[idx]
      temp = (cumsum((rev(cumsum(rev(risk * trt * score.current))))[time.index] *
                       delta/S0^2))[time.index.plus]
      w.pool[, i] = w1.pool[, i] + temp[idx]
    }
    i1 = sum(delta * Swwrr0/S0) + apply((2 * score.pool +
                                           1) * trt.pool^2 * risk.pool * ds0.pool, 2, cumsum)
    i2 = sum(delta * Swr0^2/S0^2) + apply(2 * trt.pool *
                                            risk.pool * w.pool, 2, cumsum) + i2.diff.pool
    i3 = sum(delta * Swrr0/S0) + i3.diff.pool
    i4 = sum(delta * Swr0 * Sr/S0^2) + i4.diff.pool
    v.stat = (i1 - i2) - (i3 - i4)^2 * sigma0[1, 1]
    score.stat = sum(delta * score.current * trt) - sum(delta *
                                                          Swr0/S0) + score.diff.pool
    for (i in 1:p.out) {
      idx = x.replicate[, i]
      v.stat[, i] = v.stat[idx, i]
      score.stat[, i] = score.stat[idx, i]
    }
    test = score.stat/sqrt(v.stat + 1e-08)
    if (mincut == 0) {
      mtest = apply(test, 2, max)
    }
    if (mincut > 0) {
      mtest = rep(0, p.out)
      for (i in 1:p.out) {
        test.post = test[max(effect.range) + 1, i]
        test0 = test[effect.range, i]
        test0 = test0[test0 != test.post]
        mtest[i] = max(test0)
      }
    }
    mcut = rep(0, p.out)
    i0 = rep(0, p.out)
    for (i in 1:p.out) {
      i0[i] = max(effect.range[test[effect.range, i] ==
                                 mtest[i]]) + 1
      if (is.infinite(i0[i])) i0[i]=n+1
      if (i0[i] > n)
        mcut[i] = Inf
      if (i0[i] <= n)
        mcut[i] = x.out[order.index[, i], i][i0[i]]
    }
    i.sel = which.max(mtest)
    score.current = score.current + (x.out[, i.sel] < mcut[i.sel])
    marker.sel = marker.x[i.sel]
    id.in = c(id.in, marker.sel)
    direction = c(direction, direction.x[i.sel])
    cut.value = c(cut.value, -mcut[i.sel] * direction.x[i.sel])
    num.cut[marker.sel] = num.cut[marker.sel] + 1
    if (num.cut[marker.sel] == maxnumcut) {
      id.exclude = (1:p.out)[marker.x == marker.sel]
      x.out = x.out[, -id.exclude, drop = F]
      order.index = order.index[, -id.exclude, drop = F]
      trt.pool = trt.pool[, -id.exclude, drop = F]
      risk.pool = risk.pool[, -id.exclude, drop = F]
      ds0.pool = ds0.pool[, -id.exclude, drop = F]
      i2.diff.pool = i2.diff.pool[, -id.exclude, drop = F]
      i3.diff.pool = i3.diff.pool[, -id.exclude, drop = F]
      i4.diff.pool = i4.diff.pool[, -id.exclude, drop = F]
      score.diff.pool = score.diff.pool[, -id.exclude,
                                        drop = F]
      w1.pool = w1.pool[, -id.exclude, drop = F]
      x.replicate = x.replicate[, -id.exclude, drop = F]
      direction.x = direction.x[-id.exclude]
      marker.x = marker.x[-id.exclude]
      p.out = p.out - length(id.exclude)
    }
    if (backfit == T) {
      if (count > 1) {
        x.adj = x0[, id.in]
        x.adj = -t(t(x.adj) * direction)
        cutp = -cut.value * direction
        fit = backfit.cox.interaction(x.adj, trt, y,
                                      delta, cutp, mincut = mincut)
        jmax = id.in
        cutp = -fit$cutp * direction
        maxdir = direction
        res[[count]] = cbind(jmax, cutp, maxdir)
        zvalue = c(zvalue, max(fit$zscore))
        score.current = apply((t(x.adj) < fit$cutp),
                              2, sum)
      }
      if (count == 1) {
        jmax = id.in
        cutp = cut.value
        maxdir = direction
        res[[count]] = cbind(jmax, cutp, maxdir)
        zvalue = c(zvalue, max(mtest))
      }
    }
    if (backfit == F) {
      cutp = cut.value
      maxdir = direction
      jmax = id.in
      zvalue = c(zvalue, max(mtest))
      maxsc = zvalue
      res[[count]] = cbind(jmax, cutp, maxdir, maxsc)
    }
    count = count + 1
    loop = (length(id.in) < nsteps)
  }
  return(list(res = res, maxsc = zvalue))
}


#' An internal function used in cox.interaction
#'
#' @param x the predictor matrix.
#' @param trt the treatment indicator vector. 
#' @param y the vector of the time to event response variable.
#' @param delta status indicator: 1=failure 0=alive
#' @param cutp a specific cutpoint
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' 
backfit.cox.interaction=function (x, trt, y, delta, cutp, mincut = 0) {
  n = length(y)
  p = length(x[1, ])
  id = order(y)
  y = y[id]
  delta = delta[id]
  trt = trt[id]
  x = x[id, ]
  if (mincut > 0) {
    effect.range = ceiling(n * mincut):floor(n * (1 - mincut))
  }
  if (mincut == 0) {
    effect.range = 1:n
  }
  ntime = length(unique(y))
  time.index = rep(c(1, cumsum(table(y)) + 1)[-(ntime + 1)],
                   table(y))
  num.risk = (n:1)[time.index]
  time.index.plus = rep(cumsum(table(y)), table(y))
  fit = coxph(Surv(y, delta) ~ trt, method = "breslow")
  beta = fit$coef
  sigma0 = fit$var
  risk = exp(beta * trt)
  S0 = (rev(cumsum(rev(risk))))[time.index]
  Sr = (rev(cumsum(rev(risk * trt))))[time.index]
  order.index = risk.pool = delta.pool = trt.pool = x.replicate = matrix(0,
                                                                         n, p)
  for (i in 1:p) {
    idx = order(x[, i])
    order.index[, i] = idx
    risk.pool[, i] = risk[idx]
    delta.pool[, i] = delta[idx]
    trt.pool[, i] = trt[idx]
    x.replicate[, i] = rep((1:n)[cumsum(table(x[, i]))],
                           table(x[, i]))
  }
  ds0 = (cumsum(delta/S0))[time.index.plus]
  drs0 = (cumsum(delta * Sr/S0^2))[time.index.plus]
  dss0 = (cumsum(delta/S0^2))[time.index.plus]
  ds0.pool = drs0.pool = dss0.pool = matrix(0, n, p)
  for (i in 1:p) {
    idx = order(x[, i])
    ds0.pool[, i] = ds0[idx]
    dss0.pool[, i] = dss0[idx]
    drs0.pool[, i] = drs0[idx]
  }
  score.mat = t((t(x) < cutp) * 1)
  score.tot = apply(score.mat, 1, sum)
  Swwrr0 = (rev(cumsum(rev(risk * trt^2 * score.tot^2))))[time.index]
  Swrr0 = (rev(cumsum(rev(risk * trt^2 * score.tot))))[time.index]
  Swr0 = (rev(cumsum(rev(risk * trt * score.tot))))[time.index]
  i1 = sum(delta * Swwrr0/S0)
  i2 = sum(delta * Swr0^2/S0^2)
  i3 = sum(delta * Swrr0/S0)
  i4 = sum(delta * Swr0 * Sr/S0^2)
  v.stat = (i1 - i2) - (i3 - i4)^2 * sigma0[1, 1]
  score.stat = sum(delta * score.tot * trt) - sum(delta * Swr0/S0)
  zscore = score.stat/sqrt(v.stat + 1e-08)
  subject = matrix(0, n, n)
  for (i in 1:n) {
    subject[(n - num.risk[i] + 1):n, i] = 1
  }
  i2.add = apply(trt.pool^2 * risk.pool^2 * dss0.pool, 2, cumsum)
  i3.pool = apply(trt.pool^2 * risk.pool * ds0.pool, 2, cumsum)
  i4.pool = apply(trt.pool * risk.pool * drs0.pool, 2, cumsum)
  w1.pool = matrix(0, n, p)
  for (i in 1:p) {
    idx = order.index[, i]
    subject.id = subject[idx, ]
    temp1 = apply(risk[idx] * trt[idx] * subject.id, 2, cumsum)
    temp1 = rbind(0, temp1[-n, ])
    w1.pool[, i] = diag((apply(t(temp1) * delta/S0^2, 2,
                               cumsum))[time.index.plus, ][idx, ])
  }
  score.stat.pool = apply(delta.pool * trt.pool, 2, cumsum) -
    apply(trt.pool * risk.pool * ds0.pool, 2, cumsum)
  count = 1
  loop = 1
  while (loop == 1) {
    score.mat = t((t(x) < cutp) * 1)
    score.tot = apply(score.mat, 1, sum)
    score.pool = score.pool0 = matrix(0, n, p)
    for (i in 1:p) {
      idx = order.index[, i]
      score.pool[, i] = (score.tot - score.mat[, i])[idx]
      score.pool0[, i] = (score.tot - score.mat[, i])
    }
    Swwrr0 = ((apply((risk * trt^2 * score.pool0^2)[n:1,
                                                    ], 2, cumsum))[n:1, ])[time.index, ]
    Swrr0 = ((apply((risk * trt^2 * score.pool0)[n:1, ],
                    2, cumsum))[n:1, ])[time.index, ]
    Swr0 = ((apply((risk * trt * score.pool0)[n:1, ], 2,
                   cumsum))[n:1, ])[time.index, ]
    w.pool = matrix(0, n, p)
    for (i in 1:p) {
      idx = order.index[, i]
      temp2 = (cumsum((rev(cumsum(rev(risk * trt * score.pool0[,
                                                               i]))))[time.index] * delta/S0^2))[time.index.plus]
      w.pool[, i] = w1.pool[, i] + temp2[idx]
    }
    i10 = apply(delta * Swwrr0/S0, 2, sum)
    i20 = apply(delta * Swr0^2/S0^2, 2, sum)
    i30 = apply(delta * Swrr0/S0, 2, sum)
    i40 = apply(delta * Swr0 * Sr/S0^2, 2, sum)
    i1.pool = apply((2 * score.pool + 1) * trt.pool^2 * risk.pool *
                      ds0.pool, 2, cumsum)
    i2.pool = apply(2 * trt.pool * risk.pool * w.pool, 2,
                    cumsum) + i2.add
    i1 = t(t(i1.pool) + i10)
    i2 = t(t(i2.pool) + i20)
    i3 = t(t(i3.pool) + i30)
    i4 = t(t(i4.pool) + i40)
    v.stat = (i1 - i2) - (i3 - i4)^2 * sigma0[1, 1]
    score.stat0 = apply(delta.pool * score.pool * trt.pool,
                        2, sum) - apply(delta * Swr0/S0, 2, sum)
    score.stat = t(t(score.stat.pool) + score.stat0)
    for (i in 1:p) {
      idx = x.replicate[, i]
      v.stat[, i] = v.stat[idx, i]
      score.stat[, i] = score.stat[idx, i]
    }
    test = score.stat/sqrt(v.stat + 1e-08)
    if (mincut == 0) {
      mtest = apply(test, 2, max)
    }
    if (mincut > 0) {
      mtest = rep(0, p)
      for (i in 1:p) {
        test.post = test[max(effect.range) + 1, i]
        test0 = test[effect.range, i]
        test0 = test0[test0 != test.post]
        mtest[i] = max(test0)
      }
    }
    loop = 0
    if (max(mtest) > zscore[count] + 1e-08) {
      loop = 1
      i = (1:p)[mtest == max(mtest)][1]
      i0 = max(effect.range[test[effect.range, i] == mtest[i]]) +
        1
      if (i0 > n)
        mcut = Inf
      if (i0 <= n)
        mcut = x[order.index[, i], i][i0]
      cutp[i] = mcut
      zscore = c(zscore, max(mtest))
      count = count + 1
    }
  }
  return(list(cutp = cutp, zscore = zscore))
}


###################################################################################################

###################################################################################################
#
#   CV code for AIM BATTing
#
#####################################################################################################################

#' The function for CV in aim.batting
#' @description Implements k-fold cross validation for aim.batting.
#'
#' @param y data frame containing the response
#' @param x data frame containing the predictor
#' @param censor.vec data frame giving the censor status (only for TTE data , censor=0,event=1)  - default = NULL
#' @param trt.vec data frame giving the censor status (only for TTE data , censor=0,event=1)  - default = NULL
#' @param trtref treatment reference indicator: 1=treatment, 0=control
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c"
#' @param n.boot number of bootstraps in bootstrapping step.
#' @param des.res the desired response. "larger": prefer larger response. "smaller": prefer smaller response
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#' @param mc.iter # of iterations for the MC procedure to get a stable "best number of predictors"
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering
#' @param k.fold # cross-validation folds
#' @param cv.iter Algotithm terminates after cv.iter successful iterations of cross-validation
#' @param max.iter total # iterations (including unsuccessful) allowed.
#'
#' @return "cv.aim.batting" returns a list with following entries: 
#' \item{stats.summary}{Summary of performance statistics.}
#' \item{pred.classes}{Data frame containing the predictive clases (TRUE/FALSE) for each iteration.}
#' \item{folds}{Data frame containing the fold indices (index of the fold for each row) for each iteration.}
#' \item{sig.list}{List of length cv.iter * k.fold containing the signature generated at each of the  k folds, for all iterations.}
#' \item{error.log}{List of any error messages that are returned at an iteration.}
#' \item{interplot}{Treatment*subgroup interaction plot for predictive case}
#'    

cv.aim.batting <- function(y, x, censor.vec=NULL, trt.vec=NULL, trtref=NULL, type, n.boot, des.res="larger", min.sigp.prcnt=0.20, mc.iter=1, mincut=0.1, pre.filter=NULL, filter.method=NULL, k.fold=5, cv.iter=50, max.iter=500) {  #Y#y

    y <- as.data.frame(y)
    x <- as.data.frame(x)

    if (nrow(x) != nrow(y)) {
        stop("Error: x and y have different numbers of rows.")
    }
    yvar <- names(y)
    xvars <- names(x)
    data <- cbind(y, x)

    if (!is.null(censor.vec)) {
        censor.vec <- as.data.frame(censor.vec)
        censorvar <- names(censor.vec)
        if (nrow(censor.vec) != nrow(x)) {
            stop("Error: censor.vec and x have different numbers of rows.")
        }
        data <- cbind(data, censor.vec)
    } else {
        censorvar <- NULL
    }

    if (!is.null(trt.vec)) {
        trtvar <- names(trt.vec)
        Trt <- trt.vec[, trtvar]

        if (!is.null(trtref)) {
            # If trtref is provided, create 0-1 vector for trt1 in which 1 corresponds to trtref.
            if (is.element(trtref, unique(Trt))) {
                trt1 <- rep(0,length=length(Trt))
                trt1[Trt==trtref] <- 1
            } else {
                stop("trtref must be one of the entries of trt.vec.")
            }
        } else {
            # trtref not provided, so assume Trt is 0-1 vector.
            if (setequal(union(c(0,1), unique(Trt)), c(0,1))) {
                # Elements of Trt are subset of c(0,1).
                trt1 <- Trt
            } else {
                stop("If trt.vec is not a 0-1 vector, trtref must be supplied.")
            }
        }
        trt.vec[, trtvar] <- trt1
        data <- cbind(data, trt.vec)
        trtref <- 1
    } else {
        trtvar <- NULL
    }

    if (type=="b") strata <- data[,yvar]
    if (type=="s") strata <- censor.vec[, censorvar]
    if (type=="c") strata <- NULL

    model.Rfunc <- "aim.batting.wrapper"
    model.Rfunc.args <- list(yvar=yvar, xvars=xvars, censorvar=censorvar, trtvar=trtvar, trtref=trtref, type=type,
                             n.boot=n.boot, des.res=des.res, min.sigp.prcnt=min.sigp.prcnt, mc.iter=mc.iter, mincut=mincut, pre.filter=pre.filter, filter.method=filter.method)
    # List of arguments for model.Rfunc. Must be packaged into a list to pass to kfold.cv.

    predict.Rfunc <- "pred.aim.cv"
    predict.Rfunc.args <- list(xvars=xvars)
    # List of arguments for predict.Rfunc.

    res <- kfold.cv(data=data, model.Rfunc=model.Rfunc, model.Rfunc.args=model.Rfunc.args, predict.Rfunc=predict.Rfunc,
                    predict.Rfunc.args=predict.Rfunc.args, k.fold=k.fold, cv.iter=cv.iter, strata=strata, max.iter=max.iter)

    if (length(res) > 0) {
        stats <- evaluate.cv.results(cv.data=res$cv.data, y=y, censor.vec, trt.vec=trt.vec, type=type)
        summary <- summarize.cv.stats(stats$raw.stats, trtvar, type)
        interplot=interaction.plot(data.eval=summary, type=type, main="Interaction Plot", trt.lab=c("Trt.", "Ctrl."))
        results <- list(stats.summary=summary, pred.classes=stats$pred.classes, folds=stats$folds, sig.list=res$sig.list, raw.stats=stats$raw.stats, error.log=res$error.log, interplot=interplot)
    } else {
        results <- "No successful cross-validations."
    }

    return(results)
}

###################################################################################################

#' Wrapper function for cv.aim.batting to be passed to kfold.cv.
#'
#' @param data data frame equal to cbind(y, x), where y and x are inputs to aim.batting.
#' @param args list containing all other input arguments to aim.batting except for x and y.
#'
#' @return prediction rule as returned by aim.batting.
#' 
aim.batting.wrapper <- function(data, args) {
  
  yvar=NULL; xvars=NULL; trtref=NULL; type=NULL; n.boot=0; des.res=NULL; min.sigp.prcnt=NULL; mc.iter=NULL; mincut=NULL; pre.filter=NULL; filter.method=NULL;
  
    for (name in names(args)) {
        assign(name, args[[name]])
    }
    y <- subset(data, select=yvar)
    x <- subset(data, select=xvars)

    if (!is.null(args$censorvar)) {
        censorvar <- args$censorvar
        censor.vec <- subset(data, select=censorvar)
    } else {
        censor.vec <- NULL
    }

    if (!is.null(args$trtvar)) {
        trtvar <- args$trtvar
        trt.vec <- subset(data, select=trtvar)
    } else {
        trt.vec <- NULL
    }
    res <- aim.batting(y, x, censor.vec=censor.vec, trt.vec=trt.vec, trtref=trtref, type=type, n.boot=n.boot, des.res=des.res, min.sigp.prcnt=min.sigp.prcnt, mc.iter=mc.iter, mincut=mincut, pre.filter=pre.filter, filter.method=filter.method)#Y#y
    return(res)
}

###################################################################################################
#
#  CV code for AIM-Rule BATTing
#
#####################################################################################################################

#' The function for CV in aim.rule.batting
#' @description Implements k-fold cross validation for aim.batting.
#' 
#' @param y data frame containing the response
#' @param x data frame containing the predictor
#' @param censor.vec data frame giving the censor status (only for TTE data , censor=0,event=1)  - default = NULL
#' @param trt.vec data frame giving the censor status (only for TTE data , censor=0,event=1)  - default = NULL
#' @param trtref treatment reference indicator: 1=treatment, 0=control
#' @param type data type - "c" - continuous , "b" - binary, "s" - time to event  - default = "c"
#' @param n.boot number of bootstraps in bootstrapping step.
#' @param des.res the desired response. "larger": prefer larger response. "smaller": prefer smaller response
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#' @param mc.iter # of iterations for the MC procedure to get a stable "best number of predictors"
#' @param mincut the minimum cutting proportion for the binary rule at either end. It typically is between 0 and 0.2. It is the parameter in the functions of AIM package.
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering
#' @param k.fold # cross-validation folds
#' @param cv.iter Algotithm terminates after cv.iter successful iterations of cross-validation
#' @param max.iter total # iterations (including unsuccessful) allowed.
#'
#' @return "cv.aim.batting" returns a list with following entries:
#' @return stats.summary: Summary of performance statistics
#' @return pred.classes: Data frame containing the predictive clases (TRUE/FALSE) for each iteration.
#' @return pred.classes: Data frame containing the predictive clases (TRUE/FALSE) for each iteration.
#' @return folds: Data frame containing the fold indices (index of the fold for each row) for each iteration
#' @return sig.list: List of length cv.iter * k.fold containing the signature generated at each of the  k folds, for all iterations.
#' @return error.log: List of any error messages that are returned at an iteration.
#'   
cv.aim.rule.batting <- function(y, x, censor.vec=NULL, trt.vec=NULL, trtref=NULL, type, n.boot, des.res="larger", min.sigp.prcnt=0.2, mc.iter=1, mincut=0.1, pre.filter=NULL, filter.method=NULL, k.fold=5, cv.iter=50, max.iter=500) {
  
    y <- as.data.frame(y)
    x <- as.data.frame(x)
    if (nrow(x) != nrow(y)) {
        stop("Error: x and y have different numbers of rows.")
    }
    yvar <- names(y)
    xvars <- names(x)
    data <- cbind(y, x)
    if (!is.null(censor.vec)) {
        censor.vec <- as.data.frame(censor.vec)
        censorvar <- names(censor.vec)
        if (nrow(censor.vec) != nrow(x)) {
            stop("Error: censor.vec and x have different numbers of rows.")
        }
        data <- cbind(data, censor.vec)
    } else {
        censorvar <- NULL
    }

    if (!is.null(trt.vec)) {
        trtvar <- names(trt.vec)
        Trt <- trt.vec[, trtvar]

        if (!is.null(trtref)) {
            # If trtref is provided, create 0-1 vector for trt1 in which 1 corresponds to trtref.
            if (is.element(trtref, unique(Trt))) {
                trt1 <- rep(0,length=length(Trt))
                trt1[Trt==trtref] <- 1
            } else {
                stop("trtref must be one of the entries of trt.vec.")
            }
        } else {
            # trtref not provided, so assume Trt is 0-1 vector.
            if (setequal(union(c(0,1), unique(Trt)), c(0,1))) {
                # Elements of Trt are subset of c(0,1).
                trt1 <- Trt
            } else {
                stop("If trt.vec is not a 0-1 vector, trtref must be supplied.")
            }
        }
        trt.vec[, trtvar] <- trt1
        data <- cbind(data, trt.vec)
        trtref <- 1
    } else {
        trtvar <- NULL
    }

    if (type=="b") strata <- data[,yvar]
    if (type=="s") strata <- censor.vec[, censorvar]
    if (type=="c") strata <- NULL

    model.Rfunc <- "aim.rule.batting.wrapper"
    model.Rfunc.args <- list(yvar=yvar, xvars=xvars, censorvar=censorvar, trtvar=trtvar, trtref=trtref, type=type, n.boot=n.boot,
                             des.res=des.res, min.sigp.prcnt=min.sigp.prcnt, mc.iter=mc.iter, mincut=mincut, pre.filter=pre.filter, filter.method=filter.method) #Y#y
    # List of arguments for model.Rfunc. Must be packaged into
    # a list to pass to kfold.cv.

    predict.Rfunc <- "pred.aim.cv"
    predict.Rfunc.args <- list(xvars=xvars)
    # List of arguments for predict.Rfunc.
    res <- kfold.cv(data=data, model.Rfunc=model.Rfunc, model.Rfunc.args=model.Rfunc.args, predict.Rfunc=predict.Rfunc,
                    predict.Rfunc.args=predict.Rfunc.args, k.fold=k.fold, cv.iter=cv.iter, strata=strata, max.iter=max.iter)
    if (length(res) > 0) {
        stats <- evaluate.cv.results(cv.data=res$cv.data, y=y, censor.vec=censor.vec, trt.vec=trt.vec, type=type)
        summary <- summarize.cv.stats(stats$raw.stats, trtvar, type)
        interplot=interaction.plot(data.eval=summary, type=type, main="Interaction Plot", trt.lab=c("Trt.", "Ctrl."))
        results <- list(stats.summary=summary, pred.classes=stats$pred.classes, folds=stats$folds, sig.list=res$sig.list, raw.stats=stats$raw.stats, error.log=res$error.log, interplot=interplot)

        } else {
        results <- "No successful cross-validations."
    }

    return(results)
}

###################################################################################################

#' Wrapper function for aim.rule.batting, to be passed to kfold.cv
#'
#' @param data data frame equal to cbind(y, x), where y and x are inputs to aim.rule.batting.
#' @param args list containing all other input arguments to aim.rule.batting except for x and y.
#'
#' @return prediction rule returned by aim.rule.batting.
#' 
aim.rule.batting.wrapper <- function(data, args) {
  
  yvar=NULL; xvars=NULL; trtref=NULL; type=NULL; n.boot=0; des.res=NULL; min.sigp.prcnt=NULL; mc.iter=NULL; mincut=NULL; pre.filter=NULL; filter.method=NULL;
  
    
    for (name in names(args)) {
        assign(name, args[[name]])
    }
    y <- subset(data, select=yvar)
    x <- subset(data, select=xvars)
    if (!is.null(args$censorvar)) {
        censorvar <- args$censorvar
        censor.vec <- subset(data, select=censorvar)
    } else {
        censor.vec <- NULL
    }
    if (!is.null(args$trtvar)) {
        trtvar <- args$trtvar
        trt.vec <- subset(data, select=trtvar)
    } else {
        trt.vec <- NULL
    }
    res <- aim.rule.batting(y, x, censor.vec=censor.vec, trt.vec=trt.vec, trtref=trtref, type=type, n.boot=n.boot, des.res=des.res, min.sigp.prcnt = min.sigp.prcnt, mc.iter=mc.iter, mincut=mincut, pre.filter=pre.filter, filter.method=filter.method)#Y#y
    return(res)
}

###################################################################################################

#' Make predictions for data given prediction rule in predict.rule
#'
#' @param data Data frame of form cbind(y, x), where y and x are inputs to cv.seq.batting.
#' @param predict.rule Prediction rule returned by seq.batting.
#' @param args list of the form list(xvar=xvar, yvar=yvar)
#'
#' @return The input data with an added column, a logical vector indicating the prediction for each row of data.
#' 
pred.aim.cv <-function(data, predict.rule, args) {

    xvars <- args$xvars
    aim.model <- predict.rule$aim.model
    bat.cutoff <- predict.rule$bat.cutoff
    marker.data <- data[, xvars]
    pred.aim <- index.prediction(aim.model, marker.data)
    # pred.aim is a vector of scores. Each score is the number of individual threshold conditions that are
    # true in the data for xvars.
    pred.id <- pred.aim > bat.cutoff
    pred.data <- cbind(marker.data, data.frame(pred.class=pred.id))
    # Note!!: If you don't cbind pred.class with the data in the current fold, the
    # row numbers of the current fold will be lost.
    pred.data <- subset(pred.data, select="pred.class")
    # Now that the row numbers are attached, you can get rid of marker.data.
    pred.data
}

#' Make predictions for data given prediction rule in predict.rule.
#'
#' @param x the predictor matrix.
#' @param predict.rule prediction rule returned by seq.batting.
#'
#' @return The input data with an added column, a logical vector indicating the prediction for each row of data.
#' 
pred.aim <-function(x, predict.rule) {
    aim.model <- predict.rule$aim.model
    bat.cutoff <- predict.rule$bat.cutoff
    pred.aim <- index.prediction(aim.model, x)
    # pred.aim is a vector of scores. Each score is the number of individual threshold conditions that are
    # true in the data for xvars.
    pred.id <- pred.aim > bat.cutoff
    pred.data<-data.frame(x, pred.class=pred.id)
    pred.data
}


###################################################################################################

#' Get counts of signature variables from output of cv.aim.batting.
#'
#' @param sig.list signature list output from cv.aim.batting.
#' @param xvars predictor variable names.
#'
#' @return counts of signature variables.
#' 
get.var.counts.aim <- function(sig.list, xvars) {
    sig.var.counts <- data.frame(xvars)
    sig.var.counts$count <- rep(0, length(xvars))

    for (i in 1:length(sig.list)) {
        predict.rule <- sig.list[[i]]$aim.rule
        n.rules <- nrow(predict.rule)
        if (is.null(n.rules)) {
            # Can this happen?
            var <- predict.rule["variable"]
            dir <- predict.rule["direction"]
            sig.var.counts[xvars==var, ]$count <- sig.var.counts[xvars==var, ]$count + 1
        } else {
            for(i in 1:n.rules) {
                var <- predict.rule[i, "variable"]
                dir <- predict.rule[i, "direction"]
                sig.var.counts[xvars==var, ]$count <- sig.var.counts[xvars==var, ]$count + 1
            }
        }
    }

        sig.var.counts
}

#################################################################################################################################
#
# Sequential BATTing
# Version 13.2
# Nov. 2013
#
#################################################################################################################################
#
#  Functions Shared by seq.batting and seqlr.batting
#
###################################################################################################
#
#                  Predictive BATTing Functions 
#

#' Compute score of cutoff for predictive case
#'
#' @param data input data frame.
#' @param yvar response variable name.
#' @param censorvar censoring variable name.
#' @param xvar name of predictor for which cutpoint needs to be obtained.
#' @param trtvar treatment variable name.
#' @param cutoff a specific cutpoint for which the score needs to be computed.
#' @param type "c" continuous; "s" survival; "b" binary.
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1).
#' @param dir direction of cut.
#' @param nsubj number of subjects.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#'
#' @return score (p-value of treatment*subgroup interaction) for the given cutoff.
#' 
seqlr.score.pred <- function(data, yvar, censorvar, xvar, trtvar, cutoff, type, class.wt, dir, nsubj, min.sigp.prcnt) {
    # Compute score of cutoff for predictive case.
    # score is the p-value of model.fit
    
    if (dir==">") {
        id <- data[, xvar]>cutoff
    } else {
        id <- data[, xvar]<cutoff
    }
  n.grp <- table(id) 
  sigp.prcnt<-sum(id)/nsubj
  data=cbind(data,id)
  
    if (type=="s") {
        model.text <- paste("Surv(", yvar, ", ", censorvar, ")~", trtvar,"+",trtvar, "*id",sep="")
        model.formula <- as.formula(model.text)
        model.fit <- try(coxph(model.formula, data=data), silent=TRUE)
    }
    if (type=="c") {
        model.text <- paste(yvar, "~", trtvar,"+",trtvar,"*id", sep="")
        model.formula <- as.formula(model.text)
        model.fit <- try(lm(model.formula, data=data), silent=TRUE)
    }
    if (type=="b") {
        model.text <- paste(yvar, "~", trtvar,"+",trtvar,"*id", sep="")
        model.formula <- as.formula(model.text)
        model.fit <- try(glm(model.formula, family=binomial,data=data), silent=TRUE)
    }
    # if (type=="rb") {
    #   model.text <- paste(yvar, "~", trtvar,"+",trtvar,"*id", sep="")
    #   model.formula <- as.formula(model.text)
    #   zelig.data.temp<<-data
    #   model.fit <- try(zelig(model.formula,model="relogit",data=zelig.data.temp,cite=FALSE), silent=TRUE)
    #   rm(zelig.data.temp, envir=.GlobalEnv)
    # }
    
    if (sigp.prcnt>min.sigp.prcnt & length(n.grp)>=2) {
        if (class(model.fit)!="try-error") {
            if (type=="c" | type=="b") {
                score <- try(summary(model.fit)$coefficient[4, 4],silent=TRUE)
                if(class(score)=="try-error") score=NA
            } else {
                score <- try(summary(model.fit)$coefficient[3,5],silent=TRUE)
                if (class(score)=="try-error") score=NA
            }
        } else {
            score <- NA
        }
    }
    else {
        score <- NA
    }
    score
}

#' Find cutoff for predictive case.
#'
#' @param data input data frame.
#' @param yvar response variable name.
#' @param censorvar censoring variable name.
#' @param xvar name of predictor for which cutpoint needs to be obtained.
#' @param trtvar treatment variable name.
#' @param type "c" continuous; "s" survival; "b" binary.
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1).
#' @param dir direction of cut.
#' @param nsubj number of subjects.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#'
#' @return the optimal score (p-value of subgroup*treatment interaction) for a predictor variable.
#' 
seqlr.find.cutoff.pred <- function(data, yvar, censorvar, xvar, trtvar, type, class.wt, dir, nsubj, min.sigp.prcnt) {
    # 
    cut.vec <- sort(unique(data[, xvar]))
    cut.vec <- quantile(cut.vec,  prob=seq(0.05, 0.95, 0.05),type=3)
    cut.vec=sort(unique(cut.vec))
    cut.score <- unlist(lapply(cut.vec, function(x) seqlr.score.pred(data, yvar, censorvar, xvar, trtvar, x, type, class.wt, dir, nsubj, min.sigp.prcnt)))
    cut.val <- cut.vec[which.min(cut.score)]
    cut.val
}

#' Main predictive BATTing function
#'
#' @param dataset input dataset in data frame
#' @param ids training indices
#' @param yvar response variable name
#' @param censorvar censoring variable name 1:event; 0: censor.
#' @param trtvar treatment variable name
#' @param type "c" continuous; "s" survival; "b" binary
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1)
#' @param xvar name of predictor for which cutpoint needs to be obtained
#' @param n.boot number of bootstraps for BATTing step.
#' @param des.res the desired response. "larger": prefer larger response. "smaller": prefer smaller response.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#'
#' @return a signature rule consisting of variable name, direction, optimal cutpoint and the corresponding p-value.
#' 
batting.pred <- function(dataset, ids, yvar, censorvar, trtvar, type, class.wt, xvar, n.boot, des.res, min.sigp.prcnt) {

    
    data <- dataset[ids,]
    niter <- n.boot
    nsubj<-nrow(data)

    if (type=="s") {
        model.text <- paste("Surv(", yvar, ", ", censorvar, ")~",trtvar,"+",trtvar, "*", xvar, sep="")
        model.formula <- as.formula(model.text)
        model.temp <- coxph(model.formula, data=data)
        coef.inter <- try(summary(model.temp)$coefficient[3,1],silent=TRUE)
    }
    if (type=="c") {
        model.text <- paste(yvar, "~", trtvar,"+",trtvar, "*", xvar, sep="")
        model.formula <- as.formula(model.text)
        model.temp <- lm(model.formula, data=data)
        coef.inter <- try(summary(model.temp)$coefficient[4,1],silent=TRUE)
    }
    if (type=="b") {
        model.text <- paste(yvar, "~", trtvar,"+",trtvar, "*", xvar, sep="")
        model.formula <- as.formula(model.text)
        model.temp <- glm(model.formula, family=binomial,data=data)
        coef.inter <- try(summary(model.temp)$coefficient[4,1],silent=TRUE)
    }
    # if (type=="rb") {
    #   model.text <- paste(yvar, "~", trtvar,"+",trtvar, "*", xvar, sep="")
    #   model.formula <- as.formula(model.text)
    #   zelig.data.temp<<-data
    #   model.temp <- zelig(model.formula,model="relogit",data=zelig.data.temp,cite=FALSE)
    #   rm(zelig.data.temp, envir=.GlobalEnv)
    #   coef.inter <- try(summary(model.temp)$coefficient[4,1],silent=TRUE)
    # }
    
    if (class(coef.inter)=="try-error") {
      dir <- NA
      model.pval <- NA
      cutoff.med <- NA
      cut.result <- c(xvar, dir, cutoff.med, model.pval)
      return(cut.result)   
    }
    
    
    if(des.res=="smaller"){
        if (type=="s") dir=ifelse(coef.inter>=0,">","<")
        if (type=="c") dir=ifelse(coef.inter>=0,"<",">")
        if (type=="b") dir=ifelse(coef.inter>=0,"<",">") 
    }
    
    if(des.res=="larger"){
        if (type=="s") dir=ifelse(coef.inter>=0,"<",">")
        if (type=="c") dir=ifelse(coef.inter>=0,">","<")
        if (type=="b") dir=ifelse(coef.inter>=0,">","<") 
    }
    
    
    cutoff.vec <- NULL
    
    if (!is.na(dir)) {
        
        for(i in 1:niter) {
            train.id <- sample(1:nsubj, nsubj, replace=TRUE)
            data.train <- data[train.id,]
            cutoff.temp <- seqlr.find.cutoff.pred(data.train, yvar, censorvar, xvar, trtvar, type, class.wt, dir, nsubj, min.sigp.prcnt)
            cutoff.vec <- c(cutoff.vec, cutoff.temp)
        }
        cutoff.med <- median(cutoff.vec, na.rm=TRUE)
        
        if (dir==">") {
            id <- as.numeric(data[, xvar]>cutoff.med)
        } else {
            id <- as.numeric(data[, xvar]<cutoff.med)
        }
        
        #n.id.trt <- table(data[id, trtvar])
        
        if (TRUE) {
            
            if (type=="s") {
                model.text1 <- paste("Surv(", yvar, ", ", censorvar, ")~", trtvar,"+",trtvar, "*", "id", sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(coxph(model.formula1, data=data), silent=TRUE)
            }
            
            if (type=="c") {
                model.text1 <- paste(yvar, "~", trtvar,"+",trtvar, "*", "id", sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(glm(model.formula1, data=data), silent=TRUE)
            }
            
            if (type=="b") {
                model.text1 <- paste(yvar, "~", trtvar,"+",trtvar, "*", "id", sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(glm(model.formula1, family=binomial, data=data), silent=TRUE)
            }
            
            # if (type=="rb") {
            #   model.text1 <- paste(yvar, "~", trtvar,"+",trtvar, "*", "id", sep="")
            #   model.formula1 <- as.formula(model.text1)
            #   data$id=id
            #   zelig.data.temp<<-data
            #   model.temp1 <- try(zelig(model.formula1, model="relogit",data=zelig.data.temp,cite=FALSE), silent=TRUE)
            #   rm(zelig.data.temp, envir=.GlobalEnv)
            #   data$id=NULL
            # }
            
            if (class(model.temp1)!="try-error") {
                
                if (type=="c" | type=="b") {
                    model.pval <- summary(model.temp1)$coefficient[4,4]
                } else {
                    model.pval <- summary(model.temp1)$coefficient[3,5]
                }
            } else {
                model.pval <- NA
            }
        } else {
            model.pval <- NA
        }
    } else {
        model.pval <- NA
        cutoff.med <- NA
    }
    
    cut.result <- c(xvar, dir, cutoff.med, model.pval)
    cut.result
}

###########################################################################  
#
#                  Prognostic BATTing Functions 
#
###########################################################################

#' Compute score of cutoff for prognostic case
#'
#' @param data input data frame.
#' @param yvar response variable name.
#' @param censorvar censoring variable name.
#' @param xvar name of predictor for which cutpoint needs to be obtained.
#' @param cutoff a specific cutpoint for which the score needs to be computed.
#' @param type "c" continuous; "s" survival; "b" binary.
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1).
#' @param dir direction of cut.
#' @param nsubj number of subjects.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#'
#' @return score (p-value of main effect) for the given cutoff.
#' 
seqlr.score.prog <- function(data, yvar, censorvar, xvar, cutoff, type, class.wt, dir, nsubj, min.sigp.prcnt) {
    # Compute score of cutoff for prognostic case.
    # score is the p-value of model.fit
    
    if (dir==">") {
        id <- data[, xvar]>cutoff
    }
    else {
        id <- data[, xvar]<cutoff
    }
    n.id <- table(id)
    data=cbind(data,id)
    sigp.prcnt<-sum(id)/nsubj

    if (type=="s") {
        model.text <- paste("Surv(", yvar, ", ", censorvar, ")~id", sep="")
        model.formula <- as.formula(model.text)
        model.fit <- try(coxph(model.formula, data=data), silent=TRUE)
    }
    if (type=="c") {
        model.text <- paste(yvar, "id", sep="~")
        model.formula <- as.formula(model.text)
        model.fit <- try(lm(model.formula, data=data), silent=TRUE)
    }
    if (type=="b") {

        model.text <- paste(yvar, "id", sep="~")
        model.formula <- as.formula(model.text)
        model.fit <- try(glm(model.formula, data=data, family="binomial"), silent=TRUE)
       
    }
    
    # if (type=="rb") {
    #   model.text <- paste(yvar, "id", sep="~")
    #   model.formula <- as.formula(model.text)
    #   zelig.data.temp<<-data
    #   model.fit <- try(zelig(model.formula, model="relogit",data=zelig.data.temp,cite=FALSE), silent=TRUE)
    #   rm(zelig.data.temp, envir=.GlobalEnv)
    # 
    # }
    
    if (sigp.prcnt>min.sigp.prcnt & length(n.id)>=2) {
        if (class(model.fit)!="try-error") {
            if (type=="c"|type=="b") {
                score <- try(summary(model.fit)$coefficient[2,4],silent=TRUE)
                if(class(score)=="try-error") score=NA
                
            }else {
                score <- try(summary(model.fit)$coefficient[5],silent=TRUE)
                if(class(score)=="try-error") score=NA
            }
        } else {
            score <- NA
        }
    }
    else {
        score <- NA
    }
    score
}


#' Find cutoff for prognostic case.
#'
#' @param data input data frame.
#' @param yvar response variable name.
#' @param censorvar censoring variable name.
#' @param xvar name of predictor for which cutpoint needs to be obtained.
#' @param type "c" continuous; "s" survival; "b" binary.
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1).
#' @param dir direction of cut.
#' @param nsubj number of subjects.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#'
#' @return the optimal score (p-value of main effect) for a predictor variable.
#' 
seqlr.find.cutoff.prog <- function(data, yvar, censorvar, xvar, type, class.wt, dir, nsubj, min.sigp.prcnt) {
    # Find cutoff for prognostic case.
    cut.vec <- sort(unique(data[, xvar]))
    cut.vec <- quantile(cut.vec, prob=seq(0.05, 0.95, 0.05),type=3)
    cut.vec=sort(unique(cut.vec))
    cut.score <- unlist(lapply(cut.vec, function(x) seqlr.score.prog(data, yvar, censorvar, xvar, x, type, class.wt, dir, nsubj, min.sigp.prcnt)))
    cut.val <- cut.vec[which.min(cut.score)]
    cut.val
}


#' Main prognostic BATTing function
#'
#' @param dataset input dataset in data frame
#' @param ids training indices
#' @param yvar response variable name
#' @param censorvar censoring variable name 1:event; 0: censor.
#' @param type "c" continuous; "s" survival; "b" binary
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1)
#' @param xvar name of predictor for which cutpoint needs to be obtained
#' @param n.boot number of bootstraps for BATTing step.
#' @param des.res the desired response. "larger": prefer larger response. "smaller": prefer smaller response.
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#'
#' @return a signature rule consisting of variable name, direction, optimal cutpoint and the corresponding p-value.
#' 
batting.prog <- function(dataset, ids, yvar, censorvar, type, class.wt, xvar, n.boot, des.res,min.sigp.prcnt) {
    # Main prognostic BATTing function
    data <- dataset[ids, ]
    niter <- n.boot
    nsubj<-nrow(data)

    if (type=="s") {
        model.text <- paste("Surv(", yvar, ", ", censorvar, ")~", xvar, sep="")
        model.formula <- as.formula(model.text)
        model.temp <- coxph(model.formula, data=data)
        coef.main <- try(summary(model.temp)$coefficient[1],silent=TRUE)    ##  Postive Hazard means greater risk
    }
    
    if (type=="c") {
        model.text <- paste(yvar, "~", xvar, sep="")
        model.formula <- as.formula(model.text)
        model.temp <- lm(model.formula, data=data)
        coef.main <- try(summary(model.temp)$coefficient[2,1],silent=TRUE)
    }
    
    if (type=="b") {
        model.text <- paste(yvar, "~", xvar, sep="")
        model.formula <- as.formula(model.text)
        model.temp <- glm(model.formula, family=binomial,data=data)
        coef.main <- try(summary(model.temp)$coefficient[2,1],silent=TRUE)
    }
    
    # if (type=="rb") {
    #   model.text <- paste(yvar, "~", xvar, sep="")
    #   model.formula <- as.formula(model.text)
    #   zelig.data.temp<<-data
    #   model.temp <- zelig(model.formula,model="relogit",data=zelig.data.temp,cite=FALSE)
    #   rm(zelig.data.temp, envir=.GlobalEnv)
    #   coef.main <- try(summary(model.temp)$coefficient[2,1],silent=TRUE)
    # }
    
    
    if (class(coef.main)=="try-error") {
      dir <- NA
      model.score <- NA
      cutoff.med <- NA
      cut.result <- c(xvar, dir, cutoff.med, model.score)
      return(cut.result)
    }
    
    
    if(des.res=="smaller"){
        if (type=="s") dir=ifelse(coef.main>=0,">","<")
        if (type=="c") dir=ifelse(coef.main>=0,"<",">")
        if (type=="b") dir=ifelse(coef.main>=0,"<",">") 
    }
    
    if(des.res=="larger"){
        if (type=="s") dir=ifelse(coef.main>=0,"<",">")
        if (type=="c") dir=ifelse(coef.main>=0,">","<")
        if (type=="b") dir=ifelse(coef.main>=0,">","<") 
    }
    
    
    cutoff.vec <- NULL
    
    if (!is.na(dir)) {

        for(i in 1:niter) {
            train.id <- sample(1:nsubj, nsubj, replace=TRUE)
            data.train <- data[train.id, ]
            cutoff.temp <- seqlr.find.cutoff.prog(data.train, yvar, censorvar, xvar, type, class.wt, dir, nsubj, min.sigp.prcnt)
            cutoff.vec <- c(cutoff.vec, cutoff.temp)
       
        }
        cutoff.med <- median(cutoff.vec, na.rm=TRUE)
        
        if (dir==">") {
            id <- data[, xvar]>cutoff.med
        } else {
            id <- data[, xvar]<cutoff.med
        }
        
        if (type=="s") {
            model.text1 <- paste("Surv(", yvar, ", ", censorvar, ")~id", sep="")
            model.formula1 <- as.formula(model.text1)
            model.temp1 <- try(coxph(model.formula1, data=data), silent=TRUE)
        }
        
        if (type=="b") {
            model.text1 <- paste(yvar, "id", sep="~")
            model.formula1 <- as.formula(model.text1)
            model.temp1 <- try(glm(model.formula1, family=binomial,data=data), silent=TRUE)
        }
        
        # if (type=="rb") {
        #   model.text1 <- paste(yvar, "id", sep="~")
        #   model.formula1 <- as.formula(model.text1)
        #   data$id=id
        #   zelig.data.temp<<-data
        #   model.temp1 <- try(zelig(model.formula1,model="relogit",data=zelig.data.temp,cite=FALSE), silent=TRUE)
        #   rm(zelig.data.temp, envir=.GlobalEnv)
        #   data$id=NULL
        # }
        
        if (type=="c") {
            model.text1 <- paste(yvar, "id", sep="~")
            model.formula1 <- as.formula(model.text1)
            model.temp1 <- try(lm(model.formula1, data=data), silent=TRUE)
        }
        
        if (class(model.temp1)!="try-error") {
            
            if (type=="c"|type=="b") {
                model.score <- summary(model.temp1)$coefficient[2, 4]
            } else {
                model.score <- summary(model.temp1)$coefficient[5]
            }
        } else {
            model.score <- NA
        }
    } else {
        model.score <- NA
        cutoff.med <- NA
    }
    
    cut.result <- c(xvar, dir, cutoff.med, model.score)
    cut.result
}

###############################################################################

#'  internal function used in seqlr.batting
#'
#' @param data the given dataset
#' @param rule rule is a vector of the form [x-variable, direction, cutoff, p-value]
#'
#' @return a logical variable indicating whether rules are satisfied or not. 
#' 
query.data <- function(data, rule) {
    n.rules <- nrow(rule)
    if (is.null(n.rules)) {
        var <- rule[1]
        dir <- rule[2]
        cutoff <- as.numeric(rule[3])
        if (dir=="<") {
            id <- data[, var]<cutoff
        } else {
            id <- data[, var]>cutoff
        }
        final.id <- id
    }
    else {
        rule.ids <- NULL
        for(i in 1:n.rules) {
            var <- rule[i,1]
            dir <- rule[i,2]
            cutoff <- as.numeric(rule[i,3])
            if (dir=="<") {
                id <- data[, var]<cutoff
            }
            else {
                id <- data[, var]>cutoff
            }
            rule.ids <- cbind(rule.ids, id)
        }
        final.id <- apply(as.matrix(rule.ids), 1, function(x) all(x))
        
    }
    final.id[is.na(final.id)] <- FALSE
    final.id
}


#' Perform sequential BATTing method.
#'
#' @param y data frame containing the response.
#' @param x data frame containing the predictors.
#' @param censor.vec vector containing the censor status (only for TTE data , censor=0,event=1)  - default = NULL.
#' @param trt.vec vector containing values of treatment variable ( for predictive signature). Set trt.vec to NULL for prognostic signature. 
#' @param trtref code for treatment arm.
#' @param type data type. "c" - continuous , "b" - binary, "s" - time to event : default = "c".
#' @param n.boot number of bootstraps in BATTing step.
#' @param des.res the desired response. "larger": prefer larger response. "smaller": prefer smaller response
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1)
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case.
#'
#' @return it returns a list of signature rules consisting of variable names, directions, thresholds and the loglikelihood at each step the  signatures are applied.
#' 
seqlr.batting <- function(y, x, censor.vec=NULL, trt.vec=NULL, trtref=NULL, type="c", n.boot=50,
                          des.res="larger", class.wt=c(1,1), min.sigp.prcnt=0.2, pre.filter=NULL, filter.method=NULL) {
  
    ####################  Main seqlr.BATTing Code ###############################
    # Code to reconfigure the +ve and -ve classes
    #
    ###########################################################################
    
    
    # Convert input arguments to data frames and check consistency.
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
        trt.vec <- as.data.frame(trt.vec)
        trtvar <- names(trt.vec)
        
        if (nrow(trt.vec) != nrow(x)) {
            stop("Error: trt.vec and x have different numbers of rows.")
        }
        Trt <- trt.vec[, trtvar]
        if (!is.null(trtref)) {
            # If trtref is provided, create 0-1 vector for trt1 in which 1 corresponds to trtref.
            trt1 <- rep(0,length=length(Trt))
            trt1[Trt==trtref] <- 1
            trt.vec[, trtvar] <- trt1
        } 
        
        data <- cbind(data, trt.vec)
        trtref <- NULL 
    } else {
        trtvar <- NULL
    }
    
    if (!is.null(pre.filter)){
      xvars_new=filter(data=data,type=type,yvar=yvar,xvars=xvars,censorvar=censorvar,trtvar=trtvar,n.boot=50,cv.iter=15,pre.filter=pre.filter,filter.method=filter.method)    
      del_ind=which(!(xvars %in% xvars_new))+1
      data=data[,-1*del_ind]
      x=x[,xvars_new]
      xvars_copy=xvars
      xvars=xvars_new
      marker.data=as.matrix(x)   
    }
    
    
    
    
    ll2beat <- -Inf
    dataset <- data
    data.id <- rep(TRUE, nrow(dataset))
    continue <- TRUE
    rule.string <- NULL
    var.names <- xvars
    
    while(continue) {
        if (!is.null(trt.vec)) {
            multivar.search <- lapply(var.names,  function(x) batting.pred(dataset, data.id, yvar, censorvar, trtvar, type, class.wt, x, n.boot, des.res, min.sigp.prcnt))
            mvar.result <- do.call(rbind, multivar.search)
            mvar.pvals <- as.numeric(mvar.result[, 4])
            threshold.rule.temp <- mvar.result[which.min(mvar.pvals), ]
            rule.string.temp <- rbind(rule.string, threshold.rule.temp)
            xvars.temp=rule.string.temp[,1]
            xvars.temp.fm=paste(trtvar,"*",xvars.temp,collapse="+")
     
            if (type=="s") {
                model.text1 <- paste("Surv(", yvar, ", ", censorvar, ")~", trtvar,"+",xvars.temp.fm, sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(coxph(model.formula1, data=dataset), silent=TRUE)
            }
            if (type=="c") {
                model.text1 <- paste(yvar, "~", trtvar,"+",xvars.temp.fm, sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(glm(model.formula1, data=dataset), silent=TRUE)
            }
            if (type=="b") {
                model.text1 <- paste(yvar, "~", trtvar,"+",xvars.temp.fm, sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(glm(model.formula1, family=binomial, data=dataset), silent=TRUE)
            }
            
            # if (type=="rb") {
            #   model.text1 <- paste(yvar, "~", trtvar,"+",xvars.temp.fm, sep="")
            #   model.formula1 <- as.formula(model.text1)
            #   dataset$data.id.temp=data.id.temp
            #   zelig.data.temp<<-dataset
            #   model.temp1 <- try(zelig(model.formula1, model="relogit",data=zelig.data.temp,cite=FALSE), silent=TRUE)
            #   rm(zelig.data.temp, envir=.GlobalEnv)
            #   dataset$data.id.temp=NULL
            # }
            
            if (class(model.temp1)!="try-error") {
                if (type=="c" | type=="b") {
                    min.ll <- logLik(model.temp1)[1]
                } else {
                    min.ll <- model.temp1$loglik[2]
                }
            } else {
                min.ll <- NA
            }
            
            
        } else {
            multivar.search <- lapply(var.names, function(x) batting.prog(dataset, data.id, yvar, censorvar, type, class.wt, x, n.boot, des.res, min.sigp.prcnt))
            mvar.result <- do.call(rbind, multivar.search)
            mvar.pvals <- as.numeric(mvar.result[, 4])
            threshold.rule.temp <- mvar.result[which.min(mvar.pvals), ]
            rule.string.temp <- rbind(rule.string, threshold.rule.temp)
            xvars.temp=rule.string.temp[,1]
            xvars.temp.fm=paste(xvars.temp,collapse="+")
            
            
            if (type=="s") {
                model.text1 <- paste("Surv(", yvar, ", ", censorvar, ")~",xvars.temp.fm, sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(coxph(model.formula1, data=dataset), silent=TRUE)
            }
            if (type=="c") {
                model.text1 <- paste(yvar, "~",xvars.temp.fm, sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(glm(model.formula1, data=dataset), silent=TRUE)
            }
            if (type=="b") {
                model.text1 <- paste(yvar, "~",xvars.temp.fm,  sep="")
                model.formula1 <- as.formula(model.text1)
                model.temp1 <- try(glm(model.formula1, family=binomial, data=dataset), silent=TRUE)
            }
            
            # if (type=="rb") {
            #   model.text1 <- paste(yvar, "~",xvars.temp.fm, sep="")
            #   model.formula1 <- as.formula(model.text1)
            #   dataset$data.id.temp=data.id.temp
            #   zelig.data.temp<<-dataset
            #   model.temp1 <- try(zelig(model.formula1, model="relogit",data=zelig.data.temp,cite=FALSE), silent=TRUE)
            #   rm(zelig.data.temp, envir=.GlobalEnv)
            #   dataset$data.id.temp=NULL
            # }
            
            if (class(model.temp1)!="try-error") {
                if (type=="c" | type=="b") {
                    min.ll <- logLik(model.temp1)[1]
                } else {
                    min.ll <- model.temp1$loglik[2]
                }
            } else {
                min.ll <- NA
            }                  
        }
        
        if (!is.na(min.ll)){
            obs.chisq=2*(min.ll-ll2beat)
            if (is.null(trt.vec)) pval.chisq <- ifelse(obs.chisq>0, 1-pchisq(obs.chisq, df=1), 1)
            if (!is.null(trt.vec)) pval.chisq <- ifelse(obs.chisq>0, 1-pchisq(obs.chisq, df=2), 1)
        }else{
            pval.chisq=1
        }
        
      data.id.temp<-query.data(dataset, rule.string.temp)
      data.prcnt<-sum(data.id.temp)/nrow(dataset) 
      
      
              
        if (pval.chisq<=0.05 & data.prcnt>min.sigp.prcnt) {
            threshold.rule <- threshold.rule.temp
            threshold.rule[4]=min.ll
            var.names <- var.names[var.names!=threshold.rule[1]]
            rule.string <- rbind(rule.string, threshold.rule)
            ll2beat <- min.ll
            data.id <- query.data(dataset, rule.string)
#             data.size <- sum(data.id)
#             continue <- ifelse(data.size>50, TRUE, FALSE)
            continue<-TRUE
        }else {
            continue <- FALSE
        }
    }
    
    rule.full <- rule.string
    
    
    if (!is.null(rule.full)) {
        colnames(rule.full) <- c("variable", "direction", "threshold", "LogLik")
        rule.full[, "threshold"] <- as.character(round(as.numeric(rule.full[, "threshold"]), 5))
    }
    
    
    return(rule.full)
}


#' Cross Validation for Sequential BATTing
#'
#' @param y data frame containing the response
#' @param x data frame containing the predictors
#' @param censor.vec vector giving the censor status (only for TTE data , censor=0,event=1) : default = NULL
#' @param trt.vec vector containing values of treatment variable ( for predictive signature). Set trt.vec to NULL for prognostic signature. 
#' @param trtref code for treatment arm.
#' @param type data type. "c" - continuous , "b" - binary, "s" - time to event : default = "c".
#' @param n.boot number of bootstraps in BATTing step.
#' @param des.res the desired response. "larger": prefer larger response. "smaller": prefer smaller response
#' @param class.wt vector of length 2 used to weight the accuracy score , useful when there is class imbalance in binary data defaults to c(1,1)
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff. 
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected.
#' @param filter.method NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case.
#' @param k.fold number of folds for CV.
#' @param cv.iter algorithm terminates after cv.iter successful iterations of cross-validation.
#' @param max.iter total number of iterations allowed (including unsuccessful ones).
#'
#' @return a list containing with following entries: 
#' \item{stats.summary}{Summary of performance statistics.}
#' \item{pred.classes}{Data frame containing the predictive clases (TRUE/FALSE) for each iteration.}
#' \item{folds}{Data frame containing the fold indices (index of the fold for each row) for each iteration.}
#' \item{sig.list}{List of length cv.iter * k.fold containing the signature generated at each of the  k folds, for all iterations.}
#' \item{error.log}{List of any error messages that are returned at an iteration.}
#' \item{interplot}{Treatment*subgroup interaction plot for predictive case}
#' 
cv.seqlr.batting <- function(y, x, censor.vec=NULL, trt.vec=NULL, trtref=NULL, type="c", n.boot=50, 
                             des.res="larger", class.wt=c(1, 1), min.sigp.prcnt=0.2, 
                             pre.filter=NULL, filter.method=NULL, k.fold=5, cv.iter=50, max.iter=500) {

    # Convert input arguments to data frames and check consistency.
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
        trt.vec <- as.data.frame(trt.vec)
        trtvar <- names(trt.vec)
        
        if (nrow(trt.vec) != nrow(x)) {
            stop("Error: trt.vec and x have different numbers of rows.")
        }
        Trt <- trt.vec[, trtvar]
        if (!is.null(trtref)) {
            # If trtref is provided, create 0-1 vector for trt1 in which 1 corresponds to trtref.
            trt1 <- rep(0,length=length(Trt))
            trt1[Trt==trtref] <- 1
            trt.vec[, trtvar] <- trt1
        } 
        
        data <- cbind(data, trt.vec)
        trtref <- NULL 
    } else {
        trtvar <- NULL
    }
    
    if (type=="b") strata=data[,yvar]
    if (type=="s") strata=censor.vec[, censorvar]
    if (type=="c") strata=NULL
    
    
    model.Rfunc <- "seqlr.batting.wrapper"
    # model.Rfunc.args <- make.arg.list("seq.batting")
    # Create a list containing arguments to seq.batting
    
    model.Rfunc.args <- list(yvar=yvar, xvars=xvars, censorvar=censorvar, trtvar=trtvar, trtref=trtref, 
                             type=type, n.boot=n.boot, des.res=des.res, class.wt=class.wt, 
                             min.sigp.prcnt=min.sigp.prcnt, pre.filter=pre.filter, filter.method=filter.method)
    # List of arguments for model.Rfunc. Must be packaged into 
    # a list to pass to kfold.cv
    
    predict.Rfunc <- "pred.seqlr.cv"
    # List of arguments for predict.Rfunc.
    predict.Rfunc.args <- list(yvar=yvar, xvars=xvars)
    res <- kfold.cv(data=data, model.Rfunc=model.Rfunc, model.Rfunc.args=model.Rfunc.args, predict.Rfunc=predict.Rfunc, predict.Rfunc.args=predict.Rfunc.args, k.fold=k.fold, cv.iter=cv.iter, strata=strata, max.iter=max.iter)
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



#' Wrapper function for seqlr.batting, to be passed to kfold.cv.
#'
#' @param data data frame equal to cbind(y, x, trt, censor), where y and x are inputs to seqlr.batting.
#' @param args list containing all other input arguments to seq.batting except for x and y. Also contains xvars=names(x) and yvar=names(y).
#'
#' @return  prediction rule returned by seqlr.batting.
#' 
seqlr.batting.wrapper <- function(data, args) {

  yvar <- args$yvar 
  xvars <- args$xvars
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
  
  trtref <- args$trtref
  type <- args$type
  n.boot <- args$n.boot
  des.res <- args$des.res
  class.wt <- args$class.wt
  min.sigp.prcnt <- args$min.sigp.prcnt
  pre.filter <- args$pre.filter
  filter.method = args$filter.method
  res <- seqlr.batting(y=y, x=x, censor.vec=censor.vec, trt.vec=trt.vec, trtref=trtref, type=type, n.boot=n.boot,
                     des.res=des.res, class.wt=class.wt, min.sigp.prcnt = min.sigp.prcnt, pre.filter=pre.filter, filter.method=filter.method)
  return(res)
}

#' Prediction function for CV Sequential BATTing
#' 
#' @description Assign positive and negative groups for cross-validation data given prediction rule in predict.rule.
#'
#' @param data input data frame
#' @param predict.rule Prediction rule returned by seqlr.batting.
#' @param args Prediction rule arguments
#'
#' @return a logical vector indicating the prediction for each row of data.
#' 
pred.seqlr.cv <- function(data, predict.rule, args) {

    n.rules <- nrow(predict.rule)
    if (is.null(n.rules)) {
        var <- predict.rule["variable"]
        dir <- predict.rule["direction"]
        cutoff <- as.numeric(predict.rule["threshold"])
        if (dir=="<") {
            id <- data[,var]<cutoff
        } else {
            id <- data[,var]>cutoff
        }
        pred.id <- id
    }
    else {
        rule.ids <- NULL
        for(i in 1:n.rules) {
            var <- predict.rule[i, "variable"]
            dir <- predict.rule[i, "direction"]
            cutoff <- as.numeric(predict.rule[i, "threshold"])
            if (dir=="<") {
                id <- data[,var]<cutoff
            }
            else {
                id <- data[,var]>cutoff
            }
            rule.ids <- cbind(rule.ids,id)
        }
        pred.id <- apply(as.matrix(rule.ids),1,function(x) all(x))
    }
    pred.id[is.na(pred.id)] <- FALSE
  
    pred.data <- cbind(data, data.frame(pred.class=pred.id))
    # Note!!: If you don't cbind pred.class with the data in the current fold, the 
    # row numbers of the current fold will be lost.
    pred.data <- subset(pred.data, select="pred.class")
    # Now that the row numbers are attached, you can get rid of marker.data.
    pred.data
}

#' Prediction function for Sequential BATTing
#' 
#' @description Assign positive and negative groups based on predict.rule, the output of seqlr.batting.
#'
#' @param x input predictors matrix
#' @param predict.rule Prediction rule returned by seqlr.batting.
#'
#' @return a logical vector indicating the prediction for each row of data.
#' 
pred.seqlr <- function(x, predict.rule) {
    n.rules <- nrow(predict.rule) 
    rule.ids <- NULL
    for (i in 1:n.rules) {
        var <- predict.rule[i, "variable"]
        dir <- predict.rule[i, "direction"]
        cutoff <- as.numeric(predict.rule[i, "threshold"])
        if (dir=="<") {
            id <- x[, var]<cutoff
        } else {
            id <- x[, var]>cutoff
        }
        rule.ids <- cbind(rule.ids,id)
    }
    pred.id <- apply(as.matrix(rule.ids),1,function(x) all(x))
    pred.id[is.na(pred.id)] <- FALSE
    pred.data<-data.frame(x, pred.class=pred.id)
    pred.data
}


###################################################################################################

#' Get signature variables from output of seqlr.batting.
#'
#' @param sig.list signature list returned by seqlr.batting.
#' @param xvars predictor variable names
#'
#' @return the variables included in signature rules returned by seqlr.batting
#' 
get.var.counts.seq <- function(sig.list, xvars) {
    sig.var.counts <- data.frame(xvars)
    sig.var.counts$count <- rep(0, length(xvars))
    sig.var.counts$greater.than <- rep(0, length(xvars))
    sig.var.counts$less.than <- rep(0, length(xvars))
    
    for (i in 1:length(sig.list)) {
        predict.rule <- sig.list[[i]]  
        n.rules <- nrow(predict.rule)
        if (is.null(n.rules)) {
            # Can this happen?
            var <- predict.rule["variable"]
            dir <- predict.rule["direction"]
            sig.var.counts[xvars==var, ]$count <- sig.var.counts[xvars==var, ]$count + 1
            if (dir==">") {
                sig.var.counts[xvars==var, ]$greater.than <- sig.var.counts[xvars==var, ]$greater.than + 1
            } else {
                sig.var.counts[xvars==var, ]$greater.than <- sig.var.counts[xvars==var, ]$greater.than + 1
            }
        } else {
            for(i in 1:n.rules) {
                var <- predict.rule[i, "variable"]
                dir <- predict.rule[i, "direction"]
                sig.var.counts[xvars==var, ]$count <- sig.var.counts[xvars==var, ]$count + 1
                if (dir==">") {
                    sig.var.counts[xvars==var, ]$greater.than <- sig.var.counts[xvars==var, ]$greater.than + 1
                } else {
                    sig.var.counts[xvars==var, ]$greater.than <- sig.var.counts[xvars==var, ]$greater.than + 1
                }
            }
        }    
    }
    sig.var.counts
}
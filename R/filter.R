#' Filter function for predictive/prognostic biomarker candidates for signature development
#' @description Filter function for Prognostic and preditive biomarker signature development for Exploratory Subgroup Identification in Randomized Clinical Trials 
#' 
#' @usage filter(data,
#' type="c",
#' yvar,
#' xvars,
#' censorvar=NULL,
#' trtvar=NULL,
#' trtref=1,
#' n.boot=50,
#' cv.iter=20,
#' pre.filter=length(xvars),
#' filter.method=NULL)
#' 
#' @param data input data frame 
#' @param type type of response variable: "c" continuous; "s" survival; "b" binary
#' @param yvar variable (column) name for response variable 
#' @param xvars vector of variable names for predictors (covariates)
#' @param censorvar variable name for censoring (1: event; 0: censor), default = NULL
#' @param trtvar variable name for treatment variable, default = NULL (prognostic signature)
#' @param trtref coding (in the column of trtvar) for treatment arm, default = 1 (no use for prognostic signature)
#' @param n.boot number of bootstrap for the BATTing procedure
#' @param cv.iter Algotithm terminates after cv.iter successful iterations of cross-validation, or after max.iter total iterations, whichever occurs first
#' @param pre.filter NULL (default), no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
 
#' @param filter.method NULL (default), no prefiltering; "univariate", univaraite filtering; "glmnet", glmnet filtering
#' 
#' @details The function contains two algorithms for filtering high-dimentional multivariate (prognostic/predictive) biomarker candidates via univariate fitering (used p-values of group difference for prognostic case, p-values of interaction term for predictive case); LASSO/Elastic Net method. (Tian L. et al 2012)
#'
#' @return \item{var}{a vector of filter results of variable names}
#' @references Tian L, Alizadeh A, Gentles A, Tibshirani R (2012) A Simple Method for Detecting Interactions between a Treatment and a Large Number of Covariates. J Am Stat Assoc. 2014 Oct; 109(508): 1517-1532.
#' @export 
#'
#' @examples 
#' \dontrun{
#' data(Sepsis.train)
#' 
#' yvar="survival"
#' xvars=names(Sepsis.train)[2:12]
#' trtvar="THERAPY"
#' trtref="active"
#' set.seed(123)
#' 
#' filter.res <- filter(data=Sepsis.train,
#' type="b",
#' yvar=yvar,
#' xvars=xvars,
#' trtvar=trtvar,
#' trtref=trtref,
#' pre.filter=20,
#' filter.method="univariate")
#' 
#' filter.res
#' 
#' set.seed(123)
#' filter.res <- filter(data=Sepsis.train,
#' type="b",
#' yvar=yvar,
#' xvars=xvars,
#' trtvar=trtvar,
#' trtref=trtref,
#' pre.filter="opt", 
#' filter.method="glmnet")
#' 
#' filter.res
#' }
#' 
#' 
#' @aliases filter
#' 
filter <- function(data,type="c",yvar,xvars,censorvar=NULL,trtvar=NULL,trtref=1,n.boot=50,cv.iter=20,pre.filter=length(xvars), 
                filter.method=NULL){
  if(filter.method == "glmnet"){
    filter.glmnet(data,type,yvar,xvars,censorvar,trtvar,trtref,n.boot,cv.iter,pre.filter); 
  } else if (filter.method == "univariate"){
    filter.univariate(data,type,yvar,xvars,censorvar,trtvar,trtref,pre.filter); 
  } else if (filter.method == "unicart") {
    filter.unicart(data,type,yvar,xvars,censorvar,trtvar,trtref,pre.filter)
  }else {
    warning("No pre-filtering was selected");
    return(xvars);
  }
}

#####################################################################
#' Flitering using MC glmnet
#'
#' @param data input data frame
#' @param type "c" continuous; "s" survival; "b" binary
#' @param yvar response variable name
#' @param xvars covariates variable name
#' @param censorvar censoring variable name 1:event; 0: censor.
#' @param trtvar treatment variable name 
#' @param trtref code for treatment arm
#' @param n.boot number of bootstrap for filtering 
#' @param cv.iter number of iterations required for MC glmnet filtering
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#'
#' @return variables selected after glmnet filtering
#' 
filter.glmnet <- function(data,type,yvar,xvars,censorvar,trtvar,trtref,n.boot=50,cv.iter=20,pre.filter=length(xvars)){
  pmax=pre.filter
  if (pre.filter=="opt") pmax=length(xvars)
    
#  library(glmnet)
  if (type=="s") family="cox"
  if (type=="c") family="gaussian"
  if(type=="b") family="binomial"
  
  var.sel.all=NULL

  if (!is.null(trtvar)){  ## predictive
    for (i in 1:n.boot){
      trt=1*(data[,trtvar]==trtref)
      boot_num=max(sum(trt),length(trt)-sum(trt))
      
      data_trt0=data[trt==0,]    
      boot_ind0=sample(1:dim(data_trt0)[1],size=boot_num,replace=TRUE)    
      boot_data0=data_trt0[boot_ind0,]        
      data_trt1=data[trt==1,]
      boot_ind1=sample(1:dim(data_trt1)[1],size=boot_num,replace=TRUE)
      boot_data1=data_trt1[boot_ind1,]    
      data_new=rbind(boot_data0,boot_data1)
      
      
      y=data_new[,yvar]
      x=cbind(rep(1,length(y)),scale(as.matrix(data_new[,xvars])))
      trt=1*(data_new[,trtvar]==trtref)
      trt[trt==0]=-1
      w_star=x*trt/2 
      if (type=="s") {
        censor=data_new[,censorvar]
        y=cbind(time=y,status=censor)
      }
      

      ## MCglmnet filtering
      lambda.opt=NULL
      for (j in 1:cv.iter){
        fit.cv=cv.glmnet(w_star,y,family=family,standardize=FALSE,intercept=FALSE,nfolds=5,pmax=pmax)
        lambda.opt=c(lambda.opt,fit.cv$lambda.1se)
      }
      
      lambda.unique=sort(unique(lambda.opt))
      lambda.tab=table(lambda.opt)
      lambda.med=lambda.unique[which.max(lambda.tab)]
      
      fit=glmnet(w_star,y,family=family,standardize=FALSE,intercept=FALSE,pmax=pmax)
      beta.opt=as.matrix(fit$beta)[,fit$lambda==lambda.med]
      var.sel=which(beta.opt!=0)
      var.sel.all=c(var.sel.all,var.sel)                
    }
  }else{  ## prognostic
    for (i in 1:n.boot){
      boot_num=dim(data)[1]
      boot_ind=sample(1:boot_num,size=boot_num,replace=TRUE)
      data_new=data[boot_ind,]
      
      y=data_new[,yvar]
      x=as.matrix(data_new[,xvars])
      if (type=="s") {
        censor=data_new[,censorvar]
        y=cbind(time=y,status=censor)
      }
      
      lambda.opt=NULL
      for (j in 1:cv.iter){
        fit.cv=cv.glmnet(x,y,family=family,nfolds=5,pmax=pmax)
        lambda.opt=c(lambda.opt,fit.cv$lambda.1se)
      }
      
      lambda.unique=sort(unique(lambda.opt))
      lambda.tab=table(lambda.opt)
      lambda.med=lambda.unique[which.max(lambda.tab)]
      
      fit=glmnet(x,y,family=family,pmax=pmax)
      beta.opt=as.matrix(fit$beta)[,fit$lambda==lambda.med]
      var.sel=which(beta.opt!=0)
      var.sel.all=c(var.sel.all,var.sel)
            
    }
  }
  

  
  var.tab=sort(table(var.sel.all),decreasing=TRUE)
  var.ind=as.numeric(names(var.tab)[1:min(pmax,length(var.tab))])
  if (!is.null(trtvar)) var.ind=var.ind-1
  var=xvars[var.ind]
  return(var)
}

#####################################################################
### 
#' Univariate Filtering
#'
#' @param data input data frame
#' @param type "c" continuous; "s" survival; "b" binary
#' @param yvar response variable name
#' @param xvars covariates variable name
#' @param censorvar censoring variable name 1:event; 0: censor.
#' @param trtvar treatment variable name 
#' @param trtref code for treatment arm
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#'
#' @return covariate names after univariate filtering.
#' 
filter.univariate <- function(data,type,yvar,xvars,censorvar,trtvar,trtref=1, pre.filter=length(xvars)){
  pmax=pre.filter;
  if (pmax=="opt") pmax=min(10, length(xvars))
  
  if (!is.null(trtvar)){  ## predictive
      if (type=="s"){
          res.all=NULL
          for (i in 1:length(xvars)){
              data.tmp = data.frame(time=data[,yvar], event=data[,censorvar], trt=(data[,trtvar]==trtref)*1, x=data[,xvars[i]])
              res.tmp = try(summary(coxph(Surv(time, event)~trt*x, data=data.tmp))$coefficients["trt:x","Pr(>|z|)"], silent=T)
              if(class(res.tmp)=="try-error") res.tmp=1
              res.all = c(res.all, res.tmp)
          }
      }
      
      
      if(type=="b"){
          res.all = NULL
          for (i in 1:length(xvars)){
              data.tmp = data.frame(y=data[,yvar], trt=(data[,trtvar]==trtref)*1, x=data[,xvars[i]])
              res.tmp = try(summary(glm(y~trt*x, family=binomial,data=data.tmp))$coefficients["trt:x","Pr(>|z|)"], silent=T)
              if(class(res.tmp)=="try-error") res.tmp=1
              res.all = c(res.all, res.tmp)
          }
      }
      
      if(type=="c" ){
          res.all = NULL
          for (i in 1:length(xvars)){
              data.tmp = data.frame(y=data[,yvar], trt=(data[,trtvar]==trtref)*1, x=data[,xvars[i]])
              res.tmp = try(summary(lm(y~trt*x, data=data.tmp))$coefficients["trt:x","Pr(>|t|)"], silent=T)
              if(class(res.tmp)=="try-error") res.tmp=1
              res.all = c(res.all, res.tmp)
          }
      }
        
  }else{
      
      if (type=="s"){
          res.all=NULL
          for (i in 1:length(xvars)){
              data.tmp = data.frame(time=data[,yvar], event=data[,censorvar], x=data[,xvars[i]])
              res.tmp = try(summary(coxph(Surv(time, event)~x, data=data.tmp))$coefficients["x","Pr(>|z|)"], silent=T)
              if(class(res.tmp)=="try-error") res.tmp=1
              res.all = c(res.all, res.tmp)
          }
      }
      
      
      if(type=="b"){
          res.all = NULL
          for (i in 1:length(xvars)){
              data.tmp = data.frame(y=data[,yvar], x=data[,xvars[i]])
              res.tmp = try(summary(glm(y~x, family=binomial,data=data.tmp))$coefficients["x","Pr(>|z|)"], silent=T)
              if(class(res.tmp)=="try-error") res.tmp=1
              res.all = c(res.all, res.tmp)
          }
      }
      
      if(type=="c" ){
          res.all = NULL
          for (i in 1:length(xvars)){
              data.tmp = data.frame(y=data[,yvar], x=data[,xvars[i]])
              res.tmp = try(summary(lm(y~x, data=data.tmp))$coefficients["x","Pr(>|t|)"], silent=T)
              if(class(res.tmp)=="try-error") res.tmp=1
              res.all = c(res.all, res.tmp)
          }
      }
  }
  
  if(length(res.all) != length(xvars)) stop("Not all univariate p values are properly recorded!")
  xvar.order = order(res.all)
  top.var = xvars[xvar.order][1:pmax]
  
  for (i in 1:length(xvars)){
      y.tmp = data[,yvar]
      x.tmp = data[,xvars[i]]
      table.xy = table(y.tmp, x.tmp)
      if(dim(table.xy)[1]==2 & dim(table.xy)[2]==2 & sum(table.xy==0)>0) top.var= c(top.var, xvars[i])
  }
  
  top.var = unique(top.var)
  top.var
}

#####################################################################
###### rpart filtering 
#' rpart filtering (only for prognostic case)
#'
#' @param data input data frame
#' @param type "c" continuous; "s" survival; "b" binary
#' @param yvar response variable name
#' @param xvars covariates variable name
#' @param censorvar censoring variable name 1:event; 0: censor.
#' @param trtvar treatment variable name 
#' @param trtref code for treatment arm
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#'
#' @return selected covariates after rpart filtering
#' 
filter.unicart = function (data,type,yvar,xvars,censorvar,trtvar, trtref=1, pre.filter=length(xvars)){
#    library(rpart)
    
    pmax=pre.filter
    if (pmax=="opt") pmax=min(10,length(xvars))
    
    if (!is.null(trtvar)) {
        warning("cart cannot do predictive case. No filtering performed.")
        return(xvars)
    }

    var.score=NULL
          
    for (xvar in xvars){
        if (type == "c"){
            formula.tx = as.formula(paste(yvar,"~",xvar))
            method = "anova"
        } 
        if (type == "b"){
            formula.tx = as.formula(paste(yvar,"~",xvar))
            method = "class"
        }
        if (type == "s"){
            formula.tx = as.formula(paste("Surv(",yvar,",",censorvar,")~",xvar))
            method = "exp"
        }
            
        fit = rpart(formula = formula.tx, data=data, method=method)
        if (!is.null(fit$variable.importance))
        var.score=rbind(var.score, data.frame(xvar=xvar,score=fit$variable.importance,stringsAsFactors=F))
    }
    xvar.order=order(var.score$score,decreasing=T)
    xvars.sel = var.score[xvar.order,"xvar"][1:pmax]
    xvars.sel
    
} 



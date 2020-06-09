##################### PREDICTIVE CASE ########################



########### prim train ###############

#' The main PRIM function
#'
#' @param data the input data frame
#' @param yvar name for response variable
#' @param censorvar name for censoring (1: event; 0: censor), default = NULL
#' @param trtvar name for treatment variable, default = NULL (prognostic signature)
#' @param trtref coding (in the column of trtvar) for treatment arm
#' @param xvars vector of variable names for predictors (covariates)
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response
#' @param alpha a parameter controlling the number of patients in consideration
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#' @param training.percent percentage of subjects in the initial training data
#' @param n.boot number of bootstrap for the variable selection procedure for PRIM
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case
#'
#' @return the final list of rules selected by PRIM.
prim.train <- function(
  # data info
  data, yvar, censorvar, trtvar, trtref=NULL, xvars, type, des.res,
  # control arguments
  alpha = c(0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50),
  min.sigp.prcnt = 0.2  ,  training.percent = 0.50, n.boot=0, pre.filter=NULL, filter.method=NULL) {

  data = as.data.frame(data)

  if (!is.null(trtref)){
    data[,trtvar]=(data[,trtvar]==trtref)+0
    trtref=NULL
  }

  if (!is.null(pre.filter) & type!="rb"){
    xvars_new=filter(data=data,type=type,yvar=yvar,xvars=xvars,censorvar=censorvar,trtvar=trtvar,n.boot=50,cv.iter=15,pre.filter=pre.filter,filter.method=filter.method)
    del_xvars=xvars[!(xvars %in% xvars_new)]
    data_copy=data
    data=data[,!colnames(data)%in%del_xvars]
    xvars_copy=xvars
    xvars=xvars_new
  }

  xvars.sel=xvars
  min.size.inside=floor(min.sigp.prcnt*dim(data)[1])

  if (n.boot>=1){
    n.total = dim(data)[1]
    n.sample = ceiling(n.total * 0.632)
    i=1
    j=1
    var.sel.count=matrix(0,nrow=1,ncol=length(xvars),dimnames=list(NULL,xvars))
    num.var.sel = NULL

    if (!is.null(trtvar))
    {
      while (i <= n.boot){
        data.sample = data[sample(1:n.total, size=n.sample, replace=FALSE),]
        p.collection = NULL
        sig.collection = NULL

        pidx.train.test = create.training.dataset.index(training.percent=training.percent, dim(data.sample)[1]);
        for (alpha.i in alpha){
          result.once = prim.train.pred.once(data=data.sample, yvar=yvar, censorvar=censorvar,
                                             trtvar=trtvar, xvars=xvars, type=type, des.res=des.res,
                                             alpha=alpha.i, min.size.inside=min.size.inside,
                                             pidx.train.test=pidx.train.test)

          if (length(result.once)==0) next
          p.collection = rbind(p.collection,result.once$test.result.ordered[1,])
          sig.collection = c(sig.collection,list(result.once$signature))
        }

        j=j+1
        if (is.null(p.collection) || is.null(sig.collection))  {
          if (j==3*n.boot) break
          next
        }

        sig.sel.once =sig.collection[[which.min(p.collection[,"pv1"])]]
        xvars.sel.once = sig.sel.once[,"x.nm"]
        var.sel.count[,xvars.sel.once] = var.sel.count[,xvars.sel.once] +1
        num.var.sel=c(num.var.sel,length(xvars.sel.once))
        i=i+1
      }

    }else
    {
      while (i <= n.boot){
        data.sample = data[sample(1:n.total, size=n.sample, replace=FALSE),]
        p.collection = NULL
        sig.collection = NULL

        pidx.train.test = create.training.dataset.index(training.percent=training.percent, dim(data.sample)[1]);
        for (alpha.i in alpha){
          result.once = prim.train.prog.once(data=data.sample, yvar=yvar, censorvar=censorvar,
                                             xvars=xvars, type=type, des.res=des.res,
                                             alpha=alpha.i, min.size.inside=min.size.inside,
                                             pidx.train.test=pidx.train.test)

          if (length(result.once)==0) next
          p.collection = rbind(p.collection,result.once$test.result.ordered[1,])
          sig.collection = c(sig.collection,list(result.once$signature))
        }

        j=j+1
        if (is.null(p.collection) || is.null(sig.collection))  {
          if (j==3*n.boot) break
          next
        }

        sig.sel.once =sig.collection[[which.min(p.collection[,"pv"])]]
        xvars.sel.once = sig.sel.once[,"x.nm"]
        var.sel.count[,xvars.sel.once] = var.sel.count[,xvars.sel.once] +1
        num.var.sel=c(num.var.sel,length(xvars.sel.once))
        i=i+1
      }

    }

    if (!is.null(num.var.sel) & sum(var.sel.count>0,na.rm=T))
      med.num.val.sel = min(ceiling(median(num.var.sel,na.rm=T)+mad(num.var.sel,na.rm=T)),sum(var.sel.count>0,na.rm=T))
    xvars.sel = names(sort(var.sel.count[1,]))[1:med.num.val.sel]
  }


  pidx.train.test = create.training.dataset.index(training.percent=training.percent, dim(data)[1])

  if(!is.null(trtvar))
  {
    if (!is.null(trtref)){
      data[,trtvar]=(data[,trtvar]==trtref)+0
      trtref=NULL
    }

    p.collection = NULL
    sig.collection = NULL
    for (alpha.i in alpha){
      result.once = prim.train.pred.once(data=data, yvar=yvar, censorvar=censorvar,
                                         trtvar=trtvar, xvars=xvars.sel, type=type, des.res=des.res,
                                         alpha=alpha.i, min.size.inside=min.size.inside,
                                         pidx.train.test=pidx.train.test)
      if (length(result.once)==0) next
      p.collection = rbind(p.collection,result.once$test.result.ordered[1,])
      sig.collection = c(sig.collection,list(result.once$signature))
    }



  if (is.null(p.collection) || is.null(sig.collection))  {
    return(NULL)
  }

  sig.sel.final =sig.collection[[which.min(p.collection[,"pv1"])]]

  sig.sel.final
  }
  else
  {
    p.collection = NULL
    sig.collection = NULL
    for (alpha.i in alpha){
      result.once = prim.train.prog.once(data=data, yvar=yvar, censorvar=censorvar,
                                         xvars=xvars.sel, type=type, des.res=des.res,
                                         alpha=alpha.i, min.size.inside=min.size.inside,
                                         pidx.train.test=pidx.train.test)
      if (length(result.once)==0) next
      p.collection = rbind(p.collection,result.once$test.result.ordered[1,])
      sig.collection = c(sig.collection,list(result.once$signature))
    }
    if (is.null(p.collection) || is.null(sig.collection))  {
      return(NULL)
    }

    sig.sel.final =sig.collection[[which.min(p.collection[,"pv"])]]

    sig.sel.final

  }

  }





#################### TRAIN ONCE ###############
#' Apply PRIM one time on the training data for a fixed value of alpha in predictive case
#' @description this function applies the prim procedure (peeling, pasting, and dropping operations) on training data one time.
#'
#' @param data input data frame
#' @param yvar name for response variable
#' @param censorvar name for censoring (1: event; 0: censor), default = NULL
#' @param trtvar 0-1 coded vector for treatment variable
#' @param xvars vector of variable names for predictors (covariates)
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response
#' @param alpha a parameter controlling the number of patients in consideration
#' @param min.size.inside desired number of subjects in signature positive group size for a given cutoff.
#' @param pidx.train.test training and test data index as obtained from create.training.dataset.index
#'
#' @return a list containing signature rules and test result based on the signatures.
#'
prim.train.pred.once <- function(
  # data info
  data, yvar, censorvar, trtvar, xvars,type, des.res,
  # control arguments
  alpha = 0.10,  min.size.inside = 20, pidx.train.test) {

  ## method constants
  g.str = ">=";
  l.str = "<=";


  ## create training and testing data sets
  data.train = data[pidx.train.test$pidx.train, ]; # training data
  data.test = data[pidx.train.test$pidx.test, ]; # testing data


  d.inside = data.train;
  d.outside = data.frame();


  # data structure to save results
  trace.peeling.inside.condition = data.frame(); # the most recent updated box condition
  trace.pasting.inside.condition = data.frame()
  trace.dropping.inside.condition = data.frame()

  pos.group.peeling.list = list();
  idx.pos.group.peeling.list = 1;
  pos.group.pasting.list = list();
  idx.pos.group.pasting.list = 1;
  pos.group.dropping.list = list();
  idx.pos.group.dropping.list = 1;

  # peeling
  repeat {
    one.peeling.results.ordered = one.peeling(d.inside, d.outside, xvars, alpha, min.size.inside, yvar, censorvar, trtvar, g.str, l.str, type, des.res);
    if(dim(one.peeling.results.ordered)[1] == 0) break
    new.inside.condition = one.peeling.results.ordered[1,];
    trace.peeling.inside.condition = combine.condition(trace.peeling.inside.condition, new.inside.condition)

    x.nm = new.inside.condition[1, "x.nm"];
    x.cutoff = new.inside.condition[1, "cutoff"]

    if(new.inside.condition[1, "condition"] == g.str) {
      bidx.obs.inside.to.peel = d.inside[, x.nm] >= x.cutoff;
    } else if((new.inside.condition[1, "condition"] == l.str)) {
      bidx.obs.inside.to.peel = d.inside[, x.nm] <= x.cutoff;
    }

    # NOTE: need to update outside data first
    d.outside = rbind(d.outside, d.inside[!bidx.obs.inside.to.peel, ]);
    d.inside = d.inside[bidx.obs.inside.to.peel, ];

    pos.group.peeling.list[idx.pos.group.peeling.list]=list(trace.peeling.inside.condition)
    idx.pos.group.peeling.list=idx.pos.group.peeling.list+1
  }

  if (dim(trace.peeling.inside.condition)[1]==0) return(list())

  #pasting
  trace.pasting.inside.condition = trace.peeling.inside.condition
  repeat {
    one.pasting.results.ordered = one.pasting(d.inside, d.outside, trace.pasting.inside.condition, alpha, yvar, censorvar, trtvar, g.str, l.str, type, des.res);
    if(dim(one.pasting.results.ordered)[1] == 0) break

    new.pasting.inside.condition = one.pasting.results.ordered[1,];
    temp.trace.pasting.inside.condition = combine.condition(trace.pasting.inside.condition, new.pasting.inside.condition)
    if (dim(temp.trace.pasting.inside.condition)[1]==0) break

    trace.pasting.inside.condition = combine.condition(trace.pasting.inside.condition, new.pasting.inside.condition)
    bidx.obs.outside.to.paste = query.from.condition(d.outside, trace.pasting.inside.condition, g.str, l.str)
    d.inside = rbind(d.inside, d.outside[bidx.obs.outside.to.paste,]);
    d.outside = d.outside[!bidx.obs.outside.to.paste,];


    pos.group.pasting.list[idx.pos.group.pasting.list]=list(trace.pasting.inside.condition)
    idx.pos.group.pasting.list=idx.pos.group.pasting.list+1
  }

  ## one dropping
  trace.dropping.inside.condition = trace.pasting.inside.condition
  repeat {
    one.dropping.results.ordered = one.dropping(d.inside, d.outside, trace.dropping.inside.condition, yvar, censorvar, trtvar, g.str, l.str, type, des.res);

    if(dim(one.dropping.results.ordered)[1] == 0) break

    new.dropping.inside.condition = one.dropping.results.ordered[1,];
    trace.dropping.inside.condition = combine.condition(trace.dropping.inside.condition, new.dropping.inside.condition)

    bidx.obs.outside.to.add = query.from.condition(d.outside, trace.dropping.inside.condition, g.str, l.str)
    d.inside = rbind(d.inside, d.outside[bidx.obs.outside.to.add,]);
    d.outside = d.outside[!bidx.obs.outside.to.add,];

    pos.group.dropping.list[idx.pos.group.dropping.list]=list(trace.dropping.inside.condition)
    idx.pos.group.dropping.list=idx.pos.group.dropping.list+1
  }

  pos.group.list = c(pos.group.peeling.list, pos.group.pasting.list, pos.group.dropping.list)



  ## apply candidates to testing data
  num.candidates = length (pos.group.list)
  test.result = data.frame(idx.candidate = as.numeric(), num.inside=as.numeric(), num.outside=as.numeric(),
                           pv1=as.numeric(), pv2=as.numeric(),
                           coef1=as.numeric(), coef2=as.numeric(), stringsAsFactors=F)
  pidx.testing.results = 1

  if (num.candidates==0) return (list())

  for (i in 1:num.candidates){
    idx.obs.test.inside = query.from.condition(data.test, pos.group.list[[i]], g.str, l.str)
    d.test.inside = data.test[idx.obs.test.inside,]
    d.test.outside = data.test[!idx.obs.test.inside,]
    num.inside=dim(d.test.inside)[1]
    num.outside=dim(d.test.outside)[1]

    res1 = try(pval.cal(d.test.inside, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res), silent=T)
    res2 = try(pval.cal(d.test.outside, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res), silent=T)

    if(is.na(res1$pv) || is.na(res2$pv) || class(res1)=="try-error" || class(res2)=="try-error") next

    test.result[pidx.testing.results, "idx.candidate"]=i
    test.result[pidx.testing.results, "num.inside"]=num.inside
    test.result[pidx.testing.results, "num.outside"]=num.outside
    test.result[pidx.testing.results, "pv1"]=res1$pv
    test.result[pidx.testing.results, "pv2"]=res2$pv
    test.result[pidx.testing.results, "coef1"]=res1$coef.trt
    test.result[pidx.testing.results, "coef2"]=res2$coef.trt

    pidx.testing.results = pidx.testing.results + 1
  }

  if (dim(test.result)[1]==0) return (list())

  idx.obs.qualified=test.result$pv1<test.result$pv2

  test.result = test.result[idx.obs.qualified,]
  if (dim(test.result)[1]==0) return (list())

  order.by.pv1 = order(test.result$pv1, decreasing=F)
  test.result.ordered = test.result[order.by.pv1, ]

  #final selected candidate
  idx.selected.candidate = test.result.ordered[1, "idx.candidate"]
  selected.signature = pos.group.list[[idx.selected.candidate]]
  result.final = list(test.result.ordered = test.result.ordered, signature = selected.signature, pos.group.list=pos.group.list)

  result.final

}





#' Perform peeling one time in predictive case.
#'
#' @param d.inside the dataset for subjects in consideration.
#' @param d.outside the dataset for subjects outside consideration.
#' @param xvars the vector of variable names for predictors (covariates).
#' @param alpha a parameter controlling the number of patients in consideration
#' @param min.size.inside desired number of subjects in signature positive group size for a given cutoff.
#' @param yvar the name of response variable.
#' @param censorvar the name of censoring variable (1: event; 0: censor), default = NULL).
#' @param trtvar 0-1 coded vector for treatment variable.
#' @param g.str ">=".
#' @param l.str "<=".
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#'
#' @return a data frame enlisting the signature rules (after peeling) ordered by treatment p-values in each group defined by the rules
#'
one.peeling <- function(d.inside, d.outside, xvars, alpha, min.size.inside, yvar, censorvar, trtvar, g.str, l.str, type,des.res) {
  ### iterate all possible x, peel a piece of x, compute two-sample comparison
  ### event.indicator.nm = "os.censor.0"; event.time.nm = "os";
  ### trt.nm = "trt"; control.trt.value = 0; treatment.trt.value = 1;
  ### g.str = ">="; l.str = "<=";

  one.peeling.results = data.frame(x.nm=as.character(), condition=as.character(), cutoff=as.numeric(), pv1=as.numeric(), pv2=as.numeric(),
                                   coef1=as.numeric(),coef2=as.numeric(), end.flag=as.logical(), stringsAsFactors=F);

  # check whether the d.inside after peeling has reached the limit of min number of obs
  n.obs.would.be.left = dim(d.inside)[1]*(1 - alpha);
  if(n.obs.would.be.left < min.size.inside) {
    return(one.peeling.results);
  }


  pidx.peeling.results = 1;
  num.xvars = length(xvars);
  for(i in 1:num.xvars) { # one peeling
    x.nm = xvars[i];
    cutoffs = quantile(d.inside[, x.nm], probs=c(alpha, 1-alpha), na.rm=TRUE); # PT: Added na.rm=TRUE to handle NAs
    inside.condition.strs = c(g.str, l.str);
    num.cutoffs = length(cutoffs);
    for(j in 1:num.cutoffs) { # for one of the cutoffs, compare trt vs. ctrl in d.inside.shrinked (and d.outside.extended)
      cutoff = cutoffs[j];
      inside.condition.str = inside.condition.strs[j];
      if(inside.condition.str == g.str) {
        bidx.x.inside = ifelse(d.inside[,x.nm] >= cutoff, TRUE, FALSE);
      } else {
        bidx.x.inside = ifelse(d.inside[,x.nm] <= cutoff, TRUE, FALSE);
      }
      # peel and compare
      if (dim(d.inside[!bidx.x.inside, ])[1]==0) next
      d.inside.shrinked = d.inside[bidx.x.inside, ];
      d.outside.extended = rbind(d.inside[!bidx.x.inside, ], d.outside);
      if (dim(d.inside.shrinked)[1] < min.size.inside) next

      res1 = try(pval.cal(d.inside.shrinked, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res),silent=T)
      res2 = try(pval.cal(d.outside.extended, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res),silent=T)

      # skip collecting result if not comparable
      if(is.na(res1$pv) || is.na(res2$pv)||class(res1)=="try-error"||class(res2)=="try-error") next

      # save the result of peeling of x on its one side
      one.peeling.results[pidx.peeling.results, "x.nm"] = x.nm;
      one.peeling.results[pidx.peeling.results, "condition"] = inside.condition.str;
      one.peeling.results[pidx.peeling.results, "cutoff"] = cutoff;
      one.peeling.results[pidx.peeling.results, "pv1"] = res1$pv;
      one.peeling.results[pidx.peeling.results, "pv2"] = res2$pv;
      one.peeling.results[pidx.peeling.results, "coef1"] = res1$coef.trt;
      one.peeling.results[pidx.peeling.results, "coef2"] = res2$coef.trt;
      one.peeling.results[pidx.peeling.results, "end.flag"] = FALSE
      pidx.peeling.results = pidx.peeling.results + 1;
    }
  }
  if (dim(one.peeling.results)[1]==0) return (one.peeling.results)

  # filter out unqualified peeling
  idx.obs.qualified=one.peeling.results$pv1<one.peeling.results$pv2
  one.peeling.results = one.peeling.results[idx.obs.qualified, ];
  order.by.pv1 = order(one.peeling.results$pv1, decreasing=F);
  one.peeling.results.ordered = one.peeling.results[order.by.pv1,];
  #rownames(one.peeling.results.ordered) = 1:(dim(one.peeling.results.ordered)[1]); # avoid the top one cannot be referred by [1,]
  return(one.peeling.results.ordered);
} # end of one.peeling






#' create training/testing dataset indexes.
#'
#' @param training.percent percentage of subjects in training data as mentioned in prim.train function.
#' @param n number of sbjects in the whole dataset.
#'
#' @return a list containing training and test data indices.
#'
create.training.dataset.index <- function(training.percent, n) {

  if(training.percent<1 & training.percent >0){
    pos = 1:n;
    pidx.perm = sample(pos);
    number.training.examples = floor(n * training.percent);

    # get data position indexes
    pidx.train = pidx.perm[1:number.training.examples];
    pidx.test = pidx.perm[(number.training.examples+1):n];
    return(list(pidx.test=pidx.test, pidx.train=pidx.train));
  }else if (training.percent ==1){
    pidx.train = 1:n
    pidx.test = 1:n
    return(list(pidx.test=pidx.test, pidx.train=pidx.train));
  }else{
    stop("Training set is not provided.")
  }

}





#' Calculate p-value for treatment in each subgroup in predictive case
#'
#' @param data input data frame
#' @param yvar name for response variable
#' @param censorvar name for censoring (1: event; 0: censor), default = NULL
#' @param trtvar name for treatment variable, default = NULL (prognostic signature)
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response
#'
#' @return p-value for the treatment given the dataset
#'
pval.cal <- function(data, yvar, censorvar, trtvar, type, des.res) {
  # avoid no observation in either group
  num.obs.group0 = sum(data[, trtvar] == 0);
  num.obs.group1 = sum(data[, trtvar] == 1);
  if(num.obs.group0 == 0 || num.obs.group1 == 0) {
    return(list(pv=NA, coef.trt=NA));
  }

  if (type=="s"){
    fit.formula <- as.formula(paste("Surv(", yvar, ", ", censorvar, ")~", trtvar, sep=""))
    fit.obj <- coxph(fit.formula, data=data)
    pv = summary(fit.obj)$coefficients[1,5]/2
    coef.trt=summary(fit.obj)$coefficients[1,1]
    if (des.res=="larger" & coef.trt>0) pv=1-pv
    if (des.res=="smaller" & coef.trt<0) pv=1-pv

  }

  fit.formula = as.formula(paste(yvar, "~", trtvar, sep=""))
  if (type=="b"){
    fit.obj=glm(fit.formula,family=binomial,data=data)
    pv=summary(fit.obj)$coefficients[2,4]/2
    coef.trt=summary(fit.obj)$coefficients[2,1]
    if (des.res=="larger" & coef.trt<0) pv=1-pv
    if (des.res=="smaller" & coef.trt>0) pv=1-pv
  }

  if (type=="c"){
    fit.obj=lm(fit.formula,data=data)
    pv=summary(fit.obj)$coefficients[2,4]/2
    coef.trt=summary(fit.obj)$coefficients[2,1]
    if (des.res=="larger" & coef.trt<0) pv=1-pv
    if (des.res=="smaller" & coef.trt>0) pv=1-pv
  }


  list(pv=pv, coef.trt=coef.trt)

}




#' Internal function
#'
#' @param trace.inside.condition list of signature rules
#' @param new.inside.condition new signature rule
#'
#' @return updated list of signature rules
#'
combine.condition = function(trace.inside.condition, new.inside.condition) {
  end.flag = new.inside.condition[, "end.flag"]
  new.inside.condition=new.inside.condition[,c("x.nm","condition","cutoff")]
  new.x.glcondition=paste(new.inside.condition[,"x.nm"],new.inside.condition[,"condition"])
  if(dim(trace.inside.condition)[1]==0){
    existing.x.glcondition=NULL
  }else{
    existing.x.glcondition=paste(trace.inside.condition[,"x.nm"],trace.inside.condition[,"condition"])
  }


  if (!end.flag){
    if (new.x.glcondition %in% existing.x.glcondition) {
      trace.inside.condition[existing.x.glcondition==new.x.glcondition,]=new.inside.condition
    }else{
      trace.inside.condition=rbind(trace.inside.condition,new.inside.condition)
    }
  }else{
    if (new.x.glcondition %in% existing.x.glcondition) {
      trace.inside.condition=trace.inside.condition[existing.x.glcondition!=new.x.glcondition,]
    }else{
      trace.inside.condition=trace.inside.condition
    }
  }

  trace.inside.condition
}




#' Perform pasting one time in predictive case.
#'
#' @param d.inside the dataset for subjects in consideration after peeling.
#' @param d.outside the dataset for subjects outside consideration after peeling.
#' @param trace.inside.condition list of signature rules used for d.inside.
#' @param alpha a parameter controlling the number of subjects in consideration.
#' @param yvar the name for response variable.
#' @param censorvar the name for censoring (1: event; 0: censor), default = NULL.
#' @param trtvar 0-1 coded vector for treatment variable.
#' @param g.str ">=".
#' @param l.str "<=".
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#'
#' @return a data frame enlisting the signature rules (after pasting) ordered by treatment p-values in each group defined by the rules.
#'
one.pasting = function (d.inside, d.outside, trace.inside.condition, alpha, yvar, censorvar, trtvar, g.str, l.str, type, des.res) {
  n.obs.inside = dim(d.inside)[1]
  n.obs.outside = dim(d.outside)[1]
  n.obs.to.add = ceiling(n.obs.inside * alpha)

  one.pasting.results = data.frame(x.nm=as.character(), condition=as.character(), cutoff=as.numeric(), pv1=as.numeric(),
                                   pv2=as.numeric(), coef1=as.numeric(), coef2=as.numeric(),end.flag=as.logical(),stringsAsFactors=F);

  pidx.pasting.results = 1;
  num.x.pasting=dim(trace.inside.condition)[1]

  for (i in 1:num.x.pasting){
    relax.condition=trace.inside.condition[i,]
    x.nm.relax.condition=relax.condition[,"x.nm"]
    other.condition=trace.inside.condition[-i,]
    idx.outside.by.other.condition=query.from.condition(d.outside, other.condition, g.str, l.str)
    d.outside.by.other.condition=d.outside[idx.outside.by.other.condition,]

    if(dim(d.outside.by.other.condition)[1] == 0) next

    if (relax.condition[,"condition"]==g.str){
      d.outside.by.other.condition.order = order(d.outside.by.other.condition[,x.nm.relax.condition],decreasing=T)
      d.outside.by.other.condition.ordered = d.outside.by.other.condition[d.outside.by.other.condition.order,]
      x.cutoff.relax.condition = d.outside.by.other.condition.ordered[,x.nm.relax.condition][min(n.obs.to.add,dim(d.outside.by.other.condition.ordered)[1])]
      n.obs.actually.add = sum(d.outside.by.other.condition.ordered[,x.nm.relax.condition]>= x.cutoff.relax.condition)

      end.flag = n.obs.actually.add>=dim(d.outside.by.other.condition.ordered)[1]

      d.inside.extended=rbind(d.inside, d.outside.by.other.condition.ordered[1:n.obs.actually.add,])
      d.outside.shrinked=rbind(d.outside[!idx.outside.by.other.condition,],d.outside.by.other.condition.ordered[-(1:n.obs.actually.add),])
      paste.x.extended.cutoff = d.outside.by.other.condition.ordered[n.obs.actually.add,x.nm.relax.condition]
      cutoff.inside.condition.str = g.str;
    }

    if (relax.condition[,"condition"]==l.str){
      d.outside.by.other.condition.order = order(d.outside.by.other.condition[,x.nm.relax.condition],decreasing=F)
      d.outside.by.other.condition.ordered = d.outside.by.other.condition[d.outside.by.other.condition.order,]
      x.cutoff.relax.condition = d.outside.by.other.condition.ordered[,x.nm.relax.condition][min(n.obs.to.add,dim(d.outside.by.other.condition.ordered)[1])]
      n.obs.actually.add = sum(d.outside.by.other.condition.ordered[,x.nm.relax.condition]<= x.cutoff.relax.condition)

      end.flag = n.obs.actually.add>=dim(d.outside.by.other.condition.ordered)[1]

      d.inside.extended=rbind(d.inside, d.outside.by.other.condition.ordered[1:n.obs.actually.add,])
      d.outside.shrinked=rbind(d.outside[!idx.outside.by.other.condition,],d.outside.by.other.condition.ordered[-(1:n.obs.actually.add),])
      paste.x.extended.cutoff = d.outside.by.other.condition.ordered[n.obs.actually.add,x.nm.relax.condition]
      cutoff.inside.condition.str = l.str;
    }

    res1 = try(pval.cal(d.inside.extended, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res), silent=T)
    res2 = try(pval.cal(d.outside.shrinked, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res), silent=T)

    # skip collecting result if not comparable
    if(is.na(res1$pv) || is.na(res2$pv)||class(res1)=="try-error"||class(res2)=="try-error") next

    # save the result of pasting of x on its one side
    one.pasting.results[pidx.pasting.results, "x.nm"] = x.nm.relax.condition
    one.pasting.results[pidx.pasting.results, "condition"] = cutoff.inside.condition.str
    one.pasting.results[pidx.pasting.results, "cutoff"] = paste.x.extended.cutoff;
    one.pasting.results[pidx.pasting.results, "pv1"] = res1$pv
    one.pasting.results[pidx.pasting.results, "pv2"] = res2$pv
    one.pasting.results[pidx.pasting.results, "coef1"] = res1$coef.trt
    one.pasting.results[pidx.pasting.results, "coef2"] = res2$coef.trt
    one.pasting.results[pidx.pasting.results, "end.flag"] = end.flag
    pidx.pasting.results = pidx.pasting.results + 1;

  }

  if (dim(one.pasting.results)[1]==0) return (one.pasting.results)

  idx.obs.qualified=one.pasting.results$pv1<one.pasting.results$pv2
  pv.current = pval.cal(d.inside, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type,des.res=des.res)$pv
  idx.obs.qualified = idx.obs.qualified & (one.pasting.results$pv1 < pv.current)

  one.pasting.results=one.pasting.results[idx.obs.qualified,]
  order.by.pv1 = order(one.pasting.results$pv1, decreasing=F)
  one.pasting.results.ordered = one.pasting.results[order.by.pv1, ]
  return(one.pasting.results.ordered)

}




#' An internal function inside one.pasting.
#'
#' @param d dataset for subjects in consideration.
#' @param condition signature rule in consideration.
#' @param g.str ">="
#' @param l.str "<="
#'
#' @return a vector of logical arguments indicating whether the conditions can be  satisfied for the subjects in d.
#'
query.from.condition <- function(d, condition, g.str=">=", l.str="<=") {
  idx.bound.by.condition=rep(TRUE, dim(d)[1])
  if (dim(condition)[1]==0) return(idx.bound.by.condition)

  for (i in 1:dim(condition)[1]){
    condition.i=condition[i,]
    x.nm=condition.i[,"x.nm"]
    if (condition.i[,"condition"]==g.str) {
      idx.bound.by.condition=idx.bound.by.condition & (d[,x.nm]>=condition.i[,"cutoff"])
    } else {
      idx.bound.by.condition=idx.bound.by.condition & (d[,x.nm]<=condition.i[,"cutoff"])
    }
  }
  idx.bound.by.condition
}





#' Perform dropping one time in predictive case.
#'
#' @param d.inside the dataset for subjects in consideration after pasting.
#' @param d.outside the dataset for subjects outside consideration after pasting.
#' @param trace.inside.condition list of signature rules used for d.inside.
#' @param yvar the name for response variable.
#' @param censorvar the name for censoring (1: event; 0: censor), default = NULL.
#' @param trtvar 0-1 coded vector for treatment variable.
#' @param g.str ">=".
#' @param l.str "<=".
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#'
#' @return a data frame enlisting the signature rules (after dropping) ordered by treatment p-values in each group defined by the rules.
#'
one.dropping = function (d.inside, d.outside, trace.inside.condition, yvar, censorvar, trtvar, g.str, l.str, type, des.res){
  one.dropping.results = data.frame(x.nm=as.character(), condition=as.character(), cutoff=as.numeric(), pv1=as.numeric(),
                                    pv2=as.numeric(), coef1=as.numeric(), coef2=as.numeric(),end.flag=as.logical(),stringsAsFactors=F);
  pidx.dropping.results = 1;
  num.x.dropping=dim(trace.inside.condition)[1]

  if (num.x.dropping<=1) return (one.dropping.results)

  for (i in 1:num.x.dropping){
    drop.condition=trace.inside.condition[i,]
    x.nm.drop.condition=drop.condition[,"x.nm"]
    other.condition=trace.inside.condition[-i,]
    idx.outside.by.other.condition=query.from.condition(d.outside, other.condition, g.str, l.str)
    d.outside.by.other.condition=d.outside[idx.outside.by.other.condition,]

    if(dim(d.outside.by.other.condition)[1] == 0) next

    if (drop.condition[,"condition"]==g.str){
      d.inside.extended = rbind(d.inside, d.outside.by.other.condition)
      d.outside.shrinked = d.outside[!idx.outside.by.other.condition,]
      x.extended.cutoff = min(d.inside.extended[,x.nm.drop.condition])
      cutoff.inside.condition.str = g.str;
    }

    if (drop.condition[,"condition"]==l.str){
      d.inside.extended=rbind(d.inside, d.outside.by.other.condition)
      d.outside.shrinked=d.outside[!idx.outside.by.other.condition,]
      x.extended.cutoff = max(d.inside.extended[,x.nm.drop.condition])
      cutoff.inside.condition.str = l.str;
    }

    res1 = try(pval.cal(d.inside.extended, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res),silent=T)
    res2 = try(pval.cal(d.outside.shrinked, yvar=yvar, censorvar=censorvar, trtvar=trtvar, type=type, des.res=des.res),silent=T)

    # skip collecting result if not comparable
    if(is.na(res1$pv) || is.na(res2$pv)||class(res1)=="try-error"||class(res2)=="try-error") next

    # save the result of pasting of x on its one side
    one.dropping.results[pidx.dropping.results, "x.nm"] = x.nm.drop.condition
    one.dropping.results[pidx.dropping.results, "condition"] = cutoff.inside.condition.str
    one.dropping.results[pidx.dropping.results, "cutoff"] = x.extended.cutoff;
    one.dropping.results[pidx.dropping.results, "pv1"] = res1$pv
    one.dropping.results[pidx.dropping.results, "pv2"] = res2$pv
    one.dropping.results[pidx.dropping.results, "coef1"] = res1$coef.trt
    one.dropping.results[pidx.dropping.results, "coef2"] = res2$coef.trt
    one.dropping.results[pidx.dropping.results, "end.flag"] = TRUE
    pidx.dropping.results = pidx.dropping.results + 1;

  }

  if (dim(one.dropping.results)[1]==0) return (one.dropping.results)

  idx.obs.qualified=one.dropping.results$pv1<one.dropping.results$pv2

  one.dropping.results=one.dropping.results[idx.obs.qualified,]
  order.by.pv1 = order(one.dropping.results$pv1, decreasing=F)
  one.dropping.results.ordered = one.dropping.results[order.by.pv1, ]
  return(one.dropping.results.ordered)

}




################# cv prim ######################
#' Cross-validation for PRIM
#'
#' @param data the input data frame
#' @param yvar name for response variable
#' @param censorvar name for censoring (1: event; 0: censor), default = NULL
#' @param trtvar name for treatment variable, default = NULL (prognostic signature)
#' @param trtref coding (in the column of trtvar) for treatment arm
#' @param xvars vector of variable names for predictors (covariates)
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response
#' @param alpha a parameter controlling the number of patients in consideration
#' @param min.sigp.prcnt desired proportion of signature positive group size for a given cutoff.
#' @param training.percent percentage of subjects in the initial training data
#' @param n.boot number of bootstrap for the variable selection procedure for PRIM
#' @param pre.filter NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
#' @param filter.method NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case
#' @param k.fold number of folds for CV.
#' @param cv.iter Algorithm terminates after cv.iter successful iterations of cross-validation
#' @param max.iter total number of iterations allowed (including unsuccessful ones)
#'
#' @return a list containing with following entries:
#' \item{stats.summary}{Summary of performance statistics.}
#' \item{pred.classes}{Data frame containing the predictive clases (TRUE/FALSE) for each iteration.}
#' \item{folds}{Data frame containing the fold indices (index of the fold for each row) for each iteration.}
#' \item{sig.list}{List of length cv.iter * k.fold containing the signature generated at each of the  k folds, for all iterations.}
#' \item{error.log}{List of any error messages that are returned at an iteration.}
#' \item{interplot}{Treatment*subgroup interaction plot for predictive case}
#'
prim.cv = function(data, yvar, censorvar, trtvar, trtref=NULL, xvars, type, des.res,
                        alpha = c(0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50),
                   min.sigp.prcnt = 0.2,  training.percent = 0.50, n.boot=0, pre.filter=NULL, filter.method=NULL, k.fold=5, cv.iter=50, max.iter=500) {

  data = as.data.frame(data)
  y = data[,yvar,drop=F]
  if (!is.null(censorvar)) {
    censor.vec = data[,censorvar,drop=F]
  } else {
    censor.vec = NULL
  }

  if (!is.null(trtvar))
  {
    if (!is.null(trtref)){
      data[,trtvar]=(data[,trtvar]==trtref)+0
      trtref=NULL
    }
    trt.vec=data[,trtvar,drop=F]
  } else
  {
    trt.vec<-NULL
  }

  if (type=="b" | type=="rb") strata=data[,yvar]
  if (type=="s") strata=censor.vec[, censorvar]
  if (type=="c") strata=NULL
  model.Rfunc = "prim.train.wrapper"
  model.Rfunc.args = list(yvar=yvar, xvars=xvars, censorvar=censorvar, trtvar=trtvar, trtref=trtref,
                          type=type, n.boot=n.boot, des.res=des.res, min.sigp.prcnt=min.sigp.prcnt,
                          alpha=alpha, training.percent=training.percent, pre.filter=pre.filter, filter.method=filter.method)

  predict.Rfunc = "pred.prim.cv"
  predict.Rfunc.args = list(yvar=yvar)

  res <- kfold.cv(data=data, model.Rfunc=model.Rfunc, model.Rfunc.args=model.Rfunc.args, predict.Rfunc=predict.Rfunc, predict.Rfunc.args=predict.Rfunc.args, k.fold=k.fold, cv.iter=cv.iter, strata=strata, max.iter=max.iter)
  if (length(res) > 0) {
    stats = evaluate.cv.results(cv.data=res$cv.data, y=y, censor.vec=censor.vec, trt.vec=trt.vec, type=type)
    summary = summarize.cv.stats(stats$raw.stats, trtvar, type)
    interplot=interaction.plot(data.eval=summary, type=type, main="Interaction Plot", trt.lab=c("Trt.", "Ctrl."))
    results <- list(stats.summary=summary, pred.classes=stats$pred.classes, folds=stats$folds, sig.list=res$sig.list, raw.stats=stats$raw.stats, error.log=res$error.log, interplot=interplot)

  } else {
    results <- "No successful cross-validations."
  }

  return(results)

}




#' Wrapper function for PRIM CV
#'
#' @param data input data frame
#' @param args list containing all other input arguments to prim.train except for x and y.
#'
#' @return prediction rule as returned by prim.train
#'
prim.train.wrapper=function(data, args){

  yvar<-censorvar<-trtvar<-trtref<-xvars<-type<-des.res<-alpha<-min.sigp.prcnt<-training.percent<-n.boot<-pre.filter<-filter.method<-NULL

  for (name in names(args)) {
    assign(name, args[[name]])
  }

  res=prim.train(data=data, yvar=yvar, censorvar=censorvar, trtvar=trtvar,
                      trtref=trtref, xvars=xvars,type=type, des.res=des.res,
                      alpha=alpha, min.sigp.prcnt=min.sigp.prcnt,
                      training.percent = training.percent, n.boot=n.boot,
                      pre.filter=pre.filter, filter.method=filter.method)
  res
}




#' Prediction function for PRIM CV
#'
#' @param data input data frame
#' @param predict.rule signature rules as returned by prim.train
#' @param args list of the form list(yvar=yvar)
#'
#' @return The input data with an added column, a logical vector indicating the prediction for each row of data.
#'
pred.prim.cv <- function(data, predict.rule, args) {
  pred.id = query.from.condition(d=data, condition=predict.rule, g.str=">=", l.str="<=")
  pred.data = cbind(data, data.frame(pred.class=pred.id))
  pred.data = subset(pred.data, select="pred.class")
  pred.data
}




#' Prediction function for PRIM
#'
#' @param data input data frame (only covariates)
#' @param predict.rule signature rules returned by prim.train
#'
#' @return The input data with an added column, a logical vector indicating the prediction for each row of data.
#'
pred.prim <- function(data, predict.rule) {
  pred.id = query.from.condition(d=data, condition=predict.rule, g.str=">=", l.str="<=")
  pred.data = cbind(data, data.frame(pred.class=pred.id))
  pred.data
}


#################### TRAIN ONCE ###############
#' Apply PRIM one time on the training data for a fixed value of alpha (in prognostic case)
#' @description this function applies the prim procedure (peeling, pasting, and dropping operations) on training data one time.
#'
#' @param data input data frame
#' @param yvar name for response variable
#' @param censorvar name for censoring (1: event; 0: censor), default = NULL
#' @param xvars vector of variable names for predictors (covariates)
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response
#' @param alpha a parameter controlling the number of patients in consideration
#' @param min.size.inside desired number of subjects in signature positive group size for a given cutoff.
#' @param pidx.train.test training and test data index as obtained from create.training.dataset.index
#'
#' @return a list containing signature rules and test result based on the signatures.
#'
prim.train.prog.once <- function(
  # data info
  data, yvar, censorvar, xvars,type, des.res,
  # control arguments
  alpha = 0.10,  min.size.inside = 20, pidx.train.test) {
  ### this function is the prim procedure calling peeling, pasting, and dropping operations.

  ## method constants
  g.str = ">=";
  l.str = "<=";


  ## create training and testing data sets
  data.train = data[pidx.train.test$pidx.train, ]; # training data
  data.test = data[pidx.train.test$pidx.test, ]; # testing data


  d.inside = data.train;
  d.outside = data.frame();


  # data structure to save results
  trace.peeling.inside.condition = data.frame(); # the most recent updated box condition
  trace.pasting.inside.condition = data.frame()
  trace.dropping.inside.condition = data.frame()

  pos.group.peeling.list = list();
  idx.pos.group.peeling.list = 1;
  pos.group.pasting.list = list();
  idx.pos.group.pasting.list = 1;
  pos.group.dropping.list = list();
  idx.pos.group.dropping.list = 1;

  # peeling
  repeat {
    one.peeling.results.ordered = one.peeling.prog(d.inside, d.outside, xvars, alpha, min.size.inside, yvar, censorvar, g.str, l.str, type, des.res);
    if(dim(one.peeling.results.ordered)[1] == 0) break
    new.inside.condition = one.peeling.results.ordered[1,];
    trace.peeling.inside.condition = combine.condition(trace.peeling.inside.condition, new.inside.condition)

    x.nm = new.inside.condition[1, "x.nm"];
    x.cutoff = new.inside.condition[1, "cutoff"]

    if(new.inside.condition[1, "condition"] == g.str) {
      bidx.obs.inside.to.peel = d.inside[, x.nm] >= x.cutoff;
    } else if((new.inside.condition[1, "condition"] == l.str)) {
      bidx.obs.inside.to.peel = d.inside[, x.nm] <= x.cutoff;
    }

    # NOTE: need to update outside data first
    d.outside = rbind(d.outside, d.inside[!bidx.obs.inside.to.peel, ]);
    d.inside = d.inside[bidx.obs.inside.to.peel, ];

    pos.group.peeling.list[idx.pos.group.peeling.list]=list(trace.peeling.inside.condition)
    idx.pos.group.peeling.list=idx.pos.group.peeling.list+1
  }

  if (dim(trace.peeling.inside.condition)[1]==0) return(list())

  #pasting
  trace.pasting.inside.condition = trace.peeling.inside.condition
  repeat {
    one.pasting.results.ordered = one.pasting.prog(d.inside, d.outside, trace.pasting.inside.condition, alpha, yvar, censorvar, g.str, l.str, type, des.res);
    if(dim(one.pasting.results.ordered)[1] == 0) break

    new.pasting.inside.condition = one.pasting.results.ordered[1,];
    temp.trace.pasting.inside.condition = combine.condition(trace.pasting.inside.condition, new.pasting.inside.condition)
    if (dim(temp.trace.pasting.inside.condition)[1]==0) break

    trace.pasting.inside.condition = combine.condition(trace.pasting.inside.condition, new.pasting.inside.condition)
    bidx.obs.outside.to.paste = query.from.condition(d.outside, trace.pasting.inside.condition, g.str, l.str)
    d.inside = rbind(d.inside, d.outside[bidx.obs.outside.to.paste,]);
    d.outside = d.outside[!bidx.obs.outside.to.paste,];


    pos.group.pasting.list[idx.pos.group.pasting.list]=list(trace.pasting.inside.condition)
    idx.pos.group.pasting.list=idx.pos.group.pasting.list+1
  }

  ## one dropping
  trace.dropping.inside.condition = trace.pasting.inside.condition
  repeat {
    one.dropping.results.ordered = one.dropping.prog(d.inside, d.outside, trace.dropping.inside.condition, yvar, censorvar, g.str, l.str, type, des.res);

    if(dim(one.dropping.results.ordered)[1] == 0) break

    new.dropping.inside.condition = one.dropping.results.ordered[1,];
    trace.dropping.inside.condition = combine.condition(trace.dropping.inside.condition, new.dropping.inside.condition)

    bidx.obs.outside.to.add = query.from.condition(d.outside, trace.dropping.inside.condition, g.str, l.str)
    d.inside = rbind(d.inside, d.outside[bidx.obs.outside.to.add,]);
    d.outside = d.outside[!bidx.obs.outside.to.add,];

    pos.group.dropping.list[idx.pos.group.dropping.list]=list(trace.dropping.inside.condition)
    idx.pos.group.dropping.list=idx.pos.group.dropping.list+1
  }

  pos.group.list = c(pos.group.peeling.list, pos.group.pasting.list, pos.group.dropping.list)



  ## apply candidates to testing data
  num.candidates = length (pos.group.list)
  test.result = data.frame(idx.candidate = as.numeric(), num.inside=as.numeric(), num.outside=as.numeric(),
                           pv=as.numeric(), coef=as.numeric(), stringsAsFactors=F)
  pidx.testing.results = 1

  if (num.candidates==0) return (list())

  for (i in 1:num.candidates){
    grp.id.test = query.from.condition(data.test, pos.group.list[[i]], g.str, l.str)
    d.test.inside = data.test[grp.id.test,]
    d.test.outside = data.test[!grp.id.test,]
    num.inside=dim(d.test.inside)[1]
    num.outside=dim(d.test.outside)[1]

    res = try(pval.cal.prog(data.test, yvar=yvar, censorvar=censorvar, grp.id=grp.id.test, type=type, des.res=des.res), silent=T)

    if(is.na(res$pv) || class(res)=="try-error" ) next

    test.result[pidx.testing.results, "idx.candidate"]=i
    test.result[pidx.testing.results, "num.inside"]=num.inside
    test.result[pidx.testing.results, "num.outside"]=num.outside
    test.result[pidx.testing.results, "pv"]=res$pv
    test.result[pidx.testing.results, "coef"]=res$coef.grp

    pidx.testing.results = pidx.testing.results + 1
  }

  if (dim(test.result)[1]==0) return (list())

  order.by.pv = order(test.result$pv, decreasing=F)
  test.result.ordered = test.result[order.by.pv, ]

  #final selected candidate
  idx.selected.candidate = test.result.ordered[1, "idx.candidate"]
  selected.signature = pos.group.list[[idx.selected.candidate]]
  result.final = list(test.result.ordered = test.result.ordered, signature = selected.signature, pos.group.list=pos.group.list)

  result.final

}


#' Perform peeling one time in prognostic case.
#'
#' @param d.inside the dataset for subjects in consideration.
#' @param d.outside the dataset for subjects outside consideration.
#' @param xvars the vector of variable names for predictors (covariates).
#' @param alpha a parameter controlling the number of patients in consideration
#' @param min.size.inside desired number of subjects in signature positive group size for a given cutoff.
#' @param yvar the name of response variable.
#' @param censorvar the name of censoring variable (1: event; 0: censor), default = NULL).
#' @param g.str ">=".
#' @param l.str "<=".
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#'
#' @return a data frame enlisting the signature rules (after peeling) ordered by main effect p-values in each group defined by the rules.
#'
one.peeling.prog <- function(d.inside, d.outside, xvars, alpha, min.size.inside, yvar, censorvar, g.str, l.str, type,des.res) {
  ### one.peeling for prognostic case
  ### iterate all possible x, peel a piece of x, compute two-sample comparison
  ### event.indicator.nm = "os.censor.0"; event.time.nm = "os";
  ### trt.nm = "trt"; control.trt.value = 0; treatment.trt.value = 1;
  ### g.str = ">="; l.str = "<=";


  one.peeling.results = data.frame(x.nm=as.character(), condition=as.character(), cutoff=as.numeric(), pv=as.numeric(), coef=as.numeric(), end.flag=as.logical(), stringsAsFactors=F);

  # check whether the d.inside after peeling has reached the limit of min number of obs
  n.obs.would.be.left = dim(d.inside)[1]*(1 - alpha);
  if(n.obs.would.be.left < min.size.inside) {
    return(one.peeling.results);
  }


  pidx.peeling.results = 1;
  num.xvars = length(xvars);
  for(i in 1:num.xvars) { # one peeling
    x.nm = xvars[i];
    cutoffs = quantile(d.inside[, x.nm], probs=c(alpha, 1-alpha), na.rm=TRUE); # PT: Added na.rm=TRUE to handle NAs
    inside.condition.strs = c(g.str, l.str);
    num.cutoffs = length(cutoffs);
    for(j in 1:num.cutoffs) {
      cutoff = cutoffs[j];
      inside.condition.str = inside.condition.strs[j];
      if(inside.condition.str == g.str) {
        grp.id = ifelse(d.inside[,x.nm] >= cutoff, TRUE, FALSE);
      } else {
        grp.id = ifelse(d.inside[,x.nm] <= cutoff, TRUE, FALSE);
      }
      # peel and compare
      if (dim(d.inside[!grp.id, ])[1]==0) next

      res = try(pval.cal.prog(d.inside, yvar=yvar, censorvar=censorvar, grp.id=grp.id, type=type, des.res=des.res),silent=T)
      # skip collecting result if not comparable
      if(is.na(res$pv) ||class(res)=="try-error") next

      # save the result of peeling of x on its one side
      one.peeling.results[pidx.peeling.results, "x.nm"] = x.nm;
      one.peeling.results[pidx.peeling.results, "condition"] = inside.condition.str;
      one.peeling.results[pidx.peeling.results, "cutoff"] = cutoff;
      one.peeling.results[pidx.peeling.results, "pv"] = res$pv;
      one.peeling.results[pidx.peeling.results, "coef"] = res$coef.grp;
      one.peeling.results[pidx.peeling.results, "end.flag"] = FALSE
      pidx.peeling.results = pidx.peeling.results + 1;
    }
  }
  if (dim(one.peeling.results)[1]==0) return (one.peeling.results)

  order.by.pv = order(one.peeling.results$pv, decreasing=F);
  one.peeling.results.ordered = one.peeling.results[order.by.pv,];
  return(one.peeling.results.ordered);
} # end of one.peeling.prog





#' Perform pasting one time in prognostic case.
#'
#' @param d.inside the dataset for subjects in consideration after peeling.
#' @param d.outside the dataset for subjects outside consideration after peeling.
#' @param trace.inside.condition list of signature rules used for d.inside.
#' @param alpha a parameter controlling the number of subjects in consideration.
#' @param yvar the name for response variable.
#' @param censorvar the name for censoring (1: event; 0: censor), default = NULL.
#' @param g.str ">=".
#' @param l.str "<=".
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#'
#' @return a data frame enlisting the signature rules (after pasting) ordered by main effect p-values in each group defined by the rules.
#'
one.pasting.prog = function (d.inside, d.outside, trace.inside.condition, alpha, yvar, censorvar, g.str, l.str, type, des.res) {
  ### one.pasting for prognostic case
  n.obs.inside = dim(d.inside)[1]
  n.obs.outside = dim(d.outside)[1]
  n.obs.to.add = ceiling(n.obs.inside * alpha)

  one.pasting.results = data.frame(x.nm=as.character(), condition=as.character(), cutoff=as.numeric(), pv=as.numeric(), coef=as.numeric(), end.flag=as.logical(), stringsAsFactors=F);

  pidx.pasting.results = 1;
  num.x.pasting=dim(trace.inside.condition)[1]

  for (i in 1:num.x.pasting){
    relax.condition=trace.inside.condition[i,]
    x.nm.relax.condition=relax.condition[,"x.nm"]
    other.condition=trace.inside.condition[-i,]
    idx.outside.by.other.condition=query.from.condition(d.outside, other.condition, g.str, l.str)
    d.outside.by.other.condition=d.outside[idx.outside.by.other.condition,]

    if(dim(d.outside.by.other.condition)[1] == 0) next

    if (relax.condition[,"condition"]==g.str){
      d.outside.by.other.condition.order = order(d.outside.by.other.condition[,x.nm.relax.condition],decreasing=T)
      d.outside.by.other.condition.ordered = d.outside.by.other.condition[d.outside.by.other.condition.order,]
      x.cutoff.relax.condition = d.outside.by.other.condition.ordered[,x.nm.relax.condition][min(n.obs.to.add,dim(d.outside.by.other.condition.ordered)[1])]
      n.obs.actually.add = sum(d.outside.by.other.condition.ordered[,x.nm.relax.condition]>= x.cutoff.relax.condition)

      end.flag = n.obs.actually.add>=dim(d.outside.by.other.condition.ordered)[1]

      d.inside.extended=rbind(d.inside, d.outside.by.other.condition.ordered[1:n.obs.actually.add,])
      d.outside.shrinked=rbind(d.outside[!idx.outside.by.other.condition,],d.outside.by.other.condition.ordered[-(1:n.obs.actually.add),])
      paste.x.extended.cutoff = d.outside.by.other.condition.ordered[n.obs.actually.add,x.nm.relax.condition]
      cutoff.inside.condition.str = g.str;
    }

    if (relax.condition[,"condition"]==l.str){
      d.outside.by.other.condition.order = order(d.outside.by.other.condition[,x.nm.relax.condition],decreasing=F)
      d.outside.by.other.condition.ordered = d.outside.by.other.condition[d.outside.by.other.condition.order,]
      x.cutoff.relax.condition = d.outside.by.other.condition.ordered[,x.nm.relax.condition][min(n.obs.to.add,dim(d.outside.by.other.condition.ordered)[1])]
      n.obs.actually.add = sum(d.outside.by.other.condition.ordered[,x.nm.relax.condition]<= x.cutoff.relax.condition)

      end.flag = n.obs.actually.add>=dim(d.outside.by.other.condition.ordered)[1]

      d.inside.extended=rbind(d.inside, d.outside.by.other.condition.ordered[1:n.obs.actually.add,])
      d.outside.shrinked=rbind(d.outside[!idx.outside.by.other.condition,],d.outside.by.other.condition.ordered[-(1:n.obs.actually.add),])
      paste.x.extended.cutoff = d.outside.by.other.condition.ordered[n.obs.actually.add,x.nm.relax.condition]
      cutoff.inside.condition.str = l.str;
    }

    d.complete<-rbind(d.inside.extended,d.outside.shrinked)
    grp.id<-c(rep(TRUE,nrow(d.inside.extended)),rep(FALSE,nrow(d.outside.shrinked)))

    res = try(pval.cal.prog(d.complete, yvar=yvar, censorvar=censorvar, grp.id=grp.id, type=type, des.res=des.res), silent=T)

    # skip collecting result if not comparable
    if(is.na(res$pv) ||class(res)=="try-error") next

    # save the result of pasting of x on its one side
    one.pasting.results[pidx.pasting.results, "x.nm"] = x.nm.relax.condition
    one.pasting.results[pidx.pasting.results, "condition"] = cutoff.inside.condition.str
    one.pasting.results[pidx.pasting.results, "cutoff"] = paste.x.extended.cutoff;
    one.pasting.results[pidx.pasting.results, "pv"] = res$pv
    one.pasting.results[pidx.pasting.results, "coef"] = res$coef.grp
    one.pasting.results[pidx.pasting.results, "end.flag"] = end.flag
    pidx.pasting.results = pidx.pasting.results + 1;

  }

  if (dim(one.pasting.results)[1]==0) return (one.pasting.results)

  d.current<-rbind(d.inside,d.outside)
  grp.id.current<-c(rep(TRUE,nrow(d.inside)),rep(FALSE,nrow(d.outside)))
  pv.current = pval.cal.prog(d.current, yvar=yvar, censorvar=censorvar, grp.id=grp.id.current, type=type, des.res=des.res)$pv
  idx.obs.qualified = one.pasting.results$pv < pv.current
  one.pasting.results=one.pasting.results[idx.obs.qualified,]
  order.by.pv = order(one.pasting.results$pv, decreasing=F)
  one.pasting.results.ordered = one.pasting.results[order.by.pv, ]
  return(one.pasting.results.ordered)

}




#' Perform dropping one time in prognostic case.
#'
#' @param d.inside the dataset for subjects in consideration after pasting.
#' @param d.outside the dataset for subjects outside consideration after pasting.
#' @param trace.inside.condition list of signature rules used for d.inside.
#' @param yvar the name for response variable.
#' @param censorvar the name for censoring (1: event; 0: censor), default = NULL.
#' @param g.str ">=".
#' @param l.str "<=".
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary.
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response.
#'
#' @return a data frame enlisting the signature rules (after dropping) ordered by main effect p-values in each group defined by the rules.
#'
one.dropping.prog = function (d.inside, d.outside, trace.inside.condition, yvar, censorvar, g.str, l.str, type, des.res){
  ### one.dropping for prognostic case
  one.dropping.results = data.frame(x.nm=as.character(), condition=as.character(), cutoff=as.numeric(), pv=as.numeric(), coef=as.numeric(), end.flag=as.logical(),stringsAsFactors=F);
  pidx.dropping.results = 1;
  num.x.dropping=dim(trace.inside.condition)[1]

  if (num.x.dropping<=1) return (one.dropping.results)

  for (i in 1:num.x.dropping){
    drop.condition=trace.inside.condition[i,]
    x.nm.drop.condition=drop.condition[,"x.nm"]
    other.condition=trace.inside.condition[-i,]
    idx.outside.by.other.condition=query.from.condition(d.outside, other.condition, g.str, l.str)
    d.outside.by.other.condition=d.outside[idx.outside.by.other.condition,]

    if(dim(d.outside.by.other.condition)[1] == 0) next

    if (drop.condition[,"condition"]==g.str){
      d.inside.extended = rbind(d.inside, d.outside.by.other.condition)
      d.outside.shrinked = d.outside[!idx.outside.by.other.condition,]
      x.extended.cutoff = min(d.inside.extended[,x.nm.drop.condition])
      cutoff.inside.condition.str = g.str;
    }

    if (drop.condition[,"condition"]==l.str){
      d.inside.extended=rbind(d.inside, d.outside.by.other.condition)
      d.outside.shrinked=d.outside[!idx.outside.by.other.condition,]
      x.extended.cutoff = max(d.inside.extended[,x.nm.drop.condition])
      cutoff.inside.condition.str = l.str;
    }

    d.complete<-rbind(d.inside.extended,d.outside.shrinked)
    grp.id<-c(rep(TRUE,nrow(d.inside.extended)),rep(FALSE,nrow(d.outside.shrinked)))

    res = try(pval.cal.prog(d.complete, yvar=yvar, censorvar=censorvar, grp.id=grp.id, type=type, des.res=des.res),silent=T)

    # skip collecting result if not comparable
    if(is.na(res$pv)|| class(res)=="try-error") next

    # save the result of pasting of x on its one side
    one.dropping.results[pidx.dropping.results, "x.nm"] = x.nm.drop.condition
    one.dropping.results[pidx.dropping.results, "condition"] = cutoff.inside.condition.str
    one.dropping.results[pidx.dropping.results, "cutoff"] = x.extended.cutoff;
    one.dropping.results[pidx.dropping.results, "pv"] = res$pv
    one.dropping.results[pidx.dropping.results, "coef"] = res$coef.grp
    one.dropping.results[pidx.dropping.results, "end.flag"] = TRUE
    pidx.dropping.results = pidx.dropping.results + 1;

  }

  if (dim(one.dropping.results)[1]==0) return (one.dropping.results)

  order.by.pv = order(one.dropping.results$pv, decreasing=F)
  one.dropping.results.ordered = one.dropping.results[order.by.pv, ]
  return(one.dropping.results.ordered)

}


#' Calculate p-value for treatment in each subgroup in prognostic case
#'
#' @param data input data frame
#' @param yvar name for response variable
#' @param censorvar name for censoring (1: event; 0: censor), default = NULL
#' @param grp.id subgroup id
#' @param type type of response variable: "c" continuous (default); "s" survival; "b" binary
#' @param des.res the desired response. "larger": prefer larger response (default) "smaller": prefer smaller response
#'
#' @return p-value for the main effect given the dataset
#'
pval.cal.prog <- function(data, yvar, censorvar, grp.id, type, des.res) {
  # avoid no observation in either group
  if (type=="s"){
    fit.formula <- as.formula(paste("Surv(", yvar, ", ", censorvar, ")~grp.id", sep=""))
    fit.obj <- coxph(fit.formula, data=data)
    pv = summary(fit.obj)$coefficients[5]/2
    coef.grp=summary(fit.obj)$coefficients[1]

    if (is.na(pv) || is.na(coef.grp)) # avoid NA case for pv and coefficient
    {
      pv <- 1
      coef.grp<-0 # 0 is arbitrary since anyway pv = 1 (pv is greater than every other p-value)
    }else
    {
      if (des.res=="larger" & coef.grp>0) pv=1-pv
      if (des.res=="smaller" & coef.grp<0) pv=1-pv
    }

  }

  fit.formula = as.formula(paste(yvar, "~grp.id" , sep=""))
  if (type=="b"){
    fit.obj=glm(fit.formula,family=binomial,data=data)
    if (nrow(summary(fit.obj)$coefficients)<2)
    {
      pv <- 1
      coef.grp<-0 # 0 is arbitrary since anyway pv = 1 (pv is greater than every other p-value)
    }else{
      pv=summary(fit.obj)$coefficients[2,4]/2
      coef.grp=summary(fit.obj)$coefficients[2,1]
      if (des.res=="larger" & coef.grp<0) pv=1-pv
      if (des.res=="smaller" & coef.grp>0) pv=1-pv

    }
  }

  if (type=="c"){
    fit.obj=lm(fit.formula,data=data)
    if (nrow(summary(fit.obj)$coefficients)<2)
    {
      pv <- 1
      coef.grp<-0 # 0 is arbitrary since anyway pv = 1 (pv is greater than every other p-value)
    }else{
      pv=summary(fit.obj)$coefficients[2,4]/2
      coef.grp=summary(fit.obj)$coefficients[2,1]
      if (des.res=="larger" & coef.grp<0) pv=1-pv
      if (des.res=="smaller" & coef.grp>0) pv=1-pv

    }
  }


  list(pv=pv, coef.grp=coef.grp)

}

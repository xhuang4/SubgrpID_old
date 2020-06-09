##################################################################################

######### !!! source required files before running the algorithms !!!##############

##################################################################################

SubgrpID = function(data.train, data.test=NULL, yvar, censorvar=NULL, trtvar=NULL, trtref=NULL,
                    xvars, type="c",n.boot=ifelse(method=="PRIM",0,25), des.res="larger",
                    min.sigp.prcnt=0.20, pre.filter=NULL, filter.method=NULL,k.fold=5, cv.iter=20,
                    max.iter=500, mc.iter=20, method=c("AIM.Rule"), train.percent.prim=0.5,
                    do.cv=FALSE, out.file=NULL, file.path="", plots=F, ...)
{

  #data.train: training dataframe
  #data.test: testing dataframe
  #yvar: response variable name
  #cencorvar: censoring variable name 1:event; 0: censor.
  #trtvar: treatment variable name
  #trtref: code for treatment arm
  #xvars: covariates variable name
  #type: "c" continuous; "s" survival; "b" binary
  #n.boot: number of bootstrap for batting procedure, or the variable selection procedure for PRIM
  #        For PRIM, when n.boot=0, bootstrapping for variable selection is not conducted
  #des.res: the desired response. "larger": prefer larger response. "smaller": prefer smaller response
  # min.sigp.prcnt: desired proportion of signature positive group size for a given cutoff.
  #pre.filter: NULL, no prefiltering conducted;"opt", optimized number of predictors selected; An integer: min(opt, integer) of predictors selected
  #filter.method: NULL, no prefiltering, "univariate", univaraite filtering; "glmnet", glmnet filtering, "unicart": univariate rpart filtering for prognostic case
  # k.fold: # cross-validation folds
  # cv.iter: Algotithm terminates after cv.iter successful iterations of cross-validation, or after
  # max.iter total iterations, whichever occurs first.
  # mc.iter: # of iterations for the MC procedure to get a stable "best number of predictors"
  #method: "AIM", "AIM.Rule", "Seq.BT", "PRIM", "SIDES", "optAUC"
  #train.percent.prim: IF train.percent.prim=1, all data will be used both for sub-training and sub-testing inside the PRIM
  #do.cv: whether to do cross validation or not. TRUE or FALSE
  #out.file: Name of output result files. If NULL no output file would be saved.
  #file.path: default current working directory. When specifying a dir, use "/" at the end. e.g. "TEMP/"
  #plots: whether to save plots


  # OptAUC arguments:
  # method.subgrp: method for subgroup identification in OptAUC: "BATTing" or "youden".
  # scale: logical variable indicating if scaling is required.
  # xselect: it is a logical flag. TRUE: use LASSO for variable selection; FALSE: use all predictors to calculate the composite score. Default = FALSE.
  # cv.iter.xselect: it is the number of iterations of cross-validation used for xselect.
  # k.fold.xselect: it is the number of folds of cross-validation used for xselect.
  # method.xselect: it is the method used for variable selection: CV (default) cross validated AUC/ABC; "aAUC" approximate AUC (only for prognostic case).

  if(is.null(censorvar)) {
    censor.vec.train=NULL
    censor.vec.test=NULL
  }else{
    censor.vec.train=data.train[censorvar]
    censor.vec.test=data.test[censorvar]
  }
  if(is.null(trtvar)) {
    trt.vec.train=NULL
    trt.vec.test=NULL
  }else{
    trt.vec.train=data.train[trtvar]
    trt.vec.test=data.test[trtvar]
  }


  if(method=="Seq.BT"){
    res=seqlr.batting(y=data.train[yvar], x=data.train[xvars], censor.vec=censor.vec.train,
                      trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
                      des.res=des.res, min.sigp.prcnt=min.sigp.prcnt,
                      pre.filter=pre.filter, filter.method=filter.method)
    pred.data=pred.seqlr(data.train[xvars], res)
    train.stat=evaluate.results(data.train[yvar], pred.data, censor.vec=censor.vec.train,
                                trt.vec=trt.vec.train, trtref=trtref, type=type)
    if (!is.null(data.test)){
      pred.data=pred.seqlr(data.test[xvars], res)
      test.stat=evaluate.results(data.test[yvar], pred.data, censor.vec=censor.vec.test,
                                 trt.vec=trt.vec.test, trtref=trtref, type=type)
    }

    if (do.cv){
      cv.res=cv.seqlr.batting(y=data.train[yvar],x=data.train[xvars], censor.vec=censor.vec.train,
                              trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
                              des.res=des.res,  min.sigp.prcnt=min.sigp.prcnt,
                              pre.filter=pre.filter, filter.method=filter.method,
                              k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)
    }


    if (!is.null(trtvar)){
      train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (Seq.BT, Train)", trt.lab=c("Trt.", "Ctrl."))
      if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (Seq.BT, Test)", trt.lab=c("Trt.", "Ctrl."))
    }
  }

  #############################################################################################
  if (method == "PRIM") {
    res = prim.train(data=data.train, yvar=yvar, censorvar=censorvar, trtvar=trtvar,
                     trtref=trtref, xvars=xvars, type=type, des.res=des.res,
                     alpha = c(0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50),
                     min.sigp.prcnt=min.sigp.prcnt,  training.percent = train.percent.prim,
                     n.boot=n.boot, pre.filter=pre.filter, filter.method=filter.method)
    pred.data=pred.prim(data.train[xvars], res)
    train.stat=evaluate.results(data.train[yvar], pred.data, censor.vec=censor.vec.train,
                                trt.vec=trt.vec.train, trtref=trtref, type=type)

    if (!is.null(data.test)){
      pred.data=pred.prim(data.test[xvars], res)
      test.stat=evaluate.results(data.test[yvar], pred.data, censor.vec=censor.vec.test,
                                 trt.vec=trt.vec.test, trtref=trtref, type=type)
    }

    if (do.cv){
      cv.res=prim.cv(data=data.train, yvar=yvar, censorvar=censorvar, trtvar=trtvar,
                     trtref=trtref, xvars=xvars, type=type, des.res=des.res,
                     alpha = c(0.05,0.06,0.07,0.08,0.09,0.10,0.20,0.30,0.40,0.50),
                     min.sigp.prcnt=min.sigp.prcnt,  training.percent = train.percent.prim,
                     n.boot=n.boot, pre.filter=pre.filter, filter.method=filter.method,
                     k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)

    }

    if (!is.null(trtvar)){
      train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (PRIM, Train)", trt.lab=c("Trt.", "Ctrl."))
      if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (PRIM, Test)", trt.lab=c("Trt.", "Ctrl."))
    }

  }


  #############################################################################################
  if(method=="AIM"){
    #        library(AIM)
    res=aim.batting(y=data.train[yvar],x=data.train[xvars], censor.vec=censor.vec.train,
                    trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
                    des.res=des.res, min.sigp.prcnt=min.sigp.prcnt, mc.iter=mc.iter,
                    pre.filter=pre.filter, filter.method=filter.method)
    pred.data=pred.aim(data.train[xvars], res)
    train.stat=evaluate.results(data.train[yvar], pred.data, censor.vec=censor.vec.train,
                                trt.vec=trt.vec.train, trtref=trtref, type=type)

    if (!is.null(data.test)){
      pred.data=pred.aim(data.test[xvars], res)
      test.stat=evaluate.results(data.test[yvar], pred.data, censor.vec=censor.vec.test,
                                 trt.vec=trt.vec.test, trtref=trtref, type=type)
    }
    if (do.cv){
      cv.res=cv.aim.batting(y=data.train[yvar],x=data.train[xvars], censor.vec=censor.vec.train,
                            trt.vec=trt.vec.train, trtref=trtref, type=type,
                            n.boot=n.boot, des.res=des.res, min.sigp.prcnt=min.sigp.prcnt,
                            mc.iter=mc.iter, pre.filter=pre.filter, filter.method=filter.method,
                            k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)

    }

    if (!is.null(trtvar)){
      train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (AIM, Train)", trt.lab=c("Trt.", "Ctrl."))
      if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (AIM, Test)", trt.lab=c("Trt.", "Ctrl."))
    }


  }


  #############################################################################################
  if(method=="AIM.Rule"){
    #        library(AIM)
    res=aim.rule.batting(y=data.train[yvar], x=data.train[xvars], censor.vec=censor.vec.train,
                         trt.vec=trt.vec.train, trtref=trtref, type=type, n.boot=n.boot,
                         des.res=des.res, min.sigp.prcnt=min.sigp.prcnt, mc.iter=mc.iter,
                         pre.filter=pre.filter, filter.method=filter.method)

    pred.data=pred.aim(data.train[xvars], res)
    train.stat=evaluate.results(data.train[yvar], pred.data, censor.vec=censor.vec.train,
                                trt.vec=trt.vec.train, trtref=trtref, type=type)

    if (!is.null(data.test)){
      pred.data=pred.aim(data.test[xvars], res)
      test.stat=evaluate.results(data.test[yvar], pred.data, censor.vec=censor.vec.test,
                                 trt.vec=trt.vec.test, trtref=trtref, type=type)

    }

    if (do.cv){
      cv.res=cv.aim.rule.batting(y=data.train[yvar],x=data.train[xvars], censor.vec=censor.vec.train,
                                 trt.vec=trt.vec.train, trtref=trtref, type=type,
                                 n.boot=n.boot, des.res=des.res, min.sigp.prcnt=min.sigp.prcnt,
                                 mc.iter=mc.iter, pre.filter=pre.filter, filter.method=filter.method,
                                 k.fold=k.fold, cv.iter=cv.iter, max.iter=max.iter)

    }

    if (!is.null(trtvar)){
      train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (AIM.Rule, Train)", trt.lab=c("Trt.", "Ctrl."))
      if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (AIM.Rule, Test)", trt.lab=c("Trt.", "Ctrl."))
    }


  }

  ########################################### SIDES #################################################

  if(method=="SIDES"){

    if(is.null(trtvar)) stop("SIDES only works for predictive case for now.")

    data("sidak_adj")

    # p1=function(ZE1,ZE2,ZS1=NULL,ZS2=NULL)
    # {
    #   2*(1-pnorm(abs(ZE1-ZE2)/sqrt(2)))
    # }
    #
    # p2=function(ZE1,ZE2,ZS1=NULL,ZS2=NULL)
    # {
    #   2*min(1-pnorm(ZE1),1-pnorm(ZE2))
    # }
    #
    # p3=function(ZE1,ZE2,ZS1,ZS2)
    # {
    #   max(p1(ZE1,ZE2),p2(ZE1,ZE2))
    # }
    #
    # p4=function(ZE1,ZE2,ZS1,ZS2)
    # {
    #   w=0.5 # a weight parameter balance efficacy and risk
    #
    #   2*(1-pnorm(abs(w*(ZE1-ZE2)+(1-w)*(ZS2-ZS1))/sqrt(2*(w^2+(1+w)^2))))
    # }

    M=10; split.score=p1;  maxncov=3; gamma=rep(0.5,length(maxncov));
    minSplit.prcnt=0.10; minNode=3; saf.test=NULL; maxCand=NULL;
    pm.level=0.1; pm.repeat=100; ord.bin=10;

    SIDES.args<-list(...)

    if(!is.null(SIDES.args$M)) M<-SIDES.args$M
    if(!is.null(SIDES.args$split.score)) split.score<-get(SIDES.args$split.score)
    if(!is.null(SIDES.args$maxncov)) maxncov<-SIDES.args$maxncov
    if(!is.null(SIDES.args$gamma)) gamma<-SIDES.args$gamma
    if(!is.null(SIDES.args$minSplit.prcnt)) minSplit.prcnt<-SIDES.args$minSplit.prcnt
    if(!is.null(SIDES.args$minNode)) minNode<-SIDES.args$minNode
    if(!is.null(SIDES.args$saf.test)) saf.test<-SIDES.args$saf.test
    if(!is.null(SIDES.args$maxCand)) maxCand<-SIDES.args$maxCand
    if(!is.null(SIDES.args$pm.level)) pm.level<-SIDES.args$pm.level
    if(!is.null(SIDES.args$pm.repeat)) pm.repeat<-SIDES.args$pm.repeat
    if(!is.null(SIDES.args$ord.bin)) ord.bin<-SIDES.args$ord.bin


    # # ================================================================
    # # ONE-SIDED TEST FUNCTION FOR TRT EFFECT/SAFETY
    # # ================================================================
    #
    # # Linear Regression for Continuous Response
    # lm.test=function(y,trt,id=NULL,uppertail=TRUE)
    # {
    #   if(length(unique(y))==1)
    #     return(list(t=0,p=1))
    #
    #   temp=summary(lm(y~trt))
    #   if(nrow(temp$coef)<2)
    #     return(list(t=0,p=1))
    #
    #   t=temp$coef[2,3]
    #   p=1-pt(t,temp$df[2])
    #
    #   if(!uppertail)
    #   {
    #     t=-t
    #     p=1-p
    #   }
    #   return(list(t=t,p=p))
    # }
    #
    # # Logistic Regression for Binary Response
    # logistic.test=function(y,trt,id=NULL,uppertail=TRUE)
    # {
    #   if(length(unique(y))==1)
    #     return(list(t=0,p=1))
    #   if(all(y==trt))
    #     return(list(t=Inf,p=0))
    #
    #   temp=summary(glm(y~trt,family=binomial))
    #   if(nrow(temp$coef)<2)
    #     return(list(t=0,p=1))
    #
    #   t=temp$coef[2,3]
    #   p=1-pt(t,temp$df[2])
    #   if(!uppertail)
    #   {
    #     t=-t
    #     p=1-p
    #   }
    #   return(list(t=t,p=p))
    # }
    #
    # logrank.test=function(y,trt,id=NULL,uppertail=TRUE)
    # {
    #   t=summary(coxph(y~trt))$coef[4]
    #
    #   if(is.na(t))
    #     return(list(t=0,p=1))
    #
    #   p=1-pnorm(t)
    #   if(!uppertail)
    #   {
    #     t=-t
    #     p=1-p
    #   }
    #   return(list(t=t,p=p))
    # }
    #


    if (type == "c") {
      eff.test=lm.test
    } else if ((type == "b")) {
      eff.test=logistic.test
    } else {
      eff.test=logrank.test
    }

    if (des.res == "larger") {
      uppertail = TRUE
    } else {
      uppertail = FALSE
    }


    res <- SIDES.tree(data=data.train, yvar=yvar,censorvar=censorvar, trtvar=trtvar, trtref=trtref, xvars=xvars, type=type, M=M, split.score=split.score, gamma=gamma, eff.test=eff.test, maxncov=maxncov, minSplit.prcnt=minSplit.prcnt, minNode=minNode, saf.test=saf.test, maxCand=maxCand, pm.level=pm.level, pm.repeat=pm.repeat, uppertail=uppertail, ord.bin=ord.bin)

    pred.data <- predict.sides(data.train, xvars, res)
    train.stat=evaluate.results(y=data.train[,yvar], predict.data=pred.data, censor.vec=censor.vec.train, trt.vec=trt.vec.train, trtref=trtref, type=type)

    if (!is.null(data.test)){
      pred.data=predict.sides(data.test, xvars, res)
      if (class(try(pred.data[,"pred.class"],silent = T))=="try-error") stop("no signature is found in test data")
      test.stat=evaluate.results(y=data.test[,yvar], predict.data=pred.data, censor.vec=censor.vec.test,
                                 trt.vec=trt.vec.test, trtref=trtref, type=type)
    }

    if (do.cv){

      cv.res=cv.SIDES(data=data.train, yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtref =trtref, xvars=xvars, type=type, M=M, split.score=split.score,
                      gamma=gamma, eff.test=eff.test, saf.test=saf.test, maxncov=maxncov,
                      minSplit.prcnt=minSplit.prcnt, minNode=minNode, maxCand=maxCand, pm.level=pm.level,                        pm.repeat=pm.repeat, uppertail=uppertail, ord.bin=ord.bin, k.fold=k.fold,
                      cv.iter=cv.iter, max.iter=max.iter)



    }


    if (!is.null(trtvar)){
      train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (SIDES, Train)", trt.lab=c("Trt.", "Ctrl."))
      if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (SIDES, Test)", trt.lab=c("Trt.", "Ctrl."))
    }

    best <- bestnode(res)$best
    sig <- get.sig(best)
    sig.df <- data.frame()
    for (i in 1:length(sig$var)) {
      if (sig$side[i]) {
        dir.tmp <- "<="
      } else {
        dir.tmp <- ">"
      }
      sig.tmp <- data.frame(var=xvars[sig$var[i]], dir=dir.tmp, thresh=sig$thresh[i])
      sig.df <- rbind(sig.df, sig.tmp)
    }

    res<-sig.df
  }

  ########################################### optAUC #################################################

  if(method=="optAUC"){
    # library(optAUCX)
    yvar.ref=1; method.subgrp="BATTing"; scale=TRUE; xselect=TRUE; cv.iter.xselect=20; k.fold.xselect=5; method.xselect="CV"

    args<-list(...)
    if(!is.null(args$yvar.ref)) yvar.ref<-args$yvar.ref
    if(!is.null(args$scale)) scale<-args$scale
    if(!is.null(args$xselect)) xselect<-args$xselect
    if(!is.null(args$cv.iter.xselect)) cv.iter.xselect<-args$cv.iter.xselect
    if(!is.null(args$k.fold.xselect)) k.fold.xselect<-args$k.fold.xselect
    if(!is.null(args$method.xselect)) method.xselect<-args$method.xselect
    if(!is.null(args$method.OptAUC)) method.subgrp<-args$method.OptAUC

    res<-optAUCX.fit(data=data.train, yvar=yvar, yvar.ref=yvar.ref, xvars=xvars, censorvar=censorvar, trtvar=trtvar, type=type, trtref=trtref, des.res=des.res, scale=scale, xselect=xselect, cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect, pre.filter=pre.filter, filter.method=filter.method, n.boot=n.boot, method.subgrp=method.subgrp)
    res<-res$subgrp$res
    pred.data<-pred.optAUCX(data.train[xvars],res)
    train.stat=evaluate.results(data.train[yvar], pred.data, censor.vec=censor.vec.train, trt.vec=trt.vec.train, trtref=trtref, type=type)
    if (!is.null(data.test)){
      pred.data<-pred.optAUCX(data.test[xvars],res)
      test.stat=evaluate.results(data.test[yvar], pred.data, censor.vec=censor.vec.test, trt.vec=trt.vec.test, trtref=trtref, type=type)
    }

    if (do.cv){
      cv.res=cv.optAUCX.fit(data=data.train, yvar=yvar, yvar.ref=yvar.ref, xvars=xvars, censorvar=censorvar, trtvar=trtvar, type=type, trtref=trtref, des.res=des.res, method.subgrp=method.subgrp, scale=scale, xselect=xselect, cv.iter.xselect=cv.iter.xselect, k.fold.xselect=k.fold.xselect, method.xselect=method.xselect, pre.filter=pre.filter, filter.method=filter.method, n.boot=n.boot, cv.iter.eval=cv.iter, k.fold.eval=k.fold, max.iter=max.iter)
    }

    if (!is.null(trtvar)){
      train.plot=interaction.plot(train.stat, type=type, main="Interaction Plot (OptAUC, Train)", trt.lab=c("Trt.", "Ctrl."))
      if (!is.null(data.test)) test.plot=interaction.plot(test.stat, type=type, main="Interaction Plot (OptAUC, Test)", trt.lab=c("Trt.", "Ctrl."))
    }

  }



  ############################ output files ########################################
  if (!is.null(out.file)){

    if (method=="AIM") {
      res.temp=res$aim.rule
      res.temp=rbind(res.temp,data.frame(variable="Score", direction=">", cutoff=res$bat.cutoff))
      res=res.temp
      row.names(res)=NULL
    }

    if (method=="AIM.Rule") {
      res=res$aim.rule
      row.names(res)=NULL
    }


    res.file=paste(file.path,out.file,".",method,".res.csv",sep="")
    train.pval.file=paste(file.path,out.file,".",method,".train.pval.csv",sep="")
    train.ratios.file=paste(file.path,out.file,".",method,".train.ratios.csv",sep="")
    train.gpstats.file=paste(file.path,out.file,".",method,".train.group.stats.csv",sep="")
    test.pval.file=paste(file.path,out.file,".",method,".test.pval.csv",sep="")
    test.ratios.file=paste(file.path,out.file,".",method,".test.ratios.csv",sep="")
    test.gpstats.file=paste(file.path,out.file,".",method,".test.group.stats.csv",sep="")
    cv.pval.file=paste(file.path,out.file,".",method,".cv.pval.csv",sep="")
    cv.ratios.file=paste(file.path,out.file,".",method,".cv.ratios.csv",sep="")
    cv.gpstats.file=paste(file.path,out.file,".",method,".cv.group.stats.csv",sep="")

    write.csv(res,file=res.file, row.names=FALSE)
    write.csv(train.stat$pval,file=train.pval.file, row.names=FALSE)
    if(type!="c") write.csv(train.stat$ratios,file=train.ratios.file, row.names=FALSE)
    write.csv(train.stat$group.stats,file=train.gpstats.file)

    if (!is.null(data.test)){
      write.csv(test.stat$pval,file=test.pval.file, row.names=FALSE)
      if(type!="c") write.csv(test.stat$ratios,file=test.ratios.file, row.names=FALSE)
      write.csv(test.stat$group.stats,file=test.gpstats.file)
    }

    if (do.cv){
      write.csv(cv.res$stats.summary$pvals,file=cv.pval.file, row.names=FALSE)
      if(type!="c") write.csv(cv.res$stats.summary$ratios,file=cv.ratios.file, row.names=FALSE)
      write.csv(cv.res$stats.summary$group.stats,file=cv.gpstats.file)
    }

    if (plots)
    {
      if (!is.null(trtvar))
      {
        train.plot.file=paste(file.path,out.file,".",method,".train.plot.jpg",sep="")
        ggsave(filename=train.plot.file,plot=train.plot, width=4, height=4,dpi=300)
      }

      if (!is.null(data.test))
      {
        if (!is.null(trtvar))
        {
          test.plot.file=paste(file.path,out.file,".",method,".test.plot.jpg",sep="")
          ggsave(filename=test.plot.file,plot=test.plot, width=4, height=4,dpi=300)
        }
      }

      if (do.cv)
      {
        if (!is.null(trtvar))
        {
          cv.plot.file=paste(file.path, out.file,".",method,".cv.plot.jpg",sep="")
          ggsave(filename=cv.plot.file,plot=cv.res$interplot, width=4, height=4,dpi=300)
        }
      }
    }
  }


  ###################################################
  if (is.null(trtvar)) {
    train.plot=NULL
    test.plot=NULL
  }

  if (is.null(data.test)){
    test.stat=NULL
    test.plot=NULL
  }

  if (!do.cv){
    cv.res=NULL
  }

    result=list(res=res, train.stat=train.stat, test.stat=test.stat, cv.res=cv.res,
                train.plot=train.plot, test.plot=test.plot)

  if (!is.null(out.file)){
    save(data.train, data.test, result, file = paste(file.path,out.file,".",method,".RData",sep=""))
  }

  result
}






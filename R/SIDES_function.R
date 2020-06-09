# ================================================================
# CLASS DEFINITION OF SUBGROUP OBJECTS
# ================================================================

# A single split in SIDES
SIDEScut = setClass("SIDEScut", representation(cov="numeric",factor="logical",
		side="logical", value="character"))

# A subgroup in SIDES
SIDESsub = setClass("SIDESsub",
	representation(data="data.frame",sub="vector",level="numeric",
		signif="numeric",used="vector",cut="list",n="numeric"))

# A node in recursive partition process
SIDESnode=setClass("SIDESnode", representation(node="SIDESsub",children="list",m="numeric"))

# ================================================================
# UTILITY FUNCTION LISTING POSSIBLE CUTS IN A SPLIT
# ================================================================

# Norminal variable
power.set <- function(x) {
        if(length(x) == 0)
                return(vector(mode(x), 0))
        x <- as.character(unique(x)); n <- length(x); K <- NULL
        for(m in x) K <- rbind(cbind(K, F), cbind(K, T))
        out <- apply(K, 1, function(x, s) s[x], s = x)
         out <- out[-c(1, length(out))]
         l <- length(out); i <- 1
         while (i <= l) {
            if (length(out[[i]]) >= ceiling(n/2+.5)) {out[[i]] <- NULL; i <- i-1}
		if(length(out[[i]])==n/2)
			for(j in 1:i)
				if(all(is.element(x,c(out[[i]],out[[j]]))))
					{out[[i]] <- NULL; i <- i-1; break}
            i <- i+1 
            l <- length(out) 
        }
        out=sapply(out,factor,levels=x)
	out
}      

# Continuous/Ordinal variable
ord.cut = function(x,ord.bin)
{
	x=x[!is.na(x)]
	temp=sort(unique(x))
	temp = temp[-length(temp)]
	if(length(temp) < ord.bin)
		return(list(cut=temp,g=length(temp)))

	ordseq=round(seq(0,length(x),length(x)/ord.bin))
	ordseq=ordseq[c(-1,-11)]
	temp=unique(sort(x)[ordseq])
	return(list(cut=temp,g=length(temp)))
}

# ================================================================
# ONE-SIDED TEST FUNCTION FOR TRT EFFECT/SAFETY
# ================================================================

# Linear Regression for Continuous Response
lm.test=function(y,trt,id=NULL,uppertail=TRUE)
{
	if(length(unique(y))==1)
		return(list(t=0,p=1))

	temp=summary(lm(y~trt))
	if(nrow(temp$coef)<2)
		return(list(t=0,p=1))
  
	t=temp$coef[2,3]
	p=1-pt(t,temp$df[2])
  
  if(!uppertail)
  {
    t=-t
    p=1-p
  }
	return(list(t=t,p=p))
}

# Logistic Regression for Binary Response
logistic.test=function(y,trt,id=NULL,uppertail=TRUE)
{
	if(length(unique(y))==1)
		return(list(t=0,p=1))
	if(all(y==trt))
		return(list(t=Inf,p=0))

	temp=summary(glm(y~trt,family=binomial))
	if(nrow(temp$coef)<2)
		return(list(t=0,p=1))

	t=temp$coef[2,3]
	p=1-pt(t,temp$df[2])
	if(!uppertail)
	{
	  t=-t
	  p=1-p
	}
	return(list(t=t,p=p))
}

logrank.test=function(y,trt,id=NULL,uppertail=TRUE)
{
	t=summary(coxph(y~trt))$coef[4]
  
	if(is.na(t))
	  return(list(t=0,p=1))
  
	p=1-pnorm(t)
	if(!uppertail)
	{
	  t=-t
	  p=1-p
	}
	return(list(t=t,p=p))
}

# Sidak Adjust According to the Given Table
sidak.adjust=function(p,g,type)
{
	return(1-sign(1-p)*(abs(1-p))^(sidak_adj[[type+1]][g-1]))
}

# ================================================================
# GROWING A SIDES TREE
# ================================================================

# Function to Initiate a SIDES Tree
SIDES.tree = function(data, yvar, censorvar, trtvar, trtref, xvars, type, M,split.score, gamma,eff.test,saf.test=NULL,maxncov=3,minSplit.prcnt=0.1,minNode=3,maxCand=NULL, pm.level=0.1,pm.repeat=500,
  uppertail=TRUE,ord.bin=10)
{


  if (!is.null(trtref)){
    data[,trtvar]=(data[,trtvar]==trtref)+0
  }  
  
	if(length(gamma)<maxncov)
		gamma=c(gamma,rep(gamma[length(gamma)],maxncov-length(gamma)))

	minSplit<-floor(minSplit.prcnt*nrow(data))
	
	if (type == "s")
	{
	  if (!is.null(censorvar))
	  {
	    data[,yvar]=Surv(data[,yvar], data[,censorvar])
	  } else {
	    stop("Censoring variable not found")
	  }
	}
	  
	data<-cbind(data,id=1:nrow(data))
	  

	root = SIDESnode(node=SIDESsub(data=data,sub=1:nrow(data),
    signif=eff.test(data[,yvar],data[,trtvar],id=1:nrow(data),uppertail)$p,
		level=0,used=0,cut=list(),n=nrow(data)), children = vector("list",M),m=0)


	root = branch.SIDESnode(root,yvar, trtvar, xvars, M, split.score, gamma,
	                        eff.test,saf.test,
		maxncov,minSplit,minNode,maxCand,pm.level,pm.repeat,uppertail,ord.bin)
	
	if (length(root@children)==0) stop("no signature is found")	
	
	root
}

# Function to Build a SIDES Tree Recursively
branch.SIDESnode=function(root,yvar, trtvar, xvars, M, split.score, gamma,eff.test,saf.test,
		maxncov,minSplit,minNode,maxCand,pm.level,pm.repeat,uppertail,ord.bin)
{
# root contains all subjects
  

  level=root@node@level+1
	temp=SIDES.candidate(root@node,yvar, trtvar, xvars,M,split.score, eff.test, saf.test, maxncov,
		minSplit, minNode,maxCand,uppertail,ord.bin) #outputs top M no. of candidate branches
	candidate=temp$topM
	sublist=list()
	plist=NULL

	for(i in 1:M)
	{
		if(is.null(candidate[[i]]) || (candidate[[i]]@signif > gamma[level]*root@node@signif))
			next

		newnode=SIDESnode(node=candidate[[i]],children=vector("list",M),m=0)
		root@m=root@m+1
		root@children[[root@m]]=newnode
		sublist=append(sublist,list(temp$sublist[[i]]))
		plist=c(plist,temp$plist[i])
	}
	

	selectid=permute.select(root@node@data,yvar, trtvar, xvars,sublist,plist, eff.test,pm.level,pm.repeat,uppertail)
	root@m=length(selectid)
	templist=vector("list",root@m)
	if(root@m>0)
		for(i in 1:length(selectid))
		templist[[i]]=root@children[[selectid[i]]]
	root@children=templist

	if(root@m>0)
		for(i in 1:root@m)
		root@children[[i]]=branch.SIDESnode(root@children[[i]],yvar, trtvar, xvars, M,split.score, gamma,eff.test,
			saf.test,maxncov,minSplit,minNode,maxCand,pm.level,pm.repeat,uppertail,ord.bin)

	return(root)
}

# Function to List the Promissing Subset at Current Step
SIDES.candidate = function(parent,yvar, trtvar, xvars, M,split.score, eff.test, saf.test, maxncov,
minSplit, minNode,maxCand=NULL,uppertail,ord.bin)
{
  

	sub=parent@sub
	used=parent@used
	l=parent@level
	n=parent@n
	rowid=1:n

	data=parent@data
	y=data[,yvar]
	trt=data[,trtvar]
	id=data$id
	colid=(1:ncol(data[,xvars]))[!is.element(1:ncol(data[,xvars]),used)]

	if(length(used) >= maxncov)
		return(NULL)
	
	if(is.null(maxCand))
		mtry=length(colid)

	topM=vector("list",M)
	best.score=rep(ord.bin,M)
	plist=rep(0,M)
	sublist=vector("list",M)

	for(i in sample(colid,mtry,replace=F))
	{
		x=data[,xvars][,i]

		if(is.factor(x)) 
		{
			zcut = power.set(unique(x))
			temp.g = length(unique(x))
		}
		else 
		{
			temp = ord.cut(x,ord.bin)
			zcut = temp$cut
			temp.g = temp$g
		}

		if(temp.g < 2) next

		for(j in zcut)
		{
			if(is.factor(x))
			{
				left=rowid[is.element(x , j)]
				right=rowid[!is.element(x , j)]
			}
			else
			{
				left=rowid[!is.na(x) & x <= j]
				right=rowid[!is.na(x) & x > j]
			}

			if(min(length(left),length(right))<minNode)
				next

			test1=eff.test(y[left],trt[left],id[left],uppertail)
      ZE1=test1$t
			test2=eff.test(y[right],trt[right],id[right],uppertail)
      ZE2=test2$t

			if(ZE1>ZE2)
			{
				temp.side=T
				temp.id=left
				minarm=min(sum(trt[left]),sum(1-trt[left]))
				ncand=length(left)
				temp.p=test1$p
			}
			else
			{
				temp.side=F
				temp.id=right
				minarm=min(sum(trt[right]),sum(1-trt[right]))
				ncand=length(right)
				temp.p=test2$p
			}
			if(minarm<minSplit)
				next

			if(!is.null(saf.test))
			{
				ZS1=saf.test(y[left],trt[left],id[left])
				ZS2=saf.test(y[right],trt[right],id[right])
			}
			else
				ZS1=ZS2=NA

			temp.score = sidak.adjust(split.score(ZE1,ZE2,ZS1,ZS2),temp.g,is.factor(x))
			if(temp.score < max(best.score))
			{
				k=which.max(best.score)
				best.score[k] = temp.score
				sublist[[k]] = rowid[temp.id]
				plist[k] = temp.p
				temp.cut=list(SIDEScut(cov=i,factor=is.factor(x),side=temp.side,
					value=as.character(j)))
				topM[[k]] = SIDESsub(data=data[temp.id,],sub=sub[temp.id],signif=temp.p,
					level=l+1,used=c(used,i),cut=append(parent@cut,temp.cut),n=ncand)
			}
		}
	}

	return(list(topM=topM,sublist=sublist,plist=plist))
}

# Score Functions

p1=function(ZE1,ZE2,ZS1=NULL,ZS2=NULL)
{
	2*(1-pnorm(abs(ZE1-ZE2)/sqrt(2)))
}

p2=function(ZE1,ZE2,ZS1=NULL,ZS2=NULL)
{
	2*min(1-pnorm(ZE1),1-pnorm(ZE2))
}

p3=function(ZE1,ZE2,ZS1,ZS2)
{
	max(p1(ZE1,ZE2),p2(ZE1,ZE2))
}

p4=function(ZE1,ZE2,ZS1,ZS2)
{
	w=0.5 # a weight parameter balance efficacy and risk

	2*(1-pnorm(abs(w*(ZE1-ZE2)+(1-w)*(ZS2-ZS1))/sqrt(2*(w^2+(1+w)^2))))
}

# ================================================================
# FUNCTIONS FOR SUMMARY
# ================================================================

# Print all Candidate Subgroups
print.all=function(root,cutoff=NULL,rdigit=4)
{
	print.node(root@node,cutoff,rdigit=rdigit)
	if(root@m>0)
		for(i in 1:root@m)
			print.all(root@children[[i]],cutoff,rdigit)
}

# Print a Single Subgroup: Description, Size, P-Value
print.node=function(node,cutoff=NULL,out=T,rdigit=4)
{
	if(!is.null(cutoff) && node@signif>cutoff)
		return(NULL)

	text=NULL
	
	if(node@level==0)
		text=paste(text,"Full data. ")
	else
	for(i in 1:length(node@cut))
	{
		text=paste(text,colnames(node@data$x)[node@cut[[i]]@cov])
		if(node@cut[[i]]@factor)
			if(node@cut[[i]]@side)
				text=paste(text,"is")
			else
				text=paste(text,"not")
		else	
			if(node@cut[[i]]@side)
				text=paste(text,"<=")
			else
				text=paste(text,'>')
		if(node@cut[[i]]@factor)
			text=paste(text,'{',paste(node@cut[[i]]@value,collapse=','),'}', ';')
		else
			text=paste(text,paste(node@cut[[i]]@value,collapse=','), ';')
	}
	if(out)
	{
		print(paste("Subgroup: ",text))
		if(is.null(rdigit))
			print(paste("n=",node@n, "p=",node@signif,collapse=", "))
		else
			print(paste("n=",node@n, "p=",round(node@signif,digit=rdigit),collapse=", "))
	}
	return(text)
}


# Return the Indicators of all Candidate Subgroup
candidate.labs=function(root,n,cutoff=0.025)
{
	if(root@node@signif<cutoff) 
	{
		labs=rbind(is.element(1:n,root@node@sub))
		rownames(labs)=print.node(root@node,out=F)
	}
	else labs=NULL

	if(root@m>0)
		for(i in 1:root@m)
		{
			labs=rbind(labs,candidate.labs(root@children[[i]],n,cutoff))
		}
	labs
}

# List all Candidate Subgroup Objects in a Vector
list.all=function(root)
{
	nodes=list(root@node)
	count=1

	if(root@m>0)
		for(i in 1:root@m)
		{
			temp=list.all(root@children[[i]])
			nodes=append(nodes,temp$nodes)
			count=count+temp$count
		}

	return(list(nodes=nodes,count=count))
}

# ================================================================
# FUNCTIONS FOR PERMUTATION TEST
# ================================================================

permute.test=function(data,yvar, trtvar, xvars, sublist , plist, grid, eff.test, nrepeat=500,uppertail)
{
#  #
  
  datid=1:nrow(data)
	y=data[,yvar]
	trt=data[,trtvar]
	sg.id=1:length(plist)
	sg.grid=plist %*% t(rep(1,length(grid))) <= 
		rep(1,length(plist)) %*% t(grid)

	p.permute=rep(0,length(grid))

	for(i in 1:nrepeat)
	{
		datid=sample(datid)
		y.sig=y[datid]
		trt.sig=trt[datid]

		for(j in 1:length(grid))
		{
			for(k in sg.id[sg.grid[,j]])
			{
				if(eff.test(y.sig[sublist[[k]]],trt.sig[sublist[[k]]],uppertail)$p 
				< plist[k])
				{
					p.permute[j]=p.permute[j]+1/nrepeat
					break
				}
			}
		}
	}

	return(p.permute)
}

permute.select=function(data,yvar, trtvar, xvars,sublist,plist, eff.test ,pm.level,pm.repeat,uppertail)
{
  #

	if(length(plist)==0)
		return(integer())

	if(pm.repeat==0)
		return((1:length(plist))[plist <= pm.level])

	grid=sort(unique(plist))
	
	pm.result=permute.test(data,yvar, trtvar, xvars,sublist,plist,grid, eff.test, pm.repeat,uppertail)

	if(length(grid[pm.result <= pm.level])==0)
		return(integer())

	cutoff=max(grid[pm.result <= pm.level])
	return((1:length(plist))[plist <= cutoff])
}

SIDES.permute=function(root,grid,eff.test,nrepeat=500,uppertail)
{
	data=root@node@data
	templist=list.all(root)
	sublist=vector("list",templist$count-1)
	plist=rep(1,templist$count-1)

	for(i in 2:templist$count)
	{
		sublist[[i-1]]=templist$nodes[[i]]@sub
		plist[i-1]=templist$nodes[[i]]@signif
	}

	weakerr=permute.test(data, sublist , plist, grid, eff.test, nrepeat,uppertail)
	return(list(errtable=rbind(grid,weakerr),plist=plist))
}

# ================================================================
# Cross-Validation
# ================================================================

cv.SIDES <- function(data, yvar, censorvar, trtvar, trtref, xvars, type=type, M=10, split.score=p1, 
                     gamma=c(0.5,0.5,0.5),eff.test=eff.test, saf.test=NULL, maxncov=3, minSplit.prcnt=0.1,
                     minNode=3,maxCand=NULL, pm.level=0.1, pm.repeat=500,uppertail=T, ord.bin=10, k.fold,
                     cv.iter, max.iter) {
    #
    # Perform cross-validation for SIDES.
    #
    
  #
  
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
  
  if (type=="b") strata=data[,yvar]
  if (type=="s") strata=censor.vec[, censorvar]
  if (type=="c") strata=NULL
  
    model.Rfunc <- "SIDES.wrapper"
    # model.Rfunc.args <- make.arg.list("seq.batting")
    # Create a list containing arguments to seq.batting
    
    model.Rfunc.args <- list(yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtref=trtref, xvars=xvars, type=type, M=M, split.score=p1, gamma=gamma, 
                             eff.test=eff.test, saf.test=saf.test, trtref=trtref, maxncov=maxncov, 
                             minSplit.prcnt=minSplit.prcnt, minNode=minNode,maxCand=maxCand, 
                             pm.level=pm.level, pm.repeat=pm.repeat, uppertail=uppertail, ord.bin=ord.bin)
    # List of arguments for model.Rfunc. Must be packaged into 
    # a list to pass to kfold.cv
    
    predict.Rfunc <- "predict.SIDES.cv"
    # List of arguments for predict.Rfunc.
    predict.Rfunc.args <- list(xvars=xvars)
    # Just keeping this here for consistency
    res <- kfold.cv(data=data, model.Rfunc=model.Rfunc, model.Rfunc.args=model.Rfunc.args, predict.Rfunc=predict.Rfunc, predict.Rfunc.args=predict.Rfunc.args, k.fold=k.fold, cv.iter=cv.iter, strata=strata, max.iter=max.iter)
    
    
    
    if (length(res) > 0) {
        stats <- evaluate.cv.results(cv.data=res$cv.data, y=y, censor.vec, trt.vec=trt.vec, type=type)
        

        summary = summarize.cv.stats(stats$raw.stats, trtvar, type)
        interplot=interaction.plot(data.eval=summary, type=type, main="Interaction Plot", trt.lab=c("Trt.", "Ctrl.")) 
                
        results <- list(stats.summary=summary, pred.classes=stats$pred.classes, folds=stats$folds, sig.list=res$sig.list, raw.stats=stats$raw.stats, error.log=res$error.log, interplot=interplot)
    } else {
        results <- "No successful cross-validations."
    }
    results
}

cv.stats.SIDES <- function(cv.data, y, censor.vec, trt.vec, type) {
    # Calculate performance statistics for the output of cv.SIDES
    yvar <- names(y)
    trtvar <- names(trt.vec)
    censorvar <- names(censor.vec)
    if (length(cv.data) > 0) {
        for (i in 1:length(cv.data)) {
            # y is a Surv object, which
            # evaluate.cv.results does not accept. 
            data.tmp <- cv.data[[i]]
            data.tmp.short <- subset(data.tmp, select=c(pred.class, fold))
            data.tmp <- cbind(y, trt.vec, data.tmp.short)
            if (type=="s") {
                data.tmp <- cbind(data.tmp, censor.vec)
            }
            cv.data[[i]] <- data.tmp
        }   
        stats <- evaluate.cv.results(cv.data=cv.data, y=y, censor.vec=censor.vec, trt.vec=trt.vec, type=type)
    } 
}

SIDES.wrapper <- function(data, args) {
    #
    # Wrapper function for SIDES.tree, to be passed to kfold.cv
    #
    # Input arguments:
    #     data - data frame equal to cbind(y, x), where y and x are inputs to SIDES.tree.
    #     args - list containing all other input arguments to SIDES.tree except for x and y. Also contains xvars=names(x)
    #            and yvar=names(y).
    #
    # Output: prediction rule returned by SIDES.tree.
    #
    
    # Run SIDES for Survival Data
    
  
    yvar <- args$yvar
    censorvar <- args$censorvar
    trtvar <- args$trtvar
    trtref <- args$trtref
    type <- args$type 
    xvars<-args$xvars
    M <- args$M
    split.score <- args$split.score
    gamma <- args$gamma
    eff.test <- args$eff.test
    saf.test <- args$saf.test
    maxncov <- args$maxncov
    minSplit.prcnt <- args$minSplit.prcnt
    minNode <- args$minNode
    maxCand <- args$maxCand
    pm.level <- args$pm.level
    pm.repeat <- args$pm.repeat
    uppertail <- args$uppertail
    ord.bin <- args$ord.bin
    
    
    
    
    
    ttree <- SIDES.tree(data=data, yvar=yvar, censorvar=censorvar, trtvar=trtvar, trtref=trtref, xvars=xvars, type=type, M=M, split.score=p1, gamma=gamma,
                        eff.test=eff.test, saf.test=saf.test, maxncov=maxncov, 
                        minSplit.prcnt=minSplit.prcnt, minNode=minNode,maxCand=maxCand, pm.level=pm.level,
                        pm.repeat=pm.repeat,uppertail=uppertail, ord.bin=ord.bin)
    return(ttree)
}

predict.SIDES.cv <- function(data, predict.rule, args) {
    
  
  xvars<- args$xvars
  x <- data[,xvars]
    best <- bestnode(predict.rule)$best
    sig <- get.sig(best)
    
    pred.class <- rep(TRUE, nrow(x))
    if (!is.null(sig$var)) {
        # If sig$var is NULL, return pred.class as all positive.
        for (i in 1:length(sig$var)) {
            x.tmp <- x[,sig$var[i]]
            thresh.tmp <- sig$thresh[i]
            if (sig$side[i]) {
                pred.class.tmp <- (x.tmp <= thresh.tmp)
            } else {
                pred.class.tmp <- (x.tmp > thresh.tmp)
            }
            pred.class <- (pred.class & pred.class.tmp)
        }
    }
    data <- cbind(data, pred.class)
    return(data)
}

predict.sides <- function(data, xvars, predict.rule) {
    x <- data[,xvars]    
    best <- bestnode(predict.rule)$best
    sig <- best@cut
    pred.class <- rep(TRUE, dim(x)[1])
    if (length(sig) > 0) {
        for (i in 1:length(sig)) {
            x.tmp <- x[, sig[[i]]@cov]
            if (sig[[i]]@factor) {
                if (sig[[i]]@side) {
                    pred.class.tmp <- (x.tmp == sig[[i]]@value)
                } else {
                    pred.class.tmp <- (x.tmp != sig[[i]]@value)
                }     
            } else {
                if (sig[[i]]@side) {
                    pred.class.tmp <- (x.tmp <= as.numeric(sig[[i]]@value))
                } else {
                    pred.class.tmp <- (x.tmp > as.numeric(sig[[i]]@value))
                }        
            }
            pred.class <- (pred.class & pred.class.tmp)
        }
        data <- cbind(subset(data, select=xvars), pred.class)
    }
    return(data)
}

# ================================================================
# TUNING PARAMETER: DETERMINE GAMMA
# ================================================================

# Main Function for CV
SIDES.cv=function(data,M,grid=list(seq(0.2,1,0.2),seq(0,1,0.2),seq(0,1,0.2)),nfold=5,
	filename, append=F,
	split.score=p1,eff.test=lm.test,saf.test=NULL, trtref=NULL,
	maxncov=3,minSplit,minNode=3,maxCand=NULL, pm.level=0.1,pm.repeat=500,
  uppertail=TRUE,ord.bin=10)
{
	if(append)
	{
		temp=readfile.cv(filename)
		combo=temp$combo
		cvid=temp$cvid
		nfold=length(cvid)
		
		siglist=rep(0,nrow(combo))
		siglist[1:length(temp$result)]=temp$result
		nexti=length(temp$result)+1
		best.sig=min(temp$result)
		best.comb=combo[siglist==best.sig,]
	}
	else
	{
		L=length(grid)
		g=rep(0,L)
		for(i in 1:L)
			g[i]=length(grid[[i]])

		combo=cv.comb(g,L)

		for(i in 1:L)
			combo[,i]=grid[[i]][combo[,i]]

		combo=drop.nonsense(combo)
	
		cvid=vector("list",nfold)
		newid=sample(nrow(data))
		bound=round(seq(0,nrow(data),nrow(data)/nfold))
		for(i in 1:nfold)
			cvid[[i]]=newid[(bound[i]+1):bound[i+1]]

		siglist=rep(0,nrow(combo))
		best.sig=1
		best.comb=NULL
		nexti=1

		writefile.cv(filename,cvid,combo)
	}

	for(i in nexti:nrow(combo))
	{
		for(j in 1:nfold)
		{
			temp.tree=SIDES.tree(data[-cvid[[j]],],M,split.score, combo[i,],eff.test,saf.test, trtref=trtref,
				maxncov,minSplit,minNode,maxCand, pm.level,pm.repeat,uppertail,ord.bin)
			tempbest=bestnode(temp.tree)$best
			temp.sig=predict.sig(tempbest,data[cvid[[j]],],eff.test,uppertail)
			siglist[i]=siglist[i]+temp.sig/nfold
		}
		if(siglist[i]<best.sig)
		{
			best.sig=siglist[i]
			best.comb=combo[i,]
		}

		fileio=file(filename,"a")
		writeLines(paste(siglist[i]),fileio,"\n")
		close(fileio)
		print(combo[i,])
		print(siglist[i])
	}

	return(list(table=cbind(combo,siglist),best=best.comb,signif=best.sig))
}

# Return the Candidate Subgroup of Best p-value
bestnode=function(root)
{
	best.node=root@node
	best.sig=root@node@signif

	if(root@m>0)
		for(i in 1:root@m)
		{
			temp=bestnode(root@children[[i]])
			if(temp$sig<best.sig)
			{
				best.sig=temp$sig
				best.node=temp$best
			}
		}

	return(list(best=best.node,sig=best.sig))
}

# Get the p-value for Test Set
predict.sig=function(node,newdat,eff.test,uppertail)
{

	if(node@level==0)
		return(eff.test(newdat$y,newdat$trt,newdat$id,uppertail)$p)

	cut=node@cut

	for(i in 1:node@level)
	{
		x=newdat$x[,cut[[i]]@cov]
		if(cut[[i]]@factor)
		{
			if(cut[[i]]@side)
				newdat=newdat[is.element(x,cut[[i]]@value),]
			else
				newdat=newdat[!is.element(x,cut[[i]]@value),]
		}
		else
		{
			if(cut[[i]]@side)
				newdat=newdat[x <= cut[[i]]@value,]
			else
				newdat=newdat[x > cut[[i]]@value,]
		}
	}

	return(eff.test(newdat$y,newdat$trt,newdat$id,uppertail)$p)
}

# List all Possible Combination of Gamma at 
# Different Levels
cv.comb=function(g,L)
{
	temp=rep(1,L)
	combo=NULL

	while(temp[L]<=g[L])
	{
		combo=rbind(combo,temp)
		temp[1]=temp[1]+1
		for(i in 1:(L-1))
			if(temp[i]>g[i])
			{
				temp[i]=1
				temp[i+1]=temp[i+1]+1
			}
			else break
	}

	
	return(combo)
}

# Drop Nonsense Combinations, e.g. 0,0.5,0.5
drop.nonsense=function(combo)
{
	pos=(combo==0)
	sense=apply(t(apply(pos,1,sort))==pos,1,all)
	return(combo[sense,])
}

# File IO for CV
# Save and Resume due to Long Running Time
readfile.cv=function(filename)
{
	temp=file(filename,"r")
	myline=readLines(temp,1)
	if(myline!="cvid")
		stop("invalid file format")
	newid=list()
	myline=readLines(temp,1)
	while(myline!="combo")
	{
		tempid=list(as.numeric(strsplit(myline,' ')[[1]]))
		newid=append(newid,tempid)
		myline=readLines(temp,1)
	}

	newcombo=NULL
	myline=readLines(temp,1)
	while(myline!="result")
	{
		newcombo=rbind(newcombo,as.numeric(strsplit(myline,' ')[[1]]))
		myline=readLines(temp,1)
	}
	myline=readLines(temp,-1)

	if(length(myline)==0) result=numeric(0) 
	else result=as.numeric(myline)

	close(temp)

	return(list(cvid=newid,combo=newcombo,result=result))
}

writefile.cv=function(filename,cvid,combo)
{
	temp=file(filename,"w")
	writeLines("cvid",temp,sep="\n")
	for(i in 1:length(cvid))
		writeLines(paste(cvid[[i]],collapse=" "),temp,sep="\n")
	writeLines("combo",temp,sep="\n")
	for(i in 1:nrow(combo))
		writeLines(paste(combo[i,],collapse=" "),temp,sep="\n")
	writeLines("result",temp,sep="\n")
	close(temp)
}

####### Update 1 #################

# Add a function that returns the column numbers of the split variables
cutlist=function(node)
{
	k=length(node@cut)
	out=NULL
	i=1
	while(i<=k)
	{
		out=c(out,node@cut[[i]]@cov)
		i=i+1
	}

	return(out)
}

# Return the signature for a node.
get.sig=function(node)
{
    k=length(node@cut)
    var=NULL
    thresh=NULL
    side=NULL
    i=1
    while(i<=k)
    {
        var <- c(var,node@cut[[i]]@cov)
        thresh <- c(thresh, as.numeric(node@cut[[i]]@value))
        side <- c(side, node@cut[[i]]@side)
        i=i+1
    }
    sig <- list(var=var, thresh=thresh, side=side)
    return(sig)
}
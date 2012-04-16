#!/usr/bin/env Rscript

# sink(file=devnull, type="message")
library(tools)
library(Biobase)
library(splines)
library(survival)
library(multtest)

library(getopt)
require(methods)


# library(siggenes)

#### BEGIN INLINE
# I inlined library(siggenes) in order to fix a bug in fudge2
#



ignoreThis=setClass("EBAM",representation(z="numeric",posterior="numeric",p0="numeric",local="numeric",
	mat.fdr="matrix",a0="numeric",mat.samp="matrix",vec.pos="numeric",vec.neg="numeric",
	msg="character",chip="character"))

ignoreThis=setMethod("show","EBAM",function(object) print(object))

ignoreThis=setMethod("print","EBAM",
	function(x,delta=NULL,n.digits=4){
		if(is.null(delta))
			mat.fdr<-x@mat.fdr[,1:3,drop=FALSE]
		else
			mat.fdr<-compNumber(x@z,x@posterior,x@p0,nrow(x@mat.samp),delta=delta,
				vec.pos=x@vec.pos,vec.neg=x@vec.neg)[,1:3,drop=FALSE]
		cat(x@msg[1])
		if(length(x@a0)>0)
			cat("Fudge Factor:  a0 =",round(x@a0,n.digits),"\n\n")
		print(round(mat.fdr,n.digits))
	}
)

ignoreThis=setMethod("summary","EBAM",
	function(object,delta=NULL,n.digits=4,what="both",ll=FALSE,chip="",file="",
			sep="\t",quote=FALSE,dec="."){
		if(is.null(delta))
			stop("delta must be specified.")
		if(!what%in%c("both","stats","genes"))
			stop("'what' must be either \"stats\", \"genes\" ",
				"or \"both\".")
		mat.fdr<-compNumber(object@z,object@posterior,object@p0,nrow(object@mat.samp),
			delta=delta,vec.pos=object@vec.pos,vec.neg=object@vec.neg)
		sig.genes<-which(object@z>=mat.fdr[,"CU"] | object@z<=mat.fdr[,"CL"])
		if(what%in%c("genes","both") & length(sig.genes)!=0){
			mat.sig<-cbind(Row=sig.genes,z.value=object@z[sig.genes],
				posterior=object@posterior[sig.genes],
				local.fdr=object@local[sig.genes])
			row.names(mat.sig)<-names(object@z)[sig.genes]
			mat.sig<-mat.sig[rev(order(abs(mat.sig[,"z.value"]))),,drop=FALSE]
			mat.sig<-as.data.frame(mat.sig)
			if(ll){
				if(chip=="" & object@chip==""){
					ll<-FALSE
					warning("Since the chip type is neither specified by ",
						"'chip' nor by the EBAM object,\n",
						"ll is set to FALSE.",call.=FALSE)
				}
				if(all(row.names(mat.sig)==as.character(1:nrow(mat.sig)))){
					ll<-FALSE
					warning("Since no gene names are available, it is not",
						" possible to obtain locus links.\n",
						"Thus, 'll' is set to FALSE.",call.=FALSE)
				}
			}
			if(ll){
				if(chip=="")
					chip<-object@chip
				if(chip!=object@chip & object@chip!="")
					stop("'chip' differs from the chip type of the EBAM object.")
				require(annotate)
				LL<-getLL(row.names(mat.sig),chip)
				sym<-getSYMBOL(row.names(mat.sig),chip)
				mat.sig<-data.frame(Row=mat.sig[,1],Symbol=sym,LocusLink=LL,
					mat.sig[,-1])
			} 
		}
		else
			mat.sig<-data.frame(NULL)
		list.args<-list(n.digits=n.digits,what=what,file=file,sep=sep,quote=quote,
			dec=dec,msg=object@msg,p0=object@p0,a0=object@a0)
		new("sumEBAM",row.sig.genes=sig.genes,mat.fdr=mat.fdr,mat.sig=mat.sig,
			list.args=list.args)
	}
)
		
ignoreThis=setMethod("plot","EBAM",
	function(x,y,pos.stats=NULL,sig.col=3,sig.cex=1,pch=NULL,stats.cex=0.8,main=NULL,
			xlab=NULL,ylab=NULL,y.intersp=1.3,...){
		z<-x@z
		post<-x@posterior
		if(missing(y))
			stop("No delta value has been specified.")
		if(length(y)!=1)
			stop("More than one delta value has been specified.")
		mat.fdr<-compNumber(x@z,x@posterior,x@p0,nrow(x@mat.samp),delta=y,
			vec.pos=x@vec.pos,vec.neg=x@vec.neg)
		if(is.null(main))
			main<-paste("EBAM Plot for Delta =",y)
		if(is.null(xlab))
			xlab<-"z Value"
		if(is.null(ylab))
			ylab<-"Posterior"
		if(length(sig.col)>1)
			stop("sig.col must be of length 1.")
		ids<-which(z<=mat.fdr[,4] | z>=mat.fdr[,5])
		twosided<-any(z<0)
		if(is.null(pos.stats))
			pos.stats<-ifelse(twosided,2,4)
		if(!pos.stats%in%(0:4))
			stop("pos.stats must be an integer between 0 and 4.")
		if(length(ids)==0)
			plot(z,post,main=main,xlab=xlab,ylab=ylab,pch=pch,...)
		else{
			plot(z[-ids],post[-ids],main=main,xlab=xlab,ylab=ylab,pch=pch,
				xlim=range(z),ylim=range(post),...)
			points(z[ids],post[ids],cex=sig.cex,col=sig.col,pch=pch)
		}
		abline(h=y,lty="dashed")
		if(pos.stats!=0){
			tmp<-c("Significant:","FDR:","p0:",if(length(x@a0)==1) "a0:",
				if(twosided) "Cutlow:","Cutup:")
			tmp2<-c(mat.fdr[,2],round(mat.fdr[,3],3),round(x@p0,3),
				if(length(x@a0)==1) round(x@a0,3), 
				if(twosided) round(mat.fdr[,4],3), round(mat.fdr[,5],3))
			textLegend<-paste(tmp,tmp2,sep="  ")
			where<-switch(pos.stats,"top","bottomright","bottomleft","topleft")
			legend(where,legend=textLegend,cex=stats.cex,bty="n",y.intersp=y.intersp)
		}
	}
)



			


					require(methods)

ignoreThis=setClass("FindA0",representation(mat.z="matrix",mat.posterior="matrix",mat.center="matrix",
	mat.success="matrix",mat.failure="matrix",z.norm="numeric",p0="numeric",
	mat.a0="data.frame",mat.samp="matrix",vec.a0="numeric",suggested="numeric",
	delta="numeric",df.ratio="numeric",msg="character",chip="character"))

ignoreThis=setMethod("show","FindA0",function(object) print(object))

ignoreThis=setMethod("print","FindA0",
	function(x,delta=NULL){
		if(is.null(delta))
			delta<-x@delta
		if(length(delta)>1)
			stop("In find.a0, the statistics can be computed only for one value of delta.")
		if(delta<=0 | delta>1)
			stop("delta must be between 0 and 1.")
		cat(x@msg[1])
		cat("Selection Criterion: Posterior >=",delta,"\n\n")
		if(delta==x@delta){
			print(x@mat.a0)
			sugg<-x@suggested
		}
		else{
			tmp<-makeA0mat(x@z.norm,x@mat.posterior,x@p0,x@vec.a0,
				nrow(x@mat.samp),delta=delta)
			print(tmp$tab)
			sugg<-tmp$suggest
		}
		cat("\n","Suggested Choice for a0: ",round(sugg,4),sep="")
		if(names(sugg)!="-")
			cat(" (the ",100*as.numeric(names(sugg)),"% quantile of the s values)",
				sep="")
		cat("\n\n")
	}
)

ignoreThis=setMethod("plot","FindA0",
	function(x,y,logit=TRUE,pos.legend=NULL,legend.cex=0.8,col=NULL,main=NULL,xlab=NULL,
			ylab=NULL,only.a0=FALSE,lty=1,lwd=1,y.intersp=1.1,...){
		mat.post<-x@mat.posterior
		z<-x@z.norm
		if(logit)
			mat.post<-log(mat.post/(1-mat.post))
		if(missing(y))
			y<-x@delta
		if(is.null(main))
			main<-paste("Transformed z Values vs. ",if(logit) "Logit of the ", 
				"Posterior",sep="")
		if(is.null(xlab))
			xlab<-"Transformed z Value"
		if(is.null(ylab))
			ylab<-ifelse(logit,"logit(Posterior)","Posterior")
		if(is.null(col))
			col<-1:ncol(mat.post)
		if(!length(col)%in%c(1,ncol(mat.post)))
			stop("col must be either of length 1 or equal to the number of a0 values.")
		ids<-is.finite(mat.post)
		if(any(!ids))
			warning("Some of the logit posterior probabilites are Inf. ",
				"These probabilities are not plotted.",call.=FALSE)
		ylim<-c(0,max(mat.post[ids]))
		if(is.null(pos.legend))
			pos.legend<-ifelse(any(z<0),1,4)
		if(!pos.legend%in%(0:4))
			stop("pos.legend must be an integer between 0 and 4.")
		plot(range(z),ylim,type="n",main=main,xlab=xlab,ylab=ylab)
		for(i in 1:ncol(mat.post))
			lines(z,mat.post[,i],col=col[i],lwd=lwd,lty=lty)
		h<-ifelse(logit,log(y/(1-y)),y)
		abline(h=h,lty="dashed")
		if(pos.legend!=0){
			if(y!=x@delta)
				mat.legend<-makeA0mat(x@z.norm,x@mat.posterior,x@p0,x@vec.a0,
					nrow(x@mat.samp),delta=y)$tab
			else
				mat.legend<-x@mat.a0
			textLegend<-round(mat.legend[,1],3)
			titleLegend<-"a0"
			if(!only.a0){
				textLegend<-paste(textLegend," (",mat.legend[,3],")",sep="")
				titleLegend<-"a0 (Number)"
			}
			where<-switch(pos.legend,"top","bottomright","bottomleft","topleft")
			legend(where,legend=textLegend,lty=lty,lwd=lwd,cex=legend.cex,col=col,
				bty="n",title=titleLegend,y.intersp=y.intersp)
		}
	}
)
			




			


					Rfold.cal<-function(mat,cl,unlog=TRUE,R.fold=1){
	if(R.fold<=0)
		stop("R.fold must be larger than 0.")
	if(R.fold<1)
		R.fold<-1/R.fold
	if(unlog)
		mat<-2^mat
	uni.cl<-sort(unique(cl))
	g1<-if(length(uni.cl)==2) which(cl==uni.cl[1]) else which(cl<0)
	mean.g1<-rowMeans(mat[,g1])
	mean.g1[mean.g1<=0]<-NA
	mean.g2<-rowMeans(mat[,-g1])
	mean.g2[mean.g2<=0]<-NA
	fold<-mean.g2/mean.g1
	fulfill<-numeric(length(fold))
	fulfill[fold>=R.fold | fold<=1/R.fold]<-1
	cbind(fold=fold,fulfill=fulfill)
}
	 
 

require(methods)
ignoreThis=setClass("SAM",representation(d="numeric",d.bar="numeric",vec.false="numeric",p.value="numeric",
	s="numeric",s0="numeric",mat.samp="matrix",p0="numeric",mat.fdr="matrix",q.value="numeric",
	fold="numeric",msg="character",chip="character"))

ignoreThis=setMethod("show","SAM",
	function(object){
		cat(object@msg[1])
		print(pretty.mat.fdr(object@mat.fdr[,1:5],digits=3))
	}
)

ignoreThis=setMethod("print","SAM",
	function(x,delta=NULL,n.digits=3){
		cat(x@msg[1])
		if(is.null(delta))
			print(pretty.mat.fdr(x@mat.fdr[,1:5],digits=n.digits))
		else{
			mat.fdr<-as.data.frame(stats.cal(x@d,x@d.bar,x@vec.false,x@p0,
				delta=delta))
			print(pretty.mat.fdr(mat.fdr[,1:5],digits=n.digits))
		}
	}
)



ignoreThis=setMethod("summary","SAM",
	function(object,delta=NULL,n.digits=3,what="both",ll=FALSE,chip="",file="",
			sep="\t",quote=FALSE,dec="."){
		list.args<-list(n.digits=n.digits,what=what,file=file,sep=sep,quote=quote,
			dec=dec,msg=object@msg)
		if(length(delta)!=1){
			mat.fdr<-if(is.null(delta)) object@mat.fdr
				else  stats.cal(object@d,object@d.bar,object@vec.false,
					object@p0,delta=delta)
			sig.genes<-numeric(0)
			mat.sig<-data.frame(NULL)
		}
		else{
			if(!what%in%c("both","stats","genes"))
				stop("'what' must be either \"stats\", \"genes\" or \"both\".")
			mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,
				delta=delta)
			sig.genes<-which(object@d>=mat.fdr[,"cutup"] | object@d<=mat.fdr[,"cutlow"])
			mat.sig<-data.frame(NULL)
			if(what%in%c("genes","both") & length(sig.genes)!=0){
				if(chip=="" & object@chip=="" & ll){
					ll<-FALSE
					warning("Since the chip type is neither specified by 'chip' ",
						"nor by the SAM object,","\n","ll is set to FALSE.",
						call.=FALSE)
				}
				mat.sig<-cbind(Row=(1:length(object@d))[sig.genes],
					d.value=object@d[sig.genes],stdev=object@s[sig.genes],
					rawp=object@p.value[sig.genes],
					q.value=object@q.value[sig.genes],
					R.fold=object@fold[sig.genes])
				if(length(sig.genes)>1)
					mat.sig<-mat.sig[rev(order(abs(mat.sig[,"d.value"]))),]
				mat.sig<-as.data.frame(mat.sig)
				if(ll & all(row.names(mat.sig)==as.character(1:nrow(mat.sig)))){
					ll<-FALSE
					warning("Since no gene names are available it is not",
						" possible to obtain locus links.\n",
						"'ll' is thus set to FALSE.",call.=FALSE)
				}
				if(ll){
					if(chip=="")
						chip<-object@chip
					if(chip!=object@chip & object@chip!="")
						stop("'chip' differs from the chip type of the SAM object.")
					require(annotate)
					LL<-getLL(row.names(mat.sig),chip)
					sym<-getSYMBOL(row.names(mat.sig),chip)
					mat.sig<-data.frame(Row=mat.sig[,1],Symbol=sym,
						LocusLink=LL,mat.sig[,-1])
				}
			}
			retval<-new("sumSAM",row.sig.genes=sig.genes,mat.fdr=mat.fdr,
				mat.sig=mat.sig,list.args=list.args)
		}
		new("sumSAM",row.sig.genes=sig.genes,mat.fdr=mat.fdr,mat.sig=mat.sig,
			list.args=list.args)
	}
)

ignoreThis=setMethod("plot","SAM",
	function(x,y,pos.stats=NULL,sig.col=3,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,ylab=NULL,
		pty="s",lab=c(10,10,7),pch=NULL,sig.cex=1,helplines=FALSE,...){
		if(missing(y)){
			y<-NULL
			cat("To obtain a SAM plot, delta has to be specified in plot(object,delta,...).",
				"\n")
		}	
		if(length(y)==1)
			sam.plot2(x,delta=y,pos.stats=pos.stats,sig.col=sig.col,xlim=xlim,ylim=ylim,
				main=main,xlab=xlab,ylab=ylab,pty=pty,lab=lab,pch=pch,sig.cex=sig.cex,
				...)
		else
			delta.plot(x,delta=y,helplines=helplines)
	})

ignoreThis=setMethod("identify","SAM",
	function(x,showText=TRUE,getInfo=TRUE,pos=4,cex=0.8,add.xy=numeric(2),n.digits=3,
			ask=FALSE,ll=FALSE,browse=FALSE,chip="",...){
		if(length(add.xy)!=2)
			stop("add.xy must have length 2.")
		d<-x@d
		if(is.null(names(d)))
			names(d)<-paste("G",1:length(d),sep="")
		sorted.d<-sort(d,index.return=TRUE)
		d.sort<-sorted.d$x
		mat.xy<-cbind(x@d.bar,d.sort)
		if(getInfo & chip=="" & x@chip=="" & ll){
			ll<-FALSE
			warning("Since the chip type is neither specified by 'chip' nor by the ",
				"SAM object,","\n","ll is set to FALSE.",call.=FALSE)
			}
		repeat{
			id.out<-identify(mat.xy,labels=names(d.sort),pos=TRUE,n=1,plot=FALSE)
			if(length(id.out$pos)==0)
				break
			ind<-id.out$ind
			x.ided<-mat.xy[ind,1]
			y.ided<-mat.xy[ind,2]
			if(showText)
				text(x.ided+add.xy[1],y.ided+add.xy[2],names(d.sort)[ind],pos=pos,
					cex=cex,...)
			if(any((x.ided-mat.xy[-ind,1])^2+(y.ided-mat.xy[-ind,2])^2==0)){
				same.ided<-which((x.ided-mat.xy[,1])^2+(y.ided-mat.xy[,2])^2==0)
				if(!getInfo)
					warning("There are ",length(same.ided)-1," other genes ",
						"having the same coordinates as ",names(d.sort)[ind],
						".",call.=FALSE)
			}
			else
				same.ided<-ind
			if(getInfo){
				which.d<-sorted.d$ix[same.ided]
				mat.ided<-cbind(d.value=na.exclude(d)[which.d],
					stdev=if(length(x@s)==0) numeric(0) 
						else na.exclude(x@s)[which.d],
					rawp=na.exclude(x@p.value)[which.d],
					q.value=na.exclude(x@q.value)[which.d],
					R.fold=x@fold[!is.na(d)][which.d])
				mat.ided<-data.frame(round(mat.ided,n.digits))
				if(ll){
					if(chip=="")
						chip<-x@chip
					if(chip!=x@chip & x@chip!="")
						stop("'chip' differs from the chip type of the SAM object.")
					require(annotate)
					LL<-getLL(row.names(mat.ided),chip)
					sym<-getSYMBOL(row.names(mat.ided),chip)
					mat.ided<-data.frame(Symbol=sym,LocusLink=LL,mat.ided)
					if(browse){
						for(i in 1:length(na.exclude(LL))){
							browseURL(getQuery4LL(na.exclude(LL)[i]))
							if(i<length(na.exclude(LL))){
								answer2<-readline("Next LocusLink?")
								if(tolower(answer2)%in%c("n","no"))
									break
							}
						}
					}
				}		
				print(mat.ided)
				cat("\n")
				if(ask){
					answer<-readline("Press <Return> to continue. Type 'n' to stop. ")
					if(tolower(answer)%in%c("n","no"))
						break
				}
				cat("\n")				
			}
		}
        invisible(id.out)
	}
)



add.target2href<-function(tr,split="\">"){
	tmp<-strsplit(tr,split)
	for(i in 1:length(tmp)){
		tmp1<-tmp[[i]]
		if(length(tmp1)>1)
			tr[i]<-paste(tmp1[1],paste("\" target=\"_blank\">",tmp1[-1],
				sep="",collapse=""),sep="")
	}
	tr
}



adjust.for.mt<-function(data,cl,var.equal=FALSE,eb=FALSE,wilc=FALSE){
	if(is.data.frame(cl) || is.matrix(cl))
		cl<-pairt.cl.transform(cl,ncol(data))
	ana.type<-ifelse(eb,"EBAM","SAM")
	n.cl<-length(cl)
	if(n.cl!=ncol(data))
		stop("The length of cl must be equal to the number of columns of data.")
	lev<-unique(cl)
	uni.cl<-length(lev)
	if(any(lev<0)){
		if(any(lev==0) | length(lev[lev>0])!=length(lev[lev<0]))
			stop(paste("Negative values are only allowed for the paired analysis. Or,", "\n",
				"if a paired analysis should be done: There is something wrong with the class labels."))
		uni.cl.abs<-uni.cl/2
	}
	else
		uni.cl.abs<-uni.cl
	type.mt<-NULL
	if(uni.cl==1){
		type.mt<-"pairt"
		X<-matrix(0,nrow(data),2*n.cl)
		X[,seq(1,2*n.cl-1,2)]<-as.matrix(data)
		msg<-paste(ana.type,"Analysis for the One-Class Case",
			if(wilc) "Using Wilcoxon Signed Rank Statistics","\n\n")
		if(lev!=1)
			warning("Expected class label is 1, ",lev," is thus set to 1.",call.=FALSE)
		cl.mt<-rep(c(1,0),n.cl)
	}
	if(uni.cl==2){
		msg<-paste(ana.type,"Analysis for the Two-Class Unpaired Case",
			ifelse(!wilc,paste("Assuming",ifelse(var.equal,"Equal","Unequal"),"Variances"),
			"Using Wilcoxon Rank Sums"),"\n","\n")
		if(any(table(cl)<2))
			stop("Each group must consist of at least two samples.")
		type.mt<-ifelse(var.equal,"t.equalvar","t")
		X<-as.matrix(data)
		if(!all(sort(lev)==0:1)){
			warning("Expected class labels are 0 and 1. ",sort(lev)[1]," is set to 0, and ",
				sort(lev)[2]," is set to 1.",call.=FALSE)
			cl[which(cl==sort(lev)[1])]<-0
			cl[which(cl==sort(lev)[2])]<-1
		}
		cl.mt<-as.numeric(cl)
		
	}
	if(uni.cl==2*uni.cl.abs & uni.cl==n.cl){
		msg<-paste(ana.type,"Analysis for the Two-Class Paired Case",
			if(wilc) "Using Wilcoxon Signed Rank Scores","\n\n")
		type.mt<-"pairt"
		sort.cl<-sort(cl,index=TRUE)
		if(!all(sort.cl$x==c(-uni.cl.abs:-1,1:uni.cl.abs)))
			stop(paste("The class labels must be the integers between ",-uni.cl.abs,
				" and 1, and between 1 and ",uni.cl.abs,".",sep=""))
		x<-sort.cl$ix[(uni.cl.abs+1):uni.cl]
		y<-sort.cl$ix[uni.cl.abs:1]
		index<-as.vector(rbind(x,y))
		X<-as.matrix(data[,index])	
		cl.mt<-rep(c(1,0),n.cl/2)
	}
	if(uni.cl>2 & all(lev>=0)){
		if(wilc)
			stop("Currently no rank-based SAM version for the multi-class case available.")
		msg<-paste(ana.type,"Analysis for the Multi-Class Case with",uni.cl,"Classes","\n","\n")	
		type.mt<-"f"
		if(any(table(cl)<2))
			stop("Each group must consists of at least two samples.","\n","\n")
		sort.cl<-sort(lev)
		cl.new<-cl
		if(!all(sort.cl==1:uni.cl)){
			mnc<-matrix(0,uni.cl,2)
			for(i in 1:uni.cl){
				cl.new[which(cl==sort.cl[i])]<-i
				mnc[i,]<-c(sort.cl[i],i)
			}
			warning("Expected class labels are the integers between 1 and ",uni.cl,".",
				"\n","The new class labels are thus ",paste(mnc[,2],"(was ",mnc[,1],
				ifelse(mnc[,2]!=mnc[nrow(mnc),2],"), ",")."),sep="",collapse=""),
				call.=FALSE)
		}
		X<-as.matrix(data)
		cl.mt<-as.numeric(cl.new)-1
	}
	if(is.null(type.mt))
		stop("There is something wrong with the class labels.")
	mode(X)<-"numeric"
	structure(list(X=X,cl.mt=cl.mt,type.mt=type.mt,msg=msg))
}
 
 

`args.ebam` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","summary","plot"))
		stop("'method' must be either print, summary, or plot.")
	if(method=="print")
		cat("function (x, delta = NULL, n.digits = 4)\n")
	if(method=="summary")
		cat("function(object, delta = NULL, n.digits = 4, what = \"both\",",
			"ll = FALSE,\n    chip = \"\", file = \"\", sep = \"\\t\",",
			"quote = FALSE, dec = \".\")\n")
	if(method=="plot")
		cat("function (x, y, pos.stats = NULL, sig.col = 3, sig.cex = 1,",
			"pch = NULL,\n    cexStats = 0.8, main = NULL, xlab = NULL,",
			"ylab = NULL, y.intersp = 1.3, ...)\n")
}

`args.finda0` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","plot"))
		stop("'method' must be either print or plot.")
	if(method=="print")
		cat("function (x, delta = NULL)\n")
	if(method=="plot")
		cat("function (x, y, logit = TRUE, pos.legend = NULL, cexLegend = 0.8,",
			"col = NULL,\n    main = NULL, xlab = NULL, ylab = NULL,",
			"only.a0 = FALSE, lty = 1, lwd = 1,\n    y.intersp = 1.1, ...)\n")
}

args.sam<-function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","summary","plot","identify"))
		stop("'method' must be either print, summary, plot or identify.")
	if(method=="print")
		cat("function (x, delta = NULL, n.digits = 3)\n")
	if(method=="summary")
		cat("function(object, delta = NULL, n.digits = 5, what = \"both\",",
			"ll = FALSE,\n    chip = \"\", file = \"\", sep = \"\\t\",",
			"quote = FALSE, dec = \".\")\n")
	if(method=="plot")
		cat("function (x, y, pos.stats = NULL, sig.col = 3, xlim = NULL,",
			"ylim = NULL,\n    main = NULL, xlab = NULL, ylab = NULL,",
			"pty = \"s\", lab = c(10, 10, 7),\n    pch = NULL, sig.cex = 1,",
			"helplines = TRUE,","...)\n")
	if(method=="identify")
		cat("function (x, showText = TRUE, getInfo = TRUE, pos = 4, cex = 0.8,\n",
			"   add.xy = numeric(2), n.digits = 4, ask = FALSE,",
			"ll = FALSE,\n    browse = FALSE, chip= \"\", ...)\n")
}


build.dperm<-function(X,tmp.samp,type.mt,s0,n.row,le.cl){
	B<-nrow(tmp.samp)
	max.cl<-max(tmp.samp[1,])+1
	d.perm<-matrix(0,n.row,B)
	print("hello")
	for(i in 1:B){
		tmp<-.C("get_stat_num_denum", as.double(X), as.integer(n.row), 
        	as.integer(le.cl), as.integer(tmp.samp[i,]), as.double(.mt.naNUM), 
        	t.num = double(n.row), t.denum=double(n.row),as.character(c(type.mt,"abs","y")), 
		as.integer(max.cl), PACKAGE = "multtest")
		d.perm[,i]<-as.matrix(sort(tmp$t.num/(tmp$t.denum+s0)))
	}
	d.perm
}	 
 

cat.ebam<-function(data,cl,approx=FALSE,B=100,check.levels=TRUE,check.for.NN=FALSE,lev=NULL,
		B.more=0.1,B.max=50000,n.subset=10,fast=FALSE,n.interval=139,df.ratio=3,
		df.glm=NULL,rand=NA){
	data<-as.matrix(data)
	if(any(is.na(data)))
		stop("No NAs allowed.")
	if(check.for.NN){
		if(any(data=="NN"))
			stop("No NNs allowed.")
	}
	if(ncol(data)!=length(cl))
		stop("The number of columns of data must be equal to the length of cl.")
	if(!is.null(lev)){
		for(i in 1:length(lev))
			data[data==lev[i]]<-i
	}
	if(mode(data)!="numeric")
		mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("data must contain integers.")
	n.cat<-max(data)
	if(any(!data%in%1:n.cat))
		stop("data must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%data))
		stop("Some of the values between 1 and ",n.cat," are not in data.")
	stats<-chisqClass(data,cl,n.cat,check=check.levels)
	if(is.null(n.interval))
		n.interval<-ifelse(approx,floor(sqrt(length(stats))),139)
	n.cl<-length(unique(cl))
	msg<-c("EBAM Analysis for Categorical Data\n\n",
		paste("Null Distribution:\n",
		ifelse(approx,"Approximation by ChiSquare Distribution",
		paste("Estimated Based on",B,"Permutations")),"\n\n",sep=""))
	if(approx){
		fail.out<-cat.null.approx2(stats,n.cl,n.cat,n.interval=n.interval,df.glm=df.glm)
		return(list(z=stats,ratio=fail.out$ratio,vec.pos=fail.out$vec.pos,
			vec.neg=fail.out$vec.neg,msg=msg))
	}
	out<-getSuccesses(stats,n.interval=n.interval)
	mat.samp<-setup.mat.samp(cl,"f",B=B,B.more=B.more,B.max=B.max,rand=rand)
	fail.out<-compFailure(data,mat.samp,stats,out$interval,n.subset=n.subset,
		fast=fast,n.cat=n.cat)
	ratio<-compRatio(out$center,out$success,fail.out$vec.fail,df=df.ratio,z=stats)$ratio
	if(fast)
		return(list(z=stats,ratio=ratio,success=out$success,failure=fail.out$vec.fail,
			center=out$center,mat.samp=mat.samp,msg=msg))
	else
		return(list(z=stats,ratio=ratio,vec.pos=fail.out$vec.pos/B,
			vec.neg=fail.out$vec.neg/B,mat.samp=mat.samp,msg=msg))
}

"cat.null" <-
function(data,mat.samp,d,n.subset,n.cat){
	B<-nrow(mat.samp)
	vecB<-c(seq(1,B,n.subset),B+1)
	vecB<-unique(vecB)
	n.B<-length(vecB)-1
	n.var<-nrow(data)
	vec.dperm<-vec.false<-numeric(n.var)
	d.rank<-rank(-d,ties="first")
	for(i in 1:n.B){
		tmp<-compPermStat(data,mat.samp[vecB[i]:(vecB[i+1]-1),],n.cat)
		tmp<-apply(tmp,2,sort)
		vec.dperm<-vec.dperm+rowSums(tmp)
		tmp2<-c(as.vector(tmp),d)
		vec.false<-vec.false+rank(-tmp2,ties="first")[length(tmp)+(1:n.var)]-d.rank
	}
	d.bar<-vec.dperm/B
	vec.false<-vec.false/B
	p.value<-vec.false/n.var
	return(list(d.bar=d.bar,vec.false=vec.false,p.value=p.value))
}

"cat.null.approx" <-
function(d,n.cl,n.cat){
	df<-(n.cl-1)*(n.cat-1)
	n.var<-length(d)
	d.bar<-qchisq(((1:n.var)-0.5)/n.var,df)
	p.value<-pchisq(d,df,lower.tail=FALSE)
	vec.false<-p.value*n.var
	return(list(d.bar=d.bar,p.value=p.value,vec.false=vec.false))
}



	

cat.null.approx2<-function(z,n.cl,n.cat,n.interval=NULL,df.glm=NULL){
	df.chisq<-(n.cl-1)*(n.cat-1)
	p.value<-1-pchisq(z,df.chisq)
	n.z<-length(z)
	vec.pos<-length(z)*p.value
	z.null<-dchisq(z,df.chisq)
	if(is.null(df.glm))
		df.glm<-min(df.chisq+1,5)
	z.fit<-denspr(z,n.interval=n.interval,df=df.glm)
	return(list(ratio=z.null/z.fit,vec.pos=vec.pos,vec.neg=numeric(n.z)))
}

cat.stat<-function(data,cl,B=100,approx=FALSE,check.levels=TRUE,check.for.NN=FALSE,lev=NULL,
		B.more=0.1,B.max=50000,n.subset=10,rand=NA){
	data<-as.matrix(data)
	if(any(is.na(data)))
		stop("No NAs allowed.")
	if(check.for.NN){
		if(any(data=="NN"))
			stop("No NNs allowed.")
	}
	if(ncol(data)!=length(cl))
		stop("The number of columns of data must be equal to the length of cl.")
	if(!is.null(lev)){
		for(i in 1:length(lev))
			data[data==lev[i]]<-i
	}
	if(mode(data)!="numeric")
		mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("data must contain integers.")
	n.cat<-max(data)
	if(any(!data%in%1:n.cat))
		stop("data must consist of integers between 1 and ",n.cat,".")
	if(any(!(1:n.cat)%in%data))
		stop("Some of the values between 1 and ",n.cat," are not in data.")
	stats<-chisqClass(data,cl,n.cat,check=check.levels)
	n.cl<-length(unique(cl))
	if(approx){
		null.out<-cat.null.approx(stats,n.cl,n.cat)
		mat.samp<-matrix(numeric(0))
	}
	else{
		mat.samp<-setup.mat.samp(cl,"f",B=B,B.more=B.more,B.max=B.max,
			rand=rand)
		null.out<-cat.null(data,mat.samp,stats,n.subset,n.cat)
	}
	msg<-c("SAM Analysis for Categorical Data\n\n",
		paste("Null Distribution:\n",
		ifelse(approx,"Approximation by ChiSquare Distribution",
		paste("Estimated Based on",B,"Permutations")),"\n\n",sep=""))
	structure(list(d=stats,d.bar=null.out$d.bar,p.value=null.out$p.value,
		vec.false=null.out$vec.false,discrete=FALSE,s=numeric(0),
		s0=numeric(0),mat.samp=mat.samp,msg=msg,fold=numeric(0)))
}

check.chipname<-function(chip,obj.chip,cdfname=NULL){
		if(chip=="" & obj.chip=="" & is.null(cdfname))
			stop("The chip type must be specified either by the SAM object",
				"\n or by 'chipname' or 'cdfname'.")
		if(chip=="")
			chip<-obj.chip
		if(chip!=obj.chip & obj.chip!="")
			stop("'chipname' differs from the chip type of the SAM object.")
		if(!is.null(cdfname)){
			if(!is.character(cdfname))
				stop("'cdfname' must be a character string.")
			tmp<-tolower(cdfname)
			tmp<-gsub("_","",tmp)
			tmp<-gsub("-","",tmp)
			tmp<-gsub(" ","",tmp)
			if(chip!="" & tmp!=chip)
				stop("'chipname' differs from the chip type specified by",
					" 'cdfname'.")
			chip<-tmp
		}
		chip
}



`checkA0` <-
function(s,a0=NULL,quan.a0=NULL){
	if(is.null(a0) & is.null(quan.a0))
		stop("Either a0 or quan.a0 must be specified. See the help of z.ebam.")
	if(!is.null(a0) & !is.null(quan.a0))
		stop("Only one of a0 or quan.a0 should be specified.")
	if(!is.null(a0)){
		if(length(a0)>1 | !is.numeric(a0))
			stop("a0 should be a numeric value (and not a vector of values).")
		if(a0<0)
			stop("a0 must be larger than or equal to 0.")
		return(a0)
	}
	if(length(quan.a0)>1 | !is.numeric(quan.a0))
		stop("quan.a0 must be a numeric value (and not, e.g., a vector).")
	if(quan.a0<0 | quan.a0>1)
		stop("quan.a0 must be between 0 and 1.")
	quantile(s,quan.a0)
}

`checkFUNout` <-
function(fun.out){
	tmp.names<-names(fun.out)
	if(any(!c("z","ratio")%in%tmp.names))
		stop("The output of the function specified by method must contain objects\n",
			"called z and ratio.")
	if(length(fun.out$z)!=length(fun.out$ratio))
		stop("z and ratio must have the same length.")
	exact<-ifelse(any(tmp.names=="vec.pos"),TRUE,FALSE)
	a0<-if(any(tmp.names=="a0")) fun.out$a0 else numeric(0)
	msg<-if(any(tmp.names=="msg")) fun.out$msg else ""
	if(any(tmp.names=="mat.samp")){
		mat.samp<-fun.out$mat.samp
		B<-nrow(mat.samp)
	}
	else{
		mat.samp<-matrix(numeric(0))
		B<-NA
	}
	if(!exact)
		vec.pos<-vec.neg<-numeric(0)
	else{
		vec.pos<-fun.out$vec.pos
		if(length(vec.pos)!=length(fun.out$z))
			stop("vec.pos must have the same length as z.")
		if(any(tmp.names=="vec.neg")){
			vec.neg<-fun.out$vec.neg
			if(length(vec.neg)!=length(vec.pos))
				stop("vec.neg must have the same length as vec.pos.")
		}
		else{
			twosided<-ifelse(any(fun.out$z<0),TRUE,FALSE)
			vec.neg<-if(twosided) vec.pos else numeric(length(vec.pos))
			warning("vec.neg is not specified even though vec.pos is specified.\n",
				"Since ",ifelse(twosided,"some","none")," of the z values ",
				ifelse(twosided,"are","is")," smaller than zero, ",
				"vec.neg is set to ",ifelse(twosided,"vec.pos","0"),".",
				call.=FALSE)
		}
	}
	return(list(vec.pos=vec.pos,vec.neg=vec.neg,exact=exact,B=B,mat.samp=mat.samp,a0=a0,msg=msg))
}

`checkQuantiles` <-
function(s,quan.a0=(0:10)/100,include.zero=TRUE){
	if(length(quan.a0)<2)
		stop("quan.a0 must be a vector containing at least two values.")
	if(any(quan.a0<0 | quan.a0>1))
		stop("The values in quan.a0 must be between 0 and 1.")
	if (any(round(100 * quan.a0, 10) != round(100 * quan.a0,0))){
        	warning("At least one alpha is not a percentile. Only the first two decimal digits", 
            		" are retained.")
        	quan.a0 <- signif(quan.a0, 2)
    	}
	quans<-quantile(s,quan.a0)
	if(include.zero)
		quans<-c(0,quans)
	quans
}

"chisqClass" <-
function(data,cl,n.cat,check=TRUE,pam=FALSE){
	if(any(is.na(data)))
		stop("No NAs allowed.")
	if(missing(n.cat))
		n.cat<-max(data)
	uni.cl<-sort(unique(cl))
	n.lev<-length(uni.cl)
	if(length(n.lev)>10)
		stop("cl contains more than 10 different values.")
	if(any(uni.cl!=1:n.lev))
		stop("The labels of the classes must be 1 to ",n.lev,".")
	n.obs<-ncol(data)
	n.snp<-nrow(data)
	if(length(cl)!=n.obs)
		stop("The length of cl must be equal to the number of observations.")
	CL<-matrix(0,n.obs,n.lev)
	for(i in 1:n.lev)
		CL[cl==i,i]<-1
	vec.ncl<-colSums(CL)
	listCells<-vector("list",n.cat)
	for(i in 1:n.cat)
		listCells[[i]]<-data==i
	if(check){
		for(i in 1:n.cat){
			tmp.rS<-rowSums(listCells[[i]])
			if(any(tmp.rS==0))
				stop("All variables must have the same number of levels.")
		}
	}
	if(pam){
		mat.stats<-matrix(0,n.snp,n.lev)
		list.Nobs<-vector("list",n.cat)
		list.Nexp<-vector("list",n.cat)
		for(i in 1:n.cat){
			tmp<-computeContCols(listCells[[i]],CL,vec.ncl,n.obs,pam=TRUE)
			mat.stats<-mat.stats+tmp$mat.stat
			list.Nobs[[i]]<-tmp$N.obs
			list.Nexp[[i]]<-tmp$N.exp
		}
		tmp<-t(vec.ncl)%x%rep(1,n.snp)
		out<-mat.stats-tmp
		if(!is.null(rownames(data)))
			rownames(out)<-rownames(data)
		return(list(mat.stat=out,list.Nobs=list.Nobs,list.Nexp=list.Nexp,n=n.obs))
	}
	else{
		listNexp<-lapply(listCells,computeContCols,CL=CL,vec.ncl=vec.ncl,n.obs=n.obs)
		mat.stats<-matrix(unlist(listNexp),ncol=n.cat)
		out<-rowSums(mat.stats)-n.obs
	}
	if(!is.null(rownames(data)))
		names(out)<-rownames(data)
	out
}

col2hex<-function(col){
	rgb<-col2rgb(col)[,1]
	hex<-c(0:9,LETTERS[1:6])
	paste("#",paste(hex[rgb%/%16+1],hex[rgb%%16+1],sep="",collapse=""),sep="")
}



`compFailure` <-
function(data,mat.samp,z,interval,a0=0,type.mt=NULL,n.subset=5,fast=FALSE,n.cat=NULL){
	B<-nrow(mat.samp)
	seq.samp<-unique(c(seq(1,B,ceiling(B/n.subset)),B+1))
	n.int<-length(interval)-1
	n.seq<-length(seq.samp)-1
	vec.fail<-numeric(n.int)
	n.row<-nrow(data)
	le.cl<-ncol(mat.samp)
	z.range<-range(z)
	if(!fast){
		vec.pos<-vec.neg<-numeric(n.row)
		z.rank<-rank(-abs(z),ties="first")
	}
	else
		vec.pos<-vec.neg<-NULL
	for(i in 1:n.seq){
		tmp.samp<-mat.samp[seq.samp[i]:(seq.samp[i+1]-1),]
		if(is.null(n.cat))
			z.perm<-build.dperm(data,tmp.samp,type.mt,a0,n.row,le.cl)
		else
			z.perm<-compPermStat(data,tmp.samp,n.cat)
		z.perm<-as.vector(z.perm)
		tmp<-getFailure(z.perm,z,interval,z.range=z.range,n.interval=n.int)
		vec.fail<-vec.fail+tmp
		if(!fast){
			tmp<-compFalse(z,z.perm,z.rank,n.row)
			vec.pos<-vec.pos+tmp$vec.pos
			vec.neg<-vec.neg+tmp$vec.neg
		}
	}
	return(list(vec.fail=vec.fail,vec.pos=vec.pos,vec.neg=vec.neg))
}

`compFailureMat2` <-
function(data,mat.samp,z.mat,ints,z.norm,ints.norm,FUN,vec.a0,n.chunk=1,
		n.interval=139){
	B<-nrow(mat.samp)
	seq.samp<-unique(c(seq(1,B,ceiling(B/n.chunk)),B+1))
	n.seq<-length(seq.samp)-1
	mat.failure<-mat.fail.norm<-matrix(0,n.interval,length(vec.a0))
	z.range<-apply(z.mat,2,range)
	for(j in 1:n.seq){
		tmp.samp<-mat.samp[seq.samp[j]:(seq.samp[j+1]-1),]
		tmp<-compFailureSubset(data,tmp.samp,z.mat,ints,z.norm,ints.norm,FUN,
			vec.a0,z.range,n.interval=n.interval)
		mat.failure<-mat.failure+tmp$tmp.fail
		mat.fail.norm<-mat.fail.norm+tmp$tmp.fail.norm
	}
	colnames(mat.failure)<-colnames(mat.fail.norm)<-round(vec.a0,4)
	return(list(mat.failure=mat.failure,mat.fail.norm=mat.fail.norm))
}

`compFailureSubset` <-
function(data,mat.samp,z.mat,ints,z.norm,ints.norm,FUN,vec.a0,
		z.range,n.interval=139){
	B<-nrow(mat.samp)
	mat.r<-mat.s<-matrix(0,nrow(data),B)
	n.a0<-length(vec.a0)
	tmp.fail<-tmp.fail.norm<-matrix(0,n.interval,n.a0)
	for(i in 1:B){
		fun.out<-FUN(data,mat.samp[i,])
		mat.r[,i]<-fun.out$r
		mat.s[,i]<-fun.out$s
	}
	for(i in 1:n.a0){
		tmp<-as.vector(mat.r/(mat.s+vec.a0[i]))
		tmp.fail[,i]<-getFailure(tmp,z.mat[,i],ints[[i]],z.range=z.range[,i])
		tmp.fail.norm[,i]<-getFailure(tmp,z.mat[,i],ints.norm,z.norm=z.norm,
			n.interval=n.interval)
	}
	return(list(tmp.fail=tmp.fail,tmp.fail.norm=tmp.fail.norm))
}

`compFalse` <-
function(z,z.perm,z.rank,n.row){
	ids<-which(z.perm<0)
	n.neg<-length(ids)
	if(n.neg>0)
		vec.neg<-rank(c(z.perm[ids],-abs(z)),ties="first")[n.neg+(1:n.row)]-z.rank
	else
		vec.neg<-numeric(n.row)
	ids<-which(z.perm>=0)
	vec.pos<-rank(-c(z.perm[ids],abs(z)),ties="first")[length(ids)+(1:n.row)]-z.rank
	return(list(vec.pos=vec.pos,vec.neg=vec.neg))
}

`compNumber` <-
function(z,post,p0,B,delta=0.9,vec.pos=NULL,vec.neg=NULL){
	if(any(delta<=0 | delta>1))
		stop("The delta values must be between 0 and 1.")
	z.sort<-sort(z)
	z.order<-order(z)
	post<-post[z.order]
	if(length(vec.pos)==0)
		probs<-1/(B*(1-post)/p0+1)
	else{
		vec.pos<-vec.pos[z.order]
		vec.neg<-vec.neg[z.order]
	}
	n.delta<-length(delta)
	m<-length(z)
	mat.delta<-matrix(0,n.delta,5)
	rownames(mat.delta)<-1:n.delta
	colnames(mat.delta)<-c("Delta","Number","FDR","CL","CU")
	mat.delta[,1]<-delta
	for(i in 1:n.delta){
		if(any(z.sort<0) & post[1]>=delta[i]){
			tmp<-post[z.sort<0]
			j1<-min(which(tmp<delta[i]))-1
			f1<-if(length(vec.neg)==0) sum(1/probs[1:j1]-1)/B
				else vec.neg[j1]
			mat.delta[i,4]<-z.sort[j1]
		}
		else{
			j1<-0
			f1<-0
			mat.delta[i,4]<- -Inf
		}
		if(post[m]>=delta[i]){
			tmp<-rev(post[z.sort>=0])
			j2<-min(which(tmp<delta[i]))-1
			f2<-if(length(vec.pos)==0) sum(1/probs[(m-j2+1):m]-1)/B
				else vec.pos[m-j2+1]
			mat.delta[i,5]<-z.sort[m-j2+1]
		}
		else{
			j2<-0
			f2<-0
			mat.delta[i,5]<-Inf
		}
		mat.delta[i,2]<-j1+j2
		mat.delta[i,3]<-min(1,p0*(f1+f2)/max(1,j1+j2))
	}
	mat.delta		
}

"compPermStat" <-
function(X,P,n.cat){
	if(missing(n.cat))
		n.cat<-max(X)
	n.class<-max(P)
	n.obs<-ncol(X)
	vecX<-vector("list",n.cat)
	vecP<-vector("list",n.class)
	for(i in 1:n.cat)
		vecX[[i]]<-X==i
	for(i in 1:n.class)
		vecP[[i]]<-P==i
	vecXs<-lapply(vecX,rowSums)
	vecRs<-table(P[1,])
	#vecS<-vector("list",n.cat*n.class)
	matStat<-matrix(0,nrow(X),nrow(P))
	if(n.class>2){
		for(i in 1:n.class){
			for(j in 1:n.cat){
				tmp<-vecX[[j]]%*%t(vecP[[i]])
				#tmp<-tmp*tmp
				tmp2<-vecXs[[j]]*vecRs[i]/n.obs
				#vecS[[(i-1)*n.cat+j]]<-tmp*tmp/tmp2
				matStat<-matStat+tmp*tmp/tmp2	
			}
		}
	}
	else{
		for(j in 1:n.cat){
			tmp<-vecX[[j]]%*%t(vecP[[1]])
			tmp2<-vecXs[[j]]*vecRs[1]/n.obs
			matStat<-matStat+tmp*tmp/tmp2
			tmp<-vecXs[[j]]-tmp
			tmp2<-vecXs[[j]]-tmp2
			matStat<-matStat+tmp*tmp/tmp2
		}
	}
	#matStat<-matrix(0,nrow(X),nrow(P))
	#for(i in 1:length(vecS))
	#	matStat<-matStat+vecS[[i]]
	matStat-n.obs
}

`compRatio` <-
function(center,succ,fail,df=5,z=NULL){ 
	require(splines)
	n<-succ+fail
	ids<-which(n>0)
	if(is.null(z))
		z<-center
	tmp<-ns.out<-ns(center[ids],df)
	class(tmp)<-"matrix"
	mat<-data.frame(s=succ[ids],f=fail[ids],tmp)
	glm.out<-glm(cbind(s,f)~.,data=mat,family=binomial)
	newz<-predict(ns.out,z)
    	class(newz)<-"matrix"
	probs<-predict(glm.out,data.frame(newz),type="response")
	B<-sum(fail)/sum(succ)
	r<-(1-probs)/(B*probs)
	return(list(z=z,probs=probs,ratio=r))
}

"computeContCols" <-
function(x,CL,vec.ncl,n.obs,pam=FALSE){
	x<-x%*%CL
	r<-rowSums(x)
	if(any(r==0))
		stop("All rows of data must have the same number of categories.")
	tmp<-r%*%t(vec.ncl)/n.obs
	#rowSums((x-tmp)^2/tmp)
	if(pam)
		return(list(N.obs=x,N.exp=tmp,mat.stat=x*x/tmp))
	rowSums(x*x/tmp)
}

`computeRS` <-
function(data,cl,type){
	n.row<-nrow(data)
	le.cl<-length(cl)
	max.cl<-max(cl)+1
	tmp<-.C("get_stat_num_denum",as.double(data),as.integer(n.row),as.integer(le.cl),
		as.integer(cl),as.double(.mt.naNUM),t.num=double(n.row),t.denum=double(n.row),
		as.character(c(type,"abs","y")),as.integer(max.cl),PACKAGE="multtest")
	return(list(r=tmp$t.num,s=tmp$t.denum))
}

d.null<-function(X,mat.samp,d,type.mt,s0,med=FALSE,n.subset=10){
	if(med)
		n.subset<-1
	n.samp<-nrow(mat.samp)
	seq.samp<-unique(c(seq(1,n.samp,n.subset),n.samp+1))
	n.int<-length(seq.samp)-1
	n.row<-nrow(X)
	d.mat<-mat.neg<-mat.pos<-matrix(0,n.row,n.int)
	le.cl<-ncol(mat.samp)/ifelse(type.mt=="pairt",2,1)
	if(type.mt=="pairt"){
		X<-X[,2*(1:le.cl)]-X[,2*(1:le.cl)-1]
		mat.samp<-mat.samp[,2*(1:le.cl)]
	}
	d.rank<-rank(-abs(d),ties="first")
	for(j in 1:n.int){
		tmp<-mat.samp[seq.samp[j]:(seq.samp[j+1]-1),]
		if(!is.matrix(tmp))
			tmp<-matrix(tmp,1)
		dperm.out<-build.dperm(X,tmp,type.mt,s0,n.row,le.cl)	
		d.mat[,j]<-rowSums(as.matrix(dperm.out))
		mat.pos[,j]<-rank(-c(dperm.out[dperm.out>=0],abs(d)),ties="first")[sum(dperm.out>=0)+(1:n.row)]-d.rank
		mat.neg[,j]<-rank(c(dperm.out[dperm.out<0],-abs(d)),ties="first")[sum(dperm.out<0)+(1:n.row)]-d.rank
	}
	B<-nrow(mat.samp)
	d.bar<-rowSums(as.matrix(d.mat))/B
	p.value<-(rowSums(as.matrix(mat.pos))+rowSums(as.matrix(mat.neg)))/(n.row*B)
	vec.false<-numeric(n.row)
	if(med){
		vec.false[d>=0]<-apply(as.matrix(mat.pos[d>=0,]),1,median)
		if(type.mt!="f")
			vec.false[d<0]<-apply(as.matrix(mat.neg[d<0,]),1,median)
	}
	else{
		vec.false[d>=0]<-rowSums(as.matrix(mat.pos[d>=0,]))/B
		if(type.mt!="f")
			vec.false[d<0]<-rowSums(as.matrix(mat.neg[d<0,]))/B
	}
	invisible(return(list(d.bar=d.bar,p.value=p.value,vec.false=vec.false,d.mat=d.mat,
		mat.pos=mat.pos,mat.neg=mat.neg)))
} 
 

.mt.BLIM<-2^30
.mt.naNUM<- -93074815


d.stat<-function(data,cl,var.equal=FALSE,B=100,med=FALSE,s0=NA,s.alpha=seq(0,1,0.05),
		include.zero=TRUE,n.subset=10,mat.samp=NULL,B.more=0.1,B.max=30000,
		gene.names=NULL,R.fold=1,R.unlog=TRUE,na.replace=TRUE,na.method="mean",rand=NA){
	data<-as.matrix(data)
	mode(data)<-"numeric"
	n.genes<-nrow(data)
	excluded.genes<-logical(n.genes)
	if(any(rowSums(is.na(data))>0)){
		na.out<-na.handling(data,na.replace=na.replace,na.method=na.method)
		data<-na.out$X
		excluded.genes[na.out$NA.genes]<-TRUE
		rm(na.out)
	}
	adjust.out<-adjust.for.mt(data,cl,var.equal=var.equal)
	X<-adjust.out$X
	cl.mt<-adjust.out$cl.mt
	type.mt<-adjust.out$type.mt
	msg<-adjust.out$msg
	if(type.mt=="f" || length(unique(cl))==1)
		R.fold<-0
	if(R.fold>0){
		mat.fold<-Rfold.cal(data,cl,unlog=R.unlog,R.fold=R.fold)
		n.fulfill<-sum(mat.fold[,2])
		if(n.fulfill==0)
			stop("None of the variables has a fold change larger than ",R.fold,".")
		if(n.fulfill<20)
			stop("Less than 20 variables have a fold change larger than ",R.fold,".")
		if(R.fold!=1)
			msg<-c(msg,paste("Number of variables having a fold change >=",R.fold,"or <=",
				round(1/R.fold,4),":",n.fulfill,"\n\n"))
		fold.out<-mat.fold[,1]
		if(n.fulfill<nrow(mat.fold)){
			fc.genes<-which(mat.fold[,2]==0)
			fold.out<-fold.out[-fc.genes]
			X<-X[-fc.genes,]
			excluded.genes[!excluded.genes][fc.genes]<-TRUE
		}
	}			
	else
		fold.out<-numeric(0)
	rm(data,adjust.out)
	mt1.out<-mt.teststat.num.denum(X,cl.mt,test=type.mt)
	r<-mt1.out$teststat.num
	s<-mt1.out$teststat.denum
	if(any(round(s,10)==0)){
		var0.genes<-which(round(s,10)==0)
		r<-r[-var0.genes]
		s<-s[-var0.genes]
		X<-X[-var0.genes,]
		if(R.fold>0)
			fold.out<-fold.out[-var0.genes]
		excluded.genes[!excluded.genes][var0.genes]<-TRUE
		warning("There are ",length(var0.genes)," variables with zero variance. These variables are removed,",
			"\n","and their d-values are set to NA.",call.=FALSE)
	}
	if(is.na(s0)){
		s0.out<-fudge2(r,s,alpha=s.alpha,include.zero=include.zero)
		s0<-s0.out$s.zero
		msg<-c(msg,s0.out$msg)
	}
	else
		msg<-c(msg,paste("s0 =",round(s0,4),"\n\n"))
	d<-r/(s+s0)
	#d.sort<-sort(d)
	mat.samp<-setup.mat.samp(cl.mt,type.mt,B=B,mat.samp=mat.samp,B.more=B.more,B.max=B.max,rand=rand)
	B.full<-round(ifelse(type.mt=="pairt",2^(length(cl.mt)/2),
		exp(lgamma(length(cl.mt)+1)-sum(lgamma(table(cl.mt)+1)))))
	msg<-c(msg,paste("Number of permutations:",nrow(mat.samp),
		ifelse(nrow(mat.samp)==B.full,"(complete permutation)",""),"\n\n"),
		paste(ifelse(med,"MEDIAN","MEAN"),"number of falsely called variables is computed.\n\n"))
	dnull.out<-d.null(X,mat.samp,d,type.mt,s0,med=med,n.subset=n.subset)
	d.bar<-dnull.out$d.bar
	if(nrow(X)==n.genes){
		p.value<-dnull.out$p.value
		vec.false<-dnull.out$vec.false
	}
	else{
		p.value<-vec.false<-d.new<-s.new<-rep(NA,n.genes)
		p.value[!excluded.genes]<-dnull.out$p.value
		vec.false[!excluded.genes]<-dnull.out$vec.false
		d.new[!excluded.genes]<-d
		s.new[!excluded.genes]<-s
		d<-d.new
		s<-s.new
		if(R.fold>0){
			f.new<-rep(NA,n.genes)
			f.new[!excluded.genes]<-fold.out
			fold.out<-f.new
		}
			
	}
	rm(dnull.out)
	if(!is.null(gene.names))
		names(d)<-names(p.value)<-names(vec.false)<-names(s)<-substring(gene.names,1,50)
	invisible(list(d=d,d.bar=d.bar,p.value=p.value,vec.false=vec.false,discrete=FALSE,s=s,s0=s0,
		mat.samp=mat.samp,msg=msg,fold=fold.out))	
}
 
 

delta.plot<-function(object,delta=NULL,helplines=FALSE){
	par(mfrow=c(1,2))
	mat.fdr<-if(is.null(delta)) object@mat.fdr 
		else stats.cal(object@d,object@d.bar,object@vec.false,object@p0,delta=delta)
	plot(mat.fdr[,"Delta"],100*mat.fdr[,"FDR"],main="Delta vs. FDR",xlab=expression(Delta),
		ylab="FDR (in %)",type="b")
	if(helplines){
        	segments(0,100*mat.fdr[,"FDR"],mat.fdr[,"Delta"],100*mat.fdr[,"FDR"],lty=2)
        	segments(mat.fdr[,"Delta"],100*mat.fdr[,"FDR"],mat.fdr[,"Delta"],
			-100*max(mat.fdr[,"FDR"]),lty=2)
    	}
	plot(mat.fdr[,"Delta"],mat.fdr[,"Called"],main="Delta vs. Significant Genes", 
        	xlab=expression(Delta),ylab="Number of Significant Genes",type="b")
    	if(helplines){
		segments(0,mat.fdr[,"Called"],mat.fdr[,"Delta"],mat.fdr[,"Called"],lty=2)
        	segments(mat.fdr[,"Delta"],mat.fdr[,"Called"],mat.fdr[,"Delta"],
			-max(mat.fdr[,"Called"]),lty=2)
	}
	par(mfrow=c(1,1))
}
 
 

denspr<-function(x,n.interval=NULL,df=7){
	require(splines)
	if(is.null(n.interval))
		n.interval<-floor(sqrt(length(x)))
	breaks<-seq(min(x),max(x),length=n.interval+1)
	valHist<-hist(x,breaks=breaks,plot=FALSE)
	center<-valHist$mids
	counts<-valHist$counts
	ids<-which(counts>0)
	center<-center[ids]
	tmp<-ns.out<-ns(center,df=df)
	class(tmp)<-"matrix"
	mat<-data.frame(Number=counts[ids],tmp)
	glm.out<-glm(Number~.,data=mat,family=poisson)
	scale<-sum(diff(breaks)*counts)
	newx<-predict(ns.out,x)
	class(newx)<-"matrix"
	preds<-predict(glm.out,data.frame(newx),type="response")
	preds/scale
}

	


`ebam` <-
function(x,cl,method=z.ebam,delta=.9,which.a0=NULL,p0=NA,
		p0.estimation=c("splines","interval","adhoc"),lambda=NULL,ncs.value="max",
		use.weights=FALSE,gene.names=dimnames(x)[[1]],...){
	xclass<-class(x)
	if(!xclass%in%c("FindA0","ExpressionSet","matrix","data.frame"))
		stop("x must be an object of class FindA0, ExpressionSet, matrix, or data.frame.")
	if(xclass=="FindA0"){
		out<-ebamA0(x,which.a0=which.a0)
		chip.name<-x@chip
	}
	else{
		if(missing(cl))
			stop("cl must be specified if x is a matrix, a data frame, ",
				"or an ExpressionSet object.")
		if(is(x,"ExpressionSet")){
			require(affy,quietly=TRUE)
			chip.name<-annotation(x)
			if(is.character(cl) & length(cl)<=2)
			cl<-pData(x)[,cl]
			x<-exprs(x)
		}
		else
			chip.name<-""
		if(is.factor(cl))
			cl<-as.character(cl)
		if(ncol(x)!=length(cl))
			stop("The number of columns of data must be equal to the length of cl.")
		FUN<-match.fun(method)
		out<-FUN(x,cl,...)
	}
	check.out<-checkFUNout(out)
	if(is.na(p0))
		p0<-pi0.est3(out,p0.estimation,exact=check.out$exact,lambda=lambda,
			ncs.value=ncs.value,use.weights=use.weights)
	if(!is.null(gene.names)){
		names(out$z)<-names(out$ratio)<-substring(gene.names,1,50)
		if(length(check.out$vec.pos)!=0)
			names(check.out$vec.pos)<-names(check.out$vec.neg)<-names(out$z)
	}
	posterior<-1-p0*out$ratio
	posterior[posterior<0]<-0
	mat<-compNumber(out$z,posterior,p0,check.out$B,delta=delta,vec.pos=check.out$vec.pos,
		vec.neg=check.out$vec.neg)	
	new("EBAM",z=out$z,posterior=posterior,p0=p0,local=p0*out$ratio,mat.fdr=mat,
		a0=check.out$a0,mat.samp=check.out$mat.samp,vec.pos=check.out$vec.pos,
		vec.neg=check.out$vec.neg,msg=check.out$msg,chip=chip.name)		
		
}

`ebam2excel` <-
function(object,delta,file,excel.version=1,n.digits=4,what="both",ll=FALSE,
		chip="",quote=FALSE){
	if(!is(object,"EBAM"))
		stop("object must be an object of class EBAM.")
	siggenes2excel(object,delta,file,excel.version=excel.version,n.digits=n.digits,
		what=what,ll=ll,chip=chip,quote=quote)
}

ebam2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,findA0=NULL,
		varName=NULL,ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),plotFindArgs=plotFindArguments(),
		bg.plot.adjust=FALSE,plotname=NULL,plotborder=0,tableborder=1,
		new.window=TRUE,...){
	if(!is(object,"EBAM"))
		stop("object must be an object of class EBAM.")
	if(!is.null(findA0) && !is(findA0,"FindA0"))
		stop("findA0 must be an object of class FindA0.")
	siggenes2html(object,delta,filename,addStats=addStats,addPlot=addPlot,addGenes=addGenes,
		findA0=findA0,varName=varName,ll=ll,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		chipname=chipname,cdfname=cdfname,n.digits=n.digits,bg.col=bg.col,
		text.col=text.col,link.col=link.col,plotArgs=plotArgs,plotFindArgs=plotFindArgs,
		bg.plot.adjust=bg.plot.adjust,plotname=plotname,plotborder=plotborder,
		tableborder=tableborder,new.window=new.window,...)
}

`ebamA0` <-
function(x,which.a0=NULL){
	if(is.null(which.a0)){
		a0<-x@suggested
		which.a0<-which(x@vec.a0==a0)
	}
	else{
		if(!which.a0%in%(1:length(x@vec.a0)))
			stop("which.a0 must be an integer between 1 and ",length(x@vec.a0),".")
		a0<-x@vec.a0[which.a0]
	}
	ratio<-compRatio(x@mat.center[,which.a0],x@mat.success[,which.a0],x@mat.failure[,which.a0],
		df=x@df.ratio,z=x@mat.z[,which.a0])$ratio
	list(z=x@mat.z[,which.a0],ratio=ratio,a0=a0,success=x@mat.success[,which.a0],
		failure=x@mat.failure[,which.a0],center=x@mat.center[,which.a0],
		mat.samp=x@mat.samp,msg=x@msg)
}

filterALL<-function(){
	pdat<-pData(ALL)
	subset<-intersect(grep("^B",as.character(pdat$BT)),
		which(pdat$mol %in% c("BCR/ABL","NEG")))
	eset<-ALL[,subset]
	require(genefilter)
	f1<-pOverA(0.25,log2(100))
	f2<-function(x) IQR(x)>0.5
	selected<-genefilter(eset,filterfun(f1,f2))
	esetSub<-eset[selected,]
	pdat<-pData(esetSub)
	esetSub$mol.biol<-as.character(esetSub$mol.biol)
	#newPData<-new("phenoData",pData=pdat,
	#	varLabels=phenoData(esetSub)@varLabels)
	#phenoData(esetSub)<-newPData
	esetSub
}

`find.a0` <-
function(data,cl,method=z.find,B=100,delta=0.9,quan.a0=(0:5)/5,include.zero=TRUE,
		gene.names=dimnames(data)[[1]],n.chunk=5,n.interval=139,df.ratio=NULL,
		p0.estimation=c("splines","adhoc","interval"),lambda=NULL,ncs.value="max",
		use.weights=FALSE,rand=NA,...){
	if(length(delta)>1)
		stop("For find.a0, only one value of delta is allowed.\n",
			"Use print to obtain the number of interesting genes for other values ",
			"of delta.")
	if(is(data,"ExpressionSet")){
		require(affy,quietly=TRUE)
		if(is.character(cl) & length(cl)<=2)
			cl<-pData(data)[,cl]
		data<-exprs(data)
		chip.name<-annotation(data)
	}
	else
		chip.name<-""
	if(length(cl)!=ncol(data))
		stop("The number of columns of data must be equal to the length of cl.")
	FUN<-match.fun(method)
	if(!is.na(rand))
		set.seed(rand)
	out<-FUN(data,cl,B=B,...)
	tmp.names<-names(out)
	if(any(!c("r","s","z.fun","mat.samp")%in%tmp.names))
		stop("The function specified by method must return a list consisting of\n",
			"objects called r, s, z.fun, and mat.samp.")
	r<-out$r
	s<-out$s
	z.fun<-out$z.fun
	mat.samp<-out$mat.samp
	n.genes<-length(r)
	z.norm<-if(is.null(out$z.norm)) qnorm(((1:n.genes)-0.375)/(n.genes+0.25))
		else out$z.norm
	if(substitute(method)=="z.find"){
		data<-out$data
		cl<-out$cl
	}
	msg<-if(any(names(out)=="msg")) out$msg  else "" 
	rm(out)
	vec.a0<-checkQuantiles(s,quan.a0=quan.a0,include.zero=include.zero)
	z.mat<-r/outer(s,vec.a0,"+")
	succ.out<-apply(z.mat,2,getSuccesses,n.interval=n.interval)
	tmp<-lapply(succ.out,function(x) x$success)
	mat.success<-matrix(unlist(tmp),nrow=n.interval)
	ints<-lapply(succ.out,function(x) x$interval)
	tmp<-lapply(succ.out,function(x) x$center)
	mat.center<-matrix(unlist(tmp),nrow=n.interval)
	succ.norm<-getSuccesses(z.norm,n.interval=n.interval)
	fail.mats<-compFailureMat2(data,mat.samp,z.mat,ints,z.norm,succ.norm$interval,z.fun,
		vec.a0,n.chunk=n.chunk,n.interval=n.interval)
	mat.fn<-fail.mats$mat.fail.norm
	n.a0<-length(vec.a0)
	mat.ratio<-matrix(0,n.genes,n.a0)
	if(is.null(df.ratio))
		df.ratio<-ifelse(any(r<0),5,3)
	for(i in 1:n.a0)
		mat.ratio[,i]<-compRatio(succ.norm$center,succ.norm$success,mat.fn[,i],
			df=df.ratio,z=z.norm)$ratio
	type.p0<-match.arg(p0.estimation)
	if(type.p0=="adhoc")
		vec.p0<-apply(mat.ratio,2,function(x) min(1/x))
	else
		vec.p0<-apply(mat.fn,2,pi0.est2,success=succ.norm$success,z=succ.norm$center,
			type=type.p0,lambda=lambda,ncs.value=ncs.value,use.weights=use.weights)
	mat.posterior<-1-rep(vec.p0,e=n.genes)*mat.ratio
	mat.posterior[mat.posterior<0]<-0
	names(vec.a0)<-c(if(vec.a0[1]==0) "-", quan.a0)
	a0.out<-makeA0mat(z.norm,mat.posterior,vec.p0,vec.a0,B,delta=delta)
	colnames(z.mat)<-colnames(mat.posterior)<-colnames(mat.success)<-
		colnames(mat.center)<-names(vec.p0)<-round(vec.a0,4)
	if(!is.null(gene.names))
		rownames(z.mat)<-rownames(mat.posterior)<-gene.names
	new("FindA0",mat.z=z.mat,mat.posterior=mat.posterior,mat.center=mat.center,
		mat.success=mat.success,mat.failure=fail.mats$mat.failure,
		z.norm=z.norm,p0=vec.p0,mat.a0=a0.out$tab,
		mat.samp=mat.samp,vec.a0=vec.a0,suggested=a0.out$suggest,delta=delta,
		df.ratio=df.ratio,msg=msg,chip=chip.name)
}

finda02html<-function(x,delta,file,plotArgs=plotFindArguments(),plotnames=NULL,
		bg.plot.adjust=FALSE,bg.col=col2hex("white"),n.digits=4,tableborder=1,
		plotborder=0){
	if(!is(x,"FindA0"))
		stop("findA0 must be an object of class FindA0.")
	if(length(plotnames)!=2)
		stop("plotnames must have length 2.")
	if(delta==x@delta){
		mat.a0<-x@mat.a0
		sugg<-x@suggested
	}
	else{
		tmp<-makeA0mat(x@z.norm,x@mat.posterior,x@p0,x@vec.a0,nrow(x@mat.samp),
			delta=delta)
		mat.a0<-tmp$tab
		sugg<-tmp$suggest
	}
	mat.a0[,"FDR"]<-formatSAM(mat.a0[,"FDR"],digits=n.digits)
	mat.a0<-format(mat.a0,digits=n.digits)
	cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
	cat("<h3>Specification of the Fudge Factor</h3>",sep="\n",file=file)
	codeA0<-make.tablecode(rownames(mat.a0),ll=FALSE,refseq=FALSE,symbol=FALSE,omim=FALSE,
		ug=FALSE,dataframe=mat.a0,new.window=FALSE,tableborder=tableborder,
		name1stcol="&nbsp;")
	cat("<p><font color=",bg.col," size=1> HALLO</font></p>","\n",sep="",file=file)
	cat("<style type=text/css>","p{ margin-top: 2px; margin-bottom: 2px; word-spacing: 1px}",
		"</style>",codeA0,sep="\n",file=file)
	if(!plotArgs$onlyTab){
		suf.plot<-unlist(strsplit(plotnames[1],"\\."))
		suf.plot<-suf.plot[length(suf.plot)]
		FUN<-match.fun(suf.plot)
		FUN(plotnames[2])
		if(bg.plot.adjust)
			par(bg=bg.col)
		else
			par(bg="white")
		plot(x,delta,logit=plotArgs$logit,pos.legend=plotArgs$pos.legend,
			legend.cex=plotArgs$legend.cex,col=plotArgs$col,main=plotArgs$main,
			xlab=plotArgs$xlab,ylab=plotArgs$ylab,only.a0=plotArgs$only.a0,
			lty=plotArgs$lty,lwd=plotArgs$lwd,y.intersp=plotArgs$y.intersp)
		dev.off()
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
		cat("<div style=\"text-align: center\"><img src=",
			if(plotnames[1]==plotnames[2]) "file:///",plotnames[1]," border=",
			plotborder,"></div>","\n",sep="",file=file)
		cat("<p><b>Suggested Choice:</b>"," a0 = ",round(sugg,n.digits),
			"   (Selection Criterion: Number of Genes with Posterior >= ",delta,")",
			"</p>","\n",sep="",file=file)
	}
	cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",file=file)
}
		
		

	formatSAM<-function(x,digits=3){
	x<-format.pval(x,digits=digits,eps=0)
	x[as.numeric(x)==0]<-0
	x
}

pretty.mat.fdr<-function(x,digits=3){
	x<-as.data.frame(x)
	x[,"FDR"]<-formatSAM(x[,"FDR"],digits=digits)
	if(any(colnames(x)=="False")){
		x[,"False"]<-sapply(x[,"False"], 
			function(x) if(x>10^-digits) round(x,digits) else format(x,digits=3))
		idnr<-c(1,3,5)
	}
	else
		idnr<-c(1,3)
	x[,-idnr]<-round(x[,-idnr],digits)
	x
}

pretty.mat.sig<-function(x,digits=3){
	if(any(colnames(x)=="rawp")){
		x[,"rawp"]<-formatSAM(x[,"rawp"],digits=digits)
		x[,"q.value"]<-formatSAM(x[,"q.value"],digits=digits)
	}
	else
		x[,"local.fdr"]<-formatSAM(x[,"local.fdr"],digits=digits)
	format(x,digits=digits)
}

fudge2<-function(r,s,alpha=seq(0,1,0.05),include.zero=TRUE){
	if (max(alpha) > 1 || min(alpha) < 0) 
        	stop("alpha has to be between 0 and 1")
    	if (any(round(100 * alpha, 10) != round(100 * alpha, 0))) {
        	warning("At least one alpha is not a percentile. Only the first two decimal digits",
			 " are retained.")
        alpha <- signif(alpha, 2)
    	}
	if(length(alpha)==1){
		s.zero<-quantile(s,alpha)
		msg<-paste("s0 =",round(s.zero,4)," (The",100*alpha,"% quantile of the s values.) \n \n")
		invisible(return(list(s.zero=s.zero,alpha.hat=alpha,vec.cv=NULL,msg=msg)))
	}
	fudge.quan<-quantile(s,alpha)
	if(include.zero)
		fudge.quan<-c(0,fudge.quan)
	n.alpha<-length(fudge.quan)
	d.mat<-r/outer(s,fudge.quan,"+")
	n.uni.s <- length(unique(s))
	if(n.uni.s<25)
		stop("For the computation of the fugde factor,","\n",
			"there should be at least 25 genes with differing standard deviations.")
    	n.int <- ifelse(n.uni.s > 500, 101, floor(n.uni.s/5))
    	quan <- quantile(s, seq(0, 1, le = n.int))
	quan<-unique(quan)
	n.int<-length(quan)
	int.s<-as.numeric(cut(s,quan,include.lowest=TRUE,right=FALSE))
	mad.mat<-matrix(0,n.int-1,ncol(d.mat))
	for(i in 1:(n.int-1)){
		mad.mat[i,]<-apply(as.matrix(d.mat[which(int.s==i),]),2,mad)   # added as.matrix
		# med.s<-which(d.mat[which(int.s==i),1]==median(d.mat[which(int.s==i),1]))
		# mad.mat[which(int.s==i),]<-d.mat[med.s,]
	}
	cv<-function(x){
		sd(x)/mean(x)
	}
	vec.cv<-apply(mad.mat,2,cv)
	which.min<-which(vec.cv==min(vec.cv), arr.ind=FALSE)
	cat(which.min)
	if (length(which.min) == 0) {
		msg<-"s0 = 0 \n \n"
		s.zero <- 0
		invisible(return(list(s.zero=s.zero,vec.cv=vec.cv,msg=msg)))
	}
	if(include.zero & which.min==1){
		msg<-"s0 = 0 \n \n"
		s.zero <- 0
		invisible(return(list(s.zero=s.zero,vec.cv=vec.cv,msg=msg)))
	}
	s.zero<-fudge.quan[which.min]

	if(include.zero)
		which.min<-which.min-1
	alpha.hat<-alpha[which.min]
        msg<-paste("s0 =", round(s.zero, 4), " (The", 100 * alpha.hat, 
            "% quantile of the s values.)", "\n", "\n")
	invisible(return(list(alpha.hat=alpha.hat,s.zero=s.zero,vec.cv=vec.cv,msg=msg)))
}



 
 

`getFailure` <-
function(z.perm,z,interval,z.norm=NULL,z.range=NULL,n.interval=139){
	if(is.null(z.norm))
		z.perm<-truncZ(z.perm,z.range[1],z.range[2])
	else{
		z.sort<-sort(z)
		z.perm<-approx(z.sort,z.norm,z.perm,rule=2)$y
	}
	bin<-cut(z.perm,interval,include.lowest=TRUE)
	tabulate(bin,n.interval)
}

`getSuccesses` <-
function(z,n.interval=139){
	min.int<-floor(100*min(z))/100
	max.int<-ceiling(100*max(z))/100
	interval<-seq(min.int,max.int,length=n.interval+1)
	center<-(interval[2]-interval[1])/2 + interval[-length(interval)]
	bin<-cut(z,interval,include.lowest=TRUE)
	success<-tabulate(bin,n.interval)
	return(list(success=success,interval=interval,center=center))
}

getTD4Affy<-function(genenames,cdfname){
	tmp<-paste("https://www.affymetrix.com/analysis/netaffx/fullrecord.affx?pk=",
		cdfname,":",genenames,sep="")
	td<-paste("<TD><A HREF=\"",tmp,"\">",genenames,"</A></TD>",sep="")
	td
}



`help.ebam` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","summary","plot"))
		stop("'method' must be either print, summary, or plot.")
	fp<-.find.package("siggenes")
	fp<-file.path(fp,"doc",paste(method,"ebam.html",sep="."))
	browseURL(fp)
}

`help.finda0` <-
function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method%in%c("print","plot"))
		stop("'method' must be either print or plot.")
	fp<-.find.package("siggenes")
	fp<-file.path(fp,"doc",paste(method,"finda0.html",sep="."))
	browseURL(fp)
}

help.sam<-function(method){
	if(!is.character(method))
		method<-as.character(match.call(method)[[2]])
	if(!method %in% c("print","summary","plot","identify"))
		stop("'method' must be either print, summary, plot or identify.")
	fp<-.find.package("siggenes")
	fp<-file.path(fp,"doc",paste(method,"sam.html",sep="."))
	browseURL(fp)
}



link.genes<-function(genenames,filename,ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,
		ug=TRUE,chipname="",cdfname=NULL,dataframe=NULL,title=NULL,
		bg.col="white",text.col="black",link.col="blue",tableborder=1,
		new.window=TRUE){
	tr<-make.tablecode(genenames,ll=ll,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		chipname=chipname,cdfname=cdfname,dataframe=dataframe,
		tableborder=tableborder,new.window=new.window)
	suffix<-tolower(substring(filename,nchar(filename)-4,nchar(filename)))
	if(suffix!=".html"){
		filename<-paste(filename,"html",sep=".")
		warning("Since the suffix of 'filename' is not 'html' '.html' is added",
			" to 'filename'.",call.=FALSE)
	}
	bg.col<-col2hex(bg.col)
	text.col<-col2hex(text.col)
	link.col<-col2hex(link.col)
	if(is.null(title))
		title<-"Links for a Set of Genes to Public Repositories"
	outfile<-file(filename,"w")
	cat("<html>","<head>","<title>Links to Public Repositories</title>","</head>",
		paste("<body bgcolor=",bg.col," text=",text.col," link=",link.col,">",
		sep=""),
		paste("<h1 align=center>",title,"</h1>",sep=""),
		paste("<style type=text/css>",sep=""),
		"p{ margin-top: 2px; margin-bottom: 2px;}",
		"</style>",
		tr,"</body>","</html>",sep="\n",file=outfile)
	close(outfile)
}




link.siggenes<-function(object,delta,filename,gene.names=NULL,addDataFrame=TRUE,
		ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,n.digits=3,title=NULL,bg.col="white",text.col="black",
		link.col="blue",tableborder=1,new.window=TRUE){
	if(any(c(ll,refseq,symbol,omim,ug)))
		chipname<-check.chipname(chipname,object@chip,cdfname)
	genenames<-list.siggenes(object,delta,gene.names=gene.names)
	dataframe<-if(addDataFrame) pretty.mat.sig(summary(object,delta,what="genes")@mat.sig[,-1],
		digits=n.digits) else NULL
	link.genes(genenames,filename,ll=ll,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		chipname=chipname,cdfname=cdfname,dataframe=dataframe,
		title=title,bg.col=bg.col,text.col=text.col,link.col=link.col,
		tableborder=tableborder,new.window=new.window)
}
	


		`list.siggenes` <-
function(object,delta,file="",gene.names=NULL,order=TRUE,text=NULL,append=FALSE){
	if(!class(object)%in%c("SAM","EBAM"))
		stop("object must be either a SAM or an EBAM object.")
	if(length(delta)!=1)
		stop("delta must have length 1.")
	if(delta<=0)
		stop("delta must be larger than 0.")
	if(is(object,"EBAM") & delta>1)
		stop("delta must be smaller than 1.")
	d<-if(is(object,"SAM")) object@d else object@z
	if(is.null(gene.names))
		gene.names<-names(d)
	if(is.null(gene.names))
		stop("gene.names must be specified.")
	if(length(gene.names)!=length(d))
		stop("The length of gene.names differs from the number of genes.")
	if(!is.null(names(d)) & any(gene.names!=names(d)))
		stop("Some of the gene.names differ from the gene names of the ",
			ifelse(is(object,"SAM"),"SAM","EBAM")," object.")
	rsg<-summary(object,delta=delta,what="stats")@row.sig.genes
	if(order)
		rsg<-rsg[rev(order(abs(d[rsg])))]
	sig.genes<-gene.names[rsg]
	if(file!=""){
		if(is.null(text))
			text<-paste(object@msg[1],"Delta = ",delta,"\n\n",sep="")
		if(text=="")
			text<-NULL
		cat(text,if(!is.null(text)) "\n",file=file,append=append,sep="")
		write(sig.genes,file=file,append=TRUE)
	}
	else
		sig.genes
}

make.tablecode<-function(genenames,ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,
		chipname="",cdfname=NULL,dataframe=NULL,new.window=TRUE,tableborder=1,
		name1stcol="Name"){
	require(annotate)
	if(any(c(ll,refseq,symbol,omim,ug))){
		if(chipname=="" & is.null(cdfname))
			stop("Either 'chipname' or 'cdfname' must be specified.")
	}
	if(!is.null(cdfname)){
		require(affy)
		clean<-cleancdfname(cdfname,addcdf=FALSE)
		if(chipname=="")
			chipname<-clean
		if(clean!=chipname)
			stop("'chipname' and 'cdfname' do not specify the same chip.")
		tmp<-new("AffyBatch",cdfName=cdfname,annotation=clean)
		gN<-geneNames(tmp)
		if(!all(genenames%in%gN))
			stop("Some of the 'genenames' do not specify genes of the ",
				cdfname," chip.")
		tr<-getTD4Affy(genenames,cdfname)
	}
	else
		tr<-paste("<TD>",genenames,"</TD>",sep="")
	tr<-paste(tr,"\n")
	th<-name1stcol
	if(symbol){
		sym<-getSYMBOL(genenames,chipname)
		sym[is.na(sym)]<-"&nbsp;"
		sym.cols<-paste("<TD>",sym,"</TD>\n",sep="")
		th<-c(th,"Symbol")
		tr<-paste(tr,sym.cols)
	}
	if(ll){
		LL<-getLL(genenames,chipname)
		LL[is.na(LL)]<-"&nbsp;"
		LL.cols<-getTDRows(LL,"ll")
		th<-c(th,"LocusLink")
		tr<-paste(tr,LL.cols,"\n")
	}
	if(refseq){
		tmp.refseq<-lookUp(genenames,chipname,"REFSEQ")
		nm.select<-function(x) x[substring(x,1,2)%in%c("NM","nm")]
		RefSeq<-lapply(tmp.refseq,nm.select)
		RefSeq[lapply(RefSeq,length)==0]<-"&nbsp;"
		refseq.cols<-getTDRows(RefSeq,"GB")
		th<-c(th,"RefSeq")
		tr<-paste(tr,refseq.cols,"\n")
	}
	any.na<-function(x){any(is.na(x))}
	if(omim){
		OMIM<-lookUp(genenames,chipname,"OMIM")
		OMIM[unlist(lapply(OMIM,any.na))]<-"&nbsp;"
		OMIM.cols<-getTDRows(OMIM,"omim")
		th<-c(th,"OMIM")
		tr<-paste(tr,OMIM.cols,"\n")
	}
	if(ug){
		UG<-lookUp(genenames,chipname,"UNIGENE")
		UG[lapply(UG,length)==0]<-"&nbsp;"
		# Notloesung:
		UG<-unlist(lapply(UG,function(x) x[1]))
		UG.cols<-getTDRows(UG,"ug")
		th<-c(th,"UniGene")
		tr<-paste(tr,UG.cols,"\n")
	}
	if(new.window)
		tr<-add.target2href(tr)
	if(!is.null(dataframe)){
		if(!is.data.frame(dataframe))
			stop("'dataframe' must be a data.frame.")
		if(nrow(dataframe)!=length(genenames))
			stop("'The number of rows of 'data.frame' differ from the length",
				" of 'genenames'.")
		th<-c(th,colnames(dataframe))
		for(i in 1:ncol(dataframe)){
			tmp<-paste("<TD align=\"right\">",dataframe[,i],"</TD>",sep="")
			tr<-paste(tr,tmp,"\n")
		}
	}
	tr<-paste("<TR>\n",tr,"</TR>\n")
	trth<-paste(c("<TR>\n",paste("<TH>",th,"</TH>\n",sep=""),"</TR>"),collapse=" ")
	tr<-c(trth,"\n",tr)
	tr<-c(paste("<table border=",tableborder," align=center cellpadding=5px>\n",
		sep=""),tr,"</table>\n")
	tr	
}




`makeA0mat` <-
function(z.norm,mat.post,vec.p0,vec.a0,B,delta=0.9){
	n.a0<-length(vec.a0)
	mat<-matrix(0,n.a0,2)
	colnames(mat)<-c("Number","FDR")
	for(i in 1:n.a0)
		mat[i,]<-compNumber(z.norm,mat.post[,i],vec.p0[i],B,
			delta=delta)[,2:3]
	out<-data.frame(a0=round(vec.a0,4),Quantile=names(vec.a0),round(mat,4))
	rownames(out)<-1:nrow(out)
	ids<-which.max(mat[,1])
	list(tab=out,suggest=vec.a0[ids])
}

na.handling<-function(X,na.replace=TRUE,na.method="mean"){
	if(!is.matrix(X))
		stop("X must be a matrix.")
	vec.na<-rowSums(is.na(X))
	NA.genes<-which(vec.na>0)
	if(length(NA.genes)==0)
		return(list(X=X,NA.genes=NULL))
	n.col<-ncol(X)
	if(!na.replace){
		X<-X[-NA.genes,]
		warning("There are ",length(NA.genes)," genes with at least one missing expression value.",
			"\n","All these genes are removed, and their d-values are set to NA.",call.=FALSE)
	}
	else{
		X[NA.genes,]<-na.replace.cont(X[NA.genes,],na.method=na.method)
		warning("There are ",length(NA.genes)," genes with at least one missing expression value.",
			"\n","The NAs are replaced by the gene-wise ",na.method,".",call.=FALSE)
		if(any(vec.na>=n.col-1)){
			NA.genes0<-which(vec.na==n.col)
			NA.genes1<-which(vec.na==n.col-1)
			warning(length(NA.genes0)," of the ",length(NA.genes)," genes with at least one NA have no and ",
				length(NA.genes1)," have one non-missing expression value.","\n","All these ",
				length(c(NA.genes0,NA.genes1))," genes are removed, and their d-values are set to NA.",
				call.=FALSE)
			NA.genes<-c(NA.genes0,NA.genes1)
			X<-X[-c(NA.genes0,NA.genes1),]
		}
		else
			NA.genes<-NULL
	}
	list(X=X,NA.genes=NA.genes)
}
 
 

na.replace.cont<-function(X,na.method="mean"){
	if(!na.method%in%c("mean","median"))
		stop("na.method must be either mean or median.")
	FUN<-match.fun(na.method)
	if(is.vector(X))
		X<-replace(X,which(is.na(X)),FUN(X,na.rm=TRUE))
	else{
		for (i in 1:nrow(X)) 
			X[i,]<-replace(X[i,],which(is.na(X[i,])),FUN(X[i,],na.rm=TRUE))
	}
	return(X)
} 
 

pairt.cl.transform<-function(mat,n.row){
	if(!all(dim(mat)==c(n.row,2)))
		stop("cl must be either a vector or a ncol(data) x 2 matrix.")
	if(mode(as.matrix(mat))=="character"){
		vec1<-if(!all(substring(paste(mat[,1],collapse=""),1:sum(nchar(mat[,1])),
			1:sum(nchar(mat[,1])))%in%c(as.character(0:9),".","-"))) 
			as.character(mat[,1]) else as.numeric(mat[,1])
		vec2<-if(!all(substring(paste(mat[,2],collapse=""),1:sum(nchar(mat[,2])),
			1:sum(nchar(mat[,2])))%in%c(as.character(0:9),".","-"))) 
			as.character(mat[,2]) else as.numeric(mat[,2])
		mat<-data.frame(vec1,vec2)
	}
	tab1<-table(mat[,1])
	tab2<-table(mat[,2])
	n.tab1<-length(tab1)
	n.tab2<-length(tab2)
	if(n.tab1!=2 && n.tab2!=2)
		stop("One of the columns of cl must contain 2 different values.")
	if(n.tab1!=n.row/2 && n.tab2!=n.row/2)
		stop("One of the columns of cl must contain ",n.row/2," different values.")
	if(any(table(mat[,1],mat[,2])!=rep(1,n.row)))
		stop("There is something wrong with the cl matrix.")
	vec.sign<-if(n.tab1==2) mat[,1] else mat[,2]
	vec.pair<-if(n.tab1==n.row/2) mat[,1] else mat[,2]
	sort.sign<-sort(unique(vec.sign))
	sort.pair<-sort(unique(vec.pair))
	if(!all(sort.sign==c(-1,1))){
		warning("Expected values of one of the columns of cl are -1 and 1.","\n",
			as.character(sort.sign[1])," is thus set to -1, and ", 
			as.character(sort.sign[2])," to 1.",call.=FALSE)
		vec.sign<-ifelse(vec.sign==sort.sign[1],-1,1)
	}
	if(!all(sort.pair==1:(n.row/2))){
		mnc<-matrix(0,n.row/2,2)
		vec.tmp<-numeric(length(vec.pair))
                for(i in 1:(n.row/2)){
                	vec.tmp[vec.pair==sort.pair[i]]<-i
                	mnc[i,]<-c(as.character(sort.pair[i]),i)
                }
		vec.pair<-vec.tmp
		warning("Expected values of one of the columns are the integers between 1 and ",
			n.row/2,".","\n","The new class labels are thus ",paste(mnc[,2],"(was ",
			mnc[,1],ifelse(mnc[,2]!=mnc[nrow(mnc),2],"), ",")."),sep="",collapse=""),
			call.=FALSE)
	}
	cl<-as.numeric(vec.sign)*as.numeric(vec.pair)
	cl
}
	 
 

pairt.samp<-function(n.pair){
	mat<-matrix(0,2^n.pair,n.pair)
	for(i in 1:n.pair)
		mat[,i]<-rep(rep(c(1,-1),e=2^(n.pair-i)),2^(i-1))
	mat
} 
 

pairt.samp.transform<-function(mat){
	if(!all(mat%in%c(1,-1)))
		stop("Some of the values in mat.samp are neither 1 nor -1.")
	vec.samp<-as.vector(t(mat))
	n.samp<-length(vec.samp)
	vec.new<-numeric(2*n.samp)
	vec.new[2*(1:n.samp)-1]<-vec.samp>0
	vec.new[2*(1:n.samp)]<-vec.samp<0
	mat.samp<-matrix(vec.new,nrow(mat),2*ncol(mat),byrow=TRUE)
	mat.samp
}
	 
 

pi0.est<-function(p,lambda=seq(0,.95,.05),ncs.value="max",ncs.weights=NULL){
	if(any(lambda>=1) || any(lambda<0))
		stop("lambda has to be between 0 and 1.")
	if(!ncs.value%in%c("paper","max"))
		stop("ncs.value has to be either paper oder max.")
	le.la<-length(lambda)
	ncs.pred<-ifelse(ncs.value=="paper",1,max(lambda))
	vec.p0<-numeric(le.la)
	for(i in 1:le.la)
		vec.p0[i]<-mean(p>lambda[i])/(1-lambda[i])
	if(le.la==1)
		return(list(p0=min(vec.p0,1),spline.out=NULL))
	spline.out<-smooth.spline(lambda,vec.p0,w=ncs.weights,df=3)
	p0<-min(predict(spline.out,ncs.pred)$y,1)
	if(p0<=0){
		warning("The spline based estimation of pi0 results in a non-positive value of pi0.",
			"\n","Therefore, pi0 is estimated by using lambda = 0.5.",call.=FALSE)
		p0<-vec.p0[lambda==.5]
		spline.out<-NULL
	}
	return(list(p0=p0,spline.out=spline.out))
} 
 

`pi0.est2` <-
function(fail,success,z,type,lambda=NULL,ncs.value="max",use.weights=FALSE){
	if(!type%in%c("splines","interval"))
		stop("type must be either 'splines' or 'interval' in pi0.est2.")
	if(is.null(lambda)){
		lambda<-if(type=="splines") seq(0,0.95,0.05) else 0.5
	}
	w<-if(use.weights) 1-lambda else NULL
	z.order<-order(abs(z))
	z.fail<-cumsum(fail[z.order])
	z.success<-cumsum(success[z.order])
	n.fail<-sum(fail)
	B<-n.fail/sum(success)
	if(type=="interval"){
		num<-n.fail*(1-lambda)
		id<-which.min(abs(z.fail-num))
		pi0<-B*z.success[id]/(z.fail[id])
	}
	else{
		tmp<-(n.fail-z.fail[-length(z.fail)])/n.fail
		p<-rep(c(1,tmp),success[z.order])
		pi0<-pi0.est(p,lambda=lambda,ncs.value=ncs.value,ncs.weights=w)$p0
	}	
	pi0	
}

`pi0.est3` <-
function(fun.out,p0.estimation=c("splines","interval","adhoc"),exact=TRUE,lambda=NULL,
		ncs.value="max",use.weights=FALSE){
	type.p0<-match.arg(p0.estimation)
	if(type.p0=="adhoc")
		return(min(1/fun.out$ratio))
	if(!exact)
		return(pi0.est2(fun.out$failure,fun.out$success,fun.out$center,type.p0,
			lambda=lambda,ncs.value=ncs.value,use.weights=use.weights))
	if(is.null(lambda))
		lambda<-if(type.p0=="splines") seq(0,0.95,0.05) else 0.5
	w<-if(use.weights) 1-lambda else NULL
	p<-(fun.out$vec.pos+fun.out$vec.neg)/length(fun.out$vec.pos)
	pi0.est(p,lambda=lambda,ncs.value=ncs.value,ncs.weights=w)$p0
}

plotArguments<-function(pos.stats=NULL,sig.col=3,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,
		ylab=NULL,pty="s",lab=c(10,10,7),pch=NULL,sig.cex=1,stats.cex=0.8,
		y.intersp=1.3){
	list(pos.stats=pos.stats,sig.col=sig.col,xlim=xlim,ylim=ylim,main=main,
		xlab=xlab,ylab=ylab,pty=pty,lab=lab,pch=pch,sig.cex=sig.cex,
		stats.cex=stats.cex,y.intersp=y.intersp)
}


plotFindArguments<-function(onlyTab=FALSE,logit=TRUE,pos.legend=NULL,legend.cex=0.8,col=NULL,
		main=NULL,xlab=NULL,ylab=NULL,only.a0=FALSE,lty=1,lwd=1,y.intersp=1.1){
	list(onlyTab=onlyTab,logit=logit,pos.legend=pos.legend,legend.cex=legend.cex,col=col,
		main=main,xlab=xlab,ylab=ylab,only.a0=only.a0,lty=lty,lwd=lwd,
		y.intersp=y.intersp)
}

quantiles<-function(x,prob){
    if(any(prob>1 | prob<0))
        stop("Probabilities must be between 0 and 1")
    x.sort<-sort(x)
    prob.exclude<-which(prob>0 & prob<1)
    nprob<-length(x)*prob[prob.exclude]
    quan<-numeric(length(prob))
    if(any(prob==0))
        quan[which(prob==0)]<-min(x)
    if(any(prob==1))
        quan[which(prob==1)]<-max(x)
    int<-which(nprob==round(nprob))
    not.int<-which(nprob!=round(nprob))
    if(length(not.int)>0)
        quan[prob.exclude][not.int]<-x.sort[ceiling(nprob[not.int])]
    if(length(int)>0)
        quan[prob.exclude][int]<-.5*(x.sort[nprob[int]]+x.sort[nprob[int]+1])
    return(quan)
    }

    
qvalue.cal<-function(p,p0,version=1){
	p.na<-NULL
	p.names<-names(p)
	if(any(is.na(p))){
		p.na<-which(is.na(p))
		n.p<-length(p)
		p<-as.numeric(na.exclude(p))
	}
	if(!version%in%c(1,2))
		stop("version must be either 1 (for using the FDR) or 2 (pFDR).")
	m<-length(p)
	p.sort<-sort(p)
	qvalue<-p0*m*p.sort/(1:m)
	if(version==2)
		qvalue<-qvalue/(1-(1-ifelse(p.sort!=0,p.sort,min(p.sort[p.sort!=0],10^-15)))^m)
	qvalue[m]<-min(qvalue[m],1)
	for(i in (m-1):1)
		qvalue[i]<-min(qvalue[i],qvalue[i+1])
	qvalue[order(p)]<-qvalue
	if(!is.null(p.na)){
		qnew<-rep(NA,n.p)
		qnew[-p.na]<-qvalue
		qvalue<-qnew
	}
	names(qvalue)<-p.names
	qvalue
}
	
	 
 

sam<-function(data,cl,method=d.stat,delta=NULL,n.delta=10,p0=NA,lambda=seq(0,.95,.05),
		ncs.value="max",ncs.weights=NULL,gene.names=dimnames(data)[[1]],q.version=1,...){
	if(is(data,"exprSet") | is(data,"ExpressionSet")){
		require(affy)
		chip.name<-annotation(data)
		if(is.character(cl) & length(cl)<=2)
			cl<-pData(data)[,cl]
		data<-exprs(data)
	}
	else
		chip.name<-""
	if(is.factor(cl))
		cl<-as.character(cl)
	FUN<-match.fun(method)
	d.out<-FUN(data,cl,...)
	if(is.na(p0))
		p0<-pi0.est(na.exclude(d.out$p),lambda=lambda,ncs.value=ncs.value,
			ncs.weights=ncs.weights)$p0
	mat.fdr<-stats.cal(d.out$d,d.out$d.bar,d.out$vec.false,p0,delta=delta,le.delta=n.delta)
	if(!is.null(gene.names)){
		names(d.out$d)<-names(d.out$p.value)<-names(d.out$vec.false)<-substring(gene.names,1,50)
		if(length(d.out$s)==length(d.out$d))
			names(d.out$s)<-names(d.out$d)
	}
	if(q.version%in%c(1,2))
		q.value<-qvalue.cal(d.out$p.value,p0,version=q.version)
	else
		q.value<-numeric(0)
	new("SAM",d=d.out$d,d.bar=d.out$d.bar,vec.false=d.out$vec.false,p.value=d.out$p.value,
		s=d.out$s,s0=d.out$s0,mat.samp=d.out$mat.samp,p0=p0,mat.fdr=mat.fdr,q.value=q.value,
		fold=d.out$fold,msg=d.out$msg,chip=chip.name)
} 
 

sam.dstat<-function(...){
	cat("sam.dstat has been removed.\n",
		"Instead of sam.dstat(...), please use sam(...).","\n","\n",
	    	"See the help page of sam for general arguments for a SAM analysis,","\n",
		"and the help page of d.stat for further arguments for a SAM analysis","\n",
		"with (modified) t- or F-statistics.\n",sep="")
}

sam.plot2<-function(object,delta,pos.stats=NULL,sig.col=3,xlim=NULL,ylim=NULL,main=NULL,xlab=NULL,
		ylab=NULL,pty="s",lab=c(10,10,7),pch=NULL,sig.cex=1,...){
	if(is.null(pos.stats))
		pos.stats<-ifelse(all(object@d>=0,na.rm=TRUE),2,1)
	if(!pos.stats%in%0:2)
		stop("pos.stats must be either 0 (statistics are not displayed),\n",
			"1 (stats are displayed in the upper left of the plot), or 2 (lower right).")
	if(length(sig.col)==1)
		col.down<-col.up<-sig.col
	else{
		col.down<-sig.col[1]
		col.up<-sig.col[2]
	}
	d.sort<-sort(object@d)
	d.bar<-object@d.bar
	if(is.null(xlim))
		xlim<-c(min(d.sort,d.bar),max(d.sort,d.bar))
	if(is.null(ylim))
		ylim<-c(min(d.sort,d.bar),max(d.sort,d.bar))
	if(is.null(main))
		main<-paste("SAM Plot for Delta =",delta)
	if(is.null(xlab))
		xlab<-"Expected d(i) values"
	if(is.null(ylab))
		ylab<-"Observed d(i) values"
	par.store<-list(par()$pty,par()$lab)
	on.exit(par(pty=par.store[[1]],lab=par.store[[2]]))
	par(pty=pty,lab=lab)
	mat.fdr<-stats.cal(object@d,object@d.bar,object@vec.false,object@p0,delta=delta)
	d.up<-which(d.sort>=mat.fdr[,"cutup"])
	d.down<-which(d.sort<=mat.fdr[,"cutlow"])
	if(length(c(d.up,d.down))==0)
		plot(d.bar,d.sort,main=main,xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,pch=pch,...)
	else{
		plot(d.bar[-c(d.up,d.down)],d.sort[-c(d.up,d.down)],main=main,xlab=xlab,ylab=ylab,
			xlim=xlim,ylim=ylim,pch=pch,...)
		points(d.bar[d.up],d.sort[d.up],col=col.up,cex=sig.cex,pch=pch)
		points(d.bar[d.down],d.sort[d.down],col=col.down,cex=sig.cex,pch=pch)
	}
	abline(0,1)
	abline(delta,1,lty=2)
	abline(-delta,1,lty=2)
	abline(h=mat.fdr[,"cutup"],lty=5,cex=1.5)
	abline(h=mat.fdr[,"cutlow"],lty=5,cex=1.5)
	stats<-paste(c("cutlow:","cutup:","p0:","Significant:","False:","FDR:"),
		round(mat.fdr[1,c("cutlow","cutup","p0","Called","False","FDR")],3),sep="  ")
	if(pos.stats==1)
		text(rep(xlim[1],6),seq(ylim[2],ylim[2]-(ylim[2]-ylim[1])/4,le=6),stats,
			adj=0,cex=.75)
	if(pos.stats==2)
		text(rep(xlim[2]-(xlim[2]-xlim[1])/4,6),
			seq(ylim[1],ylim[1]+(ylim[2]-ylim[1])/4,le=6),stats,adj=0,cex=.75)
}

 
 

sam.snp<-function(...){
	cat("sam.snp has been removed.\n",
		"Instead of sam.snp(...), please use sam(..., method = cat.stat).","\n","\n",
		"See the help page of sam for general arguments for a SAM analysis,","\n",
		"and the help page of cat.stat for further arguments for a SAM analysis","\n",
		"of categorical data.\n",sep="")
}
 

sam.wilc<-function(...){
	cat("sam.wilc has been removed.\n",
		"Instead of sam.wilc(...), please use sam(..., method = wilc.stat).","\n","\n",
		"See the help page of sam for general arguments for a SAM analysis,","\n",
		"and the help page of wilc.stat for further arguments for a SAM analysis","\n",
		"with Wilcoxon rank sums.\n",sep="")
}
		 
 

`sam2excel` <-
function(object,delta,file,excel.version=1,n.digits=3,what="both",ll=FALSE,
		chip="",quote=FALSE){
	if(!is(object,"SAM"))
		stop("object must be an object of class SAM.")
	siggenes2excel(object,delta,file,excel.version=excel.version,n.digits=n.digits,
		what=what,ll=ll,chip=chip,quote=quote)
}

sam2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,varName=NULL,
		ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),bg.plot.adjust=FALSE,plotname=NULL,
		plotborder=0,tableborder=1,new.window=TRUE,...){
	if(!is(object,"SAM"))
		stop("object must be an object of class SAM.")
	if(length(delta)!=1)
		stop("delta must be a numeric value.")
	siggenes2html(object,delta,filename,addStats=addStats,addPlot=addPlot,addGenes=addGenes,
		varName=varName,ll=ll,refseq=refseq,symbol=symbol,omim=omim,ug=ug,
		chipname=chipname,cdfname=cdfname,n.digits=n.digits,bg.col=bg.col,
		text.col=text.col,link.col=link.col,plotArgs=plotArgs,
		bg.plot.adjust=bg.plot.adjust,plotname=plotname,plotborder=plotborder,
		tableborder=tableborder,new.window=new.window,...)
}




	setup.mat.samp<-function(cl,type.mt,B=1000,mat.samp=NULL,B.more=0.1,B.max=50000,rand=NA){
	if(!is.na(rand))
		set.seed(rand)
	n.cl<-length(cl)
	if(!is.null(mat.samp)){
		if(type.mt=="pairt")
			mat.samp<-pairt.samp.transform(mat.samp)
		if(ncol(mat.samp)!=n.cl)
			stop("The number of columns of mat.samp must be equal to the length of cl.")
		if(type.mt=="f")
			mat.samp<-mat.samp-1		
		a.out<-apply(mat.samp,1,function(a,y=cl) all(sort(a)==sort(y)))
		if(sum(a.out)!=nrow(mat.samp))
			stop("There is something wrong with mat.samp.")
		return(mat.samp)
	}
	if(B<0 || round(B)!=B)
		stop("B must be a positive integer.")
	if(B.more<0 || B.more>1)
		stop("B.more must be between 0 and 1.")
	if(B.max<0 || round(B.max)!=B.max)
		stop("B.max must be a positive integer.")
	B.full<-round(ifelse(type.mt=="pairt",2^(n.cl/2),exp(lgamma(n.cl+1)-sum(lgamma(table(cl)+1)))))
	if(B==0)
		B<-B.full
	B.large<-ceiling(B*(1+B.more))
	if(B.large>=B.full){
		mat.samp<-mt.sample.label(cl,test=type.mt,B=B.large+1)
		cat("\n")
		return(mat.samp)
	}
	if(B.max>=B.full){
		if(type.mt=="pairt"){
			tmp<-pairt.samp.transform(pairt.samp(n.cl/2)) 
	 		cat("We're doing",2^(n.cl/2),"complete permutations \n")
		}	
		else 
			tmp<-mt.sample.label(cl,test=type.mt,B=B.full+1)
		mat.samp<-tmp[sample(B.full,B),]
		cat("and randomly select",B,"of them.\n\n")
		return(mat.samp)
	}
	#if(balanced){
	#	if(!type.mt%in%c("t","t.equalvar"))
	#		stop("The balanced option is only available for the two class unpaired cases.")
	#	n.x<-sum(cl==0)
	#	if(n.x!=n.cl/2 || n.x/2!=round(n.x/2))
	#		stop("Balanced permutations are only possible if n1=n2 are even integer.")
	#	}
	if(type.mt=="pairt"){
		tmp<-matrix(sample(c(-1,1),B*n.cl/2,TRUE),B,n.cl/2)
		mat.samp<-pairt.samp.transform(tmp)
	}
	else{
		mat.samp<-matrix(0,B,n.cl)
		for(i in 1:B) mat.samp[i,]<-sample(cl)
	}
	return(mat.samp)
}
	 
 

`siggenes2excel` <-
function(object,delta,file,excel.version=1,n.digits=3,what="both",ll=FALSE,
		chip="",quote=FALSE){
	if(!excel.version%in%c(1,2))
		stop("'excel.version must be either 1 and 2.")
	sep<-ifelse(excel.version==1,",",";")
	dec<-ifelse(excel.version==1,".",",")
	suffix<-tolower(substring(file,nchar(file)-3,nchar(file)))
	if(suffix!=".csv"){
		file<-paste(file,"csv",sep=".")
		warning("Since the suffix of 'file' is not 'csv' this suffix is added ",
			"to 'file'.",call.=FALSE)
	}
	summary(object,delta,n.digits=n.digits,what=what,ll=ll,chip=chip,
		file=file,sep=sep,quote=quote,dec=dec)
}

siggenes2html<-function(object,delta,filename,addStats=TRUE,addPlot=TRUE,addGenes=TRUE,findA0=NULL,
		varName=NULL,ll=TRUE,refseq=TRUE,symbol=TRUE,omim=TRUE,ug=TRUE,chipname="",
		cdfname=NULL,n.digits=3,bg.col="white",text.col="black",link.col="blue",
		plotArgs=plotArguments(),plotFindArgs=plotFindArguments(),
		bg.plot.adjust=FALSE,plotname=NULL,plotborder=0,tableborder=1,
		new.window=TRUE,...){
	type<-class(object)
	isSAM<-type=="SAM"
	if(length(delta)!=1)
		stop("delta must be a numeric value.")
	if(delta<=0)
		stop("delta must be larger than 0.")
	if(!isSAM && delta>=1)
		stop("delta must be smaller than 1.")
	if(any(c(ll,refseq,symbol,omim,ug))){
		tmp<-ifelse(!isSAM,"z","d")
		if(is.null(names(slot(object,tmp)))){
			ll<-refseq<-symbol<-omim<-ug<-FALSE
			warning("Since no gene names are specified by 'object'",
				" 'll', 'refseq', 'symbol',\n","'omim' and 'ug' are set",
				" to FALSE.",call.=FALSE)
		}
		if(chipname=="" & object@chip=="" & is.null(cdfname)){
			ll<-refseq<-symbol<-omim<-ug<-FALSE
			warning("Since the chip type has been specified neither",
				" by the ",type," object\n","nor by 'chipname' or",
				" 'cdfname', 'll', 'refseq', 'symbol',\n",
				"'omim' and 'ug' are set to FALSE.",call.=FALSE)
		}
	}
	if(any(c(ll,refseq,symbol,omim,ug)))
		chipname<-check.chipname(chipname,object@chip,cdfname)
	suffix<-tolower(substring(filename,nchar(filename)-4,nchar(filename)))
	if(suffix!=".html"){
		filename<-paste(filename,"html",sep=".")
		warning("Since the suffix of 'filename' is not 'html' '.html' is added",
			" to 'filename'.",call.=FALSE)
	}
	if(is.null(plotname))
		plotname<-paste(type,"plot_for_",gsub(".html","",basename(filename)),
			".png",sep="")
	file.sep<-.Platform$file.sep
	tmp<-unlist(strsplit(filename,file.sep))
	path<-if(length(tmp)>1) paste(tmp[-length(tmp)],collapse=file.sep)
		else "."	
	bg.col<-col2hex(bg.col)
	text.col<-col2hex(text.col)
	link.col<-col2hex(link.col)
	sum.out<-summary(object,delta,n.digits=n.digits,chip=chipname)
	mat.fdr<-pretty.mat.fdr(sum.out@mat.fdr,digits=n.digits)
	mat.sig<-sum.out@mat.sig
	if(!is.null(mat.sig))
		mat.sig<-pretty.mat.sig(mat.sig,digits=n.digits)
	else{
		if(addGenes){
			addGenes<-FALSE
			warning("Since there are no significant genes, 'addGenes' is set",
				" to FALSE.",call.=FALSE)
		}
	}
	msg<-object@msg
	h2<-unlist(strsplit(msg[1],"\n"))[1]
	h2<-substring(h2,ifelse(isSAM,13,14),nchar(h2))
	isSNP<-h2==" for Categorical Data"
	if(is.null(varName))
		varName<-ifelse(isSNP,"SNPs","Genes")
	if(isSAM && isSNP){
		tmp.ids<-which(colnames(mat.sig)%in%c("stdev","R.fold"))
		mat.sig<-mat.sig[,-tmp.ids]
	}
	outfile<-file(filename,"w")
	cat("<html>","<head>",paste("<title>",type," Analysis</title>",sep=""),"</head>",
		paste("<body bgcolor=",bg.col," text=",text.col," link=",link.col,
		">",sep=""),paste("<h1 align=center>",type,"Analysis </h1>"),
		"<style type=text/css",
		"h3{ margin-top: 5px; margin-bottom: 2px; padding-left: 10px; text-indent: 10px;
			word-spacing: 1px}", 
		"li{ margin-top: 10px; padding-left: 0px; text-indent: 5px; word-spacing: 1px}",
		"ul{ padding-left: 50px}",	
		"</style>",
		paste("<h2 align=center>",h2,"</h2>"),
		sep="\n",file=outfile)
	if(addPlot | !is.null(findA0)){
		suf.plot<-unlist(strsplit(plotname,"\\."))
		suf.plot<-suf.plot[length(suf.plot)]
		if(!suf.plot%in%c("jpeg","png"))
			stop("'plotname' must be either a png or a jpeg file.")
		tmp.local<-length(unlist(strsplit(plotname,file.sep)))
		plotname2<-if(tmp.local==1) paste(path,plotname,sep=file.sep) else plotname
	}
	if(!is.null(findA0)){
		pn.finda0<-gsub(type,"finda0",c(plotname,plotname2))
		finda02html(findA0,delta,outfile,plotArgs=plotFindArgs,plotnames=pn.finda0,
			bg.plot.adjust=bg.plot.adjust,bg.col=bg.col,n.digits=n.digits,
			tableborder=tableborder,plotborder=plotborder)
	}
	if(addStats){
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat(paste("<h3>General Information",
			if(!isSAM & !is.null(findA0)) " for the Actual EBAM Analysis","</h3>",
			sep=""),sep="\n",file=outfile)
		cat("<ul>",sep="\n",file=outfile)
		s0<-slot(object,ifelse(isSAM,"s0","a0"))
		if(length(s0)==1){
			s0.msg<-if(isSAM) msg[substring(msg,1,2)=="s0"] else round(s0,n.digits)
			if(length(s0.msg)==1)
				cat("<li>Fudge Factor: ",s0.msg,"\n",sep="",file=outfile)
		}
		cat(paste("<li>Prior Probability: p0 = ",
			ifelse(isSAM,mat.fdr[,"p0"],round(object@p0,n.digits)),sep=""),
			paste("<li>Statistics for Delta = ",mat.fdr[,"Delta"],":",sep=""),
			"<ul type=\"disc\">",
			paste("<li>Number of Identified ",varName,": ",
			ifelse(isSAM,mat.fdr[,"Called"],mat.fdr[,"Number"]),sep=""),
			if(isSAM)  paste("<li>Number of Falsely Called ",varName,": ",mat.fdr[,"False"],sep=""),
			paste("<li>False Discovery Rate: ",mat.fdr[,"FDR"],sep=""),
			"</ul>","</ul>",sep="\n",file=outfile)
	}
	if(addPlot){
		FUN<-match.fun(suf.plot)
		FUN(plotname2)
		if(bg.plot.adjust)
			par(bg=bg.col)
		else
			par(bg="white")
		if(isSAM)
			plot(object,delta,pos.stats=plotArgs$pos.stats,sig.col=plotArgs$sig.col,
				xlim=plotArgs$xlim,ylim=plotArgs$ylim,main=plotArgs$main,
				xlab=plotArgs$xlab,ylab=plotArgs$ylab,pty=plotArgs$pty,
				lab=plotArgs$lab,pch=plotArgs$pch,sig.cex=plotArgs$sig.cex,...)
		else
			plot(object,delta,pos.stats=plotArgs$pos.stats,sig.col=plotArgs$sig.col,
				sig.cex=plotArgs$sig.cex,pch=plotArgs$pch,stats.cex=plotArgs$stats.cex,
				main=plotArgs$main,xlab=plotArgs$xlab,ylab=plotArgs$ylab,
				y.intersp=plotArgs$y.intersp,...)
		dev.off()
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<div style=\"text-align: center\"><img src=",
			if(plotname==plotname2) "file:///",plotname," border=",
			plotborder,"></div>","\n",sep="",file=outfile)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
	}
	if(addGenes){
		has.nonames<-all(rownames(mat.sig)==as.character(1:nrow(mat.sig)))
		if(has.nonames)
			tr<-make.tablecode(as.character(mat.sig[,"Row"]),ll=FALSE,
				refseq=FALSE,symbol=FALSE,omim=FALSE,ug=FALSE,
				chipname=chipname,dataframe=mat.sig[,-1],
				new.window=new.window,tableborder=tableborder,name1stcol="Row")
		else
			tr<-make.tablecode(rownames(mat.sig),ll=ll,refseq=refseq,
				symbol=symbol,omim=omim,ug=ug,chipname=chipname,
				cdfname=cdfname,dataframe=mat.sig[,-1],
				new.window=new.window,tableborder=tableborder)
		cat("<p><font color=",bg.col," size=2> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<h3 align=center>","Identified ",varName," (Using Delta = ",
			mat.fdr[,"Delta"],")</h3>","\n",sep="",file=outfile)
		cat("<p><font color=",bg.col," size=1> HALLO</font></p>","\n",sep="",
			file=outfile)
		cat("<style type=text/css>",
			"p{ margin-top: 2px; margin-bottom: 2px; word-spacing: 1px}",
			"</style>",tr,sep="\n",file=outfile)
		if(has.nonames)
			cat("<p><font color=",bg.col," size=1> HALLO</font></p>",
				"<p align=center><b>Annotation:</b> Since no gene names 
				have been specified the row numbers of the genes are used 
				as names.",sep="\n",file=outfile)
	}
	cat("<p><font color=",bg.col," size=5> HALLO</font></p>","\n",sep="",
			file=outfile)	
	cat("</body>","</html>",sep="\n",file=outfile)
	close(outfile)
	cat("Output is stored in ",filename,".\n",sep="")
	if(addPlot)
		cat("The ",type," Plot required by the html file is stored in ",plotname2,
			".\n",sep="")
	if(!is.null(findA0))
		cat("The FindA0 plot required by the html file is stored in ",pn.finda0[2],
			".\n",sep="")
} 

stats.cal<-function(d,d.bar,vec.false,p0,delta=NULL,le.delta=10){
	d.sort<-sort(d)
	d.diff<-d.sort-d.bar
	m<-length(d.diff)
	if(is.null(delta)){
		ra.ddiff<-range(abs(d.diff))
		delta<-round(seq(max(0.1,ra.ddiff[1]),max(1,ra.ddiff[2]),le=le.delta),1)
	}
	else{
		if(any(delta<=0))
			stop("delta must be larger than 0.")
		le.delta<-length(delta)
	}
	j0<-which(d.bar==min(d.bar[d.bar>=0]))[1]
	mat.fdr<-matrix(0,le.delta,9)
	dimnames(mat.fdr)<-list(1:le.delta,
		c("Delta","p0","False","Called","FDR","cutlow","cutup","j2","j1"))
	mat.fdr[,"Delta"]<-delta
	mat.fdr[,"p0"]<-p0
	vec.order<-as.numeric(na.exclude(vec.false[order(d)]))
	for(i in 1:le.delta){
		mat.fdr[i,"j1"]<-j1<-ifelse(any(d.diff[j0:m]>=delta[i]),
			j0-1+min(which(d.diff[j0:m]>=delta[i])),m+1)
		mat.fdr[i,"cutup"]<-ifelse(j1!=m+1,d.sort[j1],Inf)
		mat.fdr[i,"j2"]<-j2<-ifelse(any(d.diff[1:(j0-1)]<= -delta[i]) & j0!=1,
			max(which(d.diff[1:(j0-1)]<=-delta[i])),0)
		mat.fdr[i,"cutlow"]<-ifelse(j2!=0,d.sort[j2],-Inf)
		mat.fdr[i,"Called"]<-m-j1+1+j2
		mat.fdr[i,"False"]<-ifelse(j1==m+1,0,vec.order[j1])+ifelse(j2==0,0,vec.order[j2])
		mat.fdr[i,"FDR"]<-min(p0*mat.fdr[i,"False"]/max(mat.fdr[i,"Called"],1),1)
	}
	mat.fdr
}	
		 
 

require(methods)

ignoreThis=setClass("sumEBAM",representation(row.sig.genes="numeric",mat.fdr="matrix",
	mat.sig="data.frame",list.args="list"))

ignoreThis=setMethod("show","sumEBAM",function(object) print(object))

ignoreThis=setMethod("print","sumEBAM",
	function(x,varNames=NULL){
		list.args<-x@list.args
		rsg<-x@row.sig.genes
		mat.fdr<-x@mat.fdr
		mat.sig<-x@mat.sig
		nd<-list.args$n.digits
		file<-list.args$file
		msg<-list.args$msg
		if(is.null(varNames)){
			h2<-unlist(strsplit(msg[1],"\n"))[1]
			h2<-substring(h2,14,nchar(h2))
			varNames<-ifelse(h2==" for Categorical Data","SNPs","Genes")
		}	
		cat(msg,"Delta: ",mat.fdr[,"Delta"],"\n",sep="",file=file)
		if(length(list.args$a0)!=0)
			cat("a0: ",round(list.args$a0,nd),"\n",sep="",file=file,append=TRUE)
		cat("p0: ",round(list.args$p0,nd),"\n","Cutlow: ",round(mat.fdr[,"CL"],nd),"\n",
			"Cutup: ",round(mat.fdr[,"CU"],nd),"\n","Identified ",varNames,": ",
			round(mat.fdr[,"Number"],nd),"\n","Estimated FDR: ",
			formatSAM(mat.fdr[,"FDR"],nd),"\n",file=file,append=TRUE,sep="")
		if(list.args$what%in%c("genes","both") & length(rsg)!=0){
			if(all(row.names(mat.sig)!=as.character(1:nrow(mat.sig))))
				mat.sig<-data.frame(mat.sig,Name=rownames(mat.sig))
			mat.sig[,"local.fdr"]<-formatSAM(mat.sig[,"local.fdr"],digits=nd)
			mat.sig<-format(mat.sig,digits=nd)
			row.names(mat.sig)<-1:nrow(mat.sig)
			cat("\n\n","Identified ",varNames," (posterior >= ",
				mat.fdr[,"Delta"],"):\n\n",file=file,append=TRUE,sep="")
			if(file=="")
				print(mat.sig)
			else{
				write.table(t(dimnames(mat.sig)[[2]]),file=file,
					sep=list.args$sep,append=TRUE,row.names=FALSE,
					col.names=FALSE,quote=list.args$quote,dec=list.args$dec)
				write.table(mat.sig,file=file,sep=list.args$sep,append=TRUE,
					row.names=FALSE,col.names=FALSE,quote=list.args$quote,
					dec=list.args$dec)
			}
		}
		if(file!="")
			cat("Output is stored in", file,"\n")
	}
)


			


					require(methods)

ignoreThis=setClass("sumSAM",representation(row.sig.genes="numeric",mat.fdr="matrix",
	mat.sig="data.frame",list.args="list"))

ignoreThis=setMethod("show","sumSAM",
	function(object) print(object)
)

ignoreThis=setMethod("print","sumSAM",
	function(x,varNames=NULL){
		list.args<-x@list.args
		file<-list.args$file
		sep<-list.args$sep
		n.digits<-list.args$n.digits
		what<-list.args$what
		quote<-list.args$quote
		dec<-list.args$dec
		msg<-list.args$msg
		row.sig.genes<-x@row.sig.genes
		mat.fdr<-x@mat.fdr
		mat.sig<-x@mat.sig
		if(length(row.sig.genes)==0){
			cat(msg,file=file)
			mat.out<-pretty.mat.fdr(mat.fdr,digits=n.digits)
			if(file=="")
				print(mat.out)
			else{
				write.table(t(dimnames(mat.out)[[2]]),file=file,sep=sep,
					append=TRUE,row.names=FALSE,col.names=FALSE,
					quote=quote,dec=dec)
				write.table(mat.out,file=file,sep=sep,append=TRUE,
					row.names=FALSE,col.names=FALSE,quote=quote,
					dec=dec)
			}
		}
		else{	
			if(is.null(varNames)){
				h2<-unlist(strsplit(msg[1],"\n"))[1]
				h2<-substring(h2,13,nchar(h2))
				varNames<-ifelse(h2==" for Categorical Data","SNPs","Genes")
			}
			cat(msg[1],file=file)
			if(what%in%c("stats","both")){
				mat.out<-pretty.mat.fdr(mat.fdr,digits=n.digits)
				cat(msg[-1],file=file,append=TRUE)
				cat("Delta: ",mat.out[,"Delta"],"\n","cutlow: ",
					mat.out[,"cutlow"],"\n","cutup: ",mat.out[,"cutup"],
					"\n","p0: ",mat.out[,"p0"],"\n","Identified ",varNames,
					": ",
					mat.out[,"Called"],"\n","Falsely Called ",varNames,": ",
					mat.out[,"False"],"\n","FDR: ",mat.out[,"FDR"],
					ifelse(what=="both","\n\n\n","\n"),file=file,sep="",
					append=TRUE)
			}
			if(what%in%c("genes","both") & length(row.sig.genes)!=0){
				if(all(rownames(mat.sig)!=as.character(1:nrow(mat.sig))))
					mat.sig<-data.frame(mat.sig,Name=rownames(mat.sig))
				mat.out<-pretty.mat.sig(mat.sig,digits=n.digits)
				rownames(mat.out)<-1:nrow(mat.out)
				cat("Identified ",varNames," (using Delta = ",
					mat.fdr[,"Delta"],"):\n\n",file=file,append=TRUE,
					sep="")
				if(file=="")
					print(mat.out)
				else{
					write.table(t(dimnames(mat.out)[[2]]),file=file,
						sep=sep,append=TRUE,row.names=FALSE,
						col.names=FALSE,quote=quote,dec=dec)
					write.table(mat.out,file=file,sep=sep,append=TRUE,
						row.names=FALSE,col.names=FALSE,
						quote=quote,dec=dec)
				}
			}
		}
		if(file!="")
			cat("Output is stored in",file,"\n")
	}
)



			`truncZ` <-
function(z.perm,z.min,z.max){
	z.perm[z.perm<z.min]<-z.min
	z.perm[z.perm>z.max]<-z.max
	z.perm
}

wilc.ebam<-function(data,cl,approx50=TRUE,ties.method=c("min","random","max"),use.offset=TRUE,
		df.glm=5,rand=NA){
	require(splines)
	if(!is.na(rand))
		set.seed(rand)
	data<-as.matrix(data)
	mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("No missing values allowed.")
	adjust.out<-adjust.for.mt(data,cl,wilc=TRUE,eb=TRUE)
	X<-adjust.out$X
	cl.mt<-adjust.out$cl.mt
	type.mt<-adjust.out$type.mt
	msg<-adjust.out$msg
	rm(data,adjust.out)
	n.row<-nrow(X)
	n.cl<-length(cl.mt)
	ties.method<-match.arg(ties.method)
	if(type.mt=="t"){
		n0<-sum(cl.mt==0)
		n1<-sum(cl.mt==1)
		W.mean<-n1*(n.cl+1)/2
		W.min<-n1*(n1+1)/2
		W.max<-n1*(2*n0+n1+1)/2
		W.var<-n1*n0*(n.cl+1)/12
		X.rank<-matrix(0,n.row,n.cl)
		for(i in 1:n.row)
			X.rank[i,]<-rank(X[i,],ties.method=ties.method)
		W<-rowSums(X.rank[,cl.mt==1])
		W.norm<-(W-W.mean)/sqrt(W.var)
		if(n0<50 & n1<50)
			approx50<-FALSE
		if(!approx50){
			W.null<-dwilcox(W-W.min,n1,n0)
			p.value<-2*pwilcox(W.mean-abs(W-W.mean)-W.min,n1,n0)
			p.value[p.value>1]<-1
			tmp<-sort(unique(W-W.min))
			f.null<-dwilcox(tmp,n1,n0)
		}

		else{
			W.null<-dnorm(W.norm)
			p.value<-2*(1-pnorm(abs(W.norm)))
			f.null<-dnorm(sort(unique(W.norm)))
		}
	}
	if(type.mt=="pairt"){
		n.cl<-n.cl/2
		X<-X[,2*(1:n.cl)-1]-X[,2*(1:n.cl)]
		W.max<-n.cl*(n.cl+1)/2
		W.mean<-W.max/2
		W.var<-n.cl*(n.cl+1)*(2*n.cl+1)/24
		if(sum(X==0)>0){
			if(check.ties)
				warning("There are ",sum(X==0)," observation pairs having a difference of zero.",
					"\n","These differences are randomly set to either 1e-06 or -1e-06.",
					call.=FALSE)
			X[X==0]<-sample(c(1e-06,-1e-06),sum(X==0),rep=TRUE)
		}
		X.rank<-matrix(0,n.row,n.cl)
		for(i in 1:n.row)
			X.rank[i,]<-rank(abs(X[i,]),ties.method=ties.method)
		W<-rowSums(X.rank*(X>0))
		W.norm<-(W-W.mean)/sqrt(W.var)
		if(n.cl<50)
			approx50<-FALSE
		if(!approx50){
			W.null<-dsignrank(W,n.cl)
			p.value<-2+psignrank(W.mean-abs(W-W.mean),n.cl)	
			p.value[p.value>1]<-1
			f.null<-dsignrank(sort(unique(W)),n.cl)
		}
		else{
			W.null<-dnorm(W.norm)
			p.value<-2*(1-pnorm(abs(W.norm)))
			f.null<-dnorm(sort(unique(W.norm)))
		}
	}
	vec.pos<-n.row*p.value/2
	tabW<-table(W)
	valW<-as.numeric(names(tabW))
	offset.value<-if(use.offset) log(f.null) else rep(0,length(valW))
	glm.out<-glm(tabW~ns(valW,df.glm)+offset(offset.value),family=poisson)
	f.W<-glm.out$fitted/n.row
	W.fitted<-f.W[as.character(W)]
	list(z=W.norm,ratio=W.null/W.fitted,vec.pos=vec.pos,vec.neg=vec.pos,msg=msg)
}
	
	
		




	wilc.stat<-function(data,cl,gene.names=NULL,R.fold=1,R.unlog=TRUE,na.replace=TRUE,na.method="mean",
		approx50=TRUE,ties.method=c("min","random","max"),check.ties=FALSE,rand=NA){
	if(!is.na(rand))
		set.seed(rand)
	data<-as.matrix(data)
	mode(data)<-"numeric"
	n.genes<-nrow(data)
	excluded.genes<-logical(n.genes)
	if(any(rowSums(is.na(data))>0)){
		na.out<-na.handling(data,na.replace=na.replace,na.method=na.method)
		data<-na.out$X
		excluded.genes[na.out$NA.genes]<-TRUE
		rm(na.out)
	}
	adjust.out<-adjust.for.mt(data,cl,wilc=TRUE)
	X<-adjust.out$X
	cl.mt<-adjust.out$cl.mt
	type.mt<-adjust.out$type.mt
	msg<-adjust.out$msg
	if(length(unique(cl))==1)
		R.fold<-0
	if(R.fold>0){
		mat.fold<-Rfold.cal(data,cl,unlog=R.unlog,R.fold=R.fold)
		n.fulfill<-sum(mat.fold[,2])
		if(n.fulfill==0)
			stop("None of the genes has a fold change larger than ",R.fold,".")
		if(n.fulfill<20)
			stop("Less than 20 genes have a fold change larger than ",R.fold,".")
		if(R.fold!=1)
			msg<-paste("Number of genes having a fold change >=",R.fold,"or <=",
				round(1/R.fold,4),":",n.fulfill,"\n\n")
		fold.out<-mat.fold[,1]
		if(n.fulfill<nrow(mat.fold)){
			fc.genes<-which(mat.fold[,2]==0)
			fold.out<-fold.out[-fc.genes]
			X<-X[-fc.genes,]
			excluded.genes[!excluded.genes][fc.genes]<-TRUE
		}
	}
	else
		fold.out<-numeric(0)
	rm(data,adjust.out)
	n.cl<-length(cl.mt)
	x2<-function(x){
		x*x
	}
	if(type.mt=="pairt"){
		n.cl<-n.cl/2
		X<-X[,2*(1:n.cl)-1]-X[,2*(1:n.cl)]
		varX<-rowSums(x2(X-rowMeans(X)))/(n.cl-1)
	}
	else{
		varX<-rowSums(x2(X[,cl.mt==0]-rowMeans(X[,cl.mt==0])))+
			rowSums(x2(X[,cl.mt==1]-rowMeans(X[,cl.mt==1])))
		varX<-varX/(n.cl-1)
	}
	if(sum(varX<10^-10)){
		var0.genes<-which(varX<10^-10)
		X<-X[-var0.genes,]
		if(R.fold>0)
			fold.out<-fold.out[-var0.genes]
		excluded.genes[!excluded.genes][var0.genes]<-TRUE
		warning("There are ",length(var0.genes)," genes with zero variance. These genes are removed,",
			"\n","and their d-values are set to NA.",call.=FALSE)
	}
	n.row<-nrow(X)
	if(check.ties){
		n.ties<-numeric(n.row)
		for(i in 1:n.row)
			n.ties[i]<-n.cl-length(unique(X[i,]))
		if(any(n.ties>0))
			warning("There are ",sum(n.ties>0)," genes with ties. The ranks of these ties",
				" are randomly assigned.",call.=FALSE)
	}
	ties.method<-match.arg(ties.method)
	if(type.mt=="t"){
		n0<-sum(cl.mt==0)
		n1<-sum(cl.mt==1)
		W.mean<-n1*(n.cl+1)/2
		W.min<-n1*(n1+1)/2
		W.max<-n1*(2*n0+n1+1)/2
		W.var<-n1*n0*(n.cl+1)/12
		X.rank<-matrix(0,n.row,n.cl)
		for(i in 1:n.row)
			X.rank[i,]<-rank(X[i,],ties.method=ties.method)
		W<-rowSums(X.rank[,cl.mt==1])
		W.norm<-(W-W.mean)/sqrt(W.var)
		if(n0<50 & n1<50)
			approx50<-FALSE
		if(!approx50){
			W.exp<-W.min+qwilcox(((1:n.row)-.5)/n.row,n1,n0)
			W.exp<-(W.exp-W.mean)/sqrt(W.var)
			p.value<-2*pwilcox(W.mean-abs(W-W.mean)-W.min,n1,n0)
			p.value[p.value>1]<-1
		}
		else{
			W.exp<-qnorm(((1:n.row)-.5)/n.row)
			p.value<-2*(1-pnorm(abs(W.norm)))
		}
	}
	if(type.mt=="pairt"){
		W.max<-n.cl*(n.cl+1)/2
		W.mean<-W.max/2
		W.var<-n.cl*(n.cl+1)*(2*n.cl+1)/24
		if(sum(X==0)>0){
			if(check.ties)
				warning("There are ",sum(X==0)," observation pairs having a difference of zero.",
					"\n","These differences are randomly set to either 1e-06 or -1e-06.",
					call.=FALSE)
			X[X==0]<-sample(c(1e-06,-1e-06),sum(X==0),rep=TRUE)
		}
		X.rank<-matrix(0,n.row,n.cl)
		for(i in 1:n.row)
			X.rank[i,]<-rank(abs(X[i,]),ties.method=ties.method)
		W<-rowSums(X.rank*(X>0))
		W.norm<-(W-W.mean)/sqrt(W.var)
		if(n.cl<50)
			approx50<-FALSE
		if(!approx50){
			W.exp<-qsignrank(((1:n.row)-.5)/n.row,n.cl)
			W.exp<-(W.exp-W.mean)/sqrt(W.var)
			p.value<-2*psignrank(W.mean-abs(W-W.mean),n.cl)
			p.value[p.value>1]<-1
		}
		else{
			W.exp<-qnorm(((1:n.row)-.5)/n.row)
			p.value<-2*(1-pnorm(abs(W.norm)))
		}
	}
	vec.false<-n.row*p.value/2
	if(n.row!=n.genes){
		W.new<-vec.new<-p.new<-rep(NA,n.genes)
		W.new[!excluded.genes]<-W.norm
		vec.new[!excluded.genes]<-vec.false
		p.new[!excluded.genes]<-p.value
		W.norm<-W.new
		vec.false<-vec.new
		p.value<-p.new
		if(R.fold>0){
			f.new<-rep(NA,n.genes)
			f.new[!excluded.genes]<-fold.out
			fold.out<-f.new
		}
	}
	if(!is.null(gene.names))
		names(W.norm)<-names(p.value)<-names(vec.false)<-substring(gene.names,1,50)
	structure(list(d=W.norm,d.bar=W.exp,p.value=p.value,vec.false=vec.false,
		discrete=ifelse(approx50,FALSE,TRUE),s=numeric(0),s0=numeric(0),
		mat.samp=matrix(numeric(0)),msg=msg,fold=fold.out))
}
	 
 

`z.ebam` <-
function(data,cl,a0=NULL,quan.a0=NULL,B=100,var.equal=FALSE,B.more=0.1,B.max=30000,
		n.subset=10,fast=FALSE,n.interval=139,df.ratio=NULL,rand=NA){
	if(!is.na(rand))
		set.seed(rand)
	out<-z.find(data,cl,B=B,var.equal=var.equal,B.more=B.more,B.max=B.max)
	a0<-checkA0(out$s,a0=a0,quan.a0=quan.a0)
	z<-out$r/(a0+out$s)
	data<-out$data
	mat.samp<-out$mat.samp
	msg<-out$msg
	type.mt<-out$type.mt
	out<-getSuccesses(z,n.interval=n.interval)
	fail.out<-compFailure(data,mat.samp,z,out$interval,a0=a0,type.mt=type.mt,
		n.subset=n.subset,fast=fast)
	if(is.null(df.ratio))
		df.ratio<-ifelse(any(z<0),5,3)
	ratio<-compRatio(out$center,out$success,fail.out$vec.fail,df=df.ratio,z=z)$ratio
	if(fast)
		return(list(z=z,ratio=ratio,a0=a0,success=out$success,failure=fail.out$vec.fail,
			center=out$center,mat.samp=mat.samp,msg=msg))
	else
		return(list(z=z,ratio=ratio,a0=a0,vec.pos=fail.out$vec.pos/B,
			vec.neg=fail.out$vec.neg/B,mat.samp=mat.samp,msg=msg))		
}

`z.find` <-
function(data,cl,B=100,var.equal=FALSE,B.more=0.1,B.max=30000){
	data<-as.matrix(data)
	mode(data)<-"numeric"
	if(any(is.na(data)))
		stop("No missing values allowed.")
	adjust.out<-adjust.for.mt(data,cl,var.equal=var.equal,eb=TRUE)
	data<-adjust.out$X
	cl<-adjust.out$cl.mt
	msg<-adjust.out$msg
	type.mt<-adjust.out$type.mt
	z.fun<-switch(EXPR=type.mt,
		t=function(data,cl) computeRS(data,cl,"t"),
		t.equalvar=function(data,cl) computeRS(data,cl,"t.equalvar"),
		pairt=function(data,cl) computeRS(data,cl,"pairt"),
		f=function(data,cl) computeRS(data,cl,"f"))
	mt.out<-mt.teststat.num.denum(data,cl,test=type.mt)
	mat.samp<-setup.mat.samp(cl,type.mt,B=B,B.more=B.more,B.max=B.max)
	if(type.mt=="pairt"){
		le.cl<-length(cl)/2
		data<-data[,2*(1:le.cl)]-data[,2*(1:le.cl)-1]
		mat.samp<-mat.samp[,2*(1:le.cl)]
	}
	n.genes<-nrow(data)
	if(type.mt=="f"){
		n.classes<-length(unique(cl))
		z.norm<-qf(((1:n.genes)-0.5)/n.genes,n.classes-1,n.genes-n.classes)
	}
	else
		z.norm<-qnorm(((1:n.genes)-0.375)/(n.genes+0.25))
	return(list(r=mt.out$teststat.num,s=mt.out$teststat.denum,z.fun=z.fun,mat.samp=mat.samp,
		msg=msg,data=data,cl=cl,z.norm=z.norm,type.mt=type.mt))
}

.OnLoad <- function(libname, pkgname){
	require(methods)
}

.onAttach <- function(libname, pkgname) {
	message("\n","Latest changes in siggenes:","\n",
		"Version 1.9.1: Added a new, much faster version of sam.snp.","\n",
		"Version 1.9.9: Added a new, completely revised version of the EBAM functions.","\n",
		"Version 1.9.27: EBAM update completed. Added new vignette.",
		"\n")
    	if(interactive() && .Platform$OS.type == "windows" && .Platform$GUI ==  "Rgui")
        	addVigs2WinMenu("siggenes")
} 
 
#################################
# END INLINE
#################################

parseReadableFile <- function(string) {
  if (string == '-') return(file("stdin"))
  f <- try(file(string, open="r"))
  if ("try-error" %in% class(f)) return(NULL)
  return(f)
}

# sink(file=NULL, "type=message")

argspec <- c("sam.R - run Significance Analysis of Microarrays (SAM) on a
data matrix.  The top row should contain the classification, and every
subsequent row should have data.  The first two columns should have annotation
data (use the -ac option to change this).
mUsage:
    sam.R [options] -i input.tab > output.tab
Options:",
             "-ac i  annotation columns (default 2)")

main <- function(argv) {
  
  #o <- Rgetopt(argv=argv, argspec=argspec, defaults=list(ac=2))
  o <- getopt( c(
	'ac', 'a', 1, "integer",
	'input', 'i', 1, "character"
  )); 

  #if (missing(o$input)) usage("Must specify one file option", argspec)

  infile <- parseReadableFile(o$input)

  header <- strsplit(readLines(con=infile, n=1), "\t")[[1]]
  cl.cols <- 1:length(header) > o$ac
  cl.string <- header[cl.cols]
  cl <- as.integer(cl.string)
  if (any(cl != cl.string)) stop("Found non-integers in class labels")

  data.df <- read.delim(infile, header=FALSE, row.names=NULL,
                        stringsAsFactors=FALSE)
  # close(infile)
  data <- as.matrix(data.df[,cl.cols])

  if(!is.numeric(data)) stop("Non-numeric data in matrix")
  if(length(cl) != ncol(data)) stop("Header length does not match data")

  sink(stderr()) # sam() is noisy and prints to stdout
  s <- attributes(sam(data, cl))
  sink(NULL)
  warnings()


  r <- data.frame(data.df[,!cl.cols],
                  s[c("d", "vec.false", "q.value", "p.value", "s")])
  r <- r[order(r[,'d'], decreasing=TRUE),]
  colnames(r) <- c(header[!cl.cols], "Score(d)", "FalseCalls", "Q-value",
                   "P-value", "StdDev(s)")

  write.table(r, stdout(), quote=FALSE, sep="\t", row.names=FALSE)
}

main()

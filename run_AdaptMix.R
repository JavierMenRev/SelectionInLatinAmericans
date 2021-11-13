## NOTE: UNDER THE NORMAL-LIK MODEL, YOU CAN RESCALE THE LRT BY ADDING MORE DRIFT; IN PARTICULAR LRT(new) = [sqrt(drift.vec.old)/sqrt(drift.vec.new)]*LRT(old) 
## THUS HALVING THE LRT IS ACCOMPLISHED BY DOUBLING drift.vec.old

## Normal lik doesn't work particularly well with alt.drift=0 -- really favors SNPs near fixation

##########################
## WD
# setwd("/Volumes/Vanuatu2/Javier/UCL/Selection-6/")

############################
## INPUT:

temp=commandArgs(TRUE)
target.pop=as.character(temp[1])

add.on='BetaBinomInferredDriftPValII'
lik.calc='betabinom'           ## 'betabinom' or 'normal' likelihoods?
alt.drift=1
if (alt.drift==1) add.on=paste(add.on,"ALT",sep='')
if (lik.calc=='normal') add.on=paste(add.on,"NormalLik",sep='')
# out.file=paste("results/",target.pop,"_",data.type,"_",add.on,"_v2.txt",sep='')
out.file=paste("results/bbadapt/",target.pop,"_",add.on,".txt",sep='') ## I just want to check the drift -- should be the same as v2

pop.vec=target.pop   ## tested pops
surrogate.vec=c("Africans","Europeans","NativeAmericans")     ## surrogate pops  (!!!! MUST BE IN ORDER OF ADMIXTURE COLUMNS !!!!)

## get CP and ID files
id.file=paste0("./data/idfiles/",target.pop,"_CANDELA_NAM_EUR_AFR_CHB_GBR_geno_0.05_mind0.05_intersect_maf0.05.idfile")
geno.filePRE="./data/CPformatData/CANDELA_NAM_EUR_AFR_CHB_GBR_geno_0.05_mind0.05_intersect_chr"
geno.filePOST=".chromopainter.haps.gz"

nsnps.perrun=400000         ## TO CONTROL RAM
chromo.vec=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22")
# chromo.vec=c("21","22")
round.val=3     ## for AIC
round.val.sel=3     ## for sel-coeff
round.val.freq=4     ## for allele frequencies (observed and expected)

## ADMIXTURE DATA or SourceFIND DATA
admixture.props.file="./results/admixture_supervised/CANDELA_NAM_EUR_AFR_CHB_GBR_geno_0.05_mind0.05_intersect_maf0.05_pruned.3.Q"
admixture.inds.file="./results/admixture_supervised/CANDELA_NAM_EUR_AFR_CHB_GBR_geno_0.05_mind0.05_intersect_maf0.05_pruned.fam"

## REMOVE INDS WITH EURO ANCESTRY FROM NON-IBERIAN POPS BASED ON SOURCEFIND?: 
#nonIberian.pops=""
# nonIberian.pops=c("NorthWestEurope","NorthEastEurope")   ## ",Italy"
nonIberian.pops=c("")   ## ",Italy"
sourcefind.props.file="./data/globetrotter_sourcefind_combined.csv"
 
##############################
## PROGRAM:

Mb=10^6
options(scipen=999)

admixture.est.lm=function(x,admix.props,genotypes) return(sum((genotypes-admix.props%*%matrix(x,ncol=1))^2))
min.ind=function(x) return((1:length(x))[x==min(x)][1])
log.binom.calc=function(n,x,p) return(lchoose(n,x)+x*log(p)+(n-x)*log(1-p))
log.betabinom.calc=function(n,x,alpha,beta.val) return(lgamma(n+1)-lgamma(x+1)-lgamma(n-x+1)+lgamma(x+alpha)+lgamma(n-x+beta.val)-lgamma(n+alpha+beta.val)+lgamma(alpha+beta.val)-lgamma(alpha)-lgamma(beta.val))
log.normal.calc=function(x,mu,sd.val) return(-log(sqrt(2*pi))-log(sd.val)-0.5*((x-mu)^2/sd.val^2))

small.val.mle=0.001
post.admix.func=function(s,geno.vec,fixed.vec,lik.calc)
{
	prob.vec=(1+s)*fixed.vec/(1+s*fixed.vec)
	prob.vec[prob.vec<small.val.mle]=small.val.mle
	prob.vec[prob.vec>(1-small.val.mle)]=1-small.val.mle
	if (lik.calc!='normal') return(-1*sum(log.binom.calc(n=2,x=geno.vec,p=prob.vec)))
	if (lik.calc=='normal') return(-1*sum(log.normal.calc(x=geno.vec,mu=2*prob.vec,sd.val=sqrt(2*prob.vec*(1-prob.vec)))))
}

insurr.func=function(s,geno.vec,fixed.vec,anc.vec,freq.val,lik.calc)
{
	prob.vec=fixed.vec+anc.vec*((1+s)*freq.val/(1+s*freq.val))
	prob.vec[prob.vec<small.val.mle]=small.val.mle
	prob.vec[prob.vec>(1-small.val.mle)]=1-small.val.mle
	if (lik.calc!='normal') return(-1*sum(log.binom.calc(n=2,x=geno.vec,p=prob.vec)))
	if (lik.calc=='normal') return(-1*sum(log.normal.calc(x=geno.vec,mu=2*prob.vec,sd.val=sqrt(2*prob.vec*(1-prob.vec)))))
}

                                   ## (I) GET ADMIXTURE PROPS/INDS:
id.mat=read.table(id.file,as.is=TRUE)
admixture.inds=read.table(admixture.inds.file,as.is=TRUE)[,2]
admixture.props.all=read.table(admixture.props.file,as.is=T)
admixture.props.all=matrix(as.matrix(admixture.props.all),ncol=dim(admixture.props.all)[2])
rownames(admixture.props.all)=admixture.inds
colnames(admixture.props.all)=surrogate.vec
## grab admixed individuals from pop.vec
is.element.vec=(is.element(id.mat[,2],pop.vec) & as.integer(id.mat[,3])==1 & is.element(id.mat[,1],rownames(admixture.props.all)))
## grab donor individuals
is.element.donor.vec=(is.element(id.mat[,2],surrogate.vec) & as.integer(id.mat[,3])==1)

                                   ## (II) GET SOURCEFIND PROPS AND REMOVE INDS FROM TARGET (IF NECESSARY):
if (nonIberian.pops[1]!="")
{
	sourcefind.props.read=read.csv(sourcefind.props.file,sep=',',header=T,as.is=T)
	sourcefind.props.all=matrix(as.matrix(sourcefind.props.read[,18:dim(sourcefind.props.read)[2]]),ncol=dim(sourcefind.props.read)[2]-17)
	rownames(sourcefind.props.all)=as.character(sourcefind.props.read[,1])
	colnames(sourcefind.props.all)=names(sourcefind.props.read)[18:dim(sourcefind.props.read)[2]]
	target.index=id.mat[,1][is.element(id.mat[,2],pop.vec)]
	match.target=match(target.index,rownames(sourcefind.props.all))
	inds.to.remove=target.index[!is.na(match.target)][apply(sourcefind.props.all[match.target[!is.na(match.target)],match(nonIberian.pops,colnames(sourcefind.props.all))],1,sum)>=0.001]
	if (length(inds.to.remove)>0) is.element.vec[match(inds.to.remove,id.mat[,1])]=F
}

                         ## (III) GET SURROGATE ALLELE FREQS:
nsites.tot=0
for (i in 1:length(chromo.vec))
{
	readfile=gzfile(paste(geno.filePRE,chromo.vec[i],geno.filePOST,sep=''),open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
	pos.vec=as.character(line2[2:length(line2)])
	nsites.tot=nsites.tot+length(pos.vec)
	close(readfile)
}
rs.vec.final=pos.vec.final=chromo.vec.final=rep(NA,nsites.tot)
surrogate.freq.all=matrix(0,nrow=length(surrogate.vec),ncol=nsites.tot)
obs.count=matrix(0,nrow=length(pop.vec),ncol=nsites.tot)
snp.count=0
for (i in 1:length(chromo.vec))
{
	readfile=gzfile(paste(geno.filePRE,chromo.vec[i],geno.filePOST,sep=''),open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
	pos.vec=as.character(line2[2:length(line2)])
	start.snp=snp.count+1
	end.snp=snp.count+length(pos.vec)
	chromo.vec.final[start.snp:end.snp]=chromo.vec[i]
	###rs.vec.final[start.snp:end.snp]=rs.vec.geno[chromo.vec.geno==chromo.vec[i]][match(overlap.set,pos.vec.i)]
	pos.vec.final[start.snp:end.snp]=pos.vec
	print(c("SURR-FREQ",i,length(chromo.vec),length(pos.vec),nsites.tot,start.snp,end.snp))
	for (j in 1:dim(id.mat)[1])
	{
		hap1=as.integer(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
		hap2=as.integer(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
		if (is.element.vec[j])
		{
			obs.count[pop.vec==id.mat[j,2],start.snp:end.snp]=obs.count[pop.vec==id.mat[j,2],start.snp:end.snp]+(hap1+hap2)
			#obs.count.squared[pop.vec==id.mat[j,2],start.snp:end.snp]=obs.count.squared[pop.vec==id.mat[j,2],start.snp:end.snp]+(hap1+hap2)^2
		}
		if (is.element.donor.vec[j]) surrogate.freq.all[surrogate.vec==id.mat[j,2],start.snp:end.snp]=surrogate.freq.all[surrogate.vec==id.mat[j,2],start.snp:end.snp]+(hap1+hap2)
	}
	close(readfile)
	snp.count=snp.count+length(pos.vec)
}
for (k in 1:length(surrogate.vec)) surrogate.freq.all[k,]=surrogate.freq.all[k,]/(2*sum(id.mat[,2][is.element.donor.vec]==surrogate.vec[k]))
small.val=0.001
surrogate.freq.all[surrogate.freq.all<small.val]=small.val
surrogate.freq.all[surrogate.freq.all>(1-small.val)]=1-small.val

                         ## (IV) INFER DRIFT:
log.normal.drift.est.func=function(x,obs.vec,mu.vec,p.vec,alt.drift)
{
	x=exp(x)
	if (alt.drift==1) var.vec=x
	if (alt.drift==0) var.vec=x*p.vec*(1-p.vec)
	return(-sum(log.normal.calc(x=obs.vec,mu=mu.vec,sd.val=sqrt(var.vec))))
}
log.betabinom.drift.est.func=function(x,n.vec,x.vec,p.vec,alt.drift)
{
	drift.vector=exp(x)
	if (alt.drift==1)
	{
		drift.vector=drift.vector/(p.vec*(1-p.vec))
		drift.upper.val=0.99999
		drift.vector[drift.vector>drift.upper.val]=drift.upper.val
	}
	return(-sum(log.betabinom.calc(n=n.vec,x=x.vec,alpha=p.vec*(1-drift.vector)/drift.vector,beta.val=(1.0-p.vec)*(1-drift.vector)/drift.vector)))
}
drift.start=0.1
drift.est=rep(NA,length(pop.vec))
for (k in 1:length(pop.vec))
{
	admixture.props.k=admixture.props.all[match(id.mat[is.element.vec,1][id.mat[is.element.vec,2]==pop.vec[k]],rownames(admixture.props.all)),]
	freq.vec=c(matrix(apply(admixture.props.k,2,mean),nrow=1)%*%surrogate.freq.all)
	freq.vec[freq.vec<small.val]=small.val
	freq.vec[freq.vec>(1.0-small.val)]=1.0-small.val
	if (lik.calc=="normal") drift.est[k]=exp(optim(par=log(drift.start),log.normal.drift.est.func,obs.vec=obs.count[k,],mu.vec=2*c(matrix(apply(admixture.props.k,2,sum),nrow=1)%*%surrogate.freq.all),p.vec=freq.vec,alt.drift=alt.drift,method="Nelder-Mead")$par)
	if (lik.calc!="normal") drift.est[k]=exp(optim(par=log(drift.start),log.betabinom.drift.est.func,n.vec=2*sum(id.mat[,2][is.element.vec]==pop.vec[k]),x.vec=obs.count[k,],p.vec=freq.vec,alt.drift=alt.drift,method="Nelder-Mead")$par)
}
drift.est

                         ## (V) TEST EACH SNP FOR SELECTION (I.E. DEVIATIONS FROM NEUTRALITY):
pval.mat=matrix(NA,nrow=length(pop.vec),ncol=nsites.tot)
for (k in 1:length(pop.vec))
{
	admixture.props.k=admixture.props.all[match(id.mat[is.element.vec,1][id.mat[is.element.vec,2]==pop.vec[k]],rownames(admixture.props.all)),]
	neutral.freq.vec=c(matrix(apply(admixture.props.k,2,mean),nrow=1)%*%surrogate.freq.all)
	neutral.freq.vec[neutral.freq.vec<small.val]=small.val
	neutral.freq.vec[neutral.freq.vec>(1.0-small.val)]=1.0-small.val
	neutral.counts.vec=2*matrix(apply(admixture.props.k,2,sum),nrow=1)%*%surrogate.freq.all
	pval.vec=rep(NA,dim(obs.count)[2])
	if (lik.calc=="normal")
	{
		if (alt.drift==1) var.vec=rep(drift.est[k],nsites.tot)
		if (alt.drift==0) var.vec=neutral.freq.vec*(1-neutral.freq.vec)*drift.est[k]
		pval.vec[obs.count[k,]<=neutral.counts.vec]=2*pnorm(obs.count[k,][obs.count[k,]<=neutral.counts.vec],mean=neutral.counts.vec[obs.count[k,]<=neutral.counts.vec],sd=sqrt(var.vec)[obs.count[k,]<=neutral.counts.vec])
		pval.vec[obs.count[k,]>neutral.counts.vec]=2*pnorm(obs.count[k,][obs.count[k,]>neutral.counts.vec],mean=neutral.counts.vec[obs.count[k,]>neutral.counts.vec],sd=sqrt(var.vec)[obs.count[k,]>neutral.counts.vec],lower.tail=FALSE)
	}
	if (lik.calc!="normal")
	{
		drift.vector=rep(drift.est[k],length(neutral.freq.vec))
		if (alt.drift==1)
		{
			drift.vector=drift.vector/(neutral.freq.vec*(1-neutral.freq.vec))
			drift.upper.val=0.99999
			drift.vector[drift.vector>drift.upper.val]=drift.upper.val
		}
		for (j in 1:dim(obs.count)[2])
		{
			if (obs.count[k,j] <= neutral.counts.vec[j]) range.tocalc=1:(obs.count[k,j]+1)
			if (obs.count[k,j] > neutral.counts.vec[j]) range.tocalc=(obs.count[k,j]+1):(2*sum(id.mat[,2][is.element.vec]==pop.vec[k])+1)
			prob.vec.all=log.betabinom.calc(n=2*sum(id.mat[,2][is.element.vec]==pop.vec[k]),x=0:(2*sum(id.mat[,2][is.element.vec]==pop.vec[k])),alpha=neutral.freq.vec[j]*(1-drift.vector[j])/drift.vector[j],beta.val=(1.0-neutral.freq.vec[j])*(1-drift.vector[j])/drift.vector[j])
			if (round(abs(sum(exp(prob.vec.all))-1.0),5)!=0)
			{
				print(c(j,sum(exp(prob.vec.all)),neutral.freq.vec[j],(1-drift.vector[j])/drift.vector[j]))
				#break
			}
			pval.vec[j]=2*sum(exp(prob.vec.all[range.tocalc]))
		}
	}
	pval.vec[pval.vec>1]=1
	pval.mat[k,]=-log10(pval.vec)
}

                         ## (VI) FIND AIC PER POP AND PRINT:
AIC.neutral=AIC.postadmix=LRT.postadmix=selfac.postadmix=exp.freq.pop.k=matrix(NA,nrow=length(pop.vec),ncol=nsites.tot)
AIC.insurr=LRT.insurr=selfac.insurr=array(NA,dim=c(length(pop.vec),length(surrogate.vec),nsites.tot))
snp.count=0
for (i in 1:length(chromo.vec))
{
	readfile=gzfile(paste(geno.filePRE,chromo.vec[i],geno.filePOST,sep=''),open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
	pos.vec=as.character(line2[2:length(line2)])
	start.snp=snp.count+1
	end.snp=snp.count+length(pos.vec)
	snp.count=snp.count+length(pos.vec)
	print(c("SCORING",i,length(chromo.vec),length(pos.vec),nsites.tot,start.snp,end.snp))
	close(readfile)
	nsites.i=length(pos.vec)
	num.runs=ceiling(nsites.i/nsnps.perrun)
	for (a in 1:num.runs)
	{
		start.snp.a=(a-1)*nsnps.perrun+1
		end.snp.a=a*nsnps.perrun
		if (end.snp.a > nsites.i) end.snp.a=nsites.i
		print(c(a,num.runs,start.snp.a,end.snp.a))
		readfile=gzfile(paste(geno.filePRE,chromo.vec[i],geno.filePOST,sep=''),open='r')
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
		geno.data.a=matrix(NA,nrow=sum(is.element.vec),ncol=end.snp.a-start.snp.a+1)
		target.count=1
		for (j in 1:dim(id.mat)[1])
		{
			hap1=as.integer(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
			hap2=as.integer(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
			if (is.element.vec[j])
			{
				geno.freq.j=(hap1+hap2)[start.snp.a:end.snp.a]
				geno.data.a[target.count,]=geno.freq.j
				target.count=target.count+1
			}
		}
		close(readfile)
		surrogate.freq=matrix(surrogate.freq.all[,start.snp:end.snp][,start.snp.a:end.snp.a],nrow=dim(surrogate.freq.all)[1])

                                           ## FIND AIC-ish SCORE:
		for (k in 1:length(pop.vec))
		{
			admixture.props.k=admixture.props.all[match(id.mat[is.element.vec,1][id.mat[is.element.vec,2]==pop.vec[k]],rownames(admixture.props.all)),]
			inds.k=(1:dim(geno.data.a)[1])[id.mat[is.element.vec,2]==pop.vec[k]]
			neutral.freqs.h=admixture.props.k%*%surrogate.freq
			neutral.freqs.h[neutral.freqs.h<small.val]=small.val
			neutral.freqs.h[neutral.freqs.h>(1-small.val)]=1-small.val
			exp.freq.pop.k[k,(start.snp:end.snp)[start.snp.a:end.snp.a]]=apply(neutral.freqs.h,2,mean)
			for (j in 1:(end.snp.a-start.snp.a+1))
			{
				if (lik.calc!='normal') prob.neutral=sum(log.binom.calc(n=2,x=geno.data.a[inds.k,j],p=neutral.freqs.h[,j]))
				if (lik.calc=='normal') prob.neutral=sum(log.normal.calc(x=geno.data.a[inds.k,j],mu=2*neutral.freqs.h[,j],sd.val=sqrt(2*neutral.freqs.h[,j]*(1-neutral.freqs.h[,j]))))
				mle.postadmix=optim(0,post.admix.func,method="Nelder-Mead",geno.vec=geno.data.a[inds.k,j],fixed.vec=neutral.freqs.h[,j],lik.calc=lik.calc)
				selfac.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=mle.postadmix$par
				prob.postadmix=-1*mle.postadmix$value
				AIC.neutral[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=-2*prob.neutral
				AIC.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=2-2*prob.postadmix
				LRT.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=pchisq(2*(prob.postadmix-prob.neutral),1,log.p=TRUE,lower.tail=FALSE)
				for (s in 1:dim(admixture.props.k)[2])
				{
					fixed.mat=c(matrix(admixture.props.k[,-s],ncol=dim(admixture.props.k)[2]-1)%*%matrix(surrogate.freq[-s,j],ncol=1))
					mle.insurr=optim(0,insurr.func,method="Nelder-Mead",geno.vec=geno.data.a[inds.k,j],fixed.vec=fixed.mat,anc.vec=admixture.props.k[,s],freq.val=surrogate.freq[s,j],lik.calc=lik.calc)
					prob.insurr=-1*mle.insurr$value
					selfac.insurr[k,s,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=mle.insurr$par
					AIC.insurr[k,s,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=2-2*prob.insurr
					LRT.insurr[k,s,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=pchisq(2*(prob.insurr-prob.neutral),1,log.p=TRUE,lower.tail=FALSE)
				}
				##AIC.diff.postadmix=exp((apply(cbind(AIC.postadmix[k,],AIC.neutral[k,]),1,min)-apply(cbind(AIC.postadmix[k,],AIC.neutral[k,]),1,max))/2)   ## models are AIC.diff times as probable as the best-fitting model to minimize the information loss by using this model (relative to the truth) -- i.e. low value loses, "1" is the best-fitting model
				##AIC.diff.postadmix.winner=apply(cbind(AIC.postadmix[k,],AIC.neutral[k,]),1,min.ind)
			}
		}
	}
}

                         ## (VII) PRINT OUT:
to.print.vec=NULL
for (k in 1:length(pop.vec)) to.print.vec=c(to.print.vec,paste("log10.pval.target.",k,sep=''),paste("obs.freq.target.",k,sep=''),paste("exp.freq.target.",k,sep=''),paste("AIC.neutral.target.",k,sep=''),paste("AIC.postadmix.target.",k,sep=''),paste("AIC.insurr.source",1:length(surrogate.vec),"target.",k,sep=''),paste("sel.postadmix.target.",k,sep=''),paste("sel.insurr.source",1:length(surrogate.vec),"target.",k,sep=''),paste("LRT.pval.postadmix.target.",k,sep=''),paste("LRT.pval.insurr.source",1:length(surrogate.vec),".target.",k,sep=''))
write(paste("drift.",pop.vec,sep=''),file=out.file,ncolumns=length(drift.est))
write(drift.est,file=out.file,ncolumns=length(drift.est),append=TRUE)
write(c("chrom","pos",to.print.vec),file=out.file,ncolumns=length(to.print.vec)+2,append=TRUE)
to.print.mat=NULL
for (k in 1:length(pop.vec)) to.print.mat=cbind(to.print.mat,round(pval.mat[k,],round.val),round(obs.count[k,]/(2*sum(id.mat[is.element.vec,2]==pop.vec[k])),round.val.freq),round(exp.freq.pop.k[k,],round.val.freq),round(AIC.neutral[k,],round.val),round(AIC.postadmix[k,],round.val),round(t(AIC.insurr[k,,]),round.val),round(selfac.postadmix[k,],round.val.sel),round(t(selfac.insurr[k,,]),round.val.sel),round(LRT.postadmix[k,],round.val),round(t(LRT.insurr[k,,]),round.val))
write.table(cbind(chromo.vec.final,pos.vec.final,to.print.mat),file=out.file,row.names=F,col.names=F,quote=F,append=TRUE)	

warnings()

drift.est
table(id.mat[is.element.vec,2])
table(id.mat[is.element.donor.vec,2])

q(save='no')
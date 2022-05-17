## NOTE: UNDER THE NORMAL-LIK MODEL, YOU CAN RESCALE THE LRT BY ADDING MORE DRIFT; IN PARTICULAR LRT(new) = [sqrt(drift.vec.old)/sqrt(drift.vec.new)]*LRT(old) 
## THUS HALVING THE LRT IS ACCOMPLISHED BY DOUBLING drift.vec.old

## Normal lik doesn't work particularly well with alt.drift=0 -- really favors SNPs near fixation

## usage:   R < run_AdaptMix.R parameter.input.file genotypes.input.filenames id.file output.file --no-save > screenoutput.out

############################
## INPUT:

usage=function()
{
	print(noquote("run using: R < run_AdaptMix.R parameter.input.file genotypes.input.filenames id.file output.file --no-save > screenoutput.out"))
	print(noquote("parameter.input.file format (NO defaults, all fields must be entered, even if not used):"))
	print(noquote("pop.vec: [pop1 pop2 ... pop_K]"))
	print(noquote("surrogate.vec: [pop_1 pop_2 ... pop_S]"))
	print(noquote("drift.maf.bins: [0.0,...,0.5]"))
	print(noquote("min.allele.freq.shift: [0.0,...,1.0]"))
}

temp=commandArgs()

param.infile=as.character(temp[2])
if (param.infile=="help"){usage();q(save='no')}

error.input.message=function(file.name)
  {
    print(paste("Something wrong with input file ",file.name,". See below. Exiting...",sep=''))
    usage()
    q(save='no')
  }
line.check=function(file.name,skip.val,match.val)
  {
    if (as.character(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[1]) != as.character(match.val))
      {
        error.input.message(file.name)
      }
    if (as.character(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[1]) == as.character(match.val))
      return(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE)[2:length(read.table(file.name,skip=skip.val,nrows=1,as.is=TRUE))])
  }

          ## line read in and checks:
pop.vec=unique(as.character(line.check(param.infile,0,"pop.vec:")))
if (length(pop.vec)<1) error.input.message(param.infile)
surrogate.vec=as.character(line.check(param.infile,1,"surrogate.vec:"))
if (length(unique(surrogate.vec))<2) error.input.message(param.infile)
drift.maf.bins=seq(0,0.5,0.05)
drift.maf.bins=as.double(line.check(param.infile,2,"drift.maf.bins:"))  ## calculating drift separately for SNPs with minor-allele expected frequencies within these bins
if ((length(drift.maf.bins)<2)||(drift.maf.bins[1]!=0) || (drift.maf.bins[length(drift.maf.bins)]!=0.5)){print(paste("Something wrong with drift.maf.bins in ",param.infile,". Values must start at 0, end at 0.5, and be sequential. Exiting....",sep=''));q(save='no')}
if (sum(drift.maf.bins[2:length(drift.maf.bins)]-drift.maf.bins[1:(length(drift.maf.bins)-1)]<0)>0){print(paste("Something wrong with drift.maf.bins in ",param.infile,". Values must start at 0, end at 0.5, and be sequential. Exiting....",sep=''));q(save='no')}
min.allele.freq.shift=as.double(line.check(param.infile,3,"min.allele.freq.shift:"))  ## minimum drift -- make sure alleles shift by at least this much to be interesting
if ((length(min.allele.freq.shift)!=1)||(min.allele.freq.shift<=0) || (min.allele.freq.shift>=1)){print(paste("Something wrong with min.allele.freq.shift in ",param.infile,". Value must be between 0 and 1. (Make very small, but above 0, to be essentially ignored.) Exiting....",sep=''));q(save='no')}

genotypes.filein=as.character(temp[3])
id.file=as.character(temp[4])
out.file <- as.character(temp[5])

                      ## PARAMETERS:
nsnps.perrun=100000         ## TO CONTROL RAM
lik.calc='betabinom'           ## 'betabinom' or 'normal' likelihoods?
alt.drift=1
round.val=3     ## for AIC
round.val.sel=3     ## for sel-coeff
round.val.freq=4     ## for allele frequencies (observed and expected)

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
if (dim(id.mat)[2]!=length(surrogate.vec)+2){print(paste("Something wrong with file ",id.file," -- should have columns equal to number of surrogate populations + 2. Exiting...",sep=''));q(save='no')}
admixture.props.all=matrix(as.matrix(id.mat[,3:dim(id.mat)[2]]),ncol=dim(id.mat)[2]-2)
rownames(admixture.props.all)=id.mat[,1]
colnames(admixture.props.all)=surrogate.vec
           ## grab admixed individuals from pop.vec:
for (i in 1:length(pop.vec)) if (sum(is.element(id.mat[,2],pop.vec[i]))==0) {print(paste("No individuals found in tested pop ",pop.vec[i],"! Exiting....",sep=''));q(save='no')}
is.element.vec=is.element(id.mat[,2],pop.vec)
           ## grab donor individuals:
for (i in 1:length(surrogate.vec)) if (sum(is.element(id.mat[,2],surrogate.vec[i]))==0) {print(paste("No individuals found in surrogate pop ",surrogate.vec[i],"! Exiting....",sep=''));q(save='no')}
is.element.donor.vec=is.element(id.mat[,2],surrogate.vec)

                                   ## (II) FIND NUMBER OF GENOTYPE FILES:
num.genofiles=0
readfile=file(genotypes.filein,open="r")
line2=1
while(!is.na(line2[1]))
{
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)
    	if (!is.na(line2[1])) num.genofiles=num.genofiles+1
}
if (num.genofiles==0)
{
	print(paste("SOMETHING WRONG WITH INPUT FILE ",genotypes.filein," (PERHAPS EMPTY?). Exiting....",sep=''))
	q(save='no')
}
close(readfile)

                         ## (III) GET TARGET AND SURROGATE ALLELE FREQS:
nsites.tot=0
readfileNAMES=file(genotypes.filein,open="r")
for (i in 1:num.genofiles)
{
	filename.i=as.character(scan(readfileNAMES,nlines=1,what='char',quiet=TRUE))
	readfile=gzfile(filename.i,open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
	pos.vec=as.character(line2[2:length(line2)])
	nsites.tot=nsites.tot+length(pos.vec)
	close(readfile)
}
close(readfileNAMES)
rs.vec.final=pos.vec.final=chromo.vec.final=rep(NA,nsites.tot)
surrogate.freq.all=surrogate.count.all=matrix(0,nrow=length(surrogate.vec),ncol=nsites.tot)
obs.count=n.count=matrix(0,nrow=length(pop.vec),ncol=nsites.tot)
snp.count=0
readfileNAMES=file(genotypes.filein,open="r")
for (i in 1:num.genofiles)
{
	print(paste("Reading in data for file ",i," of ",num.genofiles," to calculate drift and perform selection scan....",sep=''))

	filename.i=as.character(scan(readfileNAMES,nlines=1,what='char',quiet=TRUE))
	readfile=gzfile(filename.i,open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
	pos.vec=as.character(line2[2:length(line2)])
	start.snp=snp.count+1
	end.snp=snp.count+length(pos.vec)
	chromo.vec.final[start.snp:end.snp]=i
	pos.vec.final[start.snp:end.snp]=pos.vec
	for (j in 1:dim(id.mat)[1])
	{
		hap1=as.character(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
		hap2=as.character(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
		if (is.element.vec[j] || is.element.donor.vec[j])
		{
			if (sum(is.element(hap1,c("0","1","?")))!=length(hap1)){print(paste("Haplotype in file ",filename.i," contains value other than [0,1,?]. Exiting....",sep=''));q(save='no')}
			if (sum(is.element(hap2,c("0","1","?")))!=length(hap2)){print(paste("Haplotype in file ",filename.i," contains value other than [0,1,?]. Exiting....",sep=''));q(save='no')}
			hap1[hap1=="?"]=-9
			hap2[hap2=="?"]=-9
			hap1=as.integer(hap1)
			hap2=as.integer(hap2)
			tot.counts=hap1+hap2
		}
		if (is.element.vec[j])
		{
			obs.count[pop.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]=obs.count[pop.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]+tot.counts[tot.counts>=0]
			n.count[pop.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]=n.count[pop.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]+2
			#obs.count.squared[pop.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]=obs.count.squared[pop.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]+tot.counts[tot.counts>=0]^2
		}
		if (is.element.donor.vec[j])
		{
			surrogate.freq.all[surrogate.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]=surrogate.freq.all[surrogate.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]+tot.counts[tot.counts>=0]
			surrogate.count.all[surrogate.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]=surrogate.count.all[surrogate.vec==id.mat[j,2],(start.snp:end.snp)[tot.counts>=0]]+2
		}
	}
	close(readfile)
	snp.count=snp.count+length(pos.vec)
}
close(readfileNAMES)
for (k in 1:length(surrogate.vec)) surrogate.freq.all[k,]=surrogate.freq.all[k,]/surrogate.count.all[k,]
small.val=0.001
surrogate.freq.all[surrogate.freq.all<small.val]=small.val
surrogate.freq.all[surrogate.freq.all>(1-small.val)]=1-small.val
print("Finished reading in data for drift calculation and selection scan.")

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

##drift.start=0.1
drift.maf.bins[length(drift.maf.bins)]=0.51
drift.est=matrix(NA,nrow=length(pop.vec),ncol=length(drift.maf.bins)-1)
for (k in 1:length(pop.vec))
{
	admixture.props.k=admixture.props.all[match(id.mat[is.element.vec,1][id.mat[is.element.vec,2]==pop.vec[k]],rownames(admixture.props.all)),]
	freq.vec=c(matrix(apply(admixture.props.k,2,mean),nrow=1)%*%surrogate.freq.all)
	freq.vec[freq.vec<small.val]=small.val
	freq.vec[freq.vec>(1.0-small.val)]=1.0-small.val
	freq.vec.maf=freq.vec
	freq.vec.maf[freq.vec.maf>0.5]=1-freq.vec.maf[freq.vec.maf>0.5]
	for (j in 1:(length(drift.maf.bins)-1))
	{
		drift.snps.toselect=(1:length(freq.vec))[freq.vec.maf>=drift.maf.bins[j] & freq.vec.maf<drift.maf.bins[(j+1)]]
		if (length(drift.snps.toselect)==0){ print(paste("Target pop ",pop.vec[k]," does not have any expected allele frequencies between ",drift.maf.bins[j]," and ",drift.maf.bins[(j+1)],"! Should change MAF bin size for calculating drift. Exiting...",sep='')); q(save='no')}
		drift.start=min(c(min(freq.vec.maf[drift.snps.toselect][freq.vec.maf[drift.snps.toselect]>0]),0.01))
		drift.start=drift.start-drift.start/1000000
		if (length(drift.snps.toselect)>0)
		{
			if (lik.calc=="normal") drift.est[k,j]=exp(optim(par=log(drift.start),log.normal.drift.est.func,obs.vec=obs.count[k,drift.snps.toselect],mu.vec=2*c(matrix(apply(admixture.props.k,2,sum),nrow=1)%*%surrogate.freq.all[,drift.snps.toselect]),p.vec=freq.vec[drift.snps.toselect],alt.drift=alt.drift,method="Nelder-Mead")$par)
			if (lik.calc!="normal") drift.est[k,j]=exp(optim(par=log(drift.start),log.betabinom.drift.est.func,n.vec=n.count[k,drift.snps.toselect],x.vec=obs.count[k,drift.snps.toselect],p.vec=freq.vec[drift.snps.toselect],alt.drift=alt.drift,method="Nelder-Mead")$par)
		}
		if (drift.est[k,j]<(min.allele.freq.shift^2)) drift.est[k,j]=min.allele.freq.shift^2
	}
}
print("Finished calculating drift.")

                         ## (V) TEST EACH SNP FOR SELECTION (I.E. DEVIATIONS FROM NEUTRALITY):
pval.mat=matrix(NA,nrow=length(pop.vec),ncol=nsites.tot)
for (k in 1:length(pop.vec))
{
	admixture.props.k=admixture.props.all[match(id.mat[is.element.vec,1][id.mat[is.element.vec,2]==pop.vec[k]],rownames(admixture.props.all)),]
	neutral.freq.vec=c(matrix(apply(admixture.props.k,2,mean),nrow=1)%*%surrogate.freq.all)
	neutral.freq.vec[neutral.freq.vec<small.val]=small.val
	neutral.freq.vec[neutral.freq.vec>(1.0-small.val)]=1.0-small.val
	neutral.counts.vec=2*matrix(apply(admixture.props.k,2,sum),nrow=1)%*%surrogate.freq.all
	freq.vec.maf=neutral.freq.vec
	freq.vec.maf[freq.vec.maf>0.5]=1-freq.vec.maf[freq.vec.maf>0.5]
	pval.vec=rep(NA,dim(obs.count)[2])
	if (lik.calc=="normal")
	{
		var.vec=rep(NA,nsites.tot)
		if (alt.drift==1)
		{
			for (j in 1:(length(drift.maf.bins)-1)) var.vec[freq.vec.maf>=drift.maf.bins[j] & freq.vec.maf<drift.maf.bins[(j+1)]]=drift.est[k,j]
		}
		if (alt.drift==0)
		{
			for (j in 1:(length(drift.maf.bins)-1)) var.vec[freq.vec.maf>=drift.maf.bins[j] & freq.vec.maf<drift.maf.bins[(j+1)]]=neutral.freq.vec[freq.vec.maf>=drift.maf.bins[j] & freq.vec.maf<drift.maf.bins[(j+1)]]*(1-neutral.freq.vec[freq.vec.maf>=drift.maf.bins[j] & freq.vec.maf<drift.maf.bins[(j+1)]])*drift.est[k,j]
		}
		pval.vec[obs.count[k,]<=neutral.counts.vec]=2*pnorm(obs.count[k,][obs.count[k,]<=neutral.counts.vec],mean=neutral.counts.vec[obs.count[k,]<=neutral.counts.vec],sd=sqrt(var.vec)[obs.count[k,]<=neutral.counts.vec])
		pval.vec[obs.count[k,]>neutral.counts.vec]=2*pnorm(obs.count[k,][obs.count[k,]>neutral.counts.vec],mean=neutral.counts.vec[obs.count[k,]>neutral.counts.vec],sd=sqrt(var.vec)[obs.count[k,]>neutral.counts.vec],lower.tail=FALSE)
	}
	if (lik.calc!="normal")
	{
		drift.vector=rep(NA,length(neutral.freq.vec))
		for (j in 1:(length(drift.maf.bins)-1)) drift.vector[freq.vec.maf>=drift.maf.bins[j] & freq.vec.maf<drift.maf.bins[(j+1)]]=drift.est[k,j]
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
print("Finished selection test.")

                         ## (VI) FIND AIC PER POP AND PRINT:
AIC.neutral=AIC.postadmix=LRT.postadmix=selfac.postadmix=exp.freq.pop.k=matrix(NA,nrow=length(pop.vec),ncol=nsites.tot)
AIC.insurr=LRT.insurr=selfac.insurr=array(NA,dim=c(length(pop.vec),length(surrogate.vec),nsites.tot))
snp.count=0
readfileNAMES=file(genotypes.filein,open="r")
for (i in 1:num.genofiles)
{
	print(paste("Testing for post vs pre-admixture selection in file ",i," of ",num.genofiles,"....",sep=''))

	filename.i=as.character(scan(readfileNAMES,nlines=1,what='char',quiet=TRUE))
	readfile=gzfile(filename.i,open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
	pos.vec=as.character(line2[2:length(line2)])
	start.snp=snp.count+1
	end.snp=snp.count+length(pos.vec)
	snp.count=snp.count+length(pos.vec)
	close(readfile)
	nsites.i=length(pos.vec)
	num.runs=ceiling(nsites.i/nsnps.perrun)
	for (a in 1:num.runs)
	{
		start.snp.a=(a-1)*nsnps.perrun+1
		end.snp.a=a*nsnps.perrun
		if (end.snp.a > nsites.i) end.snp.a=nsites.i
		# print(c(a,num.runs,start.snp.a,end.snp.a))
		readfile=gzfile(filename.i,open='r')
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
		geno.data.a=matrix(NA,nrow=sum(is.element.vec),ncol=end.snp.a-start.snp.a+1)
		target.count=1
		for (j in 1:dim(id.mat)[1])
		{
			hap1=as.character(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
			hap2=as.character(strsplit(scan(readfile,nlines=1,what='char',quiet=TRUE),split='')[[1]])
			if (is.element.vec[j])
			{
				hap1[hap1=="?"]="-9"
				hap2[hap2=="?"]="-9"
				hap1=as.integer(hap1)
				hap2=as.integer(hap2)
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
				if (sum(geno.data.a[inds.k,j]>=0)==0)
				{
					selfac.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=NA
					AIC.neutral[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=NA
					AIC.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=NA
					LRT.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=NA
					for (s in 1:dim(admixture.props.k)[2])
					{
						selfac.insurr[k,s,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=NA
						AIC.insurr[k,s,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=NA
						LRT.insurr[k,s,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=NA
					}
				}
				if (sum(geno.data.a[inds.k,j]>=0)>0)
				{
					if (lik.calc!='normal') prob.neutral=sum(log.binom.calc(n=2,x=geno.data.a[inds.k,j][geno.data.a[inds.k,j]>=0],p=neutral.freqs.h[,j][geno.data.a[inds.k,j]>=0]))
					if (lik.calc=='normal') prob.neutral=sum(log.normal.calc(x=geno.data.a[inds.k,j][geno.data.a[inds.k,j]>=0],mu=2*neutral.freqs.h[,j][geno.data.a[inds.k,j]>=0],sd.val=sqrt(2*neutral.freqs.h[,j][geno.data.a[inds.k,j]>=0]*(1-neutral.freqs.h[,j][geno.data.a[inds.k,j]>=0]))))
					mle.postadmix=optim(0,post.admix.func,method="Nelder-Mead",geno.vec=geno.data.a[inds.k,j][geno.data.a[inds.k,j]>=0],fixed.vec=neutral.freqs.h[,j][geno.data.a[inds.k,j]>=0],lik.calc=lik.calc)
					selfac.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=mle.postadmix$par
					prob.postadmix=-1*mle.postadmix$value
					AIC.neutral[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=-2*prob.neutral
					AIC.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=2-2*prob.postadmix
					LRT.postadmix[k,(start.snp:end.snp)[start.snp.a:end.snp.a][j]]=pchisq(2*(prob.postadmix-prob.neutral),1,log.p=TRUE,lower.tail=FALSE)
					for (s in 1:dim(admixture.props.k)[2])
					{
						fixed.mat=c(matrix(admixture.props.k[,-s],ncol=dim(admixture.props.k)[2]-1)%*%matrix(surrogate.freq[-s,j],ncol=1))
						mle.insurr=optim(0,insurr.func,method="Nelder-Mead",geno.vec=geno.data.a[inds.k,j][geno.data.a[inds.k,j]>=0],fixed.vec=fixed.mat[geno.data.a[inds.k,j]>=0],anc.vec=admixture.props.k[,s][geno.data.a[inds.k,j]>=0],freq.val=surrogate.freq[s,j],lik.calc=lik.calc)
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
}
close(readfileNAMES)

                         ## (VII) PRINT OUT:
to.print.vec=NULL
for (k in 1:length(pop.vec)) to.print.vec=c(to.print.vec,paste("log10.pval.target.",k,sep=''),paste("obs.freq.target.",k,sep=''),paste("exp.freq.target.",k,sep=''),paste("AIC.neutral.target.",k,sep=''),paste("AIC.postadmix.target.",k,sep=''),paste("AIC.insurr.source",1:length(surrogate.vec),"target.",k,sep=''),paste("sel.postadmix.target.",k,sep=''),paste("sel.insurr.source",1:length(surrogate.vec),"target.",k,sep=''))
drift.maf.bins.label=paste("[",drift.maf.bins[1:(length(drift.maf.bins)-1)],",",drift.maf.bins[2:length(drift.maf.bins)],")",sep='')
write.table(rbind(c("drift.est",drift.maf.bins.label),cbind(pop.vec,round(drift.est,6))),file=out.file,row.names=FALSE,col.names=FALSE,quote=FALSE)
write(c("file","pos",to.print.vec),file=out.file,ncolumns=length(to.print.vec)+2,append=TRUE)
to.print.mat=NULL
for (k in 1:length(pop.vec)) to.print.mat=cbind(to.print.mat,round(pval.mat[k,],round.val),round(obs.count[k,]/n.count[k,],round.val.freq),round(exp.freq.pop.k[k,],round.val.freq),round(AIC.neutral[k,],round.val),round(AIC.postadmix[k,],round.val),round(t(AIC.insurr[k,,]),round.val),round(selfac.postadmix[k,],round.val.sel),round(t(selfac.insurr[k,,]),round.val.sel))
write.table(cbind(chromo.vec.final,pos.vec.final,to.print.mat),file=out.file,row.names=F,col.names=F,quote=F,append=TRUE)	

print("Finished!")



q(save='no')
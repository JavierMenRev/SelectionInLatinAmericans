## simulates A SINGLE SNP under selection
## output file in CP format, puts selected SNP at far right

## POTENTIALLY MAJOR ISSUE WITH DECIDING WHICH SOURCE "DRIFT" VALUES MATCH REAL DATA:
## Currently I am selecting "drift" by comparing the surrogate allele frequencies with:
## (1) "source" frequencies inferred using target inds (cor.truth)
## (2) "source" frequencies simulated based on your drift values (cor.sim)
## and trying to toggle drift until the comparison of (1) matches that of (2)
## But, because I sample target inds' frequencies using a binomial distribution that uses the simulated "source" frequencies in (2), we should instead of (2) be calculating
## (3) "source" frequencies inferred using simulated target inds (cor.sim)
## and try to toggle drift until the comparison of (1) matches that of (3)



## SIMULATION PROTOCOL:
##      selection pre-admixture:
## (1) sample sim source pops' allele frequencies via beta according to surrogate allele frequencies plus user-defined drift
## (2) for one pop in (1), for selected SNPs forward simulate randomly-mating pop of size N diploids for G generations with selection (incorporating dominance, etc) 
## (3) sample admixed inds' data, at selected SNPs + X (user-specified) neutral ones, via binomial according to ADMIXTURE props * sim source pops' allele frequencies
## (4) use real-data surrogate pop allele frequencies

##      selection post-admixture:
##                (Currently this simulates an equal number of generations of selection in each source pop and then combines them to make admixed pop. NOT IDEAL, but aim is to add selection while keeping admixture props of inds the same over G generations; hence can't just make big combined pop -- would have to simulate one pop per ind and add selection)
## (1) sample sim source pops' allele frequencies via beta according to surrogate allele frequencies plus drift
## (2) at X neutral snps, sample admixed inds' data via binomial according to ADMIXTURE props * sim source pops' allele frequencies
## (3) at selected snps:
##       (3a) forward simulate K (sources) randomly-mating pops of size N_k diploids for G generations with selection (incorporating dominance, etc) 
##       (3b) sample admixed inds' data via binomial according to ADMIXTURE props * sim source pops' allele frequencies (with no additional added drift)
## (4) use real-data surrogate pop allele frequencies

## usage:   R < AdaptMixSimulator.R parameter.input.file genotypes.input.filenames id.file output.file --no-save > screenoutput.out

## example: R < AdaptMixSimulator.R simexample/PEL_simulation_paramfile.txt simexample/PEL_REFs_ALLCHR_chr.txt simexample/PEL_REFs.ids.txt simexample/PEL_SIM_ALLCHR --no-save > screenoutput.out

############################
## INPUT:

usage=function()
{
	print(noquote("run using: R < AdaptMixSimulator.R parameter.input.file genotypes.input.filenames id.file output.file --no-save > screenoutput.out"))
	print(noquote("parameter.input.file format (NO defaults, all fields must be entered, in this order, even if not used):"))
	print(noquote("selection.post-admixture?: [0,1]"))
	print(noquote("sel.coeff: [0.0,....]"))
	print(noquote("sel.type: [additive,dominant,multiplicative,recessive]"))
	print(noquote("target.pop: [name]"))
	print(noquote("surrogate.pops: [name1,name2,...]"))
	print(noquote("sources.with.selection.preadmixture: [0,1]"))
	print(noquote("generations.selection.each.source: [1,...]"))
	print(noquote("pop.size.sources: [1,...]"))
	print(noquote("drift.btwn.surrogates.and.sources: [0.0,...,1.0]"))
	print(noquote("num.neutral.snps: [1,...]"))
	print(noquote("range.startfrequency.selected.snp: [0.0,...] [...,1.0]"))
	print(noquote("infer.source.freq.using.target.data: [0,1]"))
	print(noquote("divide.into.runs.ofXX.inds(to.reduce.RAM): [1,...]"))
}

                  ## ADMIXTURE/SAMPLE INFO:
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
sel.postadmix=unique(as.integer(line.check(param.infile,0,"selection.post-admixture?:")))
if (length(sel.postadmix)!=1) error.input.message(param.infile)
sel.prob.val=as.double(line.check(param.infile,1,"sel.coeff:")) 
if ((length(sel.prob.val)!=1)||is.na(sel.prob.val)){print(paste("Something wrong with selection coefficient in ",param.infile,". Exiting....",sep=''));q(save='no')}
sel.type=as.character(line.check(param.infile,2,"sel.type:"))
if (length(sel.type)!=1) error.input.message(param.infile)
target.pop=as.character(line.check(param.infile,3,"target.pop:"))
if (length(target.pop)!=1) error.input.message(param.infile)
source.pops=as.character(line.check(param.infile,4,"surrogate.pops:"))
if (length(source.pops)<1){print(paste("Something wrong with surrogate pops specified in ",param.infile,". Exiting....",sep=''));q(save='no')}
sel.simpops=as.integer(line.check(param.infile,5,"sources.with.selection.preadmixture:"))     ## simulate selection in source pops pre-admixture?
if ((length(sel.simpops)!=length(source.pops))||sum(sel.simpops==0 | sel.simpops==1,na.rm=T)!=length(sel.simpops)){print(paste("Something wrong with sources chosen to be selected as specified in ",param.infile,". Exiting....",sep=''));q(save='no')}
gen.sel.source=as.integer(line.check(param.infile,6,"generations.selection.each.source:"))  ## generations of selection in each source (pre-admixture if "sel.postadmix=0"; otherwise (meant to reflect) post-admixture)
if ((length(gen.sel.source)!=length(sel.simpops))||sum(gen.sel.source>=0,na.rm=T)!=length(gen.sel.source)){print(paste("Something wrong with number of generations of selection specified for each source specified in ",param.infile,". Exiting....",sep=''));q(save='no')}
pop.size.vec=as.double(line.check(param.infile,7,"pop.size.sources:"))    ## diploid pop size of each source pop (used for simulating selection pre or post-admixture)
if ((length(pop.size.vec)!=length(sel.simpops))||sum(pop.size.vec>0,na.rm=T)!=length(pop.size.vec)){print(paste("Something wrong with population sizes specified for sources in ",param.infile,". Exiting....",sep=''));q(save='no')}
drift.vec=as.double(line.check(param.infile,8,"drift.btwn.surrogates.and.sources:"))
if ((length(drift.vec)!=length(sel.simpops))||sum(drift.vec>0,na.rm=T)!=length(drift.vec)){print(paste("Something wrong with drift values specified in ",param.infile,". Exiting....",sep=''));q(save='no')}
num.neutral.snps=as.character(line.check(param.infile,9,"num.neutral.snps:"))    ## 'all' or lower if trying to save time
if (num.neutral.snps!="all") num.neutral.snps=as.integer(num.neutral.snps)
if ((length(num.neutral.snps)!=1)||(num.neutral.snps!="all" && num.neutral.snps<=0)){print(paste("Something wrong with number of neutral SNPs specified in ",param.infile,". Exiting....",sep=''));q(save='no')}
startfreq.selsnp.range=as.double(line.check(param.infile,10,"range.startfrequency.selected.snp:"))    ## (0.5) selected snps will have starting frequency >= first value and <= second value (i.e. with selection forcing the frequency towards 1.0 over time)
if ((length(startfreq.selsnp.range)!=2)||min(startfreq.selsnp.range)<0||max(startfreq.selsnp.range)>1){print(paste("Something wrong with starting frequency range for selected SNP specified in ",param.infile,". Exiting....",sep=''));q(save='no')}
infer.source.allelefreqs=as.integer(line.check(param.infile,11,"infer.source.freq.using.target.data:"))           ## INFER ALLELE FREQS OF SOURCES USING ONLY TARGET IND GENOTYPES AND ADMIXTURE PROPS (TO HELP DECIDE WHAT DRIFT TO USE)
inds.perrun=as.integer(line.check(param.infile,12,"divide.into.runs.ofXX.inds(to.reduce.RAM):"))           ## TO CONTROL RAM (BUT TAKES LONGER)
if ((length(inds.perrun)!=1)||inds.perrun<=0){print(paste("Something wrong with number of individuals specified (to reduce RAM) in ",param.infile,". Exiting....",sep=''));q(save='no')}
sel.prob=rep(sel.prob.val,1)   ## prob of being selected each generation (also determines number of selected SNPs to simulate)

genotypes.filein=as.character(temp[3])
id.file=as.character(temp[4])
out.filePRE=as.character(temp[5])
out.filePOSThaps=".haps"
out.filePOSTid=".idfile.txt"
out.filePOSTtruth=".truth.txt"

##############################
## PROGRAM:

Mb=10^6
options(scipen=999)
if (sel.postadmix==1)
{
	if (sum(sel.simpops==1)!=length(sel.simpops)){print("WARNING: You have specified selection post-admixture, but not to have all sources having pre-admixture selection. For post-admixture selection, this simulator works (not ideally) by simulating selection independently in each source (each for the mean of your specified \"generations.selection.each.source\" generations) and then mixing them to make admixed pop.")}
	sel.simpops=rep(1,length(source.pops))
	gen.sel.source=rep(mean(gen.sel.source),length(source.pops))
}

                                   ## (I) GET ADMIXTURE PROPS/INDS:
id.mat=read.table(id.file,as.is=TRUE)
if (dim(id.mat)[2]!=length(source.pops)+2){print(paste("Something wrong with file ",id.file," -- should have columns equal to number of source populations + 2. Exiting...",sep=''));q(save='no')}
admixture.props.all=matrix(as.matrix(id.mat[,3:dim(id.mat)[2]]),ncol=dim(id.mat)[2]-2)
rownames(admixture.props.all)=id.mat[,1]
colnames(admixture.props.all)=source.pops
for (i in 1:length(target.pop)) if (sum(is.element(id.mat[,2],target.pop[i]))==0) {print(paste("No individuals found in target pop ",target.pop[i],"! Exiting....",sep=''));q(save='no')}
for (i in 1:length(source.pops)) if (sum(is.element(id.mat[,2],source.pops[i]))==0) {print(paste("No individuals found in surrogate pop ",source.pops[i],"! Exiting....",sep=''));q(save='no')}
admixture.props=admixture.props.all[match(id.mat[is.element(id.mat[,2],target.pop) & is.element(id.mat[,1],rownames(admixture.props.all)),1],rownames(admixture.props.all)),]
admixture.props=admixture.props/apply(admixture.props,1,sum)
is.element.vec=is.element(id.mat[,2],source.pops)
is.element.target=is.element(id.mat[,2],target.pop)

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

                                   ## (III) GET SOURCE POP & TARGET ALLELE FREQS:
                                          ## get total number of snps:
nsites.vec=rep(NA,num.genofiles)
readfileNAMES=file(genotypes.filein,open="r")
for (i in 1:num.genofiles)
{
	filename.i=as.character(scan(readfileNAMES,nlines=1,what='char',quiet=TRUE))
	readfile=gzfile(filename.i,open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## ninds
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)    ## positions
	pos.vec=as.character(line2[2:length(line2)])
	nsites.vec[i]=length(pos.vec)
	close(readfile)
}
close(readfileNAMES)
                                           ## get positions and allele freqs in surrogates and targets:
chromo.vec.all=pos.vec.all=rep(NA,sum(nsites.vec))
allelefreq.all=nind.all=matrix(0,nrow=length(source.pops),ncol=sum(nsites.vec))
snp.count=0
readfileNAMES=file(genotypes.filein,open="r")
for (i in 1:num.genofiles)
{
	filename.i=as.character(scan(readfileNAMES,nlines=1,what='char',quiet=TRUE))
	print(c(i,num.genofiles,filename.i))
	readfile=gzfile(filename.i,open='r')
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nhaps
	nind=as.integer(line2[1])/2
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
	line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## positions
	chromo.vec.all[(snp.count+1):(snp.count+nsites.vec[i])]=i
	pos.vec.all[(snp.count+1):(snp.count+nsites.vec[i])]=as.double(line2[2:length(line2)])
	for (h in 1:nind)
	{
		hap1=scan(readfile,nlines=1,what='char',quiet=TRUE)
		hap2=scan(readfile,nlines=1,what='char',quiet=TRUE)
		if (is.element.vec[h])
		{
			hap1=as.integer(strsplit(hap1,split='')[[1]])
			hap2=as.integer(strsplit(hap2,split='')[[1]])
			hap.tot=hap1+hap2
			allelefreq.all[source.pops==id.mat[h,2],(snp.count+1):(snp.count+nsites.vec[i])][hap.tot>=0]=allelefreq.all[source.pops==id.mat[h,2],(snp.count+1):(snp.count+nsites.vec[i])][hap.tot>=0]+hap.tot[hap.tot>=0]
			nind.all[source.pops==id.mat[h,2],(snp.count+1):(snp.count+nsites.vec[i])][hap.tot>=0]=nind.all[source.pops==id.mat[h,2],(snp.count+1):(snp.count+nsites.vec[i])][hap.tot>=0]+2
		}
	}
	close(readfile)
	snp.count=snp.count+nsites.vec[i]
}
close(readfileNAMES)
small.val=0.0001
start.freq.mat=matrix(NA,nrow=length(source.pops),ncol=sum(nsites.vec))
for (k in 1:length(source.pops))
{
	####allelefreq.all[k,]=allelefreq.all[k,]/(2*sum(id.mat[is.element.vec,2]==source.pops[k]))
	allelefreq.all[k,]=allelefreq.all[k,]/nind.all[k,]
	allelefreq.all[k,][allelefreq.all[k,]<small.val]=small.val
	allelefreq.all[k,][allelefreq.all[k,]>(1-small.val)]=1-small.val
	if (drift.vec[k]==0) start.freq.mat[k,]=allelefreq.all[k,]
	if (drift.vec[k]>0) start.freq.mat[k,]=rbeta(sum(nsites.vec),allelefreq.all[k,]*(1-drift.vec[k])/drift.vec[k],(1-allelefreq.all[k,])*(1-drift.vec[k])/drift.vec[k])
}

if (sel.prob.val!=0)
{
                                   ## (IV) RANDOMLY PICK SNPS FOR SIMS:
                                          ## randomly pick selected SNPs:
	#if (sel.postadmix==0) mean.freq=start.freq.mat[as.integer(sample(as.character(1:length(source.pops))[sel.simpops==1],1)),]
	if (sel.postadmix==0) mean.freq=apply(matrix(start.freq.mat[sel.simpops==1,],nrow=sum(sel.simpops==1)),2,max)
	######if (sel.postadmix==0) possible.snps=(1:length(mean.freq))[mean.freq>=startfreq.selsnp.range[1] & mean.freq<=startfreq.selsnp.range[2] & apply(matrix(start.freq.mat,nrow=dim(start.freq.mat)[1]),2,min)>=startfreq.selsnp.range[1]]
	if (sel.postadmix==1)
	{
		mean.freq=rep(NA,dim(start.freq.mat)[2])
		for (j in 1:dim(start.freq.mat)[2]) mean.freq[j]=mean(matrix(start.freq.mat[,j],nrow=1)%*%t(admixture.props))
		#####mean.freq=apply(t(start.freq.mat)%*%t(admixture.props),1,mean)
	}
	possible.snps=(1:length(mean.freq))[mean.freq>=startfreq.selsnp.range[1] & mean.freq<=startfreq.selsnp.range[2] & apply(matrix(start.freq.mat,nrow=dim(start.freq.mat)[1]),2,min)>=startfreq.selsnp.range[1]]
	if (length(possible.snps)<length(sel.prob))
	{
		print("Not enough SNPs meeting selection criterion! Exiting...")
		q(save='no')
	}
	selected.snps=sort(as.integer(sample(as.character(possible.snps),length(sel.prob),replace=F)))
	##selected.snps=NULL
	##while(length(selected.snps)<length(sel.prob))
	##{
	##	print(c(length(selected.snps),length(sel.prob)))
	##	selected.snp.propose=as.integer(sample(setdiff(as.character(1:dim(start.freq.mat)[2]),as.character(selected.snps)),1,replace=F))
	##	if (mean(matrix(start.freq.mat[,selected.snp.propose],nrow=1)%*%t(admixture.props))<max.startfreq.selsnp) selected.snps=c(selected.snps,selected.snp.propose)
	##}
	##selected.snps=sort(selected.snps)

                                   ## (V) GENERATE SOURCES' SIM DATA AT SELECTED SNPS:
        start.freq.mat.selected=matrix(start.freq.mat[,selected.snps],ncol=length(selected.snps))
	for (k in 1:length(source.pops))
	{
		if (sel.simpops[k]==1)
		{
			for (s in 1:length(selected.snps))
			{
				current.pop=rbinom(pop.size.vec[k],2,start.freq.mat.selected[k,s])
				for (g in 1:gen.sel.source[k])
				{
					prob.vec=rep(1,pop.size.vec[k])
					if (sel.type=="additive") prob.vec=1+sel.prob[s]*current.pop
					if (sel.type=="recessive") prob.vec[current.pop==2]=1+sel.prob[s]
					if (sel.type=="dominant") prob.vec[current.pop>=1]=1+sel.prob[s]
					if (sel.type=="multiplicative") prob.vec=(1+sel.prob[s])^current.pop
					current.pop=sample(current.pop,pop.size.vec[k],prob=prob.vec,replace=T)
				}
				start.freq.mat.selected[k,s]=sum(current.pop)/(2*length(current.pop))
			}
		}
	}
	#print(matrix(start.freq.mat[,selected.snps],ncol=length(selected.snps)))
	#print(start.freq.mat.selected)
}
if (sel.prob.val==0) selected.snps=integer(0)

                                          ## randomly pick neutral SNPs:
if (sel.prob.val!=0)
{
	if (num.neutral.snps!='all') neutral.snps=sort(as.integer(sample(as.character(1:sum(nsites.vec))[-selected.snps],num.neutral.snps,replace=F)))
	if (num.neutral.snps=='all') neutral.snps=(1:sum(nsites.vec))[-selected.snps]
}
if (sel.prob.val==0)
{
	if (num.neutral.snps!='all') neutral.snps=sort(as.integer(sample(as.character(1:sum(nsites.vec)),num.neutral.snps,replace=F)))
	if (num.neutral.snps=='all') neutral.snps=(1:sum(nsites.vec))
}

                                   ## (VI) PRINT OUT INFO + SURROGATE DATA:
                                            ## print truth info:
write(out.filePRE,file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1)
write(paste("sel.coeff: ",sel.prob.val,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)
write(paste("selection.post-admixture?: ",sel.postadmix,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)
write(paste("sel.type: ",sel.type,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)
write(paste("target.pop: ",target.pop,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)
write(c("surrogate.pops:",source.pops),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=length(source.pops)+1,append=TRUE)
write(c("drift.btwn.surrogates.and.sources:",drift.vec),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=length(source.pops)+1,append=TRUE)
write(c("sources.with.selection.preadmixture:",sel.simpops),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=length(source.pops)+1,append=TRUE)
write(c("generations.selection.each.source:",gen.sel.source),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=length(source.pops)+1,append=TRUE)
write(c("pop.size.sources:",pop.size.vec),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=length(source.pops)+1,append=TRUE)
write(paste("num.neutral.snps: ",num.neutral.snps,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)
write(paste("range.startfrequency.selected.snp: ",startfreq.selsnp.range,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=2,append=TRUE)
write(paste("id.file: ",id.file,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)

write(paste("genotypes.file: ",genotypes.filein,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)
write(paste("id.admixture.props.file: ",id.file,sep=''),file=paste(out.filePRE,out.filePOSTtruth,sep=''),ncolumns=1,append=TRUE)
                                            ## print new id-file:
write.table(id.mat[is.element.vec,],file=paste(out.filePRE,out.filePOSTid,sep=''),quote=FALSE,row.names=F,col.names=F)
write.table(id.mat[match(rownames(admixture.props),id.mat[,1]),],file=paste(out.filePRE,out.filePOSTid,sep=''),quote=FALSE,row.names=F,col.names=F,append=T)
                                            ## print header for hap-data:
write.table(2*(sum(is.element.vec)+dim(admixture.props)[1]),file=paste(out.filePRE,out.filePOSThaps,sep=''),row.names=F,col.names=F,quote=F)  ## nhaps
write.table(length(neutral.snps)+length(selected.snps),file=paste(out.filePRE,out.filePOSThaps,sep=''),row.names=F,col.names=F,quote=F,append=TRUE)  ## nsites
if (sel.prob.val!=0) write(c("P",paste(chromo.vec.all[neutral.snps],":",pos.vec.all[neutral.snps],sep=''),paste(chromo.vec.all[selected.snps],":",pos.vec.all[selected.snps],"SEL",sep='')),file=paste(out.filePRE,out.filePOSThaps,sep=''),ncolumns=length(neutral.snps)+length(selected.snps)+1,append=TRUE)  ## positions
if (sel.prob.val==0) write(c("P",paste(chromo.vec.all[neutral.snps],":",pos.vec.all[neutral.snps],sep='')),file=paste(out.filePRE,out.filePOSThaps,sep=''),ncolumns=length(neutral.snps)+1,append=TRUE)  ## positions
                                            ## print source hap-data:
num.runs=ceiling(sum(is.element.vec)/inds.perrun)
ind.labels=(1:length(is.element.vec))[is.element.vec]
for (a in 1:num.runs)
{
	start.ind=(a-1)*inds.perrun+1
	end.ind=a*inds.perrun
	if (end.ind>sum(is.element.vec)) end.ind=sum(is.element.vec)
	source.haps=matrix(NA,nrow=2*(end.ind-start.ind+1),ncol=length(neutral.snps)+length(selected.snps))
	is.element.vec.a=rep(F,length(is.element.vec))
	is.element.vec.a[ind.labels[start.ind:end.ind]]=T
	readfileNAMES=file(genotypes.filein,open="r")
	for (i in 1:num.genofiles)
	{
		filename.i=as.character(scan(readfileNAMES,nlines=1,what='char',quiet=TRUE))
		print(c(a,num.runs,start.ind,end.ind,sum(is.element.vec.a),i,num.genofiles,filename.i))
		readfile=gzfile(filename.i,open='r')
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nhaps
		nind=as.integer(line2[1])/2
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
		line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## positions
		snps.toprint=match(paste(chromo.vec.all[c(neutral.snps,selected.snps)],":",pos.vec.all[c(neutral.snps,selected.snps)],sep=''),paste(i,":",as.double(line2[2:length(line2)]),sep=''))
		snps.toprint=snps.toprint[!is.na(snps.toprint)]
		match.snps=match(paste(i,":",as.double(line2[2:length(line2)]),sep='')[snps.toprint],paste(chromo.vec.all[c(neutral.snps,selected.snps)],":",pos.vec.all[c(neutral.snps,selected.snps)],sep=''))
		ind.count=0
		for (h in 1:nind)
		{
			hap1=scan(readfile,nlines=1,what='char',quiet=TRUE)
			hap2=scan(readfile,nlines=1,what='char',quiet=TRUE)
			if (is.element.vec.a[h])
			{
				ind.count=ind.count+1
				source.haps[2*ind.count-1,match.snps]=as.integer(strsplit(hap1,split='')[[1]])[snps.toprint]
				source.haps[2*ind.count,match.snps]=as.integer(strsplit(hap2,split='')[[1]])[snps.toprint]
			}
		}
		close(readfile)
	}
	close(readfileNAMES)
	for (i in 1:dim(source.haps)[1]) write(source.haps[i,],file=paste(out.filePRE,out.filePOSThaps,sep=''),ncolumns=dim(source.haps)[2],append=TRUE,sep='')
}

                                            ## (VII) GENERATE AND PRINT TARGET DATA:
rbeta.func=function(x,drift.val)
{
	x[x>1]=1
	x[x<0]=0
	return(rbeta(length(x),x*(1-drift.val)/drift.val,(1-x)*(1-drift.val)/drift.val))
}
rbinom.func=function(x,count.val) return(rbinom(length(x),size=count.val,prob=x))
for (i in 1:dim(admixture.props)[1])
{
                                    ## neutral snps:
	mean.freq.mat=t(start.freq.mat[,neutral.snps])%*%matrix(admixture.props[i,],ncol=1)
	####prob.mat.neutral=rbeta.func(mean.freq.mat,drift.val=drift.target)  ## !!!!!!! BECAUSE DRIFT SHOULD ALREADY HAVE BEEN INCORPORATED WHEN SAMPLING DONOR FREQS
	prob.mat.neutral=mean.freq.mat
	prob.mat.neutral[prob.mat.neutral<0]=0
	prob.mat.neutral[prob.mat.neutral>1]=1
	#print(c(i,dim(admixture.props)[1],sum(prob.mat.neutral<0),sum(prob.mat.neutral>1)))
	geno.counts.neutral=rbinom.func(prob.mat.neutral,count.val=2)
	if (sel.prob.val!=0)
	{
                                    ## selected snps:
		mean.freq.mat=t(start.freq.mat.selected)%*%matrix(admixture.props[i,],ncol=1)
		mean.freq.mat[mean.freq.mat<0]=0
		mean.freq.mat[mean.freq.mat>1]=1
		if (sel.postadmix==0)
		{
			####prob.mat.selected=rbeta.func(mean.freq.mat,drift.val=drift.target)    ## !!!!!! SAME
			prob.mat.selected=mean.freq.mat
			prob.mat.selected[prob.mat.selected<0]=0
			prob.mat.selected[prob.mat.selected>1]=1
			geno.counts.selected=rbinom.func(prob.mat.selected,count.val=2)
		}
		if (sel.postadmix==1) geno.counts.selected=rbinom.func(mean.freq.mat,count.val=2)
	}
	if (sel.prob.val==0) geno.counts.selected=NULL
	
                                     ## print-out:
	geno.i=c(geno.counts.neutral,geno.counts.selected)
	hap.i1=geno.i
	hap.i1[hap.i1==2]=1
	hap.i2=geno.i-1
	hap.i2[hap.i2<0]=0
	write(hap.i1,file=paste(out.filePRE,out.filePOSThaps,sep=''),ncolumns=length(hap.i1),append=TRUE,sep='')
	write(hap.i2,file=paste(out.filePRE,out.filePOSThaps,sep=''),ncolumns=length(hap.i2),append=TRUE,sep='')
}

                                            ## (VIII) IF SPECIFIED, INFER SOURCE ALLELE FREQUENCIES USING TARGET DATA AND ADMIXTURE PROPS:
admixture.est.lm=function(x,admix.props,allele.counts,total.counts)
{
	p.x=admix.props%*%matrix(x,ncol=1)
	loglik.val=sum(allele.counts*log(p.x)+(total.counts-allele.counts)*log(1-p.x))
	return(-loglik.val)
}
if (infer.source.allelefreqs==1)
{
	freq.mat.nlminb=matrix(NA,nrow=length(source.pops),ncol=sum(nsites.vec))
	rownames(freq.mat.nlminb)=source.pops

	snps.perrun=10000
 	snp.count=0
	readfileNAMES=file(genotypes.filein,open="r")
	for (i in 1:num.genofiles)
	{
		filename.i=as.character(scan(readfileNAMES,nlines=1,what='char',quiet=TRUE))
		print(c(i,num.genofiles,filename.i))

 		num.runs=ceiling(nsites.vec[i]/snps.perrun)
		for (a in 1:num.runs)
		{
			start.snp=(a-1)*snps.perrun+1
			end.snp=a*snps.perrun
			if (end.snp>nsites.vec[i]) end.snp=nsites.vec[i]
			
                                                    ## read in data:
			readfile=gzfile(filename.i,open='r')
			line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nhaps
			nind=as.integer(line2[1])/2
			line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## nsites
			line2=scan(readfile,nlines=1,what='char',quiet=TRUE)   ## positions

			genotypes.target=matrix(NA,nrow=sum(is.element.target),ncol=end.snp-start.snp+1)
			target.count=1
			for (h in 1:nind)
			{
				hap1=scan(readfile,nlines=1,what='char',quiet=TRUE)
				hap2=scan(readfile,nlines=1,what='char',quiet=TRUE)
				if (is.element.target[h])
				{
					hap1=as.integer(strsplit(hap1,split='')[[1]])[start.snp:end.snp]
					hap2=as.integer(strsplit(hap2,split='')[[1]])[start.snp:end.snp]
					hap.tot=hap1+hap2
					genotypes.target[target.count,]=hap.tot
					target.count=target.count+1
				}
			}
			close(readfile)

                                                     ## infer allele freqs:
			for (j in 1:(end.snp-start.snp+1)) freq.mat.nlminb[,snp.count+start.snp+j-1]=nlminb(start=rep(0.5,length(source.pops)),admixture.est.lm,lower=0,upper=1,admix.props=admixture.props[genotypes.target[,j]>=0,],allele.counts=genotypes.target[,j][genotypes.target[,j]>=0],total.counts=rep(2,sum(genotypes.target[,j]>=0)))$par
		}
		snp.count=snp.count+nsites.vec[i]
	}
	close(readfileNAMES)
}

          ## SUMMARY INFO:
                   ## allele freqs at neutral/selected snps per source:
apply(round(start.freq.mat[,neutral.snps],4),1,summary)
if (sel.prob.val!=0)
{
	print(apply(round(start.freq.mat.selected,4),1,summary))
	print(apply(round(start.freq.mat[,setdiff(possible.snps,selected.snps)],4),1,summary))
	print(apply(start.freq.mat.selected-matrix(start.freq.mat[,selected.snps],ncol=length(selected.snps)),1,summary))
}
                   ## calculate correlations and print drift.vals:
driftoutput <- NULL
print.mat=c("source","drift.val","cor.sims","cor.truth")
for (i in 1:length(source.pops)){
  row <- c(source.pops[i],as.numeric(drift.vec[i]),round(cor(start.freq.mat[i,neutral.snps],allelefreq.all[i,neutral.snps]),4))
  driftoutput <- rbind(driftoutput,row)
  if (infer.source.allelefreqs!=1) print.mat=rbind(print.mat,c(source.pops[i],drift.vec[i],round(cor(start.freq.mat[i,neutral.snps],allelefreq.all[i,neutral.snps]),4),NA))
  if (infer.source.allelefreqs==1) print.mat=rbind(print.mat,c(source.pops[i],drift.vec[i],round(cor(start.freq.mat[i,neutral.snps],allelefreq.all[i,neutral.snps]),4),round(cor(freq.mat.nlminb[i,neutral.snps],allelefreq.all[i,neutral.snps]),4)))
}
rownames(print.mat)=rep("",dim(print.mat)[1])
colnames(print.mat)=rep("",dim(print.mat)[2])
driftoutput <- as.data.frame(t(driftoutput))
colnames(driftoutput) <- source.pops
driftoutput <- driftoutput[2,]

warnings()   ## can ignore those related to nlminb()



noquote(print.mat)
	## "cor.sims" is the correlation between the simulated "source" allele frequencies and its surrogate (at neutral SNPs)
	## "cor.truth" is the correlation between the "source" allele frequencies inferred using the target and its surrogate (at neutral SNPs)
	## the idea is to toggle "drift.val" until you get "cor.sims" ~= "cor.truth"

q(save='no')

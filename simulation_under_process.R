
#AUTHOR: Yaniv Brandvain
#DESCRIPTION: This script was written to simulate recombining chromosomes in finite populations. This works by tracing the ancestry of each chromosome portion in an admixed population to the initial populaion. Each generation consists of SELECTION, RANDOM MATING (WITH RECOMBINATION) and MIGRATION and therefore accounts for population genealogy and drift. Generations are discrete. 

#In Sedghifar et al. (2015), The Spatial Mixing of Genomes in Secondary Contact Zones , this script was used to simulate admixture (under the process) along a neutral contact zone. 

#CAUTION: This code comes with no guarantees, and may contain bugs. Please let us know if you find any. 

library(multicore)

###%%%%SET THE PARAMETERS%%%%%%
Number.Of.Generations = 5 
Total.Individuals = 200000
Number.Of.Demes = 20
Deme.Size = Number.Of.Individuals/Number.Of.Demes
Migration.Rate = 0.01
N.Chrom = 1


###FUNCTIONS:

	getRuns = function(index){
	    #This returns an array in which rows are the the beggingns and end of a TRUE/FALSE vector
	    starts  = which(c(0,index)[-1]==1 & c(0,index)[-length(index)]!=1)
	    stops   = which(c(index,0)[-1]!=1 & c(index,0)[-(1+length(index))]==1)
	    if(index[length(index)] == 1 & index[(length(index)-1)] !=1){starts = c(starts,length(index))}
	    return (cbind(starts,stops))
	}


	namesList = function(these.names){
		# this is a way to get a list where each elelment has a name equal to its value, so that lapply returns a correctly named list
		all.things = as.list(these.names)
		names(all.things)=these.names
		return(all.things)
	}
	

	#get breaks between species' ancestry
	spBreaks = function(sim){
		ind.ancest = lapply( sim$inds,function(IND){ lapply(IND,function(CHR){lapply(CHR,function(PAR){
				if(nrow(PAR)==1){return(data.frame(starts=0,stops=1,sp2 = PAR$ancest > length(sim[[1]])))}
				sp1 = getRuns(as.numeric( PAR$ancest <= length(sim[[1]])))
				sp2 = getRuns(as.numeric( PAR$ancest > length(sim[[1]])))
				tmp = rbind(data.frame(sp1,sp2=rep(F,nrow(sp1))),data.frame(sp2,sp2=rep(T,nrow(sp2))))
				tmp = tmp[order(tmp$starts),]
				tmp$starts = PAR$bp[tmp$starts]
				tmp$stops = c(PAR$bp[-1],1)[tmp$stops]
				return(tmp)
		})})})
		ancest.prop = t(sapply(ind.ancest,function(IND){
			tmp = t(sapply(IND,function(CHR){ c(X1=with(CHR$X1,sum(sp2*(stops-starts))), X2 = with(CHR$X2,sum(sp2*(stops-starts))))  }))
			return(colSums(tmp)/nrow(tmp))
		}))		
		return(list(ind.ancest = ind.ancest, ancest.prop = ancest.prop))
	}


	initializeInds = function(arr,n.chr=n.chr){	
		return(lapply(seq_along(arr[,1]),function(IND){
			X=replicate(n.chr,list( X1=data.frame(bp=0,ancest = arr[IND,"X1"]), 
				X2=data.frame(bp=0,ancest = arr[IND,"X2"])  ),simplify=FALSE)
			names(X) = paste("chr",c(1:n.chr),sep="")
			return(X)
		}))
	}


#######  GOING THROUGH MEIOSIS    #######
	meiosis = function(CHR,r=1){
		start = sample(c(0,1),1)
		breaks = sort(runif(rpois(1,1))) #Place recombination events on a chromsome of 1M length. Number of rec. events drawn from Poisson distribution with mean 1.
		A= breaks[seq_along(breaks)%%2==start]
		B= breaks[seq_along(breaks)%%2!=start]
		X = rbind(
			do.call(rbind,lapply(A,function(a){
				with(CHR[[1]],c(bp=a,ancest=ancest[max(which(bp<a))]))})),
				CHR[[1]][(findInterval(CHR[[1]][,"bp"],c(-1,breaks,1.2))%%2)!=start,],
				do.call(rbind,lapply(B,function(b){with(CHR[[2]],c(bp=b,ancest=ancest[max(which(bp<b)) ]))})),
				CHR[[2]][(findInterval(CHR[[2]][,"bp"],c(-1,breaks,1.2))%%2)==start,]
			)
		return(		X[order(X$bp),]		)
	}





#######  GENOTYPING INDIVIDUALS   #######
	genotype = function(IND.CH,L,sp){
		do.call(cbind,lapply(IND.CH,function(PAR){	
			ancest.chunk = with(PAR,ancest[findInterval( L, c(bp,1) )]  )
			sp[ancest.chunk]  
		}))
	} 	
	genoAll = function(inds,QTL,sp.ids){
		tmp.ch = namesList(names(QTL))
		sp=rep(sp.ids,each=2)
		lapply(inds,function(IND){
			lapply(tmp.ch,function(CHR){
				genotype(IND[[CHR]],QTL[[CHR]]$bp,sp)
			})
		})
	}


#######    MAKING A NEW GEN     #######

	makeNewGen=function(inds,sp.ids,demes,mig,edge="M/2"){
		w_abs = rep(1,length(inds)) 
		rel.w = w_abs/mean(w_abs)
		inds.df = data.frame(w=rel.w,deme=demes,edge=ifelse(demes%in%range(demes),2,1))
		move = sign(rbinom(nrow(inds.df),1,prob=with(inds.df,ifelse(deme==min(deme),1,ifelse(deme!=max(deme),1/2,0))))-.5) * with(inds.df, rbinom(nrow(inds.df),1, c(mig/edge)))
		inds.df$new.deme = inds.df$deme + move
		inds.df$id = seq_along(inds.df$new.deme)
		new.inds.df = do.call(rbind,lapply(  unique(inds.df$deme)  , function(d){			
			#find parents
			a = with(inds.df[inds.df$new.deme==d,], sample(id, prob=w/sum(w) , replace=T))
			# MODIFY B|A TO ALLOW FOR ASSORTATIVE MATING & eg no selfing
			b = with(inds.df[inds.df$new.deme==d,], sample(id, prob=w/sum(w) , replace=T))
			return(cbind(a,b))
		}))
		lapply(seq_along(new.inds.df[,1]),function(NEW){
			NEW = lapply(list(X1 = inds[[new.inds.df[NEW,1]]],X2 = inds[[new.inds.df[NEW,2]]]),function(PAR){  lapply(PAR,meiosis)   })
			lapply(namesList(names(NEW[[1]])),function(C){lapply(NEW,function(PARNTL){PARNTL[[C]]})})
		})
	}

#This function calculates genome-wide ancetry proportions for each ind at the end of the simulation. It is used as a sanity checker in function 'forwardSim'
	getAlphas = function(inds,sp.ids){
		sp=rep(sp.ids,each=2)
		do.call(rbind,lapply(inds,function(IND){
			X=do.call(rbind,lapply(IND,function(CHR){
				do.call(rbind,lapply(CHR,function(PAR){
					sp.origin=factor(sp[PAR$ancest],levels=unique(sp))
					tapply(diff(c(PAR$bp,1)),sp.origin,sum)
				}))
			}))
			X[is.na(X)]=0
			colMeans(X)  
		}))
	}


#This function 

	forwardSim = function(n,sp.ids,deme,n.gen,n.chr,mig,SAVE=NULL,RETURN=TRUE,...){
		initial.inds.df = data.frame(matrix(1:(2*(n)),ncol=2,byrow=TRUE),sp.ids=sp.ids,deme=deme)
		initial.inds = initializeInds(arr=initial.inds.df,n.chr=n.chr)
		inds = initial.inds
		inds.demes = initial.inds.df$deme
		for(i in 1:n.gen){
			inds = makeNewGen(inds,sp.ids=sp.ids,demes=inds.demes,mig=mig)  # add in selection and assortative mating
			#plot(getAlphas(inds,sp.ids)[,1])
			print(c(i,date()))
			if(!is.null(SAVE)){ if(i%in%SAVE$gen){  
					sim = list(inds=inds, sp.ids=sp.ids,deme=deme,pars = c(mig=mig, n.gen=n.gen))
					save(sim,file=sprintf("%s_t%i_g%i.Robj",SAVE$name,SAVE$try,i))
			}}
		}
		if(RETURN){return(list(inds=inds, sp.ids=sp.ids,deme=deme,pars = c(mig=mig, n.gen=n.gen) ))}
		return(NULL)
	}		




sim1 = 	forwardSim(n=Total.Individuals, sp.ids = rep(c("A","B"),each=Total.Individuals/2),
	deme = rep(1:Number.Of.Demes,each=Deme.Size),n.chr=N.Chrom,n.gen=Number.Of.Generations,mig =Migration.Rate)
sim1.sum = spBreaks(sim1)
save(sim1, file = sprintf("m%s_tau%s.snp",Migration.Rate,Number.Of.Generations))



######	Process output and write to file  #######
outfile.snp = sprintf("sim_m%s_tau%s.snp",Migration.Rate,Number.Of.Generations)
outfile.geno = sprintf("sim_m%s_tau%s.geno",Migration.Rate,Number.Of.Generations)

#CHOOSE LOCI
    loci= sort(runif(1000))

#MAKE file.snp
    info = data.frame(cbind(1,round(loci*1000000),0,loci))
    colnames(info) = NULL
    write.csv(info,row.names=F,file=outfile.snp)


#GENOTYPE  [ASSUMES A ONE CHR GENOME WHERE WE HAVE PHASED DIPLOID DATA]
    genoInd = function(IND,loci){sapply(IND[[1]],function(CHR){as.numeric(CHR$sp2[as.numeric(cut(loci,c(CHR$starts,1)))])})}

#GENOTYPE ALL MY INDS
    my.genos = data.frame(do.call(cbind,lapply(sim1.sum$ind.ancest, genoInd, loci)))
    colnames(my.genos) = NULL

#MAKE OUR FILE [NOTE THIS TAKES A LONG TIME AND IT MIGHT BE SMARTER TO WRITE THE THING ABOVE TO FILE ONE LINE AT A TIME]
    write.table(my.genos, file = outfile.geno,row.names = FALSE,sep="")


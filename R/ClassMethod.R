setClass("kaksCodon",
	representation(
		seq_num="numeric",
		kaks="numeric",
		lod="numeric",
		q="numeric",
		ka="numeric",
		ks="numeric"
	)
)

setClass("kaksAA",
	representation(
		seq_num="numeric",
		kaks="list",
		lod="list",
		q="list",
		ka="list",
		ks="list"
	)
)
setClass("MI",
	representation(
		mi="matrix",
		p.value="matrix"
	)
)

setClass("ckaks",
	representation(
		ckaks="matrix",
		lod="matrix"
	)
)

setClass("biCompare",
	representation(
		method = "character",
		state_1="matrix",
		statistic_1="matrix",
		state_2="matrix",
		statistic_2="matrix",
		positiveSite01="ANY",
		positiveSite02="ANY"
	)
)
####################################################################
setGeneric("filterSites", 
	function (x, ...)
	standardGeneric("filterSites")
)
####################################################################
setMethod("filterSites",signature(x="kaksCodon"),
	function(x,lod_cut = 2,freq_cut=0.01){
		prange=1:length(x@kaks)
		condition=((x@kaks>1)+(x@lod>=lod_cut)+((x@ka/x@seq_num)>=freq_cut))==3
		mutation=prange[condition]
		kaks=x@kaks[condition]
		lod=x@lod[condition]
		mutFreq=sapply((x@ka[condition])/x@seq_num,function(xx)round(xx,digits=4))
		fdataframe=data.frame(mutation=mutation,kaks=kaks,lod=lod,mutFreq=mutFreq)
		names(fdataframe)=c("mutation","kaks","lod","freq")
		return(fdataframe)
	}
)

setMethod("filterSites",signature(x="kaksAA"),
	function(x,lod_cut = 2,freq_cut=0.01){
		positions=1:length(x@kaks)
		fdataframe=c()
		for(i in positions){
			positionnames=names(x@kaks[[i]])
			if(!is.null(positionnames)){
				pkaks=x@kaks[[i]]
				for(j in 1:length(pkaks)){
					mutfreq=round((x@ka[[i]][names(pkaks[j])])/x@seq_num,digits=4)
					if(pkaks[j]>1 & x@lod[[i]][names(pkaks[j])]>=lod_cut & mutfreq >= freq_cut){
						fdataframe=rbind(fdataframe,c(i,names(pkaks[j]),pkaks[j],x@lod[[i]][names(pkaks[j])],mutfreq))
					}
				}
			}
		}
		fdataframe=as.data.frame(fdataframe)
		names(fdataframe)=c("position","mutation","kaks","lod","freq")
		return(fdataframe)
	}
)
setMethod("filterSites",signature(x="MI"),
	function(x,p_cut=0.05){
		fdataframe=c()
		datanames=rownames(x@mi)
		for(i in 2:dim(x@p.value)[1]){
			for(j in 1:(i-1)){
				if(!is.na(x@p.value[i,j])){
					if(x@p.value[i,j]<=p_cut){
						fdataframe=rbind(fdataframe,c(datanames[i],datanames[j],x@mi[i,j],x@p.value[i,j]))
					}
				}
			}
		}
		fdataframe=as.data.frame(fdataframe)
		names(fdataframe)=c("position1","position2","mi","p.value")
		return(fdataframe)
}
)

setMethod("filterSites",signature(x="ckaks"),
	function(x,lod_cut=2){
		fdataframe=c()
		datanames=rownames(x@ckaks)
		for(i in 1:dim(x@lod)[1]){
			for(j in 1:dim(x@lod)[2]){
				if(!is.na(x@lod[i,j])){
					if((x@lod[i,j]>=lod_cut)&(x@ckaks[i,j]>1)){
						fdataframe=rbind(fdataframe,c(datanames[i],datanames[j],x@ckaks[i,j],x@lod[i,j]))
					}
				}
			}
		}
		if(is.null(fdataframe)) stop("No suitable results!")
		fdataframe=as.data.frame(fdataframe)
		names(fdataframe)=c("position1","position2","ckaks","lod")
		return(fdataframe)
}
)

weightnorm=function(x){
	return((x-min(x))*2/(max(x)-min(x)))
}

setMethod("plot",signature(x="MI"),
	function(x){
	filteredata=filterSites(x)
	targraph=graph.edgelist(as.matrix(filteredata[,1:2]),directed=F)
	E(targraph)$weight=filteredata$mi
	plot(targraph,vertex.label=V(targraph)$name,edge.width=weightnorm(E(targraph)$weight),layout=layout.fruchterman.reingold,vertex.label.cex=.7)
}
)

setMethod("plot",signature(x="ckaks"),
	function(x){
	filteredata=filterSites(x)
	targraph=graph.edgelist(as.matrix(filteredata[,1:2]),directed=T)
	E(targraph)$weight=filteredata$ckaks
	plot(targraph,vertex.label=V(targraph)$name,edge.width=weightnorm(E(targraph)$weight),layout=layout.fruchterman.reingold,vertex.label.cex=.7,edge.arrow.size=0.5)
}
)

setMethod("plot",signature(x="biCompare"),
	function(x,lod_cut=2,p_cut=0.05,plotUnchanged=F){
	weightnorm=function(m){
		if(length(m)==0) return(0)
		return((m-min(m))*2/(max(m)-min(m)))
	}
	result=x
	datanames01=rownames(result@state_1)
	datanames02=rownames(result@state_2)
	if(length(datanames01)==0) stop("No positive selection sites for both sequence set!")
	if(x@method=="ckaksCodon"|x@method=="ckaksAA"){
		filtertag01=((result@state_1>1)+(result@statistic_1>=lod_cut))!=2
		filtertag02=((result@state_2>1)+(result@statistic_2>=lod_cut))!=2
	}else{
		filtertag01=(result@statistic_1>p_cut)
		filtertag02=(result@statistic_2>p_cut)
	}
	resultmat=result@state_1
	resultmat02=result@state_2
	resultmat[filtertag01]=NA
	resultmat02[filtertag02]=NA
	resultmat[is.na(resultmat)]=0
	resultmat02[is.na(resultmat02)]=0
	if(plotUnchanged==F){
		unionsites=union(x@positiveSite01,x@positiveSite02)
		resultmat=resultmat[datanames01%in%unionsites,datanames01%in%unionsites]
		resultmat02=resultmat02[datanames02%in%unionsites,datanames02%in%unionsites]
	}
	namef_01=rownames(resultmat)
	namef_02=rownames(resultmat02)
	unionrnames=union(namef_01,namef_02)
	resultf=matrix(0,length(unionrnames),length(unionrnames),dimnames=list(unionrnames,unionrnames))
	resultf02=matrix(0,length(unionrnames),length(unionrnames),dimnames=list(unionrnames,unionrnames))
	for(i in unionrnames){
		for(j in unionrnames){
			if((i%in%namef_01)&(j%in%namef_01)){
				resultf[i,j]=resultmat[i,j]
			}
			if((i%in%namef_02)&(j%in%namef_02)){
				resultf02[i,j]=resultmat02[i,j]
			}
		}
	}
	if(x@method=="ckaksCodon"|x@method=="ckaksAA"){
		zz01=graph.adjacency(resultf,mode="directed",weighted=TRUE)
		zz02=graph.adjacency(resultf02,mode="directed",weighted=TRUE)
	}else{
		zz01=graph.adjacency(resultf,mode="lower",weighted=TRUE)
		zz02=graph.adjacency(resultf02,mode="lower",weighted=TRUE)
	}
	#lay=layout.svd(zz02)
	#lay=layout.kamada.kawai(zz02,kkconstSets=vcount(zz01)**10)
	lay=layout.fruchterman.reingold(zz02,repulserad=vcount(zz02)*1.5)
	par(mfrow=c(1,2),mai=c(0.05,0.05,0.05,0.05))
	nodecolvector=rep("grey",times=length(V(zz01)$name))
	nodecolvector[(V(zz01)$name)%in%x@positiveSite01]="skyblue1"
	nodecolvector[(V(zz01)$name)%in%x@positiveSite02]="lightpink1"
	plot(zz01,vertex.label=V(zz01)$name,layout=lay,edge.arrow.size=0.3,vertex.color=nodecolvector,edge.color="olivedrab4",edge.width=weightnorm(E(zz01)$weight))
	plot(zz02,vertex.label=V(zz02)$name,layout=lay,edge.arrow.size=0.3,vertex.color=nodecolvector,edge.color="olivedrab4",edge.width=weightnorm(E(zz02)$weight))
}
)

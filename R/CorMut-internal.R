#convert a vector of chars into a single string
#c2s=function (chars = c()) {
#	    return(paste(chars, collapse = ""))
#}
#convert a single string into a vector of chars
#s2c=function (string){
#	return(strsplit(string,split="")[[1]])
#}
#
translate_cm=seqinr::translate
#
fvft <-
function(x,refseq=refseq){
	refseq02=s2c(refseq)
	ft=0
	fv=0
	result=list()
	ag=c("a","g")
	tc=c("t","c")
	for(i in x){
		mono_seq=s2c(i)
		for(j in 1:length(mono_seq)){
			if(mono_seq[j]!=refseq02[j]&is.element(mono_seq[j],c("a","t","g","c"))){
				if(setequal(c(refseq02[j],mono_seq[j]),ag)|setequal(c(refseq02[j],mono_seq[j]),tc)){
					ft=ft+1
				}else{
					fv=fv+1
				}
			}
		}
	}
	L=length(mono_seq)
	sample=length(x)
	result$ft=ft/(L*sample)
	result$fv=fv/(2*L*sample)
	return(result)
}

Ntnv <-
function(x){
	x=toupper(x)
	result=list()
	TS=list()
	TS[["A"]]="G"
	TS[["T"]]="C"
	TS[["G"]]="A"
	TS[["C"]]="T"
	TV=list()
	TV[["A"]]=c("T","C")
	TV[["T"]]=c("A","G")
	TV[["G"]]=c("T","C")
	TV[["C"]]=c("A","G")
	x1=list()
	x2=list()
	for(i in 1:3 ){
		x0=x
		x0[i]=TS[[x[i]]]
		x1=append(x1,list(x0))
		x0[i]=TV[[x[i]]][1]
		x2=append(x2,list(x0))
		x0[i]=TV[[x[i]]][2]
		x2=append(x2,list(x0))
	}
	result$Nst=sum(sapply(x1,function(q)translate_cm(q)==translate_cm(x)))
	result$Nat=length(x1)-result$Nst
	result$Nsv=sum(sapply(x2,function(q)translate_cm(q)==translate_cm(x)))
	result$Nav=length(x2)-result$Nsv
	return(result)
}

Ntnv_aa <-
function(x,aa_m){
	x=toupper(x)
	result=list()
	TS=list()
	TS[["A"]]="G"
	TS[["T"]]="C"
	TS[["G"]]="A"
	TS[["C"]]="T"
	TV=list()
	TV[["A"]]=c("T","C")
	TV[["T"]]=c("A","G")
	TV[["G"]]=c("T","C")
	TV[["C"]]=c("A","G")
	x1=list()
	x2=list()
	for(i in 1:3){
		x0=x
		x0[i]=TS[[x[i]]]
		x1=append(x1,list(x0))
		x0[i]=TV[[x[i]]][1]
		x2=append(x2,list(x0))
		x0[i]=TV[[x[i]]][2]
		x2=append(x2,list(x0))
	}
	result$Nst=sum(sapply(x1,function(q)translate_cm(q)==translate_cm(x)))
	result$Nat=sum(sapply(x1,function(q)translate_cm(q)==aa_m))
	result$Nsv=sum(sapply(x2,function(q)translate_cm(q)==translate_cm(x)))
	result$Nav=sum(sapply(x2,function(q)translate_cm(q)==aa_m))
	return(result)
}

rawreplace <-
function(xx,refseq=""){
	seqf=lapply(xx,s2c)
	xb=seqf
	refseq=s2c(tolower(refseq))
	refaa=translate_cm(refseq)
	raw_set=c("r","w","v","b","m","s","h","n","y","k","d")
	Rawdict=list()
	Rawdict[["r"]]=c("a","g")
	Rawdict[["w"]]=c("a","t")
	Rawdict[["v"]]=c("a","g","c")
	Rawdict[["b"]]=c("g","c","t")
	Rawdict[["m"]]=c("a","c")
	Rawdict[["s"]]=c("g","c")
	Rawdict[["h"]]=c("a","c","t")
	Rawdict[["n"]]=c("a","g","c","t")
	Rawdict[["y"]]=c("c","t")
	Rawdict[["k"]]=c("g","t")
	Rawdict[["d"]]=c("a","g","t")
	#test

	for(i in 1:length(seqf)){
		for(j in 1:length(seqf[[i]])){
			if(j%%3==1){ tris=seqf[[i]][j:(j+2)]
			}else if(j%%3==2){ tris=seqf[[i]][(j-1):(j+1)]
			}else if(j%%3==0){ tris=seqf[[i]][(j-2):j]
			}
			#print(tris)
			if(seqf[[i]][j] %in% raw_set){
				count=0
				baked=0
				for(k in Rawdict[[seqf[[i]][j]]]){
					if(j%%3==1|j%%3==2){
						tris[j%%3]=k
					}
					if(j%%3==0){
						tris[3]=k
					}
					if(refaa[ceiling(j/3)]!=translate_cm(tris)){
						xb[[i]][j]=k
						count=count+1
					}
				}
				if(count==0){
					xb[[i]][j]=sample(setdiff(Rawdict[[seqf[[i]][j]]],refseq[j]),1)
				}
			}
		}
	}
	return(xb)
}

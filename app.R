
################################################
###web for Easykin
##by Ran Li, E mail:654163930@qq.com
##202100308


########################################################################################
#Import packages
library(shinydashboard)
library(shiny)
library(DT)
library(parallel)
library(Familias)
library('shinylogs')
#######################################################################
###functions for simulation

generate_allele<-function(freq,n=1)sample(names(freq),prob=freq,size=n,replace=TRUE)

#generate unrelated individuals		
generate_UN<-function(freq,n=1){
  genotype=list()
  alleles<-generate_allele(freq=freq,n=2*n)
  for (i in 1:n){
    genotype<-c(genotype,list(alleles[(2*i-1):(2*i)]))
  }
  genotype
}
##generate the genotypes of a child based on the genotypes of the father and mother
generate_child<-function(fa,mo,locus){
  m1<-locus$maleMutationMatrix
  m2<-locus$femaleMutationMatrix	
  allelec1<-mut2allele(MutationMatrix=m1,allele=sample(fa,1))	
  allelec2<-mut2allele(MutationMatrix=m2,allele=sample(mo,1))
  c(allelec1,allelec2)
}
mut2allele<-function(MutationMatrix,allele){
  prob<-MutationMatrix[,colnames(MutationMatrix)==allele]
  sample(x=rownames(MutationMatrix),size=1,prob=prob)
}



######
simulationPeds_genotype<-function(nonX=NULL,ped,loci,selected=NULL){
  ids<-ped$id	
  findex<-ped$findex	
  mindex<-ped$mindex
  findex[findex==0]<-NA
  mindex[mindex==0]<-NA
  UNid<-which(is.na(findex))
  pedmatrix<-data.frame(ids=ids,findex=findex,mindex=mindex)
  genotype_all<-matrix(NA,nrow=length(ids),ncol=2*length(loci),dimnames=list(ids))
  locusnames<-rep(NA,2*length(loci))
  for (i in 1:length(loci)){
    locus<-loci[[i]]
    locusname<-locus$locusname
    locusnames[c(2*i-1,2*i)]<-paste0(locusname,".",1:2)
    freq<-locus$alleles	
    #generate alleles of founders
    UN<-generate_UN(freq=freq,n=length(UNid))		
    for (j in 1:length(UNid))genotype_all[UNid[j],c(2*i-1,2*i)]<-UN[[j]]
  }
  colnames(genotype_all)<-locusnames
  genotype_all<- fill_Pedgenotype(genotype=genotype_all,pedmatrix=pedmatrix,locus=loci)
  if (is.null(selected))selected=1:nrow(genotype_all)
  as.data.frame(genotype_all[selected,])
}
fill_Pedgenotype<-function(genotype,pedmatrix,locus){
  #generate alleles of non-founders
  if (length(which(is.na(genotype[,1])))==0) return(genotype)
  idx<-which(is.na(genotype[,1]))
  for (i in idx){
    parents_id<-pedmatrix[i,2:3]
    fa<-genotype[parents_id[1,1],]
    mo<-genotype[parents_id[1,2],]
    if(NA %in% c(fa,mo))next
    for (j in seq(1,ncol(genotype),by=2)){
      fa1<-fa[c(j,j+1)]
      mo1<-mo[c(j,j+1)]
      locus0<-locus[[(j+1)/2]]
      genotype[i,c(j,j+1)]<-generate_child(fa=fa1,mo=mo1,locus=locus0)
    }
    
  }
  ##recurse until all members have a genotype
  genotype<-fill_Pedgenotype(genotype=genotype,pedmatrix=pedmatrix,locus=locus)
  
}

##functions for LR calculation
##pairwise relationships
cal_kinshipIndex2<-function(g,k1,k2,mut,freq,rare_freq=0.001){
  ##k1	null	      ibd=2	1	0
  ##k2	alternative	ibd=2	1	0
  g<-as.character(g)
  g1<-g[1:2]
  g2<-g[3:4]
  
  novel_allele<-g[!(g %in% names(freq))]													#rare allele
  if (length(novel_allele)>0){
    novel_freq<-rep(rare_freq,length(novel_allele))
    names(novel_freq)<-novel_allele
    freq<-c(freq,novel_freq)}
  
  
  shared_allele<-as.character(g1[g1 %in% g2])
  n1<-length(shared_allele)
  n2<-length(g2[g2 %in% g1])
  n3<-max(g1[1]==g1[2],g2[1]==g2[2])
  unique_alleles<-as.character(unique(g))
  f_unique<-freq[names(freq) %in% unique_alleles]
  f_share<-freq[names(freq) %in% shared_allele]
  
  if (max(n1,n2)==0) {																	##AA BB		AB CD		AA BC
    p1<-cal_k0_likelihood(k0=k1[3],g=g,mut=mut)
    p2<-cal_k0_likelihood(k0=k2[3],g=g,mut=mut)
    
  }					
  if (min(n1,n2)==1 & n3==0) {															##AB AC
    p1<-k1[2]*prod(f_unique)+4*k1[3]*prod(f_unique)*f_share
    p2<-k2[2]*prod(f_unique)+4*k2[3]*prod(f_unique)*f_share
  }
  if (min(n1,n2)==1 & n3==1) {															##AA AB
    p1<-k1[2]*prod(f_unique)*f_share+2*k1[3]*prod(f_unique)*f_share^2
    p2<-k2[2]*prod(f_unique)*f_share+2*k2[3]*prod(f_unique)*f_share^2
  }
  if (min(n1,n2)==2 & n3==0) {															##AB AB
    p1<-2*k1[1]*prod(f_unique)+k1[2]*prod(f_unique)*sum(f_unique)+4*k1[3]*prod(f_unique)^2
    p2<-2*k2[1]*prod(f_unique)+k2[2]*prod(f_unique)*sum(f_unique)+4*k2[3]*prod(f_unique)^2
  }
  if (min(n1,n2)==2 & n3==1) 	{															##AA AA
    p1<-k1[1]*f_unique^2+k1[2]*f_unique^3+k1[3]*f_unique^4
    p2<-k2[1]*f_unique^2+k2[2]*f_unique^3+k2[3]*f_unique^4
  }
  return(p1/p2)
}

cal_k0_likelihood<-function(k0,g,mut){
  if (k0!=0) return(k0)
  g<-tryCatch({as.numeric(g)},warning="none_numeric")
  if ("none_numeric" %in% g) return(10^-8)
  gF<-g[1:2];gC<-g[3:4]
  mutstep<-c(abs(gF[1]-gC[1]),abs(gF[1]-gC[2]),abs(gF[2]-gC[1]),abs(gF[2]-gC[2]))
  minstep<-which(mutstep==min(mutstep))[1]
  pid<-ifelse(minstep<=2,1,2)
  cid<-ifelse(minstep %% 2==1,1,2)
  x1=ifelse(gF[1]==gF[2],1,0.5)
  prob<-mut[rownames(mut)==as.character(gC[cid]),rownames(mut)==as.character(gF[pid])]*x1
  return(prob)
}

###
cal_pairsLR<-function(genotype,k,loci){
  LRs<-c()
  locus<-c()
  for (i in 1:length(loci)){
    locusname<-loci[[i]]$locusname
    freq<-loci[[i]]$alleles
    mut<-loci[[i]]$maleMutationMatrix
    idx<-grep(colnames(genotype),pattern=locusname,value=FALSE)
    if (length(idx)==0)next
    g<-genotype[,idx]		
    g<-t(g)
    g<-c(g[,1],g[,2])
    if (k[[1]][3] == 1)lr <- 1 else lr<-cal_kinshipIndex2(g,k1=k[[1]],k2=k[[2]],mut=mut,freq=freq)
    LRs<-c(LRs,lr)
    locus<-c(locus,locusname)
  }
  names(LRs)<-locus
  return(LRs)
}
#a

#####################################
##cal LR for a pedigree
##one references using ITO method
##otherwise using Familias
cal_pedsLR<-function(peds,loci,refID=2,genotype,theta=0,POI="POI",mut=0.002){
  
  ped1=peds[[1]]
  ped2=peds[[2]]
  affected1<-which(peds[[3]]==1)
  affected2<-which(peds[[4]]==1)
  if (length(affected1)== -2){						#pairwise relationship, not applicable at present
    s<-ped1$id[affected1]
    kvalues<-kPairs(peds[1:2],s)
    LR<-cal_pairsLR(genotype=genotype,k=kvalues,loci=loci)
  }else{
    LR<-FamiliasPosterior(pedigrees=peds[1:2],loci=loci,datamatrix=genotype,ref=refID,kinship=theta)
    LR<-LR$LRperMarker[,1]		
    
  }
  c(CLR=prod(LR),LRs=list(LR))
}

FamiliasProb<-function(peds,loci,genotypes,refID=2,theta=0){
  LR1<-FamiliasPosterior(pedigrees=peds,loci=loci,datamatrix=genotypes,ref=refID,kinship=theta)
  sum(log10(LR1$LRperMarker[,1]))
}

simulation_pedsLR_batch<-function(peds,loci,refID=2,n=10,theta=0){
  affected1<-which(peds[[3]]==1)
  affected2<-which(peds[[4]]==1)
  if (length(affected1)==1)return()
  if (length(affected1)==2){
    s<-peds[[1]]$id[affected1]
    kvalues<-kPairs(peds[1:2],s)
  }
  #detect the cores
  ncore<-detectCores(logical = F)  														# cores
  # multicore processing, only recommended when large simulations are to be processed, say over 5000.
  ##if so,  just change the condition from ncore >3000 to ncore > 1
  if (ncore > 3000){
    #ncore <- ncore -2 
    #incProgress(1/100,detail = paste0(" Multithreads processing"))
    cl <- makeCluster(getOption("cl.cores", ncore));									
    clusterExport(cl = cl, varlist=c("generate_allele","generate_UN","generate_child","mut2allele",
                                     "kPairs","relationship_distance","pedigree_finder","parent_finder","cal_k0_likelihood",
                                     "cal_kinshipIndex2","cal_pairsLR","simulationPeds_genotype","fill_Pedgenotype","cal_pedsLR"),envir=environment())
    clusterEvalQ(cl,{library(Familias)})
    
    ped1_g=parLapply(cl,X=1:n,fun=simulationPeds_genotype,ped=peds[[1]],loci=loci,selected=affected1)
    ped2_g=parLapply(cl,X=1:n,fun=simulationPeds_genotype,ped=peds[[2]],loci=loci,selected=affected2)
    if (length(affected1)==2){
      LR1<-log10(parSapply(cl=cl,X=ped1_g,FUN=cal_pairsLR,k=kvalues,loci=loci))
      LR2<-log10(parSapply(cl=cl,X=ped2_g,FUN=cal_pairsLR,k=kvalues,loci=loci))			
      LRs<-cbind(colSums(LR1),colSums(LR2))				
    }else{
      LR1<-parSapply(cl=cl,X=ped1_g,FUN=FamiliasProb,peds=peds[1:2],loci=loci,refID=refID,theta=theta)
      LR2<-parSapply(cl=cl,X=ped2_g,FUN=FamiliasProb,peds=peds[1:2],loci=loci,refID=refID,theta=theta)				
      LRs<-cbind(LR1,LR2)
    }			
    stopCluster(cl)
    
    # single core processing	
  }else{
    t1 <- Sys.time()
    withProgress(message = 'Calculation in progress',
                 detail = 'This may take a while...', value = 0, {
                   LRs<-c()
                   for (i in 1:n){
                     t2 <- Sys.time()
                     incProgress(1/n,detail = paste0(" time remaining ",ceiling((n-i)*as.numeric(difftime(t2,t1,units = "secs"))/i)," sec"))
                     ped1_g=simulationPeds_genotype(1,ped=peds[[1]],loci=loci,selected=affected1)
                     ped2_g=simulationPeds_genotype(1,ped=peds[[2]],loci=loci,selected=affected2)
                     LR1<-cal_pedsLR(peds,loci,refID=refID,genotype=ped1_g,theta=theta,POI="POI",mut=0.002)
                     LR2<-cal_pedsLR(peds,loci,refID=refID,genotype=ped2_g,theta=theta,POI="POI",mut=0.002)
                     LRs<-rbind(LRs,c(sum(log10(LR1$CLR)),sum(log10(LR2$CLR))))
                     
                   }
                 })
  }
  colnames(LRs)<-c("H1 true","H2 true")
  
  LRs
  
}

prune_pedigree<-function(peds,loci,available=NULL,POI=NULL,refID=2,n=10,theta=0,tx=0){
  #simulation_pedsLR_batch(peds,affecteds,loci,refID=2,n=10,theta=0)
  ped1<-peds[[1]]
  ped2<-peds[[2]]
  share_ids<-ped1$id[ped1$id %in% ped2$id]
  if (is.null(available))available<-share_ids
  #the person of interest, POI
  if (is.null(POI))POI=share_ids[1]
  references_ids<-available[available!=POI]
  refs_POI_combs<-list()
  refs_POI_groups<-c()
  for (i in length(references_ids):1){
    combs<-combn(references_ids,i)
    for (j in 1:ncol(combs)){
      x<-c(POI,combs[,j])
      refs_POI_combs<-c(refs_POI_combs,list(x))
      refs_POI_groups<-c(refs_POI_groups,paste(x,collapse="+"))
    }
  }
  
  #simulation and calculate LRs	
  combs_LRs<-vector("list", length(refs_POI_combs))
  t1 <- Sys.time()
  withProgress(message = 'Calculation in progress',
               detail = 'This may take a while...', value = 0, {
                 for (idx in 1:n){
                   t2 <- Sys.time()
                   incProgress(1/n,detail = paste0(" time remaining ",ceiling((n-idx)*as.numeric(difftime(t2,t1,units = "secs"))/idx)," sec"))
                   P1genotypes<-simulationPeds_genotype(ped=ped1,loci=loci)
                   P2genotypes<-simulationPeds_genotype(ped=ped2,loci=loci)
                   for (i in 1:length(refs_POI_combs)){
                     ids<-refs_POI_combs[[i]]
                     subg1<-P1genotypes[rownames(P1genotypes) %in% ids,]
                     subg2<-P2genotypes[rownames(P2genotypes) %in% ids,]
                     peds[[3]]<-as.numeric(peds[[1]]$id %in% ids)
                     peds[[4]]<-as.numeric(peds[[2]]$id %in% ids)
                     
                     LR1 <- tryCatch(
                       {cal_pedsLR(peds,loci,refID=refID,genotype=subg1,theta=theta,POI=POI,mut=0.002)},
                       warning = function(w) {"error"}, 
                       error = function(e) {"error"}, 
                       finally = {"error"})
                     if (class(LR1)=="character")	{combs_LRs[[i]]<-NA;next}	
                     LR2<-cal_pedsLR(peds,loci,refID=refID,genotype=subg2,theta=theta,POI="POI",mut=0.002)
                     LRs<-c("H1 true"=log10(LR1$CLR),"H2 true"=log10(LR2$CLR))
                     combs_LRs[[i]]<-rbind(combs_LRs[[i]],LRs)
                     
                   }
                 }
               })
  
  names(combs_LRs)<-refs_POI_groups
  return(combs_LRs)
  #estimate the system power for each combination
  
  
}
#estimate the Cotterman coefficients for two individuals
kPairs<-function(peds,s){
  kvalue<-list()
  for (i in 1:length(peds)){
    ped<-peds[[i]]
    ids<-ped$id
    sex<-ped$sex
    fa<-ped$findex;fa[fa==0]<-NA;fa<-ids[fa]		
    mo<-ped$mindex;mo[mo==0]<-NA;mo<-ids[mo]		
    pedigreedata<-data.frame(ids=ids,sex=sex,fa=fa,mo=mo)
    d<- relationship_distance(pedigreedata=pedigreedata,s)
    if (max(d)==0)k<-c(0,0,1)
    if (d[1]==1 & d[2]==1)k<-c(0,1,0)
    if (d[1]==2 & d[2]==1)k<-c(0.25,0.5,0.25)
    if (d[2]>1){k<-c(0,1/2^(d[2]-1),1-1/2^(d[2]-1))}
    names(k)<-c("k2","k1","k0")
    kvalue<-c(kvalue,list(k))
  }
  kvalue
}

#calculate the number of degrees and miosis between two individuals
relationship_distance<-function(pedigreedata,s){
  s1<-s[1]
  s2<-s[2]
  result<-c("miosisN"=0,"degreeN"=0)
  pedigree_s1<-pedigree_finder(pedigreedata,s1)
  pedigree_s2<-pedigree_finder(pedigreedata,s2)
  
  findIDX<-FALSE
  
  if (s2 %in% unlist(pedigree_s1)){						##s2 is an ancestor of s1
    for (i in 1:length(pedigree_s1)){
      parentnames<-unlist(pedigree_s1[i])
      if (s2 %in% parentnames) {result[1]=i;result[2]=i;return(result)}
    }
  }
  
  if (s1 %in% unlist(pedigree_s2)){						###s1 is an ancestor of s2
    for (i in 1:length(pedigree_s2)){
      parentnames<-unlist(pedigree_s2[i])
      if (s1 %in% parentnames) {result[1]=i;result[2]=i;return(result)}
    }
  }
  
  if (length(pedigree_s1)==0 | length(pedigree_s2)==0) return(result)
  
  for (i in 1:length(pedigree_s1)){						###collateral
    for (j in 1:length(pedigree_s2)){
      parentnames1<-pedigree_s1[[i]]
      parentnames2<-pedigree_s2[[j]]
      coreperson<-parentnames1[parentnames1 %in% parentnames2]
      t<-sum(as.numeric(parentnames1 %in% parentnames2),na.rm=TRUE)
      if (t==1){result[1]<-i+j;result[2]<-i+j;return(result)}				#share one ancestor
      if (t==2){result[1]<-i+j;result[2]<-i+j-1;return(result)}			#share two ancestor
    }	
  }	
  return(result)
  #return(c(distance=list(result),pedigree_s1=list(pedigree_s1),pedigree_s2=list(pedigree_s2)))
}


###
pedigree_finder<-function(pedigreedata,s){
  pedigree_s<-list()
  n=0
  new_s<-s
  parentnames<-c(NA)
  while(length(parentnames)!=0){						#iteratively finder parents
    if (n>50) break		
    parentnames<-c()
    for (s1 in new_s){
      x<-as.character(parent_finder(pedigreedata,s1))
      parentnames<-c(parentnames,x[!is.na(x)])
    }	
    if (length(parentnames)!=0){pedigree_s<-c(pedigree_s,list(parentnames));n=n+1}
    new_s<-parentnames	
  }
  if (n>0)names(pedigree_s)<-paste("generation",1:n,sep="_")
  return(pedigree_s)
}

###
parent_finder<-function(pedigreedata,s){
  fam_members<-pedigreedata[,1]
  parentnames<-unlist(pedigreedata[fam_members==s,3:4])
  return(parentnames)
}
#g
transform_genotype<-function(genotype,loci,rare_freq=0.001){
  new_genotype<-c()
  new_loci<-list()
  locusnames<-c()
  genotyped_loci<-colnames(genotype)
  ref_loci<-names(loci)
  for (i in 1:length(genotyped_loci)){
    locus<-genotyped_loci[i]	
    #only include loci with frequency data available
    if (locus %in% ref_loci){
      #locus<-loci[[i]]$locusname
      locusnames<-c(locusnames,locus)
      idx<-which(ref_loci==locus)
      ref_freq<-loci[[idx]]$alleles
      ref_alleles<-names(ref_freq)
      g<-genotype[,i]
      g<-strsplit(g,split=",")
      alleles<-unlist(g)
      g<-do.call(rbind,g)
      colnames(g)<-paste0(locus,".",1:2)
      new_genotype<-cbind(new_genotype,g)
      x<-loci[[idx]]
      novel_alleles<-unique(alleles[!(alleles %in% ref_alleles)])
      if (length(novel_alleles)>0){
        new_freq<-c(ref_freq,rep(rare_freq,length(novel_alleles)))
        new_freq<-new_freq/sum(new_freq)
        mut_male<-1-x$maleMutationMatrix[1,1]
        mut_female<-1-x$femaleMutationMatrix[1,1]
        x<-FamiliasLocus(frequencies=new_freq,allelenames=c(ref_alleles,novel_alleles),name=locus,MutationModel="Equal",
                         maleMutationRate=mut_male,femaleMutationRate=mut_female)
        
      }
      new_loci<-c(new_loci,list(x))
    }
    
  }
  
  new_genotype<-as.data.frame(new_genotype)
  names(new_loci)<-locusnames
  return(c(genotype=list(new_genotype),loci=list(new_loci)))
  
}
#########################################################
##system effectiveness
##plot multi Histogram
multiHist<-function(x=as.data.frame(),thresholds=c(0,0),freq=FALSE,breakNum=30,xlab="value",ylab="Density",main="distribution",legend=NULL,xlim=c(),ylim=c(),alpha=0.3,col_lines=NULL,border=NULL,barcol=NULL){
  ObjectName<-colnames(x)
  if (length(legend)==0) legend=ObjectName
  if (length(barcol)==0)cols<-rainbow(ncol(x),alpha=alpha) else cols<-barcol
  if (length(col_lines)==0)col_lines<-rainbow(ncol(x),alpha=0.9)
  
  minValue<-min(x,na.rm=TRUE)
  maxValue<-max(x,na.rm=TRUE)
  
  width<-(maxValue-minValue)/breakNum   										#the distance between each colum
  if (length(xlim)==0) xlim=c(minValue-2*width,maxValue+2*width)
  BreakArray<-c(minValue,1:(breakNum+1)*width+minValue)
  group_density<-c()
  for ( i in 1:ncol(x)){
    histResult<-hist(x[,i],breaks= BreakArray,plot=FALSE)
    group_density<-c(group_density,histResult$density)
  }
  if (length(ylim)==0) ylim<-c(0,max(group_density)*1.2)
  for ( i in 1:ncol(x)){
    
    if (i==1) {
      histResult<-hist(x[,i],freq=freq,breaks= BreakArray,plot=T,col=cols[i],xlab=xlab,ylab=ylab,main="",border=border,xlim=xlim,ylim=ylim)
      
      xseq<-seq(from=minValue,to=maxValue,by=(maxValue-minValue)/1000)
      if (length(xlim)!=0)xseq<-seq(from=min(xlim),to=max(xlim),by=(max(xlim)-min(xlim))/1000)
      yseq<-dnorm(xseq,mean(x[,i]),sd(x[,i]))
      lines(xseq,yseq,col=col_lines[i],lwd=2)
      legend(x=xlim[1],y=ylim[2],legend=strsplit(main,split = ",")[[1]],col=cols,pch=NA,xpd=TRUE,bg="gray88",box.col = "gray88",y.intersp=1.3)
      legend(x=xlim[2],y=ylim[2],legend=legend,col=cols,pch=15,xpd=TRUE,pt.cex=2)
    } else{
      histResult<-hist(x[,i],freq=freq,breaks=BreakArray,add=T,xlab="",ylab="",col=cols[i],border=border)
      
      xseq<-seq(from=minValue,to=maxValue,by=(maxValue-minValue)/1000)
      if (length(xlim)!=0)xseq<-seq(from=min(xlim),to=max(xlim),by=(max(xlim)-min(xlim))/1000)
      yseq<-dnorm(xseq,mean(x[,i]),sd(x[,i]))
      lines(xseq,yseq,col=col_lines[i],lwd=2)
    }
  }
  ###thresholds
  if (thresholds[1]!=thresholds[2]){
    lines(x=c(thresholds[1],thresholds[1]),y=c(0,ylim[2]),col="black",lwd=2,lty=2)
    text(x=thresholds[1],y=ylim[2]*0.99,labels="t1",pos=2,cex=1.3)
    lines(x=c(thresholds[2],thresholds[2]),y=c(0,ylim[2]),col="black",lwd=2,lty=2)
    text(x=thresholds[2],y=ylim[2]*0.99,labels="t2",pos=4,cex=1.3)
  }else{
    lines(x=c(thresholds[1],thresholds[1]),y=c(0,ylim[2]),col="black",lwd=2,lty=2)
    text(x=thresholds[1],y=ylim[2]*0.99,labels="t",pos=2,cex=1.3)
    
  }
}


#confusion_matrix and system power
cal_non_overlaped<-function(u1,u2,s1,s2){
  if (u1>u2) {x<-u2;u2<-u1;u1<-x;x<-s2;s2<-s1;s1<-x}	#u1 is samller
  ux1=u1
  ux2=u2
  n=0
  y=10
  while (abs(y)>0.001){
    n=n+1
    if (n>1000) {print("time out");break}
    x<-(ux1+ux2)/2
    y<-(x-u1)^2/(2*s1^2)-(x-u2)^2/(2*s2^2)+log(s1)-log(s2)
    if (y<0)ux1<-x else ux2<-x
  }
  ##area
  prob1<-1-pnorm(x,mean=u1,sd=s1)
  prob2<-pnorm(x,mean=u2,sd=s2)
  probs<-prob1+prob2
  result<-c(x,probs)
  names(result)<-c("x","probs")
  return(result)
}

#confusion matrix and system power using observations
syspower<-function(paradata=as.matrix(),t=0){
  t1<-min(t)
  t2<-max(t)
  para1<-as.numeric(as.character(paradata[,1]))
  para2<-as.numeric(as.character(paradata[,2]))
  if (mean(para1,na.rm=TRUE)>mean(para2,na.rm=TRUE)){para_x<-para2;para2<-para1;para1<-para_x} ##para1 is samller
  n_para1<-length(para1)
  n_para2<-length(para2)
  confusion_m<-matrix(NA,2,3,dimnames=list(c("neg_group","pos_group"),c("neg_group_predict","grey_zone","pos_group_predict")))
  confusion_m[1,1]<-length(which(para1<t1))
  confusion_m[1,3]<-length(which(para1>t2))
  confusion_m[1,2]<-n_para1-sum(confusion_m[1,c(1,3)],na.rm=TRUE)
  confusion_m[2,1]<-length(which(para2<t1))
  confusion_m[2,3]<-length(which(para2>t2))
  confusion_m[2,2]<-n_para2-sum(confusion_m[2,c(1,3)],na.rm=TRUE)
  sen<-confusion_m[2,3]/n_para2
  spe<-confusion_m[1,1]/n_para1
  ppv<-confusion_m[2,3]/sum(confusion_m[,3],na.rm=TRUE)
  npv<-confusion_m[1,1]/sum(confusion_m[,1],na.rm=TRUE)
  gray_rate<-sum(confusion_m[,2])/sum(n_para1+n_para2)
  sys<-sum(confusion_m[,c(1,3)],na.rm=TRUE)/(n_para1+n_para2)
  error_rate<-(confusion_m[2,1]+confusion_m[1,3])/sum(confusion_m[,c(1,3)],na.rm=TRUE)	#average error
  overlap_area<-cal_non_overlaped(u1=mean(para1),u2=mean(para2),s1=sd(para1),s2=sd(para2))
  
  sys_para<-c(sen,spe,ppv,npv,error_rate,sys)
  names(sys_para)<-c("Sen","Spe","PPV","NPV","Error rate","Effectiveness")
  result<-list(confusion_m,sys_para)
  names(result)<-c("confusion_m","sys_parameters")
  return(result)
}




#confusion matrix and system power based normal fit
syspower2<-function(paradata=as.matrix(),t=0,W.ratio=NULL){
  t1<-min(t)
  t2<-max(t)
  para1<-as.numeric(as.character(paradata[,1]))
  para2<-as.numeric(as.character(paradata[,2]))
  if (mean(para1,na.rm=TRUE)>mean(para2,na.rm=TRUE)){para_x<-para2;para2<-para1;para1<-para_x} ##para1 is smaller
  u1<-mean(para1,na.rm=TRUE)
  u2<-mean(para2,na.rm=TRUE)
  sd1<-sd(para1,na.rm=TRUE)
  sd2<-sd(para2,na.rm=TRUE)
  sen<-1-pnorm(mean=u2,sd=sd2,q=t2)
  gray2<-pnorm(mean=u2,sd=sd2,q=t2)-pnorm(mean=u2,sd=sd2,q=t1)
  FN<-pnorm(mean=u2,sd=sd2,q=t1)
  
  spe<-pnorm(mean=u1,sd=sd1,q=t1)
  gray1<-pnorm(mean=u1,sd=sd1,q=t2)-pnorm(mean=u1,sd=sd1,q=t1)
  FP<-1-pnorm(mean=u1,sd=sd1,q=t2)
  
  ppv<-sen/(1-pnorm(mean=u1,sd=sd1,q=t2)+sen)
  npv<-spe/(pnorm(mean=u2,sd=sd2,q=t1)+spe)
  gray_rate<-(gray1+gray2)/2
  overlap_area<-cal_non_overlaped(u1=u1,u2=u2,s1=sd1,s2=sd2)
  sys<-(sen+spe)/2
  sys<-1-gray_rate
  
  error_rate<-(1-pnorm(mean=u1,sd=sd1,q=t2)+pnorm(mean=u2,sd=sd2,q=t1))/2
  confusion_m<-matrix(NA,2,3,dimnames=list(c("neg_group","pos_group"),c("neg_group_predict","grey_zone","pos_group_predict")))
  confusion_m[1,1]<-spe
  confusion_m[1,2]<-gray1
  confusion_m[1,3]<-FP
  confusion_m[2,1]<-FN
  confusion_m[2,2]<-gray2
  confusion_m[2,3]<-sen
  sys_para<-c(sen,spe,ppv,npv,FP,FN,gray_rate,sys)
  names(sys_para)<-c("Sen","Spe","PPV","NPV","FPR","FNR","Inconclusive","Effectiveness")
  result<-list(confusion_m,sys_para)
  names(result)<-c("confusion_m","sys_parameters")
  ##estimate decision threshold according Marsico FL. Forensic Sci Int Genet. 2021
  if (!is.null(W.ratio)){
    WED0=1000
    for (tx in seq(min(paradata),max(paradata),0.01)){
      FN<-pnorm(mean=u2,sd=sd2,q=tx)
      FP<-1-pnorm(mean=u1,sd=sd1,q=tx)
      WED <- sqrt((W.ratio*FP)^2+FN^2)
      if (WED < WED0){
        DT=tx
        FN.out <- FN
        FP.out <- FP
        WED0 <- WED
      }
    }
    result$DT=DT
    result$error <- c("FP"=round(FP.out,8),"FN"=round(FN.out,8))
  }
  return(result)
}
syspower_batch<-function(x,t=NULL,fit=FALSE){
  if (length(t)==2) t<-matrix(t,nrow=1,ncol=2)
  t_t<-paste(t[,1],t[,2],sep=" ~ ")
  result<-c()
  
  for (i in 1:nrow(t)){
    if (fit) y<-syspower2(x,t=t[i,])	
    if (!fit) y<-syspower(x,t=t[i,])
    result<-rbind(result,y$sys_parameters)
  }
  rownames(result)<-t_t
  return(result)
}


list_freqdata<-function(x){
  loci<-as.character(unique(x[,1]))
  result<-list()
  for (locus in loci){
    subx<-x[x[,1]==locus,]
    alleles<-as.character(subx[,2])
    freqs<-as.numeric(as.character(subx[,3]))
    unique_alleles<-unique(alleles)
    unique_freqs<-c()
    for (a in unique_alleles){
      idx<-which(alleles==a)
      unique_freqs<-c(unique_freqs,freqs[idx[1]])
    }
    
    unique_freqs<-unique_freqs/sum(unique_freqs)
    names(unique_freqs)<-unique_alleles
    result<-c(result,list(unique_freqs))
  }
  names(result)<-loci
  result
}
##
#############################################################################################################################################
##pedigrees and Hypothesis

#############################################################################################################################################


inputPedigree<-function(peddata){
  subdata<-peddata
  idx <- subdata[,1] == subdata[1,1]
  ids<-subdata[idx,2]
  faid<-subdata[idx,3]
  moid<-subdata[idx,4]
  sex<-subdata[idx,5]
  sex <- sapply(sex, FUN=function(x){ifelse(x==1,"male","female")})
  h1_affected<-subdata[idx,6]
  H1<-FamiliasPedigree(id=ids, dadid=faid, momid=moid, sex=sex)
  
  idx <- !idx
  ids<-subdata[idx,2]
  faid<-subdata[idx,3]
  moid<-subdata[idx,4]
  sex<-subdata[idx,5]
  sex <- sapply(sex, FUN=function(x){ifelse(x==1,"male","female")})
  h2_affected<-subdata[idx,6]
  H2<-FamiliasPedigree(id=ids, dadid=faid, momid=moid, sex=sex)
  
  peds<-c(H1=list(H1),
          H2=list(H2),
          h1_affected=list(h1_affected),
          h2_affected=list(h2_affected))
  return(peds)
}

#poiSex <- "female"
initialize_ped0 <- function(poiSex){
  matrix(ncol=5,byrow=TRUE,data=c(
					"POI","fa","mo",poiSex,0,
					"fs_b1","fa","mo","male",0,
					"fs_b2","fa","mo","male",0,
					"fs_b3","fa","mo","male",0,
					"fs_s1","fa","mo","female",0,
					"fs_s2","fa","mo","female",0,
					"fs_s3","fa","mo","female",0,
					"fa","pgf","pgm","male",0,
					"mo","mgf","mgm","female",0,
					"pgf",NA,NA,"male",0,
					"pgm",NA,NA,"female",0,
					"mgf",NA,NA,"male",0,
					"mgm",NA,NA,"female",0,
					"puncle1","pgf","pgm","male",0,
					"puncle2","pgf","pgm","male",0,
					"puncle3","pgf","pgm","male",0,
					"paunt1","pgf","pgm","female",0,
					"paunt2","pgf","pgm","female",0,
					"paunt3","pgf","pgm","female",0,
					"muncle1","mgf","mgm","male",0,
					"muncle2","mgf","mgm","male",0,
					"muncle3","mgf","mgm","male",0,
					"maunt1","mgf","mgm","female",0,
					"maunt2","mgf","mgm","female",0,
					"maunt3","mgf","mgm","female",0,
					"hs_mo",NA,NA,"female",0,
					"hs_fa",NA,NA,"male",0,
					"phs_b1","fa","hs_mo","male",0,
					"phs_b2","fa","hs_mo","male",0,
					"phs_b3","fa","hs_mo","male",0,
					"phs_s1","fa","hs_mo","female",0,
					"phs_s2","fa","hs_mo","female",0,
					"phs_s3","fa","hs_mo","female",0,
					"mhs_b1","hs_fa","mo","male",0,
					"mhs_b2","hs_fa","mo","male",0,
					"mhs_b3","hs_fa","mo","male",0,
					"mhs_s1","hs_fa","mo","female",0,
					"mhs_s2","hs_fa","mo","female",0,
					"mhs_s3","hs_fa","mo","female",0,
					"sis_inlaw",NA,NA,"female",0,
					"nephew_b1","fs_b1","sis_inlaw","male",0,
					"nephew_b2","fs_b1","sis_inlaw","male",0,
					"nephew_b3","fs_b1","sis_inlaw","male",0,
					"niece_b1","fs_b1","sis_inlaw","female",0,
					"niece_b2","fs_b1","sis_inlaw","female",0,
					"niece_b3","fs_b1","sis_inlaw","female",0,
					"bro_inlaw",NA,NA,"male",0,
					"nephew_s1","bro_inlaw","fs_s1","male",0,
					"nephew_s2","bro_inlaw","fs_s1","male",0,
					"nephew_s3","bro_inlaw","fs_s1","male",0,
					"niece_s1","bro_inlaw","fs_s1","female",0,
					"niece_s2","bro_inlaw","fs_s1","female",0,
					"niece_s3","bro_inlaw","fs_s1","female",0,
					"spouse",NA,NA,ifelse(poiSex=="male","female","male"),0,
					"son1",ifelse(poiSex=="male","POI","spouse"),ifelse(poiSex=="male","spouse","POI"),"male",0,
					"son2",ifelse(poiSex=="male","POI","spouse"),ifelse(poiSex=="male","spouse","POI"),"male",0,
					"son3",ifelse(poiSex=="male","POI","spouse"),ifelse(poiSex=="male","spouse","POI"),"male",0,
					"daughter1",ifelse(poiSex=="male","POI","spouse"),ifelse(poiSex=="male","spouse","POI"),"female",0,
					"daughter2",ifelse(poiSex=="male","POI","spouse"),ifelse(poiSex=="male","spouse","POI"),"female",0,
					"daughter3",ifelse(poiSex=="male","POI","spouse"),ifelse(poiSex=="male","spouse","POI"),"female",0,
					"daughter_inlaw",NA,NA,"female",0,
					"gson_s1","son1","daughter_inlaw","male",0,
					"gson_s2","son1","daughter_inlaw","male",0,
					"gson_s3","son1","daughter_inlaw","male",0,
					"gdaug_s1","son1","daughter_inlaw","female",0,
					"gdaug_s2","son1","daughter_inlaw","female",0,
					"gdaug_s3","son1","daughter_inlaw","female",0,
					"son_inlaw",NA,NA,"male",0,
					"gson_d1","son_inlaw","daughter1","male",0,
					"gson_d2","son_inlaw","daughter1","male",0,
					"gson_d3","son_inlaw","daughter1","male",0,
					"gdaug_d1","son_inlaw","daughter1","female",0,
					"gdaug_d2","son_inlaw","daughter1","female",0,
					"gdaug_d3","son_inlaw","daughter1","female",0
					
					))
}
##delete undetected individuals
clean_ped0 <- function(peddata,available=c("POI","fs_b1"),plot=FALSE){
	ids <- peddata[,1]
	idx <- ids %in% available
	peddata[idx,5] <- 1
	unav_ids <- ids[which(!idx)]
	unav_ids <- unav_ids[!( unav_ids %in% peddata[,2:3])]	
	new_peddata <- peddata[!(ids %in% unav_ids),]
	if (nrow(new_peddata)< nrow(peddata))new_peddata <- clean_ped0(peddata=new_peddata,available=available)
	
	##plot
	if (plot){
		x=FamiliasPedigree(id=new_peddata[,1], dadid=new_peddata[,2], momid=new_peddata[,3], sex=new_peddata[,4])
		plot(x,affected=as.numeric(new_peddata[,5]))
	}
	return(new_peddata)
}


clean_ped <- function(peddata,available=c("POI","fs_b1"),plot=FALSE){
	ids <- peddata[,1]
	idx <- ids %in% available
	peddata[idx,5] <- 1
	unav_ids <- ids[which(!idx)]
	unav_ids <- unav_ids[!( unav_ids %in% peddata[,2:3])]	
	new_peddata <- peddata[!(ids %in% unav_ids),]
	if (nrow(new_peddata)< nrow(peddata))new_peddata <- clean_ped0(peddata=new_peddata,available=available)
	
	###delete grandparents if necessary
	available_state <- as.list(as.numeric(new_peddata[,5]))
	names(available_state) <- new_peddata[,1]
	if (sum(unlist(available_state[c("pgf","pgm")]))==0){											#if both pgf and pgm nor available
		if (sum(c("puncle1","puncle2","puncle3","paunt1","paunt2","paunt3") %in% available)==0){	#uncle and aunt are all not available
			new_peddata[new_peddata[,1]=="fa",2:3] <- NA
			new_peddata<-new_peddata[-which(new_peddata[,1] %in% c("pgf","pgm")),]					#delete grandparents
		}
	}
	if (sum(unlist(available_state[c("mgf","mgm")]))==0){											#if both mgf and mgm nor available
		if (sum(c("muncle1","muncle2","muncle3","maunt1","maunt2","maunt3") %in% available)==0){	#uncle and aunt are all not available
			new_peddata[new_peddata[,1]=="mo",2:3] <- NA
			new_peddata<-new_peddata[-which(new_peddata[,1] %in% c("mgf","mgm")),]					#delete grandparents
		}
	}
	##
	###delete fa and mo if necessary
	if (sum(is.na(c(new_peddata[new_peddata[,1]=="mo",2:3],new_peddata[new_peddata[,1]=="fa",2:3])))==4){
		if (sum(c("fs_b1","fs_b2","fs_b3","fs_s1","fs_s2","fs_s3",
				"mhs_b1","mhs_b2","mhs_b3","mhs_s1","mhs_s2","mhs_s3",
				"phs_b1","phs_b2","phs_b3","phs_s1","phs_s2","phs_s3","fa","mo",
				"nephew_b1","nephew_b2","nephew_b3","niece_b1","niece_b2","niece_b3",
				"nephew_s1","nephew_s2","nephew_s3","niece_s1","niece_s2","niece_s3"
		) %in% available)==0){
			new_peddata[new_peddata[,1]=="POI",2:3] <- NA
			new_peddata<-new_peddata[-which(new_peddata[,1] %in% c("fa","mo")),]					#delete fa and mo
		}
	}
	
	
	##plot
	if (plot){
	  ped <- FamiliasPedigree(id=new_peddata[,1], dadid=new_peddata[,2], momid=new_peddata[,3], sex=new_peddata[,4])
		plot(ped,affected=as.numeric(new_peddata[,5]))
	}
	return(new_peddata)

}

construct_H2 <- function(H1,known_parent="none"){
	if (class(H1)=="character"){
	  H2_0 <- matrix(c("POI","fa_0","mo_0",H1[4],1,
	                   "fa_0",NA,NA,"male",0,
	                   "mo_0",NA,NA,"female",0),nrow=3,byrow=TRUE)
	  return(H2_0)
	}
  
  H2 <- H1
  idx <- which(H1[,1]=="POI")
  sex <- H1[idx,4]
	H2[idx,5] <- 0
	H2[idx,1] <- "poi"
	
	H2[H1[,2]=="POI",2] <- "poi"
	H2[H1[,3]=="POI",3] <- "poi"
	if (known_parent=="none") H2_0 <- matrix(c("POI","fa_0","mo_0",sex,1,
	                 "fa_0",NA,NA,"male",0,
	                 "mo_0",NA,NA,"female",0),nrow=3,byrow=TRUE)

	if (known_parent=="fa") {
	    H2_0 <- matrix(c("POI","fa","mo_0",sex,1,
	                     "fa",NA,NA,"male",1,
	                     "mo_0",NA,NA,"female",0),nrow=3,byrow=TRUE)
	    H2[H1[,1]=="fa",1] <- "fa_0"
	    H2[H1[,2]=="fa",2] <- "fa_0"
	    H2[H1[,1]=="fa",5] <- 0
	    
	 }              
	           
	if (known_parent=="mo") {
	  H2_0 <- matrix(c("POI","fa_0","mo",sex,1,
	                   "fa_0",NA,NA,"male",0,
	                   "mo",NA,NA,"female",1),nrow=3,byrow=TRUE)

	  H2[H1[,1]=="mo",1] <- "mo_0"
	  H2[H1[,3]=="mo",3] <- "mo_0"
	  H2[H1[,1]=="mo",5] <- 0
	}     
	
	
	H2 <- rbind(H2_0,H2)
	
	#if (known_parent=="none")H2 <- rbind(H2,c("POI",NA,NA,sex,1))
	#if (known_parent=="fa")H2 <- rbind(H2,c("POI","fa",NA,sex,1))
	#if (known_parent=="mo")H2 <- rbind(H2,c("POI",NA,"mo",sex,1))
	#if (known_parent=="both")H2 <- rbind(H2,c("POI","fa","mo",sex,1))
	clean_ped0(peddata=H2,available=H2[H2[,5]=="1",1])
	
}

plot_pedigree <- function(ped){
	ids <- ped[is.na(ped[,2]) & is.na(ped[,3]),1]
	
	single_idx <- which(!(ids %in% ped[,2:3]))
	x <- FamiliasPedigree(id=ped[,1], dadid=ped[,2], momid=ped[,3], sex=ped[,4])
	x <- plot(x,affected=as.numeric(ped[,5]))
	if (length(single_idx)>0){
		n <- length(single_idx)
		points(x=max(x$x,na.rm=TRUE)-0.3*n,y=max(x$y,na.rm=TRUE)-0.2,pch=ifelse(ped[single_idx,4]=="male",16,15),cex=4,col="red")
		text(x=max(x$x,na.rm=TRUE)-0.3*n,y=max(x$y,na.rm=TRUE)-0.1,labels=ids[single_idx],pos=1)
	}
}


#############################
col_background = "teal"
box_status="primary"
#logs
#use_tracking()
##########################
ui <- dashboardPage(
  dashboardHeader(
	title = "EasyKin"
	#title = div(tags$a(img(src="Y:/tools/EasyKS/DNA.png",height=40,align="left")),style = "position: relative; top: -5px;","Easy kinship")
	#dropdownMenu(type="notifications",badgeStatus="warning",
	#			notificationItem(icon = icon("users"), status = "info",	
    #                             text =  "5 new members joined today",	
    #                            href = "https://www.r-project.org/"	
    #             ),	
    #             	
    #             notificationItem(icon = icon("warning"), status = "danger",	
    #                             text = "One error found!"	
    #             ))
	
  ),
 
  ## Sidebar menu
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "home", icon = icon("dashboard")),
	  menuItem("Kinship", tabName = "kinship", icon = icon("user-friends"),
			menuSubItem("Settings", tabName = "Settings_page"),
			menuSubItem("Pedigrees", tabName = "Pedigree_kinship"),
			menuSubItem("Real cases", tabName = "Real_cases")),		
      #menuItem("New", tabName = "New", icon = icon("th")),
	  menuItem("Manual", tabName = "manual", icon = icon("file-pdf")),
	  menuItem("Contact", tabName = "Contact", icon = icon("envelope"))
    )
  ),
  dashboardBody(
    tabItems(
	 tabItem(tabName = "home",		
		  h1("Welcome to EasyKin",aligh="center"),
		  br(),	##		 new line 
		  h2("Introduction"),
		  h4("	This is a useful online tool for forensic kinship analysis and missing person identification. There are some promissing features despite that several useful tools, such as EasyDNA, Bonaparte, Converge Software, Familias and forrel, have been developed.
			(a) This tool may be used to estimate system power before we perform kinship analysis,with available kits and references. (b) Pedigrees can be easily constructed by just clicking and including relative references.
			(c)Simple formulas, rather than complex algorithms, are applied automatically for pairwise kinship analyses, which makes it faster for simulations.(d) The system power of all possible combinations( subsets of available relative references) 
			 can be evaluated, making it easy to choose appropriate references."),
		  h4("Several parameters are calculated for the estimation of system power, including sensitivity(sen), specificity(spe), positive predictive value(PPV), negative predictive value(NPV), 
		    false positive rate(FPR), false negative rate(FNR),inconclusive and effectiveness. "),
		    h4(".        Sen: proportion of pedigrees under H1 judged as H1 true."),
			h4(".        Spe: proportion of pedigrees under H2 judged as H2 true."),
			h4(".        PPV: proportion of pedigrees correctly judged as H1 true."),
			h4(".        NPV: proportion of pedigrees correctly judged as H2 true."),
			h4(".        FPR: proportion of pedigrees under H2 judged as H1 true. "),
			h4(".        FNR: proportion of pedigrees under H1 judged as H2 true."),
			h4(".        Inconclusive: proportion of pedigrees that cannot be judged as either H1 true or H1 true."),
			h4(".        Effectiveness: proportion of pedigrees that cannot be judged as H1 true or H2 true."),
			#(c) Both singel and multiple cases can be analyzed. 
		  br(),
		  textOutput("loci_txt"),
		  h4("Format of frequency data"),
		  fluidRow( box(dataTableOutput("freq_format_table"),width=6)),
		  #br(),
		  h4("Format of genotype data"),
		  fluidRow(box(dataTableOutput("genotype_format_table"),width=6)),
		  h4("Format of STR mutation models"),
		  fluidRow(box(tableOutput("mutation_format_table"),width=6)),
		  h4("Format of user-defined pedigrees (first cousin)"),
		  fluidRow(box(tableOutput("pedigree_format_table"),width=6))
      ),
	 ###settings 
	 tabItem(tabName = "Settings_page",	
		column(width=12,	
							box(title="Kits and markers",width=12,
								fileInput('file_frequency_ped', 'Frequency data file (62 markers are available by default)',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv')),
								checkboxInput("check_strict_marker_ped","Restrict markers in following selected kit(s) ",TRUE),
											
								selectInput("select_kits_ped", "Choose the kit(s)", c("AmpFlSTR Identifiler"="AmpFlSTR Identifiler",
								                      "AmpFlSTR Sinofiler"="Sinofiler",
																			"AmpFlSTR GlobalFiler"="AmpFlSTR GlobalFiler",
																			"Huaxia"="HUAXIA",																			
																			"PowerPlex 21"="PowerPlex 21",
																			"Goldeneye 25A"="Goldeneye 25A",
																			"Goldeneye 30A"="Goldeneye 30A",
																			"Goldeneye 22NC"="22NC",
																			"Microreader 23sp"="23SP",
																			"AGCU 21+1"="21+1",
																			"CODIS 13"="CODIS_13",
																			"CODIS 20"="CODIS_20"),selected="HUAXIA",multiple=TRUE),
								#textAreaInput("Input_markerID","Input marker(s)",height=100),
								selectInput("Input_markerID","Input marker(s)",c("TH01"),multiple=TRUE),
								textOutput("lociNum_text")
							),
								
							box(title="Mutation",width=12,
								column(width=4,
									selectInput("MutationModel_ped","Model",c("Equal","Stepwise")),
									numericInput("MutationRate_ped1", "Range", value = 0.1, min = 0, max = 1,step=0.0001),
									
								),
								column(width=4,
									
									numericInput("MutationRate_ped0", "Rate1", value = 0.002, min = 0, max = 1,step=0.0001),
									
									numericInput("MutationRate_ratio", "Ratio (male/female)", value = 3, min = 1, max = 10,step=0.1),
									
									#numericInput("MutationRate_ratio", "Ratio of male to female", value = 3, min = 1, max = 10,step=0.1),
								),
								column(width=4,
									numericInput("MutationRate_ped2", "Rate2 (non_integer)", value = 10^-8, min = 0, max = 1,step=10^-9),
									fileInput('Mutation_Customized_file', 'Modeling STR individually',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))
								)
								
							),
							box(title="Simulation",width=12,
							#column(width=6,numericInput("Input_theta", "subpopulation effects", value = 0, min = 0, max = 1,step=0.01)),
							
							column(width=6,numericInput("set_seeds", "Set seeds", value = 123, min = 1, max = 10000,step=1)),
							column(width=6,numericInput("num_Simulations_ped", "Number of simulations", value = 100, min = 10, max = 5000,step=10))
							)
			)		
					
      ),
	    
	##pedigree kinship settings
      tabItem(tabName = "Pedigree_kinship",		
		fluidRow(
			column(width=4,
				
							box(title="Person of interest (POI)",width=12,
								column(width=6,selectInput("questioned_sex", "sex ", c("male" = "male" , "female" = "female" ),multiple=FALSE)),
								column(width=6,selectInput("select_known_parent", "known parent", c("father" = "fa" , "mother" = "mo" ,"none"="none"),selected="none",multiple=FALSE)),
															
							),
							#helpText("Reference samples"),
							box(title="Reference samples",width=12,
								column(width=6,
								  selectInput("select_ref_parent", "parent (s)", c("father" = "fa" , "mother" = "mo" ),multiple=TRUE),
								  selectInput("select_ref_FS", "full siblings", c("brother 1" = "fs_b1", "brother 2" = "fs_b2","brother 3" = "fs_b3",
								                                                  "sister 1" = "fs_s1", "sister 2" = "fs_s2","sister 3" = "fs_s3"),multiple=TRUE),
									selectInput("select_ref_pGG", "paternal grandparents", c("grandfather" = "pgf", "grandmother" = "pgm" ),multiple=TRUE),
									selectInput("select_ref_pAV", "paternal avuncular", c("uncle 1" = "puncle1" , "uncle 2" = "puncle2" ,"uncle 3" = "puncle3" ,
																							"aunt 1" = "paunt1", "aunt 2" = "paunt2","aunt 3" = "paunt3"),multiple=TRUE),
									selectInput("select_ref_pHS", "paternal half siblings",c("brother 1" = "phs_b1","brother 2"="phs_b2","brother 3"="phs_b3",
									                                                     "sister 1" = "phs_s1","sister 2" = "phs_s2","sister 3" = "phs_s3"),multiple=TRUE,width = "100%"),
									
									selectInput("select_ref_gchild_son", "grandchild (son)",c("grandson 1" = "gson_s1","grandson 2"="gson_s2","grandson 3"="gson_s3",
									                                                         "granddaughter 1" = "gdaug_s1","granddaughter 2" = "gdaug_s2","granddaughter 3" = "gdaug_s3"),multiple=TRUE,width = "100%"),
									selectInput("select_ref_nephew_bro", "nephew niece (bro)",c("nephew 1" = "nephew_b1","nephew 2"="nephew_b2","nephew 3"="nephew_b3",
									                                                                "niece 1" = "niece_b1","niece 2" = "niece_b2","niece 3" = "niece_b3"),multiple=TRUE,width = "100%"),
									
									
									
																		
									#selectInput("radioButton_as_FS","Sibling as: ", c("paternal half sibling" =1 , "maternal half sibling" =2 ,"full sibling"=3),selected=1),																			
									#radioButtons("radioButton_as_FS","As ", c("paternal half sibling" =1 , "maternal half sibling" =2 ,"full sibling"=3),selected=1)																			
								),
												
								column(width=6,
								  selectInput("select_ref_children", "children", c("son 1" = "son1", "son 2" = "son2","son 3" = "son3",
								                                                        "daughter 1" = "daughter1", "daughter 2" = "daughter2","daughter 3" = "daughter3"),multiple=TRUE),
								  selectInput("select_ref_spouse", "genetically unrelated", c("spouse" = "spouse", 
								                                               "mother of paternal half sibling" = "hs_mo","father of maternal half sibling" = "hs_fa",
								                                               "daughter in law" = "daughter_inlaw","son in law" = "son_inlaw",
								                                               "sister in law" = "sis_inlaw","brother in law" = "bro_inlaw"),multiple=TRUE),
								  selectInput("select_ref_mGG", "maternal grandparents", c("grandfather" = "mgf", "grandmother" = "mgm"),multiple=TRUE),
									selectInput("select_ref_mAV", "maternal avuncular", c("uncle 1" = "muncle1", "uncle 2" = "muncle2","uncle 3" = "muncle3" ,
																						"aunt 1" = "maunt1", "aunt 2" = "maunt2","aunt 3" = "maunt3"),multiple=TRUE),																					
									#selectInput("select_ref_hsmo", "mother of pHS", c("available" = "hs_mo", "not available" = "not available"),selected="not available",multiple=FALSE),
									#selectInput("select_ref_hsfa", "father of mHS", c("available" = "hs_fa", "not available" = "not available"),selected="not available",multiple=FALSE),
									#selectInput("select_ref_daug_inlaw", "daughter in law", c("available" = "daughter_inlaw", "not available" = "not available"),selected="not available",multiple=FALSE),
									#selectInput("select_ref_son_inlaw", "son in law", c("available" = "son_inlaw", "not available" = "not available"),selected="not available",multiple=FALSE),
									#selectInput("select_ref_sis_inlaw", "sister in law", c("available" = "sis_inlaw", "not available" = "not available"),selected="not available",multiple=FALSE),
									#selectInput("select_ref_bro_inlaw", "brother in law", c("available" = "bro_inlaw", "not available" = "not available"),selected="not available",multiple=FALSE),
									selectInput("select_ref_mHS", "maternal half siblings",c("brother 1" = "mhs_b1","brother 2"="mhs_b2","brother 3"="mhs_b3",
									                                                         "sister 1" = "mhs_s1","sister 2" = "mhs_s2","sister 3" = "mhs_s3"),multiple=TRUE,width = "100%"),
									selectInput("select_ref_gchild_daug", "grandchild (daughter)",c("grandson 1" = "gson_d1","grandson 2"="gson_d2","grandson 3"="gson_d3",
									                                                                "granddaughter 1" = "gdaug_d1","granddaughter 2" = "gdaug_d2","granddaughter 3" = "gdaug_d3"),multiple=TRUE,width = "100%"),
									selectInput("select_ref_nephew_sis", "nephew niece (sis)",c("nephew 1" = "nephew_s1","nephew 2"="nephew_s2","nephew 3"="nephew_s3",
									                                                            "niece 1" = "niece_s1","niece 2" = "niece_s2","niece 3" = "niece_s3"),multiple=TRUE,width = "100%"),
									
								),
							),
							box(width=12,fileInput('file_peds', 'User defined pedigrees',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
							#actionButton("nonhalf_Button", "plot pedigree"),
							box(width=12,
								actionButton("sys_ped_Button", "OK"),
							),
							textOutput("simulation_runtime_txt")
		
			),
			column(width=8,
				box(title="Hypotheses",
					width=12,
					plotOutput("plot_ped_pedigree", height = 300)
				),
				box(
					title="Likelihood ratio(LR) distribution",
					width=12,
					plotOutput("plot_KIdistribution_ped", height = 230),
					textOutput("non_overlap_txt_ped"),
					sliderInput("slider_thrshold_ped", "Threshold", -10, 10, 1,step=0.01),
					checkboxInput("check_singleT","Single threshold",FALSE,width="25%"),
					tableOutput("sys_power_table_ped"),
					textOutput("Decision_threshold"),
					downloadButton("downloadPed", "Download pedigree"),
					downloadButton("downloadKIdis", "Download LR distribution"),
					downloadButton("downloadPedLR", "Download pedigree and LR distribution")
				),
				box(title="All possible pedigrees",width=12,collapsible=TRUE,collapsed=TRUE,
					#plotOutput("plot_allped_sys"),
					actionButton("prune_ped_Button", "Prune pedigree"),
					helpText("This process will cost much more time than single pedigree."),
					dataTableOutput("allped_sys_table"))
			
			)
		
		)
	  ),
	  ##single real case
	  tabItem(tabName = "Real_cases",
				box(title="",width=4,					
					column(width=12,  fileInput('file_input_genotype', 'Input genotype data',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
					column(width=6,numericInput("NovelAlleleFreq", "Frequency of novel alleles", value = 0.01, min = 0.000001, max = 1,step=0.0001)),
					column(width=6,numericInput("prior_prob", "Prior probability", value = 0.5, min = 0, max = 1,step=0.01)),
					column(width = 6,selectInput("select_old_labels", "old labels", c(""),multiple=FALSE)),  
					column(width = 6,selectInput("select_new_labels", "new labels", c("" ),multiple=FALSE)), 
					#column(width=6, fileInput('file_peds_real', 'Input new pedigrees',accept=c('text/csv', 'text/comma-separated-values,text/plain', '.csv'))),
					tableOutput("label_old2new_table"),
					column(width = 4,actionButton("realcase_Button", "Calculate LRs")),
					column(width = 4,actionButton("clearLabel_Button", "Clear labels")),
					column(width = 4,actionButton("changeLabel_Button", "Assign labels")),
					
					
					
				),
				
				box(title="Hypotheses:",width=8,					
					plotOutput("plot_ped_pedigree_real", height = 300)
				),
				tabBox(title="",side="left",width=12,selected="Genotype",#height="500px",
				
					tabPanel("Genotype",
						
						dataTableOutput("realcase_genotype_table")
					),					
					tabPanel("LRs",
					  
						dataTableOutput("realcase_LR_table"),
						textOutput("realcase_LR_txt")
					),
					tabPanel("Pairwise LRs", tableOutput("LRdetail_table"))
					
					
				),
				tableOutput("LR_interpretation_table"),
				 tags$head(tags$style("#LR_interpretation_table{color: red;
                                 font-size: 20px;
                                 font-style: italic;
                                 }"
                         )
              )
      ),
	  #          
	  tabItem(tabName = "manual",
	         
	        box(title="manual in English",width=12,collapsible=TRUE,collapsed=TRUE,uiOutput("sraManual")),
	        box(title="manual in Chinese",width=12,collapsible=TRUE,collapsed=TRUE,uiOutput("sraManual_Chinese")),
	       
  	),
      tabItem(tabName = "Contact",
		 h3("If you have any question or suggestion, please feel free to contact us"),
		 br(),	##		 new line
         h3("Ran Li,"),		 
		 h3("E-mail: 654163930@qq.com, liran6@mail2.sysu.edu.cn"),
		 h3("Hongyu Sun,"),		 
		 h3("E-mail: sunhy@mail.sysu.edu.cn, sunhongyu2002@163.com"),
		 br(),	##		 new line
		 h4("Faculty of Forensic Medicine, Zhongshan School of Medicine, Sun Yat-sen University"),
		 h4("Guangzhou 510080, Guangdong ,People’s Republic of China"),
		 br(),	##		 new line
		 tags$b("log errors:"),
		 verbatimTextOutput(outputId = "log_errors"),
		 tags$b("All inputs:"),
		 verbatimTextOutput(outputId = "all_inputs")
      )
    )
	
  )
)

	
	

server <- function(input, output,session) {
	#################
  track_usage(storage_mode = store_null(console = FALSE))
	#initialization
	wd_path<-getwd()
	#peds<-initialize_peds()
	#peds_references<-initialize_ped_references()
	freqdata0<-read.table(file=paste0(wd_path,"/data/Example of frequency.txt"), header = TRUE, sep= "\t", stringsAsFactors = FALSE)
	genodata0<-read.table(file=paste0(wd_path,"/data/Example of genotypeA.txt"), header = TRUE, sep= "\t", stringsAsFactors = FALSE)
	#genodata0 <- matrix(c("TH01","D12S391","11,12","10,12","12,12","13,14"),nrow = 2,dimnames = c(list(c("TH01","D12S391")),list(c("loci","POI","S1"))))
	kits_set0<-read.table(file=paste0(wd_path,"/data/kits.txt"), header = TRUE, sep= "\t", stringsAsFactors = FALSE)
	mutation_model0<-read.table(file=paste0(wd_path,"/data/Example of mutation model.txt"), header = TRUE, sep= "\t", stringsAsFactors = FALSE)
	pedigree_fc<-read.table(file=paste0(wd_path,"/data/Example of ped.txt"), header = FALSE, sep= "\t", stringsAsFactors = FALSE)
	
	H1andH2 <-reactiveValues(ped=NULL)
	old2new_matrix <-reactiveValues(Record=NULL)
	
	##################
	##home page
	output$loci_txt<-renderText({
		loci<-sort(as.character(unique(freqdata0[,1])))
		locistr<-paste0(loci,collapse=", ")
		paste("Frequency data of ",length(loci)," STRs from Chinese Han population are provided, i.e., ",locistr,". Users can also upload your own data 
		as long as you follow the required formats.",sep="")}
	)
	output$freq_format_table<- renderDataTable({datatable(freqdata0[1:24,],rownames=FALSE,editable=FALSE,extensions=c("Buttons",'Scroller'),
                options = list(dom = 'Bfrtip',
                               scrollY = 200,
                               scroller = TRUE,
                               scrollX=TRUE,
							   buttons = c('csv')
                ))})
	output$genotype_format_table<- renderDataTable({datatable(genodata0,rownames=FALSE,editable=FALSE,extensions=c("Buttons",'Scroller'),
                options = list(dom = 'Bfrtip',
                               scrollY = 200,
                               scroller = TRUE,
                               scrollX=TRUE,
							   buttons = c('csv')		#, 'excel'
                ))})
	output$mutation_format_table<- renderTable(na="NA",digits=-2,colnames=TRUE,{mutation_model0})
	output$pedigree_format_table<- renderTable(na="NA",colnames=FALSE,{pedigree_fc})
	
	#################
	##pedigree_kinship
	
	refresh_locus<-reactive({		
		inFile <- input$file_frequency_ped	
		if (is.null(inFile)) freqdata<-freqdata0 else freqdata<-read.table(inFile$datapath, header = TRUE, sep= "\t", stringsAsFactors = FALSE)
		mutation_metrics<-NULL
		inFile2 <- input$Mutation_Customized_file	
		if (!is.null(inFile2)) mutation_metrics<-read.table(inFile2$datapath, header = TRUE, sep= "\t", stringsAsFactors = FALSE,row.names=1)
			
		lociallelefreq<-list_freqdata(freqdata)
		#partial markers
		if (!!input$check_strict_marker_ped=="TRUE"){
			kits<-input$select_kits_ped
			input_markers<-input$Input_markerID
			if (length(kits)==0 & length(input_markers)==0) kits<-as.character(unique(kits_set0$kits))
			submarkers<-as.character(unique(kits_set0[as.character(kits_set0$kits) %in% kits,2]))
			
			idx<-which(names(lociallelefreq) %in% c(submarkers,input_markers))
			lociallelefreq<-lociallelefreq[idx]
		}
		MutationModel_ped<-input$MutationModel_ped
		MutationRate_ped0<-input$MutationRate_ped0
		MutationRate_ped1<-input$MutationRate_ped1
		MutationRate_ped2<-input$MutationRate_ped2
		mutRatio<-input$MutationRate_ratio
		
		myloci<-list()
		for (i in 1:length(lociallelefreq)){
			freq<-lociallelefreq[[i]]
			freq<-freq/sum(freq)
			locus<-names(lociallelefreq)[i]	
			models<-MutationModel_ped
			rate1<-MutationRate_ped0
			rate2<-MutationRate_ped2
			ranges<-MutationRate_ped1
			ratios<-mutRatio
			#print(locus)
			
			##user difined models
			if (!is.null(mutation_metrics)){
				markerID<-which(rownames(mutation_metrics)==locus)
				if (length(markerID)>0){
					models<-as.character(mutation_metrics[markerID,1])
					rate1<-as.numeric(as.character(mutation_metrics[markerID,2]))
					rate2<-as.numeric(as.character(mutation_metrics[markerID,3]))
					ranges<-as.numeric(as.character(mutation_metrics[markerID,4]))
					ratios<-as.numeric(as.character(mutation_metrics[markerID,5]))
				}
			}
			#model STR markers			
			x<-FamiliasLocus(frequencies=freq,allelenames=names(freq), name=locus, 
						MutationModel=models,MutationRate=rate1,
						MutationRange=ranges,MutationRate2=rate2,
						femaleMutationRate=rate1/ratios)
			
			myloci<-c(myloci,list(x))
		}
		names(myloci)<-names(lociallelefreq)
		myloci
	})
	#update the markers in the selection box
	observe({
		inFile <- input$file_frequency_ped	
		if (is.null(inFile)) freqdata<-freqdata0 else freqdata<-read.table(inFile$datapath, header = TRUE, sep= "\t", stringsAsFactors = FALSE)
		STRs<-sort(unique(as.character(freqdata[,1])))
		updateSelectInput(session, "Input_markerID",
			label = "Input marker(s)",
			choices = STRs,
			selected = NULL #head(now_selected,1)
		)
	})
	
	#update the sample labels of pedigree
	observe({
	  ped<-H1andH2$ped
	  samplenames1 <- ped$H1[[1]]
	  samplenames2 <- ped$H2[[1]]
	  samplenames <- samplenames1[samplenames1 %in% samplenames2]
	  updateSelectInput(session, "select_old_labels",
	                    label = "old labels",
	                    choices = sort(samplenames),
	                   
	                    selected = NULL #head(now_selected,1)
	  )
	})
	#update the sample labels of input genotype
	observe({
	  geno <- refresh_real_genotype()
	  samplenames <- colnames(geno)
	  updateSelectInput(session, "select_new_labels",
	                    label = "new labels",
	                    choices = samplenames,
	                    selected = NULL #head(now_selected,1)
	  )
	})
	#
	# user defined pedigree
	observe({
	  inFile <- input$file_peds	
	  if (!is.null(inFile)) {			
	    pedData<-read.table(file=inFile$datapath,header=FALSE,sep="\t",stringsAsFactors=FALSE)
	    H1andH2$ped <- inputPedigree(peddata=pedData)
	  }
	})
	
	#refresh_peds
	observe({ 
		target_sex<-input$questioned_sex
		known_parent <- input$select_known_parent
		ids_available <- c(
		  "POI",
		  known_parent,
		  input$select_ref_pGG,
		  input$select_ref_mGG,
		  input$select_ref_pAV,
		  input$select_ref_mAV,
		  input$select_ref_pHS,
		  input$select_ref_mHS,
		  #input$select_ref_hsfa,
		  #input$select_ref_hsmo,
		  input$select_ref_gchild_son,
		  #input$select_ref_daug_inlaw,
		  input$select_ref_gchild_daug,
		  #input$select_ref_son_inlaw,
		  input$select_ref_nephew_bro,
		  #input$select_ref_sis_inlaw,
		  input$select_ref_nephew_sis,
		  #input$select_ref_bro_inlaw,
		  input$select_ref_children,
		  input$select_ref_spouse,
		  input$select_ref_parent,
		  input$select_ref_FS
		)
		
		ped0 <- initialize_ped0(target_sex)
		H1 <- clean_ped(peddata = ped0,available = ids_available)
		if (class(H1)=="character") H1<- matrix(c("POI","fa","mo",H1[4],1,
		                 "fa",NA,NA,"male",0,
		                 "mo",NA,NA,"female",0),nrow=3,byrow=TRUE)
		
		H2 <- construct_H2(H1,known_parent=known_parent)

		
		h1_affected <- as.numeric(H1[,5])
		h2_affected <- as.numeric(H2[,5])
	
		H1 <- FamiliasPedigree(id=H1[,1], dadid=H1[,2], momid=H1[,3], sex=H1[,4])
		H2 <- FamiliasPedigree(id=H2[,1], dadid=H2[,2], momid=H2[,3], sex=H2[,4])
		
		#rename
		samplenames<-levels(factor(c(H1$id,H2$id)))
		#print(samplenames)
		for (i in 2:length(H1$id)){id<-which(samplenames==H1$id[i]);H1$id[i]<-paste0("S",id)}
		for (i in 2:length(H2$id)){id<-which(samplenames==H2$id[i]);H2$id[i]<-paste0("S",id)}
		H1andH2$ped <- c(H1=list(H1),H2=list(H2),h1_affected=list(h1_affected),h2_affected=list(h2_affected))	
	})
	
	LRs_peds<-eventReactive(input$sys_ped_Button,{
		set.seed(input$set_seeds)
		loci<-refresh_locus()
		pedigrees<-H1andH2$ped
		theta<-0                      #input$Input_theta
		n<-input$num_Simulations_ped
		runTime<-system.time({LRs<-simulation_pedsLR_batch(peds=pedigrees,loci=loci,n=n,theta=theta)})
		#print(runTime)
		if (length(LRs)!=0)LRs=LRs[c(!is.infinite(LRs[,1]) & !is.infinite(LRs[,2])),]								#exclude infinite rows	
		
		c(LRs=list(LRs),runTime=list(runTime),n_ped=list(n),n_STR=list(length(loci)))
		#c(LRs=list(read.table(file="Y:/LR_case.txt",header = TRUE,sep="\t",stringsAsFactors = FALSE)),runTime=list(runTime))
		
	})
	
	LRs_allpeds<-eventReactive(input$prune_ped_Button,{	
		set.seed(input$set_seeds)
		loci<-refresh_locus()
		#pedigrees<-refresh_peds()
		pedigrees<-H1andH2$ped
		theta<-0                   ##input$Input_theta
		n<-input$num_Simulations_ped
		
		availableID<-which(pedigrees[[3]]==1)
		available<-pedigrees[[1]]$id[availableID]
	
		runTime<-system.time({
			combs_LRs<-prune_pedigree(peds=pedigrees,loci=loci,available=available,POI="POI",refID=2,n=n,theta=theta,tx=0)			
		})
		c(LRs=list(combs_LRs),runTime=list(runTime))
		
	})
	output$simulation_runtime_txt<-renderText({
		x=LRs_peds()
		x<-x[[2]]
		paste0("Finished. Time cost:   ",round(x[3],2)," seconds.")
	})
	output$plot_KIdistribution_ped<-renderPlot({
		y=LRs_peds()
		x<-y[[1]]
		
		if (length(x)==1)return()
		n_ped<-y$n_ped
		n_STR<-y$n_STR
		thresholds<-input$slider_thrshold_ped
		singelT<-input$check_singleT
		if(singelT)tx<-c(thresholds,thresholds) else tx<-sort(c(-thresholds,thresholds))
		par(mar=c(4,4,2,5))
		if (length(x)!=0)multiHist(x=x,xlab=expression("Log"[10]*"LR"),main=paste("Peds: ",n_ped, ",STRs: ",n_STR,sep=""),thresholds=tx)
	})
	output$plot_ped_pedigree<-renderPlot({
		plotpeds()
	})
	output$sys_power_table_ped<-renderTable(striped=TRUE,width="100%",digits=8,{
		x=LRs_peds()
		
		x<-x[[1]]
		if (length(x)<=1)return()		
		thresholds<-input$slider_thrshold_ped
		singelT<-input$check_singleT
		if(singelT)tx<-c(thresholds,thresholds) else tx<-sort(c(-thresholds,thresholds))
		
		sys<-syspower_batch(x=x,t=tx,fit=TRUE)
		sys
		}
	)
	output$Decision_threshold <- renderText({
	  x=LRs_peds()
	  
	  x<-x[[1]]
	  if (length(x)<=1)return()
	  sys <- syspower2(x,1,1)
	  paste0("Decision threshold(DT):",round(sys$DT,2),",    FPR:",sys$error[1],",    FNR:",sys$error[2])
	})
	output$allped_sys_table<-renderDataTable({
		thresholds<-input$slider_thrshold_ped
		singelT<-input$check_singleT
		if(singelT)tx<-c(thresholds,thresholds) else tx<-sort(c(-thresholds,thresholds))
		combs_LRs<-LRs_allpeds()
		combs_LRs<-combs_LRs[[1]]
		sys_combination<-c()
			#print(LRs)
			for (i in 1:length(combs_LRs)){
				if (length(combs_LRs[[i]])==1){sys_combination<-rbind(sys_combination,rep(NA,8));next}
				if (sum(combs_LRs[[i]])==0){sys_combination<-rbind(sys_combination,rep(NA,8));next}
				x<-syspower2(combs_LRs[[i]],t=tx)
				x<-x$sys_parameters
				sys_combination<-rbind(sys_combination,x)
			}
			rownames(sys_combination)<-names(combs_LRs)
		datatable(round(sys_combination,6),editable=FALSE,extensions=c("Buttons",'Scroller'),
                options = list(dom = 'Bfrtip',
                               scrollY = 300,
                               scroller = TRUE,
                               scrollX=TRUE,
							   buttons = c('csv')
                ))
		
	})
	
	plotpeds<-function(){
		par(mfrow=c(1,2))
		#x=refresh_peds()
		x<-H1andH2$ped
		ped1<-x$H1
		ped2<-x$H2
		if (length(which(x$h1_affected==1))<2) return()
		status1 <- sapply(x$h1_affected,function(x){ifelse(x==0,1,0)})
		status2 <- sapply(x$h2_affected,function(x){ifelse(x==0,1,0)})
		ids1 <- ped1$id;ids1[x$h1_affected==0] <- ""
		ids2 <- ped2$id;ids2[x$h2_affected==0] <- ""
		
		cols1<-c("red",rep("black",length(ped1$id)-1))
		cols2<-c("red",rep("black",length(ped2$id)-1))
		plot(ped1,mar=c(4,4,4,4),affected=x$h1_affected,angle=45, density=25,col=cols1,status=status1,id=ids1);title("H1")
		plot(x$H2,mar=c(4,4,4,4),affected=x$h2_affected,angle=45, density=25,col=cols2,status=status2,id=ids2);title("H2")
	}
	plotLRs<-function(){
		x=LRs_peds()
		x<-x[[1]]
		if (length(x)==1)return()
		thresholds<-input$slider_thrshold_ped
		par(mar=c(4,4,2,5))
		if (length(x)!=0)multiHist(x=x,xlab="Log10(LR)",main="",thresholds=c(-thresholds,thresholds))
	}
	output$downloadPed <- downloadHandler(
		filename ="pedigree.png",
		content = function(file) {
			png(file,width = 900, height = 300)
				print(plotpeds())
			dev.off()
		}
	)
	output$downloadKIdis <- downloadHandler(
		filename ="LR distribution.png",
		contentType="image/png",
		content = function(file) {
			png(file,width = 800, height = 300)
				plotLRs()
			dev.off()
		}
	)
	output$downloadPedLR <- downloadHandler(
	  filename ="Pedigree and LR distribution.png",
	  contentType="image/png",
	  content = function(file) {
	    png(file,width = 800, height = 600)
	    layout(mat = matrix(c(1,2,3,3),nrow=2,byrow = T))
	    par(mar=c(3,4,4,4))
	    print(plotpeds())
	    par(mar=c(4,4,2,5))
	    print(plotLRs())
	    dev.off()
	  }
	)
	
	output$lociNum_text<-renderText({
		x=refresh_locus()
		locusnames<-c()
		for (i in 1:length(x))locusnames<-c(locusnames,x[[i]][[1]])
		locusnames<-sort(locusnames)
		locusnames<-paste(locusnames,collapse=", ")		
		paste0("You will type ", length(x)," loci, including ",locusnames)
	})
	
	
	#####################################################################################################
	##Real cases
	refresh_real_peds<-reactive({H1andH2$ped})
	refresh_real_genotype<-reactive({	
		inFile <- input$file_input_genotype	
		if (is.null(inFile)) return() 
		filepath<-inFile$datapath
		read.table(inFile$datapath, header = TRUE, sep= "\t", stringsAsFactors = FALSE,row.names=1)
	})
	pedigree_adjusted<-function(genotype,pedigrees){
		if (length(genotype)==0) return(c(genotype=list(genotype),pedigrees=list(pedigrees)))
		samplenames<-colnames(genotype)
		idx1<-pedigrees[[1]]$id %in% samplenames
		pedigrees$h1_affected<-as.numeric(idx1)
		s<-pedigrees[[1]]$id[idx1]
		idx2<-pedigrees[[2]]$id %in% s
		pedigrees$h2_affected<-as.numeric(idx2)
		genotype<-genotype[,colnames(genotype) %in% s]
		c(genotype=list(genotype),pedigrees=list(pedigrees))
	}
	#change sample labels
	observeEvent(input$changeLabel_Button,{
	  old_label <- input$select_old_labels
	  new_label <- input$select_new_labels
	  x <- old2new_matrix$Record
	  if (is.null(x)){
	      x<- matrix(c(old_label,new_label),nrow=1)
	  }else{
	    if (old_label %in% x[,1]){
	      x[x[,1]==old_label,2] <- new_label
	    }else{
	      x <-rbind(x,c(old_label,new_label))
	    }
	  }
	  colnames(x) <- c("labels in the pedigree","labels in the input data")
	  old2new_matrix$Record <- x
	})
	#clear labels
	observeEvent(input$clearLabel_Button,{old2new_matrix$Record <- NULL})
	
	relabeled_geno <- reactive({
	  geno <- refresh_real_genotype()
	  label_matrix <- old2new_matrix$Record
	  if (is.null(label_matrix)) return(geno)
	  old_names <- colnames(geno)
	  new_names <- old_names
	  for (i in 1:length(new_names)){
	    s <- new_names[i]
	    if (s %in% label_matrix[,2]){
	      new_names[i] <- label_matrix[label_matrix[,2] == s,1]
	    }
	  }
	  colnames(geno) <- new_names
	  geno
	})
	
	#calculate LRs
	refresh_caseLR<-eventReactive(input$realcase_Button,{		
		genotype<-relabeled_geno()
		pedigrees<-refresh_real_peds()
		loci<-refresh_locus()
		#exclude samples not in the pedigree
		if (length(genotype)>0){
			new_pedG<-pedigree_adjusted(genotype,pedigrees)
			genotype<-new_pedG$genotype
			pedigrees<-new_pedG$pedigrees
		}		
		if (length(genotype)==0) return("No genotypes available for all samples")
		genotype0<-genotype
		#excluded loci with NA genotype
		na_loci<-apply(genotype,MARGIN=1,FUN=function(x)any(is.na(x) | x==""))
		genotype<-genotype[!na_loci,]
		genotype<-t(genotype)	
		rare_freq<-input$NovelAlleleFreq
		
		x=transform_genotype(genotype=genotype,loci=loci,rare_freq=rare_freq)
		genotype<-x$genotype
		loci<-x$loci
		#details
		sample_pair <- combn(1:nrow(genotype),2)
		LR_details <- vector("list",ncol(sample_pair))
		names(LR_details) <- apply(sample_pair, MARGIN=2, FUN=function(x,y=rownames(genotype)){paste(y[x],collapse = "_")})
		for (i in 1:length(LR_details)){
		  s1 <- rownames(genotype)[sample_pair[1,i]]
		  s2 <- rownames(genotype)[sample_pair[2,i]]
		  geno <- genotype[sample_pair[,i],]
  
		  kvalues<-kPairs(pedigrees[1],c(s1,s2))
		  kvalues$k2 <- c(0,0,1)
		  LR_details[[i]] <- cal_pairsLR(genotype=geno,k=kvalues,loci=loci)
		 
		}
		#
		
	
		c(LRs=list(cal_pedsLR(peds=pedigrees,loci=loci,refID=2,genotype=genotype)),
		genotype=list(genotype0),
		pedigrees=list(pedigrees),
		details=list(LR_details))
		
	})
	output$plot_ped_pedigree_real<-renderPlot({	
		pedigrees=refresh_real_peds()
		genotype=relabeled_geno()
	
		if (length(genotype)>0){
			new_pedG<-pedigree_adjusted(genotype,pedigrees)
			genotype<-new_pedG$genotype
			pedigrees<-new_pedG$pedigrees
		}
		ped1<-pedigrees$H1
		ped2<-pedigrees$H2
		x <- old2new_matrix$Record
		sample_lables1 <- ped1[[1]]
		sample_lables2 <- ped2[[1]]
		#change the sample labels
		for (i in 1:length(sample_lables1)){if (sample_lables1[i] %in% x[,1])sample_lables1[i] <- x[x[,1] == sample_lables1[i],2]}
		for (i in 1:length(sample_lables2)){if (sample_lables2[i] %in% x[,1])sample_lables2[i] <- x[x[,1] == sample_lables2[i],2]}
		cols1<-c("red",rep("black",length(ped1$id)-1))
		cols2<-c("red",rep("black",length(ped2$id)-1))
		status1 <- sapply(pedigrees$h1_affected,function(x){ifelse(x==0,1,0)})
		status2 <- sapply(pedigrees$h2_affected,function(x){ifelse(x==0,1,0)})
		
		
		par(mfrow=c(1,2))
		plot(ped1,mar=c(4,4,4,4),affected=pedigrees$h1_affected,angle=45, density=25,col=cols1,id=sample_lables1,status=status1);title("H1")
		plot(ped2,mar=c(4,4,4,4),affected=pedigrees$h2_affected,angle=45, density=25,col=cols2,id=sample_lables2,status=status2);title("H2")
	
	})
	output$label_old2new_table<-renderTable(striped=TRUE,width="100%",colnames=TRUE,align="l",{
	  old2new_matrix$Record
	  
	})
	output$realcase_genotype_table<-renderDataTable({
		datatable(refresh_real_genotype(),editable=FALSE,extensions=c("Buttons",'Scroller'),
                options = list(dom = 'Bfrtip',
                               scrollY = 300,
                               scroller = TRUE,
                               scrollX=TRUE,
							   buttons = c('csv')
                ))}
	)
	output$realcase_LR_table<-renderDataTable({		
		x<-refresh_caseLR()		
		if (length(x)==1) return(matrix(x,1,1))
		LRs<-x$LRs
		genotype<-x$genotype

		LR<-round(LRs$LRs,6)
		loci_LR<-LR[rownames(genotype)]
		loci_LR[is.na(loci_LR)]<-"No genotypes or no frequency data available"
		loci_LR<-cbind(genotype,LR=loci_LR)
				
		datatable(loci_LR,editable=FALSE,extensions=c("Buttons",'Scroller'),
                options = list(dom = 'Bfrtip',
                               scrollY = 300,
                               scroller = TRUE,
                               scrollX=TRUE,
							   buttons = c('csv')
                ))}	
	)
	output$LRdetail_table<-renderTable(striped=TRUE,width="100%",digits=4,rownames=TRUE,colnames=TRUE,align="l",{
	  x <- refresh_caseLR()$details
	  locus <- c()
	  LRs <- c()
	  for (i in 1:length(x)){locus <- c(locus,names(x[[i]]))}
	  locus <- unique(locus)
	  for (i in 1:length(x)){
	    lr <- x[[i]]
	    LRs <- cbind(LRs,lr[locus])
	    
	  }
	  CLR <- apply(LRs,2,prod,na.rm=TRUE)
	  LRs <- rbind(LRs,CLR)
	  colnames(LRs) <- names(x)
	  rownames(LRs) <- c(locus,"CLR")
	  LRs
	})
	output$LR_interpretation_table<-renderTable(striped=TRUE,width="50%",digits=4,colnames=FALSE,align="l",{
		x<-refresh_caseLR()	
		LRs<-x$LRs
		CLR<-LRs$CLR
		priorProb<-input$prior_prob
		
		support<-ifelse(CLR>1,"H1 true","H2 true")
		logCLR<-abs(log10(CLR))
		if (logCLR==0) LRweight<-"neither"
		if(logCLR>0)LRweight<-"weak or limited"
		if(logCLR>1)LRweight<-"moderate"
		if(logCLR>2)LRweight<-"strong"
		if(logCLR>3)LRweight<-"very strong"
		if(logCLR>4)LRweight<-"extremely strong"
		posteriorProb<-CLR/(CLR+((1-priorProb)/priorProb))
		data.frame(cbind(items=c("Support ","Weight of evidence ","Combined likelihood ratio (CLR)","log10(CLR)","Posterior probability "),Result=c(support,LRweight,CLR,log10(CLR),posteriorProb)))
	})
	#Open Easykin manual
	output$sraManual <- renderUI({tags$iframe(style = "height: 800px; width: 100%", src = "Easykin manual.pdf")})
	output$sraManual_Chinese <- renderUI({tags$iframe(style = "height: 800px; width: 100%", src = "Chinese manual.pdf")})
	
	#logs
	# all inputs that have changed
	output$all_inputs <- renderPrint({
	  input$.shinylogs_input
	})
	# errors
	output$log_errors <- renderPrint({
	  all_erros <- c()
	  x <- input$.shinylogs_error
	  y <- x$errors
	  if (length(y) > 0){
	   for (i in 1:length(y))if (y[[i]]$error != "")  print(paste(y[[i]][[1]],y[[i]][[2]],y[[i]][[3]],sep="     "))
	  }
	  
	  })
	
}
shiny::devmode(TRUE)
##run app
shinyApp(ui=ui, server=server)



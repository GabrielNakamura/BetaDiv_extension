#####function used to simulate metacommunity accordingly with a process of phylogenetic habitat filtering, adapted from Debastiani's function
##### original function https://github.com/vanderleidebastiani/sim.comm

#####description####
#tree= phylogenetic tree used to simulate metacommunities
#limit= lower limit of the phylogenetic filter used
#Ncomm= number of communities in metacommunity 
#Nsp= minimmunn numbeer of species in communities
#richness.equal= Logical, if TRUE all communities will present the same number of species
#filter= numeric, a single value between 0 and 0.99 that indicate the degree of phylogenetic filter in communities
#gradual= Logical, if TRUE the phylogenetic filter will be applied gradually, from the value set in filter 
          #to the value set in limit 
sim.comm_2<-function(tree,limit=0,Ncomm=30,Nsp=5,richness.equal=TRUE,filter=0.9,gradual=TRUE){
  Nsppool=length(tree$tip.label)
  disttree<-cophenetic(tree)
  simtree<-1-(disttree/max(disttree))
  gradient1<-rep(filter,Ncomm)
  if(gradual==TRUE){
    gradient1<-seq(limit,filter,length.out=Ncomm)
  }
  ResSim<-matrix(0,Ncomm,Nsppool)
  colnames(ResSim)=tree$tip.label
  rownames(ResSim)=sprintf("Comm_%.3d",1:dim(ResSim)[1])
  for (i in 1:Ncomm){
    Sp_Ref<-sample(tree$tip.label,1,replace=TRUE)
    SimilarRef<-simtree[which(colnames(simtree)==Sp_Ref),]
    QuantileSimilar<-quantile(SimilarRef,probs=gradient1[i])
    select<-unique(names(SimilarRef[SimilarRef>=QuantileSimilar]))
    SELECT<-select
    while(length(SELECT)<Nsp){
      select<-unique(names(SimilarRef[SimilarRef>=QuantileSimilar]))
      SELECT<-select
    }
    if(richness.equal==TRUE){
      if(length(SELECT)>Nsp){
        SELECT<-SELECT[sample(1:length(SELECT),size = Nsp,replace = FALSE)]
      }
    }
    ResSim[i,match(SELECT,names(ResSim[i,]))]<-1
    }
  return(ResSim)
}


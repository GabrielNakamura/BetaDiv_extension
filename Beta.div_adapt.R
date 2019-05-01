Beta.div_adapt<- function(Y, tree, trait,nperm=999, method= "distance.based"){
  if(dim(Y)[1]<=1){
    stop("\n matrix Y must be at least three communities\n")
  }
  if(sum(is.na(match(colnames(Y),tree$tip.label)))>0){
    stop("\n there are species in community matrix no present in phylogeny\n")
  }
  if(is.matrix(Y)==FALSE){
    as.matrix(Y)
  }
  P<-SYNCSA::matrix.p(comm = Y, phylodist= cophenetic(tree),notification = FALSE)$matrix.P
  X<-SYNCSA::matrix.x(comm = Y,traits = trait,scale = FALSE,notification = FALSE)$matrix.X
  if(method=="raw.data"){
  	P<- vegan::decostand(P,method = "normalize")
  	X<- vegan::decostand(X, method= "normalize")
  	BDp_all<- beta.div(Y = P, method = "none", sqrt.D = FALSE, nperm = 2)
  	BDx_all<- beta.div(Y = X, method = "none", sqrt.D = FALSE, nperm= 2)
  	BDp<- BDp_all$SStotal_BDtotal[2]
  	BDx<- BDx_all$SStotal_BDtotal[2]
  	PLCBD<- BDp_all$LCBD
  	XLCBD<- BDx_all$LCBD
  	PSCBD<- BDp_all$SCBD
  	XSCBD<- BDx_all$SCBD
  }
  BDp_all<- beta.div(Y = P, method = "percentagedifference", sqrt.D = TRUE, nperm = 2)
  BDx_all<- beta.div(Y = X, method = "percentagedifference", sqrt.D = TRUE, nperm = 2)
  BDp<- BDp_all$SStotal_BDtotal[2]
  BDx<- BDx_all$SStotal_BDtotal[2]
  PLCBD<- BDp_all$LCBD
  XLCBD<- BDx_all$LCBD
   BD_allNull<- matrix(NA, nrow = nperm, ncol= 2, dimnames= list(paste("run",1:nperm, sep = ""), c("BDp_null", "BDx_null")))
   PLCBD_allNull<- matrix(NA, nrow = nperm, ncol= ncol(Y), dimnames= list(paste("run", 1:nperm, sep=""), colnames(Y)))
   XLCBD_allNull<- matrix(NA, nrow = nperm, ncol= ncol(Y), dimnames= list(paste("run", 1:nperm, sep=""), colnames(Y)))
   for(i in 1:nperm) {
     distspp_null <- picante::taxaShuffle(dist.spp) #taxa shuffle matrix P
     distrait_null<- picante::taxaShuffle(dist.trait) #taxa shuffle matrix X
     match.namesP <- match(colnames(comm), colnames(distspp_null)) 
     match.namesX<- match(colnames(comm), colnames(distrait_null))
     matrixP_null<- SYNCSA::matrix.p(ifelse(comm>=1,1,0),as.matrix(distspp_null[match.namesP,match.namesP]))$matrix.P #matrix P taxa shuffle
     matrixX_null<- SYNCSA::matrix.p(ifelse(comm>=1,1,0),as.matrix(distrait_null[match.namesX,match.namesX]))$matrix.P #matrix X taxa shuffle
     match.namesP <- match(colnames(Y), colnames(distspp_null))
     match.namesX<- match(colnames(Y), colnames(distrait_null))
     P.null<-SYNCSA::matrix.p(comm = Y,phylodist = distspp_null[match.namesP,match.namesP],notification = FALSE)$matrix.P
     X.perm<-SYNCSA::matrix.p(comm = Y,phylodist = distrait_null[match.namesX,match.namesX],notification = FALSE)$matrix.P
     if(method="raw.data"){
       P.null<- vegan::decostand(P.null, method="normalize")
       X.null<- vegan::decostand(X.null, method="normalize")
     }
     BDp_allNull<- beta.div(Y = P.null, method = "percentagedifference", sqrt.D = TRUE, samp = 2)
     BDx_allNull<- beta.div(Y = X.null, method = "percentagedifference", sqrt.D = TRUE, samp = 2)
     BD_allNull[i,1]<- BDp_allNull$SStotal_BDtotal[2]
     BD_allNull[i,2]<- BDx_allNull$SStotal_BDtotal[2]
     PLCBD_allNull[i,]<- BDp_allNull$LCBD
     XLCBD_allNull[i,]<- BDx_allNull$LCBD
   }
   p.BDp<- (sum(ifelse(BD_allNull[,"BDp_null"] >= BDp, 1, 0)) + 1)/(nperm + 1)
   p.BDx<- (sum(ifelse(BD_allNull[,"BDx_null"] >= BDx, 1, 0)) + 1)/(nperm + 1)
   p.PLCBD<- numeric(length = ncol(Y))
   p.XLCBD<- numeric(length = ncol(Y))
   for(j in 1:ncol(Y)){
     p.PLCBD[j]<- (sum(ifelse(PLCBD_allNull[,j] >= PLCBD, 1, 0)) + 1)/(nperm + 1)
     p.XLCBD[j]<- (sum(ifelse(XLCBD_allNull[,j] >= XLCBD, 1, 0)) + 1)/(nperm + 1)
   }
   return(list(BD.obs= c(BDp, BDx), LCBD.obs= matrix(c(PLCBD, XLCBD), nrow= 2, ncol= ncol(Y), byrow= TRUE, dimnames= list(c("PLCBD","XLCBD"), colnames(Y))),
               p.BD= c(p.BDp, p.BDx), p.LCBD= matrix(c(PLCBD, XLCBD), nrow= 2, ncol= ncol(Y), byrow= TRUE, dimnames= list(c("p.PLCBD","p.XLCBD"), colnames(Y)))))
}

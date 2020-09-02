###########AIC table for scenario 1 - Power ###################################
library(apTreeshape)
source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\beta-div Legendre and Caceres 2013\\beta-div.R",local = TRUE)
source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\simul_metacomm\\simulate_metacommunity_MEE.R",local = TRUE)
#mac
#source("/Users/Gabriel/Google Drive/Fun????es/beta-div Legendre and Caceres 2013/beta-div.R", local=TRUE)
#source("/Users/Gabriel/Google Drive/Fun????es/simul_metacomm/simulate_metacommunity_MEE.R")
library("apTreeshape")
power<-1
nspp<-200
ncomm<-50
nperm<-2
runs<-999
ntree<- 1000
multiphylo<-vector(mode="list",length=999) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}

comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity(multiphylo[[i]],n.comm=ncomm,power=power,u=5)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}
matrix.p.simul<-vector(mode="list", length=length(multiphylo))
matrix.x.simul<-vector(mode="list",length=length(multiphylo))
matrix.CWM.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  matrix.p.simul[[i]] <-SYNCSA::matrix.p(ifelse(comm.simul[[i]]>=1,1,0),cophenetic(multiphylo[[i]]))$matrix.P
  matrix.x.simul[[i]]<-SYNCSA::matrix.p(ifelse(comm.simul[[i]]>=1,1,0),as.matrix(dist(trait.simul[[i]])))$matrix.P
}
LCBD_simul1P<-vector(mode="list",length=length(multiphylo))
LCBD_simul1W<-vector(mode="list",length=length(multiphylo))
LCBD_simul1X<-vector(mode="list",length=length(multiphylo))
SCBD_simul1P<-vector(mode="list",length=length(multiphylo))
SCBD_simul1W<-vector(mode="list",length=length(multiphylo))
SCBD_simul1X<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(matrix.p.simul)){
  LCBDP<-beta.div(matrix.p.simul[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDW<-beta.div(comm.simul[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDX<-beta.div(matrix.x.simul[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBD_simul1P[[i]]<-LCBDP$LCBD
  LCBD_simul1W[[i]]<-LCBDW$LCBD
  LCBD_simul1X[[i]]<-LCBDX$LCBD
}
rowSums(matrix.p.simul[[i]])
beta.div(matrix.p.simul[[1]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
#running linear models with LCBD values generated with scenario 1
glm.simul1P<-vector(mode="list", length=length(multiphylo))
glm.simul1W<-vector(mode="list", length=length(multiphylo))
glm.simul1X<-vector(mode="list", length=length(multiphylo))
parameters_model1<-matrix(NA, nrow=length(multiphylo),ncol=6, dimnames=list(paste("tree",1:length(multiphylo),sep=""),
                                                                            c("glmP","R2_P","glmW","R2_W","glmX","R2_X")))
AIC_table1<-matrix(NA,nrow=length(multiphylo),ncol=6,dimnames = list(paste("tree",1:length(multiphylo),sep=""),c(
  "AICP_poli","AICP_spoli","AICW_poli","AICW_spoli","AICX_poli","AICX_spoli")),byrow = TRUE)
for(i in 1:length(env.simul)){
  m1P_poli<-glm(LCBD_simul1P[[i]]~poly(env.simul[[i]],2))
  m1P_spoli<-glm(LCBD_simul1P[[i]]~env.simul[[i]])
  m1W_poli<-glm(LCBD_simul1W[[i]]~poly(env.simul[[i]],2))
  m1W_spoli<-glm(LCBD_simul1W[[i]]~env.simul[[i]])
  m1X_poli<-glm(LCBD_simul1X[[i]]~poly(env.simul[[i]],2))
  m1X_spoli<-glm(LCBD_simul1X[[i]]~env.simul[[i]])
  AIC_table1[i,c(1,2)]<-AIC(m1P_poli,m1P_spoli)[,2]
  AIC_table1[i,c(3,4)]<-AIC(m1W_poli,m1W_spoli)[,2]
  AIC_table1[i,c(5,6)]<-AIC(m1X_poli,m1X_spoli)[,2]
  parameters_model1[i,c(1,2)]<-c(summary.lm(m1P_poli)[5]$coefficients[3,4],summary.lm(m1P_poli)[10]$adj.r.squared)
  parameters_model1[i,c(3,4)]<-c(summary.lm(m1W_poli)[5]$coefficients[3,4],summary.lm(m1W_poli)[10]$adj.r.squared)
  parameters_model1[i,c(5,6)]<-c(summary.lm(m1X_poli)[5]$coefficients[3,4],summary.lm(m1X_poli)[10]$adj.r.squared)
}
windows()
apply(AIC_table1,2,mean) #mean values of AIC comparing models with polynomial terms and without polinomial term
apply(parameters_model1,2,mean) #mean power and R2 of polinomial models generate with scenario 1
apply(parameters_model1,2,sd)
length(which(AIC_table1[,1]<AIC_table1[,2]))/999
length(which(AIC_table1[,3]<AIC_table1[,4]))/999
length(which(AIC_table1[,5]<AIC_table1[,6]))/999
###########AIC table for scenario 2 ###################################
library(apTreeshape)
source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\beta-div Legendre and Caceres 2013\\beta-div.R",local = TRUE)
source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\simul_metacomm\\simulate_metacommunity_adapt.R",local = TRUE)
#mac
#source("/Users/Gabriel/Google Drive/Fun????es/beta-div Legendre and Caceres 2013/beta-div.R", local=TRUE)
#source("/Users/Gabriel/Google Drive/Fun????es/simul_metacomm/simulate_metacommunity_MEE.R")
library("apTreeshape")
power<-1
nspp<-200
ncomm<-50
nperm<-2
runs<-999
ntree<- 1000
multiphylo2<-vector(mode="list",length=999) 
for(i in 1:length(multiphylo)){
  multiphylo2[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}

comm.simul2<-vector(mode="list",length=length(multiphylo2))
trait.simul2<-vector(mode="list",length=length(multiphylo2))
env.simul2<-vector(mode="list",length=length(multiphylo2))
for(i in 1:length(multiphylo2)){
  simul<-simulate.metacommunity2(multiphylo2[[i]],n.comm=ncomm,power=power,u=5)
  comm.simul2[[i]]<-simul$W
  trait.simul2[[i]]<-simul$X
  env.simul2[[i]]<-simul$E
}
matrix.p.simul2<-vector(mode="list", length=length(multiphylo2))
matrix.x.simul2<-vector(mode="list",length=length(multiphylo2))
for(i in 1:length(multiphylo2)){
  matrix.p.simul2[[i]] <-SYNCSA::matrix.p(ifelse(comm.simul2[[i]]>=1,1,0),cophenetic(multiphylo2[[i]]))$matrix.P
  matrix.x.simul2[[i]]<-SYNCSA::matrix.p(ifelse(comm.simul2[[i]]>=1,1,0),as.matrix(dist(trait.simul2[[i]])))$matrix.P
}
LCBD_simul2P<-vector(mode="list",length=length(multiphylo2))
LCBD_simul2W<-vector(mode="list",length=length(multiphylo2))
LCBD_simul2X<-vector(mode="list",length=length(multiphylo))
SCBD_simul2P<-vector(mode="list",length=length(multiphylo))
SCBD_simul2W<-vector(mode="list",length=length(multiphylo))
SCBD_simul2X<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(matrix.p.simul2)){
  LCBDP2<-beta.div(matrix.p.simul2[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDW2<-beta.div(comm.simul2[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDX2<-beta.div(matrix.x.simul2[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBD_simul2P[[i]]<-LCBDP2$LCBD
  LCBD_simul2W[[i]]<-LCBDW2$LCBD
  LCBD_simul2X[[i]]<-LCBDX2$LCBD
}
beta.div(matrix.p.simul2[[1]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
#running linear models with LCBD values generated with scenario 1
glm.simul2P<-vector(mode="list", length=length(multiphylo2))
glm.simul2W<-vector(mode="list", length=length(multiphylo2))
glm.simul2X<-vector(mode="list", length=length(multiphylo2))
parameters_model2<-matrix(NA, nrow=length(multiphylo2),ncol=6, dimnames=list(paste("tree",1:length(multiphylo2),sep=""),
                                                                            c("glmP","R2_P","glmW","R2_W","glmX","R2_X")))
AIC_table2<-matrix(NA,nrow=length(multiphylo2),ncol=6,dimnames = list(paste("tree",1:length(multiphylo2),sep=""),c(
  "AICP_poli","AICP_spoli","AICW_poli","AICW_spoli","AICX_poli","AICX_spoli")),byrow = TRUE)
for(i in 1:length(env.simul2)){
  m1P_poli<-glm(LCBD_simul2P[[i]]~poly(env.simul2[[i]],2))
  m1P_spoli<-glm(LCBD_simul2P[[i]]~env.simul2[[i]])
  m1W_poli<-glm(LCBD_simul2W[[i]]~poly(env.simul2[[i]],2))
  m1W_spoli<-glm(LCBD_simul2W[[i]]~env.simul2[[i]])
  m1X_poli<-glm(LCBD_simul2X[[i]]~poly(env.simul2[[i]],2))
  m1X_spoli<-glm(LCBD_simul2X[[i]]~env.simul2[[i]])
  AIC_table2[i,c(1,2)]<-AIC(m1P_poli,m1P_spoli)[,2]
  AIC_table2[i,c(3,4)]<-AIC(m1W_poli,m1W_spoli)[,2]
  AIC_table2[i,c(5,6)]<-AIC(m1X_poli,m1X_spoli)[,2]
  parameters_model2[i,c(1,2)]<-c(summary.lm(m1P_poli)[5]$coefficients[3,4],summary.lm(m1P_poli)[10]$adj.r.squared)
  parameters_model2[i,c(3,4)]<-c(summary.lm(m1W_poli)[5]$coefficients[3,4],summary.lm(m1W_poli)[10]$adj.r.squared)
  parameters_model2[i,c(5,6)]<-c(summary.lm(m1X_poli)[5]$coefficients[3,4],summary.lm(m1X_poli)[10]$adj.r.squared)
}
windows()
apply(AIC_table2,2,mean) #mean values of AIC comparing models with polynomial terms and without polinomial term
length(which(AIC_table2[,1]<AIC_table2[,2]))/999
length(which(AIC_table2[,3]<AIC_table2[,4]))/999
length(which(AIC_table2[,5]<AIC_table2[,6]))/999
###########AIC table for scenario 3 ###################################
library(apTreeshape)
source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\beta-div Legendre and Caceres 2013\\beta-div.R",local = TRUE)
source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\simul_metacomm\\simulate_metacommunity_MEE.R",local = TRUE)
#mac
#source("/Users/Gabriel/Google Drive/Fun????es/beta-div Legendre and Caceres 2013/beta-div.R", local=TRUE)
#source("/Users/Gabriel/Google Drive/Fun????es/simul_metacomm/simulate_metacommunity_MEE.R")
library("apTreeshape")
power<-0.0001
nspp<-200
ncomm<-50
nperm<-2
runs<-999
ntree<- 1000
multiphylo3<-vector(mode="list",length=999) 
for(i in 1:length(multiphylo3)){
  multiphylo3[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}

comm.simul3<-vector(mode="list",length=length(multiphylo3))
trait.simul3<-vector(mode="list",length=length(multiphylo3))
env.simul3<-vector(mode="list",length=length(multiphylo3))
for(i in 1:length(multiphylo3)){
  simul<-simulate.metacommunity(multiphylo3[[i]],n.comm=ncomm,power=power,u=5)
  comm.simul3[[i]]<-simul$W
  trait.simul3[[i]]<-simul$X
  env.simul3[[i]]<-simul$E
}
matrix.p.simul3<-vector(mode="list", length=length(multiphylo3))
matrix.x.simul3<-vector(mode="list",length=length(multiphylo3))
for(i in 1:length(multiphylo3)){
  matrix.p.simul3[[i]] <-SYNCSA::matrix.p(ifelse(comm.simul3[[i]]>=1,1,0),cophenetic(multiphylo3[[i]]))$matrix.P
  matrix.x.simul3[[i]]<-SYNCSA::matrix.p(ifelse(comm.simul3[[i]]>=1,1,0),as.matrix(dist(trait.simul3[[i]])))$matrix.P
}
LCBD_simul3P<-vector(mode="list",length=length(multiphylo3))
LCBD_simul3W<-vector(mode="list",length=length(multiphylo2))
LCBD_simul3X<-vector(mode="list",length=length(multiphylo))
SCBD_simul3P<-vector(mode="list",length=length(multiphylo))
SCBD_simul3W<-vector(mode="list",length=length(multiphylo))
SCBD_simul3X<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(matrix.p.simul3)){
  LCBDP3<-beta.div(matrix.p.simul3[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDW3<-beta.div(comm.simul3[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDX3<-beta.div(matrix.x.simul3[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBD_simul3P[[i]]<-LCBDP3$LCBD
  LCBD_simul3W[[i]]<-LCBDW3$LCBD
  LCBD_simul3X[[i]]<-LCBDX3$LCBD
}
beta.div(matrix.p.simul3[[1]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
#running linear models with LCBD values generated with scenario 1
glm.simul3P<-vector(mode="list", length=length(multiphylo3))
glm.simul3W<-vector(mode="list", length=length(multiphylo3))
glm.simul3X<-vector(mode="list", length=length(multiphylo3))
parameters_model3<-matrix(NA, nrow=length(multiphylo3),ncol=6, dimnames=list(paste("tree",1:length(multiphylo3),sep=""),
                                                                             c("glmP","R2_P","glmW","R2_W","glmX","R2_X")))
AIC_table3<-matrix(NA,nrow=length(multiphylo2),ncol=6,dimnames = list(paste("tree",1:length(multiphylo2),sep=""),c(
  "AICP_poli","AICP_spoli","AICW_poli","AICW_spoli","AICX_poli","AICX_spoli")),byrow = TRUE)
for(i in 1:length(env.simul3)){
  m3P_poli<-glm(LCBD_simul3P[[i]]~poly(env.simul3[[i]],2))
  m3P_spoli<-glm(LCBD_simul3P[[i]]~env.simul3[[i]])
  m3W_poli<-glm(LCBD_simul3W[[i]]~poly(env.simul3[[i]],2))
  m3W_spoli<-glm(LCBD_simul3W[[i]]~env.simul3[[i]])
  m3X_poli<-glm(LCBD_simul3X[[i]]~poly(env.simul3[[i]],2))
  m3X_spoli<-glm(LCBD_simul3X[[i]]~env.simul3[[i]])
  AIC_table3[i,c(1,2)]<-AIC(m3P_poli,m3P_spoli)[,2]
  AIC_table3[i,c(3,4)]<-AIC(m3W_poli,m3W_spoli)[,2]
  AIC_table3[i,c(5,6)]<-AIC(m3X_poli,m3X_spoli)[,2]
  parameters_model3[i,c(1,2)]<-c(summary.lm(m3P_poli)[5]$coefficients[3,4],summary.lm(m3P_poli)[10]$adj.r.squared)
  parameters_model3[i,c(3,4)]<-c(summary.lm(m3W_poli)[5]$coefficients[3,4],summary.lm(m3W_poli)[10]$adj.r.squared)
  parameters_model3[i,c(5,6)]<-c(summary.lm(m3X_poli)[5]$coefficients[3,4],summary.lm(m3X_poli)[10]$adj.r.squared)
}
windows()
apply(AIC_table3,2,mean) #mean values of AIC comparing models with polynomial terms and without polinomial term
length(which(AIC_table3[,1]<AIC_table3[,2]))/999
length(which(AIC_table3[,3]<AIC_table3[,4]))/999
length(which(AIC_table3[,5]<AIC_table3[,6]))/999

###########AIC table for scenario 4 - Type I Error ###################################

source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\beta-div Legendre and Caceres 2013\\beta-div.R",local = TRUE)
source("C:\\Users\\Gabriel\\Google Drive\\Fun??es\\simul_metacomm\\simulate_metacommunity_MEE.R",local = TRUE)
#mac
source("/Users/Gabriel/Google Drive/Fun????es/beta-div Legendre and Caceres 2013/beta-div.R", local=TRUE)
source("/Users/Gabriel/Google Drive/Fun????es/simul_metacomm/simulate_metacommunity_MEE.R")

power<-1
nspp<-200
ncomm<-50
nperm<-2
runs<-999
ntree<- 1000
multiphylo4<-vector(mode="list",length=999) 
for(i in 1:length(multiphylo)){
  multiphylo4[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}

comm.simul4<-vector(mode="list",length=length(multiphylo4))
trait.simul4<-vector(mode="list",length=length(multiphylo4))
env.simul4<-vector(mode="list",length=length(multiphylo4))
for(i in 1:length(multiphylo4)){
  simul4<-simulate.metacommunity(multiphylo4[[i]],n.comm=ncomm,power=power,u=5,new.E = TRUE)
  comm.simul4[[i]]<-simul$W
  trait.simul4[[i]]<-simul$X
  env.simul4[[i]]<-simul$E
}
matrix.p.simul4<-vector(mode="list", length=length(multiphylo4))
matrix.x.simul4<-vector(mode="list",length=length(multiphylo4))
matrix.CWM.simul4<-vector(mode="list",length=length(multiphylo4))
for(i in 1:length(multiphylo4)){
  matrix.p.simul4[[i]] <-SYNCSA::matrix.p(ifelse(comm.simul4[[i]]>=1,1,0),cophenetic(multiphylo4[[i]]))$matrix.P
  matrix.x.simul4[[i]]<-SYNCSA::matrix.p(ifelse(comm.simul4[[i]]>=1,1,0),as.matrix(dist(trait.simul4[[i]])))$matrix.P
}
LCBD_simul4P<-vector(mode="list",length=length(multiphylo4))
LCBD_simul4W<-vector(mode="list",length=length(multiphylo4))
LCBD_simul4X<-vector(mode="list",length=length(multiphylo4))
SCBD_simul4P<-vector(mode="list",length=length(multiphylo4))
SCBD_simul4W<-vector(mode="list",length=length(multiphylo4))
SCBD_simul4X<-vector(mode="list",length=length(multiphylo4))
for(i in 1:length(matrix.p.simul4)){
  LCBDP<-beta.div(matrix.p.simul4[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDW<-beta.div(comm.simul4[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBDX<-beta.div(matrix.x.simul4[[i]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
  LCBD_simul4P[[i]]<-LCBDP$LCBD
  LCBD_simul4W[[i]]<-LCBDW$LCBD
  LCBD_simul4X[[i]]<-LCBDX$LCBD
}
rowSums(matrix.p.simul[[i]])
beta.div(matrix.p.simul[[1]],method = "percentagedifference",sqrt.D = TRUE,nperm= nperm)
#running linear models with LCBD values generated with scenario 1
glm.simul4P<-vector(mode="list", length=length(multiphylo4))
glm.simul4W<-vector(mode="list", length=length(multiphylo4))
glm.simul4X<-vector(mode="list", length=length(multiphylo4))
parameters_model4<-matrix(NA, nrow=length(multiphylo4),ncol=6, dimnames=list(paste("tree",1:length(multiphylo4),sep=""),
                                                                            c("glmP","R2_P","glmW","R2_W","glmX","R2_X")))
AIC_table4<-matrix(NA,nrow=length(multiphylo4),ncol=6,dimnames = list(paste("tree",1:length(multiphylo4),sep=""),c(
  "AICP_poli","AICP_spoli","AICW_poli","AICW_spoli","AICX_poli","AICX_spoli")),byrow = TRUE)
for(i in 1:length(env.simul4)){
  m4P_poli<-glm(LCBD_simul4P[[i]]~poly(env.simul4[[i]],2))
  m4P_spoli<-glm(LCBD_simul4P[[i]]~env.simul4[[i]])
  m4W_poli<-glm(LCBD_simul4W[[i]]~poly(env.simul4[[i]],2))
  m4W_spoli<-glm(LCBD_simul4W[[i]]~env.simul4[[i]])
  m4X_poli<-glm(LCBD_simul4X[[i]]~poly(env.simul4[[i]],2))
  m4X_spoli<-glm(LCBD_simul4X[[i]]~env.simul4[[i]])
  AIC_table4[i,c(1,2)]<-AIC(m4P_poli,m4P_spoli)[,2]
  AIC_table4[i,c(3,4)]<-AIC(m4W_poli,m4W_spoli)[,2]
  AIC_table4[i,c(5,6)]<-AIC(m4X_poli,m4X_spoli)[,2]
  parameters_model4[i,c(1,2)]<-c(summary.lm(m4P_poli)[5]$coefficients[3,4],summary.lm(m4P_poli)[10]$adj.r.squared)
  parameters_model4[i,c(3,4)]<-c(summary.lm(m4W_poli)[5]$coefficients[3,4],summary.lm(m4W_poli)[10]$adj.r.squared)
  parameters_model4[i,c(5,6)]<-c(summary.lm(m4X_poli)[5]$coefficients[3,4],summary.lm(m4X_poli)[10]$adj.r.squared)
}
windows()
apply(AIC_table4,2,mean) #mean values of AIC comparing models with polynomial terms and without polinomial term
apply(parameters_model4,2,mean) #mean power and R2 of polinomial models generate with scenario 1
apply(parameters_model4,2,sd)
length(which(AIC_table4[,1]<AIC_table4[,2]))/999
length(which(AIC_table4[,3]<AIC_table4[,4]))/999
length(which(AIC_table4[,5]<AIC_table4[,6]))/999

##################Aic for null models########################
multiphylo
AIC_table_nullmod<-matrix(NA,nrow=length(multiphylo),ncol=12,dimnames = list(paste("tree",1:length(multiphylo),sep=""),c(
  "AICP_poli","AICP_spoli","AICW_poli","AICW_spoli","AICX_poli","AICX_spoli",
  "AICP_poli_site","AICP_spoli_site","AICW_poli_site","AICW_spoli_site","AICX_poli_site","AICX_spoli_site")),byrow = TRUE)
for (k in 1:length(multiphylo)) {
  distspp_null <- picante::taxaShuffle(cophenetic(multiphylo[[k]]))
  distrait_null<- picante::taxaShuffle(dist(trait.simul[[k]],diag = TRUE,upper = TRUE))
  match.namesP <- match(colnames(comm.simul[[k]]), colnames(distspp_null))
  match.namesX<- match(colnames(comm.simul[[k]]), colnames(distrait_null))
  matrixP_null<- SYNCSA::matrix.p(ifelse(comm.simul[[k]]>=1,1,0),as.matrix(distspp_null[match.namesP,match.namesP]))$matrix.P
  matrixX_null<- SYNCSA::matrix.p(ifelse(comm.simul[[k]]>=1,1,0),as.matrix(distrait_null[match.namesX,match.namesX]))$matrix.P
  LCBDP_null<- beta.div(matrixP_null,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm)$LCBD
  LCBDX_null<- beta.div(matrixX_null,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm)$LCBD
  LCBDW_null<-LCBDW$LCBD[sample(1:length(LCBDW$LCBD))]
  names(LCBDW_null)<-names(LCBDW)
  modP_null_taxa_poli <- glm(LCBDP_null~poly(env.simul[[k]],2))
  modP_null_taxa <- glm(LCBDP_null~env.simul[[k]])
  modX_null_taxa_poli<- glm(LCBDX_null~poly(env.simul[[k]],2))
  modX_null_taxa<- glm(LCBDX_null~env.simul[[k]])
  modW_null_taxa_poli<- glm(LCBDW_null~poly(env.simul[[k]],2))
  modW_null_taxa<- glm(LCBDW_null~env.simul[[k]])
  LCBDP_null_site <- LCBDP$LCBD[sample(1:length(LCBDP$LCBD))]
  LCBDX_null_site <- LCBDX$LCBD[sample(1:length(LCBDX$LCBD))]
  LCBDW_null_site <- LCBDW$LCBD[sample(1:length(LCBDW$LCBD))]
  modP_null_site_poli <- glm(LCBDP_null_site~poly(env.simul[[k]],2))
  modP_null_site <- glm(LCBDP_null_site~env.simul[[k]])
  modX_null_site_poli <- glm(LCBDX_null_site~poly(env.simul[[k]],2))
  modX_null_site <- glm(LCBDX_null_site~env.simul[[k]])
  modW_null_site_poli <- glm(LCBDW_null_site~poly(env.simul[[k]],2))
  modW_null_site <- glm(LCBDW_null_site~env.simul[[k]])
  AIC_table_nullmod[k,c(1,2)]<-AIC(modP_null_taxa_poli,modP_null_taxa)[,2]
  AIC_table_nullmod[k,c(3,4)]<-AIC(modX_null_taxa_poli,modX_null_taxa)[,2]
  AIC_table_nullmod[k,c(5,6)]<-AIC(modW_null_taxa_poli,modW_null_taxa)[,2]
  AIC_table_nullmod[k,c(7,8)]<-AIC(modP_null_site_poli,modP_null_site)[,2]
  AIC_table_nullmod[k,c(9,10)]<-AIC(modX_null_site_poli,modX_null_site)[,2]
  AIC_table_nullmod[k,c(11,12)]<-AIC(modW_null_site_poli,modW_null_site)[,2]
}
apply(AIC_table_nullmod,2,mean)
length(which(AIC_table_nullmod[,1]<AIC_table_nullmod[,2]))/999 #number of glm poli with AIC lower than glm without poli for PLCBD null taxa
length(which(AIC_table_nullmod[,3]<AIC_table_nullmod[,4]))/999 #number of glm poli with AIC lower than glm without poli for XLCBD null taxa
length(which(AIC_table_nullmod[,5]<AIC_table_nullmod[,6]))/999 #number of glm poli with AIC lower than glm without poli for WLCBD null taxa
length(which(AIC_table_nullmod[,7]<AIC_table_nullmod[,8]))/999 #number of glm poli with AIC lower than glm without poli for PLCBD null site
length(which(AIC_table_nullmod[,9]<AIC_table_nullmod[,10]))/999 #number of glm poli with AIC lower than glm without poli for XLCBD null site
length(which(AIC_table_nullmod[,11]<AIC_table_nullmod[,12]))/999 #number of glm poli with AIC lower than glm without poli for WLCBD null site


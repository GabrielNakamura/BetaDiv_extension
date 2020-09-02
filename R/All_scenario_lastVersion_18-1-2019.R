##################Scenario 5 - E, P0, X0##########
setwd("/Users/gabrielnakamura/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/analysis_performance/results_08-01-2019")
source("~/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/LCBDsig_Parfunction_18-01-2019.R",local=TRUE)
source("~/Google Drive/Funcoes/simul_metacomm/simulate_metacommunity_adapt.R",local=TRUE)
source("~/Google Drive/Funcoes/beta-div Legendre and Caceres 2013/beta-div.R",local=TRUE)

power<-0.0001
nspp<-200
ncomm<-50
nperm<-2
ntree<- 999
runs<- 999
multiphylo<-vector(mode="list",length=ntree) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}
comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity2(multiphylo[[i]],n.comm=ncomm,power=power,u=5,new.X = TRUE)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}
data_simul5<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  data_simul5[[i]]<-c(list(comm.simul[[i]]),list(cophenetic(multiphylo[[i]])),list(dist(trait.simul[[i]],diag = T,upper = T)),list(env.simul[[i]]))
}
library(parallel)
ncor<-detectCores()
CL<-makeCluster(ncor)
clusterExport(CL,c("LCBD.sigPar", "data_simul5","power","nspp","ncomm","runs","beta.div","nperm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
resLCBDsig_simul5<-parLapply(cl = CL,X = data_simul5,fun = function(i) LCBD.sigPar( i ,nperm = nperm,runs = runs,family = gaussian))
stopCluster(CL)
save.image("resLCBDsig_simul5_18-1-2019.RData")
rm(list=ls())

#extracting TI error rate and power for Scenario 5
matrix_resScenario5Raw<- matrix(unlist(lapply(resLCBDsig_simul5, function(i) i$Raw_results)), 
                                nrow= 999, ncol= ncol(resLCBDsig_simul5[[1]]$Raw_results), byrow = T, 
                                dimnames= list(paste("run", 1:999, sep= ""), 
                                               colnames(resLCBDsig_simul5[[1]]$Raw_results)
                                               )
                                )

matrix_resScenario5Dist<- matrix(unlist(lapply(resLCBDsig_simul5, function(i) i$Dist_results)), 
                                 nrow= 999, ncol= ncol(resLCBDsig_simul5[[1]]$Dist_results), byrow = T, 
                                 dimnames= list(paste("run", 1:999, sep= ""), 
                                                colnames(resLCBDsig_simul5[[1]]$Dist_results)
                                                )
                                 )

p_valuesS5_Raw<- apply(matrix_resScenario5Raw[, 7:15], 2, function(i) length(which(i<=0.05))/999)
p_valuesS5_Dist<- apply(matrix_resScenario5Dist[, 7:15], 2, function(i) length(which(i<=0.05))/999)
mean_valuesS5_Raw<- apply(matrix_resScenario5Raw[, c(1:6, 16:21)], 2, function(i) mean(i))
mean_valuesS5_Dist<- apply(matrix_resScenario5Dist[, c(1:6, 16:21)], 2, function(i) mean(i))

##################Scenario 6 - E, P, X0##########
setwd("/Users/gabrielnakamura/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/analysis_performance/results_08-01-2019")
source("~/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/LCBDsig_Parfunction_18-01-2019.R",local=TRUE)
source("~/Google Drive/Funcoes/simul_metacomm/simulate_metacommunity_adapt.R",local=TRUE)
source("~/Google Drive/Funcoes/beta-div Legendre and Caceres 2013/beta-div.R",local=TRUE)

power<-1
nspp<-200
ncomm<-50
nperm<-2
ntree<-999
runs<-999
multiphylo<-vector(mode="list",length=ntree) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}
comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity2(multiphylo[[i]],n.comm=ncomm,power=power,u=5,new.X = TRUE)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}
data_simul6<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  data_simul6[[i]]<-c(list(comm.simul[[i]]),list(cophenetic(multiphylo[[i]])),list(dist(trait.simul[[i]],diag = T,upper = T)),list(env.simul[[i]]))
}
library(parallel)
ncor<-detectCores()
CL<-makeCluster(ncor)
clusterExport(CL,c("LCBD.sigPar", "data_simul6","power","nspp","ncomm","runs","beta.div","nperm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
resLCBDsig_simul6<-parLapply(cl = CL,X = data_simul6,fun = function(i) LCBD.sigPar( i ,nperm = nperm,runs = runs,family = gaussian))
stopCluster(CL)
save.image("resLCBDsig_simul6_18-1-2019.RData")
rm(list=ls())

#extracting TI error rate and power for Scenario 6
matrix_resScenario6Raw<- matrix(unlist(lapply(resLCBDsig_simul6, function(i) i$Raw_results)), 
                                nrow= 999, ncol= ncol(resLCBDsig_simul6[[1]]$Raw_results), byrow = T, 
                                dimnames= list(paste("run", 1:999, sep= ""), 
                                               colnames(resLCBDsig_simul6[[1]]$Raw_results)
                                )
)

matrix_resScenario6Dist<- matrix(unlist(lapply(resLCBDsig_simul6, function(i) i$Dist_results)), 
                                 nrow= 999, ncol= ncol(resLCBDsig_simul6[[1]]$Dist_results), byrow = T, 
                                 dimnames= list(paste("run", 1:999, sep= ""), 
                                                colnames(resLCBDsig_simul6[[1]]$Dist_results)
                                 )
)

p_valuesS6_Raw<- apply(matrix_resScenario6Raw[, 7:15], 2, function(i) length(which(i<=0.05))/999)
p_valuesS6_Dist<- apply(matrix_resScenario6Dist[, 7:15], 2, function(i) length(which(i<=0.05))/999)
mean_valuesS6_Raw<- apply(matrix_resScenario6Raw[, c(1:6, 16:21)], 2, function(i) mean(i))
mean_valuesS6_Dist<- apply(matrix_resScenario6Dist[, c(1:6, 16:21)], 2, function(i) mean(i))


##################Scenario 7 -E, P0, X#################
setwd("/Users/gabrielnakamura/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/analysis_performance/results_08-01-201")
source("/Users/gabrielnakamura/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/LCBDsig_Parfunction_18-01-2019.R",local=TRUE)
source("/Users/gabrielnakamura/Google Drive/Funcoes/simul_metacomm/simulate_metacommunity_MEE.R",local=TRUE)
source("/Users/gabrielnakamura/Google Drive/Funcoes/beta-div Legendre and Caceres 2013/beta-div.R",local=TRUE)

power<-0.0001
nspp<-200
ncomm<-50
nperm<-2
ntree<-999
runs<-999
multiphylo<-vector(mode="list",length=ntree) 
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
data_simul7<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  data_simul7[[i]]<-c(list(comm.simul[[i]]),list(cophenetic(multiphylo[[i]])),list(dist(trait.simul[[i]],diag = T,upper = T)),list(env.simul[[i]]))
}
library(parallel)
ncor<-detectCores()
CL<-makeCluster(ncor)
clusterExport(CL,c("LCBD.sigPar", "data_simul7","power","nspp","ncomm","runs","beta.div","nperm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
resLCBDsig_simul7<-parLapply(cl = CL,X = data_simul7,fun = function(i) LCBD.sigPar( i ,nperm = nperm,runs = runs,family = gaussian))
stopCluster(CL)
save.image("resLCBDsig_simul7_18-1-2019.RData")
rm(list=ls())

#extracting TI error rate and power for Scenario 7
matrix_resScenario7Raw<- matrix(unlist(lapply(resLCBDsig_simul7, function(i) i$Raw_results)), 
                                nrow= 999, ncol= ncol(resLCBDsig_simul7[[1]]$Raw_results), byrow = T, 
                                dimnames= list(paste("run", 1:999, sep= ""), 
                                               colnames(resLCBDsig_simul7[[1]]$Raw_results)
                                )
)

matrix_resScenario7Dist<- matrix(unlist(lapply(resLCBDsig_simul7, function(i) i$Dist_results)), 
                                 nrow= 999, ncol= ncol(resLCBDsig_simul7[[1]]$Dist_results), byrow = T, 
                                 dimnames= list(paste("run", 1:999, sep= ""), 
                                                colnames(resLCBDsig_simul7[[1]]$Dist_results)
                                 )
)

p_valuesS7_Raw<- apply(matrix_resScenario7Raw[, 7:15], 2, function(i) length(which(i<=0.05))/999)
p_valuesS7_Dist<- apply(matrix_resScenario7Dist[, 7:15], 2, function(i) length(which(i<=0.05))/999)
mean_valuesS7_Raw<- apply(matrix_resScenario7Raw[, c(1:6, 16:21)], 2, function(i) mean(i))
mean_valuesS7_Dist<- apply(matrix_resScenario7Dist[, c(1:6, 16:21)], 2, function(i) mean(i))


##################Scenario 8 - E,P,X##########
setwd("/Users/gabrielnakamura/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/analysis_performance/results_08-01-2019")
source("~/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/LCBDsig_Parfunction_18-01-2019.R",local=TRUE)
source("~/Google Drive/Funcoes/simul_metacomm/simulate_metacommunity_MEE.R",local=TRUE)
source("~/Google Drive/Funcoes/beta-div Legendre and Caceres 2013/beta-div.R",local=TRUE)

power<-1
nspp<-200
ncomm<-50
nperm<-2 #number of permutations to calculate the significance of LCBD component
ntree<-999
runs<-999  #number of permutations for null BD and components

multiphylo<-vector(mode="list",length=ntree) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}
comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity(multiphylo[[i]],n.comm=ncomm,power=power,u=5,sd.E = 10,sd.X = 10,noise.E = TRUE,noise.X = TRUE)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}
data_simul8<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  data_simul8[[i]]<-c(list(comm=comm.simul[[i]]),dist.phy= list(cophenetic(multiphylo[[i]])),
                      dist.trait=list(as.matrix(dist(trait.simul[[i]],diag = T,upper = T))),env=list(env.simul[[i]]))
}

library(parallel)
ncor<-detectCores()
CL<-makeCluster(ncor)
clusterExport(CL,c("LCBD.sigPar", "data_simul8","power","nspp","ncomm","runs","beta.div","nperm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
resLCBDsig_simul8<-parLapply(cl = CL,X = data_simul8,fun = function(i) LCBD.sigPar( data = i ,nperm = nperm,runs = runs,family = gaussian))
stopCluster(CL)

save.image("resLCBDsig_simul8_18-1-2019.RData")
rm(list=ls())

#extracting TI error rate and power for Scenario 8
matrix_resScenario8Raw<- matrix(unlist(lapply(resLCBDsig_simul8, function(i) i$Raw_results)), 
                                nrow= 999, ncol= ncol(resLCBDsig_simul8[[1]]$Raw_results), byrow = T, 
                                dimnames= list(paste("run", 1:999, sep= ""), 
                                               colnames(resLCBDsig_simul8[[1]]$Raw_results)
                                )
)

matrix_resScenario8Dist<- matrix(unlist(lapply(resLCBDsig_simul8, function(i) i$Dist_results)), 
                                 nrow= 999, ncol= ncol(resLCBDsig_simul8[[1]]$Dist_results), byrow = T, 
                                 dimnames= list(paste("run", 1:999, sep= ""), 
                                                colnames(resLCBDsig_simul8[[1]]$Dist_results)
                                 )
)

p_valuesS8_Raw<- apply(matrix_resScenario8Raw[, 7:15], 2, function(i) length(which(i<=0.05))/999)
p_valuesS8_Dist<- apply(matrix_resScenario8Dist[, 7:15], 2, function(i) length(which(i<=0.05))/999)
mean_valuesS8_Raw<- apply(matrix_resScenario8Raw[, c(1:6, 16:21)], 2, function(i) mean(i))
mean_valuesS8_Dist<- apply(matrix_resScenario8Dist[, c(1:6, 16:21)], 2, function(i) mean(i))

#######W0,P1,X1###############
setwd("/Users/gabrielnakamura/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/analysis_performance/results_08-01-2019")
source("~/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/LCBDsig_Parfunction_18-01-2019.R",local=TRUE)
source("~/Google Drive/Funcoes/simul_metacomm/simulate_metacommunity_MEE.R",local=TRUE)
source("~/Google Drive/Funcoes/beta-div Legendre and Caceres 2013/beta-div.R",local=TRUE)

power<-1
nspp<-200
ncomm<-50
nperm<-2 #number of permutations to calculate the significance of LCBD component
ntree<- 999
runs<- 999  #number of permutations for null BD and components
nicheBreath<- 50 #higher niche breath, less Beta diversity

multiphylo<-vector(mode="list",length=ntree) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}
comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity(multiphylo[[i]],n.comm=ncomm,power=power,u= nicheBreath, noise.E = FALSE,noise.X = FALSE)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}
data_simul8_nullW<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  data_simul8_nullW[[i]]<-c(list(comm=comm.simul[[i]]),dist.phy= list(cophenetic(multiphylo[[i]])),
                            dist.trait=list(as.matrix(dist(trait.simul[[i]],diag = T,upper = T))),env=list(env.simul[[i]]))
}

library(parallel)
ncor<-detectCores()
CL<-makeCluster(ncor)
clusterExport(CL,c("LCBD.sigPar", "data_simul8_nullW","power","nspp","ncomm","runs","beta.div","nperm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
resLCBDsig_simul8_nullW<-parLapply(cl = CL,X = data_simul8_nullW,fun = function(i) LCBD.sigPar( data = i ,nperm = nperm,runs = runs,family = gaussian))
stopCluster(CL)
save.image("resLCBDsig_simul8_NullW_18-1-2019.RData")
rm(list=ls())

#######W0,P0,X1###############
setwd("/Users/gabrielnakamura/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/analysis_performance/results_08-01-2019")
source("~/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/simulation/LCBDsig_Parfunction_18-01-2019.R",local=TRUE)
source("~/Google Drive/Funcoes/simul_metacomm/simulate_metacommunity_MEE.R",local=TRUE)
source("~/Google Drive/Funcoes/beta-div Legendre and Caceres 2013/beta-div.R",local=TRUE)

power<- 0.0001
nspp<-200
ncomm<-50
nperm<-2 #number of permutations to calculate the significance of LCBD component
ntree<- 999
runs<- 999  #number of permutations for null BD and components
nicheBreath<- 50 #higher niche breath, less Beta diversity

multiphylo<-vector(mode="list",length=ntree) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}
comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity(multiphylo[[i]],n.comm=ncomm,power=power,u= nicheBreath, noise.E = FALSE,noise.X = FALSE)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}
data_simul8_nullW_P0<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  data_simul8_nullW_P0[[i]]<-c(list(comm=comm.simul[[i]]),dist.phy= list(cophenetic(multiphylo[[i]])),
                            dist.trait=list(as.matrix(dist(trait.simul[[i]],diag = T,upper = T))),env=list(env.simul[[i]]))
}

library(parallel)
ncor<-detectCores()
CL<-makeCluster(ncor)
clusterExport(CL,c("LCBD.sigPar", "data_simul8_nullW_P0","power","nspp","ncomm","runs","beta.div","nperm"))
clusterEvalQ(CL,library(vegan))
clusterEvalQ(CL,library(ape))
clusterEvalQ(CL,library(SYNCSA))
clusterEvalQ(CL,library(PCPS))
resLCBDsig_simul8_nullW_P0<-parLapply(cl = CL,X = data_simul8_nullW_P0,fun = function(i) LCBD.sigPar( data = i ,nperm = nperm,runs = runs,family = gaussian))
stopCluster(CL)
save.image("resLCBDsig_simul8_NullW_P0_18-1-2019.RData")
rm(list=ls())



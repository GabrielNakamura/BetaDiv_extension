####testing for dissimilarity properties through criteria P3 to P6 in Legendre and De Cáceres#######

####libraries####
library(ape)
library(picante)
library(SYNCSA)
library(here)
library(adespatial)
source(here("functions", "BetaDiv_extension_22-10-2019.R"))
source(here("functions", "simul_dimensionality.R"))

###data generation#####

#assembling phylogenetic tree
tree_toy2<- "((((A:3,B:3):2,(C:3,D:3):2):5,((E:3,F:3):2,(G:3,H:3):2):5):5):10;"
tree_toy2<- read.tree(text = tree_toy2)

#community data complete
W<- matrix(c(3,7,0,0,0,0,0,0,
             7,3,0,0,0,0,0,0,
             1,1,1,1,0,0,0,0,
             0,0,0,0,1,1,1,1),
           nrow= 4, ncol= length(tree_toy2$tip.label), byrow = T, 
           dimnames= list(c("S1", "S2", "S3", "S4"), 
                                     tree_toy2$tip.label
                                     )
           )

####testing P3 - Monotonicity in changes in abbundance####

#species phylogenetic similar
WP3_eqqSpp<- W[1:2,1:2]
WP3_eqqSpp[,]<- c(100, 50, 0, 10) 
increase_abund<- seq(10, 120, 10)
resP3_eqqSpp<- numeric(length= length(increase_abund))
resP3_eqqSpp_raw<- numeric(length= length(increase_abund))
resP3_noPhy<- numeric(length= length(increase_abund))
plot(tree_toy2)
for(i in 1:length(increase_abund)){
  WP3_eqqSpp[,]<- c(100, 50, 0, increase_abund[i]) 
  resP3_eqqSpp[i]<- Beta.div_adapt(Y = WP3_eqqSpp, 
                                   dist_spp = organize.syncsa(comm = WP3_eqqSpp, phylodist = cophenetic(tree_toy2))$phylodist, 
                                   nperm = 1, method = "dist")$BDextend.obs
  resP3_eqqSpp_raw[i]<- Beta.div_adapt(Y = WP3_eqqSpp, 
                                       dist_spp = organize.syncsa(comm = WP3_eqqSpp, phylodist = cophenetic(tree_toy2))$phylodist, 
                                       nperm = 1, method = "raw")$BDextend.obs
  resP3_noPhy[i]<- beta.div(Y = WP3_eqqSpp, method = "percentdiff", sqrt.D = T, nperm = 1)$beta[2]
}
quartz()
windows()
plot(increase_abund, resP3_eqqSpp, type = "b", pch= 19, ylim= c(0.0, 0.65))
points(increase_abund, resP3_eqqSpp_raw, type = "b", pch= 2)
points(increase_abund, resP3_noPhy, type = "b", pch= 3)
legend("topleft", legend = c("sqrt Bray", "raw", "original"), pch = c(19,2, 3), lty = 1)

#species phylogenetically dissimilar

WP3_diffSpp<- W[1:2,c(1,8)]
increase_abund<- seq(10, 120, 10)
resP3_diffSpp<- numeric(length= length(increase_abund))
resP3_diffSpp_raw<- numeric(length= length(increase_abund))
resP3_noPhy<- numeric(length= length(increase_abund))
for(i in 1:length(increase_abund)){
  WP3_diffSpp[,]<- c(100, 50, 0, increase_abund[i]) 
  resP3_diffSpp[i]<- Beta.div_adapt(Y = WP3_diffSpp, 
                                   dist_spp = organize.syncsa(comm = WP3_diffSpp, phylodist = cophenetic(tree_toy2))$phylodist, 
                                   nperm = 1, method = "dist")$BDextend.obs
  resP3_diffSpp_raw[i]<- Beta.div_adapt(Y = WP3_diffSpp, 
                                       dist_spp = organize.syncsa(comm = WP3_diffSpp, phylodist = cophenetic(tree_toy2))$phylodist, 
                                       nperm = 1, method = "raw")$BDextend.obs
  resP3_noPhy[i]<- beta.div(Y = WP3_diffSpp, method = "percentdiff", sqrt.D = T, nperm = 1)$beta[2]
}
quartz()
windows()
plot(increase_abund, resP3_diffSpp, type = "b", pch= 19, ylim= c(0.0, 0.65))
points(increase_abund, resP3_diffSpp_raw, type = "b", pch= 2)
points(increase_abund, resP3_noPhy, type = "b", pch= 3)
legend("topleft", legend = c("sqrt Bray", "raw", "original"), pch = c(19,2, 3), lty = 1)

####testing P4 NS p5- doble zero assimetry and maximun value#####
#adding species
WP4_eqqSpp<- W[3:4,]
resP4_eqqSpp<- numeric(length= 4)
resP4_eqqSpp_raw<- numeric(length= 4)
resP4_noPhy<- numeric(length= 4)
for(i in 1:4){
  WP4_eqqSpp[2,i]<- 1
  resP4_eqqSpp[i]<- Beta.div_adapt(Y = WP4_eqqSpp, 
                 dist_spp = organize.syncsa(comm = WP4_eqqSpp, phylodist = cophenetic(tree_toy2))$phylodist, 
                 nperm = 1, method = "dist")$BDextend.obs
  resP4_eqqSpp_raw[i]<- Beta.div_adapt(Y = WP4_eqqSpp, 
                 dist_spp = organize.syncsa(comm = WP4_eqqSpp, phylodist = cophenetic(tree_toy2))$phylodist, 
                 nperm = 1, method = "raw")$BDextend.obs
  resP4_noPhy[i]<- beta.div(Y = WP4_eqqSpp, method = "percentdiff", sqrt.D = T, nperm = 1)$beta[2]
}

quartz()
windows()
plot(1:4, resP4_eqqSpp, type = "b", pch= 19, ylim= c(0.0, 0.8), xlab= "addition of 1 species", ylab= "Beta-div")
points(1:4, resP4_eqqSpp_raw, type = "b", pch= 2)
points(1:4, resP4_noPhy, type = "b", pch= 5)
legend("topright", legend = c("sqrt Bray", "raw", "original"), pch = c(19,2, 5), lty = 1)

#removing species - adding doble-zeros
WP4_eqqSpp<- W[3:4,]
resP4.2_eqqSpp<- numeric(length= 3)
resP4.2_eqqSpp_raw<- numeric(length= 3)
resP4.2_noPhy<- numeric(length= 3)
subst1<- 1:3
subst2<- 5:7
for(i in 1:length(subst1)){
  WP4_eqqSpp[1 ,subst1[i]]<- 0
  WP4_eqqSpp[2 ,subst2[i]]<- 0
  resP4.2_eqqSpp[i]<- Beta.div_adapt(Y = WP4_eqqSpp, 
                                   dist_spp = cophenetic(tree_toy2), nperm = 1, method = "dist")$BDextend.obs
  resP4.2_eqqSpp_raw[i]<- Beta.div_adapt(Y = WP4_eqqSpp, 
                                       dist_spp = cophenetic(tree_toy2), nperm = 1, method = "raw")$BDextend.obs
  resP4.2_noPhy[i]<- beta.div(Y = WP4_eqqSpp, method = "percentdiff", sqrt.D = T, nperm = 1)$beta[2]
}

quartz()
windows()
plot(subst1, resP4.2_eqqSpp, type = "b", pch= 19, ylim= c(0.0, 1), xlab= "removal of 1 species", ylab= "Beta-div", cex= 1.2)
points(subst1, resP4.2_eqqSpp_raw, type = "b", pch= 2)
points(subst1, resP4.2_noPhy, type = "b", pch= 5)
legend("bottomright", legend = c("sqrt Bray", "raw", "original"), pch = c(1,2, 5), lty = 1)

####P6 -nestedness and monotonicity in nestedness properties
resP5<- numeric(length= 4)
resP5_raw<- numeric(length= 4)
resP5_noPhy<- numeric(length= 4)

WP5<- W[3,]
S4<- rep(0, ncol(W))
WP5<- rbind(WP5, S4)
WP5[2,1]<- 1
for(i in 1:4){
  WP5[1,4+i]<- 1
  resP5[i]<- Beta.div_adapt(Y = WP5, dist_spp = cophenetic(tree_toy2), nperm = 1, method = "dist")$BDextend.obs
  resP5_raw[i]<- Beta.div_adapt(Y = WP5, dist_spp = cophenetic(tree_toy2), nperm = 1, method = "raw")$BDextend.obs
  resP5_noPhy[i]<- beta.div(Y = WP5, method = "percentdiff", sqrt.D = T, nperm = 1)$beta[2]
}

quartz()
windows()
plot(1:4, resP5, type = "b", pch= 19, ylim= c(0.0, 0.5), xlab= "degree of nestedness", ylab= "Beta-div", cex= 1.2)
points(1:4, resP5_raw, type = "b", pch= 2)
points(1:4, resP5_noPhy, type = "b", pch= 5)
legend("bottomright", legend = c("sqrt Bray", "raw", "original"), pch = c(1,2, 5), lty = 1)

#evaluating monotonicity
tree_toy<- "(((((A:3,B:3):2,(C:3,D:3):2):5,((E:3,F:3):2,(G:3,H:3):2):5):5):5,((((I:3,J:3):2,(K:3,L:3):2):5,((M:3,N:3):2,(O:3,P:3):2):5):5):5):10;"
tree_toy<- read.tree(text= tree_toy)
W<- matrix(c(3,7,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             7,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
             1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,
             0,0,0,0,1,1,1,1,0,0,0,0,0,0,0,0),
           nrow= 4, ncol= length(tree_toy$tip.label), byrow = T, 
           dimnames= list(c("S1", "S2", "S3", "S4"), 
                          tree_toy$tip.label
           )
)

WP5.2<- W[c(3, 4),]
comm2<- c(1, rep(0, ncol(W)-1))
WP5.2<- rbind(comm1= WP5.2[1,], comm2)
size_samp<- c(2, 4, 6, 8, 10)
resP5.2<- numeric(length= 4)
resP5.2_raw<- numeric(length= 4)
resP5.2_noPhy<- numeric(length= 4)
list_WP5.2<- vector(mode = "list", length = length(size_samp))

#sampling species at random - any clade
for(i in 1:length(size_samp)){
  samp_W<- matrix(sample(5:16, size = size_samp[i], replace = FALSE), nrow= 2, ncol= size_samp[i]/2, byrow= T)
  WP5.2[1,samp_W[1,]]<- 1
  WP5.2[2,samp_W[2,]]<- 1
  list_WP5.2[[i]]<- WP5.2
  org<- organize.syncsa(comm = WP5.2, phylodist = cophenetic(tree_toy))
  coph_dist<- org$phylodist
  WP5.2<- org$community
  #beta div calculation
  resP5.2[i]<- Beta.div_adapt(Y = WP5.2, dist_spp = coph_dist, nperm = 1, method = "dist")$BDextend.obs
  resP5.2_raw[i]<- Beta.div_adapt(Y = WP5.2, dist_spp = coph_dist, nperm = 1, method = "raw")$BDextend.obs
  resP5.2_noPhy[i]<- beta.div(Y = WP5.2, method = "percentdiff", sqrt.D = T, nperm = 1)$beta[2]
  
  #restart the process
  
  WP5.2<- W[c(3, 4),]
  comm2<- c(1, rep(0, ncol(W)-1))
  WP5.2<- rbind(comm1= WP5.2[1,], comm2)
  
}
quartz()
windows()
plot(size_samp, resP5.2, type = "b", pch= 19, ylim= c(0.0, 0.7), xlab= "degree of nestedness", ylab= "Beta-div", cex= 1.2)
points(size_samp, resP5.2_raw, type = "b", pch= 2)
points(size_samp, resP5.2_noPhy, type = "b", pch= 5)
legend("bottomright", legend = c("sqrt Bray", "raw", "original"), pch = c(19,2, 5), lty = 1)

#sampling species from different clades
size_samp2<- c(1, 2, 3, 4)
resP5.2_nestClade<- numeric(length= length(size_samp2))
resP5.2_raw_nestClade<- numeric(length= length(size_samp2))
resP5.2_noPhy<- numeric(length= 4)
list_WP5.2_nestClade<- vector(mode = "list", length = length(size_samp2))


for(i in 1:length(size_samp2)){
  samp_clade1<- sample(LETTERS[seq(5, 8)], size_samp2[i], replace= FALSE)
  samp_clade2<- sample(LETTERS[seq(9, 16)], size_samp2[i], replace= FALSE)
  WP5.2[1,samp_clade1]<- 1
  WP5.2[2,samp_clade2]<- 1
  list_WP5.2_nestClade[[i]]<- WP5.2
  org<- organize.syncsa(comm = WP5.2, phylodist = cophenetic(tree_toy))
  coph_dist<- org$phylodist
  WP5.2<- org$community
  #beta div calculation
  resP5.2_nestClade[i]<- Beta.div_adapt(Y = WP5.2, dist_spp = coph_dist, nperm = 1, method = "dist")$BDextend.obs
  resP5.2_raw_nestClade[i]<- Beta.div_adapt(Y = WP5.2, dist_spp = coph_dist, nperm = 1, method = "raw")$BDextend.obs
  resP5.2_noPhy[i]<- beta.div(Y = WP5.2, method = "percentdiff", sqrt.D = T, nperm = 1)$beta[2]
  
  #restart the process
  
  WP5.2<- W[c(3, 4),]
  comm2<- c(1, rep(0, ncol(W)-1))
  WP5.2<- rbind(comm1= WP5.2[1,], comm2)
  
}

windows()
plot(size_samp2, resP5.2_nestClade, type = "b", pch= 19, ylim= c(0.0, 1), xlab= "degree of nestedness", ylab= "Beta-div", cex= 1.2)
points(size_samp2, resP5.2_raw_nestClade, type = "b", pch= 2)
points(size_samp2, resP5.2_noPhy, type = "b", pch= 5)
legend("bottomright", legend = c("sqrt Bray", "raw", "original"), pch = c(1,2, 5), lty = 1)

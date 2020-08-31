LCBD.sigPar<-function(data,nperm,runs,family, nspp, ncomm){
  comm<-as.matrix(data[1][[1]]) #extracting community data
  dist.spp<-as.matrix(data[2][[1]]) #extracting phylogenetic distances
  dist.trait<-as.matrix(data[3][[1]]) #extracting trait distances
  envir<-as.numeric(data[4][[1]]) #extracting environmental vector 
  matrixP.obs<- SYNCSA::matrix.p(ifelse(comm>=1,1,0),dist.spp)$matrix.P #calculating matrix P for incidence data
  matrixX.obs<-SYNCSA::matrix.p(ifelse(comm>=1,1,0),dist.trait)$matrix.P #calculating matrix X for incidence data
  #LCBDP<-beta.div(matrixP.obs,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm)$LCBD 
  #LCBDX<-beta.div(matrixX.obs,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm)$LCBD
  #LCBDW<-beta.div(ifelse(comm>=1,1,0),method ="euclidean",nperm=nperm)$LCBD
  BDp_all_dist<-beta.div(matrixP.obs,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm) #calculating distance based BDp
  BDx_all_dist<-beta.div(matrixX.obs,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm) #calculating distance based BDx
  BDw_all_dist<-beta.div(ifelse(comm>=1,1,0),method ="euclidean",nperm=nperm) #calculating distance based BD
  BDp_all_raw<- beta.div(decostand(matrixP.obs, method= "normalize"), method = "none", sqrt.D = FALSE, 2) #calculating raw data BDp
  BDx_all_raw<- beta.div(decostand(matrixX.obs, method= "normalize"), method = "none", sqrt.D = FALSE, 2) #calculating raw data BDx
  BDw_all_raw<-beta.div(decostand(comm>=1,1,0, method = "normalize"), method ="none", nperm=nperm) #calculating raw data BD
  LCBDP_dist<- BDp_all_dist$LCBD #extracting distance based PLCBD component 
  LCBDX_dist<- BDx_all_dist$LCBD #extracting distance based XLCBD component
  LCBDW_dist<- BDw_all_dist$LCBD #extracting distance based LCBD component
  LCBD_all_dist<-as.matrix(cbind(LCBDP_dist,LCBDX_dist,LCBDW_dist)) #binding distance based PLCBD, XLCBD and LCBD components
  LCBDP_raw<- BDp_all_raw$LCBD #extracting raw data PLCBD component
  LCBDX_raw<- BDx_all_raw$LCBD #extracting raw data XLCBD component
  LCBDW_raw<- BDw_all_raw$LCBD #extracting raw data LCBD component
  LCBD_all_raw<-as.matrix(cbind(LCBDP_raw,LCBDX_raw,LCBDW_raw)) #binding all raw data LCBD components
  BDp_dist<- BDp_all_dist$SStotal_BDtotal[2] #extracting distance based BDp 
  BDx_dist<- BDx_all_dist$SStotal_BDtotal[2] #extracting distance based BDx
  BDw_dist<- BDw_all_dist$SStotal_BDtotal[2] #extracting distance based BD
  BD_all_dist<- as.matrix(cbind(BDp_dist,BDx_dist,BDw_dist)) #binding all distance based BDs 
  BDp_raw<- BDp_all_raw$SStotal_BDtotal[2] #extracting raw data BDp
  BDx_raw<- BDx_all_raw$SStotal_BDtotal[2] #extracting raw data BDx
  BDw_raw<- BDw_all_raw$SStotal_BDtotal[2] #extracting raw data BD
  BD_all_raw<- as.matrix(cbind(BDp_raw,BDx_raw,BDw_raw)) #binding all raw data BDs
  data_obs_dist <- as.data.frame(cbind(LCBDP_dist,LCBDX_dist,LCBDW_dist,envir)) #binding in a data frame all distance based LCBDs
  data_obs_raw <- as.data.frame(cbind(LCBDP_raw,LCBDX_raw,LCBDW_raw,envir)) #bindign in a data frame all raw data LCBDs
  mod_obsP_dist <- glm(LCBDP_dist~poly(envir,2), data = data_obs_dist, family = family) #runing a glm for distance based PLCBD
  f_obsP_dist <- summary.lm(mod_obsP_dist)$fstatistic[1] #extracting the F ratio for PLCBDxEnvir relation
  r2_obsP_dist<-summary.lm(mod_obsP_dist)$adj.r.squared #extracting the Rsquared for PLCBDxEnvir relation 
  mod_obsX_dist <- glm(LCBDX_dist~poly(envir,2), data = data_obs_dist, family = family) #runing a glm for distance based XLCBD
  f_obsX_dist <- summary.lm(mod_obsX_dist)$fstatistic[1] #extracting the Fratio for XLCBDxEnvir relation
  r2_obsX_dist<-summary.lm(mod_obsX_dist)$adj.r.squared #extracting the Rsquared for XLCBDxEnvir relation
  mod_obsW_dist <- glm(LCBDW_dist~poly(envir,2), data = data_obs_dist, family = family) #running a glm for distance based LCBD 
  f_obsW_dist <- summary.lm(mod_obsW_dist)$fstatistic[1] #extracting the Fratio for LCBDxEnvir relation
  r2_obsW_dist <- summary.lm(mod_obsW_dist)$adj.r.squared #extracting the Fratio for LCBDxEnvir relation
  mod_obsP_raw <- glm(LCBDP_raw~poly(envir,2), data = data_obs_raw, family = family) #running glm for raw data PLCBD
  f_obsP_raw <- summary.lm(mod_obsP_raw)$fstatistic[1] #extracting Fratio for raw data PLCBD 
  r2_obsP_raw<-summary.lm(mod_obsP_raw)$adj.r.squared #extracting Rsquared for raw data LCBD
  mod_obsX_raw <- glm(LCBDX_raw~poly(envir,2), data = data_obs_raw, family = family) #running glm for raw data XLCBD
  f_obsX_raw <- summary.lm(mod_obsX_raw)$fstatistic[1] #extracting Fratio for raw data XLCBD 
  r2_obsX_raw<-summary.lm(mod_obsX_raw)$adj.r.squared #extracting rsquared for raw data XLCBD
  mod_obsW_raw <- glm(LCBDW_raw~poly(envir,2), data = data_obs_raw, family = family) #fitting glm for raw data LCBD
  f_obsW_raw <- summary.lm(mod_obsW_raw)$fstatistic[1] #extracting Fratio for raw data LCBD
  r2_obsW_raw <- summary.lm(mod_obsW_raw)$adj.r.squared #extracting Rsquared for raw data LCBD
  F_null_siteP_dist <- matrix(NA, runs, 1) #matrix to receive site null model results for distance based PLCBD
  F_null_taxaP_dist <- matrix(NA, runs, 1) #matrix to receive taxa null model results for distance based PLCBD
  F_null_siteX_dist <- matrix(NA, runs, 1) #matrix to receive site null model results for distance based XLCBD
  F_null_taxaX_dist <- matrix(NA, runs, 1) #matrix to receive site null model results for distance based XLCBD
  F_null_siteW_dist <- matrix(NA, runs, 1) #matrix to receive site null model results for distance based LCBD
  F_null_taxaW_dist <- matrix(NA, runs, 1) #matrix to receive site null model results for distance based LCBD
  F_null_siteP_raw <- matrix(NA, runs, 1) #matrix to receive site null model results for raw data PLCBD
  F_null_taxaP_raw <- matrix(NA, runs, 1) #matrix to receive taxa null model results for raw based PLCBD
  F_null_siteX_raw <- matrix(NA, runs, 1) #matrix to receive site null model results for raw data XLCBD
  F_null_taxaX_raw <- matrix(NA, runs, 1) #matrix to receive taxa null model results for raw data XLCBD
  F_null_siteW_raw <- matrix(NA, runs, 1) #matrix to receive site null model results for raw based LCBD
  F_null_taxaW_raw <- matrix(NA, runs, 1) #matrix to receive taxa null model results for raw data LCBD - this doesn't make any sense
  BD_null_all_dist<- matrix(NA, nrow= runs, ncol= 3, dimnames= list(paste("run", 1:runs), 
                                                                    c("BDp_null_dist", "BDx_null_dist", "BDw_null_dist"))) #receive the null BDs with distance approach
  BD_null_all_raw<- matrix(NA, nrow= runs, ncol= 3, dimnames= list(paste("run", 1:runs), 
                                                                   c("BDp_null_raw", "BDx_null_raw", "BDw_null_raw"))) #receive the null BDs with raw approach
  #starting the null model for all metrics                                                              
  for (k in 1:runs) {
    distspp_null <- picante::taxaShuffle(dist.spp) #taxa shuffle matrix P
    distrait_null<- picante::taxaShuffle(dist.trait) #taxa shuffle matrix X
    match.namesP <- match(colnames(comm), colnames(distspp_null)) 
    match.namesX<- match(colnames(comm), colnames(distrait_null))
    matrixP_null<- SYNCSA::matrix.p(ifelse(comm>=1,1,0),as.matrix(distspp_null[match.namesP,match.namesP]))$matrix.P #matrix P taxa shuffle
    matrixX_null<- SYNCSA::matrix.p(ifelse(comm>=1,1,0),as.matrix(distrait_null[match.namesX,match.namesX]))$matrix.P #matrix X taxa shuffle
    matrixW_null<- EcoSimR::sim1(t(ifelse(comm>=1,1,0)))
    matrixW_null<- t(matrixW_null)
    rownames(matrixW_null)<- rownames(comm)
    colnames(matrixW_null)<- colnames(comm)
    #distspp_null <- picante::taxaShuffle(dist.spp)
    #distrait_null<- picante::taxaShuffle(dist.trait)
    #match.namesP <- match(colnames(comm), colnames(distspp_null))
    #match.namesX<- match(colnames(comm), colnames(distrait_null))
    #matrixP_null<- SYNCSA::matrix.p(ifelse(comm>=1,1,0),as.matrix(distspp_null[match.namesP,match.namesP]))$matrix.P
    #matrixP_null<- SYNCSA::matrix.p(comm = ifelse(comm>=1,1,0), phylodist = null_distP)$matrix.P
    #matrixX_null<- SYNCSA::matrix.p(ifelse(comm>=1,1,0),as.matrix(distrait_null[match.namesX,match.namesX]))$matrix.P
    #matrixX_null<- SYNCSA::matrix.x(comm = ifelse(comm>=1,1,0), traits = as.matrix(null_distX), 
    #                                scale = FALSE)$matrix.X
    #LCBDP_null<- beta.div(matrixP_null,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm)$LCBD
    #LCBDX_null<- beta.div(matrixX_null,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm)$LCBD
    #LCBDW_null<-LCBDW[sample(1:length(LCBDW))]
    BDp_all_nullDist<- beta.div(matrixP_null,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm) #null distance based BDp
    BDx_all_nullDist<- beta.div(matrixX_null,method ="percentagedifference",sqrt.D = TRUE,nperm= nperm) #null distance based BDx
    BDw_all_nullDist<- beta.div(matrixW_null, method ="euclidean", nperm= nperm) #null distance based BD
    BDp_all_nullRaw<- beta.div(decostand(matrixP_null, method= "normalize"),method ="none",sqrt.D = TRUE,nperm= nperm) #null raw data approach for BDp
    BDx_all_nullRaw<- beta.div(decostand(matrixX_null, method= "normalize"),method ="none",sqrt.D = TRUE,nperm= nperm) #null raw data approach for BDx
    BDw_all_nullRaw<- beta.div(decostand(matrixW_null, method= "normalize"), method ="none", nperm= nperm) #null raw data approach for BD
    LCBDW_nullDist<- BDw_all_nullDist$LCBD #WLCBD null
    LCBDP_nullDist<- BDp_all_nullDist$LCBD #PLCBD null
    LCBDX_nullDist<- BDx_all_nullDist$LCBD #XLCBD null
    LCBDW_nullRaw<- BDw_all_nullRaw$LCBD #WLCBD null
    LCBDP_nullRaw<- BDp_all_nullRaw$LCBD #PLCBD null
    LCBDX_nullRaw<- BDx_all_nullRaw$LCBD #XLCBD null
    BDp_nullDist<- BDp_all_nullDist$SStotal_BDtotal[2] #BDphy null
    BDx_nullDist<- BDx_all_nullDist$SStotal_BDtotal[2] #BDfun null
    BDw_nullDist<- BDw_all_nullDist$SStotal_BDtotal[2] #BDw null 
    BDp_nullRaw<- BDp_all_nullRaw$SStotal_BDtotal[2] #BDphy null
    BDx_nullRaw<- BDx_all_nullRaw$SStotal_BDtotal[2] #BDfun null
    BDw_nullRaw<- BDw_all_nullRaw$SStotal_BDtotal[2] #BDw null 
    LCBD_all_nullDist<-cbind(LCBDP_nullDist,LCBDX_nullDist,LCBDW_nullDist) #binding all null distance based LCBD
    LCBD_all_nullRaw<-cbind(LCBDP_nullRaw,LCBDX_nullRaw,LCBDW_nullRaw) #binding all null raw data LCBD
    data_null_taxaDist<- as.data.frame(cbind(LCBD_all_nullDist, 
                                             envir)) #binding null dist LCBD with envir 
    data_null_taxaRaw<- as.data.frame(cbind(LCBD_all_nullRaw, 
                                            envir)) #binding null raw data LCBD with envir
    BD_null_all_dist[k,]<- c(BDp_nullDist,BDx_nullDist,BDw_nullDist) #binding all distance based null BDs
    BD_null_all_raw[k,]<- c(BDp_nullRaw, BDx_nullRaw, BDw_nullRaw) #binding all raw data null BDs
    modP_null_taxaDist <- glm(LCBDP_nullDist~poly(envir,2), data = data_null_taxaDist, 
                              family = family)
    modX_null_taxaDist<- glm(LCBDX_nullDist~poly(envir,2),data=data_null_taxaDist,family=family)
    modW_null_taxaDist<- glm(LCBDW_nullDist~envir,data=data_null_taxaDist,family=family)
    modP_null_taxaRaw <- glm(LCBDP_nullRaw~poly(envir,2), data = data_null_taxaRaw, 
                             family = family)
    modX_null_taxaRaw<- glm(LCBDX_nullRaw~poly(envir,2),data=data_null_taxaRaw,family=family)
    modW_null_taxaRaw<- glm(LCBDW_nullRaw~envir,data=data_null_taxaRaw,family=family)
    F_null_taxaP_dist[k, ] <- summary.lm(modP_null_taxaDist)$fstatistic[1]
    F_null_taxaX_dist[k, ] <- summary.lm(modX_null_taxaDist)$fstatistic[1]
    F_null_taxaW_dist[k, ] <- summary.lm(modW_null_taxaDist)$fstatistic[1]
    F_null_taxaP_raw[k, ] <- summary.lm(modP_null_taxaRaw)$fstatistic[1]
    F_null_taxaX_raw[k, ] <- summary.lm(modX_null_taxaRaw)$fstatistic[1]
    F_null_taxaW_raw[k, ] <- summary.lm(modW_null_taxaRaw)$fstatistic[1]
    LCBDP_null_siteDist <- LCBD_all_dist[sample(1:dim(LCBD_all_dist)[1]),"LCBDP_dist"]
    names(LCBDP_null_siteDist)<-rownames(LCBD_all_dist)
    LCBDP_null_siteRaw <- LCBD_all_raw[sample(1:dim(LCBD_all_raw)[1]),"LCBDP_raw"]
    names(LCBDP_null_siteRaw)<-rownames(LCBD_all_raw)
    LCBDX_null_siteDist <- LCBD_all_dist[sample(1:dim(LCBD_all_dist)[1]),"LCBDX_dist"]
    names(LCBDX_null_siteDist)<-rownames(LCBD_all_dist)
    LCBDX_null_siteRaw <- LCBD_all_raw[sample(1:dim(LCBD_all_raw)[1]),"LCBDX_raw"]
    names(LCBDX_null_siteRaw)<-rownames(LCBD_all_raw)
    LCBDW_null_siteDist <- LCBD_all_dist[sample(1:dim(LCBD_all_dist)[1]),"LCBDW_dist"]
    names(LCBDW_null_siteDist)<-rownames(LCBD_all_dist)
    LCBDW_null_siteRaw <- LCBD_all_raw[sample(1:dim(LCBD_all_raw)[1]),"LCBDW_raw"]
    names(LCBDW_null_siteRaw)<-rownames(LCBD_all_raw)
    data_null_siteDist <- as.data.frame(cbind(LCBDP_null_siteDist,LCBDX_null_siteDist,LCBDW_null_siteDist,envir))
    data_null_siteRaw <- as.data.frame(cbind(LCBDP_null_siteRaw,LCBDX_null_siteRaw,LCBDW_null_siteRaw,envir))
    modP_null_siteDist <- glm(LCBDP_null_siteDist~envir, data = data_null_siteDist, 
                              family = family)
    modP_null_siteRaw <- glm(LCBDP_null_siteRaw~envir, data = data_null_siteRaw, 
                             family = family)
    modX_null_siteDist <- glm(LCBDX_null_siteDist~envir, data = data_null_siteDist, 
                              family = family)
    modX_null_siteRaw <- glm(LCBDX_null_siteRaw~envir, data = data_null_siteRaw, 
                             family = family)
    modW_null_siteDist <- glm(LCBDW_null_siteDist~envir, data = data_null_siteDist, 
                              family = family)
    modW_null_siteRaw <- glm(LCBDW_null_siteRaw~envir, data = data_null_siteRaw, 
                             family = family)
    F_null_siteP_dist[k, ] <- summary.lm(modP_null_siteDist)$fstatistic[1]
    F_null_siteX_dist[k, ] <- summary.lm(modX_null_siteDist)$fstatistic[1]
    F_null_siteW_dist[k, ] <- summary.lm(modW_null_siteDist)$fstatistic[1]
    F_null_siteP_raw[k, ] <- summary.lm(modP_null_siteRaw)$fstatistic[1]
    F_null_siteX_raw[k, ] <- summary.lm(modX_null_siteRaw)$fstatistic[1]
    F_null_siteW_raw[k, ] <- summary.lm(modW_null_siteRaw)$fstatistic[1]
  }
  #standardized effect size for BDs - this does not work for BD, makes sense only for BDp and BDx
  mean_nullDist<- apply(BD_null_all_dist,2, mean)
  mean_nullRaw<- apply(BD_null_all_raw,2, mean)
  sd_nullDist<- apply(BD_null_all_dist, 2, sd)
  sd_nullRaw<- apply(BD_null_all_raw, 2, sd)
  ses_BDpDist<- (BD_all_dist[,"BDp_dist"] - mean_nullDist[1])/(sd_nullDist[1])
  ses_BDxDist<- (BD_all_dist[,"BDx_dist"] - mean_nullDist[2])/(sd_nullDist[2])
  ses_BDwDist<- (BD_all_dist[,"BDw_dist"] - mean_nullDist[3])/(sd_nullDist[3])
  ses_BDpRaw<- (BD_all_raw[,"BDp_raw"] - mean_nullRaw[1])/(sd_nullRaw[1])
  ses_BDxRaw<- (BD_all_raw[,"BDx_raw"] - mean_nullRaw[2])/(sd_nullRaw[2])
  ses_BDwRaw<- (BD_all_raw[,"BDw_raw"] - mean_nullRaw[3])/(sd_nullRaw[3])
  p_taxaPDist <- (sum(ifelse(F_null_taxaP_dist >= f_obsP_dist, 1, 0)) + 1)/(runs + 1)
  p_sitePDist <- (sum(ifelse(F_null_siteP_dist >= f_obsP_dist, 1, 0)) + 1)/(runs + 1)
  p_taxaXDist <- (sum(ifelse(F_null_taxaX_dist >= f_obsX_dist, 1, 0)) + 1)/(runs + 1)
  p_siteXDist <- (sum(ifelse(F_null_siteX_dist >= f_obsX_dist, 1, 0)) + 1)/(runs + 1)
  p_taxaWDist <- (sum(ifelse(F_null_taxaW_dist >= f_obsW_dist, 1, 0)) + 1)/(runs + 1)
  p_siteWDist <- (sum(ifelse(F_null_siteW_dist >= f_obsW_dist, 1, 0)) + 1)/(runs + 1)
  p_taxaPRaw <- (sum(ifelse(F_null_taxaP_raw >= f_obsP_raw, 1, 0)) + 1)/(runs + 1)
  p_sitePRaw <- (sum(ifelse(F_null_siteP_raw >= f_obsP_raw, 1, 0)) + 1)/(runs + 1)
  p_taxaXRaw <- (sum(ifelse(F_null_taxaX_raw >= f_obsX_raw, 1, 0)) + 1)/(runs + 1)
  p_siteXRaw <- (sum(ifelse(F_null_siteX_raw >= f_obsX_raw, 1, 0)) + 1)/(runs + 1)
  p_taxaWRaw <- (sum(ifelse(F_null_taxaW_raw >= f_obsW_raw, 1, 0)) + 1)/(runs + 1)
  p_siteWRaw <- (sum(ifelse(F_null_siteW_raw >= f_obsW_raw, 1, 0)) + 1)/(runs + 1)
  p_BDpDist<- (sum(ifelse(BD_null_all_dist[,"BDp_null_dist"] >= BDp_dist, 1, 0)) + 1)/(runs + 1)
  p_BDxDist<- (sum(ifelse(BD_null_all_dist[,"BDx_null_dist"] >= BDx_dist, 1, 0)) + 1)/(runs + 1)
  p_BDwDist<- (sum(ifelse(BD_null_all_dist[,"BDw_null_dist"] >= BDw_dist, 1, 0)) + 1)/(runs + 1)
  p_BDpRaw<- (sum(ifelse(BD_null_all_raw[,"BDp_null_raw"] >= BDp_raw, 1, 0)) + 1)/(runs + 1)
  p_BDxRaw<- (sum(ifelse(BD_null_all_raw[,"BDx_null_raw"] >= BDx_raw, 1, 0)) + 1)/(runs + 1)
  p_BDwRaw<- (sum(ifelse(BD_null_all_raw[,"BDw_null_raw"] >= BDw_raw, 1, 0)) + 1)/(runs + 1)
  matrix_resp_Raw<- matrix(c(BD_all_raw,ses_BDpRaw,ses_BDxRaw,ses_BDwRaw,p_taxaPRaw,p_sitePRaw,p_taxaXRaw,p_siteXRaw,p_taxaWRaw,p_siteWRaw,p_BDpRaw, p_BDxRaw, p_BDwRaw,
                             f_obsP_raw,f_obsX_raw,f_obsW_raw,r2_obsP_raw,r2_obsX_raw,r2_obsW_raw),
                           nrow=1,ncol=21,dimnames=list(1,c("BDpRaw","BDxRaw","BDwRaw","sesBDpRaw","sesBDxRaw","sesBDwRaw","p_taxaPRaw","p_sitePRaw","p_taxaXRaw","p_siteXRaw",
                                                            "p_taxaWRaw","p_siteWRaw","p_BDpRaw","p_BDxRaw","p_BDwRaw","F_obsPRaw","F_obsXRaw","F_obsWRaw","R2_PRaw","R2_XRaw","R2_WRaw")))
  matrix_resp_Dist<-matrix(c(BD_all_dist,ses_BDpDist,ses_BDxDist,ses_BDwDist,p_taxaPDist,p_sitePDist,p_taxaXDist,p_siteXDist,p_taxaWDist,p_siteWDist,p_BDpDist, p_BDxDist, p_BDwDist,
                             f_obsP_dist,f_obsX_dist,f_obsW_dist,r2_obsP_dist,r2_obsX_dist,r2_obsW_dist),
                           nrow=1,ncol=21,dimnames=list(1,c("BDp","BDx","BDw","sesBDp","sesBDx","sesBDw","p_taxaP","p_siteP","p_taxaX","p_siteX","p_taxaW","p_siteW","p_BDp","p_BDx","p_BDw","F_obsP","F_obsX","F_obsW","R2_P","R2_X","R2_W")))
  return(list(Dist_results= matrix_resp_Dist, Raw_results= matrix_resp_Raw))
}

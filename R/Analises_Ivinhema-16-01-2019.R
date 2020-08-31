##################################################################################################
#An?lises de diversidade Beta - Legendre e Fuzzy Groups com dados de Riachos da Bacia do Ivinhema#
##################################################################################################
save.image("Analises_Ivinhema-16-01-2019.RData")
#Pacotes necess?rios
library(ape)
library(picante)
library(PCPS)
library(phytools)
library(vegan)
library(geiger)


#apenas verificar sinal filogen?tico nos atributos atrav?s da estat?stica K
#atributos funcionais
atributos<- read.table("clipboard",header=TRUE)
atributos

resul.K<-list()
for(i in 1:ncol(traits.Ivinhema)){
    resul.K[[i]]<-phylosignal(traits.Ivinhema[phylo.Ivinhema$tip.label,i],phy = multi2di(phylo.Ivinhema),
                checkdata= FALSE)
    resul.K
}

names(resul.K)<-c("AO","PVO","Alr","RCau","CRPED","LRBo")
resul.K


#carregando dados necess?rios para an?lises
comu.Ivinhema<-comu.Ivinhema.completa[,3:ncol(comu.Ivinhema.completa)] #comunidades de riachos do Ivinhema
comu.Ivinhema.completa<-read.table("comu_Ivinhema.txt")
head(comu.Ivinhema.completa)
names(comu.Ivinhema.completa)
head(comu.Ivinhema.completa)
write.csv(comu.Ivinhema.completa, "comu.Ivinhema.completa.csv")
comu.Ivinhema.completa$porcao
rowSums(comu.Ivinhema)
#filogenia
phylo.Ivinhema<-read.tree("filogenia.ivinhema.riachos.txt") #Filogenia para comunidades do Ivinhema
windows()
plot(phylo.Ivinhema) #plotando filogenias do Ivinhema
dist.phylo.Ivinhema<-cophenetic(phylo.Ivinhema)
comu.Ivinhema.completa$ponto
#atributos funcionais
atributos_completo<-read.table("atributos_Ivinhema.txt",header=TRUE)
traits.Ivinhema<- atributos_completo[,c(-1,-8)]
head(traits.Ivinhema)

#dist?ncia funcional
dist.functional<-dist(traits.Ivinhema, method = "euclidean")

#matriz ambiental
amb<-read.table("amb.txt",header=TRUE)
str(amb)
amb<-decostand(amb[,2:ncol(amb)],"standardize")
dist.amb<-dist(amb,"euclidian")



#calculando os pcps filogen?ticos para as comunidades
PCPS.phylo.Ivinhema<-pcps(comm = comu.Ivinhema[,colnames(dist.phylo.Ivinhema)],dist.phylo.Ivinhema)
scores(PCPS.phylo.Ivinhema)
windows()
plot(PCPS.phylo.Ivinhema$vectors[,1],PCPS.phylo.Ivinhema$vectors[,2],type="n")
text(PCPS.phylo.Ivinhema$vectors[,1],PCPS.phylo.Ivinhema$vectors[,2],
     labels=comu.Ivinhema.completa$porcao)
plot.pcps(PCPS.phylo.Ivinhema)
abline(h=0,v=0,lty=2)
arrows(x0=0,y0=0,x1=PCPS.phylo.Ivinhema$correlations[,1],y1=PCPS.phylo.Ivinhema$correlations[,2])

windows()
plot(PCPS.phylo.Ivinhema$correlations[,1],PCPS.phylo.Ivinhema$correlations[,2],type="n")
text(PCPS.phylo.Ivinhema$correlations[,1],PCPS.phylo.Ivinhema$correlations[,2],
     labels=rownames(PCPS.phylo.Ivinhema$correlations))


#calculando matriz de pondera??o filogen?tica para as comunidades
phylo.p.Ivinhema<-matrix.p(comu.Ivinhema[,rownames(dist.phylo.Ivinhema)],dist.phylo.Ivinhema)$matrix.P
phylo.p.Ivinhema #matrix de pondera??o filogen?tica para o Ivinhema


#calculando os PCPSs funcionais para as comunidades
PCPS.func.Ivinhema<-pcps(comm = comu.Ivinhema[,colnames(as.matrix(dist.functional))],
                         as.matrix(dist.functional))
scores(PCPS.func.Ivinhema)
windows()
plot(PCPS.func.Ivinhema$vectors[,1],PCPS.func.Ivinhema$vectors[,2],type="n")
text(PCPS.func.Ivinhema$vectors[,1],PCPS.func.Ivinhema$vectors[,2],
     labels=comu.Ivinhema.completa$porcao)
abline(h=0,v=0,lty=2)
arrows(x0=0,y0=0,x1=PCPS.func.Ivinhema$correlations[,1],y1=PCPS.func.Ivinhema$correlations[,2])

windows()
plot(PCPS.func.Ivinhema$correlations[,1],PCPS.func.Ivinhema$correlations[,2],type="n")
text(PCPS.func.Ivinhema$correlations[,1],PCPS.func.Ivinhema$correlations[,2],
     labels=rownames(PCPS.func.Ivinhema$correlations))


#calculando a matriz de pondera??o funcional para as comunidades
func.p.Ivinhema<-matrix.p(comu.Ivinhema[,rownames(as.matrix(dist.functional))],
                          as.matrix(dist.functional))$matrix.P

#matriz de pondera??o funcional para o Ivinhema, mesma an?lise que a de cima
func.x.Ivinhema<-matrix.x(comu.Ivinhema[,rownames(as.matrix(traits.Ivinhema))], traits.Ivinhema)

#calculando a rela??o entre estrutura funcional e filogen?tica das comunidades de peixes do rio Ivinhema
mantel(dist(cbind(PCPS.func.Ivinhema$vectors[,1],PCPS.func.Ivinhema$vectors[,2])),
       dist(cbind(PCPS.phylo.Ivinhema$vectors[,1],PCPS.phylo.Ivinhema$vectors[,2])))
windows()
plot(dist(cbind(PCPS.func.Ivinhema$vectors[,1],PCPS.func.Ivinhema$vectors[,2])),
     dist(cbind(PCPS.phylo.Ivinhema$vectors[,1],PCPS.phylo.Ivinhema$vectors[,2])))

#rela??o entre estrutura funcional das comunidades e heterogeineidade ambiental
mantel(dist(cbind(PCPS.func.Ivinhema$vectors[,1],PCPS.func.Ivinhema$vectors[,2])),dist.amb)
windows()
plot(dist.amb,dist(cbind(PCPS.func.Ivinhema$vectors[,1],PCPS.func.Ivinhema$vectors[,2])))

#rela??o entre estrutura filogen?tica das comunidades e heterogeineidade ambiental
mantel(dist(cbind(PCPS.phylo.Ivinhema$vectors[,1],PCPS.phylo.Ivinhema$vectors[,2])),dist.amb)
windows()
plot(dist(cbind(PCPS.phylo.Ivinhema$vectors[,1],PCPS.phylo.Ivinhema$vectors[,2])),dist.amb)

#Calcular a beta diversidade atrav?s do framework proposto por Legendre e C?ceres (2013)
#primeiro ler a fun??o necess?ria para isso

source(file = "beta_div_Leg.txt") #fun??o de Legendre e C?ceres para Beta diversidade
beta.div


#####beta taxonomico######
beta.Tax<-beta.div(Y = comu.Ivinhema, method = "none", save.D = TRUE,clock = TRUE)
#scaling beta Tax
identy_matrix<- diag(length(phylo.Ivinhema$tip.label))
identy_matrix<- ifelse(identy_matrix==1,0,1)
scale_W<- SYNCSA::matrix.p(comm = comu.Ivinhema, phylodist = identy_matrix)$matrix.P
beta.Tax<-beta.div(Y = comu.Ivinhema, method = "none", save.D = TRUE,clock = TRUE)
beta.Tax_dist<- beta.div(Y = comu.Ivinhema, method= "euclidean", sqrt.D = FALSE, nperm = 999)
beta.Tax_dist_scale<- beta.div(Y = scale_W, method = "percentagedifference", sqrt.D = TRUE, nperm= 999)
BDw_dist_scale<- beta.Tax_dist_scale$SStotal_BDtotal[2]
beta.Tax_raw_scale<- beta.div(Y = vegan::decostand(scale_W, method = "normalize"), method = "none", sqrt.D = TRUE, nperm= 999)
BDw_raw_scale<- beta.Tax_raw_scale$SStotal_BDtotal[2]
BDw_dist<- beta.Tax_dist$SStotal_BDtotal[2]
LCBD_dist<- beta.Tax_dist$LCBD
runs<- 999
BDw_nullDist<- matrix(NA, nrow= runs, ncol= 1, dimnames= list(paste("run", 1:runs, sep= ""),
                                                              "BDw"))
for(i in 1:runs){
  comu.Ivinhema_null <- EcoSimR::sim9(speciesData = t(comu.Ivinhema), algo = "sim9", metric = "c_score")
  BDw_all_null<- beta.div(t(comu.Ivinhema_null$Randomized.Data),method ="euclidean",sqrt.D = TRUE,nperm= 999)
  BDw_nullDist[i,]<- BDw_all_null$SStotal_BDtotal[2]
}
BDw_mean_null<- apply(BDw_nullDist,2,mean)
sdW_null<- apply(BDw_nullDist, 2, sd)
ses_BDw<- (BDw_dist - BDw_mean_null)/(sdW_null) #standardized effect size for BDp
BDw_pvalue<- length(which(BDw_nullDist>=BDw_dist))/(runs) #p value for ses_BDp
Tax_SCBD<-beta.Tax$SCBD #contribui??o de esp?cies
Tax_LCBD<-beta.Tax$LCBD #contribui??o dos locais
beta.Tax$p.LCBD #signific?ncia dos locais
which(beta.Tax$p.LCBD<=0.05) #valores menores ou iguais a 0.05
length(which(beta.Tax_dist_scale$p.LCBD<=0.05))
SCBD.tax<-cbind(names(beta.Tax$SCBD),beta.Tax$SCBD)
write.csv(SCBD.tax,"SCBD.tax.csv") #exportando a tabela de contribui??o das esp?cies para div. beta filogen?tica
LCBD.tax<-cbind(names(beta.Tax$LCBD),beta.Tax$LCBD)                               
write.csv(LCBD.tax,"LCBD.tax.csv")


windows()
par(bty="l")
plot(1:length(Tax_SCBD), sort(Tax_SCBD, decreasing= TRUE), type= "n", xlab= "rank", ylab= "Taxonomic SCBD", axes=FALSE,
     xaxs="i", yaxs="i", xlim=c(0,110), ylim=c(0,0.04))
axis(1, at=c(seq(0,110,10)))
axis(2,at=c(seq(0,0.04,0.01)))
points(1:length(Tax_SCBD), sort(Tax_SCBD, decreasing= TRUE), pch= as.numeric(graph_tax[,2]))
text((1:3)+20, c(0.033,0.030, 0.029), labels = paste(names(sort(Tax_SCBD, decreasing= TRUE))[c(1,2,3)]))
lines(1:length(Tax_SCBD), sort(Tax_SCBD, decreasing= TRUE))
dev.copy(png,'taxonomic_SCBD.png')
dev.off()

names_TaxSCBD<-names(sort(Tax_SCBD,TRUE))
dezphylo_in_tax<-match(names(dezSCBD_phylo),names_TaxSCBD)
dezfunc_in_tax<-match(names(dezSCBD_func),names_TaxSCBD)
common_phylofun_in_tax<-match(names(dezSCBD_func)[match(names(dezSCBD_phylo),names(dezSCBD_func))],names_TaxSCBD)[!is.na(match(names(dezSCBD_phylo),names(dezSCBD_func)))]

graph_tax<-matrix(NA, 107, 2, dimnames= list(c(1:length(Tax_SCBD)),c("names","Func/Tax")))
for(i in 1:107){
	graph_tax[i,1]<-names(sort(Tax_SCBD, decreasing= TRUE))[i]
}
graph_tax[dezphylo_in_tax,2]<-21
graph_tax[dezfunc_in_tax,2]<-24
graph_tax[c(-dezphylo_in_tax,-dezfunc_in_tax),2]<-19
graph_tax[common_phylofun_in_tax,2]<-22 #phylo e fun em taxonomic quadrado vazio

common_phylofun_in_tax<-match(names(dezSCBD_func)[match(names(dezSCBD_phylo),names(dezSCBD_func))],names_TaxSCBD)[!is.na(match(names(dezSCBD_phylo),names(dezSCBD_func)))]
names_TaxSCBD[c(107,28)]

#calculando por ordem a contribucao taxonomica
OCBD_tax<- read.table("SCBD_OCBD_tax.txt", header= TRUE)
OCBD_tax
resul_OCBD_tax<-sapply(levels(OCBD_tax$order), function(x){sum(OCBD_tax$SCBD[OCBD_tax$order == x])})
#calculando a contribui??o por sub-bacia
LCBD_porc_Tax<-data.frame(comu.Ivinhema.completa$porcao,Tax_LCBD)
write.csv(LCBD_porc_Tax,"LCBD_tax.csv")
LCBD_Tax_porcao<-sapply(levels(LCBD_porc_Tax$comu.Ivinhema.completa.porcao), function(x){sum(LCBD_porc_Tax$Tax_LCBD[LCBD_porc_Tax$comu.Ivinhema.completa.porcao == x])})
sort(LCBD_Tax_porcao,TRUE)
sum(LCBD_Tax_porcao)



########beta filogenetico#######
beta.Phylo<-beta.div(Y = phylo.p.Ivinhema, method = "none", save.D = TRUE,clock = TRUE)
beta.Phylo_dist<- beta.div(Y = phylo.p.Ivinhema, method = "percentagedifference", sqrt.D = TRUE, nperm = 999)
BDp_dist<- beta.Phylo$SStotal_BDtotal[2] #BDp_dist obs
PLCBD_dist<- beta.Phylo$LCBD
runs<- 999
BDp_nullDist<- matrix(NA, nrow= runs, ncol= 1, dimnames= list(paste("run", 1:runs, sep= ""),
                                                              "BDp"))
for(i in 1:runs){
  distspp_null <- picante::taxaShuffle(dist.phylo.Ivinhema)
  match.namesP <- match(colnames(comu.Ivinhema), colnames(distspp_null))
  matrixP_null<- SYNCSA::matrix.p(comu.Ivinhema, as.matrix(distspp_null[match.namesP,match.namesP]))$matrix.P
  BDp_all_null<- beta.div(matrixP_null,method ="percentagedifference",sqrt.D = TRUE,nperm= 2)
  BDp_nullDist[i,]<- BDp_all_null$SStotal_BDtotal[2]
}
BDp_mean_null<- apply(BDp_nullDist,2,mean)
sdP_null<- apply(BDp_nullDist, 2, sd)
ses_BDp<- (BDp_dist - BDp_mean_null)/(sdP_null) #standardized effect size for BDp
BDp_pvalue<- length(which(BDp_nullDist>=BDp_dist))/(runs) #p value for ses_BDp

Phylo_SCBD<-beta.Phylo$SCBD #contribui??o de esp?cies
Phylo_LCBD<-beta.Phylo$LCBD #contribui??o dos locais
Phylo_LCBD_dist<- beta.Phylo_dist$LCBD
beta.Phylo$p.LCBD #signific?ncia dos locais
beta.Phylo_dist$p.LCBD
which(beta.Phylo$p.LCBD<=0.05) #valores menores ou iguais a 0.05
taxa_signif_PLCBD.Dist<- which(beta.Phylo_dist$p.LCBD<=0.05) #valores menores ou iguais a 0.05

####taxa shuffle for PLCBD####
PLCBD_dist_null<- matrix(NA, nrow= nrow(comu.Ivinhema), ncol= 999, dimnames= list(rownames(comu.Ivinhema), paste("run", 1:999, sep="")),
                         byrow=TRUE)
for(i in 1:999){
  distphylo_null<- picante::taxaShuffle(as.matrix(dist.phylo.Ivinhema)) #taxa shuffle matrix X
  match.namesX<- match(colnames(comu.Ivinhema), colnames(as.matrix(distphylo_null)))
  matrixP_null<- SYNCSA::matrix.p(ifelse(comu.Ivinhema>=1,1,0),as.matrix(dist.phylo.Ivinhema)[match.namesX,match.namesX])$matrix.P #matrix X taxa shuffle
  BDp_nullDist<- beta.div(matrixP_null,method ="percentagedifference",sqrt.D = TRUE,nperm= 2)
  PLCBD_dist_null[,i]<- BDp_nullDist$LCBD
}
PLCBD_taxaShuff<- vector(length= length(taxa_signif_PLCBD.Dist))
for(i in 1:length(taxa_signif_PLCBD.Dist)){
  PLCBD_taxaShuff[i]<- length(which(PLCBD_dist_null[taxa_signif_PLCBD.Dist[i],]>=PLCBD_dist[taxa_signif_PLCBD.Dist[i]]))/999
}
length(which(PLCBD_taxaShuff<=0.05))


length(taxa_signif_PLCBD.Dist)
SCBD.phylo<-cbind(names(beta.Phylo$SCBD),beta.Phylo$SCBD)
write.csv(Phylo_SCBD,"SCBD.phylo.csv" ) #exportando a tabela de contribui??o das esp?cies para div. beta
                                        #filogen?tica
LCBD.phylo<-cbind(names(beta.Phylo$LCBD),beta.Phylo$LCBD)                               
write.csv(LCBD.phylo,"LCBD.phylo.csv")

#calculando por ordem a contribui??o filogen?tica
OCBD_phylo<- read.table("SCBD_OCBD_phylo.txt", header= TRUE)
OCBD_phylo
resul_OCBD_phylo<-sapply(levels(OCBD_phylo$order), function(x){sum(OCBD_phylo$SCBD[OCBD_phylo$order == x])})
#calculando por sub-porcao
LCBD_porc_Phylo<-data.frame(comu.Ivinhema.completa$porcao,Phylo_LCBD)
write.csv(LCBD_porc_Phylo,"LCBD_Phylo.csv")
LCBD_Phylo_porcao<-sapply(levels(LCBD_porc_Phylo$comu.Ivinhema.completa.porcao), function(x){sum(LCBD_porc_Phylo$Phylo_LCBD[LCBD_porc_Phylo$comu.Ivinhema.completa.porcao == x])})
sort(LCBD_Phylo_porcao,TRUE)

windows()
par(bty="l")
plot(1:length(Phylo_SCBD), sort(Phylo_SCBD, decreasing= TRUE), type= "n", xlab= "rank", ylab= "Phylogenetic SCBD", ylim=c(0,0.12),
     xlim=c(0,110),axes=FALSE, xaxs="i", yaxs="i")
axis(1, at=c(seq(0,110,10)))
axis(2,at=c(seq(0,0.12,0.02)))
points(1:length(Phylo_SCBD), sort(Phylo_SCBD, decreasing= TRUE), pch= as.numeric(graph_phylo[,2]))
text((1:3)+20, c(0.103,0.10, 0.09), labels = paste(names(sort(Phylo_SCBD, decreasing= TRUE))[c(1,2,3)]))
lines(1:length(Phylo_SCBD), sort(Phylo_SCBD, decreasing= TRUE))
dev.copy(png,"phylogenetic_SCBD.png")
dev.off()

names_PhyloSCBD<-matrix(names(sort(Phylo_SCBD, decreasing= TRUE)))
deztax_in_phylo<-match(names(dezSCBD_tax),names_PhyloSCBD[,1])
graph_phylo<-matrix(NA, 107, 2, dimnames= list(c(1:length(Func_SCBD)),c("names","Tax/Phylo")))
for(i in 1:107){
	graph_phylo[i,1]<-names(sort(Phylo_SCBD, decreasing= TRUE))[i]
}
graph_phylo[deztax_in_phylo,2]<-21 #circle tax in phylo
graph_phylo[c(-deztax_in_phylo),2]<-19 #filled circle only phylo
graph_phylo[dezfunc_in_phylo,2]<-24 #functional in phylo 


dezfunc_in_phylo<-match(names(dezSCBD_func),names_PhyloSCBD[,1]) 


#######beta funcional######
beta.Func<-beta.div(Y = func.p.Ivinhema, method = "none", save.D = TRUE,clock = TRUE)
beta.Func_dist<- beta.div(Y = func.p.Ivinhema, method= "percentagedifference", sqrt.D = TRUE, nperm = 999)
BDx_dist<- beta.Func_dist$SStotal_BDtotal[2]
XLCBD_dist<- beta.Func_dist$LCBD
BDx_nullDistMatrix<- numeric(length = 999)
for(i in 1:999){
  distrait_null<- picante::taxaShuffle(as.matrix(dist.functional)) #taxa shuffle matrix X
  match.namesX<- match(colnames(comu.Ivinhema), colnames(as.matrix(distrait_null)))
  matrixX_null<- SYNCSA::matrix.p(ifelse(comu.Ivinhema>=1,1,0),as.matrix(dist.functional)[match.namesX,match.namesX])$matrix.P #matrix X taxa shuffle
  BDx_nullDist<- beta.div(matrixX_null,method ="percentagedifference",sqrt.D = TRUE,nperm= 2)
  BDx_nullDistMatrix[i]<- BDx_nullDist$SStotal_BDtotal[2]
}

BDx_mean_null<- apply(as.matrix(BDx_nullDistMatrix),2,mean)
sdX_null<- apply(BDx_nullDist, 2, sd)
ses_BDx<- (BDx_dist - BDx_mean_null)/(sdX_null)
BDx_pvalue<- length(which(BDx_nullDistMatrix>=BDx_dist))/(999)
Func_SCBD<-beta.Func$SCBD #contribui??o de esp?cies
Func_LCBD<-beta.Func$LCBD #contribui??o dos locais
beta.Func$p.LCBD #signific?ncia dos locais
which(beta.Func$p.LCBD<=0.05) #valores menores ou iguais a 0.05
taxa_signif_XLCBD.Dist<- which(beta.Func_dist$p.LCBD<=0.05)

XLCBD_dist_null<- matrix(NA, nrow= nrow(comu.Ivinhema), ncol= 999, dimnames= list(rownames(comu.Ivinhema), paste("run", 1:999, sep="")),
                         byrow=TRUE)
for(i in 1:999){
  distrait_null<- picante::taxaShuffle(as.matrix(dist.functional)) #taxa shuffle matrix X
  match.namesX<- match(colnames(comu.Ivinhema), colnames(as.matrix(distrait_null)))
  matrixX_null<- SYNCSA::matrix.p(ifelse(comu.Ivinhema>=1,1,0),as.matrix(dist.functional)[match.namesX,match.namesX])$matrix.P #matrix X taxa shuffle
  BDx_nullDist<- beta.div(matrixX_null,method ="percentagedifference",sqrt.D = TRUE,nperm= 2)
  XLCBD_dist_null[,i]<- BDx_nullDist$LCBD
}
###calculating taxa shuffle for XLCBD####
XLCBD_taxaShuff<- vector(length= length(taxa_signif_XLCBD.Dist))
for(i in 1:length(taxa_signif_XLCBD.Dist)){
  XLCBD_taxaShuff[i]<- length(which(XLCBD_dist_null[taxa_signif_XLCBD.Dist[i],]>=XLCBD_dist[taxa_signif_XLCBD.Dist[i]]))/999
}
length(which(XLCBD_taxaShuff<=0.05))
####

names.SCBD.phylo<-colnames(phylo.p.Ivinhema[,order(beta.Phylo$SCBD,decreasing=TRUE)])

SCBD.Func<-cbind(names(beta.Func$SCBD),beta.Func$SCBD)
write.csv(Func_SCBD,"SCBD.Func.csv" ) #exportando a tabela de contribui??o das esp?cies para div. beta
                                        #funcional
LCBD.func<-cbind(names(beta.Func$LCBD),beta.Func$LCBD)                               
write.csv(LCBD.func,"LCBD.func.csv")

#calculando por ordem a contribui??o funcional
OCBD_func<- read.table("SCBD_OCBD_func.txt", header= TRUE)
OCBD_func
resul_OCBD_func<-sapply(levels(OCBD_func$order), function(x){sum(OCBD_func$SCBD[OCBD_func$order == x])})

#calculando por sub-porcao
LCBD_porc_Func<-data.frame(comu.Ivinhema.completa$porcao,Func_LCBD)
write.csv(LCBD_porc_Func, "LCBD_Func.csv")
LCBD_Func_porcao<-sapply(levels(LCBD_porc_Func$comu.Ivinhema.completa.porcao), function(x){sum(LCBD_porc_Func$Func_LCBD[LCBD_porc_Func$comu.Ivinhema.completa.porcao == x])})
sort(LCBD_Func_porcao,TRUE)

windows()
plot(beta.Func$LCBD,beta.Phylo$LCBD,pch=19) #comparando os valores de diversidade beta funcional e 
                                           #filogen?tico dos locais
text(beta.Func$LCBD,beta.Phylo$LCBD,labels=comu.Ivinhema.completa$porcao)
summary(lm(beta.Func$LCBD ~ beta.Phylo$LCBD))

windows()
par(bty="l")
plot(1:length(Func_SCBD), sort(Func_SCBD, decreasing= TRUE), type= "n", xlab= "rank", ylab= "Functional SCBD", xaxs="i",yaxs="i",
     xlim=c(0,110),ylim=c(0,0.03), axes=FALSE)
axis(1, at=c(seq(0,110,10)))
axis(2,at=c(seq(0,0.03,0.01)))
points(1:length(Func_SCBD), sort(Func_SCBD, decreasing= TRUE), pch= as.numeric(graph_func[,2]))       
text((1:3)+20, c(0.026,0.023, 0.021), labels = paste(names(sort(Func_SCBD, decreasing= TRUE))[c(1,2,3)]))
lines(1:length(Func_SCBD), sort(Func_SCBD, decreasing= TRUE))
dev.copy(png,"functional_SCBD.png")
dev.off()

dezphylo_in_funct<-match(names(dezSCBD_phylo),names_FunctSCBD[,1])
names_FunctSCBD[dezphylo_in_funct,]<-"red"
names_FunctSCBD[-dezphylo_in_funct,]<-"black"
graph_func<-matrix(NA, 107, 2, dimnames= list(c(1:length(Func_SCBD)),c("names","Tax/Func")))
for(i in 1:107){
	graph_func[i,1]<-names(sort(Func_SCBD, decreasing= TRUE))[i]
}
graph_func[deztax_in_funct,2]<-21
graph_func[c(-deztax_in_funct,-dezphylo_in_funct),2]<-19
graph_func[dezphylo_in_funct,2]<-24
graph_func[common_phylotax_in_func,2]<-22 #phylo e tax em funcional quadrado vazio



match(names(dezSCBD_tax)[match(names(dezSCBD_phylo),names(dezSCBD_tax))],names_FunctSCBD)

common_phylotax_in_func<-match(names(dezSCBD_tax)[match(names(dezSCBD_phylo),names(dezSCBD_tax))],names_FunctSCBD)[!is.na(match(names(dezSCBD_phylo),names(dezSCBD_tax)))]
common_phylo_in_tax<-match(names(dezSCBD_phylo),names(dezSCBD_tax))
graph_func[,1]<-names(sort(Func_SCBD, decreasing= TRUE))
match(names(dezSCBD_phylo),names(dezSCBD_tax))

#correla??o entre as sub-porcoes
LCBD_porcao_all<-data.frame(LCBD_Func_porcao,LCBD_Tax_porcao,LCBD_Phylo_porcao)
windows()
plot(LCBD_porcao_all)
cor(LCBD_porcao_all)
cor.test(LCBD_porcao_all$LCBD_Func_porcao,LCBD_Tax_porcao)
cor.test(LCBD_porcao_all$LCBD_Func_porcao,LCBD_Phylo_porcao)
cor.test(LCBD_porcao_all$LCBD_Phylo_porcao,LCBD_Tax_porcao)
str(SCBD_all)
SCBD_all_matrix<-matrix(c(SCBD_all$Func_SCBD, SCBD_all$Phylo_SCBD, SCBD_all$Tax_SCBD), 
                        nrow= length(SCBD_all$Func_SCBD), ncol=length(SCBD_all), byrow=FALSE,
                        dimnames=list(names(SCBD_all$Func_SCBD),names(SCBD_all)))
cor(as.data.frame(SCBD_all_matrix))
cor.test(SCBD_all$Func_SCBD,SCBD_all$Phylo_SCBD)
cor.test(SCBD_all$Func_SCBD,SCBD_all$Tax_SCBD)
cor.test(SCBD_all$Tax_SCBD,SCBD_all$Phylo_SCBD)

######################################
#grafico de importancia para as ordens
windows()
layout.show(layout(matrix(1:3,1,3)))
par(bty="l")
plot(1:length(resul_OCBD_func), sort(resul_OCBD_func, decreasing= TRUE), type= "n", xlab= "rank", ylab= "Functional OCBD", xlim=c(0,7), 
     ylim=c(0,0.5), xaxs="i", yaxs="i", axes=FALSE)
axis(1, at=c(seq(0,7,1)))
axis(2,at=c(seq(0,0.50,0.10)))
points(1:length(resul_OCBD_func), sort(resul_OCBD_func, decreasing= TRUE), pch= 19)
text((1:3)+0.5, c(0.42,0.23,0.21), labels = paste(names(sort(resul_OCBD_func, decreasing= TRUE))[c(1,2,3)]))
lines(1:length(resul_OCBD_func), sort(resul_OCBD_func, decreasing= TRUE))

windows()
par(bty="l")
plot(1:length(resul_OCBD_phylo), sort(resul_OCBD_phylo, decreasing= TRUE), type= "n", axes=FALSE, xaxs="i", yaxs="i",
     xlab= "rank", ylab= "Phylogenetic SCBD",
     xlim=c(0,7),
     ylim=c(0,0.4))
axis(1, at=c(seq(0,7,1)))
axis(2,at=c(seq(0,0.4,0.1)))
points(1:length(resul_OCBD_phylo), sort(resul_OCBD_phylo, decreasing= TRUE), pch= 19)
text((1:3)+0.5, c(0.35,0.29,0.14), labels = paste(names(sort(resul_OCBD_phylo, decreasing= TRUE))[c(1,2,3)]))
lines(1:length(resul_OCBD_phylo), sort(resul_OCBD_phylo, decreasing= TRUE))

windows()
par(bty="l")
plot(1:length(resul_OCBD_tax), sort(resul_OCBD_tax, decreasing= TRUE), type= "n", xlab= "rank", ylab= "Taxonomic SCBD",
     xaxs="i", yaxs="i",
     ylim=c(0,0.6),xlim=c(0,7))
axis(1, at=c(seq(0,7,1)))
axis(2,at=c(seq(0,0.60,0.10)))
points(1:length(resul_OCBD_tax), sort(resul_OCBD_tax, decreasing= TRUE), pch= 19)
text((1:3)+0.5, c(0.52,0.31,0.06), labels = paste(names(sort(resul_OCBD_tax, decreasing= TRUE))[c(1,2,3)]))
lines(1:length(resul_OCBD_tax), sort(resul_OCBD_tax, decreasing= TRUE))

#calculando equitabilidade da contribui??o para as esp?cies e ordens
source("C:\\Users\\Gabriel\\Documents\\MEGA\\Fun??es\\Camargo Function\\Camargo_function.R", local = TRUE)
#para as ordens
OCBD_all<-list(Func_OCBD=resul_OCBD_func, Phylo_OCBD=resul_OCBD_phylo,Tax_OCBD=resul_OCBD_tax)
eveness_OCBD<-matrix(nrow=3,ncol=1,byrow=TRUE, dimnames = list(c("func_SCBD", "phylo_OCBD", "tax_OCBD"),"OCBD"))
for(i in 1:length(OCBD_all)) eveness_OCBD[i,]<-camargo.eveness(n_spec = OCBD_all[[i]])

#para esp?cies
SCBD_all<-list(Func_SCBD=Func_SCBD, Phylo_SCBD=Phylo_SCBD,Tax_SCBD=Tax_SCBD)
eveness_SCBD<-matrix(nrow=3,ncol=1,byrow=TRUE, dimnames = list(c("func_SCBD", "phylo_SCBD", "tax_SCBD"),"SCBD"))
for(i in 1:length(SCBD_all)) eveness_SCBD[i,]<-camargo.eveness(n_spec = SCBD_all[[i]])

#gerando intervalos de confian?a para equitabilidade
source("C:\\Users\\Gabriel\\Documents\\MEGA\\Fun??es\\bootstrap_mean function\\bootstrap_mean_function.R", local=TRUE)

matrix.P.boot<-vector("list", 999)
for(i in 1:999) {
    matrix.M.boot[[i]]<- matrix.M.stand[sample(1:nrow(matrix.M.stand), replace= TRUE),]
}
names.IV.result.boot<- list(c(1:n.sample), colnames(matrix.M))
IV.result.boot<- matrix(nrow= n.sample, ncol= ncol(matrix.M), dimnames= names.IV.result.boot)


#cinco esp?cies mais importantes para cada dimens?o 
fiveSCBD_tax<-sort(Tax_SCBD, decreasing= TRUE)[1:5]
names(fiveSCBD_tax)
a<-sum(sort(Phylo_SCBD, decreasing= TRUE)[names(fiveSCBD_tax)])/sum(Phylo_SCBD) #propor??o que as 5 esp?cies mais importantes 
                                                                             #taxon?micamente contribuem para filogenia
b<-sum(sort(Func_SCBD, decreasing= TRUE)[names(fiveSCBD_tax)])/sum(Func_SCBD) #propor??o que as 5 esp?cies mais importantes 
                                                                           #taxon?micamente contribuem para filogenia


fiveSCBD_phylo<-sort(Phylo_SCBD, decreasing= TRUE)[1:5]
c<-sum(sort(Tax_SCBD, decreasing= TRUE)[names(fiveSCBD_phylo)])/sum(Tax_SCBD)
d<-sum(sort(Func_SCBD, decreasing= TRUE)[names(fiveSCBD_phylo)])/sum(Func_SCBD)

fiveSCBD_func<-sort(Func_SCBD, decreasing= TRUE)[1:5]
e<-sum(sort(Phylo_SCBD, decreasing= TRUE)[names(fiveSCBD_func)])/sum(Phylo_SCBD)
f<-sum(sort(Tax_SCBD, decreasing= TRUE)[names(fiveSCBD_func)])/sum(Tax_SCBD)

#para dez esp?cies
dezSCBD_tax<-sort(Tax_SCBD, decreasing= TRUE)[1:10]
g<-sum(sort(Phylo_SCBD, decreasing= TRUE)[names(dezSCBD_tax)])/sum(Phylo_SCBD)
h<-sum(sort(Func_SCBD, decreasing= TRUE)[names(dezSCBD_tax)])/sum(Func_SCBD)

dezSCBD_phylo<-sort(Phylo_SCBD, decreasing= TRUE)[1:10]
i<-sum(sort(Tax_SCBD, decreasing= TRUE)[names(dezSCBD_phylo)])/sum(Tax_SCBD)
j<-sum(sort(Func_SCBD, decreasing= TRUE)[names(dezSCBD_phylo)])/sum(Func_SCBD)

dezSCBD_func<-sort(Func_SCBD, decreasing= TRUE)[1:10]
k<-sum(sort(Phylo_SCBD, decreasing= TRUE)[names(dezSCBD_func)])/sum(Phylo_SCBD)
l<-sum(sort(Tax_SCBD, decreasing= TRUE)[names(dezSCBD_func)])/sum(Tax_SCBD)

partitions5_10<-matrix(c(a,b,c,d,e,f,g,h,i,j,k,l),nrow = 6, ncol = 2, dimnames = list(c("Tax|Phylo","Tax|Func","Phylo|Tax","Phylo|Func","Func|Phylo","Func|Tax"),c(paste(c(5,10),"SCBD",sep = "_"))),byrow = FALSE)


####grafico all sites and significance 
LCBD_coordinates<-read.csv("All_LCBD_coordinates.csv", header=TRUE)
source("C:\\Users\\Gabriel\\Google Drive\\Funcoes\\LCBD_function\\BetaDiv_extension_raw.R",local = TRUE)
BetaDiv_Ivinhema<-Beta.div_adapt(Y = comu.Ivinhema, tree = phylo.Ivinhema, trait = traits.Ivinhema, nperm = 999)
length(which(BetaDiv_Ivinhema$p_values_site[1,]<=0.05)) #number of sites with significant values of WLCBD
length(which(BetaDiv_Ivinhema$p_values_site[2,]<=0.05)) #number of sites with significant values of PLCBD
length(which(BetaDiv_Ivinhema$p_values_site[3,]<=0.05)) #number of sites with significant values of XLCBD
length(which(BetaDiv_Ivinhema$p_values_taxa[1,which(BetaDiv_Ivinhema$p_values_site[2,]<=0.05)]<=0.05)) #number of sites that present also significant value of PLCBD for taxa shuffle
length(which(BetaDiv_Ivinhema$p_values_taxa[2,which(BetaDiv_Ivinhema$p_values_site[3,]<=0.05)]<=0.05)) #number of sites that present also significant value of XLCBD for taxa shuffle
cor.test(BetaDiv_Ivinhema$LCBD_all[1,],BetaDiv_Ivinhema$LCBD_all[2,]) #correlation for 
cor.test(BetaDiv_Ivinhema$LCBD_all[1,],BetaDiv_Ivinhema$LCBD_all[3,])
cor.test(BetaDiv_Ivinhema$LCBD_all[2,],BetaDiv_Ivinhema$LCBD_all[3,])
length(BetaDiv_Ivinhema$p_values_site[1,])
min(BetaDiv_Ivinhema$LCBD_all[1,])
min(BetaDiv_Ivinhema$LCBD_all[2,])
min(BetaDiv_Ivinhema$LCBD_all[3,])
max(BetaDiv_Ivinhema$LCBD_all[1,])

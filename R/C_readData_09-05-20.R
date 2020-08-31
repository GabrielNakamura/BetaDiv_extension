###read data####
#spatial data
shp_Hidro<- readOGR("rios_dominio_estado_ANA.shp")
shp_Parana<- readOGR("bacia_parana.shp")

#Ivinhema data
IvinhemaData<- read.csv("IvinhemaData_28-08-2019.csv", header=T, sep=";")
comm_Ivinhema<- IvinhemaData[,-c(1:11)]
coords_pontos<- IvinhemaData[,c("LONGITUDE", "LATITUDE")]
prj <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
Ivinhema_points <- SpatialPointsDataFrame(coords = coords_pontos, data = IvinhemaData, proj4string = prj)
sub <- crop(shp_Hidro, extent(-56, -53, -23.20, -20.5)) #cutting the map to show only ivinhema river basin
quartz()
plot(sub)
axis(1) # showing the axes helps to check whether the coordinates are what you expected
axis(2)
plot(Ivinhema_points, pch=19, col= "red", add= TRUE)

#attributes
atributos_completo<-read.table("~/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/atributos_Ivinhema.txt",header=TRUE)
traits.Ivinhema<- atributos_completo[,c(-1,-8)]
dist.functional<-dist(traits.Ivinhema, method = "euclidean")

#phylogeny
phylo.Ivinhema<-read.tree("~/Google Drive/Manuscritos/Manuscrito Beta_funct_phylo/filogenia.ivinhema.riachos.txt")

#updating phylogenetic hypothesis
list_new<- c("Melanorivulus_ivinhemensis", "Coptodon_rendali", "Trichomycterus_davisi",
             "Farlowella_hahni", "Piabarchus_stramineus", "Astyanax_lacustres", 
             "Megaleporinus_macrocephalus", "Hoplias_missioneira", "Gymnorhamphichthys_britskii",
             "Brachyhypopomus_gauderio", "Curculionichthys_insperatus", "Otothyropsis_marapoama",
             "Otothyropsis_polydon", "Megaleporinus_piavussu")

list_old<- c("Melanorivulus_apiamici", "Tilapia_rendalli", "Trichomycterus_sp.", "Farlowella_amazonum",
             "Bryconamericus_stramineus", "Astyanax_altiparanae", "Leporinus_macrocephalus", "Hoplias_malabaricus",
             "Gymnorhamphichthys_hypostomus", "Brachypopomus_pinnicaudatus", "Hisonotus_insperatus", "Hisonotus_sp.1",
             "Hisonotus_sp.2", "Leporinus_obtusidens")

position<- match(list_old, phylo.Ivinhema$tip.label)
for(i in 1:length(phylo.Ivinhema$tip.label)){
  phylo.Ivinhema$tip.label[position[i]]<- list_new[i] 
}
write.tree(phy = phylo.Ivinhema, file = here::here("data", "processed", "filogenia.ivinhema.riachos.txt"))
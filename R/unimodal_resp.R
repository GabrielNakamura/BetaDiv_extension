
#####modifyed simulate metacommunity######
simulate.metacommunity <- function(tree, n.comm, power, u, niche.breadth= 5, sd.E, sd.X, new.E = FALSE, noise.E = FALSE, noise.X = FALSE){
  n.spp <- length(tree$tip.label)
  a <- rnorm(n.spp, mean = u, sd = niche.breadth)
  h <- runif(n.spp, min = 0, max = 30)
  E <- runif(n.comm, min = 0, max = 100)
  names(E) <- paste("comm", seq_len(n.comm), sep = ".")
  W <- matrix(NA, n.comm, n.spp, dimnames = list(names(E), tree$tip.label))
  if(noise.E){
    EO <- E
    E.Noise <- rnorm(n.comm, mean = 0, sd = sd.E)
    E <- E+E.Noise
    E <- scales::rescale(E, to = c(0, 100))
  }
  X <- ape::rTraitCont(ape::compute.brlen(tree, power = power), model = "BM")
  X <- scales::rescale(X, to = c(-1, 101))
  if(noise.X){
    XO <- X
    X.Noise <- rnorm(n.spp, mean = 0, sd = sd.X)
    X <- X+X.Noise
    X <- scales::rescale(X, to = c(-1, 101))
  }
  for(i in 1:n.comm){
    for (j in 1:n.spp){
      W[i, j] <- rpois(1, h[j]*exp((-((E[i]-X[j])^2))/(2*(a[j]^2))))
    }
  }
  if(noise.E){	
    E <- EO
  }
  if(noise.X){	
    X <- XO
  }
  if(new.E){
    E <- runif(n.comm, min = 0, max = 100)
    names(E) <- paste("comm", seq_len(n.comm), sep = ".")
  }
  res<-list(tree = tree, W = W, E = E, X = X, a= a,h= h)
  return(res)
}

####setting the color pallete#####
#colfunc <- colorRampPalette(c("red", "yellow"))
#heatColors<- colfunc(30)
#quartz()
#plot(rep(1,10),col=colfunc(10),pch=19,cex=3)



####simulating metacommunity with modified function - Beta +#####

power<-2
u<- 5
niche.breadth<- 5
nspp<- 20
ncomm<-30
nperm<-2 #number of permutations to calculate the significance of LCBD component
ntree<-1
runs<-500  #number of permutations for null BD and components

multiphylo<-vector(mode="list",length=ntree) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}
comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity(multiphylo[[i]],n.comm=ncomm,power=power,u= u, niche.breadth= niche.breadth, sd.E = 10,sd.X = 10,noise.E = FALSE,noise.X = FALSE)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}

tree<- multiphylo[[1]]
x<- env.simul[[1]]
X<- trait.simul[[1]]
h<- simul$h
a<- simul$a

tol<- 0.00001
cols <- heat.colors(n = 1000)
rec <- c(X[tree$tip.label], phytools::fastAnc(tree, X))
names(rec)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
rec_rec<- matrix(rec[tree$edge], nrow(tree$edge), 2)
rec_values <- rowMeans(rec_rec)
xlims <- range(rec_values) + c(-tol, tol)
breaks <- 0:1000/1000 * (xlims[2] - xlims[1]) + xlims[1]
whichColor <- function(p, cols, breaks) {
  i <- 1
  while (p >= breaks[i] && p > breaks[i + 1]) i <- i + 
      1
  cols[i]
}
colors <- sapply(rec_values, whichColor, cols = cols, breaks = breaks)
colors_spp<- colors[match(1:length(tree$tip.label),tree$edge[,2])]

#####ploting species response curves#####
quartz()
phytools::plotBranchbyTrait(tree, x= X, mode="tips", palette="heat.colors")

quartz()
plot(x = 0, y = 0, xlab = "Gradient", ylab = "Abundance", 
     xlim = c(-20, 100), ylim = c(0, max(comm.simul[[1]])), main = "unimodal species curves with power 1 phylo", 
     type = "n")
for (i in 1:nspp) {
  species <- function(x) (h[i]*exp((-((x-X[i])^2))/(2*(a[i]^2))))
  curve(species, from = -20, to = 100, add = TRUE, col= colors_spp[i])
}


####simulating metacommunity with modified function - Beta, power 0.001#####
power<-0.0001
u<- 5
nspp<- 30
ncomm<-10
nperm<-2 #number of permutations to calculate the significance of LCBD component
ntree<-1
runs<-999  #number of permutations for null BD and components

multiphylo<-vector(mode="list",length=ntree) 
for(i in 1:length(multiphylo)){
  multiphylo[[i]]<-geiger::sim.bdtree(b=0.1,d=0,n=nspp,extinct=FALSE)
}
comm.simul<-vector(mode="list",length=length(multiphylo))
trait.simul<-vector(mode="list",length=length(multiphylo))
env.simul<-vector(mode="list",length=length(multiphylo))
for(i in 1:length(multiphylo)){
  simul<-simulate.metacommunity(multiphylo[[i]],n.comm=ncomm,power=power,u= u, niche.breadth= niche.breadth, sd.E = 10,sd.X = 10,noise.E = TRUE,noise.X = TRUE)
  comm.simul[[i]]<-simul$W
  trait.simul[[i]]<-simul$X
  env.simul[[i]]<-simul$E
}


x<- env.simul[[1]]
X<- trait.simul[[1]]
h<- simul$h
a<- simul$a

#####ploting species response curves#####
tol<- 0.00001
cols <- heat.colors(n = 1000)
x <- c(X[tree$tip.label], phytools::fastAnc(tree, X))
names(x)[1:length(tree$tip.label)] <- 1:length(tree$tip.label)
XX <- matrix(x[tree$edge], nrow(tree$edge), 2)
x <- rowMeans(XX)
xlims <- range(x) + c(-tol, tol)
breaks <- 0:1000/1000 * (xlims[2] - xlims[1]) + xlims[1]
whichColor <- function(p, cols, breaks) {
  i <- 1
  while (p >= breaks[i] && p > breaks[i + 1]) i <- i + 
      1
  cols[i]
}
colors <- sapply(x, whichColor, cols = cols, breaks = breaks)


quartz()
phytools::plotBranchbyTrait(multiphylo[[1]], x= X, mode="tips", palette="heat.colors")
quartz()
plot(x = 0, y = 0, xlab = "Gradient", ylab = "Abundance", 
     xlim = c(0, 100), ylim = c(0, max(comm.simul[[1]])), main = "unimodal species curves with power 0.0001 phylo", 
     type = "n")
for (i in 1:nspp) {
  species <- function(x) (h[i]*exp((-((x-X[i])^2))/(2*(a[i]^2))))
  curve(species, from = 0, to = 100, add = TRUE, col= colors[i])
}




simulate.metacommunity2 <- function(tree, n.comm, power, u, sd.E, sd.X, new.E = FALSE, new.X= FALSE, noise.E = FALSE, noise.X = FALSE){
  n.spp <- length(tree$tip.label)
  a <- rnorm(n.spp, mean = u, sd = 10)
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
  if(new.X){
    X<-runif(n = length(tree$tip.label),min = -1,max = 101)
    names(X)<-paste("s",1:length(tree$tip.label),sep="")
  }
  res<-list(tree = tree, W = W, E = E, X = X)
  return(res)
}
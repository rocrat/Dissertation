library(compositions)


gm_mean = function(x){
  return(sum(log(x)) / length(x))
}

MFclr <- function(x){#apply to a vector
  delta <- 0.55/sum(x)# .55 times the smallest detectable value (1 read) after closure
  tdelta <- sum(x == 0) * delta
  cx <- clo(x)
  cxt <- ifelse(cx == 0, delta, cx * (1-tdelta))
  g <- gm_mean(cxt)
  ccxt <- log(cxt) - g
  return(ccxt)
}

MFalr <- function(x, ind){#apply to a vector
  browser()
  delta <- 0.55/sum(x)# .55 times the smallest detectable value (1 read) after closure
  tdelta <- sum(x == 0) * delta
  cx <- clo(x)
  D <- cx[ind]
  cxt <- ifelse(cx == 0, delta, cx * (1-tdelta))
  acxt <- log(cxt) - log(D)
  return(acxt)
}

#Matrix manipulation
Id <- matrix(0, nrow = 3, ncol = 3)
diag(Id) <- 1
j <- matrix(1, 3, 1)
N <- Id + j%*%t(j)
Jd <- matrix(1, 3, 3)
F <- cbind(Id, -1*j)
D <- 4
ID <- matrix(0, 4, 4)
diag(ID) <- 1
d <- 3
jD <- matrix(1, 4, 1)
JD <- matrix(1, 4, 4)
GD <- ID - D^-1 * JD
Hd <- Id + Jd

Ninv <- Id - 1/d * j%*%t(j)
ND <- ID + jD%*%t(jD)
NDinv <- ID - 1/D * jD%*%t(jD)
NDinv <- solve(NDinv)
NDinv%*%NDinv == NDinv
#Make a composition
X <- matrix(NA, nrow = 4, ncol = 8)
for(i in 1:ncol(X)){
  X[, i] <- rpois(4, 10)
}

Xclr <- apply(X, 2, MFtrans)
gamma <- cov(t(Xclr))

rowSums(gamma)
gamma%*%JD
GD%*%gamma
gamma
GD%*%GD

Xalr <- apply(X, 2, MFalr, ind = nrow(X))

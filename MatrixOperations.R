library(compositions)
library(magrittr)

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
F. <- cbind(Id, -1*j)
D <- 4
ID <- matrix(0, 4, 4)
diag(ID) <- 1
d <- 3
jD <- matrix(1, 4, 1)
JD <- matrix(1, 4, 4)
GD <- ID - D^-1 * JD
Hd <- Id + Jd
Hinv <- solve(Hd)

Ninv <- Id - 1/d * j%*%t(j)
ND <- ID + jD%*%t(jD)
NDinv <- ID - 1/D * jD%*%t(jD)
NDinv <- solve(NDinv)
NDinv%*%NDinv == NDinv

C <- t(X) %*% jD
Cm <- matrix(0, ncol = 8, nrow = 8)
diag(Cm) <- C^-1




#Make a composition
X <- matrix(NA, nrow = 4, ncol = 8)
for(i in 1:ncol(X)){
  X[, i] <- rpois(4, 10)
}

Xclr <- apply(X, 2, MFclr)
Xalr <- t(X) %>% #Transpose for alr function
  clo() %>% #close composition
  alr() %>% #ALR transform
  t() #Return to original position
  
gamma <- cov(t(Xclr))
sigma <- cov(t(Xalr))

gamma.new <- t(F.) %*% Hinv %*% sigma %*% Hinv %*% F.

all.equal(gamma, gamma.new)

rowSums(gamma)
gamma%*%JD
GD%*%gamma
gamma
GD%*%GD

alrPCR <- prcomp(t(Xalr))
clrPCR <- prcomp(t(Xclr))

# Closure Operation
C <- t(X) %*% jD
Cm <- matrix(0, ncol = 8, nrow = 8)
diag(Cm) <- C

Xclo <- X %*% Cm
colSums(Xclo)
t(clo(t(X)))

Xclo2 <- Cm %*% t(X)

In <- matrix(0, ncol = 8, nrow = 8)
diag(C^-1)
Xclo2 <- Xclo %*% Cm

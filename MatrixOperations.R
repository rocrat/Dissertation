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
JD <- matrix(1, 4, 4)
GD <- ID - D^-1 * JD
Hd <- Id + Jd


#Make a composition
X <- matrix(NA, nrow = 4, ncol = 8)
for(i in 1:ncol(X)){
  X[, i] <- rpois(4, 10)
}



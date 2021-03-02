## Dictionary ##
dictionary <- "D:/Program Files/RStudio/projects/LDA-implementation/F-B Model/dictionary.txt"
w <- readLines(dictionary)
V <- unlist(strsplit(w," "))
alpha_mu <- 0.05
beta_b <- 0.01 
beta_f <- 0.01
alpha_theta <- 0.09
alpha_lambda <- 0.5
## initializing ##
B <- 5                                             # number of background documents
K_B <- 2                                           # number of background topics
T <- 6                                             # number of foreground documents
K_T <- 3                                           # number of foreground topics
library(MCMCpack)
## Creating background documents matrix
n_b <- c(45,50,77,38,77)                           # number of words in each document
Background <- matrix(nrow = B, ncol = max(n_b))
## Creating foreground documents matrix
n_t <- c(94,52,101,106,112,68)                     # number of words in each document
Foreground <- matrix(nrow = T, ncol = max(n_t))
## Creating Mu-B Matrix
Mu_B <- matrix(nrow = B, ncol = K_B)
Mu_B<- rdirichlet(B, rep(alpha_mu,K_B))
## Creating Mu-T Matrix
Mu_T <- matrix(nrow = T, ncol = K_B)
Mu_T <- rdirichlet(T, rep(alpha_mu,K_B))
## Creating Phi-B Matrix
Phi_B <- matrix(nrow = K_B, ncol = length(V))
Phi_B <- rdirichlet(K_B, rep(beta_b,length(V)))
## Creating Theta-T matrix
Theta_T <- matrix(nrow = T, ncol = K_T)
Theta_T<- rdirichlet(T, rep(alpha_theta,K_T))
## Creating Phi-T matrix
Phi_F <- matrix(nrow = K_T, ncol = length(V))
Phi_F <- rdirichlet(K_T, rep(beta_f,length(V)))
## Creating type decision distribution 
lambda <- rdirichlet(T,rep(alpha_lambda,2))
## Generating background-words
for(b in 1:B){
  Mu_b <- Mu_B[b,]
  for(i in 1:n_b[b]){
    Z_b_i <-sample(c(1:K_B),1,replace = TRUE,prob = Mu_b)
    W_b_i <-sample(c(1:length(V)),1,replace = TRUE,prob = Phi_B[Z_b_i,])
    Background[b,i]<-V[W_b_i]
  }
}
## Generating foreground-words
for(t in 1:T){
  lambda_t <- lambda[t,]
  for(i in 1:n_t[t]){
    y_t_i <-sample(c(0:1),1,replace = TRUE,prob = lambda_t)
    if(y_t_i==0){
      Theta_t <- Theta_T[t,] 
      Z_t_i <-sample(c(1:K_T),1,replace = TRUE,prob = Theta_t)
      W_t_i <-sample(c(1:length(V)),1,replace = TRUE,prob = Phi_F[Z_t_i,])
      Foreground[t,i]<-V[W_t_i]
    }else{
      Mu_t <- Mu_T[t,] 
      Z_b_i <-sample(c(1:K_B),1,replace = TRUE,prob = Mu_t)
      W_b_i <-sample(c(1:length(V)),1,replace = TRUE,prob = Phi_B[Z_b_i,])
      Foreground[t,i]<-V[W_b_i]
    }
  }
}

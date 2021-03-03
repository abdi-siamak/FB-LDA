#################################################
# (C) Siamak Abdi - 2020
# Implementation of the collapsed Gibbs sampling for Foreground-Background Latent Dirichlet Allocation (FB-LDA)
# model in R, as described in "Interpreting the Public Sentiment Variations on Twitter". Shulong Tan et al.
#################################################
library(purrr)
KB = 2                                                    # number of Background topics
KF = 2                                                   # number of Foreground topics
alpha = 0.05
beta = 0.05
iteration = 5000                                      
burn = 2500
#memory.limit(size=16279)
options(scipen=999)
#################################################         # creating (topics/ percent) Matrix
na <-c()
topics <- matrix(0,nrow = W, ncol = (KB+KF))
for(k in 1:KB){
  na <-append(na,paste0("KB",k),after = length(na))
}
for(k in 1:KF){
  na <-append(na,paste0("KF",k),after = length(na))
}
colnames(topics)<- na
percent <- matrix(nrow = M, ncol = (KB+KF))
colnames(percent)<- na
#################################################         # finding duplicate words
R <- which(duplicated(gibbs[,1]) | duplicated(gibbs[,1], fromLast = TRUE))  # repeat words
U <- setdiff(R,which(duplicated(gibbs[,1])))                                # base words that have repeat
R_N <- setdiff(which(gibbs[,1]==gibbs[,1]),R)             # words that have not repeat
V <- gibbs[union(U,R_N),1]                                # Dictionary
#################################################         # row names for topics matrix
rownames(topics)<- gibbs[1:W]
#################################################         # initialize of Z & Y
for(z in 1:W){
  if(gibbs[z,5]==1){                                         # Y=1 -> Background  Y=0 -> Foreground
    gibbs[z,4] <-1                                          
    gibbs[z,3] <- rdunif(n=1, b=KB, a = 1) 
  }else if(gibbs[z,5]==2){
    gibbs[z,4] <- rdunif(n=1, b=1, a = 0)
    if(gibbs[z,4]==1){
      gibbs[z,3] <- rdunif(n=1, b=KB, a = 1)
    }else{
      gibbs[z,3] <- rdunif(n=1, b=KF, a = 1)
    }
  }
}
#################################################         # Estimated 
PHI_B <- matrix(0,nrow = KB, ncol = length(V))
PHI_F <- matrix(0,nrow = KF, ncol = length(V))
THETA_B <- matrix(0,nrow = MB, ncol = KB)
THETA_F <- matrix(0,nrow = MF, ncol = KF)
colnames(PHI_B)<- V
colnames(PHI_F)<- V

#################################################         # Theta matrix
Theta_B <- matrix(0,nrow = MB, ncol = KB)
Theta_F_1 <- matrix(0,nrow = MF, ncol = KB)
rownames(Theta_F_1)<- c((MB+1):(MB+MF))
Theta_F_2 <- matrix(0,nrow = MF, ncol = KF)
rownames(Theta_F_2)<- c((MB+1):(MB+MF))
##################################################        # Phi matrix
Phi_B <- matrix(0,nrow = KB, ncol = length(V)) 
colnames(Phi_B)<- V
Phi_F <- matrix(0,nrow = KF, ncol = length(V)) 
colnames(Phi_F)<- V
##################################################        # M_t matrix
M_t <- matrix(0,nrow = MF, ncol = 2) 
colnames(M_t)<- c("1","0")
rownames(M_t)<- c((MB+1):(MB+MF))
##################################################        # initialize Phi and Theta and N_R and M_t matrix
N_R_B <- array(0,dim = KB)
N_R_F <- array(0,dim = KF)
for(x in 1:W){
  if(gibbs[x,5]==1){
    Theta_B[as.integer(gibbs[x,2]),as.integer(gibbs[x,3])]<-Theta_B[as.integer(gibbs[x,2]),as.integer(gibbs[x,3])]+1
    Phi_B[as.integer(gibbs[x,3]),gibbs[x,1]] <- Phi_B[as.integer(gibbs[x,3]),gibbs[x,1]] + 1
    N_R_B[as.integer(gibbs[x,3])] <- N_R_B[as.integer(gibbs[x,3])] + 1
  }else if(gibbs[x,5]==2){
    if(gibbs[x,4]==1){
       Theta_F_1[gibbs[x,2],as.integer(gibbs[x,3])]<-Theta_F_1[gibbs[x,2],as.integer(gibbs[x,3])]+1
       Phi_B[as.integer(gibbs[x,3]),gibbs[x,1]] <- Phi_B[as.integer(gibbs[x,3]),gibbs[x,1]] + 1
       N_R_B[as.integer(gibbs[x,3])] <- N_R_B[as.integer(gibbs[x,3])] + 1
       M_t[gibbs[x,2],"1"]<-M_t[gibbs[x,2],"1"]+1
    }else if(gibbs[x,4]==0){
       Theta_F_2[gibbs[x,2],as.integer(gibbs[x,3])]<-Theta_F_2[gibbs[x,2],as.integer(gibbs[x,3])]+1
       Phi_F[as.integer(gibbs[x,3]),gibbs[x,1]] <- Phi_F[as.integer(gibbs[x,3]),gibbs[x,1]] + 1
       N_R_F[as.integer(gibbs[x,3])] <- N_R_F[as.integer(gibbs[x,3])] + 1
       M_t[gibbs[x,2],"0"]<-M_t[gibbs[x,2],"0"]+1
    }
  }
}
##################################################       ## iterations ##
##################################################
n_p_B <- numeric(KB)
n_r_B <- numeric(KB)
n_t_B <- numeric(KB)
n_p_F <- numeric(KF)
n_r_F <- numeric(KF)
n_t_F_2 <- numeric(KF)
n_t_F_1 <- numeric(KB)
m_t <- numeric(2)
until <- 0
for(i in 1:iteration){
  until <- until + 1
  #startTime <- Sys.time()
  Update_1 <- FALSE
  Update_2 <- FALSE
  Update_3 <- FALSE
  for(z in 1:W){  
    #startTime <- Sys.time()
    if(z==1){                                                          ### loading new data
      n_t_B <-Theta_B[as.integer(gibbs[z,2]),] 
      n_t_B[as.integer(gibbs[z,3])] <- n_t_B[as.integer(gibbs[z,3])] - 1
      n_r_B <- N_R_B[]
      n_r_B[as.integer(gibbs[z,3])] <- n_r_B[as.integer(gibbs[z,3])]  - 1
      n_p_B <-Phi_B[,gibbs[z,1]]
      n_p_B[as.integer(gibbs[z,3])] <- n_p_B[as.integer(gibbs[z,3])] - 1 
      
      n_t_B <- n_t_B + alpha
      n_r_B <- n_r_B +((length(V))*beta)
      n_p_B <- n_p_B + beta
    }else if(z!=1){
      if(gibbs[z,5]==1){
        n_t_B <-Theta_B[as.integer(gibbs[z,2]),] 
        n_t_B[as.integer(gibbs[z,3])] <- n_t_B[as.integer(gibbs[z,3])] - 1
        n_t_B <- n_t_B + alpha
      }else if(gibbs[z,5]==2){
        n_t_F_2 <-Theta_F_2[gibbs[z,2],] 
        n_t_F_1 <-Theta_F_1[gibbs[z,2],]
        if(gibbs[z,4]==0){
          n_t_F_2[as.integer(gibbs[z,3])] <- n_t_F_2[as.integer(gibbs[z,3])] - 1
        }
        if(gibbs[z,4]==1){
          n_t_F_1[as.integer(gibbs[z,3])] <- n_t_F_1[as.integer(gibbs[z,3])] - 1
        }
        n_t_F_2 <- n_t_F_2 + alpha
        n_t_F_1 <- n_t_F_1 + alpha
      }
      n_r_B <- N_R_B[]
      if(gibbs[z,4]==1){
        n_r_B[as.integer(gibbs[z,3])] <- n_r_B[as.integer(gibbs[z,3])]  - 1
      }
      n_r_F <- N_R_F[]
      if(gibbs[z,4]==0){
        n_r_F[as.integer(gibbs[z,3])] <- n_r_F[as.integer(gibbs[z,3])] - 1
      }
      n_p_B <-Phi_B[,gibbs[z,1]]
      if(gibbs[z,4]==1){
      n_p_B[as.integer(gibbs[z,3])] <- n_p_B[as.integer(gibbs[z,3])] - 1
      }
      n_p_F <-Phi_F[,gibbs[z,1]]
      if(gibbs[z,4]==0){
        n_p_F[as.integer(gibbs[z,3])] <- n_p_F[as.integer(gibbs[z,3])] - 1
      }
      if(gibbs[z,5]==2){
        m_t <- M_t[gibbs[z,2],]
        m_t[gibbs[z,4]] <- m_t[gibbs[z,4]] - 1
        
        m_t_KF <- m_t[2] + (KF*alpha) 
        m_t_KB <- m_t[1] + (KB*alpha)
        m_t_KF <- rep(m_t_KF,KF)
        m_t_KB <- rep(m_t_KB,KB)
        m_t <- m_t + alpha
        m_t_1 <- rep(m_t[1],KB)
        m_t_0 <- rep(m_t[2],KF)
      }
      n_r_B <- n_r_B +((length(V))*beta)
      n_r_F <- n_r_F +((length(V))*beta)
      n_p_B <- n_p_B + beta
      n_p_F <- n_p_F + beta
    }
    if(gibbs[z,5]==1){
      p_Z <- (n_t_B) * (n_p_B/n_r_B)
      #p_Z <- p_Z/sum(p_Z)                                   # normalizing
      prob <- sample(c(1:KB),1,prob = p_Z)                   # sampling (Draw) - new Z & new Y (prob & y_z_t)
      #prob <-which(p_Z==max(p_Z))
      y_z_t <- 1
    }else if(gibbs[z,5]==2){
      P_Z_0 <- (n_p_F/n_r_F)*(n_t_F_2/m_t_KF)*m_t_0
      P_Z_1 <- (n_p_B/n_r_B)*(n_t_F_1/m_t_KB)*m_t_1
      p_Z <- c(P_Z_0,P_Z_1)
      prob <- sample(c(1:(KB+KF)),1,prob = p_Z)
      #prob <-which(p_Z==max(p_Z))
      #if(length(prob)>1){
      #  prob <-prob[1]
      #}
      if(prob<=KF){
        y_z_t <- 0
      }else if(prob>KF){
        prob <- prob - KF
        y_z_t <- 1
      }
    }

    if((as.integer(gibbs[z,3])!=prob) | (gibbs[z,4]!=y_z_t)){                 # updating Theta / Phi / N_R / M_t 
      if(y_z_t==1){
        if(gibbs[z,5]==1){
            Update_1 <-TRUE
            Theta_B[as.integer(gibbs[z,2]),prob]<-Theta_B[as.integer(gibbs[z,2]),prob]+1
            Theta_B[as.integer(gibbs[z,2]),as.integer(gibbs[z,3])]<-Theta_B[as.integer(gibbs[z,2]),as.integer(gibbs[z,3])]-1
            Phi_B[prob,gibbs[z,1]]<-Phi_B[prob,gibbs[z,1]]+1
            Phi_B[as.integer(gibbs[z,3]),gibbs[z,1]]<-Phi_B[as.integer(gibbs[z,3]),gibbs[z,1]]-1
            N_R_B[prob]<-N_R_B[prob]+1
            N_R_B[as.integer(gibbs[z,3])]<-N_R_B[as.integer(gibbs[z,3])]-1
         }else if(gibbs[z,5]==2){
            Update_2 <-TRUE
            Theta_F_1[gibbs[z,2],prob]<-Theta_F_1[gibbs[z,2],prob]+1
            if(gibbs[z,4]==y_z_t){
              Theta_F_1[gibbs[z,2],as.integer(gibbs[z,3])]<-Theta_F_1[gibbs[z,2],as.integer(gibbs[z,3])]-1
            }else if(gibbs[z,4]!=y_z_t){
              Theta_F_2[gibbs[z,2],as.integer(gibbs[z,3])]<-Theta_F_2[gibbs[z,2],as.integer(gibbs[z,3])]-1
            }
            Phi_B[prob,gibbs[z,1]]<-Phi_B[prob,gibbs[z,1]]+1
            if(gibbs[z,4]==y_z_t){
              Phi_B[as.integer(gibbs[z,3]),gibbs[z,1]]<-Phi_B[as.integer(gibbs[z,3]),gibbs[z,1]]-1
            }else if(gibbs[z,4]!=y_z_t){
              Phi_F[as.integer(gibbs[z,3]),gibbs[z,1]]<-Phi_F[as.integer(gibbs[z,3]),gibbs[z,1]]-1
            }
            N_R_B[prob]<-N_R_B[prob]+1
            if(gibbs[z,4]==y_z_t){
              N_R_B[as.integer(gibbs[z,3])]<-N_R_B[as.integer(gibbs[z,3])]-1
            }else if(gibbs[z,4]!=y_z_t){
              N_R_F[as.integer(gibbs[z,3])]<-N_R_F[as.integer(gibbs[z,3])]-1
            }
            if(gibbs[z,4]!=y_z_t){
              M_t[gibbs[z,2],1]<-M_t[gibbs[z,2],1]+1
              M_t[gibbs[z,2],gibbs[z,4]]<-M_t[gibbs[z,2],gibbs[z,4]]-1
            }
        }
      }else if(y_z_t==0){
        Update_3 <-TRUE
        Theta_F_2[gibbs[z,2],prob]<-Theta_F_2[gibbs[z,2],prob]+1
        if(gibbs[z,4]==y_z_t){
          Theta_F_2[gibbs[z,2],as.integer(gibbs[z,3])]<-Theta_F_2[gibbs[z,2],as.integer(gibbs[z,3])]-1
        }else if(gibbs[z,4]!=y_z_t){
          Theta_F_1[gibbs[z,2],as.integer(gibbs[z,3])]<-Theta_F_1[gibbs[z,2],as.integer(gibbs[z,3])]-1
        }
        Phi_F[prob,gibbs[z,1]]<-Phi_F[prob,gibbs[z,1]]+1
        if(gibbs[z,4]==y_z_t){
          Phi_F[as.integer(gibbs[z,3]),gibbs[z,1]]<-Phi_F[as.integer(gibbs[z,3]),gibbs[z,1]]-1
        }else if(gibbs[z,4]!=y_z_t){
          Phi_B[as.integer(gibbs[z,3]),gibbs[z,1]]<-Phi_B[as.integer(gibbs[z,3]),gibbs[z,1]]-1
        }
        N_R_F[prob]<-N_R_F[prob]+1
        if(gibbs[z,4]==y_z_t){
          N_R_F[as.integer(gibbs[z,3])]<-N_R_F[as.integer(gibbs[z,3])]-1
        }else if(gibbs[z,4]!=y_z_t){
          N_R_B[as.integer(gibbs[z,3])]<-N_R_B[as.integer(gibbs[z,3])]-1
        }
        if(gibbs[z,4]!=y_z_t){
          M_t[gibbs[z,2],2]<-M_t[gibbs[z,2],2]+1
          M_t[gibbs[z,2],gibbs[z,4]]<-M_t[gibbs[z,2],gibbs[z,4]]-1
        }
       }
      gibbs[z,3] <- prob                                     ### update z
      gibbs[z,4] <- y_z_t                                    ### update y
    }
    if(i > burn){                                        # update topics matrix
      v_4 <- vector("numeric", KB+KF)
      if(y_z_t==1){
        v_4[as.integer(gibbs[z,3])]<- 1
      }else if(y_z_t==0){
        v_4[(KB+as.integer(gibbs[z,3]))]<- 1
      }
      topics[z,] <-topics[z,] + v_4
    }
    #endTime <- Sys.time() - startTime
    #cat("\r \n ", z)
  } 
  if(i > burn){
    #################################################           # estimating PHI matrix (parameter estimation)
    for(s in 1:KB){
      for(r in 1:length(V)){
        PHI_B[s,r]<-PHI_B[s,r]+(Phi_B[s,r]+beta)/(N_R_B[s]+(length(V)*beta))
      }
    }
    for(s in 1:KF){
      for(r in 1:length(V)){
        PHI_F[s,r]<-PHI_F[s,r]+(Phi_F[s,r]+beta)/(N_R_F[s]+(length(V)*beta))
      }
    }
    #################################################           # estimating THETA matrix (parameter estimation)          
    for(j in 1:MB){
      for(s in 1:KB){
        THETA_B[j,s]<-THETA_B[j,s]+(Theta_B[j,s]+alpha)/((docs_len[j])+(KB*alpha))
      }
    }
    for(j in 1:MF){
      for(s in 1:KF){
        THETA_F[j,s]<-THETA_F[j,s]+(Theta_F_2[j,s]+alpha)/((M_t[2])+(KF*alpha))
      }
    }
    #################################################
  }
  #cat("\r learning: %",(i/iteration)*100)
  #endTime <- Sys.time() - startTime
  if(until == 100){
    cat("\r \n ", i)
    until <- 0
  }
}
##################################################          # get word topic 
for(r in 1:W){                                           
  count <- names(which(topics[r,]==max(topics[r,])))
  if(length(count)>1){
    gibbs[r,6] <- sample(count,1,prob = rep((1/length(count)),length(count)))
  }else{
    gibbs[r,6] <- count
  }
}
##################################################          # Topic percentage of each document
for(p in 1:M){
  KK <- c()
  x<-subset(gibbs, gibbs[,2]==p)                            # subset of gibbs matrix
  for(k in na){
    KK[k]<-length(which(x[,6] == k))
  }
  percent[p,] <- (KK/sum(KK))*100
}
#################################################           # estimating PHI matrix (parameter estimation)
for(i in 1:KB){
  for(r in 1:length(V)){
    PHI_B[i,r]<-PHI_B[i,r]/(iteration-burn)
  }
}
for(i in 1:KF){
  for(r in 1:length(V)){
    PHI_F[i,r]<-PHI_F[i,r]/(iteration-burn)
  }
}
#################################################           # estimating THETA matrix (parameter estimation)          
for(j in 1:MB){
  for(i in 1:KB){
    THETA_B[j,i]<-THETA_B[j,i]/(iteration-burn)
  }
}
for(j in 1:MF){
  for(i in 1:KF){
    THETA_F[j,i]<-THETA_F[j,i]/(iteration-burn)
  }
}

## output matrix: gibbs, percent, topics, PHI, THETA



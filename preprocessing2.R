#knitr::opts_chunk$set(echo = TRUE)
#require(pryr)
#require(gpuR)
#listContexts()
####################################################
Path <- setwd("D:/Program Files/RStudio/projects/F-B Model/documents")
#get listing of folders
Foldernames <- list.files(getwd())
docs_len <- c()
words <-c()
Docs <-c()
ww <- c()
dd <- c()
cc <- c()
tt <- c()
n_c<- c(0)
i <- 1
j <- 1
for(f in Foldernames){
  Path <- setwd(paste0("D:/Program Files/RStudio/projects/F-B Model/documents/",f))
  Path <- setwd(paste0("D:/Program Files/RStudio/projects/F-B Model/documents/",f))
  #load files into corpus
  #get listing of .txt files in Folder
  Filenames <- list.files(getwd(),pattern="*.txt")
  Docs <- append(Docs,length(Filenames),after = length(Docs)) 
  for (d in Filenames) {
    w <- scan(d,what=" ")
    docs_len[i] <- length(w)                               # length of each documnet
    ww <- append(ww,w,after = length(ww))                  # all of the words in all datasets
    tt <- append(tt,w,after = length(tt))                  # all of the words in each dataset
    dd <- append(dd,rep(i,length(w)),after = length(dd))   # assign document name of each word
    i <- i+1
  }
  words[j] <- length(tt)                                          # number of words in each dataset
  cc <- append(cc,rep(j,length(tt)),after = length(cc))
  tt <- NULL
  n_c <- append(n_c,i-1,after = length(n_c))
  j <- j+1
}
## creating gibbs matrix
W <- length(ww)
gibbs <- matrix(nrow = W, ncol = 6)
colnames(gibbs)<- c("words","document","Z*","Y*","DataSet","word Topic")
gibbs[,1] <- ww                                         
gibbs[,2] <- dd
gibbs[,5] <- cc                                 # 1 = Background     2 = Foreground
cat("W:\n ",length(ww))
cat("WB:\n ",words[1])
cat("WF:\n ",words[2])
cat("M:\n ",i-1)
W <- length(ww)                                 # number of words in both data sets
M <- i-1                                        # number of all documents
WB <- words[1]                                  # number of words in Background data set
WF <- words[2]                                  # number of words in Foreground data set
MB <- Docs[1]                                   # number of documents in Background data set
MF <- Docs[2]                                   # number of documents in Foreground data set

###Output: gibbs, W, M, WB, WF, MB, MF

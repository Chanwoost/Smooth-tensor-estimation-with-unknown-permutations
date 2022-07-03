source("functions_asym.R")
## The following codes are for Table 3 and Table S2


## Getting the grid of the group number (k1,k2,k3)
kind = NULL
for(i in 6:9){
  for(j in 2:9){
    for(k in 3:9){
      if(i*j<k | j*k<i |i*k<j){
        print("hi")
      }else{
        kind = rbind(kind,c(i,j,k))
      }
    }
  }
}

nrow(kind)

thre = 10+ 2*(1:22)
threind =  expand.grid(1:3,thre)


## Getting the optimal group numbers (k1,k2,k3)
optimalk = as.data.frame(matrix(nrow = 680, ncol = 6))
names(optimalk) = c("model","method","k1","k2","k3","MSE")
optimalk$model = rep(c("model1","model2","model3","model4","model5"),each =136)
optimalk$method = rep(c("Borda","LSE"), 340)

## Getting the optimal threshold t by grid search
optimalt = as.data.frame(matrix(nrow = 330, ncol = 4))
names(optimalt) = c("model","mode","threshold","MSE")
optimalt$model =rep(c("model1","model2","model3","model4","model5"),each = 66)
optimalt[,2:3] = rbind(threind,threind,threind,threind,threind)



## MSE for Borda count and LSE according to (k1,k2,k3)
MSE1 = NULL
kindex = NULL
for(m in 1:5){
  for(s in 1:66){
    kvec = kind[s,]
    kindex = rbind(kindex,kvec,kvec)
    s1 = simulation_asym(50,40,30,mode = m,signal_level = 1)
    MSE1 = c(MSE1,mean((Borda2_asym(s1$observe,2,kvec)$Theta-s1$signal)^2))
    MSE1 = c(MSE1,mean((LSE_asym(s1$observe,kvec,mode = 3)-s1$signal)^2))
    print(paste("s_",s,"_ended",sep = ""))
  }
}

optimalk[,3:5] = kindex
optimalk$MSE = MSE1

## MSE for Spectral method according to threshold (t)
MSE2 = NULL
for(m in 1:5){
  for(s in 1:66){
    s1 = simulation_asym(50,40,30,mode = m,signal_level = 1)
    MSE2 = c(MSE2,mean((Spectral(s1$observe,threind[s,1],setdiff(1:3,threind[s,1]),threshold = threind[s,2])-s1$signal)^2))
    print(paste("s_",s,"_ended",sep = ""))
  }
}

optimalt$MSE = MSE2


## Getting Table S2 for Borda count and Spectral
optimalt[optimalt$model=="model2",][which.min(optimalt[optimalt$model=="model2",4]),]

optimalk[optimalk$model=="model1"&optimalk$method=="LSE",][which.min(optimalk[optimalk$model=="model1"&optimalk$method=="LSE",6]),]
optimalk[optimalk$model=="model2"&optimalk$method=="LSE",][which.min(optimalk[optimalk$model=="model2"&optimalk$method=="LSE",6]),]
optimalk[optimalk$model=="model3"&optimalk$method=="LSE",][which.min(optimalk[optimalk$model=="model3"&optimalk$method=="LSE",6]),]
optimalk[optimalk$model=="model4"&optimalk$method=="LSE",][which.min(optimalk[optimalk$model=="model4"&optimalk$method=="LSE",6]),]
optimalk[optimalk$model=="model5"&optimalk$method=="LSE",][which.min(optimalk[optimalk$model=="model5"&optimalk$method=="LSE",6]),]

kvec1 = rbind(
optimalk[optimalk$model=="model1"&optimalk$method=="Borda",][which.min(optimalk[optimalk$model=="model1"&optimalk$method=="Borda",6]),3:5],
optimalk[optimalk$model=="model2"&optimalk$method=="Borda",][which.min(optimalk[optimalk$model=="model2"&optimalk$method=="Borda",6]),3:5],
optimalk[optimalk$model=="model3"&optimalk$method=="Borda",][which.min(optimalk[optimalk$model=="model3"&optimalk$method=="Borda",6]),3:5],
optimalk[optimalk$model=="model4"&optimalk$method=="Borda",][which.min(optimalk[optimalk$model=="model4"&optimalk$method=="Borda",6]),3:5],
optimalk[optimalk$model=="model5"&optimalk$method=="Borda",][which.min(optimalk[optimalk$model=="model5"&optimalk$method=="Borda",6]),3:5])




## Extra simulation to find optimal MSE for LSE according to (k1,k2,k3) 
optimalk2 = as.data.frame(matrix(nrow = 1100, ncol = 5))
names(optimalk2) = c("model","k1","k2","k3","MSE")
optimalk2$model = rep(c("model1","model2","model3","model4","model5"),each =220)


MSE3 = NULL
kindex = NULL
for(m in 1:5){
  for(s in 1:220){
    kvec = kind[s,]
    kindex = rbind(kindex,kvec)
    s1 = simulation_asym(50,40,30,mode = m,signal_level = 1)
    MSE3 = c(MSE3,mean((LSE_asym(s1$observe,kvec,mode = 3)-s1$signal)^2))
    print(paste("s_",s,"_ended",sep = ""))
  }
}

optimalk2[,2:4] = kindex
optimalk2$MSE = MSE3

## Getting Table S2 for LSE
kvec2 = rbind(optimalk2[optimalk$model=="model1",][which.min(optimalk[optimalk$model=="model1",5]),2:4],
optimalk2[optimalk$model=="model2",][which.min(optimalk[optimalk$model=="model2",5]),2:4],
optimalk2[optimalk$model=="model3",][which.min(optimalk[optimalk$model=="model3",5]),2:4],
optimalk2[optimalk$model=="model4",][which.min(optimalk[optimalk$model=="model4",5]),2:4],
optimalk2[optimalk$model=="model5",][which.min(optimalk[optimalk$model=="model5",5]),2:4])



asym = as.data.frame(matrix(nrow = 300, ncol = 4))
names(asym) = c("sim","model","method","mse")
asym$model = rep(1:5,each = 60)
asym$sim = rep(rep(1:20,each = 3),5)
asym$method = rep(c("Borda","LSE","Spectral"),100)


# Final TAble S2 for Borda count LSE and Spectral
threind = rbind(c(1,24),c(3,48),c(1,48),c(1,28),c(1,22))
kvec1 = rbind(c(2,1,2),c(1,2,2),c(1,3,3),c(2,1,2),c(1,4,4))
kvec2 = rbind(c(6,2,3),c(8,5,8),c(6,9,6),c(9,5,6),c(7,9,3))



# Obtain Table 3, final result is available in "Data_Table3/"
mse = NULL
for(m in 1:5){
  for(s in 1:20){
    set.seed(s)
    s1 = simulation_asym(50,40,30,mode = m,signal_level = 1)
    mse = c(mse,mean((Borda2_asym(s1$observe,2,kvec1[m,])-s1$signal)^2))
    mse = c(mse,mean((LSE_asym(s1$observe,kvec2[m,],mode = 3)-s1$signal)^2))
    mse= c(mse,mean((Spectral(s1$observe,threind[m,1],setdiff(1:3,threind[m,1]),threshold = threind[m,2])-s1$signal)^2))
    print(paste("m_",m,"s_",s,"_is done",sep = ""))
  }
}

asym$mse = mse


summaryasym = summarySE(asym, measurevar="mse", groupvars=c("model","method"))
save(summaryasym, file = "summaryasym.RData")
load(file = "summaryasym.RData")
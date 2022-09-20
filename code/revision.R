Borda3 = function(A,l,k){
  d = dim(A)[1]
  #sorting
  z = HSC(A,k,k,k,sym= T)$Cs
  o1 = 1:d
  s = 1
  for (i in 1:k) {
    o1[z==i] = sample(s:(s+length(which(z==i))-1),length(which(z==i)))
    s = s + length(which(z==i))
  }
  
  As = A[o1,o1,o1]
  
  #polynomial block approximation
  est = polytensor(As,l,k)
  
  #sorting back
  invo1 = invPerm(o1)
  Theta = est[invo1,invo1,invo1]
  
  return(Theta)
}


source('functions_blk.R')


sim = 1:20
dim  = (1:10)*10

set.seed(s)


result = as.data.frame(matrix(nrow =1800 ,ncol = 7))
names(result) =  c("sim","dim","type","model","method","MSE","gpsize")
result$type = "conti"
result$method = rep(c("cluster+poly1","cluster+poly2","cluster+poly3"),600)
  



MSE = NULL
simv = NULL
dimv = NULL

result$model = rep(c("model1","model3","model5"),each = 600)
result$gpsize = 2
i = 0
for(m in c(1,3,5)){
  for(s in sim){
    for(d in dim){
      simv = c(simv,s)
      dimv = c(dimv,d)
      s1 = simulation(d,mode = m,signal_level = 1)
      
      MSE = c(MSE,mean((Borda3(s1$observe,1,k)-s1$signal)^2))
      temp = mean((Borda3(s1$observe,2,k)-s1$signal)^2)
      MSE = c(MSE,temp)
    
      MSE = c(MSE,mean((Borda3(s1$observe,3,k)-s1$signal)^2))
      i = i+1
      print(temp)
      print(i)
    }
  }
}
  
    
length(rep(simv,each = 3))
length(rep(dimv,each = 3))
length(MSE)

result$sim[1:length(rep(simv,each = 3))] = rep(simv,each = 3)
result$dim[1:length(rep(dimv,each = 3))] = rep(dimv,each = 3)
result$MSE[1:length(MSE)] = MSE

nresult = result[1:length(MSE),]
nresult$model = as.factor(nresult$model)
nresult$method = as.factor(nresult$method)
nresult$type = as.factor(nresult$type)

save(nresult, file = "revision.RData")
load('revision.RData')

summary2 = summarySE(nresult, measurevar='MSE', groupvars=c("dim","type","model","method","gpsize"))
optimal = rbind(optimal,summary2)

best2 = optimal[optimal$method %in% c('LSE','Spectral','Borda2',"cluster+poly1","cluster+poly2","cluster+poly3"),]

library(ggplot2)
p11 = ggplot(best2[best2$model =="model1"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20));p11


p13 = ggplot(best2[best2$model =="model3"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20));p13


p15 = ggplot(best2[best2$model =="model5"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20));p15



######### Robustness against non monotonicity
source('functions_asym.R')
sim = 1:20
dim  = (2:10)*10

set.seed(s)


# Hyper-parameter tuning
result = as.data.frame(matrix(nrow = 84*9*2,ncol = 5))
names(result) = c('dim','method','k1','k2','mse')


kind = NULL
for(i in 1:10){
  for(j in 1:10){
    if(i*j<j | j*j<i){
      print("hi")
    }else{
      kind = rbind(kind,c(i,j,j))
    }
  }
}

result$method = rep(c('Borda','LSE'),84*9)
result$dim = rep(dim,each = 84*2)
result[,3] = rep(kind[,1],each = 2)
result[,4] = rep(kind[,2],each = 2)

MSE1 = NULL
for(d in dim){
  for(k in 1:84){
    kvec = kind[k,]
    s1 = simulation_asym(d,d,d,mode = 6,signal_level = 1)
    MSE1 = c(MSE1,mean((Borda2_asym(s1$observe,2,kvec)$Theta-s1$signal)^2))
    MSE1 = c(MSE1,mean((LSE_asym(s1$observe,kvec,mode = 3)-s1$signal)^2))
    print(paste("d_",d,'_k',kvec[1],"_ended",sep = ""))
  }
}

#####################################
result2 = as.data.frame(matrix(nrow = 81*9,ncol = 5))
names(result2) = c('dim','method','k1','k2','mse')
kind = NULL
for(i in 11:19){
  for(j in 11:19){
    if(i*j<j | j*j<i){
      print("hi")
    }else{
      kind = rbind(kind,c(i,j,j))
    }
  }
}

result2$method = 'LSE'
result2$dim = rep(dim,each = 81)
result2[,3] = kind[,1]
result2[,4] = kind[,2]

MSE2 = NULL
for(d in dim){
  s1 = simulation_asym(d,d,d,mode = 6,signal_level = 1)
  for(k in 1:81){
    kvec = kind[k,]
    MSE2 = c(MSE2,mean((LSE_asym(s1$observe,kvec,mode = 3)-s1$signal)^2))
    print(paste("d_",d,'_k',kvec[1],"_ended",sep = ""))
  }
}


result2$mse = MSE2

totresult = rbind(result,result2)
optimalk = NULL
for (d in dim) {
  for(m in c('Borda','LSE')){
    temp = totresult[totresult$dim==d&totresult$method == m,]
    optimalk = rbind(optimalk,temp[which.min(temp$mse),])
  }  
}

#################### optimal threhsold hyper-parameters ########################################
thre = 10+ 2*(1:22)
threind =  expand.grid(1:3,thre)
optimalt = as.data.frame(matrix(nrow = 66*9, ncol = 5))
names(optimalt) = c('dim', "method","mode","threshold","mse")
optimalt$dim = rep(dim,each = 66)
optimalt[,3:4] = rbind(threind,threind,threind,threind,threind,threind,threind,threind,threind)
optimalt$method = 'Spectral'

MSE3 = NULL
for(d in dim){
  s1 = simulation_asym(d,d,d,mode = 6,signal_level = 1)
  for(s in 1:66){
    MSE3 = c(MSE3,mean((Spectral(s1$observe,threind[s,1],setdiff(1:3,threind[s,1]),threshold = threind[s,2])-s1$signal)^2))
    print(paste("s_",s,"_ended",sep = ""))
  }
}

optimalt$mse = MSE3


thre = 54 + 2*(1:10)
threind =  expand.grid(1:3,thre)
optimalt2 = as.data.frame(matrix(nrow = 30*3, ncol = 5))
names(optimalt2) = c('dim', "method","mode","threshold","mse")
optimalt2$dim = rep(dim[7:9],each = 30)
optimalt2[,3:4] = rbind(threind,threind,threind)
optimalt2$method = 'Spectral'


#### EXTRA simultation #########
MSE4 = NULL
for(d in dim[7:9]){
  s1 = simulation_asym(d,d,d,mode = 6,signal_level = 1)
  for(s in 1:30){
    MSE4 = c(MSE4,mean((Spectral(s1$observe,threind[s,1],setdiff(1:3,threind[s,1]),threshold = threind[s,2])-s1$signal)^2))
    print(paste("s_",s,"_ended",sep = ""))
  }
}

optimalt2$mse = MSE4


optimalt = rbind(optimalt,optimalt2)

optimalthr = NULL
for (d in dim) {
  temp = optimalt[optimalt$dim==d,]
  optimalthr =  rbind(optimalthr,temp[which.min(temp$mse),])
}

save(optimalk,optimalthr, file = "hypertunning.RData")


# Final result
totresult = as.data.frame(matrix(nrow = 540,ncol = 4))
names(totresult) = c('dim','rep','method','mse')
totresult$dim = rep(dim,each = 60)
totresult$rep = rep(rep(sim,each = 3),9)
totresult$method = rep(c('Borda','LSE','Spectral'),180)

MSEt = NULL
for(d in dim){
  for(s in sim){
    s1 = simulation_asym(d,d,d,mode = 6,signal_level = 1)
    k1 = optimalk[optimalk$dim==d&optimalk$method=='Borda','k1']
    k2 = optimalk[optimalk$dim==d&optimalk$method=='Borda','k2']
    MSEt = c(MSEt,mean((Borda2_asym(s1$observe,2,c(k1,k2,k2))$Theta-s1$signal)^2))
    k1 = optimalk[optimalk$dim==d&optimalk$method=='LSE','k1']
    k2 = optimalk[optimalk$dim==d&optimalk$method=='LSE','k2']
    MSEt = c(MSEt,mean((LSE_asym(s1$observe,c(k1,k2,k2),mode = 3)-s1$signal)^2))
    md = optimalthr[optimalthr$dim==d,'mode']
    threshold = optimalthr[optimalthr$dim==d,'threshold']
    MSEt = c(MSEt,mean((Spectral(s1$observe,md,setdiff(1:3,md),threshold = threshold)-s1$signal)^2))
    print(paste('d_',d,'_s_',s,sep = ''))
  }
}

totresult$mse = MSEt
save(totresult,file = 'totresult.RData')

totresult$method = as.factor(totresult$method)
totsummary = summarySE(totresult, measurevar="mse", groupvars=c("dim","method"))

library(ggplot2)
pt = ggplot(totsummary, aes(x=dim, y=mse, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=mse-se, ymax=mse+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20));pt





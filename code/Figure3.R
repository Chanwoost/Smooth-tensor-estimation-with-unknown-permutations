# Obtaining simulation results for Figure 3 and Figure S3
# simulation repetetion (sim): 1,2,...,20
# group size (k): 1,2,...,15
# Model (m): 1,2,3,4,5
# type: binary and conti

# We use slurm to obtain the results with BATCH 1,2,.. 300
# Assigning BATCH in linux
args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])
} else {
  stop()
}


# Here we only give one example of obtaining data when BATCH = 1
# Whole result datasets are available in /Data_Figure3

BATCH = 1
source("functions_blk.R")
library(rTensor)

sim = 1:20
k = 1:15

ind = matrix(nrow = 300,ncol = 2)
r = 0
for(i in sim){
  for(j in k){
    r = r+1
    ind[r,] = c(i,j)
  }
}
s = ind[BATCH,1]
k = ind[BATCH,2]

set.seed(s)

result = as.data.frame(matrix(nrow =40 ,ncol = 6))
names(result) =  c("sim","groupsize","type","model","deg","MSE")

result$sim = s
result$groupsize = k
result$type = rep(c("conti","binary"),each = 20)
result$model = rep(rep(c("model1","model2","model3","model4","model5"),each = 4),2)
result$deg = rep(0:3,10)


MSE = NULL
for(m in 1:5){
  s1 = simulation(100,mode = m,signal_level = 1)
  MSE = c(MSE,mean((Borda2(s1$observe,0,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
}
for(m in 1:5){
  s1 = simulation_bin(100,mode = m)
  MSE = c(MSE,mean((Borda2(s1$observe,0,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,1,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
  MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
}


result$MSE = MSE

save(result,file = paste("noise1sim_",s,"_gsize_",k,".RData",sep =""))



############ Summarizing the simulation results from data ########################


tot= NULL
for(BATCH in 1:300){
  s = ind[BATCH,1]
  k = ind[BATCH,2]
  
  load(paste("Codes/Data_Figure3/noise1sim_",s,"_gsize_",k,".RData",sep =""))
  tot = rbind(tot,result)
}


tot$model = as.factor(tot$model)
tot$type = as.factor(tot$type)
tot$deg = as.factor(tot$deg)

totsummary = summarySE(tot, measurevar="MSE", groupvars=c("model","groupsize","deg","type"))

############ Plotting Figure 3 and Figure S3 ########################
# Continuous case
ggplot(totsummary[totsummary$model=="model1"&totsummary$type=='conti',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model2"&totsummary$type=='conti',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model3"&totsummary$type=='conti',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model4"&totsummary$type=='conti',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model5"&totsummary$type=='conti',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)


# Binary case
ggplot(totsummary[totsummary$model=="model1"&totsummary$type=='binary',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model2"&totsummary$type=='binary',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model3"&totsummary$type=='binary',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model4"&totsummary$type=='binary',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)
ggplot(totsummary[totsummary$model=="model5"&totsummary$type=='binary',], aes(x=groupsize, y=MSE,colour = deg)) + geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.1)












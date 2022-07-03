# Obtaining simulation results for Figure 4 and Figure S1
# simulation repetetion (sim): 1,2,...,20
# tensor dimension (dim): 10,20,...,100
# Model (m): 1,2,3,4,5
# type (t): binary, conti
# possible groupsize (k): 2,..,9
# We use slurm to obtain the results with BATCH 1,2,.. 400


# Assigning BATCH in linux
args <- (commandArgs(trailingOnly=TRUE))
cat(args[1])
if(length(args) == 1){
  BATCH <- as.numeric(args[1])
} else {
  stop()
}


# Here we only give one example of obtaining data when BATCH = 1
# Whole result datasets are available in /Data_Figure4
BATCH = 1
source("functions_blk.R")
library(rTensor)



sim = 1:20
dim  = (1:10)*10

ind = matrix(nrow = 400,ncol = 3)
r = 0
for(i in sim){
  for(j in dim){
    for(l in 1:2){
      r = r+1
      ind[r,] = c(i,j,l)
    }
  }
}
s = ind[BATCH,1]
d = ind[BATCH,2]
t = ind[BATCH,3]

set.seed(s)

if(t==1){
  result = as.data.frame(matrix(nrow =205 ,ncol = 7))
  names(result) =  c("sim","dim","type","model","method","MSE","gpsize")
  result$type = "conti"
  
}else{
  result = as.data.frame(matrix(nrow =245 ,ncol = 7))
  names(result) =  c("sim","dim","type","model","method","MSE","gpsize")
  result$type = "binary"
}

result$sim = s
result$dim = d

 

if(t==1){
  MSE = NULL
  result$model = rep(c("model1","model2","model3","model4","model5"),each = 41)
  result$gpsize = rep(c(rep(2:9,each = 5),0),5)
  result$method = rep(c(rep(c("Borda0","Borda1","Borda2","Borda3","LSE"),8),"Spectral"),5)
  for(m in 1:5){
    s1 = simulation(d,mode = m,signal_level = 1)
    for(k in 2:9){
      MSE = c(MSE,mean((Borda2(s1$observe,0,k)-s1$signal)^2))
     
      MSE = c(MSE,mean((Borda2(s1$observe,1,k)-s1$signal)^2))
      
      MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
      
      MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
      
      MSE = c(MSE,mean((LSE(s1$observe,k,mode = 3)-s1$signal)^2))
    }
    

    MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2))
    
  }
}else{
  MSE = NULL
  result$model = rep(c("model1","model2","model3","model4","model5"),each = 49)
  result$gpsize = rep(c(rep(2:9,each = 6),0),5)
  result$method = rep(c(rep(c("Borda0","Borda1","Borda2","Borda3","BAL","LSE"),8),"Spectral"),5)
  for(m in 1:5){
    s1 = simulation_bin(d,mode = m)
    for(k in 2:9){
      MSE = c(MSE,mean((Borda2(s1$observe,0,k)-s1$signal)^2))
      MSE = c(MSE,mean((Borda2(s1$observe,1,k)-s1$signal)^2))
      MSE = c(MSE,mean((Borda2(s1$observe,2,k)-s1$signal)^2))
      MSE = c(MSE,mean((Borda2(s1$observe,3,k)-s1$signal)^2))
      MSE = c(MSE,mean((LSE(s1$observe,k,mode = 2)-s1$signal)^2))
      MSE = c(MSE,mean((LSE(s1$observe,k,mode = 3)-s1$signal)^2))
    }
    
    MSE = c(MSE,mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2))
  }
}


result$MSE = MSE

save(result,file = paste("Csim_",s,"_dim_",d,"_type_",t,".RData",sep =""))




############ Summarizing the simulation results from data ########################

tot= NULL
for(BATCH in 1:400){
  s = ind[BATCH,1]
  d = ind[BATCH,2]
  t = ind[BATCH,3]
  
  load(paste("Codes/Data_Figure4/Csim_",s,"_dim_",d,"_type_",t,".RData",sep =""))
  tot = rbind(tot,result)
}

tot$model = as.factor(tot$model)
tot$methodl = as.factor(tot$method)
tot$type = as.factor(tot$type)
totsummary = summarySE(tot, measurevar="MSE", groupvars=c("dim","type","model","method","gpsize"))

library(dplyr)
optimal = as.data.frame(totsummary %>% 
  group_by(dim,type,model,method) %>% 
  slice(which.min(MSE)))

#### best is for Figure S1
best = optimal[optimal$method %in% c('Borda0','Borda1','Borda2','Borda3'),]

### best 2 is for Figure 4
best2 = optimal[optimal$method %in% c('LSE','Spectral','Borda2','BAL'),]


################################# Figure S1 ############################################################
## Dimension versus MSE according to different degrees
## Conti
load('summary.RData')
s11 = ggplot(best[best$model =="model1"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s11

s12 = ggplot(best[best$model =="model2"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s12

s13 = ggplot(best[best$model =="model3"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s13

s14 = ggplot(best[best$model =="model4"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s14

s15 = ggplot(best[best$model =="model5"&best$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s15



## Binary
s21 = ggplot(best[best$model =="model1"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s21

s22 = ggplot(best[best$model =="model2"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s22

s23 = ggplot(best[best$model =="model3"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s23

s24 = ggplot(best[best$model =="model4"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s24

s25 = ggplot(best[best$model =="model5"&best$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Degree 0", "Degree 1", "Degree 2","Degree 3"));s25


################################# Figure 4 ############################################################
## Dimension versus MSE according to 3 methods: Borda count, LSE, spectral
## Conti
p11 = ggplot(best2[best2$model =="model1"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count", "LSE","Spectral"));p11

p12 = ggplot(best2[best2$model =="model2"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Borda count", "LSE","Spectral"));p12

p13 = ggplot(best2[best2$model =="model3"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c( "Borda count", "LSE","Spectral"));p13

p14 = ggplot(best2[best2$model =="model4"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Borda count", "LSE","Spectral"));p14

p15 = ggplot(best2[best2$model =="model5"&best2$type =="conti",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("Borda count", "LSE","Spectral"));p15



## Dimension versus MSE according to 4 methods: Borda count, LSE, BAL, spectral
## Binary
p21 = ggplot(best2[best2$model =="model1"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p21

p22 = ggplot(best2[best2$model =="model2"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p22

p23 = ggplot(best2[best2$model =="model3"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p23

p24 = ggplot(best2[best2$model =="model4"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p24

p25 = ggplot(best2[best2$model =="model5"&best2$type =="binary",], aes(x=dim, y=MSE, colour=method)) + 
  geom_line() + geom_point()+ geom_errorbar(aes(ymin=MSE-se, ymax=MSE+se), width=.05)+
  labs(x="Dimension (d)", y="MSE",fill = "") + theme(text = element_text(size = 20))+
  scale_colour_discrete(name = "", labels = c("BAL","Borda count", "LSE","Spectral"));p25








load("Chicago_crime_data/Chicago.RData")
source("functions_asym.R")
source("functions_blk.R")

###################### denoising chicago crime tensor ###############################
result = Borda2_asym(ltns,2,c(6,4,10))
z1 = rep(1:6,rep(floor(24/6),6))
z2 = rep(1:4,rep(floor(77/4),4))
z3 = rep(1:10,rep(floor(31/10),10))
o1 = result$order[[1]]
o2 = result$order[[2]]
o3 = result$order[[3]]
result3 = Borda2_asym(ltns,0,c(7,18,10))

d = dim(ltns)
denoised = result$Theta
denoised[denoised<0]=0
denoised_sorted = denoised[1:d[1],result$order[[2]],result$order[[3]]]


############################ Figure 7 ##################################################
library(plot.matrix)
cM = array(dim = c(d[1],4,d[3]))
result$order[[2]]

for(i in 1:4){
  temp = matrix(0,nrow = d[1],ncol = d[3])
  for(j in which(z2==i)){
    temp = temp + denoised_sorted[,j,] 
  }
  cM[,i,] = temp/length(which(z2==i))
}

## Crime Area 1
par(mar = c(8, 6, 2, 4))
plot(cM[,1,],axis.col=NULL, axis.row=NULL, xlab='', ylab='Hour',main = "Crime Area 1", 
     breaks=0:8,col =  brewer.pal(8, "OrRd"),cex.lab = 1.3,cex.main = 1.5,cex.axis = 1.2)
text(x = 1:32,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c(crimetype_map[o3,]),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 45,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 0.6)
axis(side = 2,
     ## Rotate labels perpendicular to y-axis.
     las = 2,
     ## Adjust y-axis label positions.
     at = c(4,8,12,16,20,24),
     mgp = c(3, 0.75, 0),cex = 1.7)
mtext("Crime Type", side=1, line=5,cex = 1.3)


## Crime Area 2
plot(cM[,2,],axis.col=NULL, axis.row=NULL, xlab='', ylab='Hour',main = "Crime Area 2", 
     breaks=0:8,col =  brewer.pal(8, "OrRd"),cex.lab = 1.3,cex.main = 1.5,cex.axis = 1.2)
text(x = 1:32,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c(crimetype_map[o3,]),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 45,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 0.6)
axis(side = 2,
     ## Rotate labels perpendicular to y-axis.
     las = 2,
     ## Adjust y-axis label positions.
     at = c(4,8,12,16,20,24),
     mgp = c(3, 0.75, 0),cex = 1.7)
mtext("Crime Type", side=1, line=5,cex = 1.3)


## Crime Area 3
plot(cM[,3,],axis.col=NULL, axis.row=NULL, xlab='', ylab='Hour',main = "Crime Area 3", 
     breaks=0:8,col =  brewer.pal(8, "OrRd"),cex.lab = 1.3,cex.main = 1.5,cex.axis = 1.2)
text(x = 1:32,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c(crimetype_map[o3,]),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 45,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 0.6)
axis(side = 2,
     ## Rotate labels perpendicular to y-axis.
     las = 2,
     ## Adjust y-axis label positions.
     at = c(4,8,12,16,20,24),
     mgp = c(3, 0.75, 0),cex = 1.7)
mtext("Crime Type", side=1, line=5,cex = 1.3)


## Crime Area 4
plot(cM[,4,],axis.col=NULL, axis.row=NULL, xlab='', ylab='Hour',main = "Crime Area 4", 
     breaks=0:8,col =  brewer.pal(8, "OrRd"),cex.lab = 1.3,cex.main = 1.5,cex.axis = 1.2)
text(x = 1:32,
     ## Move labels to just below bottom of chart.
     y = par("usr")[3] - 0.45,
     ## Use names from the data list.
     labels = c(crimetype_map[o3,]),
     ## Change the clipping region.
     xpd = NA,
     ## Rotate the labels by 35 degrees.
     srt = 45,
     ## Adjust the labels to almost 100% right-justified.
     adj = 0.965,
     ## Increase label size.
     cex = 0.6)
axis(side = 2,
     ## Rotate labels perpendicular to y-axis.
     las = 2,
     ## Adjust y-axis label positions.
     at = c(4,8,12,16,20,24),
     mgp = c(3, 0.75, 0),cex = 1.7)
mtext("Crime Type", side=1, line=5,cex = 1.3)




######################### Figure 6 ############################################

sort(o2[which(z2 == 1)])
sort(o2[which(z2 == 2)])
sort(o2[which(z2 == 3)])
sort(o2[which(z2 == 4)])



########################## Table S3 ############################################
crimetype_map[o3[which(z3==1)],]
crimetype_map[o3[which(z3==2)],]
crimetype_map[o3[which(z3==3)],]
crimetype_map[o3[which(z3==4)],]
crimetype_map[o3[which(z3==5)],]
crimetype_map[o3[which(z3==6)],]
crimetype_map[o3[which(z3==7)],]
crimetype_map[o3[which(z3==8)],]
crimetype_map[o3[which(z3==9)],]
crimetype_map[o3[which(z3==10)],]




######################### Table 4 #################################
test = list()
tindex = sample(1:length(ltns),length(ltns))
q = length(tns)%/%5
for(i in 1:5){
  test[[i]] = tindex[(q*(i-1)+1):(q*i)]
}


test1 = test2 = train1 = train2  = NULL
for(i in 1:5){
  test_index = test[[i]]
  train_index = setdiff(1:length(ltns),test_index)
  train_tensor = ltns
  train_tensor[test_index] = NA
  result = Borda2_asym(train_tensor,2,c(6,4,10))
  train1 = c(train1,mean((result$Theta[train_index]-ltns[train_index])^2))
  test1 = c(test1,mean((result$Theta[test_index]-ltns[test_index])^2))

  
  result2 = Borda2_asym(train_tensor,0,c(7,11,10))
  
  train2 = c(train2,mean((result2$Theta[train_index]-ltns[train_index])^2))
  test2 = c(test2,mean((result2$Theta[test_index]-ltns[test_index])^2))
}

mean(test1)
sd(test1)
mean(test2)
sd(test2)






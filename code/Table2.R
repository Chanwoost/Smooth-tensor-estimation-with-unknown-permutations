# numerical CP and Tucker rank estimation in Table 2


source("functions_blk.R")
library("rTensor")
s1 = simulation(100,mode = 1)
tensor = as.tensor(s1$signal)


appx_rank=function(tensor,thresh=99.9,step=1){
  size=dim(tensor)
  size=sort(size)
  min=which(sqrt(cumsum(svd(unfold(tensor,1:2,3)@data)$d^2))>thresh*0.01)[1]
  res=test=NULL
  r=min
  while(r<=(size[1]*size[2])){
    rk=c(r,cp(tensor,r,max_iter = 100)$norm_percent)
    res=rbind(res,rk)
    if(rk[2]>thresh) break
    r=r+step;
  }
  return(res)
}

for(m in 1:5){
  s1 = simulation(100,mode = 5)
  tensor = as.tensor(s1$signal)
  appx_rank(tensor,step = 20)
  
}


# mode 1: 1,100.000
s1 = simulation(100,mode = 1)
tensor = as.tensor(s1$signal)
cp(tensor,1,max_iter = 100)$norm_percent
tucker(tensor,rank = c(1,1,1),max_iter = 100)$norm_percent
# mode 2: 3, 99.999
s1 = simulation(100,mode = 2)
tensor = as.tensor(s1$signal)
cp(tensor,3,max_iter = 100)$norm_percent
tucker(tensor,rank = c(2,2,2),max_iter = 100)$norm_percent


# mode 3: 9, 99.99295
s1 = simulation(100,mode = 3)
tensor = as.tensor(s1$signal)
cp(tensor,9,max_iter = 1000)$norm_percent
tucker(tensor,rank = c(4,4,4),max_iter = 100)$norm_percent

# mode 4: 100, 99.90037
s1 = simulation(100,mode = 4)
tensor = as.tensor(s1$signal)
cp(tensor,100,max_iter = 100)$norm_percent
tucker(tensor,rank = c(87,87,87),max_iter = 100)$norm_percent


# mode 5: 100,99.93
s1 = simulation(100,mode = 5)
tensor = as.tensor(s1$signal)
cp(tensor,100,max_iter = 100)$norm_percent
tucker(tensor,rank = c(56,56,56),max_iter = 10)$norm_percent






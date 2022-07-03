polytensor_asym=function(tensor, l, kvec){
  d=dim(tensor)
  est=array(dim=d)
  z1=rep(1:kvec[1],rep(floor(d[1]/kvec[1]),kvec[1]))
  z1=c(z1,rep(kvec[1],d[1]-length(z1)))
  z2=rep(1:kvec[2],rep(floor(d[2]/kvec[2]),kvec[2]))
  z2=c(z2,rep(kvec[2],d[2]-length(z2)))
  z3=rep(1:kvec[3],rep(floor(d[3]/kvec[3]),kvec[3]))
  z3=c(z3,rep(kvec[3],d[3]-length(z3)))
  
  for(i in 1:kvec[1]){
    for(j in 1:kvec[2]){
      for(q in 1:kvec[3]){
        subtensor=tensor[which(z1==i),which(z2==j),which(z3==q)]
        X1=c(slice.index(subtensor,1))
        X2=c(slice.index(subtensor,2))
        X3=c(slice.index(subtensor,3))
        if(l==0){
          fit=lm(c(subtensor)~1)
          est[which(z1==i),which(z2==j),which(z3==q)]=predict(fit,as.data.frame(cbind(X1,X2,X3)))
        }else{
          fit=lm(c(subtensor)~polym(X1,X2,X3,degree=l,raw=TRUE)) 
          est[which(z1==i),which(z2==j),which(z3==q)]=predict(fit,polym(X1,X2,X3,degree=l,raw=TRUE))
        }
      }
    }
  }
  return(est)
}

# Borda count estimation
library(Matrix)
Borda2_asym = function(A,l,kvec){
  d = dim(A)
  #sorting
  o1 = order(sapply(1:d[1], function(x) sum(A[x,,],na.rm = T)))
  o2 = order(sapply(1:d[2], function(x) sum(A[,x,],na.rm = T)))
  o3 = order(sapply(1:d[3], function(x) sum(A[,,x],na.rm = T)))
  As = A[o1,o2,o3]
  
  #polynomial block approximation
  est = polytensor_asym(As,l,kvec)
  
  #sorting back
  invo1 = invPerm(o1);invo2 = invPerm(o2);invo3 = invPerm(o3)
  
  Theta = est[invo1,invo2,invo3]
  o = list(o1,o2,o3)
  result = list(Theta = Theta,order = o)
  return(result)
}



# Spectral method with threshold.
# use the soft impute
library(softImpute)

Spectral = function(A,row_idx,col_idx,threshold = NULL,lambda = 1){
  d = dim(A)[1]
  A_unfolded = unfold(as.tensor(A),row_idx = row_idx,col_idx = col_idx)@data
  if(is.null(threshold)){
    threshold = sqrt(max(dim(A_unfolded))) 
  }
  if(any(is.na(A))==T){
    # made the hardimpute
    Decomp = softImpute(A_unfolded,rank.max = d-1,lambda = lambda)
    s = max(length(which(ifelse(Decomp$d>0,Decomp$d,0)>0)),1)
  }else{
    Decomp = svd(A_unfolded)  
    s = max(length(which(ifelse(Decomp$d>threshold,Decomp$d,0)>0)),1)
  }
  
  
  D = diag(Decomp$d,s)
  Theta = Decomp$u[,1:s,drop = F] %*% D %*% t(Decomp$v[,1:s,drop = F])
  Theta = fold(Theta,row_idx = row_idx,col_idx = col_idx,dim(A))@data
  return(Theta)
}









# Least square estimation
library(rTensor)
tensor_unfold = function(tensor,dim=1){
  if (dim == 1) unfold = aperm(tensor,c(3,2,1))
  if (dim == 2) unfold = aperm(tensor,c(1,3,2))
  if (dim == 3) unfold = tensor
  unfold = apply(unfold,3,c)
  if (is.vector(unfold)) return(as.matrix(unfold)) else return(t(unfold))
}

ReNumber = function (Cs,sort = T){
  newCs <- rep(NA, length(Cs))
  if(sort==T){
    uniq <- unique(Cs)
  }else{
    uniq <- sort(unique(Cs)) 
  }
  for (i in 1:length(uniq)) {
    newCs[Cs == uniq[i]] <- i
  }
  return(newCs)
}

UpdateMus_tensor = function (x, Cs, Ds, Es) {
  Cs = ReNumber(Cs,sort = F); Ds = ReNumber(Ds,sort = F); Es = ReNumber(Es,sort = F)
  
  mus = array(NA, c(length(unique(Cs)), length(unique(Ds)), length(unique(Es))))
  d = dim(mus)
  for (k in 1:d[1]){
    for (r in 1:d[2]){
      for (l in 1:d[3]){
        mus[k,r,l] = mean(x[Cs==k,Ds==r,Es==l],na.rm = T)
      }
    }
  }
  return(mus)
}



HSC = function(x,k,l,r,nstart = 40,sym = F){
  result = list()
  u1 = svd(tensor_unfold(x,1))$u[,1:k,drop = F]
  u2 = svd(tensor_unfold(x,2))$u[,1:l,drop = F]
  u3 = svd(tensor_unfold(x,3))$u[,1:r,drop = F]
  hu1 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u2),t(u3)),ms = c(2,3))@data,1))$u[,1:k]
  hu2 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u1),t(u3)),ms = c(1,3))@data,2))$u[,1:l]
  hu3 = svd(tensor_unfold(ttl(as.tensor(x),list(t(u1),t(u2)),ms = c(1,2))@data,3))$u[,1:r]
  Y1 = hu1%*%t(hu1)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu2),t(hu3)),ms = c(2,3))@data,1)
  Y2 = hu2%*%t(hu2)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu1),t(hu3)),ms = c(1,3))@data,2)
  Y3 = hu3%*%t(hu3)%*%tensor_unfold(ttl(as.tensor(x),list(t(hu1),t(hu2)),ms = c(1,2))@data,3)
  if(sym ==F){
    result$z1  = kmeans(Y1,k,nstart = nstart)$cluster
    result$z2  = kmeans(Y2,l,nstart = nstart)$cluster
    result$z3  = kmeans(Y3,r,nstart = nstart)$cluster  
  }else{
    result$z1 =result$z2 =result$z3 = kmeans(Y1,k,nstart = nstart)$cluster
  }
  return(result)
}




Bal_correct_asym = function(A,kvec){
  d = dim(A)
  z1 = kmeans(tensor_unfold(A,1),kvec[1])$cluster
  z2 = kmeans(tensor_unfold(A,2),kvec[2])$cluster
  z3 = kmeans(tensor_unfold(A,3),kvec[3])$cluster
  E1 = matrix(nrow = d[1],ncol = kvec[1])
  E2 = matrix(nrow = d[2],ncol = kvec[2])
  E3 = matrix(nrow = d[3],ncol = kvec[3])
  prob=array(0,dim=kvec)
  for(a in 1:kvec[1]){
    for(b in 1:kvec[2]){
      for(c in 1:kvec[3]){
        prob[a,b,c]=mean(A[which(z1==a),which(z2==b),which(z3==c)])
      }
    }
  }
  
  for(i in 1:d[1]){
    for(a in 1:kvec[1]){
      E1[i,a]=sum(dbinom(A[i,,],1, prob[a,z2,z3],log=TRUE))
    }
  }
  z1=apply(E1,1,which.max)
  
  for(i in 1:d[2]){
    for(a in 1:kvec[2]){
      E2[i,a]=sum(dbinom(A[,i,],1, prob[z1,a,z3],log=TRUE))
    }
  }
  z2 =apply(E2,1,which.max)
  
  for(i in 1:d[3]){
    for(a in 1:kvec[3]){
      E3[i,a]=sum(dbinom(A[,,i],1, prob[z1,z2,a],log=TRUE))
    }
  }
  z3=apply(E3,1,which.max)
  
  return(z)
}

LSE_asym = function(A,kvec,mode = 2){
  
  if(mode==1){
    #Spectral membership estimation
    z1 = kmeans(tensor_unfold(A,1),kvec[1],nstart = 100)$cluster
    z2 = kmeans(tensor_unfold(A,2),kvec[2],nstart = 100)$cluster
    z3 = kmeans(tensor_unfold(A,3),kvec[3],nstart = 100)$cluster
  }else if (mode==2){
    # Balasubramanian estimation
    result = Bal_correct_asym(A,kvec)
    z1 =  result$z1; z2 = result$z2; z3 = result$z3
  }else if (mode==3){
    # HSC membership estimation
    result = HSC(A,kvec[1],kvec[2],kvec[3])
    z1 =  result$z1; z2 = result$z2; z3 = result$z3
  }
  z1 = ReNumber(z1); z2 = ReNumber(z2); z3 = ReNumber(z3)
  
  mu.array = UpdateMus_tensor(A,z1,z2,z3)
  # this is for the technical error
  
  Theta = mu.array[z1,z2,z3, drop=FALSE]
  return(Theta)
}







# asymmetirc simulation
af1 = function(a){
  return(a[1]*a[2]+a[3])
}
af2 = function(a){
  return(a[1]^2+a[2]+a[2]*a[3]^2)
}
af3 = function(a){
  return(a[1]/(1+exp(-3*sum(a^2))))
}
af4 = function(a){
  return(log(1+max(a)+a[1]^2+a[2]*a[3]))
}
af5=function(a){
  return(exp(-a[1]-sqrt(a[2])-a[3]^2))
}




simulation_asym = function(d1,d2,d3, mode = 1,sigma = 0.5,signal_level=5){
  tensor=array(dim=c(d1,d2,d3))
  X1=c(slice.index(tensor,1))/d1
  X2=c(slice.index(tensor,2))/d2
  X3=c(slice.index(tensor,3))/d3
  if(mode==1){
    signal = array(apply(cbind(X1,X2,X3),1,af1),dim=c(d1,d2,d3))
  }else if(mode==2){
    signal = array(apply(cbind(X1,X2,X3),1,af2),dim=c(d1,d2,d3))
  }else if(mode==3){
    signal = array(apply(cbind(X1,X2,X3),1,af3),dim=c(d1,d2,d3))
  }else if(mode==4){
    signal = array(apply(cbind(X1,X2,X3),1,af4),dim=c(d1,d2,d3))
  }else if(mode==5){
    signal = array(apply(cbind(X1,X2,X3),1,af5),dim=c(d1,d2,d3))
  }
  signal=signal_level*signal/sqrt(mean(signal^2)) 
  observe = signal+array(rnorm(d1*d2*d3,0,sigma),dim = c(d1,d2,d3))
  
  return(list(signal= signal,observe=observe))
}









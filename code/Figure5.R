source("functions_blk.R")
source("tensor_visualization.R")
############### Tensor visualization #########################################
## This code is for Figure 5 and Figure S2
##################### Figure 5 ###################################
# Model 1
s1 = simulation(80,mode = 1,signal_level = 1,symnoise = T)
p = sample(1:80,80)
plot_tensor(s1$signal)
rgl.postscript('vS1.pdf', fmt = 'pdf')
plot_tensor(s1$observe[p,p,p])
rgl.postscript('vO1.pdf', fmt = 'pdf')



plot_tensor(Borda(s1$observe,2,3));mean((Borda(s1$observe,2,3)-s1$signal)^2)
rgl.postscript('vBorda1.pdf', fmt = 'pdf')

plot_tensor(Spectral(s1$observe,1,c(2,3)));mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2) 
rgl.postscript('vSp1.pdf', fmt = 'pdf')

plot_tensor(LSE(s1$observe,9,mode  =3));mean((LSE(s1$observe,9,mode  =3)-s1$signal)^2) 
rgl.postscript('vLSE1.pdf', fmt = 'pdf')


# Model 3
s1 = simulation(80,mode = 3,signal_level = 2)
plot_tensor(s1$signal)
rgl.postscript('vS3.pdf', fmt = 'pdf')
plot_tensor(s1$observe[p,p,p])
rgl.postscript('vO3.pdf', fmt = 'pdf')

plot_tensor(Borda(s1$observe,2,2));mean((Borda(s1$observe,2,2)-s1$signal)^2)
rgl.postscript('vBorda3.pdf', fmt = 'pdf')


plot_tensor(Spectral(s1$observe,1,c(2,3)));mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2) 
rgl.postscript('vSp3.pdf', fmt = 'pdf')

plot_tensor(LSE(s1$observe,8,mode  =3));mean((LSE(s1$observe,8,mode  =3)-s1$signal)^2) 
rgl.postscript('vLSE3.pdf', fmt = 'pdf')



# Model 5
s1 = simulation(80,mode = 5,signal_level = 1)
plot_tensor(s1$signal)
rgl.postscript('vS5.pdf', fmt = 'pdf')
plot_tensor(s1$observe[p,p,p])
rgl.postscript('vO5.pdf', fmt = 'pdf')


plot_tensor(Borda(s1$observe,2,4));
mean((Borda(s1$observe,2,4)-s1$signal)^2) 
rgl.postscript('vBorda5.pdf', fmt = 'pdf')


plot_tensor(Spectral(s1$observe,1,c(2,3)));mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2) 
rgl.postscript('vSp5.pdf', fmt = 'pdf')

plot_tensor(LSE(s1$observe,14,mode  =3));mean((LSE(s1$observe,14,mode  =3)-s1$signal)^2) 




######################## Figure S2 #############################################

# Model 2
s1 = simulation(80,mode = 2,signal_level = 1,symnoise = T)
p = sample(1:80,80)
plot_tensor(s1$signal)
rgl.postscript('vS2.pdf', fmt = 'pdf')
plot_tensor(s1$observe[p,p,p])
rgl.postscript('vO2.pdf', fmt = 'pdf')


plot_tensor(Borda(s1$observe,2,2));mean((Borda(s1$observe,2,2)-s1$signal)^2)
rgl.postscript('vBorda2.pdf', fmt = 'pdf')

plot_tensor(Spectral(s1$observe,1,c(2,3)));mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2) 
rgl.postscript('vSp2.pdf', fmt = 'pdf')

plot_tensor(LSE(s1$observe,9,mode  =3));mean((LSE(s1$observe,9,mode  =3)-s1$signal)^2)
rgl.postscript('vLSE2.pdf', fmt = 'pdf')


# Model 4
s1 = simulation(80,mode = 4,signal_level = 1,symnoise = T)
p = sample(1:80,80)
plot_tensor(s1$signal)
rgl.postscript('vS4.pdf', fmt = 'pdf')
plot_tensor(s1$observe[p,p,p])
rgl.postscript('vO4.pdf', fmt = 'pdf')


plot_tensor(Borda(s1$observe,2,3));mean((Borda(s1$observe,2,3)-s1$signal)^2)
rgl.postscript('vBorda4.pdf', fmt = 'pdf')

plot_tensor(Spectral(s1$observe,1,c(2,3)));mean((Spectral(s1$observe,1,c(2,3))-s1$signal)^2) 
rgl.postscript('vSp4.pdf', fmt = 'pdf')

plot_tensor(LSE(s1$observe,8,mode  =3));mean((LSE(s1$observe,8,mode  =3)-s1$signal)^2)
rgl.postscript('vLSE4.pdf', fmt = 'pdf')



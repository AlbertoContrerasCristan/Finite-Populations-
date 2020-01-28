###############################

# Vemos si aparecen cambios


function runif(n::Int64,aa::Float64,bb::Float64)
# random generation from a Uniform(aa,bb)  distribution
# NOTE:  
unif1 = Uniform(aa,bb)
return rand(unif1,n)
end
##############################
function rexp(n1::Int64,lambda::Float64)
# random generation from a Exponential(lambda)  distribution
# NOTE:  Julia generates from Exponential()  where lambda = 
exp1 = Exponential(lambda)
return rand(exp1,n1)
end
########################
function rbeta(n::Int64,shape1::Float64,shape2::Float64)
# random generation from a Beta(shape1,shape2)  distribution 
beta1 = Beta(shape1,shape2)
return rand(beta1,n)
end
########################
function rgamma(n::Int64,shape::Float64,theta::Float64)
# random generation from a Gamma(shape,scale)  distribution.
# NOTE:  Julia generates from Gamma(shape, theta)  where theta = 1/scale
gam = Gamma(shape,1.0/theta)
return rand(gam,n)
end
####################################
function rnorm(n::Int64,mu::Float64,sig::Float64)
# random generation from a Normal(mu,sig^2) distribution
# Note: As in R,  sig is the. s.d. Not the variance parameter
norm1 = Normal(mu, sig)
return rand(norm1,n)
end
###################
function InvXi(m::Int64)
#  computes de sequence 1/xi_k = 2^{k}
#  which  is used to reduce correlation between u_i's and w_k's

kappa = 2.0

return [kappa^k for k in 1:m]
end
####################
function Nis(uv::Array{Float64,1},nu::Int64)
# determines the set of integers N_i in slice sampler
# algorithm (see Kalli, et al , 2011)
# NOTE :  if uv[i] > 0.5 returns a 0

#n <- length(uv)

Ni = [i for i in 1:nu]

for i = 1:nu

u = uv[i]
j = 0
Xhi = 0.5
while Xhi > u
j += 1		
Xhi *= 0.5
end	

Ni[i] = j	
end 

return Ni	
end
###################
function Weight0(Nn::Int64,PPr1::Float64)
# produces a sample of the weights
# (w_1, ... ,w_Nn)  in Sethuraman's representation for D.P.

##beta1 = Beta(1.0,PPr1)
##Nu = rand(beta1,Nn)

Nu = rbeta(Nn,1.0,PPr1)
Vv = 1.0.-Nu

for  i = 2:Nn
Nu[i] *= prod(Vv[1:i-1])		
end

return Nu
end
##################
function GaMix0(MM::Int64,fn::Float64,eta::Float64)
# computes mixture of gamma densities
# as in eq (13) of Escobar & West (1995)
#  THIS VERSION TO GENERATE INITIAL VALUES FOR MCMC

aux = (fn - 1.0)/(fn*(-log(eta)))
# piw = aux/(1.0+aux)
aux /= 1.0+aux
alfa = [Float64(i) for i in 1:MM]

for ii = 1:MM 
u0 = runif(1,0.0,1.0)[1]
if u0 < aux
alpha = rgamma(1,fn,-log(eta)) 
else
alpha = rgamma(1,fn-1.0,-log(eta))
end

alfa[ii] = alpha[1]
end

return mean(alfa)
end
##################
function GaMix(aa::Float64,bb::Float64,fk::Float64,fn::Float64,etta::Float64)
# computes mixture of gamma densities
# as in eq (13) of Escobar & West (1995)

aux = (aa + fk - 1.0)/(fn*(bb-log(etta)))
# piw = aux/(1.0+aux)
aux /= 1.0+aux

u2 = runif(1,0.0,1.0)[1]
if u2 < aux
alpha = rgamma(1,aa+fk,bb-log(etta)) 
else
alpha = rgamma(1,aa+fk-1.0,bb-log(etta))
end

return alpha[1]
end
################
function Betta(alfa::Float64,fn::Float64)
# simulates from beta distribution for 
# latent var eta  as in eq. (14)  Escobar & West (1995)	
	
#eta = rbeta(1,alfa+1.0,fn)			
#return eta[1]

return rbeta(1,alfa+1.0,fn)[1]
end
################
function NIG0(m::Float64,tau::Float64,S::Float64,df::Float64,Nn::Int64)
# Generates random vectors (mui,Vi) from the bi-variate prior
# Normal-InverseGamma((mui,Vi) |  m,tau,S,df)  distribution
# as in Escobar & West (1995)
beta = S/2.0
alpha = df/2.0
V = rgamma(Nn,alpha,beta)
V=1.0./V
# conditionally on Vi generate  mu_i from N(m, tau*Vi)
sig = sqrt.(tau.*V)
mu =  collect(rnorm(1,m,sig[k])[1] for k in 1:Nn)

# return a Nn*2 matrix with each row given by (mui,Vi)
return hcat(mu,V)	
end
########################
function NIG(x::Float64,m1::Float64,tau1::Float64,S1::Float64,df1::Float64)
# Generates random vectors (mui,Vi) from the bi-variate
# posterior Normal-InverseGamma( (mui,Vi) | (df+1)/2, Sx/2)    
# distribution (see Escobar & West 1995,sect 3)
Sx = (((x-m1)^2)/(1.0+tau1))  + S1
Beta = 0.5*Sx
alpha = 0.5*(df1+1.0)
# generate from Gamma(alpha,Beta) then
# transform to Inverse-Gamma(alpha,Beta)
Vi = rgamma(1,alpha,Beta)[1]
Vi=1.0/Vi
# conditionally on Vi generate from N(c,(Vi[tau/(tau+1)])^{2})

mu = (tau1*x+m1)/(1.0+tau1)
sig = sqrt((tau1*Vi)/(1.0+tau1))
mui = rnorm(1,mu,sig)
#return hcat(mui,Vi)
return  [mui Vi]
end
#########################
function NIG00(x::Array{Float64,1},m::Float64,tau::Float64,S::Float64,df::Float64,Nn1::Int64)
# Generates random vectors (mui,Vi) from the bi-variate
# POSTERIOR Normal-InverseGamma((mui,Vi) | x, m,tau,S,df)  distribution
# x = data. Nn1 = number of samples to generate
Theta = zeros(Nn1,2)
for i =1:Nn1
Theta[i,:]  = NIG(maximum(x),m,tau,S,df) 
end
return Theta

end
##########################
function dsamplerC(weights::Array{Float64,1},vTheta::Array{Float64,2},Nvect::Array{Int64,1},data1::Array{Float64,2})
# Given vectors (w_1,...w_N) ; (u_1, ... , u_n) 
# vTheta' = ((mu_1,Ka_1), ... ,(mu_N,Ka_N))  produces a 
# sample from the discrete indicators (d1, ... , dn)
nd=size(data1)[1]
dd =  [0 for i1=1:nd]

for i = 1:nd 
ww = weights[1:Nvect[i]].*InvXi(Nvect[i])

dens = [pdf(Normal(vTheta[k,1],vTheta[k,2]),data1[i])  for k=1:Nvect[i]]
Pei = ww.*dens

# re-normalize by constant
Pei /= sum(Pei)

# sample di from this discrete dist.
Pei = cumsum(Pei)
u3 = runif(1,0.0,1.0)[1]
c1 = 1
while u3>Pei[c1]
   	c1 += 1
end

dd[i] = c1 

# end for i
end

return dd
end
########################
function unifvec(dd::Array{Int64,1},Xiv::Array{Float64,1},nn::Int64)
# Given Xi sequence of weights and 
# d1, ... , dn  samples values u1, ... , un	

##univec = [0.0 for i2=1:nn]	
##univec =  [runif(1,0.0,Xiv[dd[i]])  i=1:nn]		

	
return [runif(1,0.0,Xiv[dd[i]])[1]  for i=1:nn]	
end
#########################
function buildDj(dv::Array{Int64,1},nn::Int64,j::Int64)
# given d1, ... , dn
# builds the sets Dj = {i: d_i = j}; 	

## index <- {1:nn}
ind1 = collect(i for i=1:nn if dv[i]==j)

##ind1 <- index[  dv == j  ]	
##card1 = length(ind1)
##ind2 <- index[  dv > j  ]

ind2 = collect(i for i=1:nn if dv[i]>j)

##card2 = length(ind2)
##return(list(Dj=ind1, nj=card1, mj=card2))			

return ind1,length(ind1),length(ind2)
end
#########################
function FuConThet(yy::Array{Float64,2},em::Float64,taau::Float64,sdd::Float64,S1::Float64,nThet::Array{Float64,2},enj::Int64,Dj::Array{Int64,1})
# updates the value of Theta=(miu_j,V_j)
# according to Bayes Rule	

alpha = (float(enj)+sdd+1.0)/2.0
sum1 = (em - nThet[1,1])*(em - nThet[1,1])/taau
sum3 = sum((yy[Dj,1].-nThet[1,1]).*(yy[Dj,1].-nThet[1,1]))
Beta = (sum1 + S1 + sum3)/2.0

V = rgamma(1,alpha,Beta)[1]
V = 1.0/V
sumy = sum(yy[Dj,1])
k1 = sumy + (em/taau)
k2 = float(enj) + (1.0/taau)
miu = k1/k2
sig = sqrt(V/k2)
mu = rnorm(1,miu,sig)[1]
return [mu V]
end
######################
function ThetaWeight(dv1::Array{Int64,1},dat1::Array{Float64,2},ndat::Int64,m1::Float64,ttau::Float64,sd::Float64,S1::Float64,newN::Int64,Thetta::Array{Float64,2},Nold::Int64,M::Float64)
# Given d1, ... ,  dn (dv1);  a dataset (dat1)	
# hyperpars m (m1), tau (ttau), s (sd),S (S1);  produces posterior 
# samples  (mu_1,V_1), ... (mu_newN, V_newN), (w_1, ... , w_newN)	
Nu = [0.0 for i=1:newN]
weights = [0.0 for i=1:newN]
# this counter will register number of clusters
Kl  = 0
for j1=1:newN
listDj = buildDj(dv1,ndat,j1)
 D = listDj[1]
 njj = listDj[2]	
 mjj = listDj[3]
 
 # update Thetta_j1
if  j1 > Nold
nTheta = NIG0(m1,ttau,S1,sd,1)
  if njj == 0
 	Thetta = vcat(Thetta,nTheta)
  else
 	println("redundant?")
 	Kl += 1
    nTheta = FuConThet(dat1,m1,ttau,sd,S1,nTheta,njj,D)
  	Thetta = vcat(Thetta,nTheta)
  end
else
   if njj == 0
    nTheta = NIG0(m1,ttau,S1,sd,1)
 	Thetta[j1,:] = nTheta
   else
 	Kl += 1	
 	nTheta = FuConThet(dat1,m1,ttau,sd,S1,hcat(Thetta[j1,1],Thetta[j1,2]),njj,D)
 	Thetta[j1,:] = nTheta	
   end
end

# update v_j1
Nu[j1] = rbeta(1,1.0+njj,M+mjj)[1]
# update weights_j1
if j1 > 1
weights[j1] = Nu[j1]*prod(1.0.-Nu[1:j1-1])
end
#end for j1
end

if newN <= Nold
Thetta = Thetta[1:newN,:]
end

weights[1] = Nu[1]

return Thetta,weights,Kl
end
########################
function NormMixMCMC(aa::Float64,bb::Float64,sdf::Float64,sS::Float64,mm::Float64,tauu::Float64,nMCMC::Int64,nBurn::Int64,nfreq::Int64,PP::Float64,dd::Array{Int64,1},Theta::Array{Float64,2},data::Array{Float64,2},tN::Int64)
# Given initial value Theta'(Theta_1,...Theta_Nlast) = ((mu1,Ka1) ... (mu_Nlast,Ka_Nlast))
# and selectors  d1,...,dn, performs MCMC.  
# 

sizeMCMC = Int64(nMCMC/nfreq)
print("mcmc running and to store");print(" ",sizeMCMC," posterior samples...","\n")
n = size(data)[1]
Nsim = tN-n
Rsim = 9000-1000
print("value of R-n is");print(" ",Rsim,"\n")

Ytotal = zeros(sizeMCMC)
Nlast = size(Theta)[1]
Xivec = 1.0./InvXi(Nlast)
nc1 = 0
alphk = 0.0
Nclust = 0
Nt  = 0 

# start MCMC
for j = 1:(nMCMC+nBurn) 
 if j%10000 == 0 
   println("Numb. Iter.="," ",j,"\n")
 end
 # update u1, ... , un
      uvecj = unifvec(dd,Xivec,n)
# update N_1, ... , N_n and N
      Nvecj = Nis(uvecj,n)
      nuevaN = maximum(Nvecj)
# update Theta_1, ... , Theta_nuevaN and w_1, ... , w_nuevaN and number of clusters KK
     ParsPesos = ThetaWeight(dd,data,n,mm,tauu,sdf,sS,nuevaN,Theta,Nlast,PP)
     Theta = ParsPesos[1] 
     ww = ParsPesos[2]
     Kj = ParsPesos[3]	

# sample total mass parameter     
      etta = Betta(PP,Float64(n))
##      println(etta)
      PP =  GaMix(aa,bb,Float64(Kj),Float64(n),etta) 
# update vector of selectors
      dd = dsamplerC(ww,Theta,Nvecj,data)
# update new number   of components
      Nlast = nuevaN
      Xivec = 1.0./InvXi(Nlast)
      
## collect samples if condition applies      
if  j>nBurn && j%nfreq==0
     nc1 += 1   
     alphk = hcat(alphk,PP)
     Nclust = hcat(Nclust,Kj)
     Nt = hcat(Nt,Nlast)
#  renormalizing weights so to ensure propriety of mixture 
    ww /= sum(ww)
    ww1 = cumsum(ww)
#  this aux var  will store sum of simulated Y_i values
    Ysum  = 0
# proceed to simulate and ADD y1,...,yNsim for  j-th  mixture
#  CORRECT proceed to simulate and  ADD y1,...,yRsim for  j-th  mixture
for nY = 1:Rsim
    u = runif(1,0.0,1.0)[1]
    c2 = 1
while u>ww1[c2]
   	c2 += 1
end
Ysum = Ysum +  rnorm(1,Theta[c2,1],sqrt(Theta[c2,2]))[1]
#end nY
end
Ytotal[nc1] = Ysum

# end if j>nBurn && j%nfreq==0
end

# end for j	
end
#####     MCMC  ENDS HERE  #####
println("End of MCMC---> now storing output")
print("nc1 is: ")
println(nc1)
# ouput results here 

len = length(alphk)

alphout = alphk[2:len]
KKout = Nclust[2:len]
Nout = Nt[2:len]
#  output results for assessing convergence of MCMC
outfilealfa = open("alfa.txt","w")
outfileNclus = open("Nclus.txt","w")
outfileNj = open("Nj.txt","w")
writedlm(outfilealfa,alphout)
writedlm(outfileNclus,KKout)
writedlm(outfileNj,Nout)
close(outfilealfa)
close(outfileNclus)
close(outfileNj)

##return Ytotal
### CORRECT  Ytotal here #############
#return(Ytotal)
return (float(Nsim)/float(Rsim))*Ytotal
end
######################
##### Main Program ######

using Compat, Random, Distributions, DelimitedFiles

function main1()

seed1 = 170081
Random.seed!(seed1)

## read data
n = 1000
outfile = open("muestra.txt","r")
mat1 = Matrix{Float64}(undef,n,1)
mat1[:,1] = readdlm(outfile,dims=(n,1))
close(outfile)
trueN = 100000

Data = mat1[:,1]

ndat = size(mat1)[1]
# these are hyperparameters for Gamma prior 
# on PPr = total mass par.
a = 2.0
b = 4.0

sdf = 4.0
sS = 2.0
sm = mean(Data)
stau = 15.0  

# number of posterior samples
nS = 25000
# number of initial samples to be burned
burn = 30000
#  How frequent to keep the sample,  
#   e.g. freqS <- 2 =>  every other sample
freqS = 25   

# OBTAIN INITIAL VALUES FOR MCMC
# initial value for total mass parameter
etta = rbeta(1000,1.0,Float64(ndat))
eta0 = mean(etta)
PPr = GaMix0(1000,Float64(ndat),eta0)
#PPr = 11.87

uvec = runif(ndat,0.0,0.5)
Nvec = Nis(uvec,ndat)
NN = maximum(Nvec)

Theta0 = NIG00(Data,sm,stau,sS,sdf,NN)

weight = Weight0(NN,PPr)
dvec = dsamplerC(weight,Theta0,Nvec,mat1)

#listDj = buildDj(dvec,ndat,1)
#D = listDj[1]
#njj = listDj[2]	
#mjj = listDj[3]

########## start MCMC and get sample of totales #####################

#nTheta = NIG0(sm,stau,sS,sdf,1)
#Thet = FuConThet(mat1,sm,stau,sdf,sS,nTheta,njj,D)
#nTheta1 = NIG0(sm,stau,sS,sdf,1)
#Thetta = vcat(nTheta,nTheta1)
#Nold = size(Thetta)[1]
#println(Nold)
#return ThetaWeight(dvec,mat1,ndat,sm,stau,sm,sS,3,Thetta,Nold,PPr)

STotales = NormMixMCMC(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,mat1,trueN)

###### then analize totales #
T1 = trueN-ndat
print("N-(n) =");print(" ",T1,"\n")
#cat("N-(n) =",T1,"\n")
tt = sum(Data)
print("t =");print(" ",tt,"\n")
#cat("t =",tt,"\n")

nTotales = length(STotales)
print("nTotales =");print(" ",nTotales,"\n")

Totales = STotales.+tt
print("mean =");print(" ",mean(Totales),"\n")

print("seed =");print(" ",seed1,"\n")

qw = quantile(Totales,[0.025 0.05 0.25 0.5 0.75 0.95 0.975])
print(qw)
outfile1 = open("cuantiles-muestra.txt","w")
outfile2 = open("Totales-muestra.txt","w")
writedlm(outfile1,qw)
writedlm(outfile2,Totales)
close(outfile1)
close(outfile2)



end
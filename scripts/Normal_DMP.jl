
# Libraries
using Statistics, Compat, Random, Distributions, DelimitedFiles, DataFrames, Dates
# Load

cd("$(homedir())\\Documents\\Finite-Populations-")
Base.include(Main,"$(homedir())\\Documents\\Finite-Populations-\\scripts\\helpers.jl")
data = DataFrame(readdlm("$(homedir())\\Documents\\Finite-Populations-\\data\\muestra.txt"))
ndat = nrow(data)
data_mat = Matrix{Float64}(undef,ndat,1)

data_mat[:,1] = data[1]

# Sample and True Size


trueN = 100000
# Replication

seed1 = 1234
Random.seed!(seed1)

# Hyperparameters
## Gamma prior
### on PPr = total mass par.

a, b, sdf, sS, stau = 2.0, 4.0, 4.0, 2.0, 15.0
sm = mean(data[1])

# MCMC Parameters
## number of posterior samples & number of burned
## freqS is How frequent to keep the sample,
### #   e.g. freqS <- 2 =>  every other sample
nS, burn, freqS = 25000, 30000, 25

# OBTAIN INITIAL VALUES FOR MCMC
## initial value for total mass parameter
etta = rbeta(1000,1.0,Float64(ndat))
eta0 = mean(etta)
PPr = GaMix0(1000,Float64(ndat),eta0) #PPr = 11.87


uvec = runif(ndat,0.0,0.5)
Nvec = Nis(uvec,ndat)
NN = maximum(Nvec)

Theta0 = NIG00(data[1],sm,stau,sS,sdf,NN)

weight = Weight0(NN,PPr)
dvec = dsamplerC(weight,Theta0,Nvec,data_mat)

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

STotales = NormMixMCMC(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN)

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

#outfile2 = open("Totales-muestra.txt","w")
#writedlm(outfile2,Totales)

#close(outfile2)

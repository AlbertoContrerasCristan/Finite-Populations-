# Load

cd("$(homedir())/Documents/Finite-Popúlations-")

include("$(homedir())/Documents/Finite-Populations-/scripts/helpers.jl")

##### Main Program ######

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

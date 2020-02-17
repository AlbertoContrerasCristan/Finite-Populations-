# Libraries
using Statistics, Compat, Random, Distributions, DelimitedFiles, DataFrames, Dates
# Load

cd("$(homedir())\\Documents\\Finite-Populations-")
include("$(homedir())\\Documents\\Finite-Populations-\\scripts\\helpers.jl")
data = DataFrame(readdlm("$(homedir())\\Documents\\Finite-Populations-\\data\\muestra.txt"))
ndat = nrow(data)
data_mat = Matrix{Float64}(undef,ndat,1)

data_mat[:,1] = data[1]

# Sample and True Size


trueN = 100000
# Replication

seed1 = 170081
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

R_n1 = collect(50:50:5000) # Grid de R < N
R_n2 = collect(5000:100:10000)
R_n2 = collect(10000:250:30000)
R_n2 = collect(30000:500:60000)
R_n4 = collect(60000:1000:90000)

means_short = Array{Any}(undef,length(R_n1))
runtimes_short = Array{Any}(undef,length(R_n1))

for j in 1:length(R_n1)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n1[j])
    means_short[j] = mean(out[1])
    runtimes_short[j] = out[5]

end

graph_out1 = Matrix{Float64}(R_n1, means_short,runtimes_short)

graph1 = open("graph1.csv","w")
writedlm(graph1,graph_out1)
close(graph1)

means_short2 = Array{Any}(undef,length(R_n2))
runtimes_short2 = Array{Any}(undef,length(R_n2))

for j in 1:length(R_n2)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n2[j])
    means_short2[j] = mean(out[1])
    runtimes_short2[j] = out[5]

end

graph_out2 = Matrix{Float64}(R_n2, means_short2,runtimes_short2)

graph2 = open("graph2.csv","w")
writedlm(graph2,graph_out2)
close(graph2)

means_short3 = Array{Any}(undef,length(R_n3))
runtimes_short = Array{Any}(undef,length(R_n3))

for j in 1:length(R_n3)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n3[j])
    means_short3[j] = mean(out[1])
    runtimes_short3[j] = out[5]

end

graph_out3 = Matrix{Float64}(R_n3, means_short3,runtimes_short3)

graph3 = open("graph3.csv","w")
writedlm(graph3,graph_out3)
close(graph3)

means_short4 = Array{Any}(undef,length(R_n4))
runtimes_short4 = Array{Any}(undef,length(R_n4))

for j in 1:length(R_n4)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n4[j])
    means_short4[j] = mean(out[1])
    runtimes_short4[j] = out[5]

end

graph_out4 = Matrix{Float64}(R_n4, means_short4,runtimes_short4)

graph4 = open("graph4.csv","w")
writedlm(graph4,graph_out4)
close(graph4)

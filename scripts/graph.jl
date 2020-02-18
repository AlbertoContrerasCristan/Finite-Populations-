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

R_n1 = collect(50:50:5000) # Grid de R < N
R_n2 = collect(5100:100:10000)
R_n3 = collect(10250:250:30000)
R_n4 = collect(30500:500:60000)
R_n5 = collect(61000:1000:90000)

means_short = Array{Any}(undef,length(R_n1))
runtimes_short = Array{Any}(undef,length(R_n1))

for j in 1:length(R_n1)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n1[j])
    means_short[j] = mean(out[1])
    runtimes_short[j] = out[5]

end

r_n1 = open("r_n1.csv","w")
writedlm(r_n1,R_n1)
close(r_n1)

mean1 = open("mean1.csv","w")
writedlm(mean1,means_short/1000000)
close(mean1)

runtime1 = open("runtime1.csv","w")
writedlm(runtime1,Dates.value.(runtimes_short)/-60000)
close(runtime1)

means_short2 = Array{Any}(undef,length(R_n2))
runtimes_short2 = Array{Any}(undef,length(R_n2))

for j in 1:length(R_n2)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n2[j])
    means_short2[j] = mean(out[1])
    runtimes_short2[j] = out[5]

end

r_n2 = open("r_n2.csv","w")
writedlm(r_n2,R_n2)
close(r_n2)

mean2 = open("mean2.csv","w")
writedlm(mean2,means_short2/1000000)
close(mean2)

runtime2 = open("runtime2.csv","w")
writedlm(runtime2,Dates.value.(runtimes_short2)/-60000)
close(runtime2)


means_short3 = Array{Any}(undef,length(R_n3))
runtimes_short3 = Array{Any}(undef,length(R_n3))

for j in 1:length(R_n3)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n3[j])
    means_short3[j] = mean(out[1])
    runtimes_short3[j] = out[5]

end

r_n3 = open("r_n3.csv","w")
writedlm(r_n3,R_n3)
close(r_n3)

mean3 = open("mean3.csv","w")
writedlm(mean3,means_short3/1000000)
close(mean3)

runtime3 = open("runtime3.csv","w")
writedlm(runtime3,Dates.value.(runtimes_short3)/-60000)
close(runtime3)

means_short4 = Array{Any}(undef,length(R_n4))
runtimes_short4 = Array{Any}(undef,length(R_n4))

for j in 1:length(R_n4)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n4[j])
    means_short4[j] = mean(out[1])
    runtimes_short4[j] = out[5]

end

r_n4 = open("r_n4.csv","w")
writedlm(r_n4,R_n4)
close(r_n4)

mean4 = open("mean4.csv","w")
writedlm(mean4,means_short4/1000000)
close(mean4)

runtime4 = open("runtime4.csv","w")
writedlm(runtime4,Dates.value.(runtimes_short4)/-60000)
close(runtime4)

means_short5 = Array{Any}(undef,length(R_n5))
runtimes_short5 = Array{Any}(undef,length(R_n5))

for j in 1:length(R_n5)

    out = NormMixMCMC_plot(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,data_mat,trueN, R_n5[j])
    means_short5[j] = mean(out[1])
    runtimes_short5[j] = out[5]

end

r_n5 = open("r_n5.csv","w")
writedlm(r_n5,R_n5)
close(r_n5)

mean5 = open("mean5.csv","w")
writedlm(mean5,means_short5/1000000)
close(mean5)

runtime5 = open("runtime5.csv","w")
writedlm(runtime5,Dates.value.(runtimes_short5)/-60000)
close(runtime5)

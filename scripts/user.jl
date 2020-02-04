try =NormMixMCMC(a,b,sdf,sS,sm,stau,nS,burn,freqS,PPr,dvec,Theta0,mat1,trueN)

function NormMixMCMC_user(aa::Float64,
                    bb::Float64,
                    sdf::Float64,
                    sS::Float64,
                    tauu::Float64,
                    nMCMC::Int64,
                    nBurn::Int64,
                    nfreq::Int64,
                    shape1::Float64 = 1.0,
                    PP::Float64,
                    dd::Array{Int64,1},
                    Theta::Array{Float64,2},
                    data::Array{Float64,2},
                    tN::Int64)
    # Given initial value Theta'(Theta_1,...Theta_Nlast) = ((mu1,Ka1) ... (mu_Nlast,Ka_Nlast))
    # and selectors  d1,...,dn, performs MCMC.
    init = Dates.format(now(), "HH:MM")

    print("Getting Initial Values for MCMC")
    mm = float(mean(data)) # En vez de parámetro mm
    n = size(data)[1]
    # OBTAIN INITIAL VALUES FOR MCMC
    ## initial value for total mass parameter
    ### En vez de Parámetro PP
    etta = rbeta(1000,shape1,Float64(n))
    eta0 = mean(etta)
    PPr = GaMix0(1000,Float64(n),eta0)

    ### En vez de parametros dd y Theta
    uvec = runif(n,0.0,0.5)
    Nvec = Nis(uvec,n)
    NN = maximum(Nvec)

    Theta = NIG00(data[1],mm,stau,sS,sdf,NN)

    weight = Weight0(NN,PPr)
    dd = dsamplerC(weight,Theta0,Nvec,data_mat)

    print("Initial Values Calculated")
    print("Proceding to run MCMC")


    sizeMCMC = Int64(nMCMC/nfreq)
    print("mcmc running and to store")
    print(" ",sizeMCMC," posterior samples...","\n")
    Nsim = tN-n
    Rsim = 9000-1000

    print("value of R-n is")
    print(" ",Rsim,"\n")

    Ytotal = zeros(sizeMCMC)
    Nlast = size(Theta)[1]
    Xivec = 1.0./InvXi(Nlast)
    nc1 = 0
    alphk = 0.0
    Nclust = 0
    Nt  = 0
    for j = 1:(nMCMC+nBurn) # start MCMC
        if j%10000 == 0
            println("Numb. Iter.="," ",j,"\n")
        end
        uvecj = unifvec(dd,Xivec,n) # update u1, ... , un
        Nvecj = Nis(uvecj,n) # update N_1, ... , N_n and N
        nuevaN = maximum(Nvecj)
        # update Theta_1, ... , Theta_nuevaN
        # and w_1, ... , w_nuevaN and number of clusters KK
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
            # renormalizing weights so to ensure propriety of mixture
            ww /= sum(ww)
            ww1 = cumsum(ww)
            # this aux var  will store sum of simulated Y_i values
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
            end #end nY
            Ytotal[nc1] = Ysum
        end # end if j>nBurn && j%nfreq==0
    end # end for j
    #####     MCMC  ENDS HERE  #####
end_mcmc = Dates.format(now(), "HH:MM")

    println("End of MCMC---> now storing output")
    print("nc1 is: ")
    println(nc1)

    # ouput results here
    len = length(alphk)
    alphout = alphk[2:len]
    KKout = Nclust[2:len]
    Nout = Nt[2:len]
    #  output results for assessing convergence of MCMC
    # outfilealfa = open("alfa.txt","w")
    # outfileNclus = open("Nclus.txt","w")
    # outfileNj = open("Nj.txt","w")
    # writedlm(outfilealfa,alphout)
    # writedlm(outfileNclus,KKout)
    # writedlm(outfileNj,Nout)
    # close(outfilealfa)
    # close(outfileNclus)
    # close(outfileNj)
    ## return Ytotal
    ### CORRECT  Ytotal here #############
    #return(Ytotal)

    model = DataFrame(Total = float(Nsim)/float(Rsim))*Ytotal,
                      "Total running time" = init - end_mcmc)
    return (
end

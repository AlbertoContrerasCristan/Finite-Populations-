#include <Rcpp.h>
using namespace Rcpp;


//computes de sequence 1/xi_k = 2^{k}
//which  is used to reduce correlation between u_i's and w_k's

NumericVector InvXi(int m){ 
  NumericVector kappa(m); 
  for(int j=0; j<m; j++){
    kappa[j] = 2^(j+1);
  }
  return kappa;
}

//N_i SLICE SAMPLER
//determines the set of integers N_i in slice sampler
//Algorithm (see Kalli, et al , 2011)
// if uv[i] > 0.5 returns a 0
NumericVector Nis(NumericVector uv , int nu){ 
  int n = uv.size();
  NumericVector Ni(n); //Ni = [i for i in 1:nu]
  for(int i = 0; i<nu;i++){
    double u = uv[i];
    int j = 0;
    double Xhi = 0.5;
    while (Xhi > u){
      j += 1;
      Xhi *= 0.5;  
    }//End While
    Ni[i] = j;
  }//End for
    return Ni;

//random Beta size n
NumericVector rbeta(int n, double a, double b){
  gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937);// Create RNG seed
  NumericVector beta(n);
  for(int i = 0; i < n; i ++){
    beta[i] = gsl_ran_beta(s,a,b);
  }
  return beta;
}

//WEIGHTS SAMPLE
//produces a sample of the weights
//(w_1, ... ,w_Nn)  in Sethuraman's representation for D.P.
//beta1 = Beta(1.0,PPr1)
//Nu = rand(beta1,Nn)
NumericVector Weight0(int n, double b){ 
  NumericVector Nu = rbeta(n,1,b);
  v = 1.0.-Nu;
  for(int i = 1, i < n;i++){
    Nu[i] *= R::prod(v[1:i-1]);
  }// End for  
  return Nu;
}
      
// GAMMA MIXTURES
// computes mixture of gamma densities
// as in eq (13) of Escobar & West (1995)
// THIS VERSION TO GENERATE INITIAL VALUES FOR MCMC

NumericVector runif(int n, double a, double b){
  gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937);// Create RNG seed
  NumericVector unif(n);
  for(int i = 0; i < n; i ++){
    beta[i] = gsl_ran_flat(s,a,b);
  }
  return unif;
}

NumericVector rgamma(int n, double a, double b){
  gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937);// Create RNG seed
  NumericVector gamma(n);
  for(int i = 0; i < n; i ++){
    beta[i] = 1/gsl_ran_gamma(s,a,b);
  }
  return gamma;
}


NumericVector GaMix0(int MM,double fn,double eta){
  aux = (fn - 1)/(fn*(-log(eta)));  
  aux /= (1.0 + aux);//piw = aux/(1.0+aux)
  NumericVector alfa(MM); 
  for(int i = 0; i<5;i++){
    u0 = runif(1,0,1);
    if (u0 < aux){
      alpha = rgamma(1,fn,-log(eta)); 
    }else{
      alpha = rgamma(1,fn-1,-log(eta));
    } //End If
    alfa[ii] = alpha;
  }//End for 
  return mean(alfa);
}// End function

// computes mixture of gamma densities
// as in eq (13) of Escobar & West (1995)
NumericVector GaMix(double aa,double bb, double fk,double fn,double etta){ 
  aux = (aa + fk - 1)/(fn*(bb-log(etta)));
  aux /= 1.0+aux; // piw = aux/(1.0+aux)
  u2 = runif(1,0,1.);
  if(u2 < aux){ 
    alpha = rgamma(1,aa+fk,bb-log(etta));
  }else{
    alpha = rgamma(1,aa+fk-1.0,bb-log(etta));
  }
  return alpha[0];
}
                
// simulates from beta distribution for
// latent var eta  as in eq. (14)  Escobar & West (1995)
//eta = rbeta(1,alfa+1.0,fn)
//return eta[1]

NumericVector rbeta(int n, double a, double b){
  gsl_rng *s = gsl_rng_alloc(gsl_rng_mt19937);// Create RNG seed
  NumericVector beta(n);
  for(int i = 0; i < n; i ++){
    beta[i] = 1/gsl_ran_beta(s,a,b);
  }
  return beta;
}
              
double Beta(double alfa,double fn){
  return rbeta(1,alfa+1.0,fn);
}
                  
                
//NIG FAMILIES
//Generates random vectors (mui,Vi) from the bi-variate prior
//Normal-InverseGamma((mui,Vi) |  m,tau,S,df)  distribution
// as in Escobar & West (1995)

NumericMatrix NIG0(double m,double tau,double S,double df::Float64,int Nn){
  beta = S/2;
  alpha = df/2;
  V = rgamma(Nn,alpha,beta);
  V = 1/V;
  //conditionally on Vi generate  mu_i from N(m, tau*Vi)
  sig = sqrt.(tau.*V);
  NumericVector mu(Nn);
  for(i=0; i<Nn; i++){ 
    mu[i] = rnorm(1,m,sig[i]);
  }//End for 
  //return a Nn*2 matrix with each row given by (mui,Vi)
  return NumericMatrix(mu,v);
}//End NIG0
                    
// Generates random vectors (mui,Vi) from the bi-variate
// posterior Normal-InverseGamma( (mui,Vi) | (df+1)/2, Sx/2)
// distribution (see Escobar & West 1995,sect 3)
NumericMatrix NIG(double x,double m1,double tau1, double S1, double df1){
  Sx = (((x-m1)^2)/(1+tau1))+ S1;
  Beta = 0.5*Sx;
  alpha = 0.5*(df1+1);
  // generate from Gamma(alpha,Beta) then
  // transform to Inverse-Gamma(alpha,Beta)
  Vi = rgamma(1,alpha,Beta);
  Vi=1/Vi;
  // conditionally on Vi generate from N(c,(Vi[tau/(tau+1)])^{2})
  mu = (tau1*x+m1)/(1+tau1);
  sig = sqrt((tau1*Vi)/(1+tau1));
  mui = rnorm(1,mu,sig);
  //return hcat(mui,Vi)
  return  NumericMAtrix(mui,Vi);
}//End NIG
                      
// Generates random vectors (mui,Vi) from the bi-variate
// POSTERIOR Normal-InverseGamma((mui,Vi) | x, m,tau,S,df)  distribution
// x = data. Nn1 = number of samples to generate
                      function NIG00(x::Array{Float64,1},m::Float64,tau::Float64,S::Float64,df::Float64,Nn1::Int64)
                        Theta = zeros(Nn1,2)
                        for i =1:Nn1
                          Theta[i,:]  = NIG(maximum(x),m,tau,S,df)
                          end
                          return Theta
                          end
                          
############################################################### D SAMPLER
                          
# Given vectors (w_1,...w_N) ; (u_1, ... , u_n)
# vTheta' = ((mu_1,Ka_1), ... ,(mu_N,Ka_N))  produces a
# sample from the discrete indicators (d1, ... , dn)
                          function dsamplerC(weights::Array{Float64,1},vTheta::Array{Float64,2},Nvect::Array{Int64,1},data1::Array{Float64,2})
                            nd=size(data1)[1]
                          dd =  [0 for i1=1:nd]
                          for i = 1:nd
                            ww = weights[1:Nvect[i]].*InvXi(Nvect[i])
                            dens = [pdf(Normal(vTheta[k,1],vTheta[k,2]),data1[i])  for k=1:Nvect[i]]
                          Pei = ww.*dens
                            Pei /= sum(Pei) # re-normalize by constant
                            Pei = cumsum(Pei) # sample di from this discrete dist.
                          u3 = runif(1,0.0,1.0)[1]
                          c1 = 1
                          while u3>Pei[c1]
                          c1 += 1
                          end
                            dd[i] = c1 # end for i
                              end
                              return dd
                              end
                              
# Given Xi sequence of weights and
# d1, ... , dn  samples values u1, ... , un
##univec = [0.0 for i2=1:nn]
##univec =  [runif(1,0.0,Xiv[dd[i]])  i=1:nn]
                              function unifvec(dd::Array{Int64,1},Xiv::Array{Float64,1},nn::Int64)
                                return [runif(1,0.0,Xiv[dd[i]])[1]  for i=1:nn]
                              end
                                
# given d1, ... , dn
# builds the sets Dj = {i: d_i = j};
                                function buildDj(dv::Array{Int64,1},nn::Int64,j::Int64)
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
                                  
# updates the value of Theta=(miu_j,V_j)
# according to Bayes Rule
                                  function FuConThet(yy::Array{Float64,2},em::Float64,taau::Float64,sdd::Float64,S1::Float64,nThet::Array{Float64,2},enj::Int64,Dj::Array{Int64,1})
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
                                    
# Given d1, ... ,  dn (dv1);  a dataset (dat1)
# hyperpars m (m1), tau (ttau), s (sd),S (S1);  produces posterior
# samples  (mu_1,V_1), ... (mu_newN, V_newN), (w_1, ... , w_newN)
                                    function ThetaWeight(dv1::Array{Int64,1},dat1::Array{Float64,2},ndat::Int64,m1::Float64,ttau::Float64,sd::Float64,S1::Float64,newN::Int64,Thetta::Array{Float64,2},Nold::Int64,M::Float64)
                                      Nu = [0.0 for i=1:newN]
                                    weights = [0.0 for i=1:newN]
                                    Kl  = 0 # this counter will register number of clusters
                                      for j1=1:newN
                                        listDj = buildDj(dv1,ndat,j1)
                                        D = listDj[1]
                                      njj = listDj[2]
                                      mjj = listDj[3]
                                      if  j1 > Nold # update Thetta_j1
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
                                                  Nu[j1] = rbeta(1,1.0+njj,M+mjj)[1] # update v_j1
                                                  if j1 > 1 # update weights_j1
                                                    weights[j1] = Nu[j1]*prod(1.0.-Nu[1:j1-1])
                                                    end #end for j1
                                                      end
                                                      if newN <= Nold
                                                        Thetta = Thetta[1:newN,:]
                                                      end
                                                        weights[1] = Nu[1]
                                                      return Thetta,weights,Kl
                                                        end

                                                        
                                                        function NormMixMCMC(aa::Float64,
                                                                             bb::Float64,
                                                                             sdf::Float64,
                                                                             sS::Float64,
                                                                             mm::Float64,
                                                                             tauu::Float64,
                                                                             nMCMC::Int64,
                                                                             nBurn::Int64,
                                                                             nfreq::Int64,
                                                                             PP::Float64,
                                                                             dd::Array{Int64,1},
                                                                             Theta::Array{Float64,2},
                                                                             data::Array{Float64,2},
                                                                             tN::Int64)
# Given initial value Theta'(Theta_1,...Theta_Nlast) = ((mu1,Ka1) ... (mu_Nlast,Ka_Nlast))
# and selectors  d1,...,dn, performs MCMC.
                                                        init = now()
                                                        sizeMCMC = Int64(nMCMC/nfreq)
                                                        print("mcmc running and to store")
                                                        print(" ",sizeMCMC," posterior samples...","\n")
                                                        n = size(data)[1]
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
                                                              println("End of MCMC---> now storing output")
                                                              print("nc1 is: ")
                                                              println(nc1)
# ouput results here
                                                              len = length(alphk)
                                                              alphout = alphk[2:len]
                                                            KKout = Nclust[2:len]
                                                            Nout = Nt[2:len]
                                                            finish = now()
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
##return Ytotal
### CORRECT  Ytotal here #############
#return(Ytotal)
                                                              
#runtime = init - finish
#model = Array{Any}(undef,5)
#model[2] = alphout
#model[3] = KKout
#model[4] = Nout
#model[5] = runtime, "HH:MM"
#model[1] = (float(Nsim)/float(Rsim))*Ytotal
                                                              return (float(Nsim)/float(Rsim))*Ytotal
                                                                end                                                        

                                                        
        
// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {
  return x * 2;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/

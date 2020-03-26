# Epidemiology models
# model: multi-strain SIRS

## cross immunity: slighly diff than the Gog model, sirs model fast version
cross.sirs.fast<-function(tm_strt, tm_end, tm_step, S0, I0, N, D, L, beta, cross, birthrate=0, realdata=FALSE){
  # function to integrate to the next time step
  # use multistrain SIRS model, integrate dicretely with Poisson distributions
  # input: tm_strt: starting time; tm_end: ending time; tm_step: time step
  #         S0, I0: initial states; N: population size
  #         D: infection period, day; matrix(num_obs,num_ens)
  #         L: immune period, day; matrix(num_obs,num_ens)
  #         alpha: rate from exposed to infectious; 
  #         beta: transmission matrix at time t; matrix(num_obs,num_ens)
  #         cross: cross immunity matrix, c11=c22=c33=0; 
  #         c12:infected by strain 2, provide cross-protection to strain 1
  #         NOTE: C12 may be equal to C21, 0<cij<1
  #         CROSS=cross: matrix[diff cij, num_ens]
  #         birthrate: per day per population
  #
  # output: S, I for all time steps
  cnt=1;
  # beta stores only data during the time used for the truth
  tm_strt=tm_strt-tm.range[1]+1; # adjust the index to match beta
  tm_end=tm_end-tm.range[1]+1;
  tm_vec=seq(tm_strt,tm_end,by=tm_step)
  tm_sz=length(tm_vec)+1; # including the initial conditions and results of length(tm_vec) integration steps
  Ns=dim(I0)[1]; # number of strains
  Np=dim(I0)[2]; # number of ensemble members
  # if integrating parallelly, S0, and I0 should be passed in as a vector (Np particles)
  
  S=I=newI=array(0,c(Ns,Np,tm_sz))
  S[,,1]=S0; I[,,1]=I0; # R[,1]=N-S0-I0;
  newI[,,1]=0;
  
  # also track cummulative cross-imm 
  cross.newI=array(0,c(Ns,Np,tm_sz))
  
  # cross immunity: slighly diff than the Gog model
  # Si= -beta_i*Si*Ii/N - sum(beta_i * cimm_i_j * Sj * Ij/N) 
  # note the infection term by other strain here is Sj*Ij/N (cp Si*Ij/N in Gog & Grenfell)
  
  cross.sums=rep('0',Ns);
  for(si in 1:Ns){  # compute cross imm for each strain
    for(sj in 1:Ns){
      if(si==sj) next;
      cross.sums[si]=paste0(cross.sums[si],"+cross['cimm",si,'_',sj,"',]*Einf[",sj,',]')
      
    }
  }
  
  for (t in 1:length(tm_vec)){
    cnt=cnt+1;
    
    Einf=tm_step*(beta*I[,,cnt-1]*S[,,cnt-1]/N)
    
    Ecross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      Ecross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    
    Erecov=tm_step*(1/D*I[,,cnt-1])
    Eimmloss=tm_step*(1/L*(N-S[,,cnt-1]-I[,,cnt-1])) 
    # check aphyiscality: new infection <= availabel susceptible
    Einf=pmin(Einf,S[,,cnt-1])
    Ecross=pmin(Ecross,S[,,cnt-1]-Einf)
    Erecov=pmin(Erecov,I[,,cnt-1])
    Eimmloss=pmin(Eimmloss,N-S[,,cnt-1]-I[,,cnt-1])
    smci=Einf
    smcc=Ecross
    smcr=Erecov
    smcl=Eimmloss
    
    sk1=smcl-smci-smcc
    ik1=smci-smcr
    ik1a=smci; # new infections
    ik1c=smcc; # loss of suscepbility due to cross immunity
    
    Ts1=S[,,cnt-1]+sk1/2
    Ts1[Ts1<0]=0
    Ti1=I[,,cnt-1]+ik1/2
    
    Einf=tm_step*(beta*Ti1*Ts1/N)
    Ecross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      Ecross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    Erecov=tm_step*(1/D*Ti1)
    Eimmloss=tm_step*(1/L*(N-Ts1-Ti1)) 
    # check aphyiscality: new infection <= availabel susceptible
    Einf=pmin(Einf,Ts1)
    Ecross=pmin(Ecross,Ts1-Einf)
    Erecov=pmin(Erecov,Ti1)
    Eimmloss=pmin(Eimmloss,N-Ts1-Ti1)
    smci=Einf
    smcc=Ecross
    smcr=Erecov
    smcl=Eimmloss
    
    sk2=smcl-smci-smcc
    ik2=smci-smcr
    ik2a=smci;
    ik2c=smcc; # loss of suscepbility due to cross immunity
    Ts2=S[,,cnt-1]+sk2/2
    Ts2[Ts2<0]=0
    Ti2=I[,,cnt-1]+ik2/2
    
    Einf=tm_step*(beta*Ti2*Ts2/N)
    Ecross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      Ecross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    Erecov=tm_step*(1/D*Ti2)
    Eimmloss=tm_step*(1/L*(N-Ts2-Ti2)) 
    # check aphyiscality: new infection <= availabel susceptible
    Einf=pmin(Einf,Ts2)
    Ecross=pmin(Ecross,Ts2-Einf)
    Erecov=pmin(Erecov,Ti2)
    Eimmloss=pmin(Eimmloss,N-Ts2-Ti2)
    smci=Einf
    smcc=Ecross
    smcr=Erecov
    smcl=Eimmloss
    sk3=smcl-smci-smcc
    ik3=smci-smcr
    ik3a=smci;
    ik3c=smcc; # loss of suscepbility due to cross immunity
    Ts3=S[,,cnt-1]+sk3
    Ts3[Ts3<0]=0
    Ti3=I[,,cnt-1]+ik3
    
    Einf=tm_step*(beta*Ti3*Ts3/N)
    Ecross=matrix(0,Ns,Np)
    for(si in 1:Ns){
      Ecross[si,]=eval(parse(text=eval(parse(text=paste('cross.sums[',si,']',sep='')))))
    }
    Erecov=tm_step*(1/D*Ti3)
    Eimmloss=tm_step*(1/L*(N-Ts3-Ti3))
    # check aphyiscality: new infection <= availabel susceptible
    Einf=pmin(Einf,Ts3)
    Ecross=pmin(Ecross,Ts3-Einf)
    Erecov=pmin(Erecov,Ti3)
    Eimmloss=pmin(Eimmloss,N-Ts3-Ti3)
    smci=Einf
    smcc=Ecross
    smcr=Erecov
    smcl=Eimmloss
    sk4=smcl-smci-smcc
    ik4=smci-smcr
    ik4a=smci;
    ik4c=smcc; # loss of suscepbility due to cross immunity
    
    seed=rep(.1,Ns);
    S[,,cnt]=S[,,cnt-1]+sk1/6+sk2/3+sk3/3+sk4/6-seed + tm_step*birthrate*(N-S[,,cnt-1]) # add births/deaths
    S[S<0]=0
    I[,,cnt]=I[,,cnt-1]+ik1/6+ik2/3+ik3/3+ik4/6+seed - tm_step*birthrate*I[,,cnt-1]
    newI[,,cnt]=newI[,,cnt-1]+ik1a/6+ik2a/3+ik3a/3+ik4a/6+seed;
    cross.newI[,,cnt]=cross.newI[,,cnt-1]+ik1c/6+ik2c/3+ik3c/3+ik4c/6;

  }
  
  if (realdata==FALSE){
    rec=list(S=S,I=I); 
  } else {
    rec=list(S=S,I=I,newI=newI,cross.newI=cross.newI); 
  }
  rec;
}


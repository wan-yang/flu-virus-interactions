## 9/27/18
# Script of the multi-strain SIRS - IF2 filter 
# For Yang, Lau, and Cowling 2020 Dynamic interactions of influenza viruses in Hong Kong during 1998-2018
## USE THE MIF2 (IONIDES ET AL 2015)
## DIFFERENCE CP MIF1: perturbation is applied to each filtering step & the particles are recycled each iteration
##  as opposed to only use the mean & variance
## IF2 algorithm

# logit function and its inverse
fn_logit=function(x) {
  # if(any(x==0)) x=x+1e-6; 
  # if(any(x==1)) x=x-1e-6;
  x[x==0]=1e-6;
  x[x==1]=1-1e-6
  log(x/(1-x))
}
fn_expit=function(x) {1/(1+exp(-x))}
# function to add perturbation
fn_perturb=function(theta_tt,V.perturb_tt,v.logit,v.exp,delta_iter,v.parms,SD1.trans){  # add the perturbation
  # theta_tt: parms at time tt
  # V.perturb_tt: vector of the variance for each parm for time tt
  # delta_iter: delta for iteration iter
  # step 1: transform the parms
  # v.logit=c(grep('sf',rownames(theta_tt)),grep('cimm',rownames(theta_tt))) # those need logit transform
  # v.all=1:nrow(theta_tt)
  # v.exp=v.all[!(v.all %in% v.logit)]
  num_parm=length(v.parms)
  rownames(theta_tt)=v.parms;
  theta_tt.trans=theta_tt; # num_parm, num_ens
  theta_tt.trans[v.logit,]=fn_logit(theta_tt[v.logit,])
  theta_tt.trans[v.exp,]=log(theta_tt[v.exp,])
  # step 2: add perturbation
  sd0=0; 
  if(!is.null(SD1.trans)) sd0=SD1.trans*sqrt(V.perturb_tt)
  # in case it's stucked and converge to a single value s.t. sd=0
  xperturb=t(mvrnorm(num_ens,mu=rep(0,num_parm),Sigma=diag(delta_iter^2*V.perturb_tt)) %*% diag(apply(theta_tt.trans,1,sd) + sd0 ))
  theta_tt.trans=theta_tt.trans+xperturb 
  # step 3: inverse transform
  theta_tt.new=theta_tt
  theta_tt.new[v.logit,]=fn_expit(theta_tt.trans[v.logit,])
  theta_tt.new[v.exp,]=exp(theta_tt.trans[v.exp,])
  theta_tt.new
}
# function to check data aphysicality
Fn_checkDA<-function(xnew,bound.low,bound.up){
  b.low=bound.low;
  b.up=bound.up;
  n.var=nrow(xnew); n.ens=ncol(xnew);
  for(vi in 1:n.var){
    #  Corrects if <b.low
    ug=min(xnew[vi,]);
    if (ug<b.low[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]<b.low[vi]){
          xnew[vi,jj]=b.low[vi];
        }
      }
    }
    ug=max(xnew[vi,]);
    if (ug>b.up[vi]){  
      for (jj in 1:n.ens){
        if (xnew[vi,jj]>b.up[vi]){
          xnew[vi,jj]=b.up[vi];
        }
      }
    }
  }
  xnew;
}

IF2<-function(num_ens, tmstep, param.bound, param.bound.full=param.bound.full, # also pass the full param.bound for changes
               obs_i=obs_i, obs_vars=obs_vars,tm.ini=273, tm.range=273:500,
              epi.model='cross.sirs.fast',
              IN=30,v.perturb0=0.02,v.perturb.bp=.1,regul=T){
  
  library("truncnorm"); library("tgp"); library('mvtnorm'); # for lhs
  library("MASS"); # for multivariate normal distribution
  source(paste(dir_home_code,"multistrainSIRS.R",sep="")); # SIR model
  
  # INITIALIZATION
  if (! exists('diff_L')) diff_L=FALSE  # whether to use diff L for diff strains
  if (! exists('adj_params')) adj_params=TRUE
  
  fn_epi=get(epi.model);
  
  num_obs=ncol(obs_i);  # number of strains included
  num_parm=(4)*num_obs+(num_obs^2-num_obs); # number parms: 1 for sf; 3 for the SIRS (R0, D, L); (num_obs^2-num_obs) for the cross immunity matrix
  num_var=(4+3)*num_obs+(num_obs^2-num_obs); # 4 states: S, I, newI, scaling; 3 for the SIRS (R0, D, L); 6 for the cross immunity matrix
  num_times=nrow(obs_i);  
  
  # parm bounds
  theta_low=param.bound.full[,1]; theta_up=param.bound.full[,2]
  
  # wider DA bounds
  {
    thetaDA_low=c(rep(param.bound['sf',1]*.5,num_obs),rep(param.bound['D',1]*.5,num_obs),rep(param.bound['L',1]*.5,num_obs),
                  rep(param.bound['R',1]*.5,num_obs),rep(0,(num_obs^2-num_obs))); 
    thetaDA_up=c(rep(param.bound['sf',2]*3,num_obs),rep(param.bound['D',2]*3,num_obs),rep(param.bound['L',2]*3,num_obs),
                 rep(param.bound['R',2]*3,num_obs),rep(1,(num_obs^2-num_obs)));
  }
  
  ## FOR THE IF2
  # IF Algorithmic parameters:
  # cooling: delta[1]=1, delta[100]=.1, delta=cooling^m => 
  cooling=.1^(1/100)  # 0.9772372
  delta=cooling^(1:IN-1) # delta: perturbation scale
  ll=rep(0,IN) # to save the likelihood
  
  
  THETA_POST=t(lhs(num_ens,cbind(theta_low,theta_up))) 
  # perturbation density h_n(theta|phi;delta) ~ N(phi,delta[m]^2*V[n]) 
  # v.perturb0=0.02 # for regular parms, for transformed parm with uncertainty on the order of 1 unit
  # v.perturb.bp=0.1 # break point
  V.perturb=rep(v.perturb0^2,num_times);
  
  flag=NULL; # record err
  
  # record timing of probing - potential antigenic changes
  tm.probe=array(0,c(IN,num_times,1+num_obs)); dimnames(tm.probe)[3]=list(c('time',paste0('strain',1:num_obs)))
  
  RES=NULL; # save the results from individual iteration?
  
  for(i in 1:IN) tm.probe[i,,1]=1:num_times
  
  ### This is for the particle filter regularization
  {
    ### This is for the regularization
    #  Getting the regularized component (random draw off a Epanechnikov kernel) is a bit complex.  
    # If n is the dimension of the state vector (could include the parameters as these are also adjusted), 
    # then the Epanechnikov kernel is K = (n+2)/(2C) * (1 - x^2) if |x|<1
    # and 0 otherwise.  C is the volume of a unit hypersphere in dimensions n.  
    # So for n=2 it is simply pi (i.e. pi*r^2, where r=1).  For n=6, it is pi^3/6.
    # here, for n=8 (even number), Vn(R)=Cn*R^n; Cn=pi^(n/2)/(n/2)!=pi^4/4!=pi^4/24.
    # for odd number: V2k+1=2k!(4pi)^k/(2k+1)!
    #  The optimal bandwidth for the kernel draw is
    #  h_opt=A/(N^(1/n+4)), where N is the number of particles
    #  and A = [(8/(C^(n+4)))*(2*sqrt(pi)^n)]^(1/n+4)
    # if use SEIR model:
    # C=pi^4/24;
    # A=((8/C^(8+4))*(2*sqrt(pi)^8))^(1/(8+4));
    # hopt=2*A/num_ens^(1/(8+4));
    # SIR model, 6 variables
    ## same famular as in Arulampalam et al. 2002
    # C=pi^3/6;
    # C=8/15*pi^2;  # 21 variables: S, I, sf, D, and R, and 6 cross immunity entries 
    
    # C=16/105*pi^3; # for n=7
    # A=((8/C^(6+4))*(2*sqrt(pi)^6))^(1/(6+4)); # typo?
    # hopt=2*A/num_ens^(1/(6+4));  # typo?
    nn=num_var-num_obs; # number of parameters
    C=pi^(nn/2)/gamma(nn/2+1);
    A=((8/C*(nn+4))*(2*sqrt(pi))^nn)^(-1/(nn+4)); # # change the exponent '(1/(6+4))' to (-1/(6+4))
    # hopt=A*num_ens^(-1/(6+4));  # checked right, from Density Estimation by Silverman
    hopt=2*A*num_ens^(-1/(nn+4));  # increase the band width
    xK=seq(-1,1,by=.001)
    dK=((nn+2)/2/C)*(1-xK^2); # density of the Epanechnickov Kernel
    dK=dK/sum(dK); # normolized
    # construct the cumulative density function
    cK=dK;
    for (i in 2:length(dK)){
      cK[i]=cK[i-1]+dK[i];
    }
  }
  # run IN iteration
  for (iter in 1:IN){
    
    print(paste('Iteration #',iter),quote=F)
    
    # initation for this iter 
    
    # put everything in x
    So=matrix(0,num_var+num_obs,num_ens); # S,I,newI,sf,D,R,L,each has num_obs col, 
    # and the off-diagonal entries for the cross-immunity matrix
    x=matrix(0,num_var+num_obs,num_ens);  # also store the R0 combining seasonality and scaling for that season
    xprior_mean=xpost_mean=xpost_sd=xprior_sd=NULL;
    rownames(x)=rownames(So)=c(paste0('S',1:num_obs),paste0('I',1:num_obs),paste0('newI',1:num_obs),
                               paste0('sf',1:num_obs),paste0('D',1:num_obs),paste0('L',1:num_obs),
                               paste0('R',1:num_obs),c('cimm1_2','cimm1_3','cimm2_1','cimm2_3','cimm3_1','cimm3_2'),paste0('R0tt',1:num_obs))
    v.parms=c(paste0('sf',1:num_obs),paste0('D',1:num_obs),paste0('L',1:num_obs),
              paste0('R',1:num_obs),c('cimm1_2','cimm1_3','cimm2_1','cimm2_3','cimm3_1','cimm3_2'))
    v.non.newi=rownames(x)[!grepl('newI',rownames(x))]
    v.non.i=rownames(x)[!grepl('I',rownames(x))]
    v.non.newi=v.non.newi[!grepl('R0tt',v.non.newi)]
    v.non.i=v.non.i[!grepl('R0tt',v.non.i)]
    
    # for SR
    # variables that the filter is allowed to probe
    v.SR=c(rownames(x)[grepl('S',rownames(x))])
    bounds.full=cbind(c(rep(S_low*N,num_obs),theta_low),c(rep(S_high*N,num_obs),theta_up))
    rownames(bounds.full)=v.non.i 
    
    # for transformation
    v.logit=c(v.parms[grepl('cimm',v.parms)])
    v.exp=v.parms[!(v.parms %in% v.logit)]
    
    # Initialization
    So[c(paste0('S',1:num_obs),paste0('I',1:num_obs)),]=t(lhs(num_ens,matrix(c(rep(S_low*N,num_obs),rep(0,num_obs),rep(S_high*N,num_obs),rep(100,num_obs)),num_obs*2,2)));  # per 100,000 population
    for(i in 1:num_obs){
      if(obs_i[1,i]>100){
        So[paste0('I',i),]=rnorm(num_ens,mean=obs_i[1,i],sd=obs_i[1,i]/3);
      }
    }
    So[So<0]=0;
    
    # depending on the iteration, inital parms are diff
    if(iter==1){
      So[(3*num_obs+1):(nrow(So)-num_obs),]=THETA_POST;
    } else {
      So[(3*num_obs+1):(nrow(So)-num_obs),]=fn_perturb(theta_tt = THETA_POST,V.perturb_tt=rep(V.perturb[1],num_parm),v.logit,v.exp,delta_iter = delta[iter],v.parms,SD1.trans=0)
    }
    
    # integrate forward one step
    # try using the seasonal cycle with lower-upder bounds
    R0sn.tt=t(lhs(num_ens, rect = cbind(c(R0sn_lows[1,]),c(R0sn_ups[1,]))))
    R0tt=So[paste0('R',1:num_obs),] * R0sn.tt
    So[paste0('R0tt',1:num_obs),]=R0tt;
    
    beta=R0tt / So[paste0('D',1:num_obs),];
    
    cross=So[grep('cimm',rownames(So)),]
    tcurrent=tm.ini;
    Sr_tmp=fn_epi(tm_strt=tcurrent+dt, tm_end=tcurrent+tmstep, tm_step=dt, S0=So[paste0('S',1:num_obs),],
                  I0=So[paste0('I',1:num_obs),], N=N, D=So[paste0('D',1:num_obs),], L=So[paste0('L',1:num_obs),],
                  beta=beta, cross=cross,birthrate = 1/75/365, realdata=T)
    x[paste0('S',1:num_obs),]=Sr_tmp$S[,,1+tmstep]  # tail(Sr_tmp$S,1);
    x[paste0('I',1:num_obs),]=Sr_tmp$I[,,1+tmstep];
    x[paste0('newI',1:num_obs),]=Sr_tmp$newI[,,1+tmstep]*So[paste0('sf',1:num_obs),]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
    x[(3*num_obs+1):nrow(x),]=So[(3*num_obs+1):nrow(So),];
    
    cumlike=NULL; # to record cum likelihood of the time series
    # Particle Filtering
    SD1=apply(x,1,sd);
    if(iter==1){
      SD1.trans=rep(0,length(v.parms)); names(SD1.trans)=v.parms
      SD1.trans[v.logit]=apply(fn_logit(x[v.logit,]),1,sd)
      SD1.trans[v.exp]=apply(log(x[v.exp,]),1,sd)
    }
    
    wts=nwts=rep(1/num_ens,num_ens);
    swts=matrix(1/num_ens,num_obs,num_ens); # record the weights for each strain
    for (tt in 1:num_times){
      for(si in 1:num_obs){
        swts[si,]=dnorm(x=x[paste0('newI',si),],mean=obs_i[tt,si],sd=sqrt(obs_vars[tt,si]));
      }
      # determine which strain(s) to probe based on the weights
      wts.max=apply(swts,1,max);
      wts.cut=dnorm(x=obs_i[tt,]+sqrt(obs_vars[tt,]),mean=obs_i[tt,],sd=sqrt(obs_vars[tt,])) # [8/17/18] looser cut?
      SR_sdx=NULL;
      for(si in 1:num_obs){
        if(wts.max[si]< wts.cut[si])
          SR_sdx=c(SR_sdx,si);
      }
      
      wts=wts*dmvnorm(x=t(x[paste0('newI',1:num_obs),]),mean=obs_i[tt,],sigma=diag(obs_vars[tt,],num_obs,num_obs))
      # Normalizing the weights here
      # sometimes obs_i=0, if so, density=0 for all particles, can nornalize (divide by 0)
      if (sum(wts)==0 | any(is.na(wts))){
        print(c(paste0('Iter=',iter,'; tt=',tt,', probe: '),'reset wts'),quote = F)
        nwts=wts=1/num_ens; # assign equal weights if get no information from the likelihood
        cumlike=append(cumlike,mean(wts));
        # record err: iteration, tm step
        flag=rbind(flag,c("degeneracy",iter,tt));
      } else {
        nwts=wts/sum(wts);  # Normalizing here
      }
      
      # save it first if it's the last time
      if(tt==num_times) {
        smp=sample(1:num_ens,size=num_ens,replace = T,prob=nwts)
        final.parms=x[v.parms,smp]
      }
      # Resampling with regulation
      neff=1/(sum(nwts^2));
      if(neff<num_ens/4){ # do not resample that often
        
        currx=x;   #  Getting current state and parameters
        
        ind_new=sample(x=1:num_ens, size=num_ens, replace=T, prob=nwts)
        # should record before resetting the weights!
        cumlike=append(cumlike,mean(wts)); # mean(wts[,tt])
        nwts=1/num_ens; # reset the weights to equal after resampling.
        wts=1/num_ens;
        
        currx=currx[,ind_new];
        
        
        if(regul==TRUE){
          ## add regularization noise
          SD=apply(currx,1,sd);
          
          nze=matrix(0,nrow(x),num_ens); rownames(nze)=rownames(x)
          
          for (i in v.non.newi){ # skip newI
            nze[i,]=approx(cK,xK,runif(num_ens))$y;  # Regularize noise ?
          }
          
          currx=currx+hopt*diag(SD,length(SD))%*%nze; # small jitter
        }
        
        ## check DA aphyiscality
        currx[v.non.newi,]=Fn_checkDA(currx[v.non.newi,],bound.low=c(rep(0,num_obs*2),thetaDA_low),bound.up=c(rep(N,num_obs),rep(N*.5,num_obs),thetaDA_up));

        x=currx;  
      }
      xmn=xsd=NULL;
      for(vi in 1:nrow(x)){
        xmn=c(xmn,sum(x[vi,]*nwts));
        xsd=c(xsd,sqrt(sum(nwts*(x[vi,]-sum(x[vi,]*nwts))^2)));
      }
      
      # [4/8/19] compute Re as well
      tmpRe=x[paste0('S',1:num_obs),]/N*x[paste0('R0tt',1:num_obs),]
      for(vi in 1:num_obs){
        xmn=c(xmn,sum(tmpRe[vi,]*nwts))
        xsd=c(xsd,sqrt(sum(nwts*(tmpRe[vi,]-sum(tmpRe[vi,]*nwts))^2)));
      }
      names(xmn)=names(xsd)=c(rownames(x),paste0('Re.tt',1:num_obs))
      xpost_mean=rbind(xpost_mean,xmn);
      xpost_sd=rbind(xpost_sd,xsd);
      # replace S within particles with the 10% lowest weights
      # x[1,order(nwts[,tt])[1:ceiling(S_adj*num_ens)],tt]=sample(seq(.6*N,.8*N,by=1),ceiling(S_adj*num_ens))
      
      ## pertubation of parms before perdiction
      V.perturb_tt=rep(v.perturb0^2,num_parm);
      if(!is.null(SR_sdx)){
        for(vi in SR_sdx){
          V.perturb_tt[grep(vi,v.parms)]=v.perturb.bp^2
        }
      }
      x[v.parms,]=fn_perturb(theta_tt = x[v.parms,],V.perturb_tt,v.logit,v.exp,delta_iter = delta[iter],v.parms,SD1.trans)
      
      # SR: randomly replace S
      if(adj_params==TRUE & any(SR_sdx)){
        print(c(paste0('Iter=',iter,'; tt=',tt,', probe: '),SR_sdx),quote=F);
        # record timing of probing
        tm.probe[iter,tt,1+SR_sdx]=1
        
        if(F){
          SR_vdx=c(1,4,5,6); # probe S, sf, D, R
          SR_vind=NULL;
          for(vi in SR_vdx){
            SR_vind=c(SR_vind,num_obs*(vi-1)+SR_sdx)
          }
          # probe the cross immunity matrix as well
          SR_vind=c(SR_vind,num_obs*6+1:6)
          SR_bound=matrix(c(rep(c(S_low*N,S_high*N),length(SR_sdx)),
                            rep(c(sf_low,sf_high),length(SR_sdx)),
                            rep(c(D_low,D_high),length(SR_sdx)),
                            rep(c(R_low,R_high),length(SR_sdx)),
                            rep(c(imm_low,imm_high),6)),length(SR_sdx)*length(SR_vdx)+6,2,byrow=T)
        }
        
        # probe all those related to SR_sdx strain(s) - those identified as needed
        SR_vind=NULL;
        for(vi in SR_sdx){
          SR_vind=c(SR_vind,v.SR[grepl(vi,v.SR)])
        }
        SR_vind=unique(SR_vind)
        
        SR_bound.full=bounds.full[SR_vind,,drop=F] # full range as in the prior
        
        # get the local bounds, around the mean
        tmp.mean=colMeans(tail(xpost_mean[,rownames(x)],3)) # get the mean in the past three time steps
        bounds.local=cbind(tmp.mean*.8,tmp.mean*1.2); rownames(bounds.local)=rownames(x)
        bounds.local=bounds.local[v.non.i,]
        # wider ~local
        bounds.local.wider=cbind(tmp.mean*.5,tmp.mean*2); rownames(bounds.local.wider)=rownames(x)
        # bounds.local.wider=cbind(bounds.full[,1]*1.2,bounds.full[,2]*.9)
        bounds.local.wider=bounds.local.wider[v.non.i,]
        # check the bound and make sure they don't go over the upper/lower
        bounds.local.wider[,1]=pmax(bounds.local.wider[,1],bounds.full[,1])
        bounds.local.wider[,2]=pmin(bounds.local.wider[,2],bounds.full[,2])
        SR_bound.local.wider=bounds.local.wider[SR_vind,,drop=F]
        # allow more probing of S
        v.s=rownames(SR_bound.local.wider)[grep('S',rownames(SR_bound.local.wider))]
        SR_bound.local.wider[v.s,]=SR_bound.full[v.s,,drop=F]
        
        bounds.local[,1]=pmax(bounds.local[,1],bounds.full[,1])
        bounds.local[,2]=pmin(bounds.local[,2],bounds.full[,2])
        SR_bound.local=bounds.local[SR_vind,,drop=F]
        SR_bound.local[v.s,]=bounds.full[v.s,]
        
        
        # only allow a very small % to probe the entire space! others, probe locally
        n.full=pmin(ceiling(SR_adj*.1*num_ens),100); 
        n.local.wider=pmin(ceiling(SR_adj*.1*num_ens),200);  
        n.local=ceiling(SR_adj*num_ens)-n.full-n.local.wider
        # SR_ind.full=sample(1:num_ens,n.full); # probe entire space
        SR_ind=sample(1:num_ens,ceiling(SR_adj*num_ens));
        SR_ind.full=head(SR_ind,n.full); 
        SR_ind.local.wider=SR_ind[n.full+1:n.local.wider]
        SR_ind.local=tail(SR_ind,n.local)
        x[SR_vind,SR_ind.full]=t(lhs(n.full,SR_bound.full));
        x[SR_vind,SR_ind.local.wider]=t(lhs(n.local.wider,SR_bound.local.wider));
        # SR_ind.local=sample(1:num_ens,n.local); # probe locally
        x[SR_vind,SR_ind.local]=t(lhs(n.local,SR_bound.local));
      }
      
      # integrate forward 1 step parallelly
      tcurrent=tm.ini+tmstep*tt;
      # integrate for 1 week
      # include seasonality
      R0sn.tt = t(lhs(num_ens, rect = cbind(c(R0sn_lows[WEEKS[tt],]),c(R0sn_ups[WEEKS[tt],]))))
      R0tt=x[paste0('R',1:num_obs),] * R0sn.tt
      x[paste0('R0tt',1:num_obs),]=R0tt;
      beta=R0tt / x[paste0('D',1:num_obs),];
      
      cross=x[grep('cimm',rownames(x)),]
      
      Sr_tmp=fn_epi(tm_strt=tcurrent+dt, tm_end=tcurrent+tmstep, tm_step=dt, S0=x[paste0('S',1:num_obs),],
                    I0=x[paste0('I',1:num_obs),], N=N, D=x[paste0('D',1:num_obs),], L=x[paste0('L',1:num_obs),],
                    beta=beta, cross=cross,birthrate = 1/75/365, realdata=T)
      x[paste0('S',1:num_obs),]=Sr_tmp$S[,,1+tmstep]  # tail(Sr_tmp$S,1);
      x[paste0('I',1:num_obs),]=Sr_tmp$I[,,1+tmstep];
      x[paste0('newI',1:num_obs),]=Sr_tmp$newI[,,1+tmstep]*x[paste0('sf',1:num_obs),]; #  %*% diag(So[paste0('sf',1:num_obs),ii])
      x[(3*num_obs+1):nrow(x),]=x[(3*num_obs+1):nrow(x),];
      
      
      xmn=xsd=NULL;
      for(vi in 1:nrow(x)){
        xmn=c(xmn,sum(x[vi,]*nwts));
        # xsd=c(xsd,sd(x[vi,]*nwts));
        xsd=c(xsd,sqrt(sum(nwts*(x[vi,]-sum(x[vi,]*nwts))^2)));
      }
      xprior_mean=rbind(xprior_mean,xmn);
      xprior_sd=rbind(xprior_mean,xsd);
    } # end for-loop: particle filter
    
    # SAVE RESULTS FOR THIS ITERATION
    
    THETA_POST=final.parms;
    
    if (is.null(cumlike)){
      ll[iter]=log(mean(wts)); # no resample done
    } else{
      ll[iter]=sum(log(cumlike),log(mean(wts)))
    }

    # save best fits
    {
      tstep=seq(tm.ini,num_times+tm.ini-1,by=1); # time indices
      colnames(xprior_mean)=colnames(xprior_sd)=rownames(x)
      
      # metrics for comparison
      Y=xpost_mean[,paste0('newI',1:num_obs)];  # newI for each strain
      sf=xpost_mean[,paste0('sf',1:num_obs)];
      Ytot=Y/sf;
      R0=xpost_mean[,paste0('R0tt',1:num_obs)];
      Re=xpost_mean[,paste0('S',1:num_obs)]*R0/N;
      colnames(Re)=paste0('Re',1:num_obs)
      
      Y.tot.ili=rowSums(Y); Y.tot.inf=rowSums(Ytot)
    
      # save it?
      if(saveall==T) {
        tmp=cbind(Iter=iter,time=tstep,xpost_mean,xpost_sd)
        colnames(tmp)=c('Iter','time',colnames(xpost_mean),paste0(colnames(xpost_sd),'.sd'))
        RES=rbind(RES,tmp)
      }
      
    }
    
    if(iter==1){
      best.fit=cbind(iter=1,cumlike=ll[1],time=tstep,xpost_mean,Re,newI.tot=Y.tot.ili,Inf.tot=Y.tot.inf,xpost_sd); # ,newI.tot=Y.tot.ili,Inf.tot=Y.tot.inf
    } else {
      if(ll[iter]>best.fit[1,'cumlike']){
        best.fit=cbind(iter=iter,cumlike=ll[iter],time=tstep,xpost_mean,Re,newI.tot=Y.tot.ili,Inf.tot=Y.tot.inf,xpost_sd);
      }
    }
    
  } # end iteration
  
  # save final results
  tstep=seq(tm.ini,num_times+tm.ini-1,by=1); # time indices
  colnames(xprior_mean)=colnames(xprior_sd)=rownames(x)
  
  # metrics for comparison
  Y=xpost_mean[,paste0('newI',1:num_obs)];  # newI for each strain
  sf=xpost_mean[,paste0('sf',1:num_obs)];
  Ytot=Y/sf;
  R0=xpost_mean[,paste0('R0tt',1:num_obs)];
  Re=xpost_mean[,paste0('S',1:num_obs)]*R0/N;
  colnames(Re)=paste0('Re',1:num_obs)
  xpost_mean=cbind(xpost_mean,Re)
  
  t_on=which(onsets[,'iliiso.tot']==1);
  t_end=which(ends[,'iliiso.tot']==1);
  t_peak=which(peaks[,'iliiso.tot']==1);
  parms_on=xpost_mean[t_on,c(v.non.i,paste0('Re',1:num_obs))];
  parms_end=xpost_mean[t_end,c(v.non.i,paste0('Re',1:num_obs))];
  parms_peak=xpost_mean[t_peak,c(v.non.i,paste0('Re',1:num_obs))];
  
  corr=rms=delta_sum_newI=delta_pkwk=rep(0,num_obs);
  for(i in 1:num_obs){
    corr[i]=cor(Y[,i],head(obs_i[,i],num_times));
    rms[i]=sqrt(mean((Y[,i]-head(obs_i[,i],num_times))^2));
    delta_sum_newI[i]=sum(Y[,i])-sum(head(obs_i[,i],num_times));
    delta_pkwk[i]=which(Y==max(Y[,i]))[1]-which(obs_i[,i]==max(obs_i[,i]))[1];
  }
  Y.tot.ili=rowSums(Y); Y.tot.inf=rowSums(Ytot)
  
  out1=cbind(time=tstep,xpost_mean,newI.tot=Y.tot.ili,Inf.tot=Y.tot.inf);
  out2=cbind(time=tstep,xpost_sd);
  out3=t(c(rms,corr,delta_sum_newI));
  out4=cbind(wk_on=t_on,parms_on);
  out5=cbind(wk_peak=t_peak,parms_peak);
  out6=cbind(wk_end=t_end,parms_end);
  out7=cbind(time=tstep,xprior_mean[1:nrow(xpost_mean),]);
  colnames(out3)=c(paste("rms",1:num_obs,sep=''),paste("corr",1:num_obs,sep=''),
                   paste("delta_sum_newI",1:num_obs,sep=''));
  
  cumlike=data.frame(Iter=1:IN,cumlike=ll)
  # save tm.probe too - record when probe
  tmp=NULL
  for(i in 1:IN){
    idx=which(rowSums(tm.probe[i,,-1])>0)
    if(length(idx)>1){
      tmp=rbind(tmp,data.frame(Iter=i,tm.probe[i,idx,]))
    } else if(length(idx)==1){
      tmp=rbind(tmp,data.frame(Iter=i,t(tm.probe[i,idx,])))
    }
  }
  if(adj_params==F) tmp='SR off'
  if(adj_params==T & is.null(tmp)) tmp='SR on, but no need for it'
  
  colnames(best.fit)=c('Iter','cumlike','time',colnames(xpost_mean),'newI.tot','Inf.tot',paste0(colnames(xpost_sd),'.sd'))
  
  if(saveall==T) {
    # only save the top 5 iterations?
    top.iters=cumlike[order(cumlike$cumlike,decreasing = T)[1:5],'Iter']
    print(top.iters)
    RES=data.frame(RES)
    RES=RES[RES$Iter %in% top.iters,]
    out=list(states=out1,sd=out2,cumlike=cumlike,best.fit=best.fit,tm.probe=tmp,states.all=RES); 
  } else {
    out=list(states=out1,sd=out2,cumlike=cumlike,best.fit=best.fit,tm.probe=tmp); 
  }
}

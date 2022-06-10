NMAmodel <- function(multi, one_IF, randpcons, pcons, zellner) {
  # Likelihood
  m1 <- paste0("model {for(i in 1:NHtH)  {y[i]~dnorm(delta[i],w[i]) }")

  if (multi == TRUE) {
    m2 <- paste0("y[(NHtH+1):N]~dmnorm(delta[(NHtH+1):N],PREC[,]) ")
    m4 <- paste0("delta[(NHtH+1):N]~dmnorm(mean[(NHtH+1):N],K[,])")

    # Heterogeneity
    m5 <- paste0("  for(i in 1:(N-NHtH)) {  for(j in 1:(N-NHtH)) {K[i,j]<-precision*H[i,j] }} ")
  } else {
    m2 <- m4 <- m5 <- NULL
  }

  # Parameterization of the model
  m3 <- paste0("for(i in 1:NHtH) {delta[i]~dnorm(mean[i],precision) }")


  if (one_IF == TRUE) {
    # Parameterization of the means
    m6 <- paste0(" for(i in 1:N) {   mean[i] <- d[t[i]] - d[b[i]]+ x[i]*beta   } ")
    # SSVS
    m7a <- paste0(" beta~dnorm(0,T_inv)
    #Posterior inclusion probabilities
                  ")

    m7b <- paste0("gamma~dbern(1- pow(p.cons,1/p))")

    m7c <- paste("
    T<-pow(psi,2)*equals(gamma,0) + pow(c,2)*pow(psi,2)*equals(gamma,1)
    T_inv<-pow(T,-1)
           ")

    m7 <- paste(m7a, m7b, m7c)
  } else {
    # Parameterization of the means
    m6 <- paste0("  for(i in 1:N) {mean[i] <- d[t[i]] - d[b[i]] +inprod(x[i,],beta[])} ")
    # SSVS
    if (zellner == TRUE) {
      m7 <- paste0(" beta[1:p]~dmnorm(mean.b[1:p],T_inv[,])

        for(i in 1:p){ gamma[i]~dbern(1- pow(p.cons,1/p))}

        for(i in 1:p){
          mean.b[i]<-0
          for(j in 1:p){
            #Dg[i,j]<-equals(i,j)*psi[i]*equals(gamma[i],0) + equals(i,j)*c*psi[i]*equals(gamma[i],1)
            #inverse Dg

            Dg[i,j]<-equals(i,j)* pow(psi[i],-1) *equals(gamma[i],0) + equals(i,j)*pow(c*psi[i],-1)*equals(gamma[i],1)
            #R[i,j]<-N*sigma*ZTZ[i,j] #ZTZ from r
            #R inv
            #R_inv[i,j] <- pow(N,-1)*inv.sigma*ZTZ[i,j] #inv.siga
            R_inv[i,j] <- pow(N,-1)*sigma*ZTZ[i,j]
          }
        }
        #inv.sigma~dgamma(0.0001,0.0001)
        #sigma<-1/inv.sigma

        log.sigma~dunif(-10,10)
        sigma <- exp(log.sigma)

        T_inv <- Dg %*% R_inv %*% Dg
        #T<-Dg %*% R %*% Dg
        #T_inv[1:p,1:p]<-inverse(T[,])
                   ")
    } else {
      # independent Inconsistency factors (R = Identity matrix)
      m7 <- paste0("

        for(i in 1:p){
          gamma[i]~dbern(1- pow(p.cons,1/p))
          beta[i]~dnorm(0,T_inv[i])
          T_inv[i]<-pow(psi[i],-2)*equals(gamma[i],0) + pow(c*psi[i],-2)*equals(gamma[i],1)
        }
                     ")
    }
  }

  # Priors
  m8 <- paste0("  precision<-1/pow(sd,2)
  sd~dnorm(0,1)I(0,)")

  if (randpcons == TRUE) {
    m9 <- paste0("p.cons~dbeta(157,44)")
  } else {
    m9 <- paste0("p.cons<-", pcons)
  }

  m10 <- paste0("
  for(k in 1:(ref-1)) {  d[k] ~ dnorm(0,.0001)}
  for(k in (ref+1):NT) {  d[k] ~ dnorm(0,.0001)}
  d[ref]<- 0   }") # delete this  } if NMA estimates are presented
  # NMA estimates
  # m11 <- paste0("
  # for (cc in 1:(ref-1)) {
  #   Eff.ref[cc]<- d[cc] - d[ref]
  #   predEff.ref[cc] ~ dnorm(Eff.ref[cc],precision)
  # }
  # Eff.ref[ref]=0
  # predEff.ref[ref]=0
  # for (cc in (ref+1):NT) {
  #   Eff.ref[cc]<- d[cc] - d[ref]
  #   predEff.ref[cc] ~ dnorm(Eff.ref[cc],precision)
  # }
  # for (cc in 1:(NT-1)) {
  #   for (k in (cc+1):NT) {
  #     Eff[cc,k] <- d[k] - d[cc]
  #     predEff[cc,k] ~ dnorm(Eff[cc,k],precision)
  #   }
  # }
  #
  # order<- rank(d[])
  # for(k in 1:NT) {
  #   most.effective[k]<-equals(order[k],1)
  #
  #   for(j in 1:NT) {
  #     effectiveness[k,j]<- equals(order[k],j)
  #   }
  # }
  # for(k in 1:NT) {
  #   for(j in 1:NT) {
  #     cumeffectiveness[k,j]<- sum(effectiveness[k,1:j])
  #   }
  # }
  #
  # for(k in 1:NT) {
  #   SUCRA[k]<- sum(cumeffectiveness[k,1:(NT-1)]) /(NT-1)
  # }
  #
  # for(i in 1:NHtH) {
  #   D[i]<-w[i]*(y[i]-delta[i])*(y[i]-delta[i])  ## two-arm studies
  # }
  # for(i in (NHtH+1):N) {
  #   D[i]<-(y[i]-delta[i])*(y[i]-delta[i])*PREC[i-NHtH,i-NHtH]  ## multi-arm studies
  # }
  # D.bar<-sum(D[])   }")

  paste(m1, m2, m3, m4, m5, m6, m7, m8, m9, m10) # m11)
}

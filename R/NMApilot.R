NMApilot <- function(multi, one_IF, zellner) {
  m1 <- paste0(" model { for(i in 1:NHtH)  { y[i]~dnorm(delta[i],w[i]) }  ")
  if (multi == TRUE) {
    m2 <- paste0("y[(NHtH+1):N]~dmnorm(delta[(NHtH+1):N],PREC[,])")
    m4 <- paste0("delta[(NHtH+1):N]~dmnorm(mean[(NHtH+1):N],K[,]) ")
    m5 <- paste0("  for(i in 1:(N-NHtH)) {  for(j in 1:(N-NHtH)) {K[i,j]<-precision*H[i,j] }} ")
  } else {
    m2 <- NULL
    m4 <- NULL
    m5 <- NULL
  }

  m3 <- paste0("for(i in 1:NHtH) {delta[i]~dnorm(mean[i],precision) }")

  if (one_IF == TRUE) {
    m6 <- paste0("for(i in 1:N) {mean[i] <- d[t[i]] - d[b[i]]+ x[i]*beta}
 beta~dnorm(0,0.0001)
")
  } else {
    if (zellner) {
      m6 <- paste0("for(i in 1:N) {mean[i] <- d[t[i]] - d[b[i]] +inprod(x[i,],beta[])}

 beta[1:p]~dmnorm(mean.b[1:p],T_inv[,])

 for(i in 1:p){
      mean.b[i]<-0 #mean of beta
     for(j in 1:p){
             R[i,j]<-N*sigma*ZTZ[i,j] #ZTZ from r
             }
 }

 inv.sigma~dgamma(0.0001,0.0001)
 sigma<-1/inv.sigma

 T_inv[1:p,1:p]<-inverse(R[,])")
    } else {
      m6 <- paste0("for(i in 1:N) {mean[i] <- d[t[i]] - d[b[i]] +inprod(x[i,],beta[])}
 for(i in 1:p){beta[i]~dnorm(0, 0.0001)}

 ")
    }
  }

  m7 <- paste0("precision<-1/pow(sd,2)
  sd~dnorm(0,1)I(0,)

  for(k in 1:(ref-1)) {d[k] ~ dnorm(0,.0001)}
  for(k in (ref+1):NT) {d[k] ~ dnorm(0,.0001)}
  d[ref]<- 0
}")

  paste(m1, m2, m3, m4, m5, m6, m7)
}

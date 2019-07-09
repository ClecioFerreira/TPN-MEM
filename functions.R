##-----------------------------------------------------------------------------------------
#   paquete que contiene la densidad de la normal multivariada
##-----------------------------------------------------------------------------------------
library(mnormt)

##----------------------------------------------------------------------------------------------------------------
# se definen las funciones s(gamma)=1-gamma y t(gamma)=1+gamma, se pueden observa en página 2, ec.(1) con detalle. 
# Estas funciones son parte de la  definición de la distribución two-piece, nosotros tomaremos el caso particular
# definido anteriormente para hacer la simulación.
##----------------------------------------------------------------------------------------------------------------
#fn.s <- function(gamma) if(abs(gamma) < 1) 1-gamma else NA  
#fn.t <- function(gamma) if(abs(gamma) < 1) 1+gamma else NA  
fn.s <- function(gamma) if(abs(gamma) < 1) 1+gamma else NA        # Clecio
fn.t <- function(gamma) if(abs(gamma) < 1) 1-gamma else NA        # Clecio

##-----------------------------------------------------------------------------------------
#      Definimos dnorm y rnorm, solo por comodida
##-----------------------------------------------------------------------------------------
dgenerator <- dnorm
rgenerator <- rnorm

##----------------------------------------------------------------------------------------
#   función que nos entrega distintas impresiones 
##----------------------------------------------------------------------------------------
ecat <- function(x, digits=options()$digits, tail="\n"){ 
  cat(deparse(substitute(x)), ": ", sep="")
  if(is.matrix(x)) {
    cat("\n") 
    print(x) } 
  else  cat(format(x, digits=digits), tail)
  invisible(x)
}

##-------------------------------------------------------------------------------
##         Programamos la distribución two-piece, página 2, ec.(1)
##-------------------------------------------------------------------------------
dtp <- function(x, mu=0, sigma=1, gamma) 
{ # TP density for given 'generator'  
  z <- (x-mu)/sigma
  sP <- fn.s(gamma)
  sN <- fn.t(gamma)
  f. <- ifelse(z>0, dgenerator(z/sP), dgenerator(z/sN))
  2 * f. /(sigma*(sP+sN))
}
#---

##----------------------------------------------------------------------------------------------------------
# Generamos una muestra de la distribución two-piece, esto es posible ya que conocemos la representación
# estaocástica, se puede ver esta representacion en las  Propiedades de la distribución two-piece página 3
##----------------------------------------------------------------------------------------------------------

rntp <- function(n=1, mu=0, sigma=1, gamma) 
{
  sP <- fn.s(gamma)
  sN <- fn.t(gamma)
  p <- sP/(sP + sN)
  pm <- 2*rbinom(n, 1, p) - 1
  v <- abs(rgenerator(n, 0, 1))
  x <- pm* ifelse(pm>0, v*sP, v*sN)
  xi=mu + sigma * x
return(xi)
}

rttp <- function(n=1, mu=0, sigma=1, gamma,nu=2) 
{
  sP <- fn.s(gamma)
  sN <- fn.t(gamma)
  p <- sP/(sP + sN)
  pm <- 2*rbinom(n, 1, p) - 1
  v <- abs(rt(n,nu))
  x <- pm* ifelse(pm>0, v*sP, v*sN)
  xi=mu + sigma * x
return(xi)
}

#rtp=Y
##-------------------------------------------------------------------------------------------------------
# Generamos la muestra para nuestro modelo con error en las variables, este modelo se puede observar con detalle en la 
# página 5, ec(8),(9),(10),(11) y (12)
##-------------------------------------------------------------------------------------------------------
sample.MEM <- function(n=1, alpha, beta, mu.xi, sigma2.xi, sigma2.u, Sigma.e, gamma) # componentes de theta
{
  d <- length(alpha)
  if(length(beta) != d | ncol(Sigma.e) != d) stop("dimensions mismatch")
  z  <- matrix(NA, n, d+1)
  
  #R <- chol(Sigma.e)
  for(i in 1:n) {
    xi <- rntp(1, mu.xi, sqrt(sigma2.xi), gamma)
    #xi <- rtp[i]
    x <- xi + rnorm(1, 0, sqrt(sigma2.u))
    y <- alpha + beta * xi + rmnorm(1, rep(0,d), Sigma.e)
    z[i,] <- c(x, y)
  }
  return(z)
}

sample.tMEM <- function(n=1, alpha, beta, mu.xi, sigma2.xi, sigma2.u, Sigma.e, gamma,nu) # model with xi_i distrib. t with nu gl
{
  d <- length(alpha)
  if(length(beta) != d | ncol(Sigma.e) != d) stop("dimensions mismatch")
  z  <- matrix(NA, n, d+1)
  
  #R <- chol(Sigma.e)
  for(i in 1:n) {
    xi <- mu.xi + sqrt(sigma2.xi)*rt(1, nu)
    x <- xi + rnorm(1, 0, sqrt(sigma2.u))
    y <- alpha + beta * xi + rmnorm(1, rep(0,d), Sigma.e)
    z[i,] <- c(x, y)
  }
  return(z)
}

sample.ttpMEM <- function(n=1, alpha, beta, mu.xi, sigma2.xi, sigma2.u, Sigma.e, gamma,nu=2) # model with xi_i distrib. tTP with 2 gl
{
  d <- length(alpha)
  if(length(beta) != d | ncol(Sigma.e) != d) stop("dimensions mismatch")
  z  <- matrix(NA, n, d+1)
  
  #R <- chol(Sigma.e)
  for(i in 1:n) {
    xi <- rttp(1, mu.xi, sqrt(sigma2.xi), gamma,nu)
    #xi <- rtp[i]
    x <- xi + rnorm(1, 0, sqrt(sigma2.u))
    y <- alpha + beta * xi + rmnorm(1, rep(0,d), Sigma.e)
    z[i,] <- c(x, y)
  }
  return(z)
}

sample.MEMuniv <- function(n=1, alpha, beta, mu.xi, sigma.xi, sigma.u, sigma.e, gamma) # componentes de theta
{
  d <- 1
  z  <- matrix(0, n, 2)    
  for(i in 1:n) {
    xi <- rtp(1, mu.xi, sigma.xi, gamma)
    x <- xi + rnorm(1, 0, sigma.u)
    y <- alpha + beta * xi + rnorm(1, 0, sigma.e)
    z[i,] <- c(x, y)
  }
  return(z)
}



##----------------------------------------------------------------------------------------------------------
# Ahora programamos el algoritmo EM.
##----------------------------------------------------------------------------------------------------------


EM.NTP_MEM <- function(z) # ajuste NTP-MEM por algoritmo EM
{ # asumimos z es (n=nº de filas, (d+1)=nº de columnas), z_i =(x_i,y_i\T)\T, pag5.
  d <- ncol(z) - 1
  n <- nrow(z)
  x <- z[,1]
  y <- z[,2:(d+1)]
  z.bar <- mu <- apply(z, 2,mean)       # conjunto de valores iniciales pag11.
  x.bar <- mu.xi <- mu.x <- as.numeric(mu[1])       
  y.bar <- mu.y <- as.matrix(mu[2:(d+1)])
  S <- var(z)   # \Sigma                # valor unicial definido pag.11.
  #Sigma.e <- S[2:(d+1), 2:(d+1)]        # \Sigma esta definido como una partición pag.5.
  #Sigma.e <- var(y)# diag(apply(y,2,var))  # caso especial constante
  Sigma.e <-diag(var(y)) # caso Sigma.e diagonal
  sigma.x <- sigma.xi <- sqrt(S[1,1])   # ver pag.11.
  gamma <- 0 #mean((x-x.bar)^3)/sigma.x^3  # ver pag.11. 
  sigma2.u <- (sigma.x)^2           # valor inicial arbitrario.  
  sigma2.xi <- (sigma.x)^2
  Sigma <- diag(c(sigma2.u,Sigma.e))
  beta <- matrix(S[1, 2:(d+1)]/S[1,1])    # ver pag.11
  alpha <- matrix(y.bar) - x.bar*beta
  mu=matrix(0,d+1,1)
  mu[1]=mu.xi
  mu[2:(d+1)] <- alpha + beta*mu.xi     # ver pag.5, ec(14)
  b.xi=matrix(0,d+1,1)
  b.xi[1]=sigma.xi
  b.xi[2:(d+1)]=sigma.xi*beta
  Sigma.inv <- diag(c(1/sigma2.u,1/Sigma.e)) # caso diagonal
  log0=Q_gamma(gamma,z,mu, b.xi,Sigma,Sigma.inv)
  theta0=as.vector(c(as.vector(alpha),as.vector(beta),as.vector(Sigma.e),mu.xi,sigma2.xi,sigma2.u,gamma)) 
  criterio=1
  cont=0
  while (criterio>1e-4) {
  cont=cont+1
    #cat("------- k: " , k, "\n")        #cantidad de iteraciones que hace el prgrama
    #------------------------------ paso - E ------------------------------------------------------------- 
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    b.s <- s.gamma * b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma * b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)
    #print(pi.s)
    #----------------------------------------------------------------------------------------------------
    tmp1 <- as.vector(Sigma.inv%*%b.s)             # vector temporal para definir las ec.(21) pag.8.
    aux1 <- 1+sum(as.vector(b.s)*tmp1)
    sigma2.s <- sigma.xi^2/aux1     # ver pag.8
    sigma2.s1<- sigma.xi/aux1 
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=sigma2.s1*as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    tmp2     <- as.vector(Sigma.inv%*%b.t)             # vector temporal para definir las ec.(21) pag.8
    aux2 <- 1+sum(as.vector(b.t)*tmp2)
    sigma2.t <- sigma.xi^2/aux2
    sigma2.t1<- sigma.xi/aux2
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=sigma2.t1*as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    ##---------------------------------------------------------------------------------------------------
    R.s   <- mu.si/sqrt(sigma2.s)  
    pR.s  <- pnorm(R.s)                             # lo utilizamos varias veces pag.8.
    R.t   <- mu.ti/sqrt(sigma2.t)  
    pR.t  <- pnorm(R.t)
    pf.s  <- pi.s * dmnorm(z, as.vector(mu), Sigma.s) * pR.s
    pf.t  <- pi.t * dmnorm(z, as.vector(mu), Sigma.t) * pR.t
    pi.si <- pf.s/(pf.s + pf.t)                    # ver pag.8 
    pi.ti <- pf.t/(pf.s + pf.t)                    # ver pag.8 
    R.s.x <- pmax(R.s,-37)  # Clecio
    R.t.x <- pmax(R.t,-37)  # Clecio           
    dR.s  <- sqrt(sigma2.s)*dnorm(R.s.x)/pnorm(R.s.x)
    dR.t  <- sqrt(sigma2.t)*dnorm(R.t.x)/pnorm(R.t.x)   
    ##----------------------------------------------------------------------------------------------------
    #        calcularemos las esperanzas condicionales definidas en la página 8 ec.(24)-(28)
    ##----------------------------------------------------------------------------------------------------
    Es   <- 2*pi.si - 1                           # pag.8, ec.(24)
    Ev   <- (mu.si + dR.s)*pi.si + (mu.ti + dR.t)*pi.ti   # ec.(25)
    Ev2  <- (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si + (mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti     # ec.(26)
    Evw  <- (s.gamma * (mu.si + dR.s)*pi.si -t.gamma * (mu.ti + dR.t)*pi.ti)     # ec.(27)      
    Ev2w2 <- s.gamma^2 * (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si +t.gamma^2*(mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti# ec.(28)
#    Ev2  <- (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si + (mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti   # ec.(26)  # Clecio
#    Ev2w2 <- s.gamma^2 * (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si + t.gamma^2 * (mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti# ec.(28)    # CLecio          
    ##-------------------------------------------------------------------------------------------------
    #        Paso - M
    ##-------------------------------------------------------------------------------------------------
    #gamma<-mean(Es)
    #pi.s=(1+gamma)/2    
    sigma2.xi<-mean(Ev2)   
    sigma.xi <- sqrt(sigma2.xi)
    mu.x <- x.bar-mean(Evw)                      # ver pag.11
    # solucao Reinaldo e Clecio
    mu.y <- matrix(y.bar)-beta*mean(Evw)                # ver pag.11
    resR=y-matrix(1,n,1)%*%t(matrix(mu.y))
    beta  <- t(resR)%*%matrix(Evw)/sum(Ev2w2)
    re.x <- x-mu.x
    sigma2.u <- mean(re.x^2)-2*mean(Evw*re.x)+mean(Ev2w2)
    #Sigma.e <- (var(re.y)*(n-1)/n - 1/n*matrix(betav)%*%t(colSums(Evw*re.y)) - 1/n*colSums(Evw*re.y)%*%t(betav) +mean(Ev2w2)*outer(betav,betav)) #ver pag.11 , NAO ESTRUTURADA
    #sigmav <- t(resR)%*%resR - t(resR)%*%matrix(Evw)%*%t(beta)-beta%*%t(matrix(Evw))%*%resR +sum(Ev2w2)*beta%*%t(beta) 
    #Sigma.e <-sigmav/n # Modelo nao-estruturado
    Sigma.e<- apply(resR^2,2,mean)-2*as.vector(beta)*colSums(resR*Evw)/n+ (as.vector(beta))^2*mean(Ev2w2)   # modelo com Sigma.e matriz diagonal, vetor
    mu.xi <- mu.x                             # ver pag.11
    alpha <- mu.y - beta*mu.xi  # as.matrix(y.bar) - beta*mean(x)     
    Sigma <- diag(c(sigma2.u,as.vector(Sigma.e)))
    mu[1]=mu.xi
    mu[2:(d+1)] <- alpha + beta*mu.xi     # ver pag.5, ec(14)
    b.xi[1]=sigma.xi
    b.xi[2:(d+1)]=sigma.xi*beta
    Sigma.inv <- diag(c(1/sigma2.u,1/Sigma.e)) # caso diagonal
    gammai <-optimize(Q_gamma,c(-0.999,0.999),z,mu, b.xi,Sigma,Sigma.inv)
    gamma <- as.numeric(gammai$minimum)
    theta=as.vector(c(as.vector(alpha),as.vector(beta),as.vector(Sigma.e),mu.xi,sigma2.xi,sigma2.u,gamma))  
    print(theta)
    logi=Q_gamma(gamma,z,mu, b.xi,Sigma,Sigma.inv)
	print(logi)
    criterio <- abs(logi-log0)   #step.Q
    #criterio=sqrt(sum((theta-theta0)^2))
    theta0 <- theta
	log0<-logi
    #print(theta)
    
  }
  #message("--- EM done")
  #print(cont)
  #print(c(alpha,beta,mu.xi,sigma.xi,sigma2.u,gamma))
  #Ssigma.u <- sqrt(sigma2.u)
  list(alpha=alpha, beta=beta, mu.xi=mu.xi, sigma2.xi=sigma2.xi, sigma2.u=sigma2.u, 
       Sigma.e=Sigma.e, gamma=gamma)  
}



Q_gamma<-function(gamma,z,mu, b.xi,Sigma,Sigma.inv){
    n <- nrow(z)
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    b.s <- s.gamma*b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma*b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)    
    aux1 <- sqrt(1+as.numeric(t(b.s)%*%Sigma.inv%*%b.s))
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))/aux1
    aux2 <- sqrt(1+as.numeric(t(b.t)%*%Sigma.inv%*%b.t))
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))/aux2
    ##---------------------------------------------------------------------------------------------------
    pf.s  <- pi.s*dmnorm(z, as.vector(mu), Sigma.s)*pnorm(mu.si)
    pf.t  <- pi.t*dmnorm(z, as.vector(mu), Sigma.t)*pnorm(mu.ti)
    fz <- 2*(pf.s+pf.t)
    return(-sum(log(fz)))
}


logNTP_MEM<-function(theta,z){
    n <- nrow(z)
    alpha=theta$alpha
    d=length(as.vector(alpha))
    beta=theta$beta
    mu.xi=theta$mu.xi
    sigma2.xi=theta$sigma2.xi
    sigma.xi=sqrt(sigma2.xi)
    sigma2.u=theta$sigma2.u
    Sigma.e=as.vector(theta$Sigma.e)
    gamma=theta$gamma
    #Sigma=matrix(d+1,d+1)
    #Sigma[1,1]=sigma2.u
    #Sigma[2:(d+1),2:(d+1)]=Sigma.e 
    Sigma <-diag(c(sigma2.u,Sigma.e))
    #Sigma.inv   <- pd.solve(Sigma, log.det=TRUE)
    Sigma.inv  <-diag(c(sigma2.u,Sigma.e)) # caso diagonal
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    b.xi=matrix(0,d+1,1)
    b.xi[1]=sigma.xi
    b.xi[2:(d+1)]=sigma.xi*beta
    mu=matrix(0,d+1,1)
    mu[1]=mu.xi
    mu[2:(d+1)] <- alpha + beta*mu.xi
    b.s <- s.gamma*b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma*b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)    
    aux1 <- sqrt(1+as.numeric(t(b.s)%*%Sigma.inv%*%b.s))
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))/aux1
    aux2 <- sqrt(1+as.numeric(t(b.t)%*%Sigma.inv%*%b.t))
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))/aux2
    ##---------------------------------------------------------------------------------------------------
    pf.s  <- pi.s*dmnorm(z, as.vector(mu), Sigma.s)*pnorm(mu.si)
    pf.t  <- pi.t*dmnorm(z, as.vector(mu), Sigma.t)*pnorm(mu.ti)
    fz <- 2*(pf.s+pf.t)
    return(sum(log(fz)))
}

EMV_NTP_MEM<-function(z){ # muktivariado
  d <- ncol(z) - 1
  n <- nrow(z)
  x <- z[,1]
  y <- z[,2:(d+1)]
  z.bar <- mu <- apply(z, 2,mean)       # conjunto de valores iniciales pag11.
  mu.xi <- as.numeric(mu[1])       
  mu.y <- as.matrix(mu[2:(d+1)])
  S <- var(z)   # \Sigma                # valor unicial definido pag.11.
  #Sigma.e <- S[2:(d+1), 2:(d+1)]        # \Sigma esta definido como una partición pag.5.
  sigma.e <- apply(y,2,var)  # caso especial constante
  sigma.xi <- sqrt(S[1,1])   # ver pag.11.
  sigma2.xi=sigma.xi^2
  gamma <- mean((x-mu.xi)^3)/sigma.xi^3  # ver pag.11. 
  sigma2.u <- (sigma.xi)^2/100           # valor inicial arbitrario.  
  beta <- drop(S[1, 2:(d+1)]/S[1,1])    # ver pag.11  
  alpha <- mu.y - mu.xi*beta
  theta0=c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma.e,gamma)
  print(theta0)
  LI=c(-rep(Inf,2*d+1),rep(0.01,d+2),-0.999)
  LS=c(rep(Inf,3*d+4),0.999)
  theta=optim(theta0,logNTP_MEMmax,gr = NULL, z, method = "L-BFGS-B", lower = LI, upper = LS,control=list(parscale=theta0))
  return(theta)
}

logNTP_MEMmax<-function(theta,z){
    #c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,Sigmae,gamma)
    print(theta)
    d = ncol(z)-1
    alpha=theta[1:d]
    beta=theta[(d+1):2*d]
    mu.xi=as.numeric(theta[2*d+1])
    sigma2.xi=as.numeric(theta[2*d+2])
    sigma2.u=as.numeric(theta[2*d+3])
    Sigma.e=diag(as.vector(theta[(2*d+4):(3*d+3)]))
    gamma=as.numeric(theta[3*d+4])
    Sigma=diag(c(sigma2.u,as.vector(theta[(2*d+4):(3*d+3)])))
    Sigma.inv   <- pd.solve(Sigma, log.det=TRUE)
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    mu <- c(0,alpha) + c(1, beta)*mu.xi     # ver pag.5, ec(14)
    b.xi <- sqrt(sigma2.xi)*matrix(c(1, beta))           # ver pag.5, ec(14)
    b.s <- s.gamma*b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma*b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)    
    #Sigma.inv<- pd.solve(Sigma, log.det=TRUE)    # inversa de \Sigma, se usa paquete mnormt
    tmp1     <- drop(Sigma.inv%*%b.s)             # vector temporal para definir las ec.(21) pag.8.
    mu.si    <- drop(tmp1%*%t(z-mu))/sqrt(1+sum(b.s*tmp1))
    tmp2     <- drop(Sigma.inv%*%b.t)             # vector temporal para definir las ec.(21) pag.8
    mu.ti    <- drop(tmp2%*%t(z-mu))/sqrt(1+sum(b.t*tmp2))
    ##---------------------------------------------------------------------------------------------------
    pf.s  <- pi.s*dmnorm(z, mu, Sigma.s)*pnorm(mu.si)
    pf.t  <- pi.t*dmnorm(z, mu, Sigma.t)*pnorm(mu.ti)
    fz <- 2*(pf.s+pf.t)
    return(-sum(log(fz)))
}

MIapprox<-function(theta,z){
    #c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma)
    n <- nrow(z)
    q=ncol(z)
    d=q-1
    x <- z[,1]
    y <- z[,2:q]
    alpha=theta$alpha
    beta=theta$beta
    mu.xi=theta$mu.xi
    sigma2.xi=theta$sigma2.xi
    sigma.xi=sqrt(sigma2.xi)
    sigma2.u=theta$sigma2.u
    sigma2.e=theta$Sigma.e # estah vetor
    sigma2e.inv=diag(c(1/sigma2.e))
    gamma=theta$gamma
    Sigma=diag(c(sigma2.u,sigma2.e))
    Sigma.inv   <- diag(c(1/sigma2.u,1/sigma2.e))
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    mu <- c(0,alpha) + c(1, beta)*mu.xi     # ver pag.5, ec(14)
    mu <- matrix(mu,d+1,1)
    b.xi <- sqrt(sigma2.xi)*matrix(c(1, beta),d+1,1)           # ver pag.5, ec(14)
    b.s <- s.gamma*b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma*b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)    
    #----------------------------------------------------------------------------------------------------
    tmp1 <- as.vector(Sigma.inv%*%b.s)             # vector temporal para definir las ec.(21) pag.8.
    aux1 <- 1+sum(as.vector(b.s)*tmp1)
    sigma2.s <- sigma.xi^2/aux1     # ver pag.8
    sigma2.s1<- sigma.xi/aux1 
    #mu.si    <- drop(tmp1 %*% t(z-mu))*sigma2.s1    # ver pag.8. ec.(21)
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=sigma2.s1*as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    tmp2     <- as.vector(Sigma.inv%*%b.t)             # vector temporal para definir las ec.(21) pag.8
    aux2 <- 1+sum(as.vector(b.t)*tmp2)
    sigma2.t <- sigma.xi^2/aux2
    sigma2.t1<- sigma.xi/aux2
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=sigma2.t1*as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    ##---------------------------------------------------------------------------------------------------
    R.s   <- mu.si/sqrt(sigma2.s)  
    pR.s  <- pnorm(R.s)                             # lo utilizamos varias veces pag.8.
    R.t   <- mu.ti/sqrt(sigma2.t)  
    pR.t  <- pnorm(R.t)
    pf.s  <- pi.s * dmnorm(z, as.vector(mu), Sigma.s) * pR.s
    pf.t  <- pi.t * dmnorm(z, as.vector(mu), Sigma.t) * pR.t
    pi.si <- pf.s/(pf.s + pf.t)                    # ver pag.8 
    pi.ti <- pf.t/(pf.s + pf.t)                    # ver pag.8 
    R.s.x <- pmax(R.s,-37)  # Clecio
    R.t.x <- pmax(R.t,-37)  # Clecio           
    dR.s  <- sqrt(sigma2.s)*dnorm(R.s.x)/pnorm(R.s.x)
    dR.t  <- sqrt(sigma2.t)*dnorm(R.t.x)/pnorm(R.t.x)   
    ##----------------------------------------------------------------------------------------------------
    Es   <- pi.si - pi.ti                                 # pag.8, ec.(24)
    Ev   <- (mu.si + dR.s)*pi.si + (mu.ti + dR.t)*pi.ti   # ec.(25)
    Ev2  <- (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si + (mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti     # ec.(26)
    Evw  <- s.gamma * (mu.si + dR.s)*pi.si -t.gamma * (mu.ti + dR.t)*pi.ti     # ec.(27)      
    Ev2w2 <- s.gamma^2*(mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si +t.gamma^2*(mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti# ec.(28)
    ##-------------------------------------------------------------------------------------------------
    MI=matrix(0,3*d+4,3*d+4)
    si=matrix(0,3*d+4,1)
    #c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma)
    for (i in 1:n){
	si[1:d]=as.vector(sigma2e.inv%*%matrix(y[i]-alpha-(Evw[i]+mu.xi)*beta)) # alpha
	si[(d+1):(2*d)]=as.vector(sigma2e.inv%*%matrix((Evw[i]+mu.xi)*(y[i]-alpha)-(mu.xi^2+2*mu.xi*Evw[i]+Ev2w2[i])*beta))  # beta
	si[2*d+1]=(x[i]-mu.xi-Evw[i])/sigma2.u - mu.xi*t(matrix(beta))%*%sigma2e.inv%*%matrix(beta)+t(matrix(beta))%*%sigma2e.inv%*%matrix(y[i]-alpha-beta*Evw[i])     # mu.xi
	si[2*d+2]=-0.5/sigma2.xi+0.5/(sigma2.xi^2)*Ev2[i]			# sigma2.xi
	si[2*d+3]=-0.5/sigma2.u+0.5/(sigma2.u^2)*((x[i]-mu.xi)^2-2*Evw[i]*(x[i]-mu.xi)+Ev2w2[i])		# sigma2.u
	si[(2*d+4):(3*d+3)]=-0.5/sigma2.e+0.5/(sigma2.e^2)*((y[i]-alpha)^2-2*(Evw[i]+mu.xi)*beta*(y[i]-alpha)+(mu.xi^2+2*mu.xi*Evw[i]+Ev2w2[i])*(beta^2)) # sigma2.e
	si[3*d+4]=0.5*((1+Es[i])/(1+gamma)-(1-Es[i])/(1-gamma))		# gamma
	aux=si%*%t(si)
	MI=MI+aux
}
    return(MI)
}


MIapprox2<-function(theta,z){
    #c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma)
    n <- nrow(z)
    q=ncol(z)
    d=q-1
    x <- z[,1]
    y <- z[,2:q]
    alpha=theta$alpha
    beta=theta$beta
    mu.xi=theta$mu.xi
    sigma2.xi=theta$sigma2.xi
    sigma.xi=sqrt(sigma2.xi)
    sigma2.u=theta$sigma2.u
    sigma2.e=theta$Sigma.e # estah vetor
    sigma2e.inv=diag(c(1/sigma2.e))
    gamma=theta$gamma
    Sigma=diag(c(sigma2.u,sigma2.e))
    Sigma.inv   <- diag(c(1/sigma2.u,1/sigma2.e))
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    mu <- c(0,alpha) + c(1, beta)*mu.xi     # ver pag.5, ec(14)
    mu <- matrix(mu,d+1,1)
    b.xi <- sqrt(sigma2.xi)*matrix(c(1, beta),d+1,1)           # ver pag.5, ec(14)
    b.s <- s.gamma*b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma*b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)    
    #----------------------------------------------------------------------------------------------------
    tmp1 <- as.vector(Sigma.inv%*%b.s)             # vector temporal para definir las ec.(21) pag.8.
    aux1 <- 1+sum(as.vector(b.s)*tmp1)
    sigma2.s <- sigma.xi^2/aux1     # ver pag.8
    sigma2.s1<- sigma.xi/aux1 
    #mu.si    <- drop(tmp1 %*% t(z-mu))*sigma2.s1    # ver pag.8. ec.(21)
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=sigma2.s1*as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    tmp2     <- as.vector(Sigma.inv%*%b.t)             # vector temporal para definir las ec.(21) pag.8
    aux2 <- 1+sum(as.vector(b.t)*tmp2)
    sigma2.t <- sigma.xi^2/aux2
    sigma2.t1<- sigma.xi/aux2
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=sigma2.t1*as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    ##---------------------------------------------------------------------------------------------------
    R.s   <- mu.si/sqrt(sigma2.s)  
    pR.s  <- pnorm(R.s)                             # lo utilizamos varias veces pag.8.
    R.t   <- mu.ti/sqrt(sigma2.t)  
    pR.t  <- pnorm(R.t)
    pf.s  <- pi.s * dmnorm(z, as.vector(mu), Sigma.s) * pR.s
    pf.t  <- pi.t * dmnorm(z, as.vector(mu), Sigma.t) * pR.t
    pi.si <- pf.s/(pf.s + pf.t)                    # ver pag.8 
    pi.ti <- pf.t/(pf.s + pf.t)                    # ver pag.8 
    R.s.x <- pmax(R.s,-37)  # Clecio
    R.t.x <- pmax(R.t,-37)  # Clecio           
    dR.s  <- sqrt(sigma2.s)*dnorm(R.s.x)/pnorm(R.s.x)
    dR.t  <- sqrt(sigma2.t)*dnorm(R.t.x)/pnorm(R.t.x)   
    ##----------------------------------------------------------------------------------------------------
    Es   <- pi.si - pi.ti                                 # pag.8, ec.(24)
    Ev   <- (mu.si + dR.s)*pi.si + (mu.ti + dR.t)*pi.ti   # ec.(25)
    Ev2  <- (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si + (mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti     # ec.(26)
    Evw  <- s.gamma * (mu.si + dR.s)*pi.si -t.gamma * (mu.ti + dR.t)*pi.ti     # ec.(27)      
    Ev2w2 <- s.gamma^2*(mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si +t.gamma^2*(mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti# ec.(28)
    ##-------------------------------------------------------------------------------------------------
    MI=matrix(0,3*d+4,3*d+4)
    si=matrix(0,3*d+4,1)
    mu.y=alpha+beta*mu.xi	
    #c(alpha,beta,sigma2.e,mu.xi,sigma2.xi,sigma2.u,gamma)
    for (i in 1:n){
	si[1:d]=as.vector(sigma2e.inv%*%matrix(y[i,]-mu.y-Evw[i]*beta)) # alpha
	#si[(d+1):(2*d)]=as.vector(sigma2e.inv%*%matrix(mu.xi*(y[i,]-mu.y)+Evw[i]*(y[i,]-mu.y-mu.xi*beta)+Ev2w2[i]*beta))  # beta
      si[(d+1):(2*d)]=as.vector(sigma2e.inv%*%matrix((mu.xi+Evw[i])*(y[i,]-mu.y)-(Evw[i]*mu.xi+Ev2w2[i])*beta))  # beta
      si[2*d+1]=(x[i]-mu.xi-Evw[i])/sigma2.u - mu.xi*t(matrix(beta))%*%sigma2e.inv%*%matrix(beta)+t(matrix(beta))%*%sigma2e.inv%*%matrix(y[i,]-alpha-beta*Evw[i]) # mu.xi
	si[2*d+2]=-0.5/sigma2.xi+0.5/(sigma2.xi^2)*Ev2[i]			# sigma2.xi
	si[2*d+3]=-0.5/sigma2.u+0.5/(sigma2.u^2)*((x[i]-mu.xi)^2-2*Evw[i]*(x[i]-mu.xi)+Ev2w2[i])		# sigma2.u
      si[(2*d+4):(3*d+3)]=-0.5/sigma2.e+0.5/(sigma2.e^2)*((y[i,]-mu.y)^2-2*Evw[i]*beta*(y[i,]-mu.y)+Ev2w2[i]*(beta^2)) # sigma2.e
	si[3*d+4]=0.5*((1+Es[i])/(1+gamma)-(1-Es[i])/(1-gamma))		# gamma
	aux=si%*%t(si)
	MI=MI+aux
}
    return(MI)
}



#############################################################################################

EM.NTP_MEMuniv <- function(z) # ajuste NTP-MEM por algoritmo EM
{ # asumimos z es (n=nº de filas, (d+1)=nº de columnas), z_i =(x_i,y_i\T)\T, pag5.
  d <- ncol(z) - 1
  n <- nrow(z)
  x <- z[,1]
  y <- z[,2]
  z.bar <- mu <- apply(z, 2,mean)       # conjunto de valores iniciales pag11.
  x.bar <- mu.xi <- mu.x <- as.numeric(mu[1])       
  y.bar <- mu.y <- mu[2]
  S=var(z)
  sigma2.e <- var(y)
  sigma.x <- sigma.xi <- sd(x)   # ver pag.11.
  sigma2.xi=sigma.xi^2
  gamma <-0# mean((x-x.bar)^3)/sigma.x^3  # ver pag.11. 
  sigma2.u <- (sigma.x)^2/100           # valor inicial arbitrario.  
  Sigma <- diag(c(sigma2.u, sigma2.e))    # ver pag.5.
  beta <- as.numeric(drop(S[1,2]/S[1,1]))    # ver pag.11
  alpha <- y.bar - x.bar*beta
  mu <- matrix(c(0,alpha) + c(1, beta)*mu.xi,2,1)     # ver pag.5, ec(14)
  b.xi <- sigma.xi*matrix(c(1, beta),2,1)           # ver pag.5, ec(14)
  Sigma.inv   <- diag(c(1/sigma2.u,1/sigma2.e))
  log0=-Q_gamma(gamma,z,mu, b.xi,Sigma,Sigma.inv)
  #print(log0)
  #print(c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma))
  criterio=1
  cont=0
  while (criterio>1e-4) {
  cont=cont+1
    #cat("------- k: " , k, "\n")        #cantidad de iteraciones que hace el prgrama
    #------------------------------ paso - E ------------------------------------------------------------- 
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    b.s <- s.gamma * b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma * b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)
    #print(pi.s)
    #----------------------------------------------------------------------------------------------------
    tmp1 <- as.vector(Sigma.inv%*%b.s)             # vector temporal para definir las ec.(21) pag.8.
    aux1 <- 1+sum(as.vector(b.s)*tmp1)
    sigma2.s <- sigma.xi^2/aux1     # ver pag.8
    sigma2.s1<- sigma.xi/aux1 
    #mu.si    <- drop(tmp1 %*% t(z-mu))*sigma2.s1    # ver pag.8. ec.(21)
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=sigma2.s1*as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    tmp2     <- as.vector(Sigma.inv%*%b.t)             # vector temporal para definir las ec.(21) pag.8
    aux2 <- 1+sum(as.vector(b.t)*tmp2)
    sigma2.t <- sigma.xi^2/aux2
    sigma2.t1<- sigma.xi/aux2
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=sigma2.t1*as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    
    ##---------------------------------------------------------------------------------------------------
    R.s   <- mu.si/sqrt(sigma2.s)  
    pR.s  <- pnorm(R.s)                             # lo utilizamos varias veces pag.8.
    R.t   <- mu.ti/sqrt(sigma2.t)  
    pR.t  <- pnorm(R.t)
    pf.s  <- pi.s * dmnorm(z, as.vector(mu), Sigma.s) * pR.s
    pf.t  <- pi.t * dmnorm(z, as.vector(mu), Sigma.t) * pR.t
    pi.si <- pf.s/(pf.s + pf.t)                    # ver pag.8 
    pi.ti <- pf.t/(pf.s + pf.t)                    # ver pag.8 
    R.s.x <- pmax(R.s,-37)  # Clecio
    R.t.x <- pmax(R.t,-37)  # Clecio           
    dR.s  <- sqrt(sigma2.s)*dnorm(R.s.x)/pnorm(R.s.x)
    dR.t  <- sqrt(sigma2.t)*dnorm(R.t.x)/pnorm(R.t.x)   
    ##----------------------------------------------------------------------------------------------------
    #        calcularemos las esperanzas condicionales definidas en la página 8 ec.(24)-(28)
    ##----------------------------------------------------------------------------------------------------
    Es   <- pi.si - pi.ti                                 # pag.8, ec.(24)
    Ev   <- (mu.si + dR.s)*pi.si + (mu.ti + dR.t)*pi.ti   # ec.(25)
    Ev2  <- (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si + (mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti     # ec.(26)
    Evw  <- s.gamma * (mu.si + dR.s)*pi.si -t.gamma * (mu.ti + dR.t)*pi.ti     # ec.(27)      
    Ev2w2 <- s.gamma^2*(mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si +t.gamma^2*(mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti# ec.(28)
    ##-------------------------------------------------------------------------------------------------
    #        Paso - M
    ##-------------------------------------------------------------------------------------------------
    sigma2.xi<-mean(Ev2)   
    sigma.xi <- sqrt(sigma2.xi)
    mu.x <- x.bar-mean(Evw) 
    re.x <-x-mu.x 
    mu.y <- y.bar-beta*mean(Evw)                # ver pag.11
    re.y=y-mu.y
    beta  <- sum(re.y*Evw)/sum(Ev2w2) 
    sigma2.u <- mean(re.x^2)-2*mean(Evw*re.x)+mean(Ev2w2)
    sigma2.e <- mean(re.y^2) - 2*mean(Evw*re.y)*beta +mean(Ev2w2)*(beta^2) #ver pag.11
    mu.xi <- mu.x                             # ver pag.11
    alpha <- y.bar - beta*x.bar   
    mu <- matrix(c(0,alpha) + c(1, beta)*mu.xi,2,1)     # ver pag.5, ec(14)
    b.xi <- sigma.xi*matrix(c(1, beta),2,1)           # ver pag.5, ec(14)
    Sigma=diag(c(sigma2.u,sigma2.e))
    Sigma.inv   <- diag(c(1/sigma2.u,1/sigma2.e))
    #gammai <-optimize(Q_gamma,c(-0.999,0.999),z,mu, b.xi,Sigma,Sigma.inv)
    #gamma <- as.numeric(gammai$minimum)
    #gamma <- mean(Es)
    logi=-Q_gamma(gamma,z,mu, b.xi,Sigma,Sigma.inv)
    criterio <- abs(logi-log0) 
    log0 <- logi
    #print(logi)
    #print(c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma))
  }
  message("--- EM done")
  print(cont)
  print(logi)
  #print(c(alpha,beta,mu.xi,sigma.xi,sigma2.u,gamma))
  #Ssigma.u <- sqrt(sigma2.u)
  list(alpha=alpha, beta=beta, mu.xi=mu.xi, sigma2.xi=sigma2.xi, sigma2.u=sigma2.u, 
       sigma2.e=sigma2.e, gamma=gamma)  
}



logNTP_MEMuniv<-function(theta,z){
    n=nrow(z)
    d=1
    alpha=theta$alpha
    beta=theta$beta
    mu.xi=theta$mu.xi
    sigma2.xi=theta$sigma2.xi
    sigma2.u=theta$sigma2.u
    sigma2.e=theta$sigma2.e
    gamma=theta$gamma
    Sigma=diag(c(sigma2.u,sigma2.e))
    Sigma.inv   <- diag(c(1/sigma2.u,1/sigma2.e))
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    mu <- c(0,alpha) + c(1, beta)*mu.xi     # ver pag.5, ec(14)
    b.xi <- sqrt(sigma2.xi)*matrix(c(1, beta))           # ver pag.5, ec(14)
    b.s <- s.gamma*b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma*b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)    
    #----------------------------------------------------------------------------------------------------
    aux1 <- sqrt(1+as.numeric(t(b.s)%*%Sigma.inv%*%b.s))
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))/aux1
    aux2 <- sqrt(1+as.numeric(t(b.t)%*%Sigma.inv%*%b.t))
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))/aux2
    ##---------------------------------------------------------------------------------------------------
    pf.s  <- pi.s*dmnorm(z, mu, Sigma.s)*pnorm(mu.si)
    pf.t  <- pi.t*dmnorm(z, mu, Sigma.t)*pnorm(mu.ti)
    fz <- 2*(pf.s+pf.t)
    return(sum(log(fz)))
}

EMV_NTP_MEMuniv<-function(z){ # univariado
  d=1
  n <- nrow(z)
  x <- z[,1]
  y <- z[,2]
  z.bar <- mu <- apply(z, 2,mean)       # conjunto de valores iniciales pag11.
  mu.xi <- as.numeric(mu[1])       
  mu.y <- as.numeric(mu[2])
  S <- var(z)   # \Sigma                # valor unicial definido pag.11.
  #Sigma.e <- S[2:(d+1), 2:(d+1)]        # \Sigma esta definido como una partición pag.5.
  sigma2.e <- var(y)  # caso especial constante
  sigma.xi <- sqrt(S[1,1])   # ver pag.11.
  sigma2.xi=sigma.xi^2
  gamma <- mean((x-mu.xi)^3)/sigma.xi^3  # ver pag.11. 
  sigma2.u <- (sigma.xi)^2/100           # valor inicial arbitrario.  
  beta <- drop(S[1,2]/S[1,1])    # ver pag.11
  alpha <- mu.y - mu.xi*beta
  theta0=c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma)
  print(theta0)
  LI=c(-rep(Inf,2*d+1),rep(0.01,d+2),-0.999)
  LS=c(rep(Inf,3*d+4),0.999)
  theta=optim(theta0,logNTP_MEMmax,gr = NULL, z, method = "L-BFGS-B", lower = LI, upper = LS,control=list(parscale=theta0))
  return(theta)
}


MIapprox_univ<-function(theta,z){
    #c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma)
    n <- nrow(z)
    x <- z[,1]
    y <- z[,2]
    alpha=theta$alpha
    beta=theta$beta
    mu.xi=theta$mu.xi
    sigma2.xi=theta$sigma2.xi
    sigma2.u=theta$sigma2.u
    sigma2.e=theta$sigma2.e
    gamma=theta$gamma
    Sigma=diag(c(sigma2.u,sigma2.e))
    Sigma.inv   <- diag(c(1/sigma2.u,1/sigma2.e))
    s.gamma <- fn.s(gamma)          
    t.gamma <- fn.t(gamma)
    mu <- c(0,alpha) + c(1, beta)*mu.xi     # ver pag.5, ec(14)
    b.xi <- sqrt(sigma2.xi)*matrix(c(1, beta))           # ver pag.5, ec(14)
    b.s <- s.gamma*b.xi                   # ver pag.5, ec(14)
    b.t <- -t.gamma*b.xi
    Sigma.s <- Sigma + b.s%*%t(b.s)      # ver pag.5, ec(14)
    Sigma.t <- Sigma + b.t%*%t(b.t)
    pi.s <- s.gamma/(s.gamma+t.gamma)       # ver pag.5, ec(14)
    pi.t <- t.gamma/(s.gamma+t.gamma)    
    #----------------------------------------------------------------------------------------------------
    tmp1 <- as.vector(Sigma.inv%*%b.s)             # vector temporal para definir las ec.(21) pag.8.
    aux1 <- 1+sum(as.vector(b.s)*tmp1)
    sigma2.s <- sigma.xi^2/aux1     # ver pag.8
    sigma2.s1<- sigma.xi/aux1 
    #mu.si    <- drop(tmp1 %*% t(z-mu))*sigma2.s1    # ver pag.8. ec.(21)
    mu.si=rep(0,n)
    for (i in 1:n) mu.si[i]=sigma2.s1*as.numeric(t(b.s)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    tmp2     <- as.vector(Sigma.inv%*%b.t)             # vector temporal para definir las ec.(21) pag.8
    aux2 <- 1+sum(as.vector(b.t)*tmp2)
    sigma2.t <- sigma.xi^2/aux2
    sigma2.t1<- sigma.xi/aux2
    mu.ti=rep(0,n)
    for (i in 1:n) mu.ti[i]=sigma2.t1*as.numeric(t(b.t)%*%Sigma.inv%*%(matrix(z[i,])-mu))
    ##---------------------------------------------------------------------------------------------------
    R.s   <- mu.si/sqrt(sigma2.s)  
    pR.s  <- pnorm(R.s)                             # lo utilizamos varias veces pag.8.
    R.t   <- mu.ti/sqrt(sigma2.t)  
    pR.t  <- pnorm(R.t)
    pf.s  <- pi.s * dmnorm(z, as.vector(mu), Sigma.s) * pR.s
    pf.t  <- pi.t * dmnorm(z, as.vector(mu), Sigma.t) * pR.t
    pi.si <- pf.s/(pf.s + pf.t)                    # ver pag.8 
    pi.ti <- pf.t/(pf.s + pf.t)                    # ver pag.8 
    R.s.x <- pmax(R.s,-37)  # Clecio
    R.t.x <- pmax(R.t,-37)  # Clecio           
    dR.s  <- sqrt(sigma2.s)*dnorm(R.s.x)/pnorm(R.s.x)
    dR.t  <- sqrt(sigma2.t)*dnorm(R.t.x)/pnorm(R.t.x)   
    ##----------------------------------------------------------------------------------------------------
    Es   <- pi.si - pi.ti                                 # pag.8, ec.(24)
    Ev   <- (mu.si + dR.s)*pi.si + (mu.ti + dR.t)*pi.ti   # ec.(25)
    Ev2  <- (mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si + (mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti     # ec.(26)
    Evw  <- s.gamma * (mu.si + dR.s)*pi.si -t.gamma * (mu.ti + dR.t)*pi.ti     # ec.(27)      
    Ev2w2 <- s.gamma^2*(mu.si^2 + sigma2.s + mu.si*dR.s)*pi.si +t.gamma^2*(mu.ti^2 + sigma2.t + mu.ti*dR.t)*pi.ti# ec.(28)
    ##-------------------------------------------------------------------------------------------------
    MI=matrix(0,7,7)
    si=matrix(0,7,1)
    #c(alpha,beta,mu.xi,sigma2.xi,sigma2.u,sigma2.e,gamma)
    for (i in 1:n){
	si[1]=(y[i]-alpha-(Evw[i]+mu.xi)*beta)/sigma2.e # alpha
	si[2]=((Evw[i]+mu.xi)*(y[i]-alpha)-(mu.xi^2+2*mu.xi*Evw[i]+Ev2w2[i])*beta)/sigma2.e  # beta
	si[3]=(x[i]-mu.xi-Evw[i])/sigma2.u - mu.xi*(beta^2)/sigma2.e+beta*(y[i]-alpha-beta*Evw[i])/sigma2.e      # mu.xi
	si[4]=-0.5/sigma2.xi+0.5/(sigma2.xi^2)*Ev2[i]			# sigma2.xi
	si[5]=-0.5/sigma2.u+0.5/(sigma2.u^2)*((x[i]-mu.xi)^2-2*Evw[i]*(x[i]-mu.xi)+Ev2w2[i])		# sigma2.u
	si[6]=-0.5/sigma2.e+0.5/(sigma2.e^2)*((y[i]-alpha)^2-2*(Evw[i]+mu.xi)*beta*(y[i]-alpha)+(mu.xi^2+2*mu.xi*Evw[i]+Ev2w2[i])*(beta^2))		# sigma2.e
	si[7]=0.5*((1+Es[i])/(1+gamma)-(1-Es[i])/(1-gamma))		# gamma
	aux=si%*%t(si)
	MI=MI+aux
}
    return(MI)
}













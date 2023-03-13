###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o para gerar dados com fra??o de cura e depend?ncia entre tempos de falha e censura informativa.                                ###
### Tempos de falha seguem distribui??o Weibull e tempos de censura informativa tamb?m seguem distribui??o Weibull.                      ###
###--------------------------------------------------------------------------------------------------------------------------------------###
#simula.teoria_unif <- function(alpha.c, alpha.t,  lambda.c, lambda.t, alphaa, w,  X_cure,X_C, beta_cure,beta_C, ident)
#{
#  n <- nrow(X_cure);  t <- rep(NA,n);  c <- rep(NA,n);  e <- rep(NA,n);  y <- rep(NA,n);  u <- runif(n);  v <- runif(n)  
#  cens <- rep(NA,n)
#  
#  gamma.t <- lambda.t
#  gamma.c <- lambda.c*exp(X_C%*%beta_C + alphaa*w)  
#  s <- runif(n,0,10) # administrative censoring 
#  
#  Theta <- exp(X_cure%*%beta_cure + w)
#  
#  for(i in 1:n)
#  {
#    if ( u[i] <= (exp(-Theta[i]))){
#      t[i] <- Inf
#    }
#    if ( u[i] > (exp(-Theta[i]))){
#      t[i] <- (-log((log(u[i])/Theta[i])+1)/(gamma.t))^(1/alpha.t) 
#    }
#    c[i] <- (-log(v[i])/(gamma.c[i]))^(1/alpha.c)     
#    y[i] <- min(t[i], c[i],s[i])
#    e[i] <- as.numeric(t[i] < min(c[i],s[i]))
#    if(t[i] == min(t[i], c[i], s[i]))
#    {
#      cens[i] <- 1
#    }
#    if(c[i] == min(t[i], c[i], s[i]))
#    {
#      cens[i] <- 2
#    }
#    if(s[i] == min(t[i], c[i], s[i]))
#    {
#      cens[i] <- 3
#    }
#    
#  }  
#  data <- as.data.frame(cbind(u, v, t, c, y, e,X_cure,X_C, cens,ident))
#  names(data) <- c("u","v","t","c","time","event", "int","x1_cure","x_c1","x_c2", "cens", "ident")  
#  return(data)
#}
###--------------------------------------------------------------------------------------------------------------------------------------###
### funcao para gerar da log condicional completa da fragilidade, wk (sendo k o indice de grupo)                                         ###
###--------------------------------------------------------------------------------------------------------------------------------------###
#tambem muda pois acrescent frac cura
log_cond_wi <- function(w_k_aux,risco_a_T, risco_a_C, delta_T,delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta){

  w_k <- log(w_k_aux/(1-w_k_aux)) 
  
  if (ncol(t(t(X_cure))) == 1){
    pred_cure <- exp((t(t(X_cure))*beta_cure)+w_k)
  }
  if (ncol(t(t(X_cure))) != 1){
    pred_cure <- exp((X_cure%*%beta_cure)+w_k)
  }
  if (ncol(t(t(X_C))) == 1){
    exp_C <- exp((t(t(X_C))*beta_C)+alpha*w_k)
  }
  if (ncol(t(t(X_C))) != 1){
    exp_C <- exp((X_C%*%beta_C)+alpha*w_k)
  }
  
  log_vero_w <- sum(delta_T*(w_k) - pred_cure*(1- exp(-risco_a_T)) + delta_C*alpha*w_k - risco_a_C*exp_C)-((w_k^2)/(2*theta)) - log(w_k_aux*(1-w_k_aux))   
  
  return(log_vero_w)
}

###-----------------------------------------------------------------------------------------------------------------------------------------
# funcao usada no arms para a log condicional completa da fragilidade 
support_wi <-  function(w_k_aux,risco_a_T, risco_a_C, delta_T,delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta){(w_k_aux>0)*(w_k_aux<1)}  
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o log-verossimilhan?a, parte da falha, usada no optim para estimar os par?metros da distribui??o Weibull                        ### 
### e da fra??o de cura, considerando covari?veis e fragilidade na fra??o de cura                                                        ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_param_T <-  function(par,t,delta_T, X_cure, bi){
  p <- ncol(X_cure)
  n <- nrow(X_cure)
  
  if(ncol(t(t(X_cure))) == 1){
    pred <- X_cure*par[3]
    pred_cure <- exp(pred)*rowMeans(exp(bi))
  }
  
  else{
    pred <- X_cure%*%par[3:(p+2)]
    pred_cure <- exp(pred)*rowMeans(exp(bi))
  }
  
  a <- par[1]
  b <- exp(par[2])
  exp_T <- (t^a)*b
  
  U <- -sum(delta_T*(pred + rowMeans(bi) + log(a) + (a-1)*log(t) + log(b) - exp_T)- pred_cure*(1-exp(-exp_T)))
}

# mudam (s? entra frac de cura nois tempos de falha)
#modelo_param_T <-  function(par,t,delta_T, X_cure, bi){
#  
#  pred <- X_cure%*%par[3:4]
#  pred_cure <- exp(pred)*rowMeans(exp(bi))
#  a <- par[1]
#  b <- exp(par[2])
#  exp_T <- (t^a)*b
#  
#  U <- -sum(delta_T*(pred + rowMeans(bi) + log(a) + (a-1)*log(t) + log(b) - exp_T)- pred_cure*(1-exp(-exp_T)))
#}    
###--------------------------------------------------------------------------------------------------------------------------------------###
###------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o log-verossimilhan?a, parte da censura informativa, usada no optim para estimar os par?metros da distribui??o Weibull        ###
### considerando covari?veis na fun??o taxa de falha de C                                                                              ###
###------------------------------------------------------------------------------------------------------------------------------------###
#Essa ? igual ao depcens
modelo_param_C <-  function(par,t,delta_C, X_C, bi){
  
  q <- ncol(X_C)
  n <- nrow(X_C)
  
  if(ncol(t(t(X_C))) == 1){
    pred_C <- X_C*par[4]
  }
  
  else{
    pred_C <- X_C%*%par[4:(q+3)]
  }
  a <- par[1]
  b <- exp(par[2])
  exp_C <- (t^a)*b*exp(pred_C)
  
  U <- -sum(delta_C*(log(a) + (a-1)*log(t) + log(b) + pred_C + par[3]*rowMeans(bi))- exp_C*rowMeans(exp(par[3]*bi)))  
}
###--------------------------------------------------------------------------------------------------------------------------------------###
###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??es com as derivadas de primeira ordem    
###--------------------------------------------------------------------------------------------------------------------------------------###
#aqui tbm muda primeiro parametros de falha
Esp.Deriv.Prim.Ordem <-  function(t,delta_T, delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta, alpha_T, alpha_C, lambda_T, lambda_C, w_k_grupo, ident){
  
  p <- ncol(X_cure)
  q <- ncol(X_C)
  
  wk = w_k_grupo[ident,]   
  m <- nrow(w_k_grupo)
  
  num_param <- length(beta_cure)+2+length(beta_C)+2+length(alpha)+length(theta)
  deriv1 <- matrix(NA,num_param,ncol(wk)) 
  #deriv1 <- matrix(NA,10,ncol(wk))  
  if (ncol(t(t(X_cure))) == 1){
    pred_cure <- as.vector(exp(t(t(X_cure))*beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_cure))) != 1){
    pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_C))) == 1){
    pred_C <- as.vector(exp(t(t(X_C))*beta_C))*exp(alpha*wk)
  }
  if (ncol(t(t(X_C))) != 1){
    pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  }
  
  #pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  #pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  risco_a_T <- (t^alpha_T)*lambda_T
  risco_a_C <- (t^alpha_C)*lambda_C
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  
  for (i in 1:p){
    deriv1[i,] <- colSums(X_cure[,i]*(Delta_t - pred_cure*(1-exp(-risco_a_T)))) #derivada do vetor beta_T
  }
  
  #derivada de beta_cure_0
  #deriv1[1,] <- colSums(X_cure[,1]*(Delta_t - pred_cure*(1-exp(-risco_a_T)))) 
  
  #derivada de beta_cure_1
  #deriv1[2,] <- colSums(X_cure[,2]*(Delta_t - pred_cure*(1-exp(-risco_a_T)))) 
  
  #derivada de alpha.t
  deriv1[1+p,] <- colSums(Delta_t*(alpha_T^(-1) + log(t) - risco_a_T*log(t))- pred_cure*exp(-risco_a_T)*risco_a_T*log(t))
  
  #derivada de lambda.t
  deriv1[2+p,] <- colSums(Delta_t*(lambda_T^(-1) - t^alpha_T) - pred_cure*exp(-risco_a_T)*(t^alpha_T))
  
  
  for (i in 1:q){
    deriv1[2+p+i,] <- colSums(X_C[,i]*(Delta_c - risco_a_C*pred_C))  #derivada do vetor beta_C
  }
  
  #derivada de beta_C_1
  #deriv1[5,] <- colSums(X_C[,1]*(Delta_c - risco_a_C*pred_C)) 
  
  #derivada de beta_C_2
  #deriv1[6,] <- colSums(X_C[,2]*(Delta_c - risco_a_C*pred_C)) 
  
  #derivada de alpha
  deriv1[2+p+q+1,] <- colSums(wk*(Delta_c - risco_a_C*pred_C)) 
  
  #derivada de alpha.c
  deriv1[2+p+q+2,] <- colSums(Delta_c*(alpha_C^(-1) + log(t)) - risco_a_C*log(t)*pred_C)
  
  #derivada de lambda.c
  deriv1[2+p+q+3,] <- colSums(Delta_c*(lambda_C^(-1)) - (t^alpha_C)*pred_C)
  
  
  #derivada de theta
  deriv1[2+p+q+4,] <- -0.5*(m*theta^(-1) - (theta^(-2))*colSums(w_k_grupo^(2)))
  
  aux <- deriv1[,1]%*%t(deriv1[,1])
  
  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }
  
  return(aux/ncol(wk))
}
###--------------------------------------------------------------------------------------------------------------------------------------###



###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o com as derivadas parciais de segunda ordem                                                                                    ###
###--------------------------------------------------------------------------------------------------------------------------------------###
#idem, segue 1
Esp.DerivParciais <-  function(t,delta_T, delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta, alpha_T, alpha_C, lambda_T, lambda_C, w_k_grupo, ident){

  p <- ncol(X_cure)
  q <- ncol(X_C)
  
  num_param <- length(beta_cure)+2+length(beta_C)+2+length(alpha)+length(theta)
  deriv2 <- matrix(0,num_param,num_param)
  
  wk = w_k_grupo[ident,]   
  m <- nrow(w_k_grupo)
  if (ncol(t(t(X_cure))) == 1){
    pred_cure <- as.vector(exp(t(t(X_cure))*beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_cure))) != 1){
    pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  }
  if (ncol(t(t(X_C))) == 1){
    pred_C <- as.vector(exp(t(t(X_C))*beta_C))*exp(alpha*wk)
  }
  if (ncol(t(t(X_C))) != 1){
    pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  }
  risco_a_T <- (t^alpha_T)*lambda_T
  risco_a_C <- (t^alpha_C)*lambda_C
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  
  # ============================== #
  
  for (i in 1:p){
    for (u in 1:p) {
    deriv2[u,i] <- mean(colSums( - X_cure[,u]*X_cure[,i]*pred_cure*(1-exp(-risco_a_T))))   #derivada do vetor beta_T da derivada de alpha.t
    deriv2[i,u] <- deriv2[u,i]
    }
  }
  
  #derivada de beta_cura0 da derivada de beta_cura0
  #deriv2[1,1] <- mean(colSums( - X_cure[,1]*X_cure[,1]*pred_cure*(1-exp(-risco_a_T)))) 
  #
  ##derivada de beta_cura1 da derivada de beta_cura1
  #deriv2[2,2] <- mean(colSums( - X_cure[,2]*X_cure[,2]*pred_cure*(1-exp(-risco_a_T)))) 
  #
  ##derivada de beta_cura1 da derivada de beta_cura0
  #deriv2[1,2] <- mean(colSums( - X_cure[,1]*X_cure[,2]*pred_cure*(1-exp(-risco_a_T)))) 
  #deriv2[2,1] <- deriv2[1,2]
  
  # ============================== #
  
  for (i in 1:p){
      deriv2[i,p+1] <- mean(colSums( - X_cure[,i]*pred_cure*exp(-risco_a_T)*log(t)*risco_a_T))   #derivada de alpha.t da derivada de beta_cura0
      deriv2[p+1,i] <- deriv2[i,p+1]
  }
  
  #derivada de alpha.t da derivada de beta_cura0
  #deriv2[1,3] <- mean(colSums( - X_cure[,1]*pred_cure*exp(-risco_a_T)*log(t)*risco_a_T)) 
  #deriv2[3,1] <- deriv2[1,3]
  #
  ##derivada de alpha.t da derivada de beta_cura1
  #deriv2[2,3] <- mean(colSums( - X_cure[,2]*pred_cure*exp(-risco_a_T)*log(t)*risco_a_T)) 
  #deriv2[3,2] <- deriv2[2,3]
  
  # ============================== #
  
  for (i in 1:p){
      deriv2[i,p+2] <- mean(colSums( - X_cure[,i]*pred_cure*exp(-risco_a_T)*(t^alpha_T)))    #derivada de lambda.t da derivada de beta_cura0
      deriv2[p+2,i] <- deriv2[i,p+2]
  }
  
  #derivada de lambda.t da derivada de beta_cura0
  #deriv2[1,4] <- mean(colSums( - X_cure[,1]*pred_cure*exp(-risco_a_T)*(t^alpha_T))) 
  #deriv2[4,1] <- deriv2[1,4]
  #
  ##derivada de lambda.t da derivada de beta_cura1
  #deriv2[2,4] <- mean(colSums( - X_cure[,2]*pred_cure*exp(-risco_a_T)*(t^alpha_T))) 
  #deriv2[4,2] <- deriv2[2,4]
  
  # ============================== #
  
  #derivada de alpha.t da derivada de alpha.t
  deriv2[p+1,p+1] <- mean(colSums(Delta_t*(-alpha_T^(-2) - risco_a_T*log(t)*log(t)) + (pred_cure*exp(-risco_a_T)*(-risco_a_T)*log(t)*log(t)*(-risco_a_T))))
  
  #derivada de lambda.t da derivada de alpha.t
  deriv2[p+1,p+2] <- mean(colSums(Delta_t*(-(t^alpha_T)*log(t)) + (pred_cure*exp(-risco_a_T)*(-risco_a_T)*(-t^alpha_T)*log(t))))
  deriv2[p+2,p+1] <- deriv2[p+1,p+2]
  
  #derivada de lambda.t da derivada de lambda.t
  deriv2[p+2,p+2] <- mean(colSums(Delta_t*(-lambda_T^(-2)) + pred_cure*exp(-risco_a_T)*(t^alpha_T)*(t^alpha_T)))
  
  # ============================== #
  
  for (i in 1:q){
    for (u in 1:q) {
      deriv2[p+2+u,p+2+i] <- mean(colSums(- X_C[,u]*X_C[,i]*risco_a_C*pred_C))   #derivada do vetor beta_T da derivada de alpha.t
      deriv2[p+2+i,p+2+u] <- deriv2[p+2+u,p+2+i]
    }
  }
  
  ##derivada de beta_C_1 da derivada de beta_C_1
  #deriv2[5,5] <- mean(colSums(- X_C[,1]*X_C[,1]*risco_a_C*pred_C))
  #
  ##derivada de beta_C_2 da derivada de beta_C_2
  #deriv2[6,6] <- mean(colSums(- X_C[,2]*X_C[,2]*risco_a_C*pred_C))
  #
  ##derivada de beta_C_2 da derivada de beta_C_1
  #deriv2[5,6] <- mean(colSums(- X_C[,1]*X_C[,2]*risco_a_C*pred_C))
  #deriv2[6,5] <- deriv2[5,6]    
  
  # ============================== #
  
  for (i in 1:q){
      deriv2[p+2+i,p+q+3] <- mean(colSums(- X_C[,i]*wk*risco_a_C*pred_C))   #derivada de alpha da derivada de beta_C_1
      deriv2[p+q+3,p+2+i] <- deriv2[p+2+i,p+q+3]
    }
  
  #derivada de alpha da derivada de beta_C_1
  #deriv2[5,7] <- mean(colSums(- X_C[,1]*wk*risco_a_C*pred_C))
  #deriv2[7,5] <- deriv2[5,7]
  #
  ##derivada de alpha da derivada de beta_C_2
  #deriv2[6,7] <- mean(colSums(- X_C[,2]*wk*risco_a_C*pred_C))
  #deriv2[7,6] <- deriv2[6,7]
  
  # ============================== #
  
  for (i in 1:q){
    deriv2[p+2+i,p+q+4] <- mean(colSums(- X_C[,i]*risco_a_C*log(t)*pred_C))   #derivada de alpha da derivada de beta_C_1
    deriv2[p+q+4,p+2+i] <- deriv2[p+2+i,p+q+4]
  }
  
  #derivada de alpha.c da derivada de beta_C_1
  #deriv2[5,8] <- mean(colSums(- X_C[,1]*risco_a_C*log(t)*pred_C))
  #deriv2[8,5] <-  deriv2[5,8]
  #
  ##derivada de alpha.c da derivada de beta_C_2
  #deriv2[6,8] <- mean(colSums(- X_C[,2]*risco_a_C*log(t)*pred_C))
  #deriv2[8,6] <- deriv2[6,8]
  
  # ============================== #
  
  for (i in 1:q){
    deriv2[p+2+i,p+q+5] <- mean(colSums(- X_C[,i]*(t^alpha_C)*pred_C))   #derivada de alpha da derivada de beta_C_1
    deriv2[p+q+5,p+2+i] <- deriv2[p+2+i,p+q+5]
  }
  
  ##derivada de lambda.c da derivada de beta_C_1
  #deriv2[5,9] <- mean(colSums(- X_C[,1]*(t^alpha_C)*pred_C))
  #deriv2[9,5] <- deriv2[5,9]
  #
  ##derivada de lambda.c da derivada de beta_C_2
  #deriv2[6,9] <- mean(colSums(- X_C[,2]*(t^alpha_C)*pred_C))
  #deriv2[9,6] <-  deriv2[6,9]
  
  # ============================== #
  #derivada de alpha da derivada de alpha
  deriv2[p+q+3,p+q+3] <- mean(colSums(- risco_a_C*pred_C*wk*wk))
  #
  ##derivada de alpha.c da derivada de alpha  
  deriv2[p+q+3,p+q+4] <- mean(colSums( - risco_a_C*log(t)*pred_C*wk)) 
  deriv2[p+q+4,p+q+3] <-  deriv2[p+q+3,p+q+4]
  #
  ##derivada de lambda.c da derivada de alpha 
  deriv2[p+q+3,p+q+5] <- mean(colSums(- (t^alpha_C)*pred_C*wk))
  deriv2[p+q+5,p+q+3] <-  deriv2[p+q+3,p+q+5]
  #
  ##derivada de alpha.c da derivada de alpha.c
  deriv2[p+q+4,p+q+4] <- mean(colSums(Delta_c*(-alpha_C^(-2)) - risco_a_C*pred_C*log(t)*log(t))) 
  #
  ##derivada de lambda.c da derivada de alpha.c
  deriv2[p+q+4,p+q+5] <- mean(colSums( - (t^alpha_C)*log(t)*pred_C)) 
  deriv2[p+q+5,p+q+4] <- deriv2[p+q+4,p+q+5]
  #
  ##derivada de lambda.c da derivada de lambda.c
  deriv2[p+q+5,p+q+5] <- mean(colSums(Delta_c*(-lambda_C^(-2)))) 
  #
  ##derivada de theta da derivada de theta
  deriv2[p+q+6,p+q+6] <- mean(-0.5*(-m*(theta^(-2)) + (theta^(-3))*2*colSums(w_k_grupo^(2))))
  
  #derivada de alpha da derivada de alpha
  #deriv2[7,7] <- mean(colSums(- risco_a_C*pred_C*wk*wk))
  #
  ##derivada de alpha.c da derivada de alpha  
  #deriv2[7,8] <- mean(colSums( - risco_a_C*log(t)*pred_C*wk)) 
  #deriv2[8,7] <-  deriv2[7,8]
  #
  ##derivada de lambda.c da derivada de alpha 
  #deriv2[7,9] <- mean(colSums(- (t^alpha_C)*pred_C*wk))
  #deriv2[9,7] <-  deriv2[7,9]
  #
  ##derivada de alpha.c da derivada de alpha.c
  #deriv2[8,8] <- mean(colSums(Delta_c*(-alpha_C^(-2)) - risco_a_C*pred_C*log(t)*log(t))) 
  #
  ##derivada de lambda.c da derivada de alpha.c
  #deriv2[8,9] <- mean(colSums( - (t^alpha_C)*log(t)*pred_C)) 
  #deriv2[9,8] <- deriv2[8,9]
  #
  ##derivada de lambda.c da derivada de lambda.c
  #deriv2[9,9] <- mean(colSums(Delta_c*(-lambda_C^(-2)))) 
  #
  ##derivada de theta da derivada de theta
  #deriv2[10,10] <- mean(-0.5*(-m*(theta^(-2)) + (theta^(-3))*2*colSums(w_k_grupo^(2))))
  
  # ============================== #
  
  return((as.matrix(deriv2)))
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o para o c?lculo dos crit?rios
###--------------------------------------------------------------------------------------------------------------------------------------###
#como muda logver aqui tbm muda
crit_weibull <- function(X_cure,X_C,bi,n,alpha_T,lambda_T,beta_cura,alpha_C,lambda_C,beta_C,alpha,delta_t,delta_c,time,m,ident) {
  w <- bi                                                                       # pega a matriz de fragilidades salva na ultima iteracao
  L <- ncol(w)
  num_param <- length(c(alpha,beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C)) # numero de parametros, parametros salvos em fit (fit <- c((out[s,]+out[s-1,]+out[s-2,])/3))
  log_vero <- matrix(NA,n,L)                                                    #log-verossimilhan?a, L numero de replicas monte carlo de w
  
  pred <- X_cure%*%beta_cura
  a_T <- alpha_T
  b_T <- lambda_T
  exp_T <- (time^a_T)*b_T
  
  pred_C <- X_C%*%beta_C     
  a_C <- alpha_C
  b_C <- lambda_C
  exp_C <- (time^a_C)*b_C*exp(pred_C)
  
  for ( l in 1:L){
       log_vero[,l] <- delta_t*(pred + w[,l] + log(a_T) + (a_T-1)*log(time) + log(b_T) - exp_T) - exp(pred)*exp(w[,l])*(1-exp(-exp_T))
                      + delta_c*(log(a_C) + (a_C-1)*log(time) + log(b_C) + pred_C + alpha*w[,l])- exp_C*exp(alpha*w[,l])  
  }
  
  vero <- exp(log_vero)
  vero_grupo <- matrix(NA,m,L)  #m eh o numero de grupos/cluster
  
  for (i in 1:m){
    vero_grupo[i,] <- matrixStats::colProds(vero, rows = which(ident==i))  # verossimilhanca para cada grupo
  }
  
  log_lik <- sum(log(rowSums(vero_grupo)/L))
  
  AIC <- 2*(-log_lik + num_param)
  
  BIC <- 2*(-log_lik + 0.5*log(n)*num_param)
  
  HQ <- 2*(-log_lik + log(log(n))*num_param)
  
  return(cbind(AIC,BIC,HQ))
}
###--------------------------------------------------------------------------------------------------------------------------------------###




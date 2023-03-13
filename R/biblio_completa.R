###--------------------------------------------------------------------------------------------------------------------------------------###
### funcao para gerar da log condicional completa da fragilidade, wk (sendo k o indice de grupo)                                         ###
###--------------------------------------------------------------------------------------------------------------------------------------###
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
###-----------------------------------------------------------------------------------------------------------------------------------------
support_wi <-  function(w_k_aux,risco_a_T, risco_a_C, delta_T,delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta){(w_k_aux>0)*(w_k_aux<1)}

###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o log-verossimilhan?a, parte da falha, usada no optim para estimar os par?metros da distribui??o Weibull                        ###
### e da fra??o de cura, considerando covari?veis e fragilidade na fra??o de cura                                                        ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_param_T_Weib <-  function(par,t,delta_T, X_cure, bi){
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

###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o log-verossimilhan?a utilizada no optim() para estimar os parametros beta_cure e lambdas_j dos tempos de falha T                 ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_param_T_MEP <-  function(param_T,delta_t, xi_T,X_cure,num_int_T,bi){   #lambda_T=param_T
  colcura <- ncol(X_cure)
  n <- length(X_cure[,1])
  Delta <- delta_t%*%t(rep(1,num_int_T))
  lambda <- t((param_T[1:num_int_T])%*%t(rep(1,n)))  #estima os lambdas_j

  if (ncol(t(t(X_cure))) == 1){
    x_cure_b <- t(t(X_cure))*param_T[(num_int_T+1):(num_int_T+colcura)]   #estima os betas
  }
  if (ncol(t(t(X_cure))) != 1){
    x_cure_b <- X_cure%*%param_T[(num_int_T+1):(num_int_T+colcura)]   #estima os betas
  }

  #x_cure_b <- X_cure%*%param_T[(num_int_T+1):(num_int_T+colcura)]   #estima os betas
  prop <- exp(x_cure_b)*rowMeans(exp(bi))

  aux <- apply(exp(lambda)*xi_T,1,sum)
  soma <- aux%*%t(rep(1,num_int_T))

  aux2 <- as.vector(prop)*(1- exp(-apply(exp(lambda)*xi_T,1,sum)))
  aux3 <- aux2%*%t(rep(1,num_int_T))

  mean_bi <- rowMeans(bi)%*%t(rep(1,num_int_T))

  U <- -sum(Delta*(as.vector(x_cure_b) + lambda  + mean_bi - soma) - aux3)

}

###------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o log-verossimilhan?a, parte da censura informativa, usada no optim para estimar os par?metros da distribui??o Weibull        ###
### considerando covari?veis na fun??o taxa de falha de C                                                                              ###
###------------------------------------------------------------------------------------------------------------------------------------###
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
### Fun??o log-verossimilhan?a utilizada no multiroot() para estimar os parametros beta_C e alpha tempos de censura informativa C        ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_C_MEP <-  function(beta_C, X_C,delta_c,risco_a_C,bi,n_intMC){
  q <- ncol(X_C)
  n <- nrow(X_C)

  if(ncol(t(t(X_C))) == 1){
    w_kl_beta_C <- risco_a_C*exp((X_C[,]*beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]*beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
  }

  else{
    w_kl_beta_C <- risco_a_C*exp((X_C[,]%*%beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]%*%beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
  }

  w_kl_beta_C_num <- cbind(matrix(w_kl_beta_C, nrow = n, ncol = q))*X_C
  U_C_1 <- colSums((X_C*delta_c - w_kl_beta_C_num))

  U_C_alpha <- sum(delta_c*(rowSums(bi[,])/n_intMC) - w_kl_beta_C_alpha_num)

  c(U_C_1 = U_C_1,U_C_alpha=U_C_alpha)
}

###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o que particiona o eixo do tempo em intervalos, tamb?m retorna a indicadora dos intervalos (observa??es de cada intervao)       ###
###--------------------------------------------------------------------------------------------------------------------------------------###
time.grid <- function(time, event, n.int=NULL)
{
  o <- order(time)  #ordeno os dados
  time <- time[o]
  event <- event[o]
  time.aux <- unique(time[event==1])
  if(is.null(n.int))
  {
    n.int <- length(time.aux)
  }

  m <- length(time.aux)
  if(n.int > m)
  {
    a <- c(0,unique(time[event==1]))
    a[length(a)] <- Inf
  }
  else
  {
    b <- min(m,n.int)
    k1 <- trunc(m/b)
    r <- m-b*k1
    k2 <- k1+1
    idf1 <- seq(k1,(b-r)*k1, k1)
    idf2 <- sort(seq(m,max(idf1),-k2))
    idf <- unique(c(idf1,idf2))
    a_inf <- c(0,time.aux[idf])
    a_inf[length(a_inf)] <- Inf
    a_s_inf  <- c(0,time.aux[idf])
  }
  saida <- list(a_inf,a_s_inf)

  return(saida)
}

###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o com o c?lculo das derivadas de primeira ordem. Retorna a esperan?a da matriz de drivadas de primeira ordem.                   ###
### C?lculo ? feito como no artigo de Louis et al (1982)                                                                                 ###
###--------------------------------------------------------------------------------------------------------------------------------------###
Esp_Deriv_Prim_Ordem_MEP <-  function( X_cure, X_C,delta_T,delta_C,id_T, beta_cure, beta_C, alpha, theta,nu_ind_j_T,nu_ind_j_C, lambda_T_j,lambda_C_j, deltaij_T,deltaij_C, w_k_grupo, ident){

  nc <- nrow(X_cure)
  p <- ncol(X_cure)
  q <- ncol(X_C)

  lambdalengthJ <- length(lambda_T_j)

  n <- nrow(X_C)
  wk = w_k_grupo[ident,]
  m <- nrow(w_k_grupo)
  deriv1 <- matrix(NA,(p+q+2+length(lambda_T_j)+length(lambda_C_j)),ncol(wk))  #matriz que com as derivadas de primeira ordem,

  if(ncol(t(t(X_cure))) == 1){
    pred_cure <- as.vector(exp(t(t(X_cure))*beta_cure))*exp(wk)
  }

  else{
    pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  }

  if(ncol(t(t(X_C))) == 1){
    pred_C <- as.vector(exp(t(t(X_C))*beta_C))*exp(alpha*wk)
  }

  else{
    pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  }


  #pred_cure <- as.vector(exp(X_cure%*%beta_cure))*exp(wk)
  #pred_C <- as.vector(exp(X_C%*%beta_C))*exp(alpha*wk)
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  lambda_c <- t(lambda_C_j%*%t(rep(1,n)))
  risco_a_T <- apply(t(as.vector(lambda_T_j)*t(deltaij_T)),1,sum)
  risco_a_C <- apply(t(as.vector(lambda_C_j)*t(deltaij_C)),1,sum)

  for (i in 1:p){
    deriv1[i,] <- colSums(X_cure[,i]*(Delta_t - pred_cure*(1-exp(-risco_a_T))))  #derivadas de beta_T's
  }

  #derivada de beta_cure_0
  #deriv1[1,] <- colSums(X_cure[,1]*(Delta_t - pred_cure*(1-exp(-risco_a_T))))

  #derivada de beta_cure_1
  #deriv1[2,] <- colSums(X_cure[,2]*(Delta_t - pred_cure*(1-exp(-risco_a_T))))

  #derivada de lambda_T_j
  for ( j in 1:length(lambda_T_j)){
    deriv1[(j+p),] <- nu_ind_j_T[j]*(lambda_T_j[j]^(-1)) - sum(delta_T[id_T==j]*deltaij_T[id_T==j,j]) - colSums(pred_cure*exp(-risco_a_T)*deltaij_T[,j])
  }

  for (i in 1:q){
    deriv1[(length(lambda_T_j)+p+i),] <- colSums(X_C[,i]*(Delta_c - risco_a_C*pred_C))  #derivadas de beta_C's
  }

  #derivada de beta_C_1
  #deriv1[length(lambda_T_j)+3,] <- colSums(X_C[,1]*(Delta_c - risco_a_C*pred_C))

  #derivada de beta_C_2
  #deriv1[length(lambda_T_j)+4,] <- colSums(X_C[,2]*(Delta_c - risco_a_C*pred_C))

  #derivada de alpha (12)
  deriv1[(length(lambda_T_j)+p+q+1),] <- colSums(wk*(Delta_c - risco_a_C*pred_C))

  #derivada de lambda_C_j
  for ( j in 1:length(lambda_C_j)){
    deriv1[(length(lambda_T_j)+p+q+j+1),] <- nu_ind_j_C[j]*(lambda_C_j[j]^(-1)) - colSums(deltaij_C[,j]*pred_C)
  }

  #derivada de tau (17)
  deriv1[length(lambda_C_j)+length(lambda_T_j)+p+q+2,] <- -0.5*(m*theta^(-1) - (theta^(-2))*colSums(w_k_grupo^(2)))

  aux <- deriv1[,1]%*%t(deriv1[,1])

  for( i in 2:ncol(wk)){
    aux <- aux + deriv1[,i]%*%t(deriv1[,i])
  }

  return(aux/ncol(wk))
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o com o c?lculo das derivadas de segunda ordem. Retorna a esperan?a da matriz de drivadas de segunda ordem.                     ###
### C?lculo ? feito como no artigo de Louis et al (1982)                                                                                 ###
###--------------------------------------------------------------------------------------------------------------------------------------###
Esp_DerivParciais_MEP  <-  function( X_cure, X_C,delta_T,delta_C, beta_cure, beta_C, alpha, theta,nu_ind_j_T,nu_ind_j_C, lambda_T_j,lambda_C_j, deltaij_T,deltaij_C, w_k_grupo, ident){
  nc <- nrow(X_cure)
  p <- ncol(X_cure)
  q <- ncol(X_C)

  n <- nrow(X_C)
  wk = w_k_grupo[ident,]
  m <- nrow(w_k_grupo)
  deriv2 <- matrix(0,(p+q+2+length(lambda_T_j)+length(lambda_C_j)),(p+q+2+length(lambda_T_j)+length(lambda_C_j)))

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
  Delta_t <- delta_T%*%t(rep(1,ncol(wk)))
  Delta_c <- delta_C%*%t(rep(1,ncol(wk)))
  lambda_t <- t(lambda_T_j%*%t(rep(1,n)))
  lambda_c <- t(lambda_C_j%*%t(rep(1,n)))
  risco_a_T <- apply(t(as.vector(lambda_T_j)*t(deltaij_T)),1,sum)
  risco_a_C <- apply(t(as.vector(lambda_C_j)*t(deltaij_C)),1,sum)

  #derivada de beta_T_1 da derivada de beta_T_1
  for (i in 1:p){
    for (u in 1:p) {
      deriv2[u,i] <- mean(colSums( - X_cure[,u]*X_cure[,i]*pred_cure*(1-exp(-risco_a_T))))   #derivada do vetor beta_T da derivada de alpha.t
      deriv2[i,u] <- deriv2[u,i]
    }
  }

  #derivada de beta_T_1 da derivada de beta_T_1
  #deriv2[1,1] <- - mean(colSums(X_cure[,1]*X_cure[,1]*pred_cure*(1-exp(-risco_a_T))))
  #
  ##derivada de beta_T_2 da derivada de beta_T_2
  #deriv2[2,2] <- - mean(colSums(X_cure[,2]*X_cure[,2]*pred_cure*(1-exp(-risco_a_T))))
  #
  ##derivada de beta_T_2 da derivada de beta_T_1
  #deriv2[1,2] <- - mean(colSums(X_cure[,1]*X_cure[,2]*pred_cure*(1-exp(-risco_a_T))))
  #deriv2[2,1] <- deriv2[1,2]

  #derivada de lambda_T_j da derivada de beta_T_1
  for (i in 1:p) {
    for (j in 1:length(lambda_T_j)){
      deriv2[i,(j+p)] <- - mean(colSums(X_cure[,i]*pred_cure*exp(-risco_a_T)*deltaij_T[,j]))
      deriv2[(j+p),i] <- deriv2[i,(j+p)]
    }
  }

  #for (j in 1:length(lambda_T_j)){
  #  deriv2[1,(j+2)] <- - mean(colSums(X_cure[,1]*pred_cure*exp(-risco_a_T)*deltaij_T[,j]))
  #  deriv2[(j+2),1] <- deriv2[1,(j+2)]
  #}
  #
  ##derivada de lambda_T_j da derivada de beta_T_2
  #for (j in 1:length(lambda_T_j)){
  #  deriv2[2,(j+2)] <- - mean(colSums(X_cure[,2]*pred_cure*exp(-risco_a_T)*deltaij_T[,j]))
  #  deriv2[(j+2),2] <- deriv2[2,(j+2)]
  #}

  #derivada de lambda_T_j da derivada de lambda_T_j
  for (j in 1:length(lambda_T_j)){
    deriv2[(j+p),(j+p)] <- - nu_ind_j_T[j]*(lambda_T_j[j]^(-2)) + mean(colSums(pred_cure*exp(-risco_a_T)*deltaij_T[,j]*deltaij_T[,j]))
  }

  #derivada de beta_C _1 da derivada de beta_C_1

  #derivada de beta_C_1 da derivada de beta_C_1

  for (i in 1:q){
    for (u in 1:q) {
      deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+u)] <- - mean(colSums(X_C[,i]*X_C[,u]*risco_a_C*pred_C))
      deriv2[(length(lambda_T_j)+p+u),(length(lambda_T_j)+p+i)] <- deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+u)]
    }
  }

  #derivada de beta_C_1 da derivada de beta_C_1
  #deriv2[(length(lambda_T_j)+3),(length(lambda_T_j)+3)] <- - mean(colSums(X_C[,1]*X_C[,1]*risco_a_C*pred_C))

  #derivada de beta_C_2 da derivada de beta_C_2
  #deriv2[(length(lambda_T_j)+4),(length(lambda_T_j)+4)] <- - mean(colSums(X_C[,2]*X_C[,2]*risco_a_C*pred_C))

  # ========= ASSUMINDO QUE O CODIGO DA LINHA 337 JÁ FAZ ESSES CASOS !!!!!!!! ============
  #derivada de beta_C_2 da derivada de beta_C_1
  #for (j in 1:q){
  #  for (i in 2:(q)){
  #    if (ncol(t(t(X_C))) == 1){
  #      next
  #    }
  #    deriv2[(length(lambda_T_j)+p+j),(length(lambda_T_j)+p+1+i)] <- - mean(colSums(X_C[,j]*X_C[,i+1]*risco_a_C*pred_C))
  #    deriv2[(length(lambda_T_j)+p+1+i),(length(lambda_T_j)+p+j)] <-  deriv2[(length(lambda_T_j)+p+j),(length(lambda_T_j)+p+1+i)]   #derivada de beta_C_i da derivada de beta_C_(i-1)
  #  }
  #}

  #derivada de beta_C_2 da derivada de beta_C_1
  #deriv2[(length(lambda_T_j)+3),(length(lambda_T_j)+4)] <- - mean(colSums(X_C[,1]*X_C[,2]*risco_a_C*pred_C))
  #deriv2[(length(lambda_T_j)+4),(length(lambda_T_j)+3)] <- deriv2[(length(lambda_T_j)+3),(length(lambda_T_j)+4)]
  for (i in 1:q){
    deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+q+1)]  <- - mean(colSums(X_C[,i]*wk*risco_a_C*pred_C)) #derivada de alpha de beta_C
    deriv2[(length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+i)]  <- deriv2[(length(lambda_T_j)+p+i),(length(lambda_T_j)+p+q+1)]
  }

  #derivada de alpha da derivada de beta_C_1
  #deriv2[(length(lambda_T_j)+3),(length(lambda_T_j)+5)] <- - mean(colSums(X_C[,1]*wk*risco_a_C*pred_C))
  #deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+3)] <- deriv2[(length(lambda_T_j)+3),(length(lambda_T_j)+5)]

  #derivada de alpha da derivada de beta_C_2
  #deriv2[(length(lambda_T_j)+4),(length(lambda_T_j)+5)] <- - mean(colSums(X_C[,2]*wk*risco_a_C*pred_C))
  #deriv2[(length(lambda_T_j)+5),(length(lambda_T_j)+4)] <- deriv2[(length(lambda_T_j)+4),(length(lambda_T_j)+5)]

  #derivada de alpha da derivada de alpha
  deriv2[(length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+1)] <- - mean(colSums(wk*wk*risco_a_C*pred_C))

  #derivada de lambda_C_j da derivada de lambda_C_j
  for (j in 1:length(lambda_C_j)){
    deriv2[(j+length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)] <- - nu_ind_j_C[j]*(lambda_C_j[j]^(-2))
  }

  #derivada de beta_C_1 da derivada de lambda_C_j
  for (i in 1:q) {
    for (j in 1:length(lambda_C_j)){
      deriv2[(length(lambda_T_j)+p+i),(j+length(lambda_T_j)+p+q+1)] <- - mean(colSums(deltaij_C[,j]*pred_C*X_C[,i]))
      deriv2[(j+length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+i)] <- deriv2[(length(lambda_T_j)+p+i),(j+length(lambda_T_j)+p+q+1)]
    }
  }

  #derivada de beta_C_2 da derivada de lambda_C_j
  #for (j in 1:length(lambda_C_j)){
  #  deriv2[(length(lambda_T_j)+4),(j+length(lambda_T_j)+5)] <- -mean(colSums(deltaij_C[,j]*pred_C*X_C[,2]))
  #  deriv2[(j+length(lambda_T_j)+5),(length(lambda_T_j)+4)] <- deriv2[(length(lambda_T_j)+4),(j+length(lambda_T_j)+5)]
  #}

  #derivada de alpha da derivada de lambda_C_j
  for (j in 1:length(lambda_C_j)){
    deriv2[(length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)] <- - mean(colSums(deltaij_C[,j]*pred_C*wk))
    deriv2[(j+length(lambda_T_j)+p+q+1),(length(lambda_T_j)+p+q+1)] <- deriv2[(length(lambda_T_j)+p+q+1),(j+length(lambda_T_j)+p+q+1)]
  }

  #derivada de theta da derivada de theta
  deriv2[(length(lambda_T_j)+length(lambda_C_j)+p+q+2),(length(lambda_T_j)+length(lambda_C_j)+p+q+2)] <-  -0.5*(-m*(theta^(-2)) + (theta^(-3))*2*mean(colSums(w_k_grupo^(2))))

  return((as.matrix(deriv2)))
}
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o com o c?lculo dos crit?rios de sela??o de modelos.                                                                            ###
###--------------------------------------------------------------------------------------------------------------------------------------###
crit_MEP <- function(X_cure,X_C,bi,n,beta_cure,lambda_T_j,beta_C,alpha,lambda_C_j,delta_t,delta_c,risco_a_T,risco_a_C,id_T,id_C,m,ident) {
  w <- bi                                                           # pega a matriz de fragilidades salva na ultima itera??o
  L <- ncol(w)
  num_param <- length(c(beta_cure,lambda_T_j,beta_C,alpha,lambda_C_j)) # numero de parametros, parametros salvos em fit (fit <- c((out[s,]+out[s-1,]+out[s-2,])/3))
  log_vero <- matrix(NA,n,L)                                        # log-verossimilhan?a, L numero de replicas monte carlo de w


  for (l in 1:L){
    if (ncol(t(t(X_cure))) == 1 && ncol(t(t(X_C))) == 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure*beta_cure + w[,l] - risco_a_T) - exp(X_cure*beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C*beta_C + alpha*w[,l]) - risco_a_C*exp(X_C*beta_C + alpha*w[,l])
    }

    if (ncol(t(t(X_cure))) == 1 && ncol(t(t(X_C))) != 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure*beta_cure + w[,l] - risco_a_T) - exp(X_cure*beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C%*%beta_C + alpha*w[,l]) - risco_a_C*exp(X_C%*%beta_C + alpha*w[,l])
    }

    if (ncol(t(t(X_cure))) != 1 && ncol(t(t(X_C))) == 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure%*%beta_cure + w[,l] - risco_a_T) - exp(X_cure%*%beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C*beta_C + alpha*w[,l]) - risco_a_C*exp(X_C*beta_C + alpha*w[,l])

    }
    if (ncol(t(t(X_cure))) != 1 && ncol(t(t(X_C))) != 1){

      log_vero[,l] <-  delta_t*(t(t(log(lambda_T_j[id_T]))) + X_cure%*%beta_cure + w[,l] - risco_a_T) - exp(X_cure%*%beta_cure + w[,l])*(1-exp(-risco_a_T))
      + delta_c*(t(t(log(lambda_C_j[id_C]))) + X_C%*%beta_C + alpha*w[,l]) - risco_a_C*exp(X_C%*%beta_C + alpha*w[,l])
    }
  }

  vero <- exp(log_vero)
  vero_grupo <- matrix(NA,m,L)  #m eh o numero de grupos/cluster

  for (i in 1:m){
    vero_grupo[i,] <- matrixStats::colProds(vero, rows = which(ident==i))  # verossimilhan?a para cada grupo
  }

  log_lik <- sum(log(rowSums(vero_grupo)/L))

  AIC <- 2*(-log_lik + num_param)
  BIC <- 2*(-log_lik + 0.5*log(n)*num_param)
  HQ <- 2*(-log_lik + log(log(n))*num_param)

  return(cbind(AIC,BIC,HQ))
}

###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??es com as derivadas de primeira ordem
###--------------------------------------------------------------------------------------------------------------------------------------###
#aqui tbm muda primeiro parametros de falha
Esp_Deriv_Prim_Ordem_Weib <-  function(t,delta_T, delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta, alpha_T, alpha_C, lambda_T, lambda_C, w_k_grupo, ident){

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
Esp_DerivParciais_Weib <-  function(t,delta_T, delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta, alpha_T, alpha_C, lambda_T, lambda_C, w_k_grupo, ident){

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

model_MEP_dep <-  function(formula, data, delta_t, delta_c, ident, Num_intervals){

  #INICIO FUNÇÃO
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- stats::model.matrix(formula, data = mf, rhs = 2)
  X_cura <- Z[,-1]
  #X_cura <- cbind(int = 1, X_cura) #Adicionando a coluna correspondente a beta_0
  X_cura <- as.matrix(X_cura)
  X_cura <- cbind(1,X_cura)
  X_C <- X[,-1]
  X_C <- as.matrix(X_C)
  Xlabels <- colnames(X_cura)
  Zlabels <- colnames(X_C)
  time <- stats::model.response(mf)
  # usando 't' ao invés de time
  t <- stats::model.response(mf)
  q <- ncol(X_C)
  p <- ncol(X_cura)
  n <- nrow(data)
  m <- max(ident)

  ###--------------------------------------------------------------------------------------
  # chute inicial para os parametros
  n <- 400
  beta_cura <- rep(0.1,p)
  beta_C <- rep(0.1,q)
  theta <- 1
  alpha <- 0
  param <- c(alpha, beta_cura,beta_C,theta)
  risco_a_T <- rep(1,n)
  risco_a_C <- rep(1,n)
  bmax <- 5
  lambda_T_j <- rep(0.01,bmax)
  ###----------------------------------------------------------------------------------------------------
  # Especificacoes do algorimo EMMC
  maxit <- 200 #numero maximo de iteracoes
  eps1= rep(1e-7, length(param))
  eps2= rep(1e-8, length(param))
  n_intMCs = c(rep(10,20),rep(20,30),rep(50,50),rep(75,50),rep(100,30),rep(125,20)) #numero de replicas da fragilidade em cada iteracao

  ## Iniciando os objetos utilizados
  out =  matrix(NA,maxit+1,length(param))
  dif =  matrix(NA,maxit+1,length(param))
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE

  while (continue == TRUE) {

    count = rep(0,maxit+1)
    out[s,] =  c(alpha, beta_cura,beta_C,theta)
    n_intMC = n_intMCs[s]
    m <- max(ident)
    w_chapeu_grupo <- matrix(NA, m, n_intMC)

    for ( k in 1:m){  #k eh grupo
      w_trans <- arms(0.5, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,risco_a_T=risco_a_T[ident==k], risco_a_C=risco_a_C[ident==k], delta_T=delta_t[ident==k],delta_C=delta_c[ident==k], X_cure=X_cura[ident==k,],X_C=X_C[ident==k,], beta_cure=beta_cura,beta_C=beta_C,alpha=alpha, theta=theta)
      w_auxi <- log(w_trans)-log(1-w_trans)
      w_chapeu_grupo[k,] <- w_auxi
    }
    bi = w_chapeu_grupo[ident,]

    theta <- mean(w_chapeu_grupo^2)
    #print(theta)
    ###----------------------------------------------------------------------------------------
    a <- time.grid(t,delta_t,bmax); a <- a[[1]]; #grade dos tempos de falha
    c <- time.grid(t,delta_c,bmax); c <- c[[1]]; #grade dos tempos de censura
    b <- length(a)-1
    d <- length(c)-1
    a[length(a)] <- max(t)
    c[length(c)] <- max(t)
    num_int_T <- bmax
    num_int_C <- bmax

    id_T <- as.numeric(cut(t,a))  # id sempre ? com grade com inf na ?ltima casela
    nu_T <- tapply(delta_t,id_T ,sum)

    xi.falha <- matrix(0,n,b)
    xi.cens <- matrix(0,n,d)

    for(i in 1:n){
      for(j in 1:b){
        xi.falha[i,j] <- (min(t[i], a[j+1])-a[j])*((t[i]-a[j])>0)
      }
    }
    id.falha <- as.numeric(cut(t,a))
    id.cens  <- as.numeric(cut(t,c))
    for(i in 1:n){
      for(j in 1:d){
        xi.cens[i,j] <- (min(t[i], c[j+1])-c[j])*((t[i]-c[j])>0)
      }
    }

    param_T <- optim(par = c(lambda_T_j, beta_cura), fn = modelo_param_T_MEP, control = list(maxit = 5000), method = c( "BFGS"), delta_t=delta_t, xi_T=xi.falha,X_cure=X_cura,num_int_T=bmax, bi=bi)
    param_Ts <- param_T$par
    lambda_T_j <- exp(param_Ts[1:num_int_T])
    beta_cura  <- param_Ts[(num_int_T+1):(num_int_T+p)]
    risco_a_T <- apply(t(lambda_T_j*t(xi.falha)),1,sum)
    ###--------------------------------------------------------------------------------

    pred_linear_C <- exp(X_C%*%beta_C)*rowMeans(exp(alpha*bi[,]))
    id_C <- as.numeric(cut(t,c))  # id sempre ? com grade com inf na ?ltima casela
    nu_C <- tapply(delta_c,id_C ,sum)

    xi.cens_pred <- matrix(0,n,d)
    for(i in 1:n){
      for(j in 1:d){
        xi.cens_pred[i,j] <- (min(t[i], c[j+1])-c[j])*((t[i]-c[j])>0)*pred_linear_C[i]
      }
    }

    lambda_C_j <- nu_C/apply(xi.cens_pred,2,sum)
    risco_a_C <- apply(t(as.vector(lambda_C_j)*t(xi.cens)),1,sum)

    S_C <- multiroot(f = modelo_C_MEP, start = c(beta_C,alpha), X_C=X_C,delta_c=delta_c,risco_a_C=risco_a_C,bi=bi,n_intMC=n_intMC)
    betas_C <- S_C$root
    beta_C <- betas_C[1:q]
    alpha <- betas_C[q+1]
    ###-----------------------------------------------------------------------------------
    # criterio de parada
    out[s+1,]= c(alpha, beta_cura,beta_C,theta)
    print(out[s+1,])
    if (s>2) {
      dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)

      for (z in 1: length(param)){
        if (dif[s,z]<eps2[z]) {
          count[s] = count[s] + 1
        }
      }
    }
    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param)))
    }
  } ## fim do EMMC


  param_est <- c((out[s,]+out[s-1,]+out[s-2,])/3)

  ###-----------------------------------------------------------------------------------
  # c?lculo da esperan?a das derivadas para o c?lculo do Erro Padr?o segundo Louis et al (1982)
  Esp_deriv_ordem1 <- Esp_Deriv_Prim_Ordem_MEP( X_cure=X_cura, X_C=X_C,delta_T=delta_t,delta_C=delta_c,id_T=id_T, beta_cure=beta_cura, beta_C=beta_C, alpha=alpha,
                                                theta=theta,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j,
                                                deltaij_T=xi.falha,deltaij_C=xi.cens, w_k_grupo=w_chapeu_grupo, ident=ident)

  Esp_deriv_ordem2 <- Esp_DerivParciais_MEP( X_cure=X_cura, X_C=X_C,delta_T=delta_t,delta_C=delta_c, beta_cure=beta_cura, beta_C=beta_C, alpha=alpha,
                                             theta=theta,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j,
                                             deltaij_T=xi.falha,deltaij_C=xi.cens, w_k_grupo=w_chapeu_grupo, ident=ident)

  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  ErroPadrao <- sqrt(diag(Var))

  ###-----------------------------------------------------------------------------------
  # c?lculo dos crit?rios de sele??o de modelos
  criterios <- crit_MEP(X_cure=X_cura,X_C,bi,n,beta_cure=beta_cura,lambda_T_j,beta_C,alpha,lambda_C_j,delta_t=delta_t,delta_c=delta_c,risco_a_T,risco_a_C,id_T,id_C,m,ident)

  # Ajustando a classe/lista
  fit <- param_est
  fit <- list(fit=fit)
  fit$stde <- ErroPadrao
  fit$crit <- criterios
  #fit$pvalue <- c(p_value_alpha, p_value_t, p_value_c)
  fit$n <- n
  fit$p <- p
  fit$q <- q
  fit$call <- match.call()
  fit$formula <- stats::formula(Terms)
  fit$terms <- stats::terms.formula(formula)
  fit$labels1 <- Zlabels
  fit$labels2 <- Xlabels
  fit$risco_a_T <- risco_a_T
  fit$risco_a_C <- risco_a_C
  fit$bi <- bi
  fit$X_cura <- X_cura
  fit$X_C <- X_C
  fit$t <- t
  class(fit) <- "dcensoring"
  return(fit)
}

model_Weibull_dep <-  function(formula, data, delta_t, delta_c, ident){

  #INICIO FUNÇÃO
  formula <- Formula::Formula(formula)
  mf <- stats::model.frame(formula=formula, data=data)
  Terms <- stats::terms(mf)
  Z <- stats::model.matrix(formula, data = mf, rhs = 1)
  X <- stats::model.matrix(formula, data = mf, rhs = 2)
  X_cura <- Z[,-1]
  #X_cura <- cbind(int = 1, X_cura) #Adicionando a coluna correspondente a beta_0
  X_cura <- as.matrix(X_cura)
  X_cura <- cbind(1,X_cura)
  X_C <- X[,-1]
  X_C <- as.matrix(X_C)
  Xlabels <- colnames(X_cura)
  Zlabels <- colnames(X_C)
  time <- stats::model.response(mf)
  # usando 't' ao invés de time
  t <- stats::model.response(mf)
  q <- ncol(X_C)
  p <- ncol(X_cura)
  n <- nrow(data)
  m <- max(ident)

  #model_Weibull_dep(time ~ x1_cure + int | x_c1 + x_c2, data=dados2, delta_t=delta_t, delta_c=delta_c, ident=dados2$ident)


  #  o <- order(dados$time)  #ordeno os dados
  #  dados <- dados[o,]
  #  ident <- dados$ident
  #  X_cura<- cbind(dados$int,dados$x1_cure)
  #  X_C <- cbind(dados$x_c1, dados$x_c2)
  #  p <- ncol(X_cura)
  #  q <- ncol(X_C)
  #  t <- dados$time

  #  delta_t <- ifelse(dados$cens==1,1,0)
  #  delta_c <- ifelse(dados$cens==2,1,0)

  ###--------------------------------------------------------------------------------------
  # chute inicial
  beta_cura <- rep(0.1,p)
  beta_C <- rep(0.1,q)
  alpha_T <- 1
  lambda_T <- 1
  alpha_C <- 1
  lambda_C <- 1
  theta <- 1
  alpha <- 0
  param <- c(alpha, beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C,theta)
  risco_a_T <- rep(1,n)
  risco_a_C <- rep(1,n)

  ###----------------------------------------------------------------------------------------------------
  # Especificacoes do algorimo EMMC
  maxit <- 200 #numero maximo de iteracoes
  eps1= rep(1e-7, length(param))  #criterio de parada
  eps2= rep(1e-8, length(param))
  n_intMCs = c(rep(10,20),rep(20,30), rep(50,50), rep(75,50), rep(100,30),  rep(125,20)) #numero de replicas da fragilidade em cada iteracao

  ## Iniciando os objetos utilizados
  out =  matrix(NA,maxit+1,length(param))
  dif =  matrix(NA,maxit+1,length(param))
  final = length(param)
  count = rep(0,maxit+1)
  s=1
  continue=TRUE

  while (continue == TRUE) {

    count = rep(0,maxit+1)
    out[s,] =  c(alpha, beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C,theta)
    n_intMC = n_intMCs[s]
    m <- max(ident)
    w_chapeu_grupo <- matrix(NA, m, n_intMC)

    for ( k in 1:m){  #k eh grupo
      w_trans <- arms(0.4, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,risco_a_T=risco_a_T[ident==k], risco_a_C=risco_a_C[ident==k], delta_T=delta_t[ident==k],delta_C=delta_c[ident==k], X_cure=X_cura[ident==k,],X_C=X_C[ident==k,], beta_cure=beta_cura,beta_C=beta_C,alpha=alpha, theta=theta)
      w_auxi <- log(w_trans)-log(1-w_trans)
      w_chapeu_grupo[k,] <- w_auxi
    }
    bi = w_chapeu_grupo[ident,]

    theta <- mean(w_chapeu_grupo^2)
    ###----------------------------------------------------------------------------------------
    param_T <- optim(par = c(alpha_T,lambda_T, beta_cura), fn = modelo_param_T_Weib, control = list(maxit = 5000), method = c( "Nelder-Mead"),t=t,delta_T=delta_t, X_cure=X_cura, bi=bi)
    par_T <- param_T$par
    alpha_T <- par_T[1]
    lambda_T <- exp(par_T[2])
    beta_cura <- par_T[3:(2+p)] #AQUI ============================
    risco_a_T <- (t^alpha_T)*(lambda_T)

    ###--------------------------------------------------------------------------------
    param_C <- optim(par = c(alpha_C,lambda_C,alpha,beta_C), fn = modelo_param_C, control = list(maxit = 5000), method = c( "Nelder-Mead"),t=t,delta_C=delta_c, X_C=X_C, bi=bi)
    par_C <- param_C$par
    alpha_C <- par_C[1]
    lambda_C <- exp(par_C[2])
    alpha <- par_C[3]
    beta_C <- par_C[4:(3+q)]
    risco_a_C <- (t^alpha_C)*(lambda_C)

    ###-----------------------------------------------------------------------------------
    # criterio de parada
    out[s+1,]= c(alpha, beta_cura,beta_C,alpha_T,alpha_C,lambda_T,lambda_C,theta)
    print(out[s+1,])
    if (s>2) {
      dif[s,] <- (abs(out[s+1,]-out[s,]))/(abs(out[s,])-eps1)

      for (z in 1: length(param)){
        if (dif[s,z]<eps2[z]) {
          count[s] = count[s] + 1
        }
      }
    }

    s=s+1
    if (s>3) {
      continue=((s <= maxit) & (count[s] < length(param))&(count[s-1] < length(param))&(count[s-2] < length(param)))
    }
  } ## fim do EMMC

  param_est <- c((out[s,]+out[s-1,]+out[s-2,])/3)


  Esp_deriv_ordem1 <- Esp_Deriv_Prim_Ordem_Weib(t,delta_T=delta_t, delta_C=delta_c, X_cure=X_cura, X_C=X_C,
                                                beta_cure=beta_cura, beta_C=beta_C, alpha=alpha, theta=theta,
                                                alpha_T=alpha_T, alpha_C=alpha_C, lambda_T=lambda_T, lambda_C=lambda_C,
                                                w_k_grupo=w_chapeu_grupo, ident=ident)


  Esp_deriv_ordem2 <- Esp_DerivParciais_Weib(t,delta_T=delta_t, delta_C=delta_c, X_cure=X_cura, X_C=X_C,
                                             beta_cure=beta_cura, beta_C=beta_C, alpha=alpha, theta=theta,
                                             alpha_T=alpha_T, alpha_C=alpha_C, lambda_T=lambda_T, lambda_C=lambda_C,
                                             w_k_grupo=w_chapeu_grupo, ident=ident)


  InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
  Var <- solve(InfFisher)
  ErroPadrao <- sqrt(diag(Var))


  # calulo dos criterios
  criterios <- crit_weibull(X_cure=X_cura,X_C,bi,n,alpha_T,lambda_T,beta_cura,alpha_C,lambda_C,beta_C,alpha,delta_t=delta_t,delta_c=delta_c,time=t,m,ident)

  # Ajustando a classe/lista
  fit <- param_est
  fit <- list(fit=fit)
  fit$stde <- ErroPadrao
  fit$crit <- criterios
  #fit$pvalue <- c(p_value_alpha, p_value_t, p_value_c)
  fit$n <- n
  fit$p <- p
  fit$q <- q
  fit$call <- match.call()
  fit$formula <- stats::formula(Terms)
  fit$terms <- stats::terms.formula(formula)
  fit$labels1 <- Zlabels
  fit$labels2 <- Xlabels
  fit$risco_a_T <- risco_a_T
  fit$risco_a_C <- risco_a_C
  fit$bi <- bi
  fit$X_cura <- X_cura
  fit$X_C <- X_C
  fit$t <- t
  class(fit) <- "dcensoring"
  return(fit)
}

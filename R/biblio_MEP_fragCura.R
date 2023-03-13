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

###-----------------------------------------------------------------------------------------------
# funcao usada no arms para a log condicional completa da fragilidade 
support_wi <-  function(w_k_aux,risco_a_T, risco_a_C, delta_T,delta_C, X_cure, X_C, beta_cure, beta_C, alpha, theta){(w_k_aux>0)*(w_k_aux<1)}  
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o log-verossimilhan?a utilizada no optim() para estimar os parametros beta_cure e lambdas_j dos tempos de falha T                 ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_param_T <-  function(param_T,delta.t, xi_T,X_cure,num_int_T,bi){   #lambda_T=param_T
  colcura <- ncol(X_cure)
  n <- length(X_cure[,1])
  Delta <- delta.t%*%t(rep(1,num_int_T))
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

#modelo_param_T <-  function(param_T,delta.t, xi_T,X_cure,num_int_T,bi){   #lambda_T=param_T
#  
#  n <- length(X_cure[,1])
#  Delta <- delta.t%*%t(rep(1,num_int_T))
#  lambda <- t((param_T[1:num_int_T])%*%t(rep(1,n)))  #estima os lambdas_j
#  
#  x_cure_b <- X_cure%*%param_T[(num_int_T+1):(num_int_T+2)]   #estima os betas 
#  prop <- exp(x_cure_b)*rowMeans(exp(bi))
#  
#  aux <- apply(exp(lambda)*xi_T,1,sum)
#  soma <- aux%*%t(rep(1,num_int_T))
#  
#  aux2 <- as.vector(prop)*(1- exp(-apply(exp(lambda)*xi_T,1,sum)))
#  aux3 <- aux2%*%t(rep(1,num_int_T))
#  
#  mean_bi <- rowMeans(bi)%*%t(rep(1,num_int_T))
#  
#  U <- -sum(Delta*(as.vector(x_cure_b) + lambda  + mean_bi - soma) - aux3)
#  
#}      
###--------------------------------------------------------------------------------------------------------------------------------------###


###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o log-verossimilhan?a utilizada no multiroot() para estimar os parametros beta_C e alpha tempos de censura informativa C        ###
###--------------------------------------------------------------------------------------------------------------------------------------###
modelo_C_MEP <-  function(beta_C, X_C,delta.c,risco_a_C,bi,n_intMC){
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
  U_C_1 <- colSums((X_C*delta.c - w_kl_beta_C_num))
  
  U_C_alpha <- sum(delta.c*(rowSums(bi[,])/n_intMC) - w_kl_beta_C_alpha_num)
  
  c(U_C_1 = U_C_1,U_C_alpha=U_C_alpha)
}

#modelo_C_MEP <-  function(beta_C, X_C,delta.c,risco_a_C,bi,n_intMC){
#  q <- ncol(X_C)
#  n <- nrow(X_C)
#  
#  if(ncol(t(t(X_C))) == 1){
#    w_kl_beta_C <- risco_a_C*exp((X_C[,]*beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
#    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]*beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
#  }
#  
#  else{
#    w_kl_beta_C <- risco_a_C*exp((X_C[,]%*%beta_C[1:q]))*(rowSums(exp(beta_C[q+1]*bi[,]))/n_intMC)
#    w_kl_beta_C_alpha_num <- risco_a_C*exp(X_C[,]%*%beta_C[1:q])*(rowSums(bi[,]*exp(beta_C[q+1]*bi[,]))/n_intMC)
#  }
#  
#  w_kl_beta_C_num <- cbind(matrix(w_kl_beta_C, nrow = n, ncol = q))*X_C
#  U_C_1 <- colSums((X_C*delta.c - w_kl_beta_C_num))
#  
#  U_C_alpha <- sum(delta.c*(rowSums(bi[,])/n_intMC) - w_kl_beta_C_alpha_num)
#  
#  c(U_C_1 = U_C_1,U_C_alpha=U_C_alpha)
#}
###--------------------------------------------------------------------------------------------------------------------------------------###


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


###--------------------------------------------------------------------------------------------------------------------------------------###
### Fun??o com o c?lculo das derivadas de primeira ordem. Retorna a esperan?a da matriz de drivadas de primeira ordem.                   ###
### C?lculo ? feito como no artigo de Louis et al (1982)                                                                                 ###
###--------------------------------------------------------------------------------------------------------------------------------------###  
Esp.Deriv.Prim.Ordem <-  function( X_cure, X_C,delta_T,delta_C,id_T, beta_cure, beta_C, alpha, theta,nu_ind_j_T,nu_ind_j_C, lambda_T_j,lambda_C_j, deltaij_T,deltaij_C, w_k_grupo, ident){
  
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
Esp.DerivParciais  <-  function( X_cure, X_C,delta_T,delta_C, beta_cure, beta_C, alpha, theta,nu_ind_j_T,nu_ind_j_C, lambda_T_j,lambda_C_j, deltaij_T,deltaij_C, w_k_grupo, ident){
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

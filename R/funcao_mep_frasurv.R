###-------------------------------------------------------------------------------------------------------------###
### Programa para ajustar o modelo com fra??o de cura  e censura dependente                                     ###
### esse modelo ajustado neste programa considera censura ? direita, dependente e distribui??o para os          ###
### tempos de falha e censura dependente como sendo MEP.                                                        ###
### para o modelo de fracao de cura e coniderado modelo teoria unificada                                        ###
### com a abordagem utilizada ? abordagem frequentista                                                         ###
###-------------------------------------------------------------------------------------------------------------###

rm(list=ls(all=TRUE)) 

set.seed(123456)

source("C:/Users/Gabriel-PC/Desktop/progamation/survfracur/R/biblio_MEP_fragCura.R",local=TRUE)
require(dlm)
require(survival)
require(MfUSampler)
require(rootSolve)
library(tidyverse)

# LIBS DEPCENS
library(Formula)
library(survival)
library(stats)
library(rootSolve)
library(dlm)
library(matrixStats)
library(graphics)


dados <- read.table("C:/Users/Gabriel-PC/Desktop/progamation/Bolsa Silvana/simulacoes/data/dados2.txt", head=TRUE)

dados <- dados %>%
  mutate(
    delta_t = ifelse(dados$cens==1,1,0),
    delta_c = ifelse(dados$cens==2,1,0)
  )

delta.t = ifelse(dados$cens==1,1,0)
delta.c = ifelse(dados$cens==2,1,0)

#DADOS FUNÇÃO!!
delta_t = ifelse(dados$cens==1,1,0)
delta_c = ifelse(dados$cens==2,1,0)
delta.t=delta_t
delta.c=delta_c
ident=dados$ident
data = dados
formula = time ~ x1_cure + x_c2 | x_c1 + x_c2 + x1_cure
#DADOS FUNÇÃO!!

#=======================================================================================
########################### Inicio do Estudo Monte Carlo ###############################
#=======================================================================================



model_MEP_dep <-  function(formula, data, delta.t, delta.c, ident, Num_intervals){
  
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
    w_trans <- arms(0.5, myldens=log_cond_wi, indFunc=support_wi,n.sample=n_intMC,risco_a_T=risco_a_T[ident==k], risco_a_C=risco_a_C[ident==k], delta_T=delta.t[ident==k],delta_C=delta.c[ident==k], X_cure=X_cura[ident==k,],X_C=X_C[ident==k,], beta_cure=beta_cura,beta_C=beta_C,alpha=alpha, theta=theta) 
    w_auxi <- log(w_trans)-log(1-w_trans) 
    w_chapeu_grupo[k,] <- w_auxi
  }
  bi = w_chapeu_grupo[ident,]
  
  theta <- mean(w_chapeu_grupo^2)
  #print(theta)
  ###----------------------------------------------------------------------------------------
  a <- time.grid(t,delta.t,bmax); a <- a[[1]]; #grade dos tempos de falha
  c <- time.grid(t,delta.c,bmax); c <- c[[1]]; #grade dos tempos de censura
  b <- length(a)-1
  d <- length(c)-1
  a[length(a)] <- max(t)
  c[length(c)] <- max(t)
  num_int_T <- bmax
  num_int_C <- bmax
  
  id_T <- as.numeric(cut(t,a))  # id sempre ? com grade com inf na ?ltima casela
  nu_T <- tapply(delta.t,id_T ,sum)
  
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
  
  param_T <- optim(par = c(lambda_T_j, beta_cura), fn = modelo_param_T, control = list(maxit = 5000), method = c( "BFGS"), delta.t=delta.t, xi_T=xi.falha,X_cure=X_cura,num_int_T=bmax, bi=bi)
  param_Ts <- param_T$par
  lambda_T_j <- exp(param_Ts[1:num_int_T])
  beta_cura  <- param_Ts[(num_int_T+1):(num_int_T+p)]
  risco_a_T <- apply(t(lambda_T_j*t(xi.falha)),1,sum)
  ###--------------------------------------------------------------------------------
  
  pred_linear_C <- exp(X_C%*%beta_C)*rowMeans(exp(alpha*bi[,]))
  id_C <- as.numeric(cut(t,c))  # id sempre ? com grade com inf na ?ltima casela
  nu_C <- tapply(delta.c,id_C ,sum)
  
  xi.cens_pred <- matrix(0,n,d)
  for(i in 1:n){                        
    for(j in 1:d){
      xi.cens_pred[i,j] <- (min(t[i], c[j+1])-c[j])*((t[i]-c[j])>0)*pred_linear_C[i]
    }
  }
  
  lambda_C_j <- nu_C/apply(xi.cens_pred,2,sum)
  risco_a_C <- apply(t(as.vector(lambda_C_j)*t(xi.cens)),1,sum)

  S_C <- multiroot(f = modelo_C_MEP, start = c(beta_C,alpha), X_C=X_C,delta.c=delta.c,risco_a_C=risco_a_C,bi=bi,n_intMC=n_intMC)
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
Esp_deriv_ordem1 <- Esp.Deriv.Prim.Ordem( X_cure=X_cura, X_C=X_C,delta_T=delta.t,delta_C=delta.c,id_T=id_T, beta_cure=beta_cura, beta_C=beta_C, alpha=alpha,
                                          theta=theta,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j, 
                                          deltaij_T=xi.falha,deltaij_C=xi.cens, w_k_grupo=w_chapeu_grupo, ident=ident)

Esp_deriv_ordem2 <- Esp.DerivParciais( X_cure=X_cura, X_C=X_C,delta_T=delta.t,delta_C=delta.c, beta_cure=beta_cura, beta_C=beta_C, alpha=alpha,
                                       theta=theta,nu_ind_j_T=nu_T,nu_ind_j_C=nu_C, lambda_T_j=lambda_T_j,lambda_C_j=lambda_C_j, 
                                       deltaij_T=xi.falha,deltaij_C=xi.cens, w_k_grupo=w_chapeu_grupo, ident=ident)

InfFisher <- (-Esp_deriv_ordem2) - Esp_deriv_ordem1
Var <- solve(InfFisher)
ErroPadrao <- sqrt(diag(Var)) 

###-----------------------------------------------------------------------------------
# c?lculo dos crit?rios de sele??o de modelos
criterios <- crit_MEP(X_cure=X_cura,X_C,bi,n,beta_cure=beta_cura,lambda_T_j,beta_C,alpha,lambda_C_j,delta_t=delta.t,delta_c=delta.c,risco_a_T,risco_a_C,id_T,id_C,m,ident)

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


fit_frasurv_MEP <- model_MEP_dep(formula = time ~ x1_cure + x_c1 + x_c2 | x_c1 + x_c2, 
                              data = dados,
                              delta.t = dados$delta_t,
                              delta.c = dados$delta_c,
                              ident = dados$ident)


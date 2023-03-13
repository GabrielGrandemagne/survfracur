###-------------------------------------------------------------------------------------------------------------###
### Programa para ajustar o modelo com fra??o de cura  e censura dependente                                     ###
### esse modelo ajustado neste programa considera censura ? direita, dependente e distribui??o para os          ###
### tempos de falha como sendo Weibull, com parametrizacao dada pelo livro de Ibrahim, 159                      ###
### para o modelo de fracao de cura e coniderado modelo teoria unificada                                        ###
### com a abordagem utilizada ? abordagem frequentista                                                          ###
###-------------------------------------------------------------------------------------------------------------###

rm(list=ls(all=TRUE)) 

wd <- getwd() 
set.seed(123456)

#source("biblio_Weibull_fragCura_DOISSSS.R",local=TRUE)
source("C:/Users/Gabriel-PC/Desktop/progamation/survfracur/R/biblio_Weibull_fragCura.R",local=TRUE)
require(dlm)
require(survival)
require(rootSolve)
library(Formula)
library(survival)
library(rootSolve)
library(dlm)
library(matrixStats)
library(stats)
library(graphics)

###-------------------------------------------------------------------------------------------------
### Definir as pastas,
#modelos <- c("data", "Ajuste")  
#
#dir.create("simulacoes")
#wd.sim <- rep(NA,length(modelos))
#wd.mod <- paste(wd, "/modelos/", sep="")
#
#for(i in 1:length(modelos))
#{
#  dir.create( paste(wd, "/simulacoes/", modelos[i], sep="") )
#  wd.sim[i] <- paste(wd, "/simulacoes/", modelos[i], sep="")
#}     
#
#parametros_out <- file.path(wd.sim[2], "parametros_Weibull.txt")
#Erro_Padrao    <- file.path(wd.sim[2], "ErroPadrao_Weibull.txt")
###------------------------------------------------------------------------------------------
# especificacoes dos dados
#n <- 400
#alpha.c <- 2
#alpha.t <- 1
#lambda.c <- exp(-2)
#lambda.t <- exp(0)
#alphaa <- 1
#beta_cure <- c(-0.6, 0.4)
#beta_C <- c(0.1,-0.1)
#m <- 40 
#ident <- rep(c(1:m),10)
###--------------------------------------------------------------------------------
# gerar os dados e salvar (utilizar os mesmos dados para todos modelos
###---------------------------------------------------------------------------------
#replica <- 500 
#
#for(iteracao in 1:replica){ 
#  
#  X_C <- cbind(rnorm(n),rbinom(n,1,0.5))
#  X_cure <- cbind(1, rbinom(n,1,0.5))
#  w_aux <- rnorm(m,0,1)
#  w <- w_aux[ident]
#  
#  dados <- simula.teoria_unif(alpha.c, alpha.t,  lambda.c, lambda.t, alphaa, w,  X_cure, X_C, beta_cure,beta_C, ident)
#  
#	
#	summary.data <- matrix(ncol=7, nrow=1)
#	data.labels <- paste("dados",1:replica, ".txt", sep="")
#	data.local  <- file.path(wd.sim[1], data.labels[iteracao])
#	summary.local <- file.path(wd.sim[1], "summary.data.txt")
#	summary.data[1] <- with(dados, cor(t,c, method="pearson"))
#	summary.data[2] <- with(dados, cor(t,c, method="kendall"))
#	summary.data[3] <- with(dados, cor(t,c, method="spearman"))
#	summary.data[4] <- with(dados, mean(event))
#	summary.data[5:7] <- tabulate(dados$cens)/length(dados$cens)
#	write.table(summary.data, file =summary.local, row.names=FALSE, col.names=FALSE, append=TRUE);  
#	write.table(dados, file = data.local, row.names=FALSE, col.names=TRUE, append=FALSE);
#}

#=======================================================================================
########################### Inicio do Estudo Monte Carlo ###############################
#=======================================================================================
  #data.local  <- file.path(wd.sim[1],data.name)
  #dados <- read.table(data.local, head=TRUE)

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
  param_T <- optim(par = c(alpha_T,lambda_T, beta_cura), fn = modelo_param_T, control = list(maxit = 5000), method = c( "Nelder-Mead"),t=t,delta_T=delta_t, X_cure=X_cura, bi=bi)
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


Esp_deriv_ordem1 <- Esp.Deriv.Prim.Ordem(t,delta_T=delta_t, delta_C=delta_c, X_cure=X_cura, X_C=X_C,
                                         beta_cure=beta_cura, beta_C=beta_C, alpha=alpha, theta=theta,
                                         alpha_T=alpha_T, alpha_C=alpha_C, lambda_T=lambda_T, lambda_C=lambda_C,
                                         w_k_grupo=w_chapeu_grupo, ident=ident)


Esp_deriv_ordem2 <- Esp.DerivParciais(t,delta_T=delta_t, delta_C=delta_c, X_cure=X_cura, X_C=X_C,
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

dados <- read.table("C:/Users/Gabriel-PC/Desktop/progamation/Bolsa Silvana/simulacoes/data/dados2.txt", head=TRUE)

dados <- dados %>%
  mutate(
    delta_t = ifelse(dados$cens==1,1,0),
    delta_c = ifelse(dados$cens==2,1,0)
  )


fit_frasurv2 <- model_Weibull_dep(formula = time ~ x1_cure + x_c2 + x_c1 + x_c1 | x_c1 + x_c2 + x1_cure, 
                                 data = dados,
                                 delta_t = dados$delta_t,
                                 delta_c = dados$delta_c,
                                 ident = dados$ident)

#DADOS FUNÇÃO!!
#delta_t = ifelse(dados$cens==1,1,0)
#delta_c = ifelse(dados$cens==2,1,0)
#delta_t=delta_t
#delta_c=delta_c
#ident=dados$ident
#data = dados
#formula = time ~ x1_cure + x_c2 + x_c1 | x_c1 + x_c2
#DADOS FUNÇÃO!!



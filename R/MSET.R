library(tidyverse)
library(bigrquery)
library(pbmcapply)
library(deconvolveR)
library(reshape2)
library(pbmcapply)
library(REBayes)
library(yarrr)

source('epsest.func.R')

# ###hyperparameters
# number of iterations
iter <- 1:500
# proportion of non-nulls 
p <- 0.05
# Number of hypothesis 
m <- 5000
# Number of historical stages for estimation of non-null proportion 
R <- 5

# utility function that calculates data-driven likelihood ratio statistic in mSPRT for gaussian prior
lambda <- function(tau, obs_subset,n){
  if(n == 1){
    s = obs_subset^2
  }
  else{
    s = apply(obs_subset,2,sum)^2
  }
  result = sqrt(1/(n*tau^2+1))*exp(tau^2*s/2/(n*tau^2+1))
  return(result)
}

# main function for calculating fdr and mdr 
simulate <- function(iter,mu,p,m,max_stages,alpha=0.05, tau = 2, mu_fixed){
  # Generate new observations independently sampling from ground truth sampled from 
  theta <- sample(c(0,1),m, prob = c(1-p,p),replace=TRUE)
  if(mu_fixed == TRUE){
    mu <- mu * rep(1, m)
  }
  if(mu_fixed == FALSE){
    mu <- runif(m, min= 1, max = 3)  #random mu for each coordinate 
  }
  true_mu <-  theta * mu
  obs <- matrix(rep(true_mu,max_stages) + rnorm(m*max_stages), nrow = max_stages, byrow = TRUE)
  
  # Initialization 
  T_or <- matrix(NA,nrow=max_stages,ncol=m)
  T_dd <- matrix(NA,nrow=max_stages,ncol=m)
  #Always valid p-value likelihood ratio 
  LR_or_always = matrix(0, nrow = max_stages, ncol = m)
  LR_dd_always = matrix(0, nrow = max_stages, ncol = m)
  
  ## Calculate oracle test statistic in term of probability ratio 
  f0 <- dnorm(obs, mean = 0, sd = 1)
  f1 <- matrix(NA, nrow = max_stages, ncol = m)
  for(stage in 1:max_stages){
    f1[stage,] = dnorm(obs[stage,], mean = mu, sd = 1)
  }
  T_or <- (1-p)*apply(f0, 2, cumprod)/ ((1-p)*apply(f0, 2, cumprod) + p*apply(f1, 2, cumprod))
  
  #NPMLE in REBayes estimation of alternative density
  p_Est_list = c()
  for(r in 1:R){
    obs_past <- true_mu + rnorm(m)
    p_Est_list <- c(p_Est_list,epsest.func(obs_past,0,1)) 
  }
  p_Est = mean(p_Est_list)
  obs_past <- true_mu + rnorm(m)
  gl_res <- GLmix(obs_past, sigma = 1/sqrt(R))
  # p_Est <- epsest.func(obs_past,0,1)  # Use historical data to estimate proportion of non-nulls
  gl_summary <- tibble(x=gl_res$x,y=gl_res$y) %>% 
    mutate(x=round(x,2)) %>% group_by(x) %>% 
    summarise(y=sum(y))%>% round(2) %>% filter(y!=0) %>% mutate(y=y/sum(y)) %>% arrange(abs(x))
  thresh = gl_summary$x[which(cumsum(gl_summary$y) > 1-p_Est)[1]]  # data driven threshold for boundary between null and alternative
  # #########
  # non data-driven method if a threshold is determined ex-ante, possibly from minimum detectable threshold
  thresh = 1
  p_Est <- gl_summary %>% filter(abs(x)>thresh) %>% summarise(non_null=sum(y)) %>% .$non_null
  # #########
  non_null_mus <- gl_summary %>% filter(abs(x)>=abs(thresh)) %>% .$x
  non_null_props <- gl_summary %>% filter(abs(x)>=abs(thresh)) %>% mutate(y=y/sum(y)) %>% .$y
  f1_dd <- matrix(0,nrow=max_stages,ncol=m)
  if(length(non_null_mus) > 0){
    for (i in 1:length(non_null_mus)){
      f1_dd <- f1_dd+non_null_props[i]*dnorm(obs,non_null_mus[i],1)
    }
  }
  T_dd <- (1-p_Est)*apply(f0,2,cumprod)/((1-p_Est)*apply(f0,2,cumprod)+p_Est*apply(f1_dd,2,cumprod))
  
  for(stage in 1:max_stages){
    ## oracle procedure using always valid p-value and BH
    LR_or_always[stage,] = f1[stage,]/f0[stage,]
    ## data driven procedure using always valid p-value and BH
    LR_dd_always[stage,] = lambda(tau, obs[1:stage,],stage)
  }
  
  dec_or <- rep(0,m)
  sorted_T_or_subset <- sort(T_or[max_stages,], index.return = TRUE)
  dec_or[sorted_T_or_subset$ix[which(cummean(sorted_T_or_subset$x)<=alpha)]] <- 1
  dec_dd <- rep(0,m)
  sorted_T_dd_subset <- sort(T_dd[max_stages,], index.return = TRUE)
  dec_dd[sorted_T_dd_subset$ix[which(cummean(sorted_T_dd_subset$x)<=alpha)]] <- 1
  
  dec_or_always = rep(0,m)
  LR_active = apply(LR_or_always, 2, cumprod)
  pval_always = 1/apply(LR_active,2, max)
  dec_or_always[which(p.adjust(pval_always,method='BH') <= alpha)] <- 1
  dec_dd_always = rep(0,m)
  pval_always_dd = 1/apply(LR_dd_always,2,max)
  dec_dd_always[which(p.adjust(pval_always_dd,method='BH') <= alpha)] <- 1
  
  fpr_or <- sum((1-theta)*dec_or)/max(sum(dec_or),1)
  mdr_or <- sum(theta*(1-dec_or))/max(sum(theta),1)
  fpr_dd <- sum((1-theta)*dec_dd)/max(sum(dec_dd),1)
  mdr_dd <- sum(theta*(1-dec_dd))/max(sum(theta),1)
  fpr_or_always <- sum((1-theta)*dec_or_always)/max(sum(dec_or_always),1)
  mdr_or_always <- sum(theta*(1-dec_or_always))/max(sum(theta),1)
  fpr_dd_always <- sum((1-theta)*dec_dd_always)/max(sum(dec_dd_always),1)
  mdr_dd_always <- sum(theta*(1-dec_dd_always))/max(sum(theta),1)
  
  if(mu_fixed == TRUE){
    mu = mean(mu)
    return(tibble(iter,m,p,mu,max_stages,alpha,tau,fpr_or,mdr_or,fpr_dd,mdr_dd,fpr_or_always,mdr_or_always,fpr_dd_always,mdr_dd_always))
  }
  if(mu_fixed == FALSE){
    return(tibble(iter,m,p,max_stages,alpha,tau,fpr_or,mdr_or,fpr_dd,mdr_dd,fpr_or_always,mdr_or_always,fpr_dd_always,mdr_dd_always))
  }
}

#Setting 1: constant signal proportion=0.05, vary mu from 1 to 3
mu <- seq(1,3,0.2)
grid <- expand.grid(iter, mu)
res <- pbmcmapply(simulate,
                  iter=grid[,1],
                  mu=grid[,2],
                  p = 0.05,
                  MoreArgs=list(m=m,max_stages=5, tau = 2, mu_fixed = TRUE),
                  mc.cores=90)
res <- as_tibble(t(res)) %>% unnest
summary_res <- res %>% group_by(m,p,mu,max_stages,alpha,tau) %>%
  summarise(fdr_or=mean(fpr_or),mdr_or=mean(mdr_or),
            fdr_dd=mean(fpr_dd),mdr_dd=mean(mdr_dd),
            fdr_or_always=mean(fpr_or_always), mdr_or_always=mean(mdr_or_always),
            fdr_dd_always=mean(fpr_dd_always), mdr_dd_always=mean(mdr_dd_always))
write_csv(summary_res, 'setting1.csv')


#Setting 2: Vary signal proportion from 0.01 to 0.1 with stepsize 0.01, constant mu at 2
p <- seq(0.01, 0.1, 0.01)
grid <- expand.grid(iter, p)
res <- pbmcmapply(simulate,
                  iter=grid[,1],
                  mu=2,
                  p = grid[,2],
                  MoreArgs=list(m=m,max_stages=5, tau = 2, mu_fixed = TRUE),
                  mc.cores=90)
res <- as_tibble(t(res)) %>% unnest
summary_res <- res %>% group_by(m,p,max_stages,alpha,tau) %>%
  summarise(fdr_or=mean(fpr_or),mdr_or=mean(mdr_or),
            fdr_dd=mean(fpr_dd),mdr_dd=mean(mdr_dd),
            fdr_or_always=mean(fpr_or_always), mdr_or_always=mean(mdr_or_always),
            fdr_dd_always=mean(fpr_dd_always), mdr_dd_always=mean(mdr_dd_always),
            total_obs_or=mean(total_obs_or),total_obs_dd=mean(total_obs_dd),total_obs_or_always=mean(total_obs_or_always),total_obs_dd_always=mean(total_obs_dd_always))
write_csv(summary_res, 'setting2.csv')

# Setting 3: Vary signal proportion from 0.01 to 0.1, random mu from 1 to 3 for each coordinate
p <- seq(0.01, 0.1, 0.01)
grid <- expand.grid(iter, p)
res <- pbmcmapply(simulate,
                  iter=grid[,1],
                  mu=2,
                  p = grid[,2],
                  MoreArgs=list(m=m,max_stages=5, tau = 2, mu_fixed = FALSE),
                  mc.cores=90)
res <- as_tibble(t(res)) %>% unnest
summary_res <- res %>% group_by(m,p,max_stages,alpha,tau) %>%
  summarise(fdr_or=mean(fpr_or),mdr_or=mean(mdr_or),
            fdr_dd=mean(fpr_dd),mdr_dd=mean(mdr_dd),
            fdr_or_always=mean(fpr_or_always), mdr_or_always=mean(mdr_or_always),
            fdr_dd_always=mean(fpr_dd_always), mdr_dd_always=mean(mdr_dd_always),
            total_obs_or=mean(total_obs_or),total_obs_dd=mean(total_obs_dd),total_obs_or_always=mean(total_obs_or_always),total_obs_dd_always=mean(total_obs_dd_always))
write_csv(summary_res, 'setting3.csv')
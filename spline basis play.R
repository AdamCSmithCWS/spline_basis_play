### comparing different basis function calculations


library(rstanarm)
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)
library(tidyverse)
library(GGally)
library(shinystan)
source("functions/GAM_basis_function.R")
library(tidybayes)
library(loo)
library(mgcv)

y = -20:20
y2 = sample(y,size = length(y))
#i = y*0.04 + (y^2)*-0.01 + y2*-0.04 + (y2^2)*0.01 + rnorm(41,0,0.01)
i = y*0.04 + (y^2)*-0.01 + rnorm(41,0,0.01)


dat = data.frame(y = y,i = i)
#ks = seq(-20,20,length = 7)
#M = stan_gamm4(i ~ s(y2,k = 7) + s(y,k = 7),data = dat)
M = stan_gamm4(i ~ 1 + s(y,k = 8, bs = "tp"),data = dat,
               chains=1, iter=200,
               warmup=100,
               cores = 1)

Mj = jagam(i ~ 1 + s(y,k = 8, bs = "tp"),data = dat,
           file = "temp_jag.del")
head(M$x)
rst_X = data.frame(matrix(M$x[,-c(1)],nrow = 41,ncol = 7))
names(rst_X) <- paste("rst",1:7,sep = "_")

head(Mj$jags.data$X)
j_X = data.frame(matrix(Mj$jags.data$X[,-c(1)],nrow = 41,ncol = 7))
names(j_X) <- paste("mj",1:7,sep = "_")



# plot_nonlinear(M)
# rstan::get_stanmodel(M$stanfit)
# head(M$x)


Mbr = brm(bf(i ~ 1 + s(y,k = 8, bs = "tp")),data = dat,
          chains=1, iter=200,
          warmup=100,
          cores = 1)

#stanmod = stancode(Mbr)

st_dat = standata(Mbr)
#st_dat$Xs[,1] #linear predictor


brm_X = data.frame(st_dat$Zs_1_1)
names(brm_X) <- paste("brm",1:6,sep = "_")
brm_X[,"brm_7"] <- as.numeric(st_dat$Xs)


paired_mat = bind_cols(rst_X,j_X,brm_X)

pp = ggpairs(data = paired_mat)

print(pp)




hb = gam.basis.func(orig.preds = dat$y,nknots = 7,
                    sm_name = "home_brew",
                    standardize = "z")
hb_X = data.frame(hb[["home_brew_basis"]])
names(hb_X) <- paste("hb",1:7)
# paired_mat2 = bind_cols(brm_X,hb_X)
# 
# pp2 = ggpairs(data = paired_mat2)
# 
# print(pp2)


dat_brm <- list(y = y,
               x = i,
               Xb = as.matrix(brm_X),
               N = 41,
               nknots = 7)


dat_hb <- list(y = y,
                x = i,
                Xb = hb[["home_brew_basis"]],
               N = 41,
               nknots = 7)


dat_rst <- list(y = y,
               x = i,
               Xb = as.matrix(rst_X),
               N = 41,
               nknots = 7)

dat_j <- list(y = y,
                x = i,
                Xb = as.matrix(j_X),
                N = 41,
                nknots = 7)

modl_file = "models/simple_GAM.Stan"


modl = stan_model(file=modl_file)

parms = c("betas",
          "mu",
          "sd_beta",
          "sigma",
          "alpha",
          "log_lik")
## run sampler on model, data
stanfit_rst <- sampling(modl,
                               data=dat_rst,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=2000,
                               warmup=1000,
                               cores = 4,
                               pars = c(parms),
                               control = list(adapt_delta = 0.95,
                                              max_treedepth = 15))

stanfit_j <- sampling(modl,
                        data=dat_j,
                        verbose=TRUE, refresh=100,
                        chains=4, iter=2000,
                        warmup=1000,
                        cores = 4,
                        pars = c(parms),
                        control = list(adapt_delta = 0.95,
                                       max_treedepth = 15))


## run sampler on model, data
stanfit_hb <- sampling(modl,
                         data=dat_hb,
                         verbose=TRUE, refresh=100,
                         chains=4, iter=2000,
                         warmup=1000,
                         cores = 4,
                         pars = c(parms),
                         control = list(adapt_delta = 0.95,
                                        max_treedepth = 15))

stanfit_brm <- sampling(modl,
                       data=dat_brm,
                       verbose=TRUE, refresh=100,
                       chains=4, iter=2000,
                       warmup=1000,
                       cores = 4,
                       pars = c(parms),
                       control = list(adapt_delta = 0.95,
                                      max_treedepth = 15))



loo_rst = loo(stanfit_rst)
loo_hb = loo(stanfit_hb)
loo_j = loo(stanfit_j)
loo_brm = loo(stanfit_brm)


#launch_shinystan(stanfit_rst)

mu_hb <- gather_draws(stanfit_hb,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "hb")%>% 
  bind_cols(.,dat)



mu_brm <- gather_draws(stanfit_brm,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "brm")%>% 
  bind_cols(.,dat)



mu_rst <- gather_draws(stanfit_rst,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "rst")%>% 
  bind_cols(.,dat)


mu_j <- gather_draws(stanfit_j,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "j")%>% 
  bind_cols(.,dat)


mu = bind_rows(mu_hb,mu_rst,mu_j,mu_brm)



x11()
mupl = ggplot(data = mu,aes(x = n,y = mu))+
  geom_line(aes(colour = model))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),alpha = 0.2)+
  geom_point(aes(x = n,y = i),alpha = 0.3)
print(mupl)


loo_rst
loo_hb 
loo_j 
loo_brm 



get_elapsed_time(stanfit_hb)
get_elapsed_time(stanfit_rst)




b_hb <- gather_draws(stanfit_hb,betas[n]) %>% 
  group_by(n) %>% 
  summarise(b = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "hb") 

sdb_hb <- gather_draws(stanfit_hb,sd_beta) %>% 
  summarise(sdb = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "hb") 



b_rst <- gather_draws(stanfit_rst,betas[n]) %>% 
  group_by(n) %>% 
  summarise(b = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "rst")

sdb_rst <- gather_draws(stanfit_rst,sd_beta) %>% 
  summarise(sdb = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "rst") 


b = bind_rows(b_hb,b_rst)
sdb = bind_rows(sdb_hb,sdb_rst)
sdb
loo_rst
loo_hb



x11()
bpl = ggplot(data = b,aes(x = n,y = b))+
  geom_line(aes(colour = model))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),alpha = 0.2)
print(bpl)


launch_shinystan(stanfit_hb)
launch_shinystan(stanfit_rst)





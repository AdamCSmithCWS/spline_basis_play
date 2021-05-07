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


y = -20:20
y2 = sample(y,size = length(y))
#i = y*0.04 + (y^2)*-0.01 + y2*-0.04 + (y2^2)*0.01 + rnorm(41,0,0.01)
i = y*0.04 + (y^2)*-0.01 + rnorm(41,0,0.01)


dat = data.frame(y = y,i = i)
#ks = seq(-20,20,length = 7)
#M = stan_gamm4(i ~ s(y2,k = 7) + s(y,k = 7),data = dat)
M = stan_gamm4(i ~ 1 + s(y,k = 7),data = dat)


# plot_nonlinear(M)
# rstan::get_stanmodel(M$stanfit)
# head(M$x)
# 
# 
# Mbr = brm(bf(i ~ 1 + s(y,k = 7)),data = dat)
# 
# stanmod = stancode(Mbr)
# st_dat = standata(Mbr)
# st_dat$Xs[,1]
# plot(st_dat$Xs[,1],M$x[,7])
# abline(0,1)
# 
# brm_X = data.frame(st_dat$Zs_1_1)
# names(brm_X) <- paste("brm",1:5,sep = "_")
# brm_X[,"brm_6"] <- as.numeric(st_dat$Xs)
rst_X = data.frame(matrix(M$x[,-1],nrow = 41,ncol = 6))
names(rst_X) <- paste("rst",1:6,sep = "_")

# paired_mat = bind_cols(brm_X,rst_X)
# 
# pp = ggpairs(data = paired_mat)
# 
# print(pp)
# 
# ### they're identical, just ordered differently


hb = gam.basis.func(orig.preds = dat$y,nknots = 6,
                    sm_name = "home_brew",
                    standardize = "z")
hb_X = data.frame(hb[["home_brew_basis"]])
names(hb_X) <- paste("hb",1:6)
# paired_mat2 = bind_cols(brm_X,hb_X)
# 
# pp2 = ggpairs(data = paired_mat2)
# 
# print(pp2)


dat_hb <- list(y = y,
                x = i,
                Xb = hb[["home_brew_basis"]],
               N = 41,
               nknots = 6)


dat_rst <- list(y = y,
               x = i,
               Xb = as.matrix(rst_X),
               N = 41,
               nknots = 6)

modl_file = "models/simple_GAM.Stan"


modl = stan_model(file=modl_file)

parms = c("betas",
          "mu",
          "sd_beta",
          "sigma",
          "alpha")
## run sampler on model, data
stanfit_rst <- sampling(modl,
                               data=dat_rst,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=2000,
                               warmup=1000,
                               cores = 4,
                               pars = c(parms),
                               control = list(adapt_delta = 0.9,
                                              max_treedepth = 15))



## run sampler on model, data
stanfit_hb <- sampling(modl,
                         data=dat_hb,
                         verbose=TRUE, refresh=100,
                         chains=4, iter=2000,
                         warmup=1000,
                         cores = 4,
                         pars = c(parms),
                         control = list(adapt_delta = 0.9,
                                        max_treedepth = 15))




#launch_shinystan(stanfit_rst)
launch_shinystan(stanfit_hb)
launch_shinystan(stanfit_rst)

mu_hb <- gather_draws(stanfit_hb,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "hb")




mu_rst <- gather_draws(stanfit_rst,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "rst")

mu = bind_rows(mu_hb,mu_rst)



x11()
mupl = ggplot(data = mu,aes(x = n,y = mu))+
  geom_line(aes(colour = model))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),alpha = 0.2)
print(mupl)








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

x11()
bpl = ggplot(data = b,aes(x = n,y = b))+
  geom_line(aes(colour = model))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),alpha = 0.2)
print(bpl)







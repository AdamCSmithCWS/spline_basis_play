### comparing different basis function calculations


library(rstan)
library(rstanarm)
library(brms)
library(tidyverse)
library(GGally)
source("functions/GAM_basis_function.R")

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
head(M$x)


Mbr = brm(bf(i ~ 1 + s(y,k = 7)),data = dat)

stanmod = stancode(Mbr)
st_dat = standata(Mbr)
st_dat$Xs[,1]
plot(st_dat$Xs[,1],M$x[,7])
abline(0,1)

brm_X = data.frame(st_dat$Zs_1_1)
names(brm_X) <- paste("brm",1:5,sep = "_")
brm_X[,"brm_6"] <- as.numeric(st_dat$Xs)
rst_X = data.frame(matrix(M$x[,-1],nrow = 41,ncol = 6))
names(rst_X) <- paste("rst",1:6,sep = "_")

paired_mat = bind_cols(brm_X,rst_X)

pp = ggpairs(data = paired_mat)

print(pp)

### they're identical, just ordered differently


hb = gam.basis.func(orig.preds = dat$y,nknots = 6,
                    sm_name = "home_brew",
                    standardize = "z")
hb_X = data.frame(hb[["home_brew_basis"]])
names(hb_X) <- paste("hb",1:6)
paired_mat2 = bind_cols(brm_X,hb_X)

pp2 = ggpairs(data = paired_mat2)

print(pp2)


dat_hb <- list(y = y,
                x = i,
                Xb = hb[["home_brew_basis"]],
               N = 41,
               nknots = 6)


dat_brm <- list(y = y,
               x = i,
               Xb = as.matrix(brm_X),
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
stanfit_brms <- sampling(modl,
                               data=dat_brm,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=2000,
                               warmup=1000,
                               cores = 4,
                               pars = c(parms),
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 10))



## run sampler on model, data
stanfit_hb <- sampling(modl,
                         data=dat_hb,
                         verbose=TRUE, refresh=100,
                         chains=4, iter=2000,
                         warmup=1000,
                         cores = 4,
                         pars = c(parms),
                         control = list(adapt_delta = 0.8,
                                        max_treedepth = 10))







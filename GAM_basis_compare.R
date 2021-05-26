### comparing different basis function calculations


library(rstanarm)
library(brms)
library(rstan)
rstan_options(auto_write = TRUE)
library(tidyverse)
library(GGally)
library(shinystan)
library(tidybayes)
library(loo)
library(mgcv)


#super simple quadratic true relationship
x = -20:20
y = 2 + x*0.04 + (x^2)*-0.01 + rnorm(length(x),0,0.01)
dat = data.frame(y = y,x = x)


#  rstanarm basis ---------------------------------------------------------

## just to get the basis function
## minimal (and obviously insufficient) fit
## Is there a way to get the basis without fitting it?
M = stan_gamm4(y ~ s(x,k = 8, bs = "cr"),data = dat,
               chains=1, iter=600,
               warmup=300,
               cores = 1)
# drop the constant term re-order columns to match brms (not necessary, but tidy),
rstanarm_basis = matrix(M$x[,c(2,7,6,5,4,3,8)],nrow = 41,ncol = 7)

modl = get_stancode(M$stanfit)
cat(modl,file = "tmp_stancode.stan")


#  brms basis ------------------------------------------------------------
## just to get the basis function
## minimal (and obviously insufficient) fit
## Is there a way to get the basis without fitting it?

brms_code = make_stancode(bf(y ~ s(x,k = 8, bs = "ps")),data = dat)
Mbr = brm(bf(y ~ s(x,k = 8, bs = "ds")),data = dat,
          chains=1, iter=200,
          warmup=100,
          cores = 1)

st_dat = standata(Mbr)
## brms separates the linear term from the rest of the basis, this just adds it back to match rstanarm
brms_basis <- matrix(c(st_dat$Zs_1_1,st_dat$Xs),nrow = 41,ncol = 7)
brms_code2 <- stancode(Mbr)
# they are identical
all(rstanarm_basis == brms_basis)



#  jagam basis ------------------------------------------------------------

Mj = jagam(y ~ s(x,k = 8, bs = "tp"),data = dat,
           file = "temp_jag.del")
head(Mj$jags.data$X)

# drop the constant term 
jagam_basis = matrix(Mj$jags.data$X[,-c(1)],nrow = 41,ncol = 7)




# mgcv smooth con ---------------------------------------------------------

## specific arguments here generate the same smooth as rstanarm
sm = smoothCon(s(x,k = 8, bs = "tp"),data = dat,
               absorb.cons=TRUE,#this drops the constant
               diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 


mgcv_basis = sm[[1]]$X#[,-7]

# gamm4 basis -------------------------------------------------------------

#  g4 = gamm4(y ~ s(x,k = 8, bs = "ps"),data = dat)
# g4_basis = matrix(g4$gam$X[,-c(1)],nrow = 41,ncol = 7)
# # 


# convert rstan basis to dataframe for ggpairs plot below
rstanarm_X = data.frame(rstanarm_basis) 
rstanarm_X$pkg = "rstanarm"
# convert jagam basis to dataframe for ggpairs plot below
jagam_X = data.frame(jagam_basis)
jagam_X$pkg = "jagam"

# convert jagam basis to dataframe for ggpairs plot below
mgcv_X = data.frame(mgcv_basis)
mgcv_X$pkg = "mgcv"


# paired_mat = bind_cols(rstanarm_X,jagam_X)
# 
# pp = ggpairs(data = paired_mat)
# 
# png(filename = "pairs_plot.png",
#     res = 150,
#     width = 480*2,
#     height = 480*2)
# print(pp)
# dev.off()


## stack paired mat for plotting the basis
stacked_basis = bind_rows(rstanarm_X,jagam_X,mgcv_X) %>% 
  mutate(n = rep(1:41,3)) %>% 
  pivot_longer(.,
               cols = starts_with("X"),
               names_to = "knot",
               values_to = "basis")

bas_gg = ggplot(data = stacked_basis,aes(x = n,y = basis,colour = knot))+
  geom_line()+
  facet_wrap(~pkg,ncol = 3)
print(bas_gg)




dat_mgcv <- list(y = y,
               Xb = mgcv_basis,
               N = 41,
               nknots = 7)


dat_rstanarm <- list(y = y,
               Xb = rstanarm_basis,
               N = 41,
               nknots = 7)

dat_jagam <- list(y = y,
                Xb = jagam_basis,
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
stanfit_rstanarm <- sampling(modl,
                               data=dat_rstanarm,
                               verbose=TRUE, refresh=100,
                               chains=4, iter=2000,
                               warmup=1000,
                               cores = 4,
                               pars = c(parms),
                               control = list(adapt_delta = 0.95,
                                              max_treedepth = 15))

stanfit_jagam <- sampling(modl,
                        data=dat_jagam,
                        verbose=TRUE, refresh=100,
                        chains=4, iter=2000,
                        warmup=1000,
                        cores = 4,
                        pars = c(parms),
                        control = list(adapt_delta = 0.95,
                                       max_treedepth = 15))


stanfit_mgcv <- sampling(modl,
                       data=dat_mgcv,
                       verbose=TRUE, refresh=100,
                       chains=4, iter=2000,
                       warmup=1000,
                       cores = 4,
                       pars = c(parms),
                       control = list(adapt_delta = 0.95,
                                      max_treedepth = 15))



loo_rstanarm = loo(stanfit_rstanarm)
loo_jagam = loo(stanfit_jagam)
loo_mgcv = loo(stanfit_mgcv)


loo_rstanarm
loo_mgcv 
loo_jagam




mu_mgcv <- gather_draws(stanfit_mgcv,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "mgcv",
         n = n-21)%>% 
  bind_cols(.,dat)



mu_rstanarm <- gather_draws(stanfit_rstanarm,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "rstanarm",
         n = n-21)%>% 
  bind_cols(.,dat)


mu_jagam <- gather_draws(stanfit_jagam,mu[n]) %>% 
  group_by(n) %>% 
  summarise(mu = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "jagam",
         n = n-21)%>% 
  bind_cols(.,dat)


mu = bind_rows(mu_rstanarm,mu_jagam,mu_mgcv)



x11()
mupl = ggplot(data = mu,aes(x = n,y = mu))+
  geom_line(aes(colour = model))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),alpha = 0.2)+
  geom_point(data = dat,aes(x = x,y = y),alpha = 0.3,inherit.aes = FALSE)
print(mupl)

png(filename = "predictions_plot.png",
    res = 150,
    width = 480*2,
    height = 480*2)
print(mupl)
dev.off()




get_elapsed_time(stanfit_jagam)
get_elapsed_time(stanfit_rstanarm)
get_elapsed_time(stanfit_mgcv)




# b_hb <- gather_draws(stanfit_hb,betas[n]) %>% 
#   group_by(n) %>% 
#   summarise(b = mean(.value),
#             lci = quantile(.value,0.025),
#             uci = quantile(.value,0.975)) %>% 
#   mutate(model = "hb") 
# 
# sdb_hb <- gather_draws(stanfit_hb,sd_beta) %>% 
#   summarise(sdb = mean(.value),
#             lci = quantile(.value,0.025),
#             uci = quantile(.value,0.975)) %>% 
#   mutate(model = "hb") 



b_rstanarm <- gather_draws(stanfit_rstanarm,betas[n]) %>% 
  group_by(n) %>% 
  summarise(b = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "rstanarm")

b_sum_rstanarm <- gather_draws(stanfit_rstanarm,betas[n]) %>% 
  group_by(.draw) %>% 
  summarise(bb = sum(.value)) %>% 
  ungroup() %>% 
  summarise(b_sum = mean(bb),
            b_sumlci = quantile(bb,0.025),
            b_sumuci = quantile(bb,0.975)) %>% 
  mutate(model = "rstanarm")


sdb_rstanarm <- gather_draws(stanfit_rstanarm,sd_beta) %>% 
  summarise(sdb = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "rstanarm") 



b_jagam <- gather_draws(stanfit_jagam,betas[n]) %>% 
  group_by(n) %>% 
  summarise(b = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "jagam")

sdb_jagam <- gather_draws(stanfit_jagam,sd_beta) %>% 
  summarise(sdb = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(model = "jagam") 



b = bind_rows(b_jagam,b_rstanarm)
sdb = bind_rows(sdb_jagam,sdb_rstanarm)
sdb


x11()
bpl = ggplot(data = b,aes(x = n,y = b))+
  geom_line(aes(colour = model))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = model),alpha = 0.2)
print(bpl)


launch_shinystan(stanfit_hb)
launch_shinystan(stanfit_rstanarm)





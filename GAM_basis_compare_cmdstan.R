### comparing different basis function calculations


library(cmdstanr)
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





# mgcv smooth con ---------------------------------------------------------

## specific arguments here generate the same smooth as rstanarm
sm = smoothCon(s(x,k = 8, bs = "tp"),data = dat,
               absorb.cons=TRUE,#this drops the constant
               diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 


mgcv_basis = sm[[1]]$X#[,-7]



dat_mgcv <- list(y = y,
               Xb = mgcv_basis,
               N = 41,
               nknots = 7)



modl_file = paste0(getwd(),"/models/simple_GAM.Stan")


modl = cmdstan_model(modl_file)

parms = c("betas",
          "mu",
          "sd_beta",
          "sigma",
          "alpha",
          "log_lik")
## run sampler on model, data

stanfit_mgcv <- modl$sample(data=dat_mgcv,
                            refresh=100,
                       chains=4, iter_sampling =2000,
                       iter_warmup=1000,
                       parallel_chains = 4,
                       max_treedepth = 15,
                       adapt_delta = 0.95)




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





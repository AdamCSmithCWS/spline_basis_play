### generating trajectories with known trends, fluctuations, and change-points


library(tidyverse)
source("Functions/gam_basis_function_mgcv.R")
library(rstan)
rstan_options(auto_write = TRUE)
library(shinystan)
library(tidybayes)


years = 1966:2019
nyears = length(years)
routes = 1:10
nroutes = 10
mean_count = log(5)
sd_route = 0.1
sd_ann_fluct = log(1.15) # average annual fluctuation = 15% or -13%
sd_ann_fluct_route = sd_ann_fluct*0.05 #5% variation among routes
t = c(0,0.03,-0.03)
cp = c(1985,2001)
sdt = rep(0.002,3)

betas_alphas = data.frame(route = routes,
                   alpha = rnorm(nroutes,mean_count,sd_route),
                   beta1 = rnorm(nroutes,t[1],sdt[1]),
                   beta2 = rnorm(nroutes,t[2],sdt[2]),
                   beta3 = rnorm(nroutes,t[3],sdt[3]))
year_flucts = rep(0,nyears)
year_flucts[2:(nyears-1)] = rep(c(-1,1),length.out = nyears-2)*abs(rnorm(nyears-2,0,sd_ann_fluct))
## above, ensures year effects in first and last year = 0

count_df = data.frame(year = rep(years[1],nroutes),
                      route = routes,
                      elambda = betas_alphas$alpha)
count_df$lambda <- exp(count_df$elambda)

count_df$count <- rpois(nroutes,count_df$lambda)

for(y in 2:nyears){
  yn = years[y]
  tmp = data.frame(year = rep(yn,nroutes),
                   route = routes,
                   elambda = NA)
  if(yn < cp[1]){bb = 1}
  if(yn >= cp[1] & yn < cp[2]){bb = 2}
  if(yn >= cp[2]){bb = 3}
  yeff = year_flucts[y]
for(rt in routes){
  
  lyr_elambda = count_df[which(count_df$year == (yn-1) &
                          count_df$route == rt),"elambda"]
    tmp[rt,"elambda"] <- lyr_elambda + (t[bb]*lyr_elambda) 
    tmp[rt,"lambda"] <- exp(tmp[rt,"elambda"]+ rnorm(1,yeff,sd_ann_fluct_route))
    tmp[rt,"count"] <- rpois(1,tmp[rt,"lambda"])
  
}#rt
  
  count_df <- bind_rows(count_df,tmp)
}#y


count_df$yr <- count_df$year-(min(count_df$year)-1)


gam_dat = gam_basis(orig.preds = years,
                    nknots = 25,
                    sm_name = "year")


stan_data = list(count = count_df$count,
                 year = count_df$yr,
                 route = count_df$route,
                 nyears = nyears,
                 ncounts = nrow(count_df),
                 nroutes = nroutes,
                 nknots_year = gam_dat$nknots_year,
                 year_basispred = gam_dat$year_basis)



mod.file = "models/simple_GAM_route.stan"

## compile model
slope_model = stan_model(file=mod.file)

## run sampler on model, data
slope_stanfit <- sampling(slope_model,
                               data=stan_data,
                               verbose=TRUE, refresh=100,
                               chains=3,
                          iter=2000,
                               warmup=1000,
                          cores = 3,
                               #pars = parms,
                               control = list(adapt_delta = 0.8,
                                              max_treedepth = 14))



#launch_shinystan(slope_stanfit)
source("Functions/posterior_summary_functions.R")


years_df = unique(count_df[,c("yr","year")]) %>% 
  rename(y = yr)



sdBETA_samples = posterior_samples(fit = slope_stanfit,
                     parm = "sdBETA")
BETA_samples = posterior_samples(fit = slope_stanfit,
                                 parm = "BETA",
                                 dims = "k")

nsmooth_samples = posterior_samples(fit = slope_stanfit,
                                    parm = "nsmooth",
                                    dims = "y")
sdyear_samples = posterior_samples(fit = slope_stanfit,
                                   parm = "sdyear")

n_samples = posterior_samples(fit = slope_stanfit,
                                    parm = "n",
                                    dims = "y")


yeareff_samples = posterior_samples(fit = slope_stanfit,
                              parm = "yeareffect",
                              dims = "y")

yeareffect = posterior_sums(yeareff_samples,dims = "y") %>% 
  mutate(tr_y = year_flucts)



sdBETA = posterior_sums(sdBETA_samples)
sdyear = posterior_sums(sdyear_samples)
nsmooth = posterior_sums(nsmooth_samples,
                         dims = "y") %>% 
  left_join(.,years_df,by = "y")


BETA = posterior_sums(BETA_samples,
                         dims = "k")


ns = posterior_sums(n_samples,
                         dims = "y") %>% 
  left_join(.,years_df,by = "y")

# plotting of route-level trajectories and estimates ----------------------


overplot = ggplot(data = count_df,aes(x = year,y = lambda,colour = route,group = route))+
  geom_line()+
  geom_line(data = ns,size = 2,aes(x = year,y = median),colour = "darkorange",alpha = 0.8,inherit.aes = FALSE)+
  geom_ribbon(data = nsmooth,aes(x = year,ymin = lci,ymax = uci),alpha = 0.2,inherit.aes = FALSE)+
  geom_line(data = nsmooth,size = 2,aes(x = year,y = median),alpha = 0.8,inherit.aes = FALSE)

print(overplot)




beta_plot <- ggplot(data = BETA,aes(x = k,y = mean))+
  geom_errorbar(aes(ymin = lci,ymax = uci),width = 0,alpha = 0.2)+
  geom_point()

print(beta_plot)

sdBETA

year_plot <- ggplot(data = yeareffect,aes(x = y,y = mean))+
  geom_errorbar(mapping = aes(ymin = lci,ymax = uci),width = 0,alpha = 0.4)+
  geom_point()+
  geom_point(data = yeareffect,mapping = aes(x = y,y = tr_y),colour = "darkorange",alpha = 0.8,inherit.aes = FALSE)+
  coord_flip()

print(year_plot)







# explroe basis function --------------------------------------------------

basis = as.data.frame(stan_data$year_basispred)


stacked_basis = basis %>% 
  mutate(n = 1:50) %>% 
  pivot_longer(.,
               cols = starts_with("V"),
               names_to = "knot",
               values_to = "basis")

bas_gg = ggplot(data = stacked_basis,aes(x = n,y = basis,colour = knot))+
  geom_line()+
  facet_wrap(~pkg,ncol = 3)
print(bas_gg)

























































### generating trajectories with known trends, fluctuations, and change-points


library(tidyverse)
source("Functions/gam_basis_function_mgcv.R")
# library(rstan)
# rstan_options(auto_write = TRUE)
library(shinystan)
#library(tidybayes)
library(cmdstanr)

years = 1966:2019
nyears = length(years)
routes = 1:10
nroutes = 10
mean_count = log(5)
sd_route = 0.1
sd_ann_fluct = log(1.3) # average annual fluctuation = 15% or -13%
sd_ann_fluct_route = sd_ann_fluct*0.05 #5% variation among routes
t = c(0,0.02,0,-0.02,0.00)
cp = c(1986,1992,2001,2010)
yr_spans = as.integer(cut(years,breaks = c(-Inf,cp,Inf),ordered_result = TRUE))
sdt = rep(0.002,3)

betas_alphas = data.frame(route = routes,
                   alpha = rnorm(nroutes,mean_count,sd_route))
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
  
bb <- yr_spans[y]
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
                    nknots = 20,
                    sm_name = "year",
                    basis = "tp",
                    compare_rstanarm = FALSE,
                    diag.pen = TRUE)

# explroe basis function --------------------------------------------------
basis = as.data.frame(gam_dat$year_basis)
# basis_rstanarm = as.data.frame(gam_dat$rstanarm_basis)
# any(basis != basis_rstanarm)


stacked_basis = basis %>% 
  mutate(n = 1:54) %>% 
  pivot_longer(.,
               cols = starts_with("V"),
               names_to = "knot",
               values_to = "basis")

labls <- stacked_basis %>% filter(n == 54)
bas_gg = ggplot(data = stacked_basis,aes(x = n,y = basis,colour = knot))+
  geom_line()+
  theme(legend.position = "none")+
  scale_x_continuous(limits = c(NA,60))+
  scale_colour_viridis_d()+
  ggrepel::geom_label_repel(data = labls,aes(label = knot),nudge_x = 5)
print(bas_gg)


sdB_pr <- signif(sd(gam_dat$year_basis),2)#prior on the sd of t-distribution prior for sdBETA


#save(list = "count_df", file = "temp_data.RData")

#load("temp_data.RData")

sdB_pr <- 2.5 

stan_data = list(count = count_df$count,
                 year = count_df$yr,
                 route = count_df$route,
                 nyears = nyears,
                 ncounts = nrow(count_df),
                 nroutes = nroutes,
                 sdB_pr = sdB_pr,
                 nknots_year = gam_dat$nknots_year,
                 year_basispred = gam_dat$year_basis)




mod.file = "models/simple_GAM_route.stan"

## compile model
slope_model = cmdstan_model(stan_file=mod.file)


slope_stanfit <- slope_model$sample(data=stan_data,
                               refresh=100,
                               chains=4, iter_sampling =1000,
                               iter_warmup=1000,
                               parallel_chains = 4,
                               max_treedepth = 14,
                               adapt_delta = 0.8)


## run sampler on model, data
# slope_stanfit <- sampling(slope_model,
#                                data=stan_data,
#                                verbose=TRUE, refresh=100,
#                                chains=3,
#                           iter=1000,
#                                warmup=500,
#                           cores = 3,
#                                #pars = parms,
#                                control = list(adapt_delta = 0.8,
#                                               max_treedepth = 14))


shinyexp = FALSE
if(shinyexp){
  
slope_stanfit$save_output_files(dir = "output", basename = "temp", timestamp = FALSE, random = FALSE)
  csv_files <- dir("output/",pattern = "temp",full.names = TRUE)
  
  slope_rstan = rstan::read_stan_csv(csv_files)
  
  
  launch_shinystan(slope_rstan)
  
}
source("Functions/posterior_summary_functions.R")


years_df = unique(count_df[,c("yr","year")]) %>% 
  rename(y = yr)



sdBETA_samples = posterior_samples(fit = slope_stanfit,
                     parm = "sdBETA")
hist(sdBETA_samples$.value)



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
  geom_hline(yintercept = 0,alpha = 0.2)+
  geom_point(data = yeareffect,mapping = aes(x = y,y = tr_y),colour = "darkorange",alpha = 0.8,inherit.aes = FALSE)+
  coord_flip()

print(year_plot)






























































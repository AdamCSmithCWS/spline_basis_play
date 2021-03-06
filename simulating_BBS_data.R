### simulating BBS data using the core structure of the true data for a 
### selected species

library(bbsBayes)

library(tidyverse)
library(patchwork)
library(cmdstanr)
funct_dir <- "Functions/" # alternate = "C:/myFunctions/myFunctions/R/"
source(paste0(funct_dir,"get_basemap_function.R"))
source(paste0(funct_dir,"GAM_basis_function_mgcv.R"))
source(paste0(funct_dir,"neighbours_define.R"))
source(paste0(funct_dir,"posterior_summary_functions.R"))
source(paste0(funct_dir,"trend_function.R"))
source(paste0(funct_dir,"prepare-jags-data-alt.R")) # alternate version of the bbsBayes function that retains some extra information



# load BBS data -----------------------------------------------------------

strat_sel = "bbs_usgs"

fulldat <- stratify(by = strat_sel)

strat_df <- get_composite_regions(strat_sel)


# gather spatial info -----------------------------------------------------

strata_map = get_basemap(strat_sel,
                         append_area_weights = FALSE)



species = "Pacific Wren"
species_file = gsub(pattern = "([[:punct:]]|[[:blank:]])","",species)


model = "gamye"

jags_data = prepare_jags_data(strat_data = fulldat,
                              species_to_run = species,
                              model = model,
                              min_max_route_years = 2,
                              #min_year = first_year,
                              n_knots = 14)

years_df <- unique(data.frame(y = jags_data$year,
                              year = jags_data$r_year)) %>% 
  arrange(y)

str_link <- unique(data.frame(strat = jags_data$strat,
                              strat_name = jags_data$strat_name))
str_link <- str_link %>% arrange(strat)

realized_strata_map <- inner_join(strata_map,str_link,by = c("ST_12" = "strat_name")) %>% 
  arrange(strat)
# generate neighbourhoods -------------------------------------------------
buff_dist <- 10000




### re-arrange GEOBUGS formated nb_info into appropriate format for Stan model
car_stan_dat <- neighbours_define(real_strata_map = realized_strata_map,
                                  strat_link_fill = buff_dist,
                                  species = species,
                                  voronoi = FALSE,
                                  plot_dir = "maps/")




# fit the GAMYE spatial model in Stan -------------------------------------


stan_data = jags_data[c("ncounts",
                        "nstrata",
                        "count",
                        "strat",
                        "year",
                        "firstyr",
                        "nonzeroweight")]
stan_data[["observer"]] <- as.integer(factor(jags_data$ObsN))
stan_data[["nobservers"]] <- max(stan_data[["observer"]])

stan_data[["route"]] <- as.integer(factor(jags_data$route))
stan_data[["nroutes"]] <- max(stan_data[["route"]])
strat_route <- unique(data.frame(strat = stan_data$strat,
                                 route = stan_data$route))
strat_route <- strat_route %>% arrange(strat)

nroutes_strata <- as.integer(table(strat_route$strat))
maxnroutes_strata <- max(nroutes_strata)
stan_data[["nroutes_strata"]] <- nroutes_strata
stan_data[["maxnroutes_strata"]] <- maxnroutes_strata

rte_mat <- matrix(data = 0,
                  nrow = stan_data$nstrata,
                  ncol = maxnroutes_strata)
for(i in 1:stan_data$nstrata){
  rte_mat[i,1:nroutes_strata[i]] <- strat_route[which(strat_route$strat == i),"route"]
}

stan_data[["rte_mat"]] <- rte_mat

stan_data[["nyears"]] <- max(jags_data$year)
stan_data[["nknots_year"]] <- jags_data$nknots
stan_data[["year_basispred"]] <- jags_data$X.basis

stan_data[["N_edges"]] = car_stan_dat$N_edges
stan_data[["node1"]] = car_stan_dat$node1
stan_data[["node2"]] = car_stan_dat$node2




ncounts = stan_data$ncounts
nroutes = stan_data$nroutes
nstrata = stan_data$nstrata
nyears = stan_data$nyears
nknots_year = stan_data$nknots_year
nobservers = stan_data$nobservers



# # cmdStanR ----------------------------------------------------------------
mod.file = "models/gamye_iCAR.stan"

library(cmdstanr)
## compile model
slope_model <- cmdstan_model(mod.file)

init_def <- function(){ list(noise_raw = rnorm(ncounts,0,0.1),
                             strata_raw = rnorm(nstrata,0,0.1),
                             STRATA = 0,
                             eta = 0,
                             yeareffect_raw = matrix(rnorm(nstrata*nyears,0,0.1),nrow = nstrata,ncol = nyears),
                             obs_raw = rnorm(nobservers,0,0.1),
                             rte_raw = rnorm(nroutes,0,0.1),
                             sdnoise = 0.2,
                             sdobs = 0.1,
                             sdrte = 0.2,
                             sdbeta = runif(nknots_year,0.01,0.1),
                             sdBETA = 0.1,
                             sdyear = runif(nstrata,0.01,0.1),
                             nu = 10,
                             BETA_raw = rnorm(nknots_year,0,0.1),
                             beta_raw = matrix(rnorm(nknots_year*nstrata,0,0.01),nrow = nstrata,ncol = nknots_year))}



## run sampler on model, data
# data_file <- "tmpdata/tmp_data.json"
# write_stan_json(stan_data,file = data_file)
slope_stanfit <- slope_model$sample(
  data=stan_data,
  refresh=25,
  chains=3, iter_sampling=500,
  iter_warmup=800,
  parallel_chains = 3,
  #pars = parms,
  adapt_delta = 0.8,
  max_treedepth = 14,
  seed = 123,
  init = init_def,
  output_dir = "output",
  output_basename = paste0(species_file,"_cmdStan_out_",nyears))

csv_files <- dir("output/",pattern = paste0(species_file,"_cmdStan_out_",nyears),full.names = TRUE)



# save(list = c("slope_stanfit","stan_data","jags_data","model"),
#      file = paste0("output/cmdStan_",species_file,"ten_yr__gamye_iCAR.RData"))
# 
# 
# 
# slope_stanfit$save_object(file = paste0("output/cmdStan_file",
#                                         species_file,"_",nyears,"_gamye_iCAR.RDS"))
# 
# # export to csv and read in as rstan --------------------------------------
# 
# slope_stanfit$save_output_files(dir = "output",
#                                 basename = paste0(species_file,"_cmdStan_out_",nyears))
# 
# #sl_rstan <- As.mcmc.list(read_stan_csv(csv_files))
# if(shiny_explore){
  sl_rstan <- rstan::read_stan_csv(csv_files)
#   launch_shinystan(as.shinystan(sl_rstan))
#   
#   loo_stan = loo(sl_rstan)
# }
# # mod_sum = slope_stanfit$summary(variables = "n")



# Extract parameter values from fitted model ------------------------------

# sdyear[s]
  sdyear = posterior_samples(sl_rstan,
                                     parm = "sdyear",
                                     dims = "s") %>%
  posterior_sums(.,dims = "s") %>% 
  select(s,mean)

# route
  rte = posterior_samples(sl_rstan,
                                     parm = "rte",
                          dims = "r")%>% 
    posterior_sums(dims = "r") %>% 
    select(mean)
# sdobs
  obs = posterior_samples(sl_rstan,
                              parm = "obs",
                          dims = "o")%>% 
    posterior_sums(dims = "o") %>% 
    select(mean)
# sdstrata
  strata = posterior_samples(sl_rstan,
                            parm = "strata",
                            dims = "s")%>% 
    posterior_sums(dims = "s") %>% 
    select(mean)
# STRATA
  STRATA = posterior_samples(sl_rstan,
                            parm = "STRATA")%>% 
    posterior_sums() %>% 
    select(mean) %>% 
    as.numeric()

  sdnoise = posterior_samples(sl_rstan,
                             parm = "sdnoise")%>% 
    posterior_sums() %>% 
    select(mean) %>% 
    as.numeric()


# # Construct neighbourhood matrix ------------------------------------------
# simluate a series of slope
nb_mat = matrix(0,nrow = stan_data$nstrata,ncol = stan_data$nstrata)

  for(sr in 1:stan_data$nstrata){

    sc = car_stan_dat$node2[which(car_stan_dat$node1 == sr)]
    nb_mat[sr,sc] <- 1
    nb_mat[sc,sr] <- 1
  }
  
  
# 
#   var_mat = diag(x = 1,nrow = stan_data$nstrata,ncol = stan_data$nstrata)
#   var_mat = var_mat*0.1
#   cov_mat1 = nb_mat*0.5#trend_sd_btw_neighbours
# cov_mat = cov_mat1+var_mat
# stratum lat-long centres ------------------------------------------------

strat_centres = suppressWarnings(sf::st_centroid(realized_strata_map))
  strat_coords = data.frame(sf::st_coordinates(strat_centres))
  
 strat_coords$x_s = as.numeric(scale(strat_coords$X)) 
 strat_coords$y_s = as.numeric(scale(strat_coords$Y))  
 

# generate spatially correlated piecewise linear trends -------------------
  years = years_df$year 
  scale_trend_var_lat = 0.02 # 0.01 = ~range of 8%/year in trends values based on latitude
  scale_trend_var_long = 0.005 # 0.01 = ~range of 2%/year in trends values based on latitude
  
  true_trends = c(0.02,0.07,0,-0.07,0.03)
  true_breakpoints = c(1984,1993,2001,2010)
  trend_spans = as.integer(cut(years,breaks = c(-Inf,true_breakpoints,Inf),ordered_result = TRUE))
  ### trend spans is an integer indicator variable
  ### that identifies which of the true_trends is in effect in each year
 
    trend_mat = matrix(NA,nrow = length(true_trends),ncol = stan_data$nstrata)
  
  # loop through strata to fill in the trends based on lat and long
  for(ss in 1:stan_data$nstrata){
    for(tt in 1:length(true_trends)){
      lat_ef <- strat_coords[ss,"y_s"]*scale_trend_var_lat
      long_ef <- strat_coords[ss,"x_s"]*scale_trend_var_long
      
      trend_mat[tt,ss] <- true_trends[tt]+lat_ef+long_ef
      
    }
  }
  
  
 ### trend_mat is a matrix of true trends for each stratum and year-span
  
  

# balanced sample of routes, strata, years, w empir. obs turnover --------

rt_df = unique(data.frame(route = stan_data$route,
                          strat = stan_data$strat)) %>% 
    arrange(route)# complete list of routes

  event_df_log = as.data.frame(expand_grid(rt_df,years_df)) # complete list of all possible sampling events
  ## 
  if(nrow(event_df_log) != (stan_data$nyears*stan_data$nroutes)){
    stop("ERROR sampling event matrix structure is wrong")
  }
 
  midyr = floor(stan_data$nyears/2)
  ## TRUE route trajectories
  ## adding intercepts, routes, years
  

  
  for(rt in 1:stan_data$nroutes){
    
    wrty = which(event_df_log$route == rt & event_df_log$y == midyr)
    
    rf = rte$mean[rt]
    st = rt_df[rt,"strat"]
    sf = strata$mean[st]
    
    event_df_log[wrty,"true_traj"] <- rf + sf 
    
    
    for(y in (midyr+1):stan_data$nyears){ #mvoing forwards in time from midyr
    wrty = which(event_df_log$route == rt & event_df_log$y == y)
    wrtly = which(event_df_log$route == rt & event_df_log$y == (y-1))
    tr_tmp = trend_mat[trend_spans[y],st] #trend for that stratum and year
    event_df_log[wrty,"true_traj"] <-  event_df_log[wrtly,"true_traj"]+tr_tmp
        }#y
    
    for(y in (midyr-1):1){ #moving backwards in time (so y+1 and trend*-1)
      wrty = which(event_df_log$route == rt & event_df_log$y == y)
      wrtly = which(event_df_log$route == rt & event_df_log$y == (y+1))
      tr_tmp = (-1*trend_mat[trend_spans[y],st]) #negative of trend for that stratum and year
      event_df_log[wrty,"true_traj"] <-  event_df_log[wrtly,"true_traj"]+tr_tmp
    }#y
    
    
  }#r
  
  event_df = event_df_log
   event_df$route <- factor(event_df$route)
  plot_dim = ceiling(sqrt(max(rt_df$strat)))
  
  vis_plot = ggplot(data = event_df,aes(x = y,y = true_traj,colour = route))+
    geom_line()+
    theme(legend.position = "none")+
    facet_wrap(~strat,nrow = plot_dim,ncol = plot_dim)
  
  
  print(vis_plot)
  
  

# add yearly fluctuations -------------------------------------------------

  yeareff_df <- expand.grid(y = 1:stan_data$nyears,
                       strat = 1:stan_data$nstrata)
  
  for(st in 1:stan_data$nstrata){
    wst = which(yeareff_df$strat == st)
    yeareff_df[wst,"year_eff"] <- rnorm(stan_data$nyears,0,sdyear$mean[st])
  }
  
  for(rt in 1:stan_data$nroutes){
    for(y in 1:stan_data$nyears){
      
    wrty = which(event_df_log$route == rt & event_df_log$y == y)
    
    st = rt_df[rt,"strat"]
    wyeff = which(yeareff_df$strat == st & yeareff_df$y == y)
    
    
    event_df_log[wrty,"true_traj_yr"] <- event_df_log[wrty,"true_traj"]+yeareff_df[wyeff,"year_eff"] 

    }#y
    
    
  }#rt
  
  
  event_df = event_df_log
  event_df$route <- factor(event_df$route)
  plot_dim = ceiling(sqrt(max(rt_df$strat)))
  
  vis_plot = ggplot(data = event_df,aes(x = y,y = true_traj_yr,colour = route))+
    geom_line()+
    theme(legend.position = "none")+
    facet_wrap(~strat,nrow = plot_dim,ncol = plot_dim)
  
  
  print(vis_plot)
  

# add observer effects ----------------------------------------------------

 #add a route by year by observer
  obs_rt_df = unique(data.frame(route = stan_data$route,
                         y = stan_data$year,
                         obs = stan_data$observer)) ## not sure why the unique function is necessary
  
  
  event_wobs = event_df_log %>% dplyr::select(., route,strat,y,year) %>% 
    left_join(obs_rt_df,by = c("route","y"))
  
  ### filling in observers for sampling events that didn't happen
  ### this assumes that missing data were gathered by the last observer
  ### to survey the route
  for(rt in 1:stan_data$nroutes){
    
    wrw = which(event_wobs$route == rt)
    tmp = event_wobs[wrw,"obs"]
    frst = min(which(!is.na(tmp)))
    tmp[1:(frst-1)] <- tmp[frst]
    if(frst < length(tmp)){
    for(j in (frst+1):length(tmp)){
      tmp[j] <- ifelse(is.na(tmp[j]),tmp[j-1],tmp[j])
    }
    }
    event_wobs[wrw,"obs"] <- tmp
    event_wobs[wrw,"obs_eff"] <- obs$mean[tmp]
    
    rm("tmp")
  }
  
  
  noise = rnorm(nrow(event_df_log),0,sdnoise)
  
  #event_df_log
  event_df_log <- left_join(event_df_log,event_wobs,by = c("route","strat","y","year")) %>% 
    mutate(true_traj_yr_obs = true_traj_yr+obs_eff,
           noise_eff = rnorm(nrow(event_df_log),0,sdnoise/2),
           true_traj_yr_obs_noise = true_traj_yr_obs+noise_eff)
 
  
  event_df = event_df_log
  event_df$route <- factor(event_df$route)
  plot_dim = ceiling(sqrt(max(rt_df$strat)))
  
  vis_plot = ggplot(data = event_df,aes(x = y,y = exp(true_traj_yr_obs_noise),colour = route))+
    geom_line()+
    theme(legend.position = "none")+
    facet_wrap(~strat,nrow = plot_dim,ncol = plot_dim)
  
  
  print(vis_plot)
  

  stan_data_sim <- stan_data
  for(i in 1:stan_data_sim$ncounts){
    rt = stan_data_sim$route[i]
    yr = stan_data_sim$year[i]
    
    lmbd <- exp(event_df_log[which(event_df_log$route == rt &
                                   event_df_log$y == yr),"true_traj_yr_obs_noise"])
    smc <- rpois(1,lambda = lmbd)
    if(is.na(smc)){break}
    stan_data_sim$count[i] <- smc
    
  }
  
  







# fit simulated data ------------------------------------------------------

  
  ## run sampler on model, data
  # data_file <- "tmpdata/tmp_data.json"
  # write_stan_json(stan_data,file = data_file)
  slope_stanfit_sim <- slope_model$sample(
    data=stan_data_sim,
    refresh=25,
    chains=3, iter_sampling=500,
    iter_warmup=800,
    parallel_chains = 3,
    #pars = parms,
    adapt_delta = 0.8,
    max_treedepth = 14,
    seed = 123,
    init = init_def,
    output_dir = "output",
    output_basename = paste0(species_file,"_simulated_cmdStan_out_"))
  
  csv_files_sim <- dir("output/",pattern = paste0(species_file,"_simulated_cmdStan_out_"),full.names = TRUE)
  


  
  # compare simple estimates with true ---------------------------------------
  
  
  
  sl_rstan_sim <- rstan::read_stan_csv(csv_files_sim)
  #   launch_shinystan(as.shinystan(sl_rstan))
  #   
  #   loo_stan = loo(sl_rstan)
  # }
  # # mod_sum = slope_stanfit$summary(variables = "n")
  
  
  
  # Extract parameter values from fitted model ------------------------------
  
  
  beta_sim = posterior_samples(sl_rstan_sim,
                                 parm = "beta",
                                 dims = c("s","k")) %>% 
    posterior_sums(dims = c("s","k")) %>% 
    left_join(.,str_link,by = c("s" = "strat"))
  
  nfac = ceiling(sqrt(max(str_link$strat)))
  beta_comp = ggplot(data = beta_sim,aes(x = k,y = mean))+
    geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.5,width = 0)+
    geom_point()+
    facet_wrap(~strat_name,nrow = nfac,ncol = nfac)
  print(beta_comp)
  

# sdbeta and sdBETA ------------------------------------------------------------------

  
  sdbeta_sim = posterior_samples(sl_rstan_sim,
                                 parm = "sdbeta")%>% 
    posterior_sums() 
  
  sdbeta_sim 
  
  sdBETA_sim = posterior_samples(sl_rstan_sim,
                                 parm = "sdBETA")%>% 
    posterior_sums() 
  
  sdBETA_sim 
  
  # sdyear[s]
  sdyear_sim = posterior_samples(sl_rstan_sim,
                             parm = "sdyear",
                             dims = "s") %>% 
    posterior_sums(dims = "s")
  
  sdyear_comp = ggplot(data = sdyear_sim,aes(x = s,y = mean))+
    geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.5,width = 0)+
    geom_point()+
    geom_point(data = sdyear,aes(x = s,y = mean),colour = "blue",alpha = 0.5,shape = 2)
  print(sdyear_comp)
  
  
  # route
  rte_sim = posterior_samples(sl_rstan_sim,
                          parm = "rte",
                          dims = "r")%>% 
    posterior_sums(dims = "r") 
  rte = rename(rte,true_mean = mean)
  rte_sim = bind_cols(rte_sim,rte)
  
  rte_comp = ggplot(data = rte_sim,aes(x = true_mean,y = mean))+
    geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.5,width = 0)+
    geom_point()
  
  print(rte_comp)

  
  # sdobs
  obs_sim = posterior_samples(sl_rstan_sim,
                          parm = "obs",
                          dims = "o")%>% 
    posterior_sums(dims = "o")
  
  
  obs = rename(obs,true_mean = mean)
  obs_sim = bind_cols(obs_sim,obs)
  
  obs_comp = ggplot(data = obs_sim,aes(x = true_mean,y = mean))+
    geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.2,width = 0)+
    geom_point()
  
  print(obs_comp)

  
  
  # sdstrata
  strata_sim = posterior_samples(sl_rstan_sim,
                             parm = "strata",
                             dims = "s")%>% 
    posterior_sums(dims = "s") 
  strata = rename(strata,true_mean = mean)
  strata_sim = bind_cols(strata_sim,strata)
  
  strata_comp = ggplot(data = strata_sim,aes(x = true_mean,y = mean))+
    geom_errorbar(aes(ymin = lci,ymax = uci),alpha = 0.2,width = 0)+
    geom_point()
  
  print(strata_comp)

  
  # STRATA
  STRATA_sim = posterior_samples(sl_rstan_sim,
                             parm = "STRATA")%>% 
    posterior_sums() 
  
  STRATA_sim
  STRATA
  
  sdnoise_sim = posterior_samples(sl_rstan_sim,
                              parm = "sdnoise")%>% 
    posterior_sums() 
  sdnoise_sim
  sdnoise
  
  
  
  

# calculate trajectories and trends ---------------------------------------

  n_samples = posterior_samples(sl_rstan_sim,
                                parm = "n",
                                dims = c("s","y"))  %>% 
    left_join(.,years_df,by = "y") %>% 
    left_join(.,str_link,by = c("s" = "strat"))
  
  inds = n_samples %>% 
    posterior_sums(.,
                   dims = c("strat_name","year"))

  nsmooth_samples = posterior_samples(sl_rstan_sim,
                                     parm = "nsmooth",
                                     dims = c("s","y")) %>% 
    left_join(.,years_df,by = "y")  %>% 
    left_join(.,str_link,by = c("s" = "strat"))
  
  
  inds_smooth = nsmooth_samples %>% 
    posterior_sums(.,
                   dims = c("strat_name","year"))
  
  
  

# compare estimate trends with true ---------------------------------------

  
## calculate the estimated trends for each trend span and for the full time-series
  trend_intervals = c(min(years_df$year),true_breakpoints,max(years_df$year))
  
  trends = NULL
  for(yy in 1:length(true_trends)){
    y1 = trend_intervals[yy]
    y2 = trend_intervals[yy+1]
    
    trends_tmp = estimate_trend(nsmooth_samples,
                     scale = "strat_name",
                     start_year = y1,
                     end_year = y2)
    
    trends <- bind_rows(trends,trends_tmp)
    
    
    
  }
  
  trends = mutate(trends,period = as.character(start_year))
  # full time series trend
  trends_tmp <- estimate_trend(nsmooth_samples,
                       scale = "strat_name",
                       start_year = trend_intervals[1],
                       end_year = trend_intervals[length(trend_intervals)]) 
  trends = bind_rows(trends,trends_tmp)
  trends = mutate(trends,period = as.character(paste0(start_year,"-",end_year)))
  

# join the true trends in a data frame ------------------------------------

  str_tr_link <- NULL
  for(sp in 1:length(true_trends)){
   tmp <- data.frame(s = 1:max(str_link$strat),
                     period = as.character(paste0(trend_intervals[sp],"-",trend_intervals[sp+1])))
    for(st in 1:nrow(tmp)){
      tmp[st,"true_trend"] <- 100*(exp(trend_mat[sp,st])-1)
    }
   str_tr_link <- bind_rows(str_tr_link,tmp)
  }
  tmp <- data.frame(s = 1:max(str_link$strat),
                    period = as.character(paste0(trend_intervals[1],"-",trend_intervals[length(trend_intervals)])))
  for(st in 1:nrow(tmp)){
    tmp2 = mean(trend_mat[,st][trend_spans])
    tmp[st,"true_trend"] <- 100*(exp(tmp2)-1)
  }
  str_tr_link <- bind_rows(str_tr_link,tmp)
  str_tr_link <- left_join(str_tr_link,str_link,by = c("s" = "strat"))
  
  
  ## join true trends and estimated trends
  trends <- trends %>% left_join(.,str_tr_link,by = c("strat_name","period"))

#plot a comparison plot of true and estimates
nfac = ceiling(sqrt(max(str_link$strat)))

trends_comp = ggplot(data = trends,aes(x = s,y = mean_trend))+
  geom_errorbar(aes(ymin = lci_trend,ymax = uci_trend),alpha = 0.2,width = 0)+
  geom_point()+
  geom_point(aes(x = s,y = true_trend),colour = "blue",shape = 2)+
  facet_wrap(~period,nrow = 2,ncol = 3)

print(trends_comp)




# plot the hyperparameter trajectory ----------------------

Y_pred_samples = posterior_samples(fit = sl_rstan_sim,
                                   parm = "Y_pred",
                                   dims = "y") %>% 
  left_join(.,years_df,by = "y")

Y_pred = posterior_sums(Y_pred_samples,
                        dims = "year")

hyper_plot = ggplot(data = Y_pred,aes(x = year,y = exp(mean)))+
  geom_ribbon(aes(ymin = exp(lci),ymax = exp(uci)),alpha = 0.3)+
  geom_line()

print(hyper_plot)




# map the true and estimated trends ---------------------------------------



maps_out = vector(mode = "list",length = length(unique(trends$period)))
names(maps_out) <- unique(trends$period)

for(pp in names(maps_out)){
  
  tmp = vector(mode = "list",2)
  
  ttmp = trends %>% filter(period == pp)
  
  tmp[[1]] <- map_trends(trends = ttmp,
                         trend_col = "true_trend",
                         species = paste("True trends"))
  tmp[[2]] <- map_trends(trends = ttmp,
                         trend_col = "mean_trend",
                         species = paste("Estimated trends"))
  
  maps_out[[pp]] <- tmp
  
  
}





# plot the trajectories ---------------------------------------------------

# STratum level indices ---------------------------------------------------
first_year = min(years_df$year)

ind_sm <- nsmooth_samples %>% group_by(s,strat_name,year) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975)) %>% 
  mutate(type = "smooth")

ind_full <- n_samples %>% group_by(s,strat_name,year) %>% 
  summarise(mean = mean(.value),
            lci = quantile(.value,0.025),
            uci = quantile(.value,0.975))%>% 
  mutate(type = "full")


inds <- bind_rows(ind_sm,ind_full)

nstrata <- max(inds$s)
nf <- ceiling(sqrt(nstrata))

nonzero_w = data.frame(strat = 1:stan_data_sim$nstrata,
                       weights = stan_data_sim$nonzeroweight)

raw <- data.frame(count = stan_data_sim$count,
                  strat = stan_data_sim$strat,
                  year = stan_data_sim$year+(first_year-1)) %>% 
  left_join(.,str_link,by = "strat")




raw_full <- data.frame(count = exp(event_df_log$true_traj_yr_obs_noise),
                  strat = event_df_log$strat,
                  year = event_df_log$year) %>% 
  left_join(.,str_link,by = "strat")



raw_means <- raw %>% group_by(year,strat,strat_name) %>% 
  summarise(mean_count = mean(count),
            uqrt_count = quantile(count,0.75),
            n_surveys = n())%>% 
  left_join(.,nonzero_w,by = "strat") %>% 
  mutate(mean_count = mean_count*weights,
         uqrt_count = uqrt_count*weights)


raw_means_full <- raw_full %>% group_by(year,strat,strat_name) %>% 
  summarise(mean_count = mean(count),
            uqrt_count = quantile(count,0.75),
            n_surveys = n()) %>% 
  left_join(.,nonzero_w,by = "strat")%>% 
  mutate(mean_count = mean_count*weights,
         uqrt_count = uqrt_count*weights)



ind_fac <- ggplot(data = inds,aes(x = year,y = mean))+
  geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
  geom_line(aes(colour = type))+
  geom_point(data = raw_means,aes(x = year,y = mean_count),alpha = 0.2,inherit.aes = FALSE)+
  scale_y_continuous(limits = c(0,NA))+
  facet_wrap(~strat_name,nrow = nf,scales = "free")

#print(ind_fac)


# continental indices -----------------------------------------------------
# 
# a_weights <- as.data.frame(realized_strata_map) %>% 
#   rename(s = strat) %>% 
#   mutate(area = AREA/sum(AREA)) %>% 
#   select(s,area)
# 
# 
# I_sm <- nsmooth_samples %>% left_join(.,a_weights,by = "s") %>% 
#   mutate(.value = .value*area) %>% 
#   group_by(.draw,year) %>% 
#   summarise(.value = sum(.value)) %>% 
#   group_by(year) %>% 
#   summarise(mean = mean(.value),
#             lci = quantile(.value,0.025),
#             uci = quantile(.value,0.975)) %>% 
#   mutate(type = "smooth")
# 
# I_full <- n_samples %>% left_join(.,a_weights,by = "s") %>% 
#   mutate(.value = .value*area) %>% 
#   group_by(.draw,year) %>% 
#   summarise(.value = sum(.value)) %>% 
#   group_by(year) %>% 
#   summarise(mean = mean(.value),
#             lci = quantile(.value,0.025),
#             uci = quantile(.value,0.975)) %>%
#   mutate(type = "full")
# 
# 
# Is <- bind_rows(I_sm,I_full)
# 
# 
# I_plot <- ggplot(data = Is,aes(x = year,y = mean))+
#   geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.2)+
#   geom_line(aes(colour = type))+
#   labs(title = paste(species,"survey wide trajectory"))+
#   scale_y_continuous(limits = c(0,NA))
# 
# pdf(file = paste0("figures/",species_file,"trajectories.pdf"))
# 
# print(I_plot)






pdf(paste0("Figures/",species_file,"comparison_plots3normalregulrize.pdf"))
print(beta_comp)

print(strata_comp)
print(obs_comp)
print(rte_comp)
print(sdyear_comp)


print(hyper_plot)
print(trends_comp)


for(pp in names(maps_out)){
 tmp = maps_out[[pp]]
  print(tmp[[1]]+tmp[[2]]) 
  
}
for(st in str_link$strat_name){
  indst = filter(inds,strat_name == st)
  raw_mt = filter(raw_means,strat_name == st)
  raw_mt_full = filter(raw_means_full,strat_name == st)
  nsur = filter(raw,strat_name == st)
  ind_fac <- ggplot(data = indst,aes(x = year,y = mean))+
    geom_ribbon(aes(ymin = lci,ymax = uci,fill = type),alpha = 0.1)+
    geom_line(aes(colour = type))+
    geom_point(data = raw_mt,
               aes(x = year,y = mean_count),
               alpha = 0.2,colour = "blue",
               size = 1,inherit.aes = FALSE)+
    geom_point(data = raw_mt_full,
               aes(x = year,y = mean_count),
               alpha = 0.7,colour = "darkorange",
               size = 1,inherit.aes = FALSE)+
    geom_dotplot(data = nsur,aes(x = year),
                 inherit.aes = FALSE,binwidth = 1,
                 colour = grey(0.5),
                 fill = grey(0.5),
                 alpha = 0.1,
                 method = "histodot",dotsize = 0.5)+
    scale_y_continuous(limits = c(0,NA))+
    theme_classic()+
    labs(title = st)
  print(ind_fac)
}
dev.off()












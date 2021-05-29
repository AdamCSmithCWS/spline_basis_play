### simulating BBS data using the core structure of the true data for a 
### selected species

library(bbsBayes)

library(tidyverse)
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


species_file = gsub(pattern = "([[:punct:]]|[[:blank:]])","",species)

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
  
  true_trends = c(0,-0.07,0,-0.07,0)
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
  slope_stanfit <- slope_model$sample(
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
  
  csv_files <- dir("output/",pattern = paste0(species_file,"_simulated_cmdStan_out_"),full.names = TRUE)
  





























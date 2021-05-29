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
#   sl_rstan <- rstan::read_stan_csv(csv_files)
#   launch_shinystan(as.shinystan(sl_rstan))
#   
#   loo_stan = loo(sl_rstan)
# }
# # mod_sum = slope_stanfit$summary(variables = "n")


#select the stratum with the largest number of neighbours as the central
# stratum. This stratum will be assigned the base trajectory from which 
# all other trajectories will vary

core_stratum = as.integer(which.max(table(car_stan_dat$node1)))


































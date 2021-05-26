### a function to calculate the basis function(s) for a GAM component of a stan model

gam_basis <- function(orig.preds = dts[,"yr"],
                           nknots = 6,#number of internal knots
                           predpoints = NULL,
                           npredpoints = 100,
                           sm_name = "",
                      basis = "tp",
                      compare_rstanarm = FALSE,
                      diag.pen = TRUE){
  
  if(any(is.na(orig.preds) == T)){
    stop("This GAM formulation cannot handle missing values in the predictor")
  }
  
  # if(compare_brms){
  #   require(brms)
  # }
  # require(mgcv)
    
    dat = data.frame(x = orig.preds,
                     y = rnorm(length(orig.preds),0,0.1))
    if(is.null(predpoints)){
    predpoints = seq(min(orig.preds),max(orig.preds),length.out = npredpoints)
    }else{
      npredpoints <- length(predpoints)
    }
    
    dat_pred = data.frame(x = predpoints,
                          y = rnorm(length(predpoints),0,0.1))
    

    M = mgcv::smoothCon(mgcv::s(x,k = nknots+1, bs = basis),
                       data = dat,
                       absorb.cons=TRUE,#this drops the constant
                       diagonal.penalty=diag.pen) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 
    
    gamx.basis = M[[1]]$X
    
    M_pred = mgcv::smoothCon(mgcv::s(x,k = nknots+1, bs = basis),
                             data = dat_pred,
                             absorb.cons=TRUE,#this drops the constant
                             diagonal.penalty=diag.pen) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 
    
    gamx.basispred = M_pred[[1]]$X
    
    
    
    if(compare_rstanarm){
      
      require(rstanarm)
      M = stan_gamm4(y ~ s(x,k = nknots+1, bs = basis),data = dat,
                     chains=1, iter=600,
                     warmup=300,
                     cores = 1)
      # drop the constant term re-order columns to match brms (not necessary, but tidy),
      rstanarm_basis = matrix(M$x[,-1],nrow = length(orig.preds),ncol = nknots)
      
      
      
      
      if(any(rowSums(rstanarm_basis) != rowSums(gamx.basis))){print("WARNING brms basis differs from output")}
      
    }
    
    
    

    
    
   outlist <- list(gamx.basis = gamx.basis,
                  gamx.basispred = gamx.basispred,
                  orig.preds = orig.preds,
                  predpoints = predpoints,
                  nknots = nknots,
                  npredpoints = npredpoints,
                  M = M,
                  M_pred = M_pred)
  names(outlist) <- c(paste0(sm_name,"_basis"),
                      paste0(sm_name,"_basispred"),
                      "original_predictor_values",
                      paste0(sm_name,"_visualized_predictor_values"),
                      paste0("nknots_",sm_name),
                      paste0("npredpoints_",sm_name),
                      paste0(sm_name,"_smoothCon"),
                      paste0(sm_name,"_smoothCon_pred"))
  if(compare_rstanarm){
    outlist[["rstanarm_basis"]] <- rstanarm_basis
  }
  
  
  return(outlist)
  
  
  
}

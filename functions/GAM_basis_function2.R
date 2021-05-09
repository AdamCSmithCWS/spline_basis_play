### a function to calculate the basis function(s) for a GAM component of a JAGS or Stan model
### orig.preds = predictor variable to be modeled as a smooth, column from the dataset being modeled, must match with the indexing used in teh main likelihood definition (i.e., one value for each statistical unit)
### nknots = number of knots in the GAM smooth, (careful! you need to know what you're doing :-) 
### random = logical; TRUE if the GAM parameters are to be estimated as random effects (e.g., separate smooths within strata, centered on a mean smooth across all strata)
### npredpoints = number of points at which to teh user wants to visualize predictions for the smooth
### standardize = text; "z" if the predictor should be standardized (centered and re-scaled) to a mean = 0 and SD = 1, "range" if the user prefers a standard value wit hmean = mid-point of the range, and SD ~ 1 depending on the shape of the distribution (this approach can be helpful if predictor values are strongly skewed)
### even_gaps = T; logical, if TRUE the teh knots are spaced evenly along the range fot he predictor, if FALSE the data quantiles are used
### sm_name = text, text string to uniquely identify the predictor variable and its related parameters in the smooth


gam.basis.func2 <- function(orig.preds = dts[,"yr"],
                           nknots = 6,#number of internal knots
                           standardize = "z",
                           predpoints = NULL,
                           npredpoints = 100,
                           sm_name = ""){
  
  if(any(is.na(orig.preds) == T)){
    stop("This GAM formulation cannot handle missing values in the predictor")
  }
  
  
  require(mgcv)
    
    dat = data.frame(x = orig.preds,
                     y = orig.preds*0.1+rnorm(length(orig.preds),0,0.1))
    if(is.null(predpoints)){
    predpoints = seq(min(orig.preds),max(orig.preds),length.out = npredpoints)
    }else{
      npredpoints <- length(predpoints)
    }
    
    dat_pred = data.frame(x = predpoints,
                          y = predpoints*0.1+rnorm(length(predpoints),0,0.1))
    

    M = smoothCon(s(x,k = nknots, bs = "tp"),data = dat,
                       absorb.cons=TRUE,#this drops the constant
                       diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 
    
    
    gamx.basis = M[[1]]$X
    
    
    M_pred = smoothCon(s(x,k = nknots, bs = "tp"),data = dat_pred,
                       absorb.cons=TRUE,#this drops the constant
                       diagonal.penalty=TRUE) ## If TRUE then the smooth is reparameterized to turn the penalty into an identity matrix, with the final diagonal elements zeroed (corresponding to the penalty nullspace). 
    
    gamx.basispred = M_pred[[1]]$X
    
    
    
    
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
  
  
  return(outlist)
  
  
  
}
#############################################################################################################################
#START BY RUNNING CHECKS ON INPUT RASTER LAYERS
#############################################################################################################################

#' Checks input raster layers
#'
#' This function checks the dimensions of raster layers before running the power analysis simulations
#' @param occ  Run this function requires loading in the occupancy raster stack
#' @keywords 
#' @export
#' @examples

check.inputs <- function(occ){
  
  # Make sure number of sites specified is equal to number of sites loaded in
  if (load.sites == TRUE & all.loaded.sites == TRUE & n.sites != nrow(sites)) {
    cat('\n',"n.sites does not equal number of sites loaded in......check nrow(sites) or all.loaded.sites <- TRUE")
  }
  
  # Check dimensions of raster layers and print warning if they don't match
  if (ncol(occ) != ncol(det.method1) & ncol(det.method2) & ncol(det.method3) & ncol(parks) & ncol(remote)) {
    cat('\n',"Warning: Raster layers have unequal number of columns")
  }
  
  # Check that the number of layers in he determinstic disturbance raster stack is equal to the the length of the monitoring program (Tmax)
  if (model.disturbance == TRUE & disturbance.type == "deterministic"){
    if (nlayers(disturbance) != Tmax){
      cat('\n',"Warning: Length of Determinstic disturbance time series not equal to Tmax")}
  }
  
  #Check the number of species is consistent in the occupancy and detectability raster stacks
  if (nrow(species.list) != nlayers(occ) & nrow(species.occ) & nlayers(det.method1) & nlayers(det.method2) & nlayers(det.method3)) {
    cat('\n',"Warning: Number of species in species/list and/or species.occ not equal to number of raster layers in occ stack")
  }
  
  # Check monitoring occurs in last year of simulations 
  if (Tmax != max(s.years)) {
    cat('\n',"Warning: Monitoring must occur in final year (last element of s.years must equal Tmax")}
  
  # If power is to be estimated at a park level, make sure there's more than one park
  if (park.level == TRUE & length(n.park) == 1) {
    cat('\n',"Warning: Park level analysis selected but only one park identified in raster layer")
  }
  
  # Check there is at least one monitoring site
  if (n.sites == 0) {
    cat('\n',"Warning: No monitoring sites identified")
  }
  
  # Make sure you can't load sites and select sites with maximum occupancy
  if (load.sites == TRUE & all.loaded.sites == FALSE & new.site.selection == "maxocc") {
    cat('\n',"Warning: If new.site.selection == maxocc, set load.sites == FALSE")
  }
  
  # Make sure you can't choose two ways to select sites
  if (load.sites == TRUE & all.loaded.sites == FALSE & new.site.selection == "stratified") {
    cat('\n',"Warning: If new.site.selection == stratified, set load.sites == FALSE")
  }
  
  # Make sure you can't choose two ways to select sites
  if (load.sites == TRUE & all.loaded.sites == FALSE & new.site.selection == "manual") {
    cat('\n',"Warning: If new.site.selection == manual, set load.sites == FALSE")
  }
}


#' Site selection function
#'
#' This function positions sites in the landscape for monitoring. Sites can be: 1) pre-selected by loading a file of the XY coordinates; 2) randomly selected at the start of each simulation; 3) randomly stratified across environmental layers at the start of each simulation, or; 3) positioned on cells 
#' with the highest expected species richness. If sites are selected randomly, they can also be divided between remote and non-remote areas given the ratio R 
#' @param sites  Contains the XY coordinates of sites to be surveyed. The header of the x-coordinate should be an x, the header for the y-coordinate should be a y. If sites are to be selected randomly, dummy XY coordinates must still be provided. 
#' @param n.sites  The number of sites to be monitored. If sites are pre-determined, this must be equal to the number of rows in the site coordinate matrix      
#' @param R  The ratio of remote to non-remote sites. If R is 1, all sites will be in remote areas. If R is 0, no sites will be in remote areas, Setting R=0.6 means that 60 percent of randomly selected sites will be in remote areas
#' @param load.sites  Set to TRUE only monitor at fixed coordinates in the landscape, FALSE if sites are to be randomly selected, randomly stratified across environmental layers, or positioned on cells with the highest expected species richness
#' @param all.loaded.sites  Set to TRUE to monitor all of the XY coordinates loaded into the program, FALSE to monitor more than, or a subset of sites, loaded into the program. If n.sites is less than the number of loaded sites, a subset is selected randomly during each simulation. If n.sites is greater than the number of loaded sites, the remaining sites are positioned randomly throughout the landscape at the start of each simulation   
#' @param new.site.selection  Set to 'random' to select sites randomly throughout the landscape, 'stratified' to randomly select equal number of sites with environmental strata, or 'maxocc' to position sites on cells with the highest relative species richness
#' @param plot.sites  Set to TRUE to plot the location of sites for the first simulation, FALSE to not plot the location of sites
#' @keywords 
#' @export
#' @examples
#' #Load in fixed site coordinates to simulate monitoring at each site
#' load.sites <- TRUE
#' all.loaded.sites <- TRUE
#' new.site.selection <- 'random'
#' R <- 1
#' n.sites <- 150
#' plot.sites <- TRUE
#' xy.sites <- select.sites(sites, n.sites, R, all.loaded.sites, load.sites, new.site.selection, plot.sites) 
#' 
#' #Load in coordinates of 150 sites, but simulate monitoring at 200 sites
#' #Select the remaining 50 sites randomly throughout the landscape 
#' load.sites <- TRUE
#' all.loaded.sites <- FALSE
#' new.site.selection <- 'random'
#' R <- 1
#' n.sites <- 100
#' plot.sites <- TRUE
#' xy.sites <- select.sites(sites, n.sites, R, all.loaded.sites, load.sites, new.site.selection, plot.sites)
#' 
#' #Select 150 sites randomly stratified across environmental layers at the start of each simulation 
#' load.sites <- FALSE
#' all.loaded.sites <- FALSE
#' new.site.selection <- 'stratified'
#' R <- 1
#' n.sites <- 150
#' plot.sites <- TRUE
#' xy.sites <- select.sites(sites, n.sites, R, all.loaded.sites, load.sites, new.site.selection, plot.sites)
#' 
#' #Position 150 sites for monitoring on cells with the highest relative species richness
#' load.sites <- FALSE
#' all.loaded.sites <- FALSE
#' new.site.selection <- 'maxocc'
#' R <- 1
#' n.sites <- 150
#' plot.sites <- TRUE
#' xy.sites <- select.sites(sites, n.sites, R, all.loaded.sites, load.sites, new.site.selection, plot.sites)

select.sites <-  function(sites, n.sites, R, all.loaded.sites, load.sites, new.site.selection, plot.sites) { 
  
  #Identify and then separate remote and non-remote sites
  rem.ras <- acc.ras <- remote
  rem.ras[rem.ras==0] <- NA
  acc.ras[acc.ras==1] <- NA
  Rem <- remote[cellFromXY(remote, cbind(sites$x,sites$y))] #Extract remote ID from remote layer
  n.rem <- floor(n.sites*R) #Calculate number of remote sites
  n.acc <- n.sites - n.rem #Calculate number of non-remote sites
  XY.rem <- cbind(sites$x,sites$y)[which(Rem==1),] #Separate remote and non-remote sites
  XY.acc <- cbind(sites$x,sites$y)[which(Rem==0),]
  
  #Select all pre-selected sites
  if  (all.loaded.sites == TRUE) {
    XY.sites <- cbind(sites$x,sites$y)
  }
  
  #Select a random subset of pre-specified sites given the number of remote and non-remote sites is less then what is available
  if (load.sites == TRUE & all.loaded.sites == FALSE & n.rem <= nrow(XY.rem)  &  n.acc <= nrow(XY.acc)) {
    XY.rem <- XY.rem[sample(nrow(XY.rem), n.rem, replace=FALSE), ] #Randomly select remote and non-remote sites
    XY.acc <- XY.acc[sample(nrow(XY.acc), n.acc, replace=FALSE), ]
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  #The number of remote sites to survey is less than what is available. But the number of non-remote sites for survey is greater than what is available
  if (load.sites == TRUE & all.loaded.sites == FALSE & n.rem <= nrow(XY.rem)  &  n.acc > nrow(XY.acc)) {
    XY.rem <- XY.rem[sample(nrow(XY.rem), n.rem, replace=FALSE), ] #Randomly select remote and non-remote sites
    XY.acc <- rbind(XY.acc, sampleRandom(acc.ras,(n.acc-nrow(XY.acc)),na.rm=TRUE, xy=TRUE)[,1:2])
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  #The number of remote sites to survey is greater than what is available. But the number of 4WD sites for survey is less than what is available
  if (load.sites == TRUE & all.loaded.sites == FALSE & n.rem > nrow(XY.rem)  &  n.acc <= nrow(XY.acc)) {
    XY.acc <- XY.acc[sample(nrow(XY.acc), n.acc, replace=FALSE), ] #Randomly select remote and non-remote sites
    XY.rem <- rbind(XY.rem, sampleRandom(rem.ras,(n.rem-nrow(XY.rem)),na.rm=TRUE, xy=TRUE)[,1:2])
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  #The number of remote and 4WD sites to be surveyed is greater than the number of pre-existing sites
  if (load.sites == TRUE & all.loaded.sites == FALSE & n.rem > nrow(XY.rem)  &  n.acc > nrow(XY.acc)) {
    XY.acc <- rbind(XY.acc, sampleRandom(acc.ras, (n.acc-nrow(XY.acc)), na.rm=TRUE, xy=TRUE)[,1:2]) #Randomly select remote and non-remote sites
    XY.rem <- rbind(XY.rem, sampleRandom(rem.ras, (n.rem-nrow(XY.rem)), na.rm=TRUE, xy=TRUE)[,1:2])
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  if (load.sites == FALSE & new.site.selection == "random") {
    XY.acc <- XY.rem <- NULL
    if (n.acc > 0) {XY.acc <- sampleRandom(acc.ras, n.acc, na.rm=TRUE, xy=TRUE)[,1:2]} #Randomly select remote and non-remote sites
    if (n.rem > 0){XY.rem <- sampleRandom(rem.ras, n.rem, na.rm=TRUE, xy=TRUE)[,1:2]}
    XY.sites <- rbind(XY.rem, XY.acc) #Combine remote and non-remote site
  }
  
  if (load.sites == FALSE & new.site.selection == "stratified") {
    classes <- unique(values(stratify))
    classes <- na.omit(classes)
    strata <- length(classes)
    XY.sites <- sampleStratified(stratify, size=floor(n.sites/strata), na.rm=TRUE, xy=TRUE)[,2:3] 
    if (nrow(XY.sites) < n.sites) {
      extra <- sampleRandom(stratify, (n.sites-nrow(XY.sites)), na.rm=TRUE, xy=TRUE)[,1:2]
      XY.sites <- rbind(XY.sites, extra)
    }
  }
  
  if (load.sites == FALSE & new.site.selection == "manual") {
    X <- sum(occ)
    plot(X)
    cat('\n',"Click on the map to select sites then press Esc") 
    XY.sites <- "click"(X, n=Inf, id=FALSE, xy=TRUE, type="p", col="red")
    XY.sites <- XY.sites[,1:2]
    rm(X)
  }
  
  if (load.sites == FALSE & new.site.selection == "maxocc") {
    X <- sum(occ)
    test2 <- unlist(sort(as.vector(X), decreasing = TRUE, index.return=TRUE)$x)[n.sites]
    XY.sites <- xyFromCell(X, Which(X >= test2, cells=TRUE,na.rm=TRUE))
    rm(X)
  }
  
  if (plot.sites == TRUE) {
    plot(remote)
    points(XY.sites, col="red", pch=19)
  }
  
  rm(acc.ras, rem.ras)
  return(XY.sites)
}


#' Significance level function
#'
#' This function returns a critical value needed to calculate confidence intervals around the trend parameter given the Type I error rate and whether it is a one or two-tailed significance test.
#' @param two.tailed Set to TRUE for a two-tailed test, FALSE for a one-tailed test. 
#' @param alpha Set the type one error rate. Choose between 0.1, 0.05 and 0.01
#' @keywords 
#' @export
#' @examples
#' #Return the critical value for a two-tailed test with a significance level of 0.05
#' two.tailed <- TRUE
#' alpha <- 0.05
#' sig.test(two.tailed, alpha)
#' 
#' #Return the critical value for a one-tailed test with a significance level of 0.1
#' two.tailed <- FALSE
#' alpha <- 0.1
#' sig.test(two.tailed, alpha)

sig.test <- function(two.tailed, alpha) {
  values <- matrix(c(1.28, 1.645, 1.65, 1.96, 2.33, 2.58), ncol=3, nrow=2)
  if (two.tailed == TRUE) {ind <- 2} else {ind <- 1}
  if (alpha == 0.1) {value <- values[ind,1]}
  if (alpha == 0.05) {value <- values[ind,2]}
  if (alpha == 0.01) {value <- values[ind,3]}
  return(value)
}

#' Stochastic disturbance function
#'
#' This function models the incidence of a stochastic disturbance at monitoring sites given the probability of a cell being disturbed as a function of time since the last disturbance.
#' It returns a vector specifying whether each site is disturbed (1) or not (0) in a given year jj
#' @param time.hist An array specifying the time since a disturbance at each monitoring site prior to simulations
#' @param Pr.dist A vector specifying the probability of a cell being disturbed given the time since the last disturbance.
#' @param n.sites The number of sites monitored during each simulation
#' @keywords 
#' @export
#' @examples

stochastic.disturbance <- function(time.hist, Pr.dist, n.sites){
  burn <- Pr.dist[time.hist]
  burn <- ifelse(runif(n.sites)<burn,1,0)
  return(burn)
}

#' Refit occupancy raster layers function
#'
#' This function is only called if stochastic disturbances are modelled at sites. It re-maps occupancy for disturbance-sensitive species given the simulated disturbance events and the relationship between distrubances and occupancy
#' It returns a new raster stack with the updated occupancy raster layers for each species
#' @param occ.new A raster stack of occupancy maps. Raster layers for fire sensitive species are updated depending on the simulated fire history at time jj
#' @param layers An array of covararite values at each site. These are used to update the occ.new raster stack   
#' @param time.fire A vector specifying the time since fire at each site in year jj given the simulated fire history
#' @param fire.freq A vector specifying the number of fires at sites during a 15 year moving window
#' @param jj The year of the monitoring program
#' @param dist.occ.time An integer indentifying which layer in the covariate raster stack is the time since the last disturbance layer
#' @param dist.occ.freq An integer indentifying which layer in the covariate raster stack is the disturbance frequency layer
#' @param n.species Number of species included in power simulations
#' @keywords 
#' @export
#' @examples

refit.occ <- function(occ.new, layers, time.fire, fire.freq, jj, dist.occ.time, dist.occ.freq, n.species) {
  for (i in 1:n.species) {
    covr <- as.numeric(species.occ[i,-1])
    covr.int <- covr[1]
    covr.pred <- covr[-1]
    if (covr.pred[dist.occ.time] | covr.pred[dist.occ.freq] !=0) { #Re-map occupancy for fire-sensitive species only
      time.scale <- array(scale(time.fire[,jj])[,1])
      freq.scale <- array(scale(fire.freq[,jj])[,1])
      
      layers[,dist.occ.time] <- time.scale #Replace disturbance layers
      layers[,dist.occ.freq] <- freq.scale
      
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit.occ <- covr.int + rowSums(update)
      
      occ.new[,i,jj] <- exp(logit.occ)/(1+exp(logit.occ))
    }
  }
  return(occ.new)
}

#' Refit detectability raster layers function
#'
#' This function is called only if a stochastic disturbance is modelled at sites. It re-fits the detectability raster layers for disturbance sensitive species given simulated disturbance events and the disturbance history
#' A new raster stack is returned for each detection method
#' @param det.method1.new An array of detectability estimates for each species at sites over time using method 1. 
#' @param det.method2.new An array of detectability estimates for each species at sites over time using method 2.  
#' @param det.method3.new An array of detectability estimates for each species at sites over time using method 3. 
#' @param det.method4.new An array of detectability estimates for each species at sites over time using method 4. 
#' @param layers An array of covararites at each site used to update the detectability arrays given the simulated disturbances 
#' @param time.fire A vector specifying the time since a disturbance at each site in year jj given the simulated disturbances
#' @param fire.freq A vector specifying the number of disturbances at sites during a 15 year moving window
#' @param jj The year of the monitoring program
#' @param dist.occ.time An integer indentifying which layer in the covariate raster stack is the time since the last disturbance layer
#' @param dist.occ.freq An integer indentifying which layer in the covariate raster stack is the disturbance frequency layer
#' @param n.species Number of species included in power simulations
#' @keywords 
#' @export
#' @examples

refit.det <- function(det.method1.new, det.method2.new, det.method3.new, det.method4.new, layers, time.fire, fire.freq, jj, dist.occ.time, dist.occ.freq, n.species) {
  
  time.scale <- array(scale(time.fire)[,jj])
  freq.scale <- array(scale(fire.freq)[,jj])
  layers[,dist.occ.time] <- time.scale 
  layers[,dist.occ.freq] <- freq.scale
  obj <- list(det.method1.new, det.method2.new, det.method3.new, det.method4.new)
  
  for (i in 1:n.species) {
    covr <- as.numeric(species.det[i,-1])
    covr.int <- covr[1]
    covr.method <- covr[2:4]
    covr.pred <- covr[-c(1:4)]
    
    #Species detected using 1 method
    if ((sum(species.list[i,-1]) == 1) & (covr.pred[dist.occ.time] | covr.pred[dist.occ.freq] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      #Method 1
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + rowSums(update)
      obj[[pos]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
    #Species detected using 2 methods
    if ((sum(species.list[i,-1]) == 2) & (covr.pred[dist.occ.time] | covr.pred[dist.occ.freq] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      #Method 1
      update <- sweep(layers,MARGIN=2,covr.pred,'*') 
      logit <- covr.int + rowSums(update)
      obj[[pos[1]]][,i,jj] <- exp(logit)/(1+exp(logit))
      #Method 2
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + covr.method[1] + rowSums(update)
      obj[[pos[2]]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
    #Species detected using 3 methods
    if ((sum(species.list[i,-1]) == 3) & (covr.pred[dist.occ.time] | covr.pred[dist.occ.freq] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      #Method 1
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + rowSums(update)
      obj[[pos[1]]][,i,jj] <- exp(logit)/(1+exp(logit))
      #Method 2
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + covr.method[1] + rowSums(update)
      obj[[pos[2]]][,i,jj] <- exp(logit)/(1+exp(logit))
      #Method 3
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + covr.method[2] + rowSums(update)
      obj[[pos[3]]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
    #Species detected using 4 methods
    if ((sum(species.list[i,-1]) == 4) & (covr.pred[dist.occ.time] | covr.pred[dist.occ.freq] != 0)){ 
      pos <- which(species.list[i,-1]==1) 
      #Method 1
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + rowSums(update)
      obj[[pos[1]]][,i,jj] <- exp(logit)/(1+exp(logit))
      #Method 2
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + covr.method[2] + rowSums(update)
      obj[[pos[2]]][,i,jj] <- exp(logit)/(1+exp(logit))
      #Method 3
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + covr.method[3] + rowSums(update)
      obj[[pos[3]]][,i,jj] <- exp(logit)/(1+exp(logit))
      #Method 4
      update <- sweep(layers,MARGIN=2,covr.pred,'*')
      logit <- covr.int + covr.method[4] + rowSums(update)
      obj[[pos[4]]][,i,jj] <- exp(logit)/(1+exp(logit))
    }
    
  }
  return(obj)
}


#' Fit occupancy model with one detection method
#'
#' This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected using only one detection method
#' A trend in occupancy is estimated from the simulated detection histories. Confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
#' A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
#' intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from within each park-level management unit. 
#' @param method The detection method relevant to species ss
#' @param repeats The number of repeat visits for the specified detection method
#' @param s.years A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)
#' @param n.sites The number of sites in which to simulate monitoring 
#' @param xy.sites The XY coordinates of monitored sites, plus the park site ID
#' @param park.ID A vector identifying which park each site is in
#' @param park.level Set to TRUE if power is estimated within park-level management unit, FALSE to estimate power only across the landscape
#' @param powcnt A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape
#' @param fail A vector that keeps track of how many times the occupancy model could not be fitted to the simulated detection histories from all sites in the landscape
#' @param pow.park A vector that keeps track of how many times a significant trend in occupancy is detected within each park
#' @param fail.park A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories generated from within each park
#' @param value The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test
#' @param ss An index to loop through each species 
#' @param two.tailed Set to TRUE if conducting a two-tailed test, FALSE to conduct a one-tailed test
#' @param n.park The number of parks in which to estimate power
#' @keywords 
#' @export
#' @examples

fit.occ.1method <- function(method, repeats, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats == 1) {
    y <- aperm(method[,ss,,s.years],c(1,2))
    dim(y)<- c(n.sites*length(s.years), repeats)
  } else {
    y <- aperm(method[,ss,,s.years],c(1,3,2))
    dim(y)<- c(n.sites*length(s.years), repeats)
  }
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = NULL)
  
  inits<- rbind(c(0,0,0),c(-1,-1,-1),c(1,1,1))
  for (vv in 1:nrow(inits)) {
    mod<- try(occu(~1 ~Time, umf, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
    if (!is(mod,"try-error")){break}
  }
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        
        inits<- rbind(c(0,0,0),c(-1,-1,-1),c(1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~1 ~Time, y.park, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}

#' Fit occupancy model with two detection methods
#'
#' This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected with two detection methods
#' A trend in occupancy is estimated and confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
#' A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
#' intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from each regional level management unit. 
#' @param method1 The first detection method relevant to species ss
#' @param method2 The second detection method relevant to species ss
#' @param repeats1 The number of repeat visits for the first detection method
#' @param repeats2 The number of repeat visits for the second detection method
#' @param s.years A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)
#' @param n.sites The number of sites monitored
#' @param xy.sites The XY coordinates of monitored sites
#' @param park.ID A vector identifying the park that each site is positioned within
#' @param park.level Set to TRUE if power is estimated within park-level management units, FALSE if power is only estimated across the landscape
#' @param powcnt A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape
#' @param fail A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a landscape level in unmarked
#' @param pow.park A vector that keeps track of how many times a significant trend in occupancy is detected within each park
#' @param fail.park A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a park level in unmarked
#' @param value The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test
#' @param ss An index to loop through each species 
#' @param two.tailed Set to TRUE if conducting a two-tailed test, FALSE otherwise
#' @param n.park The number of parks in which to estimate power
#' @keywords 
#' @export
#' @examples

#Function to fit occupancy model when only 2 methods are used to detect species s
fit.occ.2method <- function(method1, method2, repeats1, repeats2, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats1 == 1) {
    y1 <- aperm(method1[,ss,,s.years],c(1,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  } else {
    y1 <- aperm(method1[,ss,,s.years],c(1,3,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  }
  if (repeats2 == 1) {
    y2 <- aperm(method2[,ss,,s.years],c(1,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  } else {
    y2 <- aperm(method2[,ss,,s.years],c(1,3,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  }
  
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  
  m1.code <- matrix(factor(1),nrow=nrow(xy.sites)*length(s.years), ncol=repeats1)
  m2.code <- matrix(factor(2),nrow=nrow(xy.sites)*length(s.years), ncol=repeats2)
  obsCovs <- data.frame(as.matrix(cbind(m1.code,m2.code)))
  colnames(obsCovs) <- c(rep("M1",repeats1),rep("M2",repeats2))
  
  umf <- unmarkedFrameOccu(y = cbind(y1, y2), siteCovs = siteCovs, obsCovs = list(METHOD = obsCovs))
  
  inits<- rbind(c(0,0,0,0),c(-1,-1,-1,-1),c(1,1,1,1))
  for (vv in 1:nrow(inits)) {
    mod<- try(occu(~METHOD ~Time, umf, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
    if (!is(mod,"try-error")){break}
  }
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        
        inits<- rbind(c(0,0,0,0),c(-1,-1,-1,-1),c(1,1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~METHOD ~Time, y.park, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}

#' Fit occupancy model with three detection methods
#'
#' This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected using 3 detection methods
#' A trend in occupancy is estimated and confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
#' A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
#' intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from each regional level management unit. 
#' @param method1 The first detection method relevant to species ss
#' @param method2 The second detection method relevant to species ss
#' @param method3 The third detection method relevant to species ss
#' @param repeats1 The number of repeat visits for the first detection method
#' @param repeats2 The number of repeat visits for the second detection method
#' @param repeats3 The number of repeat visits for the third detection method
#' @param s.years A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)
#' @param n.sites The number of sites monitored
#' @param xy.sites The XY coordinates of monitored sites
#' @param park.ID A vector specifying that location of each site with sub-level parks
#' @param park.level Set to TRUE is power is estimated within regional level management unit, FALSE otherwise
#' @param powcnt A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape
#' @param fail A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a landscape level in unmarked
#' @param pow.park A vector that keeps track of how many times a significant trend in occupancy is detected within each park
#' @param fail.park A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a park level in unmarked
#' @param value The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test
#' @param ss An index to loop through each species 
#' @param two.tailed Set to TRUE if conducting a two-tailed test, FALSE otherwise
#' @param n.park The number of parks in which to estimate power
#' @keywords 
#' @export
#' @examples

#Function to fit occupancy model when 3 methods are used to detect species s
fit.occ.3method <- function(method1, method2, method3, repeats1, repeats2, repeats3, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats1 == 1) {
    y1 <- aperm(method1[,ss,,s.years],c(1,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  } else {
    y1 <- aperm(method1[,ss,,s.years],c(1,3,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  }
  if (repeats2 == 1) {
    y2 <- aperm(method2[,ss,,s.years],c(1,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  } else {
    y2 <- aperm(method2[,ss,,s.years],c(1,3,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  }
  if (repeats3 == 1) {
    y3 <- aperm(method3[,ss,,s.years],c(1,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  } else {
    y3 <- aperm(method3[,ss,,s.years],c(1,3,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  }
  
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  
  m1.code <- matrix(factor(1),nrow=nrow(xy.sites)*length(s.years), ncol=repeats1)
  m2.code <- matrix(factor(2),nrow=nrow(xy.sites)*length(s.years), ncol=repeats2)
  m3.code <- matrix(factor(3),nrow=nrow(xy.sites)*length(s.years), ncol=repeats3)
  obsCovs <- data.frame(as.matrix(cbind(m1.code, m2.code, m3.code)))
  colnames(obsCovs) <- c(rep("M1",repeats1), rep("M2",repeats2), rep("M3",repeats3))
  
  umf <- unmarkedFrameOccu(y = cbind(y1, y2, y3), siteCovs = siteCovs, obsCovs = list(METHOD = obsCovs))
  
  inits<- rbind(c(0,0,0,0,0),c(-1,-1,-1,-1,-1),c(1,1,1,1,1))
  for (vv in 1:nrow(inits)) {
    mod<- try(occu(~METHOD ~Time, umf, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
    if (!is(mod,"try-error")){break}
  }
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        inits<- rbind(c(0,0,0,0,0),c(-1,-1,-1,-1,-1),c(1,1,1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~METHOD ~Time, y.park, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}


#' Fit occupancy model with four detection methods
#'
#' This function fits an occupancy model to simulated detection histories using the package unmarked for species that are detected using 4 detection methods
#' A trend in occupancy is estimated and confidence intervals are calculated depending on the Type I error rate. A one-tailed or two-tailed significance test is then conducted on the trend parameter
#' A one-tailed test looks to see if the upper or lower confidence interval is greater than or less than zero. A two-tailed test assesses whether both the upper and lower confidence
#' intervals have the same sign (i.e. are both positive or negative). If park.power = TRUE, model fitting is repeated on sites from each regional level management unit. 
#' @param method1 The first detection method relevant to species ss
#' @param method2 The second detection method relevant to species ss
#' @param method3 The third detection method relevant to species ss
#' @param method4 The fourth detection method relevant to species ss
#' @param repeats1 The number of repeat visits for the first detection method
#' @param repeats2 The number of repeat visits for the second detection method
#' @param repeats3 The number of repeat visits for the third detection method
#' @param repeats4 The number of repeat visits for the fourth detection method
#' @param s.years A vector specifying the years that monitoring occurs. Note, monitoring must be done in the final year (i.e. Tmax)
#' @param n.sites The number of sites monitored
#' @param xy.sites The XY coordinates of monitored sites
#' @param park.ID A vector specifying that location of each site with sub-level parks
#' @param park.level Set to TRUE is power is estimated within regional level management unit, FALSE otherwise
#' @param powcnt A vector that keeps track of how many times a significant trend in occupancy is detected across the landscape
#' @param fail A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a landscape level in unmarked
#' @param pow.park A vector that keeps track of how many times a significant trend in occupancy is detected within each park
#' @param fail.park A vector that keeps track of how many times the occupancy model could not be fitted to simulated detection histories at a park level in unmarked
#' @param value The critical value used to calculate confidence intervals around the trend parameter, depending on the Type I error rate and a one-tailed or two-tailed test
#' @param ss An index to loop through each species 
#' @param two.tailed Set to TRUE if conducting a two-tailed test, FALSE otherwise
#' @param n.park The number of parks in which to estimate power
#' @keywords 
#' @export
#' @examples

fit.occ.4method <- function(method1, method2, method3, method4, repeats1, repeats2, repeats3, repeats4, s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park) {
  
  if (repeats1 == 1) {
    y1 <- aperm(method1[,ss,,s.years],c(1,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  } else {
    y1 <- aperm(method1[,ss,,s.years],c(1,3,2))
    dim(y1)<- c(n.sites*length(s.years), repeats1)
  }
  if (repeats2 == 1) {
    y2 <- aperm(method2[,ss,,s.years],c(1,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  } else {
    y2 <- aperm(method2[,ss,,s.years],c(1,3,2))
    dim(y2)<- c(n.sites*length(s.years), repeats2)
  }
  if (repeats3 == 1) {
    y3 <- aperm(method3[,ss,,s.years],c(1,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  } else {
    y3 <- aperm(method3[,ss,,s.years],c(1,3,2))
    dim(y3)<- c(n.sites*length(s.years), repeats3)
  }
  if (repeats4 == 1) {
    y4 <- aperm(method4[,ss,,s.years],c(1,2))
    dim(y4)<- c(n.sites*length(s.years), repeats4)
  } else {
    y4 <- aperm(method4[,ss,,s.years],c(1,3,2))
    dim(y4)<- c(n.sites*length(s.years), repeats4)
  }
  
  siteCovs <- data.frame(matrix(rep(s.years,each=nrow(xy.sites)),nrow(xy.sites)*length(s.years),1))
  park.ID.year <- rep(park.ID, length(s.years))
  siteCovs <- cbind(siteCovs, park.ID.year)
  colnames(siteCovs) <- c("Time", "Park")
  
  m1.code <- matrix(factor(1),nrow=nrow(xy.sites)*length(s.years), ncol=repeats1)
  m2.code <- matrix(factor(2),nrow=nrow(xy.sites)*length(s.years), ncol=repeats2)
  m3.code <- matrix(factor(3),nrow=nrow(xy.sites)*length(s.years), ncol=repeats3)
  m4.code <- matrix(factor(4),nrow=nrow(xy.sites)*length(s.years), ncol=repeats4)
  obsCovs <- data.frame(as.matrix(cbind(m1.code, m2.code, m3.code, m4.code)))
  
  colnames(obsCovs) <- c(rep("M1",repeats1), rep("M2",repeats2), rep("M3",repeats3), rep("M4",repeats4))
  
  umf <- unmarkedFrameOccu(y = cbind(y1, y2, y3, y4), siteCovs = siteCovs, obsCovs = list(METHOD = obsCovs))
  
  inits<- rbind(c(0,0,0,0,0,0),c(-1,-1,-1,-1,-1,-1),c(1,1,1,1,1,1))
  for (vv in 1:nrow(inits)) {
    mod<- try(occu(~METHOD ~Time, umf, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
    if (!is(mod,"try-error")){break}
  }
  
  if (is(mod,"try-error")) {fail[ss] <- fail[ss] + 1} else {
    if (!is.nan(SE(mod)[2])) {
      Time.mean <- coef(mod)[2]
      Time.CI <- SE(mod)[2] * value
      if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
      if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
        powcnt[ss] <- powcnt[ss] + 1}
    } else {fail[ss] <- fail[ss] + 1}}
  
  if (park.level == TRUE) {
    for (pp in n.park) {
      park.site <- which(xy.sites[,3] == pp)
      if (length(park.site) == 0) {pow.park[pp,ss] <- pow.park[pp,ss]}
      if (length(park.site) > 0) {
        
        y.park <- umf[which(park.ID.year == pp),]
        inits<- rbind(c(0,0,0,0,0,0),c(-1,-1,-1,-1,-1,-1),c(1,1,1,1,1,1))
        for (vv in 1:nrow(inits)) {
          mod<- try(occu(~METHOD ~Time, y.park, control=list(maxit=10000), starts=inits[vv,]), silent=TRUE)
          if (!is(mod,"try-error")){break}
        }
        
        if (is(mod,"try-error")) {fail.park[pp,ss] <- fail.park[pp,ss] + 1} else {
          if (!is.nan(SE(mod)[2])) {
            Time.mean <- coef(mod)[2]
            Time.CI <- SE(mod)[2] * value
            
            if (two.tailed == TRUE & sign(Time.mean + Time.CI) == sign(Time.mean - Time.CI)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
            if (two.tailed == FALSE & (Time.mean + Time.CI < 0 | Time.mean - Time.CI > 0)) { #Add on-tailed test option here
              pow.park[pp,ss] <- pow.park[pp,ss] + 1}
          } else {fail.park[pp,ss] <- fail.park[pp,ss] + 1}}
      }
    }
  }			
  return(list(powcnt, fail, pow.park, fail.park))
}

#' Run power analysis
#'
#' This function runs the main power analysis. It starts by creating empty vectors to record the number of times a significant trend is detected in the simulated detection histories
#' It then simulates fire at sites if model.fire = TRUE, and models a decline in occupancy over time depending on the specified effect size. 
#' Detection histories are simulated at each site given the occupancy status at each site, the number of survey methods used to detect a species and the detection probability of each relevent method.
#' Simulated detection histories are loaded into the package unmarked, and a trend in occupancy is modelled over time. The proportion of times a significant trend is detected is recorded.
#' This process is repeated for each effect size.
#' @param effect.size An integer specifying the proportional reduction in occupancy between the start and end of the monitoring program
#' @param nsims An integer specifying the number of simulations 
#' @param alpha The Type I error rate. Must be either 0.1, 0.05 or 0.01
#' @param Tmax An integer specifying the length of the monitoring program
#' @param s.years A vector specifying the years in which monitoring occurs. Note, the final year or monitoring must be equal to Tmax
#' @param trend Set to 'increasing' to model an increasing trend in occupancy or 'decreasing' to model a decreasing trend
#' @param sites An array containing the XY coordinates of pre-specified monitoring sites
#' @param model.disturbance Set to TRUE to model a disturbance, FALSE otherwise
#' @param disturbance.type Set to 'deterministic' or 'stochastic' depending on the type of disturbance
#' @param species.list An array containing the name of each species in the first column, and a zero or one describing how many methods are used to detect the species
#' @param dist.occ.time An integer specifying which layer in the raster covariate stack represents time since the last disturbance
#' @param dist.occ.freq An integer specifying which layer in the raster covariate stack represents the number of disturbances 
#' @param park.level Set to TRUE if power is to be estimated within sub-level management units
#' @param xy.sites The XY coordinates of monitoring sites
#' @param R Ratio of remote to non-remote sites
#' @param two.tailed Set to TRUE to conduct a two-tailed significance test, FALSE if conducting a one-tailed test
#' @param plot.sites Set to TRUE to plot the location of sites at the start of each simulation, FALSE otherwise
#' @param n.sites The number of sites to be monitored.  
#' @param n.species The number of species for which monitoring is simulated. Is equal to the number of layers in the occ raster stack
#' @param n.park The number of nested management units in which to estimate power
#' @param all.loaded.sites Set to TRUE if all of the loaded sites are to be monitored, FALSE otherwise
#' @param loaded.sites Set to TRUE if monitoring is simulated at XY-coordinates provided, FALSE otherwise
#' @param new.site.selection Set to 'random' to select new sites randomly, 'stratified' to randomly position sites within environmental strata, or 'maxocc' to position sites on cells with the highest relative species richness
#' @param occ.time An array containing the occupancy value for each species at monitoring sites over time 
#' @param det.method1.time An array containing the detectability values of method 1 for each species at monitoring sites over time
#' @param det.method2.time An array containing the detectability values of method 2 for each species at monitoring sites over time
#' @param det.method3.time An array containing the detectability values of method 3 for each species at monitoring sites over time
#' @param det.method4.time An array containing the detectability values of method 4 for each species at monitoring sites over time
#' @param n.method The number of repeat visits to a site in a given survey year. Each element corresponds to each of the four detection methods
#' @param occ Occupancy raster stack with the number of layers equal to the number of species
#' @param det1.method1 An array that records simulated detection histories at sites using method 1
#' @param det1.method2 An array that records simulated detection histories at sites using method 2
#' @param det1.method3 An array that records simulated detection histories at sites using method 3
#' @param det1.method4 An array that records simulated detection histories at sites using method 4
#' @param stratify A raster layer of the environmental strata in which to randomly stratify monitoring sites. Each strata should have its own unique integer 
#' @param remote A raster layer identifying remote and non-remote areas. Remote areas should be identified with a 1, non-remote a 0
#' @param parks A raster layer identifying nested areas within the landscape in which to estimate power. Each sub-unit should be identified with a separate integer
#' @param Pr.dist A vector specifying the probability of a disturbance given the time since the last disturbance
#' @param disturbance A raster stack of the disturbance history durign a proceeding time period. Disturbned cells are givena  value of 1, undisturbed cells a value of 0.
#' @keywords 
#' @export
#' @examples

run.power <- function(effect.size, nsims, alpha, Tmax, s.years, trend, sites, model.disturbance, disturbance.type, species.list, dist.occ.time, dist.occ.freq, 
                      park.level, xy.sites, R, two.tailed, plot.sites, n.sites, n.species, n.park, all.loaded.sites, load.sites, new.site.selection, occ.time, 
                      det.method1.time, det.method2.time, det.method3.time, det.method4.time, n.method, occ, det.method1, det.method2, det.method3, det.method4, 
                      stratify, remote, parks, Pr.dist, disturbance) {
  
  if (n.method[1] == 0) {n.method[1] <- 1}
  if (n.method[2] == 0) {n.method[2] <- 1}
  if (n.method[3] == 0) {n.method[3] <- 1}
  if (n.method[4] == 0) {n.method[4] <- 1}
  
  det1.method1 <- array(NA, dim=c(n.sites,n.species,n.method[1],Tmax))
  det1.method2 <- array(NA, dim=c(n.sites,n.species,n.method[2],Tmax))
  det1.method3 <- array(NA, dim=c(n.sites,n.species,n.method[3],Tmax)) 
  det1.method4 <- array(NA, dim=c(n.sites,n.species,n.method[4],Tmax))
  det.method1.new <- det.method2.new <- det.method3.new <- det.method4.new <- occ.new <- array(0, dim=c(n.sites,n.species,Tmax))
  combined.effect<- matrix(NA,nsims,n.species)
  fr.hist <- time.hist <-veg.ID <- park.ID <- rep(NA,n.sites)
  time.fire <- array(NA, dim=c(n.sites,Tmax))
  fire.freq <- array(NA, dim=c(n.sites,Tmax))
  
  #Set up vectors and matrices to record the number of instances when we detect a change in occupancy
  powcnt <-rep(0,n.species) 
  pow.park <- matrix(0,ncol=n.species,nrow=length(n.park)) 
  fail <- rep(0,n.species) 
  fail.park <- matrix(0,ncol=n.species,nrow=length(n.park)) 
  
  for (ii in 1:nsims){ #Loop through number of simulations from 1 to nsims
    cat('\n',"Effect size = ", effect.size, " Simulation = ", ii) 
    
    #Select survey sites
    if (all.loaded.sites == FALSE & new.site.selection != "maxocc" & new.site.selection != "manual") {
      cat('\n',"Selecting survey sites.....") 
      xy.sites <- select.sites(sites=sites, n.sites=n.sites, R=R, all.loaded.sites=all.loaded.sites, load.sites=load.sites, new.site.selection=new.site.selection, plot.sites=plot.sites) #Returns matrix of XY coordinates plus description for which park each point belongs in 	
      
      #Extract occ and det values from selected sites
      cat('\n',"Extracting values from selected sites.....") 
      occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)]
      
      for (ss in 1:n.species){
        if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]}
        if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]}
        if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]}
        if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]}
      }
    }
    
    #Extract park ID, veg type and fire histories from selected sites
    park.ID <- parks[cellFromXY(parks, xy.sites)]
    xy.sites <- cbind(xy.sites,park.ID)
    
    #Make a copy of the occ and det values at each site so they can be manipulated by fire and the effect size over time
    det.method1.new[] <- det.method1.time
    det.method2.new[] <- det.method2.time
    det.method3.new[] <- det.method3.time
    det.method4.new[] <- det.method4.time
    occ.new[] <- occ.time
    
    #Model fire at monitoring sites if the model.fire == TRUE
    if (model.disturbance == TRUE & disturbance.type == "stochastic") { #If we specified a stochastic disturbance event
      cat('\n',"Modelling stochastic disturbance and updating occ and det layers.....")
      fire.ID <- matrix(NA, nrow = n.sites, ncol=nlayers(fire.hist) + Tmax)
      veg.ID <- stratify[cellFromXY(stratify, xy.sites)]
      layers <- covariates[cellFromXY(covariates, xy.sites)]
      fire.ID[,1:nlayers(fire.hist)] <- fire.hist[cellFromXY(fire.hist, xy.sites)]
      
      time.hist <- apply(fire.ID[,1:nlayers(fire.hist)], 1, function(x) (nlayers(fire.hist)+1)-max(which(x==1))) #time.hist <- round(fire.hist[cellFromXY(fire.hist, xy.sites)]) #time.fire <- 
      fr.hist <- rowSums(fire.ID[,1:nlayers(fire.hist)]) #fire.freq <- rowSums(fire.ID[,1:nlayers(fire)])
      
      for (jj in 1:Tmax) { #For each simulation loop ii, loop through time from 1 to Tmax
        #burn <- fire.point.model(time.hist, time.fire, veg.ID, jj) #Simulate whether monitoring sites burn
        burn <- stochastic.disturbance(time.hist, Pr.dist, n.sites)
        fire.ID[,nlayers(fire.hist) + jj] <- burn #Record fire history at sites for year jj
        if (jj==1) {
          time.fire[,jj] <- ifelse(burn == 1, 1, time.hist+1)
          fire.freq[,jj] <- rowSums(fire.ID[, (jj+1):(jj+nlayers(fire.hist))]) #Sum number of fires in 15 year moving window
        } 
        if (jj>1) {
          time.fire[,jj] <- ifelse(burn == 1, 1, time.fire[,jj-1]+1)
          fire.freq[,jj] <- rowSums(fire.ID[, (jj+1):(jj+nlayers(fire.hist))]) #Sum number of fires in 15 year moving window
        }	
        if (jj %in% s.years) { 
          occ.new <- refit.occ(occ.new, layers, time.fire, fire.freq, jj, dist.occ.time, dist.occ.freq, n.species) 
          det.out <- refit.det(det.method1.new=det.method1.new, det.method2.new=det.method2.new, det.method3.new=det.method3.new, det.method4.new=det.method4.new, layers, time.fire, fire.freq, jj, dist.occ.time, dist.occ.freq, n.species)
        }
      }
      det.method1.new[] <- det.out[[1]]
      det.method2.new[] <- det.out[[2]]
      det.method3.new[] <- det.out[[3]]
      det.method4.new[] <- det.out[[4]]
    }
    
    if (model.disturbance == TRUE & disturbance.type == "deterministic") { #If we specified a stochastic disturbance event
      dist.prop <- matrix(0,n.sites,Tmax)
      for (i in 1:Tmax) {
        dist.prop[,i] <- disturbance[[i]][cellFromXY(disturbance[[i]], xy.sites)]
      } 
      dist.copy <- dist.prop
      if (any(dist.prop<0)) { #If the determinstic disturbance has a negative effect on occupancy
        dist.copy[dist.copy<0] <- dist.copy[dist.copy<0]*-1
        for (i in 1:Tmax) {occ.new[,,i] <- occ.new[,,i]*dist.copy[,i]}
      }
      
      if (all(dist.prop>0)) { #If the deterministic disturbance has a positive effect on occupancy
        dist.copy[dist.copy==1] <- 0
        for (i in 1:Tmax) {occ.new[,,i] <- occ.new[,,i] + (1-occ.new[,,i])*dist.copy[,i]}
      }
    } 
    
    if (trend == "decreasing") { #Model decreasing trend
      effect.time <- 1-(effect.size/Tmax*c(1:Tmax))
      for (i in 1:Tmax) {occ.new[,,i] <- occ.new[,,i]*effect.time[i]}} 
    if (trend == "increasing") { #Model increasing trend
      effect.time <- effect.size/Tmax*c(1:Tmax)
      for (i in 1:Tmax) {occ.new[,,i] <- occ.new[,,i] + ((1-occ.new[,,i])*effect.time[i])}} 
    
    #Calculate the combined effect size
    magnitude.change <- -(occ.time[1,,1]-occ.new[1,,Tmax])
    
    if (any(magnitude.change<0)) {
      combined.effect[ii,which(magnitude.change<0)] <- (occ.time[1,which(magnitude.change<0),1]-occ.new[1,which(magnitude.change<0),Tmax])/occ.time[1,which(magnitude.change<0),1]
    }
    if (any(magnitude.change>0)) {
      combined.effect[ii,which(magnitude.change>0)] <- -((occ.new[1,which(magnitude.change>0),Tmax]) - occ.time[1,which(magnitude.change>0),1])/(1 - occ.time[1,which(magnitude.change>0),1])
    }
   
    for (i in 1:Tmax) {
      occ.new[,,i] <- ifelse(runif(ncell(occ.new[,,i])) < occ.new[,,i],1,0)
    }
    
    cat('\n',"Simulating monitoring at sites.....")
    for (jj in 1:n.method[1]){
      for (tt in 1:Tmax) { 
        det1.method1[,,jj,tt] <- ifelse(runif(ncell(det.method1.new[,,tt])) < det.method1.new[,,tt] & occ.new[,,tt] == 1, 1, 0)
      }}
    
    for (jj in 1:n.method[2]){
      for (tt in 1:Tmax) { 
        det1.method2[,,jj,tt] <- ifelse(runif(ncell(det.method2.new[,,tt])) < det.method2.new[,,tt] & occ.new[,,tt] == 1,1,0)
      }
    }
    
    for (jj in 1:n.method[3]){
      for (tt in 1:Tmax) { 
        det1.method3[,,jj,tt] <- ifelse(runif(ncell(det.method3.new[,,tt])) < det.method3.new[,,tt] & occ.new[,,tt] == 1 ,1,0)
      }
    }
    
    for (jj in 1:n.method[4]){
      for (tt in 1:Tmax) { 
        det1.method4[,,jj,tt] <- ifelse(runif(ncell(det.method4.new[,,tt])) < det.method4.new[,,tt] & occ.new[,,tt] == 1 ,1,0)
      }
    }
    
    #Combine the results of each detection method for each species, create detection histories and fit occupancy model 
    cat('\n',"Fitting occ model.....")
    value <- sig.test(two.tailed, alpha)
    methods <- list(det1.method1,det1.method2,det1.method3,det1.method4)
    
    for (ss in 1:n.species) {
      
      #Fit occ model for species detected using method 1
      if (sum(species.list[ss,-1]) == 1) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.1method(method=methods[[pos]], repeats=n.method[pos], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      #Fit occ model for species detected using 2 methods
      if (sum(species.list[ss,-1]) == 2) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.2method(method1=methods[[pos[1]]], method2=methods[[pos[2]]], repeats1=n.method[pos[1]], repeats2=n.method[pos[2]], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      #Fit occ model for species detected using 3 methods
      if (sum(species.list[ss,-1]) == 3) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.3method(method1=methods[[pos[1]]], method2=methods[[pos[2]]], method3=methods[[pos[3]]], repeats1=n.method[pos[1]], repeats2=n.method[pos[2]], repeats3=n.method[pos[3]], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      #Fit occ model for species detected using 4 methods
      if (sum(species.list[ss,-1]) == 4) {
        pos <- which(species.list[ss,-1]==1)
        occ.fit <- fit.occ.4method(method1=methods[[pos[1]]], method2=methods[[pos[2]]], method3=methods[[pos[3]]], method4=methods[[pos[4]]], repeats1=n.method[pos[1]], repeats2=n.method[pos[2]], repeats3=n.method[pos[3]], repeats4=n.method[pos[4]], s.years, n.sites, xy.sites, park.ID, park.level, powcnt, fail, pow.park, fail.park, value, ss, two.tailed, n.park)
      }
      
      powcnt <- occ.fit[[1]]
      fail <- occ.fit[[2]]
      pow.park <- occ.fit[[3]]
      fail.park <- occ.fit[[4]]
    } #End species loop
    
  } #End sims loop (ii)
  
  #Calculate pwr from simulations for each species
  combined.effect <- colMeans(combined.effect)
  output <- rbind(powcnt,fail,pow.park,fail.park,combined.effect) #If running in parallel
  return(output)
}

#' Plot results of power analysis
#'
#' This function plots results of the power analysis. It first unpacks the list generated by the foreach function and estimates power given the number of simulations 
#' and the number of cores used. Power is plotted for each species on a single figure. Users can then scroll through separate plots for power at a park level if it park.power = TRUE
#' @param pwr An array containing the results of the power analysis
#' @param n.species The number of species in which to estimate power
#' @param nsims The number of simulations across each core
#' @param effect.size A vector specifying the proportional decline in occupancy
#' @param n.park The number of nested management units in which to estimate power
#' @param species.list A vector specifying the name of each species analysed
#' @param n.cores A matrix specifying the name of each species analysed and the relevant detection method
#' @keywords 
#' @export
#' @examples

plot.power <- function(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) {
  Results <- array(dim=c(length(n.park)+1,n.species,length(effect.size)))
  new.effect <- matrix(NA, length(effect.size), n.species)
  pwr_all <- pwr
  for (i in 1:length(effect.size)) {
    pwr1 <- pwr_all[,i]
    dim(pwr1) <- c(length(pwr1)/n.species, n.species)
    results <- matrix(NA,length(n.park)+1,n.species)
    results[1,] <- pwr1[1,]/((nsims*n.cores)-pwr1[2,])
    results[2:(length(n.park)+1),] <- pwr1[3:(length(n.park)+2),]/((nsims*n.cores)-pwr1[(length(n.park)+3):(length(n.park)*2+2),])
    new.effect[i,] <- pwr1[nrow(pwr1),]/n.cores
    Results[,,i] <- results
  }
  
  if (park.level==TRUE) { 
    for (v in 1:(length(n.park)+1)) {
      cl <- rainbow(n.species)
      par(mfcol=c(1,1), mar=c(0.1,0.1,0.1,0.1), oma=c(4,4,4,4), mai=c(0.1,0.1,0.1,1),xpd=NA)
      plot(Results[v,1,]~new.effect[,1], type="l",ylim=c(0,1), xlim=c(-1,1), lwd=1, main="", col="grey", ylab="Statistical power", xlab="Effect size",cex.lab=1.2)
      for (i in 1:n.species) {lines(Results[v,i,]~new.effect[,i], type="l",ylim=c(0,1), xlim=c(-1,1), lwd=2, col=cl[i], lty=1)}
      legend("topright", c(as.character(species.list[,1])), inset=c(-0.5,0), lwd=rep(2,n.species), cex=0.7, col=cl[1:n.species])
      mtext("Increasing", side = 1, outer=TRUE, adj=0.15, line = 1.3, cex=0.9)
      mtext("Decreasing", side = 1, outer=TRUE, adj=0.60, line = 1.3, cex=0.9)
      if (v==1) {mtext("Power: landscape level", side = 3, outer=TRUE, adj=0.78, line = 1, cex=1.3)} 
      else {mtext(paste("Power: Park ", v-1), side = 3, outer=TRUE, adj=0.78, line = 1, cex=1.3)}
      readline(prompt="Press [enter] to continue.....")
    }
  } else {
    cl <- rainbow(n.species)
    par(mfcol=c(1,1), mar=c(0.1,0.1,0.1,0.1), oma=c(4,4,4,4), mai=c(0.1,0.1,0.1,1),xpd=NA)
    plot(Results[1,1,]~new.effect[,1], type="l",ylim=c(0,1), xlim=c(0,1), lwd=2, main="", col="grey", ylab="Statistical power", xlab="Effect size",cex.lab=1.2)
    mtext("Landscape level (all sites)", side = 3, outer=TRUE, adj=0.78, line = 1, cex=1.3)
    for (i in 1:n.species) {lines(Results[1,i,]~new.effect[,i], type="l",ylim=c(0,1), xlim=c(0,1), lwd=2, col=cl[i], lty=1)}
    legend("topright", c(as.character(species.list[,1])), inset=c(-0.5,0), lwd=rep(2,n.species), cex=0.7, col=cl[1:n.species])
    mtext("Increasing", side = 1, outer=TRUE, adj=0.15, line = 1.3, cex=0.9)
    mtext("Decreasing", side = 1, outer=TRUE, adj=0.60, line = 1.3, cex=0.9)
  }
  
  cat('\n',"##########################################################") 
  cat('\n',"OUTPUT SUMMARY.....") 
  cat('\n',"##########################################################")  
  cat('\n',"Number of species = ", n.species)
  cat('\n',"Number of parks = ", length(unique(n.park)))
  cat('\n',"Disturbance modelled at sites = ", model.disturbance)
  cat('\n',"Disturbance type = ", disturbance.type)
  cat('\n',"Analysis based on loaded sites = ", load.sites)
  cat('\n',"All loaded sites monitored = ", all.loaded.sites)
  cat('\n',"Number of sites monitored = ", n.sites)
  cat('\n',"Ratio of remote to non-remote sites = ", R)
  cat('\n',"Time horizon = ", Tmax, " years")
  cat('\n',"Survey years = ", s.years)
  cat('\n',"Number of repeat visits (method 1) = ", n.method[1])
  cat('\n',"Number of repeat visits (method 2) = ", n.method[2])
  cat('\n',"Number of repeat visits (method 3) = ", n.method[3])
  cat('\n',"Number of repeat visits (method 4) = ", n.method[4])
  cat('\n',"Direction of trend = ", trend) 
  cat('\n',"Effect size(s) = ", effect.size) 
  cat('\n',"Two tailed test = ", two.tailed)
  cat('\n',"Type I error rate = ", alpha)
  cat('\n',"Number of simulations = ", nsims*n.cores)
  
  return(Results) 
}  





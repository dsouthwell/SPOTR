
###################################################################################################################
#SPATIAL POWER ANALYSIS MODELLING IN R (SPOTR)
###################################################################################################################

# Darren Southwell
# Email: darren.southwell@unimelb.edu.au

# This package calculates the statistical power to detect simulated occupancy trends for multiple species in spatially explicit landscapes
# The package requires an occupancy raster layer for each species loaded into R as a raster stack.
# Simulations also require detectability raster layers for each species. Up to four independent detection methods can be simulated for each species separately.
# Users must specify the length of the monitoring program, the years in which monitoring occurs, and the number of repeat visits to sites for each detection method.
# Pre-determined monitoring sites can be loaded into the package for analysis.
# Alternatively, monitoring sites can be selected randomly at the start of each simulation, be stratified across environmental layers, or be placed on cells with the highest relative species richness.
# Occupancy and detectability raster layers can remain static over time or be dynamic in response to deterministic or stochastic disturbance events.
# To model a deterministic disturbance, users must load in a time series of raster layers specifying where the disturbance will occur and its proportional effect on occupancy
# To model a stochastic disturbance, users must load in a history of disturbances, a function specifying the relationship between a disturbance event and time since a disturbance, and functions relating occupancy and detectability with time since a disturbance for each species.
# Users simulate either an increasing or decreasing trend in occupancy, the magnitude of which is given by the effect size.
# Simulations can be run in series or in parallel. Users choose the alpha significance level and either a one-tailed or two-tailed significance test
# After running the simulations, the package returns the proportion of times (i.e. statistical power) a significant trend in occupancy was detected in the simulated detection histories
# Statistical power can be estimated across the landscape, and/or within smaller-level, nested management units

###################################################################################################################
#STEP 1 - Load package and set working directory
###################################################################################################################

rm(list=ls()) #Remove any existing files in the global environment
require(devtools)
install_github("dsouthwell/SPOTR") #Install the SPOTR package from GitHub
require(SPOTR) #Loading this package will automatically load the other required packages
#setwd("~/workspace/") #Set the working directory if estimating power for your own example

###################################################################################################################
#STEP 2 - Load in occupancy and detectability raster layers for each species
###################################################################################################################

# Simulations require raster layers of occupancy and detectability to be loaded into R for each species. These raster layers need to be built beforehand and
# loaded into R as raster stacks. All cells should be between 0 and 1. Cells not considered for monitoring should be set to NA.
# Occupancy values can be constant or vary across space in response to environmental covariates
# The package contains sample occupancy raster layers for 10 species. To view these, type 'occ', or load in your own below from your designated workspace.

#occ <- stack("occ.tif") #Load in species occupancy raster stack
#plot(occ[[1]]) #Plot species 1 in occupancy raster stack

# Simulations also require detectability raster layers for each species. All cells should be between 0 and 1. Detectability layers should have the same dimensions as the occupancy layers.
# Detectability values can be constant or vary across space in response to environmental covariates.
# Any combination of 4 separate detection methods are allowed for each species. For example, species A might be detected with methods 1 and 2, while species B might be detected with methods 3 and 4, and so fourth.
# Importantly, 4 raster stacks (one for each detection method) must be loaded into R regardless of whether they are all used to detect each species.
# If a particular method does not apply to a species, all cells in the correspnding raster layer should be set to zero.
# The package contains example detectability raster layers for 10 species. Type det.method1 to view the first the detectability estimates for the first method for each species. If all cells equal zero, then method 1 is not used to detect the species

#det.method1 <- stack("Method1.tif") #Load in species detectability layers for method 1
#det.method2 <- stack("Method2.tif") #Load in species detectability layers for method 2
#det.method3 <- stack("Method3.tif") #Load in species detectability layers for method 3
#det.method4 <- stack("Method4.tif") #Load in species detectability layers for method 4

###################################################################################################################
#STEP 3 - Load a table specifying which detection method applies to each species
###################################################################################################################

# Simulations require that the user loads in a table (species.list) which specifies which detection method is relevant to each species.
# The rows in this table lists each species, the columns refer to each of the four detection methods. A '1' in the table means that detection method is relevant to the species, a '0' means that it is not relevant.
# Type 'species.list' to see the example table for the 10 species. The number of rows equals the number of species, the number of columns equals the number of detection methods.
# To run your own example, load in a 'species.list' table as a csv file.

#species.list <- read.csv("workspace/Species_list.csv",stringsAsFactors=FALSE) #Load in a table specifying which methods are relevant to each species

###################################################################################################################
#STEP 4 - Decide whether or not to model deterministic or stochastic disturbances
###################################################################################################################

# The package can simulate disturbances in the landscape over time
# To model disturbances, set model.disturbance <- TRUE, otherwise set model.disturbance <- FALSE

model.disturbance <- TRUE #Choose whether disturbances are simulated at monitoring sites

# Two types of disturbances can be simulated - 1) Determistic disturbances, or 2) Stochastic disturbances.
# For deterministic disturbances, the location of the disturbance and its effect on occupancy is known beforehand during each time step.
# For stochastic disturbances, the incidence of a disturbance is probabilistic and is simulated when the package is run, but the effect of a disturbance on occupancy and detectability is known.

disturbance.type <- "stochastic" #Specify a 'stochastic' or 'deterministic' disturbance

# For deterministic disturbances, a raster stack must be loaded into R, with a raster layer for each year of the simulation period specifying the location and effect of the disturbance.
# Cells not subject to the deterministic disturbance must have values of 1, while cells disturbed should contain a value specifying the proportional change in occupancy due to the disturbance
# To look at an example of a deterministic disturbance, type 'disturbance'. Or alterantively, load in your own deterministic disturbance raster stack

#disturbance <- stack("~/workspace/disturbance.tif") #Load in deterministic disturbance event
#plot(disturbance) #Plot the disturbance stack to see where the deterministic disturbance occurs

# For stochastic disturbances, the user must load in a raster stack of the disturbance history (fire.hist).
# The disturbance history stack should contain a raster layer for each preceding year. Cells should indicate whether a disturbance occurred (1) or not (0) in each of the preceding years.
# The package then uses this stack to calculate a 'time since a disturbance' layer and 'number of disturbances' layer for the preceding period
# Load in the disturbance history raster stack or type 'fire.hist' to look at an example.

#fire.hist <- stack("~/workspace/fire.hist.tif") #Load in disturbance history raster stack

# For stochastic disturbances, the user must also load a raster stack of the environmental covariates that influence occupancy and detectability (covariates).
# The covariates stack is used to update occupancy and detectability raster layers in response to stochastic disturbance events
# There can be any number of layers in the covariates stack - it depends on what is thought to influence species occupancy and detectability
# Load in a covariates stack below or type 'covariates' to look at an example.
# In the example covariates stack, the raster layers correspond to: 1=Elevation, 2=Distance to creek, 3=Time since last fire, 4=Fire frequency, 5=Maximum temperature,
# 6=Fire extent, 7=Fire patchiness, 8=Soil type, 9=Terrain ruggedness, 10=Annual rainfall, 11=Projected foliage cover

#covariates <- stack("~/workspace/covariates.tif") #Load in covariate raster stack

# For stochastic disturbances, a vector (Pr.dist) must be specified relating the probability of a disturbance with the time since the last disturbance.
# This is used to simulate the incidence of disturbances at monitoring sites given the disturbance history at sites
# This function should be a vector with length equal to the number of layers in the fire.hist stack. All values should lie between 0 and 1.
# Below, the first element specifies the probability of a disturbance in one year since a disturbance, the second element is the probability of a disturbance event after 2 years without a disturbance, and so forth.
# An example function is specified below.

Pr.dist <- c(0.48,0.53,0.42,0.40,0.39,0.26,0.32,0.25,0.28,0.30,0.26,0.21,0.22,0.26,0.05) #Probability of a disturbance for each time period since the last disturbance

# For stochastic disturbances, users must specify how each species will respond to disturbances in terms of occupancy and detectability.
# Users must load in a table (species.occ) that defines how species occupancy is related to environmental covariate raster stack.
# In this table, the species are listed by rows, and the environmental covariates are listed by columns. An extra column is added for the intercept term.
# The elements of the species.occ table therefore specifies the species-specific relationships to environmental covariates
# Note, the number of columns is equal to the number layers in the covariates raster layer (not inlcuding the intercept or species column). The order of covariates must also match the order of layers the covariate raster stack.
# If a particular covariate is not relevant to a species, then it should contain a 0 value.
# Load in a species.occ table or type 'species.occ' to look at an example.

#species.occ <- read.csv("workspace/species.occ.csv", header=TRUE,stringsAsFactors=FALSE) #Load in statistical models for each species

# The package must be told which layers are dynamic (i.e. are updated by simulating the stochastic disturbance events)
# Dist.occ.time identifies the 'time since disturbance' layer in the covariates raster stack.
# Dist.occ.freq identifies the 'number of disturbances' layer in the covariates raster stack.
# Specify them below or enter the following for the example dataset.

dist.occ.time <- 3 #Records that the time since distrubance layer is the 3rd layer in the covariates raster stack
dist.occ.freq <- 4 #Records that the number of disturbance layer is the 4rd layer in the covariates raster stack

# If simulating stochastic disturbances, a table (species.det) must also be specified that relates detectability with environmental covariates for each species.
# Note, this table has 4 intercept columns for the 4 possible detection methods.
# Once again, the number and order of environmental covariates should match the number of layers in the 'covariates' raster stack
# Load in the species.det table or type 'species.det' to look at an example.

#species.det <- read.csv("workspace/species.det.csv", header=TRUE,stringsAsFactors=FALSE)

###################################################################################################################
#STEP 4 - Define HOW MANY sites to monitor
###################################################################################################################

# Users can either: 1) load in existing monitoring sites; 2) select sites at random; 3) stratify sites across environmental layers; 4) position sites on cells with the highest relative species richness
# If assessing existing monitoring sites, set load.sites <- TRUE, FALSE otherwise.

load.sites <- TRUE #Set to TRUE to simulate monitoring at known, fixed sites in the landscape, FALSE otherwise.

# If users do wish to simulate monitoring at fixed sites, they must provide the XY coordinates of these sites in a csv file.
# Importantly the X coordinate should have the heading POINT_X and the Y coordinate should have the heading POINT_Y
# Type 'sites' to see an example.

#sites<-read.csv("workspace/sites.csv",head=TRUE, sep=",") #Load in the XY coordinates of fixed monitoring sites. Note, the column header should be POINT_X and POINT_Y.

# Users might wish to monitor only a subset of the loaded sites, or alternatively, all of the loaded sites plus others selected randomly throughout the landscape
# Users should set all.loaded.sites <- TRUE to ONLY monitor the loaded sites, or set all.loaded.sites <- FALSE if they wish to monitor more or less than what is specified.

all.loaded.sites <- FALSE #Set to TRUE to monitor all loaded sites, otherwise FALSE to selecte a subset at random during each simulation

#Specify the number of sites to simulate monitoring.
#Importantly, if load.sites <- TRUE and all.loaded.sites <- TRUE, then n.sites must equal the number of rows in the sites table.

n.sites <- 100 #Enter the number of sites you wish to simulate monitoring.

###################################################################################################################
#STEP 5 - Decide WHERE to monitor
###################################################################################################################

# Instead of simulating monitoring at fixed sites in the landscape, sites can be positioned randomly in a number of ways.
# To randomly sample within enviromental layers, users must load in a raster layer (stratify) specifying the environmental layers in which to sample.
# Each environmental layer should be identified with a different integer value.
# The package automatically distributes an equal number of sites randomly within each of the layers
# Type 'stratify' to look at an example.

#stratify <- raster("~/workspace/stratify.tif") #Load in environmental strata

# Alternatively, users can randomly select a proportion of sites in 'remote' and 'non-remote' areas.
# For example, a remote area might include all cells that exceed a certain distance from roads etc.
# In this raster layer, 1's must be used to indicate remote cells, 0s to indicate non-remote cells.
# Load in a remote layer or type 'remote' to look at an example...

#remote <- raster("~/workspace/remote.tif") #Load in remote layer

# Users can then define the proportion of sites that are randomly selected in remote or non-remote areas
# Users not wanting to sepearte sites in remote and non-remote areas can make all cells in the remote layer equal to 1, and set the ratio equal to 1.

R <- 1 #the ratio of remote to non-remote sites. For example, R=0.4 means that 40% of sites will be in remote areas and 60% will be in non-remote areas.

# If users have selected load.sites <- FALSE, they must now specify how to position the new sites in the landscape during each simulation.
# There are 3 options
    # 1) "random" - n.sites are selected at random at the start of each simulation
    # 2) "stratified" - an equal number of sites are randomly selected within each layer in the 'stratified' raster layer
    # 3) "maxocc" - n.sites are positioned on the cells with the highest relative species richness; that is, the highest summed occupancy values

new.site.selection <- "random" #Define how new sites are selected

# The package doesn't plot the location of sites during each step of the simulation.
# However, users can plot the site locations for the first simulation to check whether sites are positioned as expected.

plot.sites <- TRUE #Set to TRUE to plot sites at the start of the first simulation

###################################################################################################################
#Step 6: Define WHEN to monitor
###################################################################################################################

# Users can now specify the length, timing and intensity of monitoring
# First, the length of the monitoring program must be specified.

Tmax <- 10 #Define length of monitoring program

# The time steps in which monitoring occurs must also be specified.
# It is assumed that all sites are monitored each survey year
# Note, surveys must occur on the last year; that is, be equal to Tmax.
# In the example below, monitoring is simulated in the 1st, 5th and 10th years of a 10 year monitoring program

s.years <- c(1,5,Tmax) #the years when monitoring occurs

# Specify the number of repeat visits at a site during a survey year for each of the four detection methods
# The first element in the vector below (n.method) corresponds to to the number of repeat visits using detection method 1.
# The second element in the vector below (n.method) corresponds to to the number of repeat visits using detection method 2, and so forth.
# Note, the number of repeat visits can be equal to 0 or should exceed 1 (i.e. not be equal to 1).
# In the example below, monitoring is conducted for 3 days/nights using method 1, 3 days/nights using method 2, 3 days/nights using method 3 and 5 days/nights using method 4

n.method <- c(3,3,3,5) #Number of repeat visits to sites corresponding to methods 1,2,3,4

###################################################################################################################
#Step 7: Set up parameters for power analysis
###################################################################################################################

# Users must now decide whether calculate power just across the landscape or also within smaller-level, nested management units

park.level <- TRUE #Set to TRUE to estimate power at a landscape and within smaller nested units, FALSE to just estimate power at a landscape-level

# If power is to be calculated within smaller management units (park.level == TRUE), users must idenitify where these parks are within the landscape.
# A parks layer should be loaded into R as a raster layer. Cells within each park should be identified by a unique integers, starting at 1.
# Type 'parks' for an example.

#parks <- raster("~/workspace/parks.tif") #Load in parks layer

# There is an option to model either a decreasing or increasing trend in occupancy through time.
# Set trend <- "decreasing" to model a decreasing trend, or trend <- "increasing" to model an increasing trend

trend <- "decreasing" #Choose from "increasing" or "decreasing"

# A vector effect.size must be defined, which specifies the effect sizes to be simulated.
# All elements of the vector must be between 0 and 1. There is no restriction on the length of the vector (i.e. any number of effect sizes can be tested)
# If trend == decreasing, the effect sizes should be interpreted as the proportional reduction in initial occupancy of a cell
# If trend == increasing, the effect sizes should be intepreted as the proportional increase in the difference between initial occupancy of a cell and 1.
# In the below example, 4 effect sizes are simulated - 10%, 40%, 60% and 90%

effect.size <- c(0.1,0.4,0.6,0.9) #Decide on the proportional change in occupancy at Tmax. This can be a vector ranging from a small changes in occupancy to very large changes

# Users must set the significance level of the statistical test. Choose from 3 options -0.01, 0.05 or 0.1.

alpha <- 0.05 #Set significance level

# Decide on a one-tailed or two-tailed significance test. A one-tailed test should be selected if the direction of change in occupancy is known.

two.tailed <- TRUE #Set to TRUE for a two tailed test, FALSE for a one-tailed test

# Define the number of simulations from which to calculate statistical power. Power is the proportion of these simulations in which there is a significant trend in the trend parameter

nsims <- 5 #Define the number of simulations used to calculate power


###################################################################################################################
#Step 8: Create empty arrays
###################################################################################################################

# The code below creates empty arrays to store the results of simulations. This is done outside of the main loop to speed up computations.

n.species <- nlayers(occ) #Calculates the number of species
n.park <- unique(parks) #Calculates the number of parks
det.method1.time <- det.method2.time <- det.method3.time <- det.method4.time <- occ.time <- array(0, dim=c(n.sites,n.species,Tmax)) #Create empty arrays

###################################################################################################################
#Step 9: Select monitoring sites
###################################################################################################################

# This section selectes monitoring sites for the first simulation given the information provided above. Users can then check to see if sites are being positioned in the landscape as expected
# Occupancy and detectability values are then extracted from the sites for further processing.

xy.sites <- select.sites(sites, n.sites, R, all.loaded.sites, load.sites, new.site.selection, plot.sites) #Define monitoring sites

park.ID <- parks[cellFromXY(parks, xy.sites)] #Extract parks values at monitoring sites
xy.sites <- cbind(xy.sites,park.ID) #combined XY coordinates of sites with park values
occ.time[,,s.years] <- occ[cellFromXY(occ, xy.sites)] #Extract occupancy values from monitoring sites
for (ss in 1:n.species){
  if (species.list[ss,2] == 1) {det.method1.time[,ss,s.years] <- det.method1[[ss]][cellFromXY(det.method1[[ss]], xy.sites)]} # Extract detectability values from sites using method 1
  if (species.list[ss,3] == 1) {det.method2.time[,ss,s.years] <- det.method2[[ss]][cellFromXY(det.method2[[ss]], xy.sites)]} # Extract detectability values from sites using method 2
  if (species.list[ss,4] == 1) {det.method3.time[,ss,s.years] <- det.method3[[ss]][cellFromXY(det.method3[[ss]], xy.sites)]} # Extract detectability values from sites using method 3
  if (species.list[ss,5] == 1) {det.method4.time[,ss,s.years] <- det.method4[[ss]][cellFromXY(det.method4[[ss]], xy.sites)]} # Extract detectability values from sites using method 4
}

###################################################################################################################
#Step 10: Check input layers and tables and run simulations
###################################################################################################################

# The last step before running the power analysis is to check that all raster layers have the same dimensions.
# A warning message will be displayed if there is a mismatch between raster stacks/layers.
# Users should run the function below. If no warnings are returned, the simulations are ready to run.

check.inputs(occ) #This function checks for any mismatches bewteen raster layers. If nothing is returned, everything's OK

# This type of analysis can take a long time to run on a desktop depending on the size of the raster layers, the number of species, whether disturbances are modelled etc.
# Simulations can be sped up by running them across multiple cores.
# If running in parallel, select the number of cores to run. Set to 1 if running in series

n.cores <- detectCores()-1 #Define the number of cores to run in parallel. If you don't want to run in parallel, set to 1
cl <- makeCluster(n.cores) #initiate clusters
registerDoParallel(cl) #initiate clusters

# Now, run the simulations by copying this function into R.
# NOTE - TO RUN IN PARALLEL USING N.CORES, MAKE SURE ITS SET TO "%DOPAR%"
# IF RUNNING IN SERIES, SET TO %DO%"
# Also, you cannot track progress if running in parallel. If running in series, a text output will keep track of the simulations

start.time<-proc.time() #Record start time
pwr <- foreach(i = 1:n.cores,.combine='+',.packages=c("SPOTR"))%dopar%{ #Run the power analysis. Set %dopar% to run in parallel, %do% otherwise.
  sapply(effect.size, run.power,
         nsims=nsims,
         alpha=alpha,
         Tmax=Tmax,
         s.years=s.years,
         trend=trend,
         sites=sites,
         model.disturbance=model.disturbance,
         disturbance.type=disturbance.type,
         species.list=species.list,
         dist.occ.time=dist.occ.time,
         dist.occ.freq=dist.occ.freq,
         park.level=park.level,
         xy.sites=xy.sites,
         R=R,
         two.tailed=two.tailed,
         plot.sites=plot.sites,
         n.sites=n.sites,
         n.species=n.species,
         n.park=n.park,
         all.loaded.sites=all.loaded.sites,
         load.sites=load.sites,
         new.site.selection=new.site.selection,
         occ.time=occ.time,
         det.method1.time=det.method1.time,
         det.method2.time=det.method2.time,
         det.method3.time=det.method3.time,
         det.method4.time=det.method4.time,
         n.method=n.method,
         occ=occ,
         det.method1=det.method1,
         det.method2=det.method2,
         det.method3=det.method3,
         det.method4=det.method4,
         stratify=stratify,
         remote=remote,
         parks=parks,
         Pr.dist=Pr.dist,
         disturbance=disturbance)
}
stopCluster(cl)
time.elapsed<-proc.time()-start.time #Record end time
time.elapsed #Print elapsed time

###################################################################################################################
#Step 11: Plot results
###################################################################################################################

#After running the simulations, calculate power and plot with respect to the effect size

Results <- plot.power(pwr, n.species, nsims, effect.size, n.park, species.list, n.cores) #Calls a function that finds the results and calculates power

###################################################################################################################
#Step 12: Save results
###################################################################################################################

setwd("~/workspace/Results")
save(pwr, file = paste("Results",".RData", sep = ""))



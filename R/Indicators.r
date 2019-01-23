# Indicators

# Large Fish Indicator
# The proportion (by weight) of the community that is larger than a threshold length
# For NS the threshold length is 40cm and suggested Ecological Quality Objective
# is an LFI >= 0.3.
# For CS, an LFI of > 0.4 has been proposed (see Interpreting the large fish indicator
# for the Celtic Sea)

# Also, only sizes above 10g should be included (or some other min g)
# species is the row number of the species - not the name
LFI <- function(model, Lt=40, species=NULL, minw = 10, maxw = max(model$w), minl = NULL, maxl = NULL)
{

	# Sort out which species we are looking at
	all.species <- 1:nrow(model$param$species)
	if (is.null(species))
		species <- all.species
	if (!all(species %in% all.species))
		stop("Incorrectly specified species")

    nsp <- length(species)

    # minl is a single value, so is minw - expand into vector
    # Use either W range or L range - cannot use both
    if(!is.null(minl))
        minw <- model$param$species$a * minl ^ model$param$species$b
    else
        minw <- rep(minw,nsp)
    # same for maxw
    if(!is.null(maxl))
        maxw <- model$param$species$a * maxl ^ model$param$species$b
    else
        maxw <- rep(maxw,nsp)
    # From this point on, minw and maxw are vectors of length nsp

	# Calculate the proportion by weight of the community that is larger than threshold L
	# So for each species get the weight at threshold length
	# Sum up total weight greater than threshold weight
	# Compare to total weight of all size classes
	Wt <-	model$param$species$a * Lt ^ model$param$species$b
	# Size classes greater than the threshold size
	Wtmat <- aaply(Wt[species], 1, function(x) x <= model$w)
	Wminmat <- aaply(minw, 1, function(x) x < model$w)
    Wmaxmat <- aaply(maxw, 1, function(x) x > model$w)
    Wtmat <- Wtmat & Wminmat & Wmaxmat
    Wallmat <- Wminmat & Wmaxmat

    # Abundances in those size classes greater than threshold size
	N_greater_Wt <- sweep(model$N[,species,],c(2,3),Wtmat,"*")
	# Biomasses in those size classes greater than Wt
	B_greater_Wt <- rowSums(sweep(N_greater_Wt,3,model$w*model$dw,"*"),dims=2)
	# Total Biomass in size classes greater than Wt for the species selected
	total_B_greater_Wt <- apply(B_greater_Wt,1,sum)

	# Total Biomass in size classes > winc for species selected
	N_all <- sweep(model$N[,species,],c(2,3),Wallmat,"*")
	B_all <- rowSums(sweep(N_all,3,model$w*model$dw,"*"),dims=2)
    total_B <- apply(B_all,1,sum)
    # Get indicator
	LFI <- total_B_greater_Wt / total_B
	return(LFI)
}

#\item Slope of the community size spectrum.
#\item Mean weight of all individuals in the community.
#\item Mean maximum weight
#\item Proportion of large fish in the community, aka the Large Fish Indicator (LFI).
#\item Number of collapsed species.



#The slope of the community size spectrum was calculated as the slope
# obtained from performing a linear regression of log($totalN_w$) against log$w$
#where $totalN_w$ is the total abundance of individuals at size $w$.
# meantsteps is number of rows of to average over (not actual years or timesteps)
CommunitySlope <- function(model, meantsteps=NULL, plotfit=FALSE, species=NULL, minw=10, maxw = max(model$w), minl = NULL, maxl = NULL)
{

 	# Sort out which species we are looking at
	all.species <- 1:nrow(model$param$species)
	if (is.null(species))
		species <- all.species
	if (!all(species %in% all.species))
		stop("Incorrectly specified species")

    nsp <- length(species)

    # minl is a single value, so is minw - expand into vector
    # Use either W range or L range - cannot use both
    if(!is.null(minl))
        minw <- model$param$species$a * minl ^ model$param$species$b
    else
        minw <- rep(minw,nsp)
    # same for maxw
    if(!is.null(maxl))
        maxw <- model$param$species$a * maxl ^ model$param$species$b
    else
        maxw <- rep(maxw,nsp)
    # From this point on, minw and maxw are vectors of length nsp

    # matrix of sizes within range
	Wminmat <- aaply(minw, 1, function(x) x < model$w)
    Wmaxmat <- aaply(maxw, 1, function(x) x > model$w)
    Wallmat <- Wminmat & Wmaxmat

    # Not taking mean of last tsteps, just last timestep
    nt <- dim(model$N)[1]
    if (is.null(meantsteps))
        Nt <- model$N[nt,species,]
    else
        Nt <- apply(model$N[(nt-meantsteps+1):nt,species,,drop=FALSE],c(2,3),mean)

    # Only include sizes within range
	Nt <- Nt * Wallmat
    
    Nt <- apply(Nt,2,sum)
    # hack to remove 0
    #Nt[Nt<=0] <- 1e-9
    Nt[Nt<=0] <- NA

    fit <- lm(log(Nt) ~ log(model$w))
    if(plotfit)
    {
        plot(x=log(model$w),y=log(Nt),xlab="log w", ylab="log total abundance")
        abline(fit)
    }
    out <- c(slope = fit$coefficients[[2]], r2 = summary(fit)$r.squared)
    return(out)
}

# Mean weight of all individuals in the community
MeanWeight <- function(model, species=NULL, minw=10, maxw = max(model$w), minl = NULL, maxl = NULL)
{

#browser()
 	# Sort out which species we are looking at
	all.species <- 1:nrow(model$param$species)
	if (is.null(species))
		species <- all.species
	if (!all(species %in% all.species))
		stop("Incorrectly specified species")

    nsp <- length(species)

    # minl is a single value, so is minw - expand into vector
    # Use either W range or L range - cannot use both
    if(!is.null(minl))
        minw <- model$param$species$a * minl ^ model$param$species$b
    else
        minw <- rep(minw,nsp)
    # same for maxw
    if(!is.null(maxl))
        maxw <- model$param$species$a * maxl ^ model$param$species$b
    else
        maxw <- rep(maxw,nsp)
    # From this point on, minw and maxw are vectors of length nsp

    # matrix of sizes within range
	Wminmat <- aaply(minw, 1, function(x) x < model$w)
    Wmaxmat <- aaply(maxw, 1, function(x) x > model$w)
    Wallmat <- Wminmat & Wmaxmat
    #winc <- model$w > minw
 
    # Get total weight and abundance of of community
	N_all <- sweep(model$N[,species,],c(2,3),Wallmat,"*")
	B_all <- rowSums(sweep(N_all,3,model$w*model$dw,"*"),dims=2)
    total_B <- apply(B_all,1,sum)
    # total_N <- rowSums(apply(N_all,c(1,3),sum))
    total_N <- rowSums(sweep(apply(N_all,c(1,3),sum),2,model$dw,"*"))
    
    #  Btotal <- apply(sweep(model$N[,species,winc],3,model$w[winc]*model$dw[winc],"*"),1,sum)
    # Get total abundance in community
    #  Ntotal <- apply(sweep(model$N[,species,winc],3,model$dw[winc],"*"),1,sum)
    # Get mean weight
    return(total_B / total_N)
}

# Max weight taken to be Winf
# i.e. what is the mean Winf of the community?
MeanMaxWeight <- function(model, species=NULL, minw=10, maxw = max(model$w), minl = NULL, maxl = NULL)
{


  # sum of (Winf_i * total_B_i) / sum of total_B_i

 	# Sort out which species we are looking at
	all.species <- 1:nrow(model$param$species)
	if (is.null(species))
		species <- all.species
	if (!all(species %in% all.species))
		stop("Incorrectly specified species")

    nsp <- length(species)

    # minl is a single value, so is minw - expand into vector
    # Use either W range or L range - cannot use both
    if(!is.null(minl))
        minw <- model$param$species$a * minl ^ model$param$species$b
    else
        minw <- rep(minw,nsp)
    # same for maxw
    if(!is.null(maxl))
        maxw <- model$param$species$a * maxl ^ model$param$species$b
    else
        maxw <- rep(maxw,nsp)
    # From this point on, minw and maxw are vectors of length nsp

    # matrix of sizes within range
	Wminmat <- aaply(minw, 1, function(x) x < model$w)
    Wmaxmat <- aaply(maxw, 1, function(x) x > model$w)
    Wallmat <- Wminmat & Wmaxmat

    #  winc <- model$w > minw
	N_all <- sweep(model$N[,species,],c(2,3),Wallmat,"*")
	B_all <- rowSums(sweep(N_all,3,model$w*model$dw,"*"),dims=2)
#    WmaxB <- apply(sweep(B_all[,species],2,model$param$species$Winf[species],"*"),1,sum)
    WmaxB <- apply(sweep(B_all,2,model$param$species$Winf[species],"*"),1,sum)
    total_biomass <- apply(B_all,1,sum)

    #  total_biomass <- apply(sweep(model$N[,species,],3,model$w*model$dw,"*"),1,sum)
    #  total_biomass_species <- apply(sweep(model$N[,species,],3,model$w*model$dw,"*"),c(1,2),sum)
    #  total_biomass_species <- apply(sweep(model$N[,species,winc],3,model$w[winc]*model$dw[winc],"*"),c(1,2),sum)

    #  WmaxB <- apply(sweep(total_biomass_species,2,model$param$species$Winf[species],"*"),1,sum)
    #  total_biomass <- apply(sweep(model$N[,species,winc],3,model$w[winc]*model$dw[winc],"*"),1,sum)
  return(WmaxB / total_biomass)
}

# Collapsed species indicator
# Number of collapsed species and who...?
# Collapsed species is one whose biomass is less than 10% of unfished biomass.
# Ignore first x timesteps
Collapse <- function(model, species=NULL, ignoretime = 48, extinctratio = 0.1, unfishedbiomass)
{
  # Look through time series, if anyone goes below threshhold at any point = collapse
 	all.species <- 1:nrow(model$param$species)
	if (is.null(species))
		species <- all.species
	if (!all(species %in% all.species))
		stop("Incorrectly specified species")

  total_biomass_species <- apply(sweep(model$N[,species,],3,model$w*model$dw,"*"),c(1,2),sum)
  relbiomass <- sweep(total_biomass_species,2,unfishedbiomass,"/")
  collapse <- apply(relbiomass[-(1:ignoretime),],2, function(x) any(x < extinctratio))
  names(collapse) <- model$param$species$species[species]
  # When collapsed
  #  collapse <- apply(relbiomass[-(1:ignoretime),],2, function(x) min(which(x < extinctratio)))

  return(collapse)
}
# Or do we want who collapsed and when?


#*******************************************************************************
# Economic and fleet based indicators
#*******************************************************************************
# Net Present Value - NPV

# Calculates fishing mortality through time, by species, size and gear
FGear <- function(model)
{
  # Number of time steps we saved data at
  tsave <- dim(model$F)[1]
  # Fmortality by time, species, w and gear
  Fg <- array(NA,dim=c(tsave, model$param$nspp, length(model$w), ncol(model$param$Q)),
              dimnames = list(time = 1:tsave*model$param$isave*model$param$dt, species=model$param$species$species, size=signif(model$w,3), gear=dimnames(model$param$Q)$gear))
  # Have to be careful with number of steps saved not equalling number of effort steps - use isave
  # This for loop is quite manky - can we clear it up?
  for (i in 1:tsave)
    Fg[i,,,] <- sweep(model$selectivity,c(1,3),sweep(model$param$Q,2,model$effort[i*model$param$isave,],"*"),"*")
  return(Fg)
}

# Calculate yield is through time, by species, size and gear
YieldGear <- function(model)
{
  Fg <- FGear(model)
  Yg <- sweep(Fg,c(1,2,3),sweep(model$N,3,model$w*model$dw,"*"),"*")
  return(Yg)
}

# Calculate revenue through time, by species, size and gear
# For the moment assume that price is same for all gears and stationary through time
Revenue <- function(model,price)
{
  # Check price dims
  Yg <- YieldGear(model)
  r <- sweep(Yg,c(2,3),price, "*")
  return(r)
}

# fg <- Fgear(NSmodel)
# yg <- YieldGear(NSmodel)
# r <- Revenue(NSmodel,price)

# Calculates the discounted revenue by species and gear through time (not by size...)
# timerange is a vector of length 2 of time (in years) that you are calculating NPV over
DiscountedRevenue <- function(model,price,timerange=NULL,discount_rate=0.05)
{
  Rg <- Revenue(model,price)
  # total revenue by gear
  total_revenue <- apply(Rg,c(1,2,4),sum)  

  # window the team here
    if (!is.null(timerange))
      #total_revenue <- total_revenue[as.character(seq(from=timerange[1],to=timerange[2],by=model$param$isave)),,,drop=FALSE]
      total_revenue <- total_revenue[as.character(seq(from=timerange[1],to=timerange[2],by=model$param$isave * model$param$dt)),,,drop=FALSE]      

    time <- as.numeric(dimnames(total_revenue)$time)

    drfactor <- (1+discount_rate) ^ (time-1)
    # get the discounted revenue by species and gear through time
    discounted_revenue <- sweep(total_revenue,1,drfactor,"/")
    return(discounted_revenue)
}

#dr <- DiscountedRevenue(NSmodel,price)
# dr <- DiscountedRevenue(NSmodel,price, timerange=c(80,90))


# Net present value - discounted revenue summed over time
NPV <- function(model,price,timerange=NULL,discount_rate=0.05)
{
  dr <- DiscountedRevenue(model,price,timerange,discount_rate)
  npv <- apply(dr,c(2,3),sum)
  return(npv)
}

# NPV(NSmodel,price, timerange=c(80,90))
# NPV(NSmodel,price)

#-----------------------------------------------------------------------
# Get the slope at every point in time
CommunitySlopeTime <- function(model, species=NULL, minw=10, maxw = max(model$w), minl = NULL, maxl = NULL)
{

 	# Sort out which species we are looking at
	all.species <- 1:nrow(model$param$species)
	if (is.null(species))
		species <- all.species
	if (!all(species %in% all.species))
		stop("Incorrectly specified species")

    nsp <- length(species)

    # minl is a single value, so is minw - expand into vector
    # Use either W range or L range - cannot use both
    if(!is.null(minl))
        minw <- model$param$species$a * minl ^ model$param$species$b
    else
        minw <- rep(minw,nsp)
    # same for maxw
    if(!is.null(maxl))
        maxw <- model$param$species$a * maxl ^ model$param$species$b
    else
        maxw <- rep(maxw,nsp)
    # From this point on, minw and maxw are vectors of length nsp

    # matrix of sizes within range
    Wminmat <- aaply(minw, 1, function(x) x < model$w)
    Wmaxmat <- aaply(maxw, 1, function(x) x > model$w)
    Wallmat <- Wminmat & Wmaxmat

    # Not taking mean of last tsteps, just last timestep
    #nt <- dim(model$N)[1]
    #if (is.null(meantsteps))
    #    Nt <- model$N[nt,species,]
    #else
    #    Nt <- apply(model$N[(nt-meantsteps+1):nt,species,,drop=FALSE],c(2,3),mean)
    Nt <- model$N

    # Only include sizes within range
    #	Nt <- Nt * Wallmat
    Nt <- sweep(Nt[,species,],c(2,3),Wallmat,"*")
    # Sum over to get total N    
    Nt <- apply(Nt,c(1,3),sum)
    # hack to remove 0
    #Nt[Nt<=0] <- 1e-9
    Nt[Nt<=0] <- NA

    # finally, do the fits
    slope <- rep(NA,dim(Nt)[1])
    r2 <- rep(NA,dim(Nt)[1])
    for (i in 1:dim(Nt)[1]){
	fit <- lm(log(Nt[i,]) ~ log(model$w))
	slope[i] <- fit$coefficients[[2]]
	r2[i] = summary(fit)$r.squared
    }

    #fit <- lm(log(Nt) ~ log(model$w))
    #out <- c(slope = fit$coefficients[[2]], r2 = summary(fit)$r.squared)
    out <- data.frame(time = 1:dim(Nt)[1], slope = slope, r2 = r2)
    return(out)
}


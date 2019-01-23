# Summary funcs from size-based model

summarySpecies <- function(x, price = NA, baseBiomass = NA, baseSSB = NA, NPVTimeRange = 10, discount_rate=0.05,meanTimeSteps = 10,minw=0){
    #browser()
    dim1 <- dim(x$N)[1]
    # Biological measures
    Nhat <- apply(x$N[(dim1-meanTimeSteps+1):dim1,,,drop=FALSE],c(2,3),mean)
    SSB <- apply(sweep(x$psi * Nhat,2,x$w * x$dw,"*"),1,sum) /1e6
    # only include > minw in these calcs
    wInclude <- x$w >= minw
    biomass <- apply(sweep(Nhat,2,x$w * x$dw * wInclude,"*"),1,sum) /1e6
    relativeBiomass <- biomass / baseBiomass
    relativeSSB <- SSB / baseSSB
    meanRDD <- apply(x$RDD[(dim1 - meanTimeSteps+1):dim1,,drop=FALSE],2,mean)
    meanRDI <- apply(x$RDI[(dim1 - meanTimeSteps+1):dim1,,drop=FALSE],2,mean)
    relativeRDD <- meanRDD / x$param$species$R0
    relativeRDI <- meanRDI / x$param$species$R0

    # Economic measures
    yield <- apply(x$Yield[(dim(x$Yield)[1]-meanTimeSteps+1):dim(x$Yield)[1],,drop=FALSE],2,mean)/1e6
    #browser()
    revenue <- Revenue(x,price) # Revenue by time, species, size and gear
    meanRevenue <- apply(revenue[(dim(revenue)[1]-meanTimeSteps+1):dim(revenue)[1],,,,drop=FALSE],c(2,3,4),mean) # mean over time
    totalMeanRevenue <- apply(meanRevenue,1,sum) # total over gears and sizes
    NPVTimeRange <- c((x$param$tmax-NPVTimeRange),x$param$tmax)
    npv <- NPV(x,price,timerange=NPVTimeRange,discount_rate=discount_rate)

    op <- data.frame(species = x$param$species$species,
                        SSB = SSB,
                        biomass = biomass,
                        relativeBiomass = relativeBiomass,
                        relativeSSB = relativeSSB,
                        meanRDD = meanRDD,
                        meanRDI = meanRDI,
                        relativeRDD = relativeRDD,
                        relativeRDI = relativeRDI,
                        yield = yield,
                        revenue = as.numeric(totalMeanRevenue),
                        npv = as.numeric(npv[,1])
                        #effort = x$effort[1],
                        #q = x$param$Q[,1]
			)
    return(op)
}

summaryCommunity <- function(x, Lt = 40, price = NA, NPVTimeRange = 10, discount_rate=0.05,meanTimeSteps = 10, species = NULL, minw = 10, maxw = max(x$w), minl=NULL, maxl=NULL){
    #browser()

    dim1 <- dim(x$N)[1]
    #browser()
    # Biological measures
    lfi <- LFI(x, Lt=Lt, species=species, minw = minw, maxw = maxw, minl = minl, maxl = maxl)
    lfi <- mean(lfi[(length(lfi) - meanTimeSteps + 1) : length(lfi)])
    cs <- CommunitySlope(x, meantsteps=meanTimeSteps, species=species, minw=minw, maxw = maxw, minl = minl, maxl = maxl)["slope"]
    mw <- MeanWeight(x, species=species, minw=minw, maxw = maxw, minl = minl, maxl = maxl)
    mw <- mean(mw[(length(mw) - meanTimeSteps + 1) : length(mw)])
    mmw <- MeanMaxWeight(x, species=species, minw=minw, maxw = maxw, minl = minl, maxl = maxl)
    mmw <- mean(mmw[(length(mmw) - meanTimeSteps + 1) : length(mmw)])

    # Economic
    yield <- apply(x$Yield[(dim(x$Yield)[1]-meanTimeSteps+1):dim(x$Yield)[1],,drop=FALSE],2,mean)/1e6
    revenue <- Revenue(x,price) # Revenue by time, species, size and gear
    meanRevenue <- apply(revenue[(dim(revenue)[1]-meanTimeSteps+1):dim(revenue)[1],,,,drop=FALSE],c(2,3,4),mean) # mean over time
    totalMeanRevenue <- apply(meanRevenue,1,sum) # total over gears and sizes
    NPVTimeRange <- c((x$param$tmax-NPVTimeRange),x$param$tmax)
    npv <- NPV(x,price,timerange=NPVTimeRange,discount_rate=discount_rate)



    op <- data.frame(lfi = lfi, cs = cs, mw = mw, mmw = mmw,
                     yield = sum(yield), revenue = sum(meanRevenue), npv = sum(npv))
    return(op)
}

#-----------------------------------------------------------------------------
# Summarise whole time series


summarySpeciesTime <- function(x, baseBiomass = NA, baseSSB = NA, minw=0){
    # only include > minw in these calcs
    wInclude <- x$w >= minw

    # Biological measures
    #sweep(x$N,c(2,3),x$psi,"*")
    #sweep(sweep(x$N,c(2,3),x$psi,"*"), 3, x$w * x$dw, "*")
    SSB <- apply(sweep(sweep(x$N,c(2,3),x$psi,"*"), 3, x$w * x$dw * wInclude, "*"), c(1,2), sum) / 1e6
    biomass <- apply(sweep(x$N,3,x$w * x$dw * wInclude, "*"), c(1,2), sum) / 1e6
    relativeBiomass <- sweep(biomass,2, baseBiomass, "/")
    relativeSSB <- sweep(SSB,2, baseSSB, "/")
    RDD <- x$RDD
    RDI <- x$RDI
    relativeRDD <- sweep(RDD, 2, x$param$species$R0, "/")
    relativeRDI <- sweep(RDI, 2, x$param$species$R0, "/")

    # Economic measures
    #yield <- apply(x$Yield[(dim(x$Yield)[1]-meanTimeSteps+1):dim(x$Yield)[1],,drop=FALSE],2,mean)/1e6
    yield <- x$Yield

    # Leave economics alone for now

op <- data.frame(time = 1:dim(SSB)[1],
		 species = rep(x$param$species$species,each=dim(SSB)[1]),
		 SSB = c(SSB),
		 biomass = c(biomass),
		 relativeBiomass = c(relativeBiomass),
		 relativeSSB = c(relativeSSB),
		 RDD = c(RDD),
		 RDI = c(RDI),
		 relativeRDD = c(relativeRDD),
		 relativeRDI = c(relativeRDI),
		 yield = c(yield)
		 )
    return(op)
}


summaryCommunityTime <- function(x, Lt = 40, species = NULL, minw = 10, maxw = max(x$w), minl=NULL, maxl=NULL){
    # Biological measures
    lfi <- LFI(x, Lt=Lt, species=species, minw = minw, maxw = maxw, minl = minl, maxl = maxl)
    cs <- CommunitySlopeTime(x, species=species, minw=minw, maxw = maxw, minl = minl, maxl = maxl)["slope"]
    mw <- MeanWeight(x, species=species, minw=minw, maxw = maxw, minl = minl, maxl = maxl)
    mmw <- MeanMaxWeight(x, species=species, minw=minw, maxw = maxw, minl = minl, maxl = maxl)

    # Economic
    yield <- apply(x$Yield, 1, sum)

    #revenue <- Revenue(x,price) # Revenue by time, species, size and gear
    #meanRevenue <- apply(revenue[(dim(revenue)[1]-meanTimeSteps+1):dim(revenue)[1],,,,drop=FALSE],c(2,3,4),mean) # mean over time
    #totalMeanRevenue <- apply(meanRevenue,1,sum) # total over gears and sizes
    #NPVTimeRange <- c((x$param$tmax-NPVTimeRange),x$param$tmax)
    #npv <- NPV(x,price,timerange=NPVTimeRange,discount_rate=discount_rate)

    op <- data.frame(time = 1:dim(x$N)[1], lfi = lfi, cs = cs, mw = mw, mmw = mmw, yield = yield)
    return(op)
}



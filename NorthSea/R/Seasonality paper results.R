## R code for the multispecies size spectrum model and simulations in Datta & Blanchard 2016 Canadian Journal of Fisheries and Aquatic Sciences.
## This code is a modification of the original multispecies size spectrum model published by Blanchard et al. 2014 (Journal of Applied Ecology), which has been now generalised as an R package 'mizer' Scott, Blanchard & Andersen 2014 (mizer).
## Modifications to the model equations can be found in SizeBasedModel.R involved adding a seasonal plankton bloom and seasonal fish reproduction process

# Preliminary stuff
rm(list=ls())

# setwd("~/MyFolder/R code for size spectrum model (figshare)/NorthSea") # Mac
# setwd("C:/MyFolder/R code for size spectrum model (figshare)/NorthSea") # Windows

tic <- proc.time() # start timer

source("R/paramNorthSeaModel.r") # parameters for North Sea
source("../R/SizeBasedModel.r") # functions for model equations integrating over time
source("../R/SelectivityFuncs.r") # functions for fishing parameters used
source("../R/plots.r") # code for making different plots
source("../R/Indicators.r")
source("../R/summaryFuncs.r")

library(plyr) 

#------------------------------------------------------------------------------
# Set up the calibrated base model from Blanchard et al. 2014
fileSpecies <- "input/nsea_params.csv"
load("input/interactionmatrix_Schoener_twostage2D.RData")
param <- paramNSModel(fileSpecies=fileSpecies,
                      theta = theta,
                      kap=1e11)

param$tmax <- 500 # making sure what time is in years (500 normal, 150 for weekly timestep)
param$dt = 1/12 # time step (years) (1/12 for monthly, 1/52 for weekly)
juliabase <- Setup(param) # set up parameters and functions
ss_frac = 0.7 # proportion of resource spectrum present outside of peak
p_peak_s = 30 # severity of peak (30 lasts around 50 weeks)
p_peakt = 0.4 # proportion of year where bloom peak occurs  
# R species order: Sprat, Sandeel, N.pout, Herring, Dab, Whiting, Sole, Grey Gurnard, Plaice, Haddock, Cod, Saithe
# w_inf order: Sprat, Sandeel, N.pout, Dab, Herring, Grey Gurnard, Sole, Whiting, Plaice, Haddock, Cod, Saithe
s_order = c(1:3, 5, 4, 8, 7, 6, 9:12)
r_peak_s = c(3.6047, 2.9994, 1.944, 0.40493, 1.141, 1.4257, 4.973, 0.795, 4.0495, 3.6567, 5.4732, 1.951) # reproduction peakinesses
r_peakt = c(0.54765, 0.064324, 0.3574, 0.85759, 0.38245, 0.37164, 0.48564, 0.50213, 0.22454, 0.41806, 0.31236, 0.33333) # time of peak

#### SCENARIOS TO SIMULATE

# (1) Now set up and project the (non-seasonal) base case model
# no seasonality
r_peak = array(0, 12)
p_peak = 0
baseModel <- Project(juliabase)

# (2) ADD SEASONALITY
r_peak = r_peak_s
p_peak = p_peak_s
baseModel2 <- Project(juliabase)

# (3) check effects of seasonal plankton alone (not specifically needed for plots)
# plankton blooms, no reproduction peaks
# r_peak = array(0, 12)
# p_peak = p_peak_s
# baseModel3 <- Project(juliabase)

# (4) check effects of seasonal spawning alone (not specifically needed for plots)
# plankton blooms, no reproduction peaks
# r_peak = r_peak_s
# p_peak = 0
# baseModel4 <- Project(juliabase)
#
#
#
# FISHING

# Time-varying fishing mortality scenario projections

# Set up parameters based on base model but with correct number of time steps
fishParam <- baseModel$param
fishParam$tmax <- 50
fishModel <- Setup(fishParam,ContinueCalculation=TRUE,initialcommunity=baseModel)
endt = length(fishModel$effort[,1])
# Set the new effort - double this after 12 months
fishModel$effort[((1/param$dt)+1):endt,] <- fishModel$effort[((1/param$dt)+1):endt,]*2
r_peak = array(0, 12)
p_peak = 0
fishModel <- Project(fishModel)
# plotResults(fishModel)

# Do same as above with seasonal model
fishModel2 <- Setup(fishParam,ContinueCalculation=TRUE,initialcommunity=baseModel2)
# Set the new effort - double this after 12 months
fishModel2$effort[((1/param$dt)+1):endt,] <- fishModel2$effort[((1/param$dt)+1):endt,]*2
r_peak = r_peak_s
p_peak = p_peak_s
fishModel2 <- Project(fishModel2)
# plotResults(fishModel2)

## Halving fishing effort for non-seasonal system
unfishModel <- Setup(fishParam,ContinueCalculation=TRUE,initialcommunity=baseModel)
# Set the new effort - reduce this after 12 months
unfishModel$effort[((1/param$dt)+1):endt,] <- unfishModel$effort[((1/param$dt)+1):endt,]*0.5
r_peak = array(0, 12)
p_peak = 0
unfishModel <- Project(unfishModel)
# plotResults(unfishModel)

# Do same as above with seasonal model
unfishModel2 <- Setup(fishParam,ContinueCalculation=TRUE,initialcommunity=baseModel2)
# Set the new effort - reduce this after 12 months
unfishModel2$effort[((1/param$dt)+1):endt,] <-unfishModel2$effort[((1/param$dt)+1):endt,]*0.5
r_peak = r_peak_s
p_peak = p_peak_s
unfishModel2 <- Project(unfishModel2)
# plotResults(unfishModel2)
#
#
# # SAVE model results for later (if wanted - warning, files will be >100MB)
# save(baseModel,file="output/baseModel.RData")
# save(baseModel2,file="output/baseModel2.RData")
# save(baseModel3,file="output/baseModel3.RData")
# save(baseModel4,file="output/baseModel4.RData")
# save(fishModel,file="output/fishModel.RData")
# save(fishModel2,file="output/fishModel2.RData")
# save(unfishModel,file="output/unfishModel.RData")
# save(unfishModel2,file="output/unfishModel2.RData")
# 

##---- ONCE RESULTS ARE DONE, COMMENT OUT LINES 43 TO 122 AND RUN FOLLOWING 8 LINES  
# load(file="output/baseModel.RData")
# load(file="output/baseModel2.RData")
# load(file="output/baseModel3.RData")
# load(file="output/baseModel4.RData")
# load(file="output/fishModel.RData")
# load(file="output/fishModel2.RData")
# load(file="output/unfishModel.RData")
# load(file="output/unfishModel2.RData")

# SHORTHAND FOR CONSTANTS
laststep = length(baseModel$N[,1,1])
dt = param$dt # define dt
months = laststep - (1/dt) + round(seq(1,12)/(12*dt)) # indices in a year for monthly snapshots
weeks = laststep - (1/dt) + round(seq(1,52)/(52*dt)) # indices in a year for monthly snapshots
w = baseModel$w
growth_st = baseModel$gg
pred_st = baseModel$M2
natmort_st = baseModel$M2background
growth = baseModel2$gg
pred = baseModel2$M2
natmort = baseModel2$M2background
cspectrum<-apply(baseModel$N,c(1,3),sum)
cspectrum2<-apply(baseModel2$N,c(1,3),sum)

# Figures for paper

ts<-seq(param$dt,param$tmax,param$dt)

to_save = 2 # 0 for one big PDF, 1 for separate EPS files, else just in plot window

# Figure 1 - figure of spawning periods you made

# save figures in a pdf file to look at later

if (to_save == 0) pdf("paperfigs.pdf")


# Figure 2 - Effects of seasonality on the community size spectrum - not sure we need this one? is the  species by species one  better?

par(mfrow=c(1,1),mai=c(0,0,0,0),omi=c(1,1,1,1))

# (a) Non-seasonal and seasonal community spectrum
if (to_save == 1){
  setEPS()
  postscript("communities.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
plot(w,cspectrum2[laststep,],log="xy",type="l",lwd=4,ylab="Number density",xlab="Body mass (g)",xlim=c(1e-4,1e5),ylim=c(1e-5,1e18),col="grey",
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
points(w,cspectrum[laststep,],type="l",lty=2,lwd=4)
points(baseModel$wFull,baseModel$nPP[laststep,],type="l",lty=3,lwd=2)
if (to_save == 1) dev.off()

# (b) Seasonal speciues spectra
if (to_save == 1){
  setEPS()
  postscript("comm_seasonal.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
plot(w,cspectrum2[laststep,],log="xy",typ="l",lwd=3,ylab="Number density",xlab="Body mass (g)",ylim=c(1e-5,1e15),
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
for(i in 1:12) {
  lines(w,baseModel2$N[laststep,i,],log="xy",lty=2,xlim=c(1e-5,1e5),ylim=c(1e-3,1e25),col="grey")
  bla = which(baseModel2$N[laststep,i,] < 1e-6)
  points(w[bla[1]-1],baseModel2$N[laststep,i,bla[1]-1],pch=15, col="black")
}
if (to_save == 1) dev.off()


# Figure 3 - Seasonal species spectra relative to non-seasonal ones

if (to_save == 1){
  setEPS()
  postscript("species_scaled.eps", horizontal = FALSE, onefile = FALSE, paper = "special", height = 10, width = 8)
}
par(mfrow=c(4,3),mai=c(0,0,0,0),omi=c(1,1,1,1))
specs<-baseModel$param$species$species
cl <- c("blue","blue","deepskyblue","deepskyblue","red","red","yellow","yellow","green",
        "green","green4","green4")
for (ispec in s_order) {
  bla = which(baseModel$N[laststep,ispec,] < 1) # set limit
  pick = bla[1]
  plot(w[1:pick],baseModel$N[laststep,ispec,1:pick]/baseModel$N[laststep,ispec,1:pick],log="xy",typ="l",lty=1,lwd=2,xaxt="n",yaxt="n",ylim=c(5e-1,5))
  for (i in seq(1,12,2)) points(w[1:pick],baseModel2$N[months[i],ispec,1:pick]/baseModel$N[laststep,ispec,1:pick],typ="l",col=cl[i])
  points(w[1:pick],baseModel2$N[months[5],ispec,1:pick]/baseModel$N[laststep,ispec,1:pick],typ="l",col=cl[5]) # re-do bloom lines
  lines(w[1:pick],baseModel$N[laststep,ispec,1:pick]/baseModel$N[laststep,ispec,1:pick],typ="l",lty=1,lwd=2) # re-do horizontal line
  mtext(specs[ispec],3,line=-1,cex=0.8)
  # add axes to outer plots    
  if (is.element(ispec,c(1,5,7,10)))	axis(side=2, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  if (is.element(ispec,c(10,11,12)))	axis(side=1, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  if (ispec == 1) legend("topright", legend = c("Jan","Mar","May","Jul","Sep","Nov"), col=cl[seq(1,12,2)], lty=1, ncol=2, bty="n") # optional legend
}
# add axes labels on outside
mtext("Body mass (g)",1,outer=T,line=3)
mtext(expression(paste("Seasonal density / non-seasonal density ",sep="")),2,outer=T,line=3)
if (to_save == 1) dev.off()

# TRACK THE DENSITY OF A COHORT THROUGH TIME 

# LOOK AT GROWTH AND SURVIVAL
T = 15/dt; # number of time steps you want to follow cohort for
t_start = laststep - T; # setting the start time for tracking the cohort
z_st = array(0, c(12, T+1)); # row vector for following cohort weight
dc_st = array(0, c(12, T+1)); # vector for cohort survival
z = array(0, c(12, T+1)); # row vector for following cohort weight
dc = array(0, c(12, T+1)); # vector for cohort survival

# NEWBORNS OVER LIFETIME
z_st[,1] = w[1]; # log weight initially (newborn)
dc_st[,1] = baseModel$N[t_start,,1]; # initial population in spectrum
z[,1] = w[1]; # same for seasonal system
dc[,1] = baseModel2$N[t_start,,1];

for (q in seq(1,12)){ 
  for (t in seq(1,T)){ # within time period you're interested in
    zy = max(which(z[q,t] - w >= 0)) # weight bin of cohort from last time step (this will probably need updating for your code)
    z[q,t+1] = z[q,t]+dt*growth[t_start+t-1,q,zy] # using growth rate in that bin to update to z(t-t_start+1)
    dc[q,t+1] = dc[q,t]*exp(-dt*(pred[t_start+t-1,q,zy]+natmort[t_start+t-1,30+zy])) # updating amount surviving using death rate
    zy = max(which(z_st[q,t] - w >= 0)) # weight bin of cohort from last time step (this will probably need updating for your code)
    z_st[q,t+1] = z_st[q,t]+dt*growth_st[t_start+t-1,q,zy] # using growth rate in that bin to update to z(t-t_start+1)
    dc_st[q,t+1] = dc_st[q,t]*exp(-dt*(pred_st[t_start+t-1,q,zy]+natmort_st[t_start+t-1,30+zy])) # updating amount surviving using death rate
  }
}

specs<-baseModel$param$species$species

# GROWTH

if (to_save == 1){
  setEPS()
  postscript("growthcurves.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
par(mfrow=c(4,3),mai=c(0,0,0,0),omi=c(1,1,1,1))
for (ispec in s_order) {
  plot(dt*seq(0,T),z_st[ispec,],type="l",lty=2,lwd=2,xaxt="n",yaxt="n",ylim=c(range(z[ispec,],z_st[ispec,])),col="grey")
  lines(dt*seq(0,T),z[ispec,],type="l",lty=1,lwd=2)
  mtext(specs[ispec],3,line=-1,cex=0.8)
  
  # add axes to outer plots    
  if (is.element(ispec,c(1,5,7,10))) {
    axis(side=2, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  }	else {
    axis(side=4, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)   
  }
  if (is.element(ispec,c(10,11,12)))	axis(side=1, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
}
# add axes labels on outside
mtext("Time (years)",1,outer=T,line=3)
mtext(expression(paste("Mass (g)",sep="")),2,outer=T,line=3)
if (to_save == 1) dev.off()

# SURVIVAL - not in paper
if (to_save == 1){
  setEPS()
  postscript("survival.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
par(mfrow=c(4,3))
for (ispec in 1:12) {
  plot(dt*seq(0,T),dc_st[ispec,],type="l",lty=2,lwd=2,xaxt="n",yaxt="n",ylim=c(1e-4,max(dc[ispec,],dc_st[ispec,])))
  lines(dt*seq(0,T),dc[ispec,],type="l",lwd=2,col="grey")
  mtext(specs[ispec],3,line=-1,cex=0.8)
  # add axes to outer plots    
  if (is.element(ispec,c(1,5,9)))	axis(side=2, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
  if (is.element(ispec,c(length(param$species$species)-3,length(param$species$species)-2,length(param$species$species)-1,length(param$species$species))))	axis(side=1, at = NULL, labels = TRUE, tick = TRUE,cex.axis=0.8,hadj=TRUE)
}
# add axes labels on outside
mtext("Time (years)",1,outer=T,line=3)
mtext(expression(paste("Biomass remaining",sep="")),2,outer=T,line=3)
if (to_save == 1) dev.off()

# Figure 5 - spawning within a year

if (to_save == 1){
  setEPS()
  postscript("mat_growth.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
par(mfrow=c(1,1))
pick = 7 # choose species (7 = sole)
t_pick = 10 # year to start 
t_v=seq(t_pick/dt,t_pick/dt+1/dt) # 13 time steps
plot(seq(0,1,dt), z[pick,t_v],type="l",lty=1,lwd=4,ylab="Weight (g)",xlab="Time (years)",
     ylim=c(min(z_st[pick,t_v],z[pick,t_v]),max(z_st[pick,t_v],z[pick,t_v])),
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
lines(seq(0,1,dt),z_st[pick,t_v],type="l",lty=2,lwd=2,col="grey")
lines(seq(0,1,1/1000),min(z_st[pick,t_v],z[pick,t_v])+(max(z_st[pick,t_v],z[pick,t_v])-min(z_st[pick,t_v],z[pick,t_v+12]))*
        exp(r_peak_s[pick]*cos(2*pi*(seq(0,1,1/1000)-r_peakt[pick])))/exp(r_peak_s[pick]),type="l",lty=3,lwd=2)
if (to_save == 1) dev.off()


# HOW WILL WE SHOW THE EFFECTS OF FISHING ON THE SIZE SPECTRUM?

ns<-summaryCommunityTime(fishModel,minw=1,maxw=10^4)
s<-summaryCommunityTime(fishModel2,minw=1,maxw=10^4)
ufns<-summaryCommunityTime(unfishModel,minw=1,maxw=10^4)
ufs<-summaryCommunityTime(unfishModel2,minw=1,maxw=10^4)
dfns <- summaryCommunityTime(baseModel,minw=1,maxw=10^4)
dfs <- summaryCommunityTime(baseModel2,minw=1,maxw=10^4)

## FIGURE 6

# (a) slope
select = (length(dfns$yield)-length(s$yield)+1):length(dfns$yield)
if (to_save == 1){
  setEPS()
  postscript("slope.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
par(mfrow=c(1,1), mai=c(0.5,1,0.5,1), oma=c(3,1,0,0))
plot(s$time*dt,ns$slope+1,typ="l",lwd=2,ylim=c(-2,-1.4),xlim=c(0,10),ylab="Numerical Size Spectrum Slope",xlab="Time (years)",bty="n",
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
points(s$time*dt,s$slope+1,typ="l",lwd=2,col="grey")
points(s$time*dt,ufns$slope+1,typ="l",lwd=2,lty=2)
points(s$time*dt,ufs$slope+1,typ="l",col="grey",lwd=2,lty=2)
points(s$time*dt,dfns$slope[select]+1,typ="l",lwd=2,lty=4)
points(s$time*dt,dfs$slope[select]+1,typ="l",lwd=2,lty=4,col="grey")
mtext("Time (years)",1,2.5)
if (to_save == 1) dev.off()

# (b) yields
if (to_save == 1){
  setEPS()
  postscript("fishing.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
plot(ns$time*dt,ns$yield/1e12,typ="l",lwd=2,xlim=c(0,10),ylab="Yield (Mt per year)",xlab="Time (years)",ylim=c(1.5,7.5),bty="n",
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2)
points(s$time*dt,s$yield/1e12,typ="l",lwd=2,col="grey")
points(s$time*dt,ufns$yield/1e12,typ="l",lwd=2,lty=2)
points(s$time*dt,ufs$yield/1e12,typ="l",lwd=2,col="grey",lty=2)
points(s$time*dt,dfns$yield[select]/1e12,typ="l",lwd=2,lty=4)
points(s$time*dt,dfs$yield[select]/1e12,typ="l",lwd=2,lty=4,col="grey")
legend("topright", legend = c("Double, non-seasonal","Double, seasonal","Half, non-seasonal","Half, seasonal","Normal, non-seasonal",
                              "Normal, seasonal"), col=c("black","grey","black","grey","black","grey"), lty=c(1,1,2,2,4,4), bty="n") # optional legend
if (to_save == 1) dev.off()

## Figure 7 - total desnity in different portions of baseline
if (to_save == 1){
  setEPS()
  postscript("abundance.eps", horizontal = FALSE, onefile = FALSE, paper = "special")
}
par(mfrow=c(2,1), mai=c(0.5,1,0.5,1), oma=c(3,1,0,0))
plot(ns$time*dt,apply(cspectrum[select,1:60], 1, sum),typ="l",lwd=2,ylab="Total density",bty="n",
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2, xlim=c(0, 10), main="< 40g")
points(ns$time*dt,apply(cspectrum2[select,1:60], 1, sum),typ="l",lwd=2,lty=2,col="grey")
plot(ns$time*dt,apply(cspectrum[select,61:100], 1, sum),typ="l",lwd=2,ylab="Total density",bty="n",
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=1.2, xlim=c(0, 10), main=">= 40g")
points(ns$time*dt,apply(cspectrum2[select,61:100], 1, sum),typ="l",lwd=2,lty=2,col="grey")
mtext("Time (years)",1,2.5)
legend("topright", legend = c("Non-seasonal","Seasonal"), col=c("black", "grey"), lty=c(1,2), ncol=2, bty="n") # optional legend
if (to_save == 1) dev.off()

if (to_save == 0) dev.off() # close pdf

toc = proc.time() - tic # time taken to run it


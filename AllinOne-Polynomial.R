#Uncertainty in Renewability
rm(list=ls())
require(ggplot2)
require(DiceKriging)
require(reshape2)
directory=setwd(getwd())
data <- read.csv(paste(directory,"/data.csv",sep=""))
source('/Users/eansar2/Desktop/Geostatistics/Chapter4/Code/filled.contour3.R')
#Kriging instead of Polynomial 
data$DynamicDateTime001<-as.numeric(substr(data$DynamicDateTime001,1,4))-2012

Ra.fit <- with(data,
               lm(DynamicDateTime001 ~ ReservoirDip + Length + 
                    Height + Permeability + porosity   + I(Height^2) +
                    I(Length^2) + I(porosity^2) + I(Permeability^2) + I(ReservoirDip^2)
                  + Height:Length + Height:porosity  + Height:Permeability 
                  + Height:ReservoirDip + Length:porosity + Length:Permeability
                  + Length:ReservoirDip + porosity:Permeability + porosity:ReservoirDip
                  + Permeability:ReservoirDip))
Factors = intersect(labels(terms(Ra.fit)),colnames(data))
d <- 5; 
design.fact <- with(data,data.frame(x1=ReservoirDip, x2=Length, x3=Height,
                                    x4=Permeability,x5=porosity))
y <- with(data,DynamicDateTime001)
# kriging model 1 : matern5_2 covariance structure, no trend, no nugget effect

m1 <- km(design=design.fact, response=y)

# kriging model 2 : matern5_2 covariance structure,
# linear trend + interactions, no nugget effect
m2 <- km(design=design.fact, response=y)
# graphics
n.grid <- 10
x1.grid <-  with(data,seq(min(ReservoirDip),max(ReservoirDip),length=n.grid))
x2.grid <-  with(data,seq(min(Length),max(Length),length=n.grid))
x3.grid <-  with(data,seq(min(Height),max(Height),length=n.grid))
x4.grid <-  with(data,seq(min(Permeability),max(Permeability),length=n.grid))
x5.grid <-  with(data,seq(min(porosity),max(porosity),length=n.grid))
design.grid <- expand.grid(x1=x1.grid, x2=x2.grid, x3=x3.grid, x4=x4.grid,x5=x5.grid)
response.grid <- with(data,DynamicDateTime001)
predicted.values.model1 <- predict(m1, design.grid, "UK")$mean #UK is OK with no trend.
predicted.values.model2 <- predict(m2, design.grid, "SK")$mean


#------------
#Height-Perm (not needed)
#------------
colnames(design.grid)=Factors
response.grid = predict(Ra.fit,design.grid)
mean.response = matrix(rep(0,length(x3.grid)*length(x4.grid)),
                       c(length(x3.grid),length(x4.grid)))

for (i in 1:length(x3.grid)) {
  for (j in 1:length(x4.grid)){
    
    mean.response[i,j] = mean(response.grid[which(design.grid[,3]==x3.grid[i] & 
                                                    design.grid[,4]==x4.grid[j])])
  }
  
}

pdf("Height_Perm.pdf")
par(mfrow=c(1,1),pty="s")
adjust=0.001
xlim=c(range(x3.grid)[1]-0.01*range(x3.grid)[1],range(x3.grid)[2]+0.01*range(x3.grid)[2])
ylim=c(range(x4.grid)[1]-.01*range(x4.grid)[1],range(x4.grid)[2]+.12*range(x4.grid)[2])
filled.contour3(x3.grid, x4.grid, mean.response,
                color.palette=topo.colors,xlim=xlim,ylim=ylim,xlab="Height"
                ,ylab="Permeability",yaxs='i')
contour(x3.grid, x4.grid, mean.response, 5,labcex=1,add=TRUE,yaxs='i')
points(design.fact[,3], design.fact[,4], pch=19, cex=1.5, col="blue",yaxs='i')

dev.off()

#------------
#Dip-Length
#------------
colnames(design.grid)=Factors
response.grid = predict(Ra.fit,design.grid)
mean.response = matrix(rep(0,length(x1.grid)*length(x2.grid)),
                       c(length(x1.grid),length(x2.grid)))

for (i in 1:length(x1.grid)) {
  for (j in 1:length(x2.grid)){
    
    mean.response[i,j] = mean(response.grid[which(design.grid[,1]==x1.grid[i] & 
                                                    design.grid[,2]==x2.grid[j])])
  }
  
}

pdf("Dip_Length.pdf")
par(mfrow=c(1,1),pty="s")
adjust=0.001
xlim=c(range(x1.grid)[1]-0.01*range(x1.grid)[1],range(x1.grid)[2]+0.01*range(x1.grid)[2])
ylim=c(range(x2.grid)[1]-.01*range(x2.grid)[1],range(x2.grid)[2]+.12*range(x2.grid)[2])
filled.contour3(x1.grid, x2.grid, mean.response,
                color.palette=topo.colors,xlim=xlim,ylim=ylim,xlab="Dip"
                ,ylab="Length",yaxs='i')
contour(x1.grid, x2.grid, mean.response, 5,labcex=1,add=TRUE,yaxs='i')
points(design.fact[,1], design.fact[,2], pch=19, cex=1.5, col="blue",yaxs='i')

dev.off()

#------------
#Dip-Thickness
#------------
colnames(design.grid)=Factors
response.grid = predict(Ra.fit,design.grid)
mean.response = matrix(rep(0,length(x1.grid)*length(x3.grid)),
                       c(length(x1.grid),length(x3.grid)))

for (i in 1:length(x1.grid)) {
  for (j in 1:length(x3.grid)){
    
    mean.response[i,j] = mean(response.grid[which(design.grid[,1]==x1.grid[i] & 
                                                    design.grid[,3]==x3.grid[j])])
  }
  
}

pdf("Dip_Thickness.pdf")
par(mfrow=c(1,1),pty="s")
adjust=0.001
xlim=c(range(x1.grid)[1]-0.01*range(x1.grid)[1],range(x1.grid)[2]+0.01*range(x1.grid)[2])
ylim=c(range(x3.grid)[1]-.01*range(x3.grid)[1],range(x3.grid)[2]+.12*range(x3.grid)[2])
filled.contour3(x1.grid, x3.grid, mean.response,
                color.palette=topo.colors,xlim=xlim,ylim=ylim,xlab="Dip"
                ,ylab="Thickness",yaxs='i')
contour(x1.grid, x3.grid, mean.response, 5,labcex=1,add=TRUE,yaxs='i')
points(design.fact[,1], design.fact[,3], pch=19, cex=1.5, col="blue",yaxs='i')

dev.off()

#------------
#Poro-Thickness
#------------
colnames(design.grid)=Factors
response.grid = predict(Ra.fit,design.grid)
mean.response = matrix(rep(0,length(x5.grid)*length(x3.grid)),
                       c(length(x5.grid),length(x3.grid)))

for (i in 1:length(x5.grid)) {
  for (j in 1:length(x3.grid)){
    
    mean.response[i,j] = mean(response.grid[which(design.grid[,5]==x5.grid[i] & 
                                                    design.grid[,3]==x3.grid[j])])
  }
  
}

pdf("Porosity_Thickness.pdf")
par(mfrow=c(1,1),pty="s")
adjust=0.001
xlim=c(range(x5.grid)[1]-0.01*range(x5.grid)[1],range(x5.grid)[2]+0.01*range(x5.grid)[2])
ylim=c(range(x3.grid)[1]-.01*range(x3.grid)[1],range(x3.grid)[2]+.12*range(x3.grid)[2])
filled.contour3(x5.grid, x3.grid, mean.response,
                color.palette=topo.colors,xlim=xlim,ylim=ylim,xlab="Porosity"
                ,ylab="Thickness",yaxs='i')
contour(x5.grid, x3.grid, mean.response, 5,labcex=1,add=TRUE,yaxs='i')
points(design.fact[,5], design.fact[,3], pch=19, cex=1.5, col="blue",yaxs='i')

dev.off()

#------------
#Porosity-Permeability 
#------------
colnames(design.grid)=Factors
response.grid = predict(Ra.fit,design.grid)
mean.response = matrix(rep(0,length(x5.grid)*length(x4.grid)),
                       c(length(x5.grid),length(x4.grid)))

for (i in 1:length(x5.grid)) {
  for (j in 1:length(x4.grid)){
    
    mean.response[i,j] = mean(response.grid[which(design.grid[,5]==x5.grid[i] & 
                                                    design.grid[,4]==x4.grid[j])])
  }
  
}

pdf("Porosity_Permeability.pdf")
par(mfrow=c(1,1),pty="s")
adjust=0.001
xlim=c(range(x5.grid)[1]-0.01*range(x5.grid)[1],range(x5.grid)[2]+0.01*range(x5.grid)[2])
ylim=c(range(x4.grid)[1]-.01*range(x4.grid)[1],range(x4.grid)[2]+.12*range(x4.grid)[2])
filled.contour3(x5.grid, x4.grid, mean.response,
                color.palette=topo.colors,xlim=xlim,ylim=ylim,xlab="Porosity",
                ylab="Permeability",yaxs='i')
contour(x5.grid, x4.grid, mean.response, 5,labcex=1,add=TRUE,yaxs='i')
points(design.fact[,5], design.fact[,4], pch=19, cex=1.5, col="blue",yaxs='i')

dev.off()

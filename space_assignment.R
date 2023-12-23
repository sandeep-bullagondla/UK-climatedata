setwd("C:/Users/Sandeep kumar/OneDrive/Desktop/R")
temp<-read.csv("MaxTemp.csv")

meta<-read.csv("metadata.csv")

mutate(temp, date=as.Date(Date,origin = "2020-01-01",format="%Y-%m-%d"))
time<-as.POSIXct(temp$Date,format="%Y-%m-%d")

library(knitr)
kable(summary(temp))
library(lubridate)
temp$Date=ymd(temp$Date)
library(ggplot2)
library(tidyverse)
library(maps)
UK <- map_data(map = "world", region = "UK") # changed map to "world"


library(stringr)
library(tidyverse)

meltData <- melt(temp,id.vars = "Date")

library(ggplot2)
ggplot(meltData, aes(factor(variable), value))+
 geom_boxplot() + facet_wrap(~variable, scale="free")
  
a<-ggplot() + 
    geom_polygon(data = UK, aes(x = long, y = lat, group = group), 
                           fill = 'gray90', color = 'black') + 
    geom_point(data=sept12,aes(x=Longitude,y=Latitude,colour=value,size=Elevation))+
    coord_fixed(ratio = 1.3, xlim = c(-10,3),ylim = c(50, 59))+
  labs(x="Longitude",y="Latitude",title="(a) Twenty location points on UK map")
a+b


library(reshape)
sept12<-melt(temp[256,],id.vars = "Date",variable_name ="Location",value.name="Maximum Temperature")
sept12<-merge(sept12,meta,on="Location")
predicts<-sept12[c(2,6,14),]
sept12<-sept12[-c(2,6,14),]


library(geoR)
set.seed(356)
#for sept 12
gdata1 <- as.geodata(sept12,coords.col=4:5,data.col=3,borders=as.matrix(uk_cords)) # selecting SST
dup <- dup.coords(gdata1)
gdata1 <- jitterDupCoords(gdata1,max=0.1,min=0.05) # jittering duplicates (geoR can't handle different data values in the same position)



gdata3 <- as.geodata(predicts,coords.col=4:5,data.col=3)
dup <- dup.coords(gdata3)
gdata3 <- jitterDupCoords(gdata3,max=0.1,min=0.05) # jittering duplicates (geoR can't handle different data values in the same position)

plot(gdata1)

parana
plot(variog4(gdata1))
plot(variog4(gdata1,trend = "1st"))
library(gstat)
var<-variogram(value~1,
               locations=~Longitude+Latitude,
               data=sept12,alpha=c(0,45,90,135))
plot(var,ylim=c(0,6))

sample_vario<-variog(gdata1,trend="1st",option='bin')
plot(sample_vario)

sample_vario$n
par(mar=c(4,4,2,2), mfrow=c(1,3))
#plotting variogram using hawkins and cressie estimator
plot(variog(gdata1, trend="1st",option='cloud',
            estimator.type='modulus'),pch = 19)


sample.vario <- variog(gdata1,option='bin', trend="1st",ini=inits,max.dist=4,
                       estimator.type='modulus', bin.cloud = TRUE)

plot(sample.vario, pch = 19)
plot(sample.vario, bin.cloud = TRUE)

var.default<-variofit(sample.vario)

var.mat1<-variofit(sample.vario,kappa = 1,fix.kappa = T)

var.mat1.5<-variofit(sample.vario,ini=inits,kappa = 1.5,fix.kappa = T)


xv.ml <- xvalid(parana, model = var.mat1)
par(mfcol = c(5,2), mar = c(4,4,1,1))
plot(xv.ml)


par(mfrow=c(1,1))
plot(sample.vario,max.dist = 4)
lines(ml_matern)
lines(var.mat1)
lines(var.mat1.5)
lines(reml_exponential,lty=2)

pred.grid <- expand.grid(seq(min(gdata1$coords[,1]),max(gdata1$coords[,1]),by=0.1),
                         seq(min(gdata1$coords[,2]),max(gdata1$coords[,2]),by=0.1))
#predictions of data by using vari.mat model
preds <- krige.conv(gdata1,
                    loc=pred.grid, krige=krige.control(obj.model=var.mat1.5))
min(preds$predict)
image(preds, col = viridis::viridis(100), zlim = c(0,max(c(preds$predict))),
      coords.data =gdata3[1]$coords , main = 'Mean', xlab = 'x', ylab = 'y',
      x.leg = c(1,4), y.leg = c(51, 52))
par(mfrow=c(1,2))
image(preds, values = preds$krige.var, col = heat.colors(100)[100:1],
      zlim = c(0,max(c(preds$krige.var))), coords.data = gdata3[1]$coords,
      main = 'Variance', xlab = 'x', ylab = 'y', x.leg = c(1, 4), y.leg = c(51, 51.5))

par(mfrow=c(1,1))
points(gdata1, pch = 24,col = c("violet","orange","green","red"))

inits<-c(4,0.5)
#matern
ml_matern<-likfit(gdata1,cov.model = "matern",kappa = 1.5,trend="1st",ini=inits,fix.nugget = F,
                  fix.kappa = T)

reml_matern<-likfit(gdata1,cov.model = "matern",kappa = 1.5,trend="1st",ini=inits,
                    fix.nugget = F,
                    lik.method = 'REML')

summary(ml_matern)
summary(reml_matern)
plot(variog(gdata1))
xv.matern<-xvalid(gdata1,model = reml_matern)
par(mfcol=c(5,2),mar=c(4,4,1,1))
plot(xv.matern)

#grid of coordinates in sampled data on 0.5 degree
pred.grid <- expand.grid(seq(min(gdata1$coords[,1]),max(gdata1$coords[,1]),by=0.1),
                         seq(min(gdata1$coords[,2]),max(gdata1$coords[,2]),by=0.1))
#predictions for matern model 1st trend
kc <- krige.conv(gdata1, loc = pred.grid, krige = krige.control(obj.model = reml_matern))

## krige.conv: model with mean given by a 1st order polynomial on the coordinates
## krige.conv: Kriging performed using global neighbourhood
par(mfrow=c(1,1))
#Image for mean of predictions with reml_matern1 model
image(kc, col = viridis::viridis(100), zlim = c(0,max(c(kc$predict))),
      coords.data = gdata1[1]$coords, main = 'Mean with actual locations',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))


image(kc, col = viridis::viridis(100), zlim = c(0,max(c(kc$predict))),
      coords.data = gdata3[1]$coords, main = 'Mean with prediction locations',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))

#Image for variance of predictions with reml_matern1 model
image(kc, values = kc$krige.var, col =  heat.colors(100)[100:1],
      zlim = c(0,max(c(kc$krige.var))), coords.data = gdata1[1]$coords,
      main = 'Variance with actual locations',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))

image(kc, values = kc$krige.var, col =  heat.colors(100)[100:1],
      zlim = c(0,max(c(kc$krige.var))), coords.data = gdata3[1]$coords,
      main = 'Variance with predicted locations',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))

predicts



#exponential
ml_exponential<-likfit(gdata1,cov.model = "exponential",ini=inits,trend="1st",fix.nugget = F)

reml_exponential<-likfit(gdata1,cov.model = "exponential",trend="1st",ini=inits,fix.nugget = F,
                         lik.method = 'REML')

summary(ml_exponential)

xv.exponential<-xvalid(gdata1,model = reml_exponential)
par(mfcol=c(5,2),mar=c(4,4,1,1))
plot(xv.exponential)


kc_exp <- krige.conv(gdata1, loc = pred.grid, krige = krige.control(obj.model = reml_exponential))
## krige.conv: model with mean given by a 1st order polynomial on the coordinates
## krige.conv: Kriging performed using global neighbourhood
par(mfrow=c(1,1))
#Image for mean of predictions with reml_matern1 model
image(kc_exp, col = viridis::viridis(100), zlim = c(0,max(c(kc_exp$predict))),
      coords.data = gdata3[1]$coords, main = 'Mean-exponential',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))

min(kc_exp$predict)
gdata1
#Image for variance of predictions with reml_matern1 model
image(kc_exp, values = kc_exp$krige.var, col = heat.colors(100)[100:1],
      zlim = c(0,max(c(kc_exp$krige.var))), coords.data = gdata3[1]$coords,
      main = 'Variance-Matern REML 1st trend',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))

predicts
max(temp[,-1])
min(temp[,-1])
ggplot(temp)+aes(x=Date,y=Temp)



#for 26th april
aprl26<-melt(temp[117,],id.vars = "Date",variable_name ="Location",value.name="Maximum Temperature")
aprl26<-merge(aprl26,meta,on="Location")

set.seed(356)
gdata4 <- as.geodata(aprl26,coords.col=4:5,data.col=3) # selecting SST
dup <- dup.coords(gdata4)
gdata4 <- jitterDupCoords(gdata4,max=0.1,min=0.05) # jittering duplicates (geoR can't handle different data values in the same position)


reml_matern_aprl25<-likfit(gdata4,cov.model = 'matern',kappa = 1.5,ini=inits,
                           trend = "1st",fix.nugget = F)
  
pred.grid1 <- expand.grid(seq(min(gdata4$coords[,1]),max(gdata4$coords[,1]),by=0.1),
                         seq(min(gdata4$coords[,2]),max(gdata4$coords[,2]),by=0.1))
#predictions for matern model 1st trend
kc_a <- krige.conv(gdata4, loc = pred.grid1, krige = krige.control(obj.model = reml_matern_aprl25))

## krige.conv: model with mean given by a 1st order polynomial on the coordinates
## krige.conv: Kriging performed using global neighbourhood
par(mfrow=c(1,1))
#Image for mean of predictions with reml_matern1 model
image(kc_a, col = viridis::viridis(100), zlim = c(0,max(c(kc_a$predict))),
      coords.data = gdata4[1]$coords, main = 'Mean for April 26',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))


#Image for variance of predictions with reml_matern1 model
image(kc_a, values = kc_a$krige.var, col =  heat.colors(100)[100:1],
      zlim = c(0,max(c(kc_a$krige.var))), coords.data = gdata4[1]$coords,
      main = 'Variance for April 26',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5)) 

#for dec 5th
dec5<-melt(temp[366,],id.vars = "Date",variable_name ="Location",value.name="Maximum Temperature")
dec5<-merge(dec5,meta,on="Location")


gdata5 <- as.geodata(dec5,coords.col=4:5,data.col=3) # selecting SST
dup <- dup.coords(gdata5)
gdata5 <- jitterDupCoords(gdata5,max=0.1,min=0.05) # jittering duplicates (geoR can't handle different data values in the same position)


reml_matern_dec5<-likfit(gdata5,cov.model = 'matern',kappa = 1.5,ini=inits,
                           trend = "1st",fix.nugget = F)

pred.grid2 <- expand.grid(seq(min(gdata5$coords[,1]),max(gdata5$coords[,1]),by=0.1),
                          seq(min(gdata5$coords[,2]),max(gdata5$coords[,2]),by=0.1))
#predictions for matern model 1st trend
kc_d <- krige.conv(gdata5, loc = pred.grid2, krige = krige.control(obj.model = reml_matern_dec5))

## krige.conv: model with mean given by a 1st order polynomial on the coordinates
## krige.conv: Kriging performed using global neighbourhood
par(mfrow=c(1,1))
#Image for mean of predictions with reml_matern1 model
image(kc_d, col = viridis::viridis(100), zlim = c(-2,max(c(kc_d$predict))),
      coords.data = gdata5[1]$coords, main = 'Mean for Dec 31',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))


#Image for variance of predictions with reml_matern1 model
image(kc_d, values = kc_d$krige.var, col =  heat.colors(100)[100:1],
      zlim = c(0,max(c(kc_d$krige.var))), coords.data = gdata5[1]$coords,
      main = 'Variance for Dec 31',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))
  
temp[351,]

#for july 5th
jul5<-melt(temp[187,],id.vars = "Date",variable_name ="Location",value.name="Maximum Temperature")
jul5<-merge(jul5,meta,on="Location")

set.seed(356)
gdata6 <- as.geodata(jul5,coords.col=4:5,data.col=3) # selecting SST
dup <- dup.coords(gdata6)
gdata6 <- jitterDupCoords(gdata6,max=0.1,min=0.05) # jittering duplicates (geoR can't handle different data values in the same position)


reml_matern_jul5<-likfit(gdata4,cov.model = 'matern',kappa = 1.5,ini=inits,
                           trend = "1st",fix.nugget = F)

pred.grid1 <- expand.grid(seq(min(gdata4$coords[,1]),max(gdata4$coords[,1]),by=0.1),
                          seq(min(gdata4$coords[,2]),max(gdata4$coords[,2]),by=0.1))
#predictions for matern model 1st trend
kc_ju <- krige.conv(gdata6, loc = pred.grid1, krige = krige.control(obj.model = reml_matern_jul5))

## krige.conv: model with mean given by a 1st order polynomial on the coordinates
## krige.conv: Kriging performed using global neighbourhood
par(mfrow=c(1,1))
#Image for mean of predictions with reml_matern1 model
image(kc_ju, col = viridis::viridis(100), zlim = c(0,max(c(kc_ju$predict))),
      coords.data = gdata6[1]$coords, main = 'Mean for July 5',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))


#Image for variance of predictions with reml_matern1 model
image(kc_ju, values = kc_ju$krige.var, col =  heat.colors(100)[100:1],
      zlim = c(0,max(c(kc_ju$krige.var))), coords.data = gdata6[1]$coords,
      main = 'Variance for July 5',
      xlab = 'X Coord', ylab = 'Y Coords',
      x.leg = c(1,4), y.leg = c(51, 51.5))

time_series<-ts(temp,start = 1,frequency = 1)
time(time_series)

Date<-temp$Date
library(ggplot2)
b<-ggplot()+geom_line(aes(x=Date,y=yev))+
  geom_line(aes(x=Date,y=mor,colour="Morecambe"))+
  geom_line(aes(x=Date,y=sheff,colour="Sheffield"))+
  geom_line(aes(x=Date,y=lon,colour="London"))+
  geom_line(aes(x=Date,y=hig,colour="High Wycombe"))+
  geom_line(aes(x=Date,y=ler,colour="Lerwick"))+
  geom_line(aes(x=Date,y=cam,colour="Camborne"))+
  geom_line(aes(x=Date,y=loss,colour="Lossiemouth"))+
  geom_line(aes(x=Date,y=kin,colour="Kinross"))+
  labs(y="Maximum Temperatures",x="Date",
       title = "(b) Maximum temperatures in a day of some locations in 2020")
  

library(patchwork)
a+b

diff_yev<-(diff(yev))

c<-ggplot()+geom_line(aes(x=Date,y=yev),colour="blue")+
  labs(y="maximum temperatue",title="a) Yeovilton time series till October")

d<-ggplot()+geom_line(aes(x=Date[1:304],y=diff_yev),colour="blue")+
  labs(x="Date",y="Difference",title="b) 1st diffrencing of Yeovilton time series")
  
c+d
yev<-time_series[,"Yeovilton"]
mor<-time_series[,"Morecambe"]
sheff<-time_series[,"Sheffield"]
lon<-time_series[,"London"]
hig<-time_series[,"High_Wycombe"]
ler<-time_series[,"Lerwick"]
cam<-time_series[,"Camborne"]
loss<-time_series[,"Lossiemouth"]
kin<-time_series[,"Kinross"]
#fitting auto arma with p and q values
library(forecast)
ggplot()+geom_line(aes(x=time(time_series),y=yev))+
  geom_line(aes(x=time(time_series),y=predicts,colour="Predicts"))
  


model1<-auto.arima(yev, max.p = 4, max.q = 4, max.d = 4, seasonal = FALSE)

tsdiag(model1)

pre<-forecast(model3,h=7)

model2<-auto.arima(yev, trace = T)

summary(model2)
predicts<-fitted(model1)

model3<-arima(yev,order=c(1,1,2))
model3
predict(model1,start=305,end=312)

time_series[,"Yeovilton"][305:312]

as.data.frame(pre)

ggplot(meltData, aes(factor(variable), value))+
  geom_boxplot() + facet_wrap(~variable, scale="free")+xlab("Location")+
  ggtitle("Boxplot for outliers")

summary(ml_matern)
summary(reml_matern)
summary(ml_exponential)
summary(reml_exponential)

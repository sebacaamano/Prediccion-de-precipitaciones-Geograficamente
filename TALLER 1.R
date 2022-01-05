setwd("C:/Users/sebas/Universidad/Seminario I/Proyecto/lluvia de 2010-2020")
library(dplyr)
library(geoR)
library(readr)
library(spdep)
library(akima)
library(automap)
# promedio=read.csv("promedio.txt", sep=" ")
# minimo=read.csv("minimo.txt", sep=" ")
# maximo = read.csv("maximo.txt", sep=" ")
# cordena = read.csv("cor9.txt", sep=" ")


# p1=data.frame(promedio)
# p2=data.frame(minimo)
# p3=data.frame(maximo)
# p4=data.frame(cordena)


# regiones=merge(p4,p1)
# regiones=merge(regiones,p3)
# regiones=merge(regiones,p4)
# write.table(regiones,"precipitaciones.txt", sep=" ", 
#            col.names=TRUE, row.names=FALSE, quote=TRUE, na="NA", dec=".")


lluvia=read.table("precipitaciones.txt", header=TRUE)
X=lluvia$Longitud
Y=lluvia$Latitud
lluviaP=lluvia$Promedio


summary(lluviaP)

# Graficos
par(mfrow = c(1,1))
boxplot(lluviaP,
        main = "Promedio de Precipitaciones, en Verano",
        xlab = "Milimetros",
        ylab = "Precipitacion",
        col = "orange",
        border = "brown",
        horizontal = TRUE,
        notch = TRUE
)
par(mfrow = c(1,2))
hist(lluviaP,
     main = "Histograma de Precipitaciones, en Verano",
     xlab = "Milimetros",
     ylab = "Frecuencia",
     col = "orange",
     xlim = c(0,250),
     ylim = c(0,30))
qqnorm(lluviaP)
qqline(lluviaP, col = 2)

# Prueba de normalidad

shapiro.test(lluviaP)

logdataO<- log(lluviaP)
hist(logdataO)
shapiro.test(logdataO)

par(mfrow = c(2,2))
hist(X)
qqnorm(X)
qqline(X, col = 2)

hist(Y)
qqnorm(Y)
qqline(Y, col = 2)

# interpolacion para la lluviaP
par(mfrow=c(1,1))

int.lluviaP <- interp.new(Y, X, lluviaP)
image(int.lluviaP, ylab="X (latitud)", xlab="Y (longitud)", 
      main ="Interpolación para lluvia promedio")

par(mfrow=c(1,1))
# Curvas de nivel para la lluviaP
contour(int.lluviaP, ylab="X (latitud)", xlab="Y (longitud)", 
        main = "Curvas de Nivel para lluvia promedio")

# Combinamos las cordenadas con la lluviaP
obj <- cbind(X, Y, lluviaP)

# Transformamos los datos a cordenadas espaciales
lluviaP.geo <- as.geodata(obj, coords.col = 1:2, data.col = 3)
plot(lluviaP.geo)
summary(lluviaP.geo)
maxd = summary(lluviaP.geo)$distances.summary[[2]] #distancia maxima

vario.b = variog(lluviaP.geo, lambda=0, estimator.type = "modulus", max.dist = 4)
plot(vario.b)
#eyefit(vario.b)

par(mfrow = c(1,1))

vario.4 <- variog4(lluviaP.geo, max.dist = maxd/2)
plot(vario.4, omnidirectional = T, legend = F)


vario.wls = variofit(vario.b, cov.model = "exponential", ini = c(0.46,2.18), 
                     nugget = 0.02 , weights = "cressie")
summary(vario.wls)#minimoscuadradosponderados
datgeo=lluviaP.geo
vario.ml = likfit(datgeo, trend = "cte", lambda=0, cov.model = "exponential", 
                  ini = c(0.46,2.18), nugget = 0.02, lik.method = "ML", 
                  fix.psiA = FALSE, fix.psiR = FALSE)#max.verosimilitud
summary(vario.ml)#max.verosimilitud

vario.reml = likfit(datgeo, trend = "cte", lambda=0, cov.model = "exponential", 
                    ini = c(0.46,2.18), nugget = 0.02 , lik.method = "RML", 
                    fix.psiA = FALSE, fix.psiR = FALSE)
summary(vario.reml)#m.v.restricta

plot(vario.b, main= "Variograma empirico",ylab="Semivarianza", xlab="Distancia")
lines(vario.ml, max.dist=, pch= 14, col="blue") #max.verosimilitud
lines(vario.reml, lwd= 2,  max.dist=4, pch= 14, col="red") #m.v.restricta
lines(vario.wls, lty= 2, lwd= 2,  max.dist=4, pch= 14, col="purple") #minimos cuadrados pesados
legend(x = "bottomright" ,legend= c("MaxV", "MaxVRe", "MCP"), col=c("blue","red","purple"),lty =c(1,1,1),lwd=2)

######## Kroos Validación

xv.wls = xvalid(datgeo, model = vario.wls) 
summary(xv.wls)#MCP

xv.reml = xvalid(datgeo, model = vario.reml)
summary(xv.reml)#MaxVRe

xv.ml = xvalid(datgeo, model = vario.ml)
summary(xv.ml)#MaxV


#######KRIGING#####
#leemos la data nuevamente
Data=read.table("precipitaciones.txt", header=TRUE)
Data=data.frame(Data[,1],Data[,2],Data[,3])

Temperatura.geo <- as.geodata(Data, coords.col = 1:2, data.col = 3)
delta = 0.5
min1 <- summary(Temperatura.geo)$coords.summary[1,1]-delta
min2 <- summary(Temperatura.geo)$coords.summary[1,2]-delta
max1 <- summary(Temperatura.geo)$coords.summary[2,1]+delta
max2 <- summary(Temperatura.geo)$coords.summary[2,2]+delta
loci <- expand.grid(seq(min1, max1, by=0.05), seq(min2, max2, by=0.05))



### Geo Data ###
datagrid= as.geodata(Data, coords.col = 1:2, data.col = 3)
coords=datagrid$coords

######## Modelo 1 Kriging Ordinario
KrigeMVR= krige.conv(datagrid, coords= datagrid$coords, data= datagrid$data, locations= loci, 
                     krige=krige.control(type.krige = "OK", trend.d = "cte", trend.l = "cte",
                                         cov.model = "exponential", cov.pars = c(0.2374,1.048),
                                         nugget = 0.0087, kappa = 0))


persp(KrigeMVR, locations= loci, values = KrigeMVR$predict)
points(datagrid, pt.divide = "quintile", xlab= "Coord X",ylab= "Coord Y",main="Gráfico de datos espaciales, Kriging Ordinario")
contour(KrigeMVR, add=T)

par(mfrow=c(1,1))
image(KrigeMVR, filled=TRUE, main="Kriging Ordinario")
points(datgeo$coords,pch=20)
contour(KrigeMVR, add=T)

#muestra lo q predice respecto las coordenadas
Rmvord=cbind(loci$Var1,loci$Var2,KrigeMVR$predict);Rmvord


###### MODELO 2 Kriging Simple

KrigSimple = krige.conv(datagrid, coords = datagrid$coords, data = datagrid$data, locations = loci,
                        krige = krige.control(type.krige = "SK", trend.d = "cte", trend.l = "cte",
                                              cov.model = "exponential", cov.pars = c(0.2374,1.048), 
                                              nugget = 0.0087, kappa = 0, beta=4.3213)) 

persp(KrigSimple, locations = loci, values=KrigSimple$predict)

points(datagrid, pt.divide = "quintile", xlab="Coord X", ylab= "Coord Y",main="Gráfico de datos espaciales, Kriging Simple")
contour(KrigSimple, add=T)

par(mfrow=c(1,1))
image(KrigSimple, filled=TRUE, main="Kriging Simple")
points(datgeo$coords,pch=20)
contour(KrigSimple, add=T)

#muestra lo q predice respecto las coordenadas
Rmvsimp=cbind(loci$Var1,loci$Var2,KrigSimple$predict);Rmvsimp




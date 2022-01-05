setwd("C:/Users/sebas/Universidad/Seminario I/Proyecto/lluvia de 2010-2020")
library(dplyr)
library(geoR)
library(readr)
library(spdep)
library(akima)
library(automap)
#######KRIGING#####
#leemos la data nuevamente
Data=read.table("precipitaciones.txt", header=TRUE)
Data=data.frame(Data[,1],Data[,2],Data[,3])

Temperatura.geo <- as.geodata(Data, coords.col = 1:2, data.col = 3)
delta = 6
min1 <- summary(Temperatura.geo)$coords.summary[1,1]-delta
min2 <- summary(Temperatura.geo)$coords.summary[1,2]-delta
max1 <- summary(Temperatura.geo)$coords.summary[2,1]+delta
max2 <- summary(Temperatura.geo)$coords.summary[2,2]+delta
loci <- expand.grid(seq(min1-delta, max1+delta, by=0.05), seq(min2-delta, max2+delta, by=0.05))



### Geo Data ###
datagrid= as.geodata(Data, coords.col = 1:2, data.col = 3)
coords=datagrid$coords

######## Modelo 1 Kriging Ordinario
KrigeMVR= krige.conv(datagrid, coords= datagrid$coords, data= datagrid$data, locations= loci, 
                     krige=krige.control(type.krige = "OK", trend.d = "cte", trend.l = "cte",
                                         cov.model = "exponential", cov.pars = c(0.2374,1.048),
                                         nugget = 0.0087, kappa = 0))


points(datagrid, pt.divide = "quintile", xlab= "Coord X",ylab= "Coord Y",main="Gráfico de datos espaciales, Kriging Ordinario")
contour(KrigeMVR, add=T)

###### MODELO 2 Kriging Simple

KrigSimple = krige.conv(datagrid, coords = datagrid$coords, data = datagrid$data, locations = loci,
                        krige = krige.control(type.krige = "SK", trend.d = "cte", trend.l = "cte",
                                              cov.model = "exponential", cov.pars = c(0.2374,1.048), 
                                              nugget = 0.0087, kappa = 0, beta=4.3213)) 


points(datagrid, pt.divide = "quintile", xlab="Coord X", ylab= "Coord Y",main="Gráfico de datos espaciales, Kriging Simple")
contour(KrigSimple, add=T)


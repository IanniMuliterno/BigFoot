#library(kableExtra) # MEMO
#require(map)
#require(purrr)
#require(xtable)
library(akima)
require(geoR)  
require(moments)
require(sp)
require(gstat)
library(readxl)
library(tmaptools)
library(rgeos)
library(raster)
library(leaflet)
library(tibble)
library(mvtnorm)
library(spTest)

bf <- read_excel("C:/Users/tuchi/Desktop/BigFoot-master/BigFoot-master/bigfoot USA.xlsx")

#Removendo Alaska
#Alaska tecnicamente faz parte dos EUA, mas fica muito distante do resto do pais
#comentar a linha abaixo para incluir o Alaska
bf <- bf[-2,] 

# Garantindo q os centroides sao coletados para os estados por ordem alfabetica
states <- sort(bf$State) 

#carregando shapefile dos EUA 
usa <- raster::getData(country="USA",level=1) #comentar esta linha depois de rodar pela primeira vez (tempo)
# visualizar o mapa retirado do GADM
#plot(usa)
#library(zoom)
#zoom()
#usa <- usa[usa@data$NAME_1 != "Alaska",]# removendo duplicata de um estado no shapefile do GADM
usa_db <- usa@data


#encontrando os estados sem avistamento do sasquatch
usa_db$NAME_1[which( usa_db$NAME_1 %in% states == FALSE)]
#which( states %in% usa_db$NAME_1 == FALSE)

# incluindo na analise os dois estados sem registro de avistamento
#bf <- rbind(bf,c("District of Columbia" ,0,NA))
#bf <- rbind(bf,c("Hawaii" ,0,NA))


#ordenando banco por ordem alfabetica dos estados
bf <- bf[order(bf$State),]


states <- sort(bf$State)  # atualizando a lista com todos os estados


#coletando centroides 


centros <- data.frame(lati=0,long=0)

for(i in 1:length(states)){
  if(states[i] %in% usa_db$NAME_1){  
  centros[i,]<-coordinates(gCentroid(usa[usa_db$NAME_1==states[i],]))
  }
i} 

bf <- cbind(bf,centros)

centros <- SpatialPoints(centros)
crs.geo <- CRS("+init=EPSG:32633")  
proj4string(centros) <- crs.geo 
centros<-SpatialPointsDataFrame(centros,tibble(attr1 = bf$State))

#visualizacao

plot(centros)

# visualizar dados interpolados

bf.interp=interp(bf$lati,bf$long,log(bf$Sightings))
image(bf.interp,col=terrain.colors(100))
title(main="Interpolated data")
#leaflet() %>% addTiles() %>% addMarkers(data = centros, 
#                                        popup = ~attr1)

#library(tmap)
# Funciona melhor sem considerar o Alaska

#tm_shape(usa) +
#  tm_fill("NAME_1", thres.poly = 0) 

#conversao necessaria tibble -> matrix -> numeric
Isaw <- bf$Sightings
Isaw <- as.matrix(Isaw)
Isaw <- as.numeric(Isaw)
hist(Isaw)
boxplot(Isaw)
classes <- factor(cut(Isaw, breaks=nclass.Sturges(Isaw))) #quebra em classes
output <- as.data.frame(table(classes))
#Add cumFreq and proportions
output <- transform(output, cumFreq = cumsum(Freq), relative = prop.table(Freq))
#xtable(output)

#tabela de estatisticas descritivas
dv <- sd(Isaw)
CV = sd(Isaw)*100/mean(Isaw)
summary(Isaw)
tst<-as.data.frame.vector(as.table(summary(Isaw)))
name <- c(row.names(tst),"dv","CV","skewness","kurtosis")
value <- c(tst$`as.table(summary(Isaw))`,dv,CV,skewness(Isaw),kurtosis(Isaw))
tst<-data.frame(value)
row.names(tst)<- name
tst

bf2 <- bf[,-c(1,3)] # retirando colunas 'obsoletas'
coordinates(bf2) <- cbind(bf2$lati,bf2$long)
bf2$Sightings <- as.numeric(bf2$Sightings)

d <- as.geodata(bf2)
#normalidade da vari?vel resposta
mean(d[["data"]])
sd(d[["data"]])
ks.test(d[["data"]], "pnorm", 5, 8)
mean(log(d[["data"]]))
sd(log(d[["data"]]))
ks.test(log(d[["data"]]), "pnorm", 4, 1.2)
qqnorm(log(d[["data"]]))
abline(4,1,col="red")

d[["data"]] <- log(d[["data"]])
#voltar a variavel resposta sem transforma??es
#d[["data"]] <- P2$data
# estudo de tendencia direcional
points(d,l=1,pt.div="equal",col=gray(seq(1,0,l=11)),main="Post-Plot") 
points(d,pt.div="quartile",col=c("yellow","green","red","blue"),main="Post-Plot")
plot(d)
# ??? rever como preparar a entrada da spplot e bubble

#d2 <- SpatialPointsDataFrame(coords=d$coords,data=bf)
#d2 <- SpatialPointsDataFrame(coords=d$coords,data=as.data.frame(d2@data[,2]))
#spplot(d2,"Avistamento", do.log=T)
#bubble(d2, "Avistamento", dolog=TRUE, key.space="bottom")  # não ajuda em nada nesse caso
max(dist(d[[1]]))/2
x11()
Variog <- variog4(d, uvec=seq(0,18,l=11)  )
plot(Variog, legend = T)

#visualizando separadamente as 4 direcoes (util tambem definir 1 modelo e 
#adiciona-lo neste plot para cada direcao a fim checar anisotropia visualmente)

x11()
par(mfrow = c(2, 2))
plot(variog(d,uvec=seq(0,18,l=9), 
            dir=0, 
            tol=pi/8 ),main="N-S")
plot(variog(d,uvec=seq(0,18,l=9), 
            dir=pi/2, 
            tol=pi/8 ),main="E-W")
plot(variog(d,uvec=seq(0,18,l=9), 
        dir=3*pi/4, 
        tol=pi/8 ),main="SE-NW")
plot(variog(d,uvec=seq(0,18,l=9), 
            dir=pi/4, 
            tol=pi/8 ),main="SW-NE")

par(mfrow = c(1, 1))

alcances <- c(Variog$`0`$max.dist,Variog$`45`$max.dist,Variog$`90`$max.dist,Variog$`135`$max.dist)

#### Rodando o teste de anisotropia

mydt <- cbind(d[[1]],d[["data"]])

tol <- .5
mylags <- rbind(c(1, 0), c(0, 1), c(1, 1), c(-1, 1))  # representa as direcoes
myA <- rbind(c(1, -1, 0, 0), c(0, 0, 1, -1))
my.xlims <-  c(min(mydt[,1])-tol, max(mydt[,1])+tol)
my.ylims <-  c(min(mydt[,2])-tol, max(mydt[,2])+tol)
xlen <-  my.xlims[2]-my.xlims[1]
ylen <-  my.ylims[2]-my.ylims[1]
my.grid.spacing <-  c(xlen/16, ylen/12)

myh <-  0.70
myh.sb <-  0.85

tr.guan <- GuanTestUnif(spdata = mydt, lagmat = mylags, A = myA, df = 2, h = myh,
                        kernel = "norm", truncation = 1.5, xlims = my.xlims, ylims = my.ylims, 
                      grid.spacing = my.grid.spacing, window.dims = c(8,5), subblock.h = myh.sb)
tr.guan$p.value # segundo a GuanTestUnif com a configuração menos problematica possivel, o H0 e rejeitado por muito pouco

#arranjar os dados em formato matricial para aplicar median polish
bf.mat<-tapply(d[["data"]],list(factor(bf2$lati),factor(bf2$long)),function(x)x)

# necessário remover NA
bf.med<-medpolish(bf.mat,na.rm=T)
bf.trend.m<-bf2$Sightings-bf.med$resid

bf.trend<-na.omit(as.vector(bf.trend.m))

bf.res<-na.omit(as.vector(bf.med$res))
bf.res.g<-as.geodata(cbind(d[[1]],bf.res))
plot(bf.res)
#summary(bf.res)

#bf.finalvar<-variog(bf.res.g,max.dist=18,estimator.type="modulus")
#plot(bf.finalvar,xlab="Distance",ylab="Gamma",pch=20)

#bf.variog.d<-variog4(bf.res.g,max.dist=18,estimator.type="modulus")
#plot(bf.variog.d,lty=7,lwd=2)

#Remover tendência usando regressão
reg1<-lm(d[["data"]]~bf2$lati+bf2$long)#regressão nas coordenadas
summary(reg1) # talvez a longitude nao devesse estar inclusa

#reg2<-lm(d[["data"]]~bf2$lati)
#summary(reg2) # ou talvez devesse


bf.res<-reg1$res #extrair res e plotar
summary(bf.res)
bf.reg.g<-as.geodata(cbind(bf2@coords,bf.res)) 
plot(bf.reg.g)
#plot(d) #esse plot permite comparar os dados originais (onde apenas o log é aplicado a var. resposta)
#com os dados pós tentativa de correção de tendência. Aparentemente corrigir tendencia nao foi interessante (o semivariograma que o diga)

# o tipo de estimador foi alterado (classico -> Cressie), mas não causa mudança relevante
bf.res.reg.variog2<-variog4(bf.reg.g,max.dist=18,estimator.type="modulus",trend = "2nd")
x11()
plot(bf.res.reg.variog2)
title(main="Linear trend Residuals Directional Variograms")

# omnidirecional

d.var <- variog(d, uvec=seq(0,25,l=11),estimator.type="classical",pairs.min=30)
# visualizacao da semivariograma experimental ondimensional de Matheron 
plot(d.var, main= '')
d.var$u
d.var$v

#############################################################################
###### AJUSTE DO MODELO EXPONENCIAL POR ML   
#############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma') 
exp.ml<-likfit(d,ini=c(1.5,15),lambda=1, lik.method= "ML", cov.model="exp") 
#exp.ml 
lines(exp.ml,col="blue")
summary(exp.ml)


#############################################################################
###### AJUSTE DO MODELO EXPONENCIAL POR ML   
#############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma') 
wave.ml<-likfit(d,ini=c(1.5,15),lambda=1, lik.method= "ML", cov.model="wave") 
#exp.ml 
lines(exp.ml,col="blue")
summary(exp.ml)


#############################################################################
###### AJUSTE DO MODELO GAUSSIANO POR ML  
#############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma') 
gaus.ml<-likfit(d,ini=c(2,15),lambda=1, lik.method= "ML", cov.model="gaus",kappa=2.5) 
#gaus.ml 
lines(gaus.ml,col="blue")
summary(gaus.ml)


#############################################################################
###### AJUSTE DO MODELO Matern k=0.7 POR ML 
############################################################################ 

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma dmatern0.7.ml modificado') 
dmatern07.ml<-likfit(d,ini=c(0.5,0.1),lambda=1, method= "ML", cov.model="matern", kappa=0.7) 
dmatern07.ml
lines(dmatern07.ml,col="black")
summary(dmatern07.ml)

#############################################################################
###### AJUSTE DO MODELO Matern k=1 POR ML   
##############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma dmatern1.ml modificado') 
dmatern1.ml<-likfit(d,ini=c(0.5,0.1),lambda=1, method= "ML", cov.model="matern", kappa=1) 
dmatern1.ml
lines(dmatern1.ml,col="black")
summary(dmatern1.ml)

#############################################################################
###### AJUSTE DO MODELO Matern k=1.5 POR ML 
#############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma dmatern15.ml modificado') 
dmatern15.ml<-likfit(d,ini=c(0.5,0.08),lambda=1, method= "ML", cov.model="matern", kappa=1.5) 
dmatern15.ml
lines(dmatern15.ml,col="black")
summary(dmatern15.ml)


# atribuindo modelo ao semivariograma por OLS
# WAVE

plot(d.var,xlab='Distancia',ylab='semivariancia',main='Semivariograma ') 
dwave.ols <- variofit(d.var,ini=c(1,5),weights= "equal",cov.model="wave") 
#  Mostra o gr?fico das semivari?ncias com o modelo ajusto e o r?tulo.
lines(dwave.ols,col="blue")
summary(dwave.ols) 
# GAUSSIANO

plot(d.var,xlab='Distancia',ylab='semivariancia',main='Semivariograma ') 
dgaus.ols <- variofit(d.var,ini=c(1,15),weights= "equal",cov.model="gaus") 
#  Mostra o gr?fico das semivari?ncias com o modelo ajusto e o r?tulo.
lines(dgaus.ols,col="blue")
summary(dgaus.ols) 



# criando borda
points(d,l=1,pt.div="equal",col=gray(seq(1,0,l=11)),main="Post-Plot") 

contorno <- locator()
borda <- cbind(contorno$x,contorno$y)
polygon(borda)
#partindo para o processo de Krigagem Ordinaria
shift <- 5
gr<-expand.grid(x=seq(min(bf2$lati)-shift,max(bf2$lati)+shift,by=1),
                y=seq(min(bf2$long)-shift,max(bf2$long)+shift,by=1))
gi<-polygrid(gr,bor=borda)
plot(gr)
poly<-polygrid(gr,bor=borda)
points(poly,pch="+",col=2)
KC1=krige.control(obj=dwave.ols ,lam=1)
d.k1=krige.conv(d,loc=gi,krige=KC1)

#install.packages("classInt")
require(classInt)

#mapa com 5 classes de mesma amplitude
valores=d.k1$predict
interval_equal1=classIntervals(valores,5,style="equal",intervalClosure="right")
interval_equal1$brks
image(d.k1,loc=gr,border=borda,col=gray(seq(1,0,l=5)),breaks=c(interval_equal1$brks),main="") #automaticamente recebendo os valores dos intervalos


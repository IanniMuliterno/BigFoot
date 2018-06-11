
library(akima)
require(geoR)  
require(moments)
require(sp)
require(gstat)
library(readxl)
library(rgeos)
library(raster)
library(leaflet)
library(tibble)
library(spTest)

#put the path to your BigFoot data below

bf <- read_excel("~/GitHub/BigFoot/bigfoot USA.xlsx")

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
round(tst,3)

bf2 <- bf[,-c(1,3)] # retirando colunas 'obsoletas'
coordinates(bf2) <- cbind(bf2$lati,bf2$long)
bf2$Sightings <- as.numeric(bf2$Sightings)


d <- as.geodata(bf2)
#normalidade da vari?vel resposta
par(mfrow = c(1, 2))
mean(d[["data"]])
sd(d[["data"]])
qqnorm(bf2$Sightings,main="Normal Q-Q Plot for Original data")
ks.test(d[["data"]], "pnorm", 5, 8)
mean(log(d[["data"]]))
sd(log(d[["data"]]))
ks.test(log(d[["data"]]), "pnorm", 4, 1.2)
qqnorm(log(d[["data"]]),main="Normal Q-Q Plot for transformed data")
abline(4,1,col="red")
par(mfrow = c(1, 1))

d[["data"]] <-log(d[["data"]])
#voltar a variavel resposta sem transformacoes
#d[["data"]] <- P2$data
# estudo de tendencia direcional
points(d,l=1,pt.div="equal",col=gray(seq(1,0,l=11)),main="Post-Plot") 
points(d,pt.div="quartile",col=c("yellow","green","red","blue"),main="Post-Plot")
plot(d)
#max(dist(d[[1]]))/2
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
tr.guan$p.value 
# segundo a GuanTestUnif com a configuração menos problematica possivel,
# o H0 e rejeitado por muito pouco

#tentando melhorar a estética do semivariograma
#foram utilizadas Median Polish e regressão.
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

               
# o fato de não parecer existir um patamar no semivariograma inicial
# leva a crer que existe tendência, porem o median polish devolve zeros
               
               
#Remover tendência usando regressão
reg1<-lm(d[["data"]]~bf2$lati+bf2$long)#regressão nas coordenadas
summary(reg1) # talvez a longitude nao devesse estar inclusa # |
                                                             # |
 #reg2<-lm(d[["data"]]~bf2$lati)                             # |
#summary(reg2) # ou talvez devesse                           # |


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
# Nenhuma das duas tentativas devolveu algo melhor que o semivariograma dos dados transformados (log), voltamos a eles então
               
               
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
###### AJUSTE DO MODELO WAVE POR ML   
#############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma') 
wave.ml<-likfit(d,ini=c(1.5,5),lambda=1, lik.method= "ML", cov.model="wave") 
#exp.ml 
lines(wave.ml,col="blue")
summary(wave.ml)


#############################################################################
###### AJUSTE DO MODELO GAUSSIANO POR ML  
#############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma') 
gaus.ml<-likfit(d,ini=c(2,15),lambda=1, lik.method= "ML", cov.model="gaus") 
lines(gaus.ml,col="blue")
summary(gaus.ml)


#############################################################################
###### AJUSTE DO MODELO Matern k=0.7 POR ML 
############################################################################ 

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma dmatern0.7.ml modificado') 
dmatern07.ml<-likfit(d,ini=c(2,15),lambda=1, method= "ML", cov.model="matern", kappa=0.7) 
#dmatern07.ml
lines(dmatern07.ml,col="red")
summary(dmatern07.ml)

#############################################################################
###### AJUSTE DO MODELO Matern k=1 POR ML   
##############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma dmatern1.ml modificado') 
dmatern1.ml<-likfit(d,ini=c(2,15),lambda=1, method= "ML", cov.model="matern", kappa=1) 
#dmatern1.ml
lines(dmatern1.ml,col="black")
summary(dmatern1.ml)

#############################################################################
###### AJUSTE DO MODELO Matern k=1.5 POR ML 
#############################################################################

plot(d.var,xlab='Distancia',ylab='Semivariancia',main='Semivariograma dmatern15.ml modificado') 
dmatern15.ml<-likfit(d,ini=c(2,15),lambda=1, method= "ML", cov.model="matern", kappa=1.5) 
dmatern15.ml
lines(dmatern15.ml,col="black")
summary(dmatern15.ml)

#############################################################################
# Validacao Cruzada e Erro Absoluto do modelo 
#wave
vcwave=xvalid(d,model=wave.ml)
EAwave=sum(abs(vcwave$predicted-vcwave$data))
EAwave
summary(vcwave)
#exp
vcexp=xvalid(d,model=exp.ml)
EAexp=sum(abs(vcexp$predicted-vcexp$data))
EAexp
summary(vcexp)
#gaussiano
vcgaus=xvalid(d,model=gaus.ml)
EAgaus=sum(abs(vcgaus$predicted-vcgaus$data))
EAgaus
summary(vcexp)
#mat 0.7
vcmat07=xvalid(d,model=dmatern07.ml)
EAmat07=sum(abs(vcmat07$predicted-vcmat07$data))
EAmat07
summary(vcmat07)
#mat 1
vcmat1=xvalid(d,model=dmatern1.ml)
EAmat1=sum(abs(vcmat1$predicted-vcmat1$data))
EAmat1
summary(vcmat1)
#mat 1.5
vcmat15=xvalid(d,model=dmatern15.ml)
EAmat15=sum(abs(vcmat15$predicted-vcmat15$data))
EAmat15
summary(vcmat15)

modelos <- c("wave","exp","gaus","Mat_0.7","Mat_1","Mat_1.5")
EA_store <- c(EAwave,EAexp,EAgaus,EAmat07,EAmat1,EAmat15)
AIC_store <- c(wave.ml$AIC,exp.ml$AIC,gaus.ml$AIC,dmatern07.ml$AIC,
               dmatern1.ml$AIC,dmatern15.ml$AIC)
BIC_store <- c(wave.ml$BIC,exp.ml$BIC,gaus.ml$BIC,dmatern07.ml$BIC,
               dmatern1.ml$BIC,dmatern15.ml$BIC)
rbind(modelos,round(EA_store,4),round(AIC_store,4),round(BIC_store,4))

# Wave se saiu melhor
# mas o Matern com kappa 1.5 e o gaussiano podem ser considerados 
par(mfrow = c(1, 3))
# envelope com wave
env.wave<-variog.model.env(d,model.pars= wave.ml, obj.variog= d.var)  
plot(d.var,env= env.wave,xlab='distância',ylab='semivariancia',main='Envelope - Wave') 

# envelope com Mat 1.5
env.mat15<-variog.model.env(d,model.pars= dmatern15.ml, obj.variog= d.var)  
plot(d.var,env= env.mat15,xlab='distância',ylab='semivariancia',main='Envelope - Matern k=1.5') 

# envelope com gaus
env.gaus<-variog.model.env(d,model.pars= gaus.ml, obj.variog= d.var)  
plot(d.var,env= env.gaus,xlab='distância',ylab='semivariancia',main='Envelope - Gaussiano') 
par(mfrow = c(1, 1))


# criando borda
points(d,l=1,pt.div="equal",col=gray(seq(1,0,l=11)),main="Post-Plot") 

#COMENTADO PARA NAO ACABAR CRIANDO OUTRA BORDA POR CIMA DA ATUAL
#OS COMANDOS DEVEM SER RODADOS SE O SCRIPT ESTIVER SENDO EXECUTADO
#PELA PRIMEIRA VEZ

#contorno <- locator()
#borda <- cbind(contorno$x,contorno$y)
#write.table(borda,"border.txt")
polygon(borda)
#partindo para o processo de Krigagem Ordinaria
shift <- 5
gr<-expand.grid(x=seq(min(bf2$lati)-shift,max(bf2$lati)+shift,by=1),
                y=seq(min(bf2$long)-shift,max(bf2$long)+shift,by=1))
gi<-polygrid(gr,bor=borda)
plot(gr)
poly<-polygrid(gr,bor=borda)
points(poly,pch="+",col=2)
#d[["data"]] <- log(bf2$Sightings)
#d[["data"]] <- bf2$Sightings
#considerar usar ksline no lugar do krige.control
bf.UK <-ksline(d,locations=gi,cov.pars=dmatern15.ml$cov.pars,
       nugget=dmatern15.ml$nugget,cov.model="matern",kappa = 1.5,
       m0="kt",trend=1,ktedata=NULL,ktelocations=NULL)
# ou este, lembrando que esta configurado para krigagem ordinaria
#KC1=krige.control(obj=dmatern15.ml ,lam=1)
#d.k1=krige.conv(d,loc=gi,krige=KC1)
#install.packages("classInt")
require(classInt)

#mapa com 20 classes de mesma amplitude (40 para o mapa do desvio padrao)
valores <- bf.UK$predict
#valores=d.k1$predict

# para voltar a usar krig.cov e trabalhar com seu respectivo mapa retire os 
# os comentarios das linhas 387, 388 e 394 e substitua bf.UK por d.k1 nos comandos ab
interval_equal1=classIntervals(valores,5,style="equal",intervalClosure="right")
#interval_equal1$brks

image(bf.UK,loc=gr,border=borda,col=gray(seq(1,0,l=20)),
      val=valores,main="Mapa temático") #automaticamente recebendo os valores dos intervalos
#exp(bf.UK$)

#mapa do desvio padrao
image(bf.UK,loc=gr,border=borda,col=gray(seq(0,1,l=40)),
      val=sqrt(bf.UK$krige.var),main="Mapa do desvio padrão") #automaticamente recebendo os valores dos intervalos

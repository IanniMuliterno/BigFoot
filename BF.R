#library(kableExtra) # MEMO
#require(map)
#require(purrr)
#require(xtable)
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

bf <- read_excel("C:/Users/dapproplan2/Downloads/bigfoot USA.xlsx")

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

#normalidade da variável resposta
mean(d[["data"]])
sd(d[["data"]])
ks.test(d[["data"]], "pnorm", 5, 8)
mean(log(d[["data"]]))
sd(log(d[["data"]]))
ks.test(log(d[["data"]]), "pnorm", 4, 1.2)
qqnorm(log(d[["data"]]))
abline(4,1,col="red")

d[["data"]] <- log(d[["data"]])
#voltar a variavel resposta sem transformações
#d[["data"]] <- P2$data
# estudo de tendencia direcional
points(d,l=1,pt.div="equal",col=gray(seq(1,0,l=11)),main="Post-Plot") 
points(d,pt.div="quartile",col=c("yellow","green","red","blue"),main="Post-Plot")
plot(d)
# ??? rever como preparar a entrada da spplot e bubble

#d2 <- SpatialPointsDataFrame(coords=d$coords,data=bf)
#d2 <- SpatialPointsDataFrame(coords=d$coords,data=as.data.frame(d2@data[,2]))
#spplot(d2,"Avistamento", do.log=T)
#bubble(d2, "Avistamento", dolog=TRUE, key.space="bottom")  # nÃ£o ajuda em nada nesse caso
max(dist(d[[1]]))/2
x11()
Variog <- variog4(d, uvec=seq(0,max(dist(d[[1]]))/2,l=9))
plot(variog4(d, uvec=seq(0,max(dist(d[[1]]))/2,l=9)), legend = T)
alcances <- c(Variog$`0`$max.dist,Variog$`45`$max.dist,Variog$`90`$max.dist,Variog$`135`$max.dist)


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



d.var <- variog(d, uvec=seq(0,25,l=11),estimator.type="classical",pairs.min=30)
# visualizacao da semivariograma experimental ondimensional de Matheron 
plot(d.var, main= '')
d.var


plot(d.var,xlab='Distancia',ylab='semivariancia',main='Semivariograma ') 
dwave.ols <- variofit(d.var,ini=c(1,5),weights= "equal",cov.model="wave") 
#  Mostra o gr?fico das semivari?ncias com o modelo ajusto e o r?tulo.
lines(dwave.ols,col="blue")
summary(dwave.ols) 


plot(d.var,xlab='Distancia',ylab='semivariancia',main='Semivariograma ') 
dgaus.ols <- variofit(d.var,ini=c(1,5),weights= "equal",cov.model="gaus") 
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


##Establecer directorio y cargar librerias##
setwd()
library(adephylo)
library(ape)
library(BAMMtools)
library(betapart)
library(Biostrings)
library(CommEcol)
library(devtools)
library(dismo)
library(fossil)
library(geiger)
library(ggmap)
library(ggplot2)
library(gridExtra)
library(hisse)
library(kableExtra)
library(lattice)
library(latticeExtra)
library(mapdata)
library(maps)
library(maptools)
library(MASS)
library(matrixcalc)
library(Matrix)
library(mvMORPH)
library(OUwie)
library(PDcalc)
library(pez)
library(phylobase)
library(PhyloMeasures)
library(phyloregion)
library(phylosignal)
library(phytools)
library(picante)
library(raster)
library(rasterVis)
library(RColorBrewer)
library(reshape2)
library(rgbif)
library(rgdal)
library(rgeos)
library(sf)
library(shapefiles)
library(sp)
library(vegan)

###Para hacer un grid cuadrado de medio grado#####

shp<-shapefile("faja volcanica.shp")
gridsqu=st_make_grid(shp, square=TRUE,cellsize=0.25)
gridsqu=st_sf(geom=gridsqu)
st_write(gridsqu,"cuadro25.shp")

##Obtencion de la PAM###

map<-readOGR("cuadro25b.shp")
locs<-read.csv("anfibios_FVT.csv")
de cada uno de los campos del shapefile
coordinates(locs)<-~LON+LAT
crs(locs)<-crs(map)
class(locs)
class(map)
ovr<-over(locs,map)
frame<-data.frame(locs$Especie,locs$LAT,locs$LON,ovr$FID)
mat<-create.matrix(frame,tax.name = "locs.Especie", locality= "ovr.FID",abund=F)
write.csv(mat,"PAM25.csv")

Endemicidad
data<-read.csv("Pam500.csv")
matriz=data.matrix(data, rownames.force = NA)
sparse_comm <- matriz
dense_comm <- as.matrix(sparse_comm) 
object.size(dense_comm)
dense_comm
Calcular el indice de endemismo ponderado
endm <- weighted_endemism(dense_comm)
Convertir el vector de salida en una lista y guardarla
lista=as.list(endm)
lapply(lista, write, "endemismo cuadro.txt", append=TRUE, ncolumns=1000)
shp<-shapefile("cuadro 0.5.shp")
shp$cell<- as.numeric(shp$FID)
Endemismo<-as.data.frame(read.csv("endemismo cuadro.csv"))
shp2<-merge(shp, Endemismo, by.x='cell', by.y="cell")
writeSpatialShape(shp2,"Endemismo FAO.shp")
data<-read.csv("PAM_FAOendemismo.csv")
Cargar el arbol filogenetico
tree<-chronoMPL(read.tree("Phylo.tre"))

##Verificar que los nombres coincidan
Clean_tree<- match.phylo.comm(phy= tree, comm = data)$phy
Clean_data<- match.phylo.comm(phy = tree, comm = data)$comm
filoendemismo=phyloendemism(Clean_data, Clean_tree, weighted = TRUE)
Convertir el vector de salida en una lista y guardarla
lista=as.list(filoendemismo)
lapply(lista, write, "filoendemismo fao.txt", append=TRUE, ncolumns=1000)
shp<-shapefile("fao recorte.shp")
shp$cell<- as.numeric(shp$cell)
filoEndemismo2<-as.data.frame(read.csv("filoendemismo_fao2.csv"))
Se genera el shape de diversidad y se guarda el resultado 
shp2<-merge(shp, filoEndemismo2, by.x='cell', by.y="cell")
writeSpatialShape(shp2,"filoEndemismo FAO.shp")

###Calculo de diversidad alfa y beta###
data<-read.csv("PAM1000.csv")[,-(1)]
tree<-chronoMPL(read.tree("phylo.tre"))


Clean_tree<- match.phylo.comm(phy= tree, comm = data)$phy
Clean_data<- match.phylo.comm(phy = tree, comm = data)$comm


PD<-pd(samp= Clean_data, tree = Clean_tree, include.root = T)


test<-cor.test(PD$PD,PD$SR)
test
plot(PD$PD, PD$SR, xlab = "Phylogenetic Diversity", ylab = "Species Richness", pch = 16)

PDestandar=ses.pd(data, tree, null.model = "taxa.labels", runs = 99)
PDcorregida=PDestandar$pd.obs-PDestandar$pd.rand.mean
PDcorregida2=PDcorregida/PDestandar$pd.rand.sd
write.csv(PDcorregida2,file="PDcorregida.csv")

PDcorregida3=read.csv("PDcorregida.csv")


PDI<-pd.query(Clean_tree, Clean_data, standardize= T, null.model="uniform", reps=1000)
NRI<-mpd.query(Clean_tree, Clean_data, standardize=T, null.model="uniform", reps=1000)
NTI<-mntd.query(Clean_tree, Clean_data,standardize=T, null.model="uniform", reps=1000)


frame<-data.frame(PD$SR,PD$PD,PDI,NTI,NRI)
write.csv(frame, "alfa.csv")


shp<-shapefile("cuadro 0.5c.shp")


shp$cell<- as.numeric(shp$cell)
DIV<-as.data.frame(read.csv("alfa.csv "))
Se genera el shape de diversidad y se escribe/guarda el resultado 
shp2<-merge(shp, DIV, by.x='cell', by.y="cell")
writeSpatialShape(shp2,"div_alfa.shp")

###Calcular las matrices de diversidad beta###

Beta.taxo<-beta.pair(data, index.family = "sorensen")

Beta.taxo.sor<-as.matrix(Beta.taxo$beta.sor)
write.csv(Beta.taxo.sor, "Betataxo_sor_index.csv")

Beta.taxo.sim<-as.matrix(Beta.taxo$beta.sim)
write.csv(Beta.taxo.sim, "Betataxo_turn_index.csv")

Beta.taxo.nees<-as.matrix(Beta.taxo$beta.sne)
write.csv(Beta.taxo.nees, "Betataxo_ness_index.csv")

###Filogenetica####

beta_part<-phylo.beta.pair(Clean_data,Clean_tree,index.family = "sorensen")

Beta.filo.sor<-as.matrix(beta_part$phylo.beta.sor)
write.csv(Beta.filo.sor, "phyloBetataxo_sor_index.csv")

Beta.filo.turn<-as.matrix(beta_part$phylo.beta.sim)
write.csv(Beta.taxo.sim, "phyloBetataxo_turn_index.csv")

Beta.filo.nees<-as.matrix(beta_part$phylo.beta.sne)
write.csv(Beta.taxo.nees, "phyloBetataxo_ness_index.csv")

###Graficando beta###

mitotero<-function(data, radius, phylotree, phylobeta=F, index="sorensen"){
  mean_turnover<-numeric(length(data[,1]))
  mean_nestedness<-numeric(length(data[,1]))
  mean_beta<-numeric(length(data[,1]))
  for(i in 1:length(data[,1])){
    adj<-select.window(xf=data[i,1], yf=data[i,2], radius, xydata=data)[,-c(1,2)]
    if(phylobeta==F){
      ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-beta.pair(adj, index.family=index))
    }else if(phylobeta==T){
      ifelse(sum(nrow(adj))==0 || ncol(adj)==0, res<-0 , res<-phylo.beta.pair(adj, phylotree, index.family=index))
    }
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_turnover[i]<-0 , mean_turnover[i]<-mean(as.matrix(res[[1]])[2:length(as.matrix(res[[1]])[,1]),1],na.rm=TRUE) )
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_nestedness[i]<-0 , mean_nestedness[i]<-mean(as.matrix(res[[2]])[2:length(as.matrix(res[[2]])[,1]),1],na.rm=TRUE) )
    ifelse(sum(nrow(adj))==0 || ncol(adj)==0, mean_beta[i]<-0 , mean_beta[i]<-mean(as.matrix(res[[3]])[2:length(as.matrix(res[[3]])[,1]),1],na.rm=TRUE) )  
  }
  return(data.frame(cell=row.names(data), mean_turnover, mean_nestedness, mean_beta))
}
data <- read.table(choose.files(), row.names=1, head=T)
phylo <- read.tree("D:/Phylo.tre")

results <- mitotero(data,radius=1, index="sorensen")
Y filogenéticos
resultsphylo<- mitotero(data, radius=1, phylotree=phylo, phylobeta=T, index="sorensen")

frame<-data.frame(results$mean_beta, results$mean_turnover, results$mean_nestedness, resultsphylo$mean_beta, resultsphylo$mean_turnover, resultsphylo$mean_nestedness)
write.csv(frame,"DIV_BETA.csv ")
shp<-shapefile("D:/Cuencas FAO.shp")
shp$cell<- as.numeric(shp$cell)
DIV<-as.data.frame(read.csv("DIV_BETA.csv"))
shp2<-merge(shp, DIV, by.x='cell', by.y="cell")
writeSpatialShape(shp2,"DIV_beta.shp")

###Generar el centroide por celda
mapa=st_read("cuadro125gradosc.shp")
centroide=st_centroid(mapa)
coordenadas=st_coordinates(centroide)
write.csv(coordenadas,'centroide125.csv')
###Calcular la distancia geografica entre
los centroides

distanciaH=read.csv("centroide125.csv")
pts_sf <- st_as_sf(distanciaH, coords = c("lon", "lat"), crs = 4326)
matriz_distancia = st_distance(pts_sf)
matriz_km=matriz_distancia/1000

Mantelbetataxo=mantel(Beta.taxo$beta.sor, matriz_km, method = "pearson", permutations = 999)
print(Mantelbetataxo)

mantel_corrBetatotal <- mantel.correlog(Beta.taxo$beta.sor, matriz_km, nperm = 999)
plot(mantel_corrBetatotal)


###Calcular el factor de inflación de varianza para eliminar las variables correlacionadas###


library(usdm)
library(dplyr)
library(hier.part)
library("ggplot2")

##Cargar los datos para las regresiones####
T_resDiv=read.csv("recortadas.csv")


env <- T_resDiv[,4:11]

vif_results <- vifstep(env, th = 10)

vif_cor=vifcor(env,th=10)


##Cargar los datos para las regresiones##

T_resDiv=read.csv("recortadas.csv")

###Verificar que las variables cumplan el criterio de normalidad###

NormalDgeo=shapiro.test(T_resDiv$Betaphy.ness)
print(NormalDgeo)

##Graficar la forma de la distribución###
plot(density(T_resDiv$Beta.tax), main="Density Plot")

Grafica ern formato QQ para ver cuales son los puntos que desvian###
qqnorm(T_resDiv$Beta.tax, main = "Normal Q-Q Plot")
qqline(T_resDiv$Beta.tax, col = "red", lwd = 2)

#Primero, del dataframe que contiene los datos para las 

env <- T_resDiv[,3:10] 

###Asumiendo una distribucion gamma

gofs=all.regs(T_resDiv$Beta.tax, env, fam = "Gamma", gof = "RMSPE",
              print.vars = TRUE)

### una distribucion Gaussiana

gofsg=all.regs(T_resDiv$Beta.tax, env, fam = "gaussian", gof = "RMSPE",
               print.vars = TRUE)



prob=rand.hp(T_resDiv$Betaphy.ness, env, fam = "gaussian",
             gof = "logLik", num.reps = 999)$Iprobs

prob$sig95

part=partition(gofsg, pcan = 8, var.names = names(T_resDiv[,3,10]))
part$I.perc

part$gfs
part$IJ

gofsg$variable.combination
gofsg$gof



write.csv(part$I.perc,file="varianzaBetatax.csv")


betataxgraph=read.csv("varianzaBetatax.csv")

ggplot(betataxgraph,aes(fill=Variable, y=Variance,x=Estimator))+geom_bar(position="stack",stat="identity")+ scale_fill_manual(values=c("BIO 6"="chocolate","BIO 14"="red4","BIO 4D"="orangered","BIO 15D"="chocolate1","BIO 17D"="orange2","slope"="darkgoldenrod1","D Geo"="gold","soil temp"="khaki1"))+theme_bw()+
  theme(axis.text.x = element_text(angle = 90, size = 8, colour = "black", vjust = 0.5, hjust = 1), 
        axis.title.y = element_text(size = 8), legend.title = element_text(size = 8), 
        legend.text = element_text(size = 8, colour = "black"), 
        axis.text.y = element_text(colour = "black", size = 8))








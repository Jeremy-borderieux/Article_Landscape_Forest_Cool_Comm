####  GIS data  ####
GRECO_poly<-read_sf(file.path("Data","spatial_polygon"),layer="greco_l93")
eucoreg_poly<-read_sf(file.path("Data","spatial_polygon"),layer="map_fig4_boundaries")


#### all NFI datas ####


ListZip <- dir(path=file.path("Data","data_NFI"),pattern = "[1234567890]-fr.zip")
ListZip<-file.path(file.path("Data","data_NFI"),ListZip)
ListFichArb <- grep(pattern = "arbres_foret", unlist(lapply(ListZip, FUN = unzip, list = TRUE)), value = T)
ListFichCouvert <- grep(pattern = "couverts_foret", unlist(lapply(ListZip, FUN = unzip, list = TRUE)), value = T)
ListFichPlacette <- grep(pattern = "placettes_foret_20", unlist(lapply(ListZip, FUN = unzip, list = TRUE)), value = T)
ListFichEcologie <- grep(pattern = "ecologie_20", unlist(lapply(ListZip, FUN = unzip, list = TRUE)), value = T)


ListArbresVivants <- mapply(function(x, y) read.table(unz(x, y),sep=";",header=T), ListZip, ListFichArb, SIMPLIFY = F)
ListPlacette <- mapply(function(x, y) read.table(unz(x, y), header = TRUE, sep = ";", dec = ".", quote = "", encoding = "UTF-8",na.strings = c("","NULL")), ListZip, ListFichPlacette, SIMPLIFY = F)
ListCouvert <- mapply(function(x, y) read.table(unz(x, y), header = TRUE, sep = ";", dec = ".", quote = "", encoding = "UTF-8",na.strings = c("","NULL")), ListZip[grep("2005", ListZip, invert = T)], ListFichCouvert, SIMPLIFY = F)
ListEcologie <- mapply(function(x, y) read.table(unz(x, y), header = TRUE, sep = ";", dec = ".", quote = "", encoding = "UTF-8",na.strings = c("","NULL")), ListZip, ListFichEcologie, SIMPLIFY = F)

## dendrometrics variable
arbresVivants <- rbindlist(ListArbresVivants, fill=T); setkey(arbresVivants, "idp")
## location and general informations
placette <- rbindlist(ListPlacette, fill=T)
names(placette)[names(placette)=="tm2"]<-"tm_2"
# canopy cover variable
couvert <- rbindlist(ListCouvert, fill=T)
# ecological vairable
ecologie <- rbindlist(ListEcologie, fill=T)


## dendrometrics variable computation from individual trees 

setkey(arbresVivants, "espar");
arbresVivants <- droplevels(arbresVivants[c("2", "3", "4", "5", "6", "7", "9"), espar:=paste(0,espar,sep="")])

# see NFI protocole
arbresVivants[, dimess := cut(c13, breaks = c(0, 70.5, 117.5, 164.5, 1000), labels = c("PB", "BM", "GB", "TGB"), right = FALSE)]
arbresVivants[, htot := ifelse(is.na(htot), mean(htot, na.rm = TRUE), htot), by = c("idp", "espar", "dimess")]
arbresVivants[, ir5 := ifelse(is.na(ir5), mean(ir5, na.rm = TRUE), ir5), by = c("idp", "espar", "dimess")]

arbresVivants[,basal_area:=c13*c13*w/(40000*pi)]


dendro <- arbresVivants[,.(basal_area=sum(basal_area), Ntot=sum(w)),  by = idp]


## exact Elevation (because the cooridnates are blurred)
alti <- fread(file.path("Data","elevation_NFI.csv"))
alti <- alti[idp!="NULL"]
colnames(alti)<-c("idp","alti")
alti[,idp:=as.integer(idp)]


placette[ , greco := substr(ser, 1, 1)]
placette[,annee:=floor(idp/100000)+2005]
placette<-merge(placette,ecologie[,c("idp","dateeco")],by="idp")

rm(list=grep("List[:alnum:]*",ls(),value = T))



placette[,campagne := floor(idp/100000)+2005]


## normalization of the plot classification , .g. deforested
placette[,campagne := floor(idp/100000)+2005]
placette[,gest_struct:=sfo]
placette[,tmp:=switchv(as.character(sver),`NA`="0",`0`="0",`X`="0",`2`="1",`3`="4",`4`="2",`5`="3",`6`="1")]
placette$gest_struct[placette$campagne>2013]<-placette$tmp[placette$campagne>2013]
placette[,gest_struct:=switchv(as.character(gest_struct),`NA`=NA,`0`="debois",`1`="FR",`2`="FIR",`3`="TSF",`4`="T")]
placette[,gest_struct:=as.factor(gest_struct)]

### canopy cover vairable computation


couvert[,couverttot:=sum(tcl),by=.(idp)]


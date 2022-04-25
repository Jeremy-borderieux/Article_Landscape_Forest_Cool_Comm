#### Packages ####
for(pack in c("data.table","broman","stringr","ggplot2","ggspatial","sf","foreach","RColorBrewer","smoothr",
              "doParallel","purrr","raster","dplyr","car","sjPlot","mgcv","ggpubr","lubridate","patchwork","rcartocolor")){
  if(!pack %in%installed.packages()[,"Package"] )install.packages(pack)
  library(pack,character.only=TRUE)}

# We Use the Rtree package to improve calculation speed during the pairing, however this package is not mandatory
# 
#devtools::install_github("hunzikp/rtree")
library(rtree)


#### Loading the survey and covariable database  ####

## loading the NFI: survey location, canopy cover, dendrometric variable
setwd("~/Article_Landscape_Forest_Cool_Comm")
source("script_data_loading.R")
plots_data<-placette
plots_data[,ident:=as.character(idp)]# idp is the plot id of NFI plots

## climate_data contains the extraction of several climatic variable, calculated with the french meteorological station
climate_data<-readRDS(file.path("Data","floristic_climate_data","climatic_database.RData"))

## loading the survey database, which contain the harmonized floristic surveys of the NFI
surveys_data<-readRDS(file.path("Data","floristic_climate_data","floristic_survey_database.RData"))


## Loading the thermal optimum data from ClimPlant, with other traits from french Phytosociology databases and Ellenberg
sp_it_climplant<-fread(file.path("Data","CIT_data","climplant_names_trait.csv"))

# merging the traits with the surveys
surveys_data<-surveys_data[ident%in%plots_data$ident,]
surveys_data<-merge(surveys_data,sp_it_climplant[,c("lb_nom_final","YearMeanMean","YearMeanMedian","YearMean95","YearMean05","YearMaxMean","Area","indFor_Freq","indFor_Chytry","N_ellenberg","R_ellenberg","L_Ellenberg" ,"vi_pH" ,"vi_CN", "azote","topt_picq" , "topt" , "pHopt" , "CNopt"  , "Nopt" , "STopt","Li")],by.x="nom_espece",by.y="lb_nom_final",all.x=T)
surveys_data[,idp:=as.numeric(ident)]

# aggregating the survey to have a Community Inferred Temperature at the plot level (community)
cit_climplant<-surveys_data[tree!=1 & !nom_espece%in% c("Sambucus nigra","Sambucus racemosa","Ligustrum vulgare","Crataegus monogyna","Crataegus laevigata"),## remove of the tree species and the main woody species
                          .(cit_climplant=mean(YearMeanMean,na.rm=T),median_climplant=mean(YearMeanMedian,na.rm=T),
                            cit_tmax_climplant=mean(YearMaxMean,na.rm=T),
                            cit_climplant_05=mean(YearMean05,na.rm=T),
                            cit_climplant_95=mean(YearMean95,na.rm=T),
                            cit_ecoplant=mean(topt,na.rm=T),
                            cit_ecoplant_picq=mean(topt_picq,na.rm=T),
                            n_sp_climplant=sum(!is.na(YearMeanMean)),
                            mean_area=mean(Area,na.rm=T),
                            freq_for_mean=mean(indFor_Freq,na.rm=T),
                            indfor_Chytry=mean(indFor_Chytry,na.rm=T),
                            mean_azote=mean(azote,na.rm=T),
                            mean_N=mean(Nopt,na.rm=T),
                            mean_R=mean(R_ellenberg,na.rm=T),
                            mean_pH=mean(pHopt,na.rm=T),
                            mean_CN=mean(vi_CN,na.rm=T),
                            mean_L=mean(Li,na.rm=T)),by=idp]




plots_data$ident<-as.character(plots_data$idp)

## get elevation and climate data for the plot
plots_data<-merge(plots_data,climate_data,by="ident")
plots_data$idp<-as.numeric(plots_data$idp)
plots_data<-merge(plots_data,alti,by="idp",all.x=T)
plots_data[,alti:=ifelse(is.na(alti),elevation,alti)]

## get the CIt data
plots_data<-merge(plots_data,cit_climplant,by="idp",all.x=T)


## Landscape metrics calculated with a 20m meter resolution forest map of France (IGN), see main text
land_metrics_1000<-fread(file.path("Data","Landscape_data","Landscape_metrics_survey_1000m.csv"))
distance_edge<-fread(file.path("Data","Landscape_data","distance_edge_patch_metrics_2005_2011_nfi.csv"))

plots_data<-merge(plots_data,land_metrics_1000,by.x="ident",by.y="plot_id",all.x=T)
plots_data<-merge(plots_data,distance_edge,by="idp",all.x=T)

# pland= percentage of forested in a 1km buffer around the point
plots_data[,forest_type:=ifelse(pland_1000m>80,"Forested",ifelse(pland_1000m<30 ,"Not_Forested","nc"))]
plots_data<-plots_data[!is.na(forest_type),]# few plots didn't matched with the forest map at all



## get the dendrometrics and canopy cover variable
plots_data<-merge(plots_data,dendro,by="idp",all.x=T)
plots_data<-merge(plots_data,couvert[,.(couverttot =couverttot [1]),by=idp],by="idp",all.x=T)

plots_data<-plots_data[!is.na(tmoy_v2b_6186_13),]# some plots are outsid of the climatic grid
plots_data<-plots_data[!is.na(cit_climplant),]# some plots lacks any species with a thermal optimum
plots_data<-plots_data[gest_struct!="debois",]## the nfi classify deforested plots, we remove them



plots_data<-plots_data[!greco%in%c("J","H","K"),]# this study focuses the temperate biome
#we remove plots accoridng to the classification of Ecological regions of france = https://fr.wikipedia.org/wiki/Sylvoécorégion
plots_data<-plots_data[n_sp_climplant>=5,]# we want plots we enough species with a know thermal optimum

plots_data[,idp:=as.character(idp)]# the ID is transformed to a character, for consistency

# Create a spatial object
plots_data_sf<-st_as_sf(plots_data,coords=c("xl93","yl93"),crs=st_crs(2154))




# now that we have the final selection of NFI surveys, lets move on to the geographical pairing

#### functions for geographical pairing ####


get_distances<-function(ifn_sf,clust_factor,clust_level_1,clust_level_2,distance_tresh=6,use_rtree=FALSE){
  
  coord_clust_1<-ifn_sf[clust_factor==clust_level_1,]
  coord_clust_2<-ifn_sf[clust_factor==clust_level_2,]
  coord_clust_1<-coord_clust_1[,"idp"]
  coord_clust_2<-coord_clust_2[,"idp"]
  
  if(use_rtree){## Rtree fasten the process by exluding plots too far way from each other
 
  rtree_clust_1<-RTree(st_coordinates(coord_clust_1))
  
  close_neig<-withinDistance(rtree_clust_1, st_coordinates(coord_clust_2), distance_tresh*1000)
  
  coord_clust_1<-coord_clust_1[unique(unlist(close_neig)),]
  coord_clust_2<-coord_clust_2[,]
  }
  ## the function st _distance return a matrix of distances, transformed into a data.table
  dist_mat<-st_distance(coord_clust_2,coord_clust_1,tolerance=5000)
  units(dist_mat)<-"km"
  colnames(dist_mat)<-coord_clust_1$idp
  rownames(dist_mat)<-coord_clust_2$idp

  all_duet<-expand.grid(coord_clust_2$idp,coord_clust_1$idp,stringsAsFactors = FALSE)
  
  all_duet<-data.table(all_duet)
  all_duet[,dist_duet:=as.numeric(unlist(dist_mat))]
  all_duet_2_res<-all_duet[dist_duet<distance_tresh,]
  
  names(all_duet_2_res)<-c("idp2","idp1","dist_duet")
  
  all_duet_2_res<-all_duet_2_res[order(dist_duet),]
  all_duet_2_res<-all_duet_2_res[,c(2,1,3)]

  return(all_duet_2_res)
}

delete_dt_row <- function(DT, del.idxs) { ## useful function to fasten get_duet_optimal
  # delete row of potential instead of re-writing the datable into memory
  keep.idxs <- setdiff(DT[, .I], del.idxs)
  cols = names(DT)
  DT.subset <- data.table(DT[[1]][keep.idxs])
  setnames(DT.subset, cols[1])
  for (col in cols[2:length(cols)]) {
    DT.subset[, (col) := DT[[col]][keep.idxs]]
    DT[, (col) := NULL]
  }
  return(DT.subset)
}

get_centroid_duets<-function(sf_2rows){ res<-st_as_sf( st_centroid(st_cast(st_combine(sf_2rows),"LINESTRING")),crs=st_crs(2154))
  res$duet_id<-sf_2rows$duet_id[1]
  res}


get_duet_optimal<-function(sf_ifn,dist_table,distance_tresh=5,differences_to_compute=NULL,differences_tresh=if(is.null(differences_to_compute)) NULL else rep(0,length(differences_to_compute)),ident_is_num=T){
  
  dist_table<-dist_table[dist_duet<=distance_tresh,]# remove plots too far for each other:5km
  
  ## data preprocessing to get differences betwenn covariable
  dist_table[,duet_id:=1:nrow(dist_table)]
  subsest_idp1<-merge(dist_table[,c("idp1","duet_id")],as.data.table(sf_ifn),by.x="idp1",by.y="idp")
  subsest_idp2<-merge(dist_table[,c("idp2","duet_id")],as.data.table(sf_ifn),by.x="idp2",by.y="idp")
  subsest_idp1<-subsest_idp1[order(duet_id),]
  subsest_idp2<-subsest_idp2[order(duet_id),]
  
  ## compute the differences
  if(length(differences_to_compute)!=length(differences_tresh))stop("Length of columns to compute differences is different to the number of difference tresholds")
  for(col in differences_to_compute)dist_table[,paste0(col,"_dif")]<-subsest_idp1[[col]]-subsest_idp2[[col]]
  
  ## check for a differences criterion, here the elevation threshold only
  foreach(col = differences_to_compute ,trsh = differences_tresh)%do%{
    new_col<-paste0(col,"_dif")
    if(trsh!=0)dist_table<-dist_table[abs(dist_table[[new_col]])<= trsh,]
    
  }
  
  
  # this empty data.frame will keep the selected pairs (duets), and mimic the structure of dist_table
  kept_duet<-data.frame(matrix(NA,nrow = nrow(dist_table),ncol=ncol(dist_table)))
  
  colnames(kept_duet)<-colnames(dist_table)
  
  incr<-1
  while(nrow(dist_table)!=0){
    
    cat(paste("->",incr))
    
    # we use the function table to count the number od time each plot_id (called idp) are found, as it is also the number 
    # of neighbor it has
    nb_neig_idp1_dt<-as.data.table(table(dist_table$idp1))
    colnames(nb_neig_idp1_dt)<-c("idp1","nb_neig1")
    nb_neig_idp2_dt<-as.data.table(table(dist_table$idp2))
    colnames(nb_neig_idp2_dt)<-c("idp2","nb_neig2")
    
    # we match the number of neigboh with dist_table, that containt all the potential pairs
    dist_table[,nb_neig1:=nb_neig_idp1_dt[match(dist_table$idp1,idp1),2]]
    dist_table[,nb_neig2:=nb_neig_idp2_dt[match(dist_table$idp2,idp2),2]]
    
    ## this line check the minimum number of neigbor a plot in the pair has
    dist_table[,priority_low_neig:=pmin(nb_neig1,nb_neig2)]
    ## this line check the maximum number of neigbor a plot in the pair has
    dist_table[,priority_max_neig:=pmax(nb_neig1,nb_neig2)]
    
    # priority is given to the pair containing a plot with the lower minimal neighbor count
    ## then to the pair with the lower maximal count, then to the closest pair
    setorder(dist_table,priority_low_neig,priority_max_neig,dist_duet)
    
    # fast way to remove those columns now that the pairs have been prioritize
    dist_table[,nb_neig1:=NULL]
    dist_table[,nb_neig2:=NULL]
    dist_table[,priority_low_neig:=NULL]  
    dist_table[,priority_max_neig:=NULL]
    
    # the first row of dist_table is the selected one
    kept_duet[incr,]<-as.data.frame(dist_table[1,])
    
    incr<-incr+1
    kept_idp1<-dist_table[1,idp1]
    kept_idp2<-dist_table[1,idp2]
    
    ## we remove every potential pair with a plot that has just been selected
    dist_table<-delete_dt_row(dist_table,(1:nrow(dist_table))[dist_table[,idp1==kept_idp1]])
    dist_table<-delete_dt_row(dist_table,(1:nrow(dist_table))[dist_table[,idp2==kept_idp2]])
    
    
    
  }
  
  colnames(kept_duet)<-colnames(dist_table)
  
  kept_duet<-data.table(kept_duet)
  kept_duet<-kept_duet[!is.na(idp1),]## remove the empty rows, kept_duet was too large by design
  kept_duet<-kept_duet[order(dist_duet),]
  kept_duet[,duet_id:=1:nrow(kept_duet)]
  
  return(kept_duet)
}

#### geographical pairing ####

# get_duet_optimal and get_distances
# the results of these functions are saved in a .RData if you don't want to run the functions

# get_distances return all the potential pairs and the distance between the plots
duet_distances<-get_distances(plots_data_sf, # spatial object of interest (points)
                              plots_data_sf$forest_type, # the factor used to cluster the ponit in two classes
                              "Forested","Not_Forested", # the two classes
                              5.5,# distance threshold km
                              use_rtree = T) # use the package Rtree to faste nthe computation

# from all the potential pairs, get_duet_optimal return the maximum number of pairs, the constrain being that a plot can only be used once
duets<-get_duet_optimal(plots_data_sf,# same spatial object
                        duet_distances,# list of potential pairs and their distance, created with get_distances()
                        distance_tresh=5,# maximum distance of two points within a pair 
                        # compute pairwise differences within a pair of various variable
                        # here we compute forested- not_forested  (the order was set in get_distances())
                        differences_to_compute=c("alti","cit_climplant","tmoy_v2b_9015_13","mean_L","mean_pH","mean_pH","mean_N","annee","pland_1000m","couverttot","basal_area","Ntot","cit_ecoplant","cit_ecoplant_picq","distance_edge"),
                        differences_tresh =   c(50,rep(0,14)))# vector of same length as differences_to_compute, will exclude pairs with an absolute value of the difference higher than the threshold (0=no treshold)

summary(duets)

# load the saved computation of distances and pairing
list_distances<-readRDS(file.path("Data","saved_calculation","list_distances.RData"))
list_duets<-readRDS(file.path("Data","saved_calculation","list_duets.RData"))

duet_distances<-list_distances$`80`
duets<-list_duets$`80` # 80-100 %forested landscape in the buffer

summary(duets)


#### Analysis of the pairing , producing the results table and figure ####


melt_duets<-melt(duets,c(3:ncol(duets)))
melt_duets<-subset(melt_duets,select=- variable)


plots_data_duets_info<-merge(plots_data,melt_duets,by.x="idp",by.y="value")
plots_data_duets_info$clust_studied<-as.factor(plots_data_duets_info$forest_type)  # change the cluster studied here
plots_data_duets_info$clust_studied<-relevel(plots_data_duets_info$clust_studied,"Not_Forested")
plots_data_duets_info[,bioifn:=as.factor(ifelse(greco%in% c("D","E","G","H","I","K"),"montagne",ifelse(greco%in% c("A","B","C","F"),"plaine","mediterranee")))]
plots_data_duets_info$bioifn<-relevel(plots_data_duets_info$bioifn,"plaine")
plots_data_duets_info<-plots_data_duets_info[order(duet_id,clust_studied),]


plots_data_duets_info$dendro_info_duets<-plots_data_duets_info[,rep(sum(!is.na(basal_area)),2),by=duet_id]$V1
plots_data_duets_info$couvert_info_duets<-plots_data_duets_info[,rep(sum(!is.na(couverttot)),2),by=duet_id]$V1
plots_data_duets_info$edge_info_duets<-plots_data_duets_info[,rep(sum(!is.na(distance_edge)),2),by=duet_id]$V1


## figure 2: histogram of the differences within a pair
duets[,Forested_landscape:=ifelse(cit_climplant_dif<0,"Cooler CIT","Warmer CIT")]
fig_2_out<-ggplot(duets,aes(x= cit_climplant_dif,fill= Forested_landscape))+theme_bw()+
  geom_histogram(color="black",binwidth = 0.25,origin=0)+
  geom_vline(xintercept = mean(duets$cit_climplant_dif),col="black",lwd=1,linetype=2)+geom_segment(aes(x=0,xend=0,y=0,yend=254),col="grey20",lwd=1)+
  scale_fill_manual(values=c("lightblue3","lightcoral"))+
  ylab("Number of plot pairs")+xlab("∆CIT: Forested CIT  -  Non forested CIT ")+labs(fill = "Forested \nlandscape")+scale_y_continuous(limits=c(0,300))+scale_x_continuous(labels = paste0(c(-4,-2,0,2,4) , "°C"))+
  annotate("text", x = -1.25, y = 295, label = paste0("Mean ∆CIT",round(mean(duets$cit_climplant_dif),2),"°C"))


reso<-2
tiff(file.path("results_and_figures","fig_2_hist.tif"),width=725*reso,height=450*reso,res=100*reso)

fig_2_out

dev.off()

## small numbers for result core text main text
mean(duets$cit_climplant_dif)
sd(duets$cit_climplant_dif)
sum(duets$cit_climplant_dif<0)/nrow(duets)
sum(duets$cit_climplant_dif< (-1))/nrow(duets)

wilcox.test(duets$cit_climplant_dif)




## linear model for analysis
#number of pairs for the main and the submodels
n_full_model<-nrow(duets)
n_cover_model<-nrow(duets[!is.na(couverttot_dif),])
n_edge_model<-nrow(duets[!is.na(distance_edge_dif),])

scope<-c("cit_climplant_dif","tmoy_v2b_9015_13_dif","mean_L_dif","mean_pH_dif","mean_N_dif","annee_dif","alti_dif")

lm_duets_pair<-lm(cit_climplant_dif~1,data=duets[,..scope],y=T,x=T)
lm_duets_pair_2<-lm(cit_climplant_dif~.,data=duets[,..scope],y=T,x=T)

step_pair<-step(lm_duets_pair,direction="both",scope=list(upper=lm_duets_pair_2,lower=lm_duets_pair))
summary(step_pair)
## all the variable are included, we run the submodels with all the variable, and the addition of canopy cover or distance to the edge
lm_duets_pair<-lm(cit_climplant_dif~alti_dif+tmoy_v2b_9015_13_dif+mean_L_dif+mean_pH_dif+mean_N_dif+annee_dif,data=duets,y=T,x=T)
lm_cover<-lm(cit_climplant_dif~alti_dif+tmoy_v2b_9015_13_dif+mean_L_dif+mean_pH_dif+mean_N_dif+annee_dif+couverttot_dif,data=duets,y=T,x=T)
lm_edge<-lm(cit_climplant_dif~alti_dif+tmoy_v2b_9015_13_dif+mean_L_dif+mean_pH_dif+mean_N_dif+annee_dif+distance_edge_dif,data=duets,y=T,x=T)


summ<-summary(lm_duets_pair)
summ
coef(summ)

## checking the residuals
vif(step_pair)
vif(lm_duets_pair)
res_pairs<-residuals(lm_duets_pair)
hist(res_pairs,nc=40,col="lightgrey")
duets$residuals<-res_pairs


ggplot(duets,aes(y=cit_climplant_dif,x=residuals+coef(lm_duets_pair)[1]))+theme_bw()+geom_point(size=0.5)+geom_smooth(method="lm")
ggplot(duets,aes(x=dist_duet,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth()
ggplot(duets,aes(x=alti_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth()
ggplot(duets,aes(x=tmoy_v2b_9015_13_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth(method="lm")
ggplot(duets,aes(x=mean_L_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth()
ggplot(duets,aes(x=mean_pH_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth()
ggplot(duets,aes(x=mean_N_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth(method="lm")
ggplot(duets,aes(x=annee_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth()
ggplot(duets,aes(x=basal_area_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth()
ggplot(duets,aes(x=couverttot_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth()
ggplot(duets,aes(x=distance_edge_dif,y=residuals))+theme_bw()+geom_point(size=0.5)+geom_smooth(method="lm")


## some mean of differents environementale variables, for F and NF plots, to show the sampling balance
plots_data_duets_info[,mean(cit_climplant),by=.(clust_studied)]

plots_data_duets_info[,mean(tmoy_v2b_9015_13),by=.(clust_studied)]
plots_data_duets_info[,quantile(tmoy_v2b_9015_13,probs=c(0.05,0.95)),by=.(clust_studied)]
plots_data_duets_info[,mean(prec_v2b_9015_4_9),by=.(clust_studied)]

plots_data_duets_info[,mean(tmoy_wc_7100_13),by=.(clust_studied)]
plots_data_duets_info[,mean(pp_wc_7100_13),by=.(clust_studied)]

plots_data_duets_info[,mean(alti),by=.(clust_studied)]

plots_data_duets_info[,mean(pland_1000m),by=.(clust_studied)]
plots_data_duets_info[,mean(ca_1000m),by=.(clust_studied)]


plots_data_duets_info[dendro_info_duets==2,mean(basal_area),by=.(clust_studied)]
plots_data_duets_info[dendro_info_duets==2,mean(Ntot),by=.(clust_studied)]

plots_data_duets_info[couvert_info_duets==2,mean(couverttot),by=.(clust_studied)]

plots_data_duets_info[dendro_info_duets==2,quantile(basal_area,probs=c(0.05,0.95)),by=.(clust_studied)]
plots_data_duets_info[dendro_info_duets==2,quantile(Ntot,probs=c(0.05,0.95)),by=.(clust_studied)]

plots_data_duets_info[couvert_info_duets==2,quantile(couverttot,probs=c(0.05,0.95)),by=.(clust_studied)]


plots_data_duets_info[,mean(mean_N),by=.(clust_studied)]
plots_data_duets_info[,mean(mean_pH),by=.(clust_studied)]
plots_data_duets_info[,mean(mean_CN,na.rm=T),by=.(clust_studied)]
plots_data_duets_info[,mean(mean_L,na.rm=T),by=.(clust_studied)]
plots_data_duets_info[,mean(n_sp_climplant,na.rm=T),by=.(clust_studied)]

plots_data_duets_info[,mean(n_sp_climplant,na.rm=T),by=.(clust_studied)]


## table1 for mat_met
summary_data<-function(x,round){
  m<-round(mean(x),round)
  q<-round(quantile(x,c(0.05,0.95)),round)
  
  return(paste(m," (",q[1],"-",q[2],")",sep="",collaspe=""))
}



table_1_dsummary<-plots_data_duets_info[,.(n=.N,tmoy=summary_data(tmoy_v2b_9015_13,2),alti=summary_data(alti,0),ca=summary_data(ca_1000m,0),pland=summary_data(pland_1000m,1),N_ell=summary_data(mean_N,1),pH_vi=summary_data(mean_pH,1)),by=.(clust_studied)]
table_1_dsummary<-cbind(table_1_dsummary,plots_data_duets_info[dendro_info_duets==2,.(gtot=summary_data(basal_area,1)),by=.(clust_studied)][,2])
table_1_dsummary<-cbind(table_1_dsummary,plots_data_duets_info[couvert_info_duets==2,.(couvert=summary_data(couverttot,1)),by=.(clust_studied)][,2])
table_1_dsummary
name<-data.table(`Landscape of the plot`=c("Not forested","Forested"))
final<-cbind(name,table_1_dsummary)
final$clust_studied<-NULL
write.table(final,file.path("results_and_figures","table_1_mat_met.csv"),sep=";",row.names = F)


## main text mat met, nb of species with topt and indicator value
n_sp_biodind<-surveys_data[ident%in%plots_data_duets_info$idp  & tree==0,]
length(unique(n_sp_biodind[!is.na(YearMeanMean),nom_espece]))
length(unique(n_sp_biodind[!is.na(pHopt),nom_espece]))
length(unique(n_sp_biodind[!is.na(Nopt),nom_espece]))
length(unique(n_sp_biodind[!is.na(Li),nom_espece]))


n_sp_biodind[,.(n_oc_topt=sum(!is.na(YearMeanMean))/.N,
                n_oc_pH=sum(!is.na(pHopt))/.N,
                n_oc_N=sum(!is.na(Nopt))/.N,
                n_oc_L=sum(!is.na(Li))/.N),]

summary(n_sp_biodind[ ,.(n_oc_topt=sum(!is.na(YearMeanMean)),
                                         prop_topt=sum(!is.na(YearMeanMean))/.N,
                                         n_oc_pH=sum(!is.na(pHopt)),
                                         n_oc_N=sum(!is.na(Nopt)),
                                         n_oc_L=sum(!is.na(Li))),by=ident])# 14 sp in climplant on average

## table 2 : results

make_result_table<-function(lm){
  if(is.null(lm$x))stop( "use lm(...,x=TRUE)")
  summ<-summary(lm)
  
  table_model<-data.table(coef(summ)[,c(1,2,4)])
  table_model[,`P-value`:=ifelse(`Pr(>|t|)`<10e-4,"<10-4",as.character(signif(`Pr(>|t|)`,2)))]
  table_model<-cbind(data.table(estimates=names(coef(lm))),table_model)
  
  mean_pred<-apply(lm$x, 2, mean)
  table_model[,Effect_size:=mean_pred*Estimate]
  table_model$`Pr(>|t|)`<-NULL
  
  table_model[,Estimate:=signif(Estimate,3)]
  table_model[,`Std. Error`:=signif(`Std. Error`,2)]
  table_model[,Effect_size:=signif(Effect_size,2)]
  
  table_model[,`R²`:=signif(summ$r.squared,3)]
  return(table_model)
  
}


table_model<-make_result_table(lm_duets_pair)
table_model$`P-value`<-str_replace_all(table_model$`P-value`,"[.]",",")



write.table(table_model,file.path("results_and_figures","table_2_resu_model.csv"),dec = ",",sep=";",row.names = F)

table_suppl_cover<-make_result_table(lm_cover)
table_suppl_cover$`P-value`<-str_replace_all(table_suppl_cover$`P-value`,"[.]",",")
table_suppl_cover$`R²`<-NULL

table_suppl_edge<-make_result_table(lm_edge)
table_suppl_edge$`P-value`<-str_replace_all(table_suppl_edge$`P-value`,"[.]",",")
table_suppl_edge$`R²`<-NULL

write.table(table_suppl_cover,file.path("results_and_figures","table_supplementary_1_resu_model.csv"),dec = ",",sep=";",row.names = F)
write.table(table_suppl_edge,file.path("results_and_figures","table_supplementary_2_resu_model.csv"),dec = ",",sep=";",row.names = F)


#### changing scale analysis, fig 3 ####

buffers<-c(30,35,40,45,50,55,60,65,70,75,80)
for(i in buffers ){
  iplus<-i+20
  plots_data[,paste0("forest_type_",i,"_",iplus)]<-ifelse(plots_data_sf$pland_1000m> i & plots_data_sf$pland_1000m<= iplus ,"Forested",ifelse(plots_data_sf$pland_1000m<30  ,"Not_Forested","nc"))
  
  plots_data_sf[,paste0("forest_type_",i,"_",iplus)]<-ifelse(plots_data_sf$pland_1000m> i & plots_data_sf$pland_1000m<= iplus ,"Forested",ifelse(plots_data_sf$pland_1000m<30  ,"Not_Forested","nc"))
  
  
}


### the two foreach loop reproduce the calcultation of the distances and the duets pairing for a lot of different
# percentages of forested landscape, the results are saved in the same .RData we saw earlier if you don't want to run these lines

cluster_duets_vary<-makeCluster(4)
registerDoParallel(cluster_duets_vary)

list_distances<-foreach(buff=buffers,.combine = list,.multicombine = T,.packages = c("data.table","stringr","foreach","broman","sf","rtree"))%dopar%{
  buff_plus<-buff+20
  distances<-get_distances(plots_data_sf,as.factor(unlist(st_drop_geometry(plots_data_sf[,paste0("forest_type_",buff,"_",buff_plus)]))),"Forested","Not_Forested",6,use_rtree=T)
  cat(paste0(buff_plus," - "))
  distances
  
  
}

names(list_distances)<-buffers



list_duets<-foreach(dist=list_distances,.combine = list,.multicombine = T,.packages = c("data.table","stringr","foreach","broman"))%dopar%{
  
  
  res<-get_duet_optimal(plots_data_sf,dist,distance_tresh=5,differences_to_compute=c("alti","cit_climplant","tmoy_v2b_9015_13","mean_L","mean_pH","mean_pH","mean_N","annee","pland_1000m","couverttot","basal_area","Ntot","cit_ecoplant","cit_ecoplant_picq"),differences_tresh =   c(50,rep(0,13)))
  res
  
}

names(list_duets)<-buffers

stopCluster(cluster_duets_vary)


#saveRDS(list_distances,file.path("Data","saved_calculation","list_distances.RData"))
#saveRDS(list_duets,file.path("Data","saved_calculation","list_duets.RData"))

list_distances<-readRDS(file.path("Data","saved_calculation","list_distances.RData"))
list_duets<-readRDS(file.path("Data","saved_calculation","list_duets.RData"))



### this function reproduce the analysis for each new sets of duets
analyse_one_duet_subset<-function(x,buff_min){
  lm_duets<-lm(cit_climplant_dif~alti_dif+tmoy_v2b_9015_13_dif+mean_L_dif+mean_pH_dif+mean_N_dif+annee_dif,data=x,y=T,x=T)
  topt_dif<-mean(x$cit_climplant_dif)
  var_landscape_dif<-mean(x$pland_1000m_dif)
  resu<-make_result_table(lm_duets)
  dt_eff_size<-t(data.table(resu$Effect_size))
  colnames(dt_eff_size)<-paste0("Eff_size_",resu$estimates)
  
  whole_buffer<-mean(x$cit_climplant_dif)/mean(x$pland_1000m_dif)
  forested_ls_buffer<-coef(summary(lm_duets))[1,1]/mean(x$pland_1000m_dif)
  
  
  
  resu<-data.table(buff_min=buff_min,topt_dif=topt_dif,pland_dif=var_landscape_dif,n_duets=nrow(x),r_squared=resu$`R²`[1],p_value_landscape=resu[1,"P-value"],whole_buffer=whole_buffer,forested_ls_buffer=forested_ls_buffer)
  resu<-cbind(resu,dt_eff_size)
  return(resu)
  
}



res_dif<-foreach(buff=buffers,duets_different_buffer=list_duets,.combine = rbind)%do%{
  
  multi_scale_results<-analyse_one_duet_subset(duets_different_buffer,buff)
  
  multi_scale_results
}

res_dif[,buffer_width:=paste0("[",buff_min," : ",buff_min+20,"]")]
res_dif[,coef_thermo:=topt_dif/pland_dif]
res_dif[,coef_land:=`Eff_size_(Intercept)`/pland_dif]
summary(res_dif)

### figure 3 : with different scale



fig_3_out<-ggplot(res_dif,aes(x=buffer_width))+theme_bw()+
  geom_line(aes(y=Eff_size_mean_N_dif,color="Δ Ellenberg N",group=1),lwd=1.25)+geom_point(aes(y=Eff_size_mean_N_dif,fill="Δ Ellenberg N"),color="black",size=3,shape=24)+
  geom_line(aes(y=Eff_size_mean_pH_dif,color="Δ Bioindicated pH",group=1),lwd=1.25)+geom_point(aes(y=Eff_size_mean_pH_dif,fill="Δ Bioindicated pH"),color="black",size=3,shape=25)+
  geom_line(aes(y=`Eff_size_(Intercept)`,color="Forested landscape",group=1),lwd=1.25)+ geom_point(aes(y=`Eff_size_(Intercept)`,fill="Forested landscape"),color="black",size=3,shape=25)+
  geom_line(aes(y=topt_dif,colour="Δ CIT",group=1),lwd=1.25)+ geom_point(aes(y=topt_dif,fill="Δ CIT"),size=3)+
  geom_line(aes(y=Eff_size_alti_dif,colour="Marginal effects",group=1),lwd=1,linetype=2)+
  geom_line(aes(y=Eff_size_mean_L_dif,colour="Marginal effects",group=1),lwd=1,linetype=2)+
  geom_line(aes(y=Eff_size_annee_dif,colour="Marginal effects",group=1),lwd=1,linetype=2)+
  geom_line(aes(y=Eff_size_tmoy_v2b_9015_13_dif,colour="Marginal effects",group=1),lwd=1,linetype=2)+
  geom_hline(yintercept = 0,lwd=1.25,lty=2,col="grey40")+
  scale_color_manual(values=c("chocolate2", "black","cornflowerblue" ,"chartreuse4","grey80"))+
  scale_fill_manual(values=c("chocolate2", "black","cornflowerblue" ,"chartreuse4","grey80"))+
  guides(fill = FALSE)+
  geom_text(aes(y=pmin(Eff_size_mean_pH_dif,topt_dif),label=paste0("n= ",n_duets)),nudge_y=-0.035 ,nudge_x=-0.05,angle=25)+
  labs(x="Forest amount in the 1km buffer to consider a landscape forested",y="Δ CIT",colour="Contribution of \neach effects to \nΔ CIT ")+
  theme(axis.text.x = element_text(size=12, angle=45,colour="black", hjust = 1),
        axis.text.y = element_text(size=12,color="black"),
        axis.title.x = element_text(size=13,color="black"),
        axis.title.y = element_text(size=13,color="black"))

fig_3_out

reso<-2
tiff(file.path("results_and_figures","fig_3.tif"),width =955*reso ,height =570*reso ,res =100*reso)

fig_3_out

dev.off()


#### maps for representation, fig 1#### 


 plots_data_duets_info_sf<-st_as_sf(plots_data_duets_info,coords=c("xl93","yl93"),crs=st_crs(2154))#
## getting the centroids of each pairs
get_centroid_duets<-function(sf_2rows){ res<-st_as_sf( st_centroid(st_cast(st_combine(sf_2rows),"LINESTRING")),crs=st_crs(2154))
 res$duet_id<-sf_2rows$duet_id[1]
 res}

 all_centroid<-foreach(duet=sort(unique(plots_data_duets_info_sf$duet_id)),.combine=rbind)%do%{get_centroid_duets(plots_data_duets_info_sf[plots_data_duets_info_sf$duet_id==duet,])}


## loading France forest cover at 1km scale for representation
rasterize_forest_2<-raster(file.path("Data","spatial_raster","forest_cover_france.tif"))
## Using the french ecological region as a border for basemap
GRECO_temperate_forest<-GRECO_poly[!GRECO_poly$GRECO%in%c("-","H","K","J"),]
rasterize_forest_2<-mask(rasterize_forest_2,as_Spatial( st_union(GRECO_temperate_forest)))


## maps, using ggspatial and opens street map basemaps
raster_dt<-data.table(xyFromCell(rasterize_forest_2,1:ncell(rasterize_forest_2)) )
raster_dt$forested<-getValues(rasterize_forest_2)


out_map_1<-ggplot(all_centroid)+theme_bw()+  annotation_map_tile(type="thunderforestlandscape",zoomin = 0) + 
  geom_raster(raster_dt,mapping=aes(x=x,y=y,fill=forested))+scale_fill_gradient(low="cornsilk",high ="chartreuse4" ,na.value = "transparent")+
  geom_sf(size=0.9)+geom_sf(data=st_union(GRECO_temperate_forest),color="black",fill=NA)+ labs(x="",y="")+guides(fill = FALSE)+
  annotation_north_arrow(location = "tl")+annotation_scale(location = "bl", text_cex = 1.3)

out_map_2<-ggplot(all_centroid)+theme_bw()+  annotation_map_tile(type="cartolight",zoomin = 0) + 
  geom_raster(raster_dt,mapping=aes(x=x,y=y,fill=forested))+scale_fill_gradient(low="cornsilk",high ="chartreuse4" ,na.value = "transparent")+
  geom_sf(size=0.9)+geom_sf(data=st_union(GRECO_temperate_forest),color="black",fill=NA)+ labs(x="",y="")+guides(fill = FALSE)+
 annotation_north_arrow(location = "tl")+annotation_scale(location = "bl", text_cex = 1.3)


## changing the scale, allowing the figure 1 to be smaller
out_map_1_zoomed<-out_map_1+theme(axis.text = element_text(size=13))
out_map_2_zoomed<-out_map_2+theme(axis.text = element_text(size=13))


## this function is used to create the fig1.b zoomed landscape, to illustrate a pair

Visualize_duets_landscape<-function(dt_2_points,raster_forest,buffers=c(250,1000,2500,5000),use_osm=F,style_map="osm",extent_x_change=0,extent_y_change=0,legend.position="right"){
  
  dt_2_points[,`Plot classification`:=ifelse(clust_studied=="Not_Forested","Not forested","Forested")]
  
  points_2_sf<-st_as_sf(dt_2_points,coords=c("xl93","yl93"),crs=st_crs(2154))
  
  sf_buffer<-NULL
  
  for(radius in buffers ){
    bu<-st_buffer(points_2_sf,radius)
    bu$radius<-radius
    sf_buffer<-rbind(sf_buffer,bu)
    
  }
  ext<-bbox(as_Spatial( st_buffer(sf_buffer,max(buffers)*2.5)))
  ext<-extent(ext)
  ext<-extent(ext[1]-(ext[2]-ext[1])*extent_x_change,ext[2]+(ext[2]-ext[1])*extent_x_change,ext[3] - (ext[4]-ext[3])*extent_y_change, ext[4] +(ext[4]-ext[3])*extent_y_change)
  ext<-ext-1000


  forest<-raster_forest
  val<-getValues(forest)
  xy<-data.table(xyFromCell(forest,1:ncell(forest)))
  xy<-cbind(xy,val)
  if(use_osm) xy[,val:=ifelse(val==0,NA,1)]
  p<-ggplot(data=xy)

  if(use_osm) forest_color<-c("chartreuse4") else forest_color<-c("cornsilk","chartreuse4")
  
  if(!use_osm){ p<-p+geom_raster(aes(x=x,y=y,fill=as.factor(val)),show.legend=F)+scale_fill_manual(values=forest_color,na.value = "transparent")+theme_bw()
  p<-p+geom_sf(data=points_2_sf,aes(color=`Plot classification`),size=2.3)+scale_color_manual(values=c("tan4","darkgoldenrod1"))
  p<-p+geom_sf(data=sf_buffer,aes(color=`Plot classification`),size=1.2,show.legend=F,fill=NA,alpha=0.2)+guides(color = guide_legend(override.aes = list(size = 4)))+
    p<-p+annotation_scale(location = "bl")+ coord_sf(xlim=c(ext[1],ext[2]),ylim=c(ext[3],ext[4]),datum = st_crs(4326))+labs(x="",y="")+theme(legend.position =legend.position,legend.direction = ifelse(legend.position%in%c("bottom","top"),"horizontal","vertical") )
  p
  }else{
    p2<-ggplot(points_2_sf)+annotation_map_tile(type=style_map,zoom =   13)+geom_raster(data=xy,aes(x=x,y=y,fill=as.factor(val)),show.legend=F)+scale_fill_manual(values=forest_color,na.value = "transparent")+theme_bw()+
      geom_sf(mapping=aes(color=`Plot classification`),size=2.3)+scale_color_manual(values=c("tan4","darkgoldenrod1"))+
      geom_sf(data=sf_buffer,aes(color=`Plot classification`),size=1.2,show.legend=F,fill=NA,alpha=0.2)+guides(color = guide_legend(override.aes = list(size = 4)))+
      annotation_scale(location = "bl",text_cex =1.3)+ coord_sf(xlim=c(ext[1],ext[2]),ylim=c(ext[3],ext[4]),datum = st_crs(4326))+labs(x="",y="")+theme(legend.position =legend.position,legend.direction = ifelse(legend.position%in%c("bottom","top"),"horizontal","vertical"),axis.text = element_blank(),axis.ticks = element_blank() )
    p2}
  
}


# the raster is large, we use a small portion of it for reproducibility
sample_raster_fig_1<-raster(file.path("Data","spatial_raster","fig_1_landscape.tif"))

# we use a selected pair for display, as real coordinates of NFI are private
dt_fig1_b<-data.table(xl93=c(517606,513803),yl93=c(6787715,6789970),clust_studied=c("Not_Forested","Forested"),duet_id=1)

out_final_map3_zoom<-Visualize_duets_landscape(dt_fig1_b,sample_raster_fig_1,buffers=c(1000),
                          use_osm = T,style_map = "osm",
                          extent_x_change = -0.05,extent_y_change = -0.15,
                          legend.position = "bottom")

out_final_map3_zoom<-out_final_map3_zoom+theme(legend.text = element_text(size=14),legend.title = element_text(size=16,face="bold"))

square_inset<-st_buffer(get_centroid_duets(st_as_sf(dt_fig1_b,coords = c("xl93","yl93"),crs=st_crs(2154))),endCapStyle = "SQUARE",dist=10000)

## arrange fig1.a and fig.b
out_map_inset_zommed<-ggarrange(out_map_2_zoomed+geom_sf(data=square_inset,fill=NA,lwd=1,col="brown"),out_final_map3_zoom,nrow=2,align = "v",heights = c(1.5,1),hjust = 0,labels = c("(a)    ","(b)    "),font.label = list(size = 26, color = "black", face = "bold", family = NULL))


tiff(file.path("results_and_figures","map_fig_1_ggarange_smaller.tif"),width =650*reso /2,height =1150*reso /2,res =100*reso/2)
out_map_inset_zommed

dev.off()



#### Mapping the effect across Europe, fig 4####


## calculation of a pland (percentage of forest in a 1km buffer around the point) factor
## how much pland to buffer 0.1°c of bioindicated temperature
whole_buffer<-mean(duets$cit_climplant_dif)/mean(duets$pland_1000m_dif)

##used for later raster computation
w_raster_eur<-focalWeight(diverse_tile,1000,type="circle")

## the computation described on the whole Europe, with forest cover stemming from Copernicus land monitoring service
## is too costly to  host or reproduce
## we load the result of the computation
## the result is quartile describing the destribution of the percentage of forest around forests cells (see fig 4.b)
## for each land systems tiles, those tiles are 9250meter in resolution
dt_lus_act<-readRDS(file.path("Data","saved_calculation","Europe_forest_cover_landsystems.RData"))
refr_lus<-fread(file.path("Data","saved_calculation","Land_systems_legend.csv"))
colnames(refr_lus)<-c("code","LUS",colnames(refr_lus)[3:ncol(refr_lus)])
refr_lus$order<-1:nrow(refr_lus)

## the raster of forest cover is large, we provide a subset of it to reproduce the figure
diverse_tile<-raster(file.path("Data","spatial_raster","M.cropland-Forest.tif"))
cropland_tile<-raster(file.path("Data","spatial_raster","Cropland.tif"))
forest_tile<-raster(file.path("Data","spatial_raster","Forest.tif"))


## get inset is a function used to create a smaller map, of one landsystem tile, and an histogram of the pland value
## open street map is also used for the background

get_inset<-function(x,y,size,rast,style_map="thunderforestlandscape",factor_pland ,title="",save_raster=FALSE){
  ext<-(extent(x,x,y,y)+size)
  
  tmp_rast<-reclassify(rast,matrix(c(0,1,2,254,255,1,NA,NA,1,1),ncol=2,nrow=5))  
  tmp_rast<-compute_pland_raster_na(tmp_rast,w_raster_eur)
  tmp_rast<-crop(tmp_rast,ext)
  tmp_dt<-data.table(xyFromCell(tmp_rast,1:ncell(tmp_rast)) )
  tmp_dt[,forest_1000m:=getValues(tmp_rast)]
  
  
  
  poly_forest<-st_union(st_as_sf(rasterToPolygons(tmp_rast)))
  poly_forest<- smoothr::smooth(poly_forest, method = "ksmooth")
  
  dummy_sf<-st_as_sf(data.frame(x=0,y=0),coords=c("x","y"),crs=st_crs(tmp_rast))[-1,]
  
  palette_cool_looking<- list(scale_fill_carto_c(palette = "ArmyRose",direction=-1,limits = c(0,100)))
  
  gg_res<-ggplot(dummy_sf)+annotation_map_tile(type=style_map,zoom =   13)+geom_tile(data=tmp_dt[!is.na(forest_1000m),],aes(x=x,y=y,fill=forest_1000m),show.legend=F)+
    theme_bw()+labs(x="",y="")+scale_x_continuous(expand=c(0,0)) +scale_y_continuous(expand=c(0,0))+theme(panel.border = element_rect(fill=NA,size=1),axis.text = element_blank(),axis.ticks = element_blank())+
    palette_cool_looking+geom_sf(data=poly_forest,color="grey10",fill=NA,size=0.5)+labs(title=title)+theme( plot.title = element_text(hjust=0.5))+ annotation_scale(location = "bl")
  
  histo_res<-ggplot(tmp_dt)+geom_histogram(aes(forest_1000m,y=stat(count)/sum(count),fill=..x..),color="black",show.legend = F,breaks=seq(from=0,to=100,length.out = 20))+theme_bw()+ theme(axis.text.y = element_blank(),axis.ticks.y = element_blank())+scale_y_continuous(labels=scales::percent)+theme(aspect.ratio = 1/3,axis.title.x = element_text(size=7))+palette_cool_looking+ylab("")+xlab("Forest cover in a 1km buffer (%)")+scale_x_continuous(sec.axis = sec_axis(~ . *factor_pland, name = "CIT cooled (C°)"))
  print(q_05(tmp_dt[,forest_1000m],na.rm=T))
  print(q_95(tmp_dt[,forest_1000m],na.rm=T))
  print("- - - -")
  print(q_10(tmp_dt[,forest_1000m],na.rm=T))
  print(q_90(tmp_dt[,forest_1000m],na.rm=T))

  return(gg_res + histo_res + plot_layout(ncol=1,nrow=2,heights = c(3,1)))
  
}

## perfome a fast computation of the percentage of forest around each forest cell
compute_pland_raster_na<-function(rast_na,w){
  rast_na<-raster::focal(rast_na,w,sum,na.rm=T,pad=T,padValue=1,NAonly=T)
  rast_na<-100-(rast_na*100)
  rast_na<-reclassify(rast_na,matrix(c(NA,100),ncol=2))
  rast_na<-reclassify(rast_na,matrix(c(0,NA),ncol=2))
  return(rast_na)
  
}
q_05<-function(x,na.rm=T)quantile(x,na.rm=T,probs=0.05)
q_95<-function(x,na.rm=T)quantile(x,na.rm=T,probs=0.95)
q_10<-function(x,na.rm=T)quantile(x,na.rm=T,probs=0.1)
q_90<-function(x,na.rm=T)quantile(x,na.rm=T,probs=0.9)

## fig 4.a
out_map_lus<-ggplot(eucoreg_poly)+theme_bw()+  annotation_map_tile(type="cartolight",zoomin = 1) + labs(x="",y="")+ 
  geom_tile(dt_lus_act,mapping=aes(x=x,y=y,fill=fill2))+scale_fill_identity(na.value="transparent")+geom_sf(fill=NA,color="grey10")+
  annotation_scale(location = "bl") +theme(panel.border = element_rect(fill = NA,size=1))  

out_map_lus<-out_map_lus+geom_tile(dt_lus_act[extrapo=="Extrapo",],mapping=aes(x=x,y=y),fill="white",alpha=0.25)


## suppl figure of the land systems
suppl_map_LUS<-ggplot(eucoreg_poly)+theme_bw()+  annotation_map_tile(type="cartolight",zoomin = 1) + labs(x="",y="")+ 
  geom_tile(dt_lus_act[!is.na(Q10),],mapping=aes(x=x,y=y,fill=LUS))+geom_sf(fill=NA,color="grey10")+
  annotation_scale(location = "bl") +theme(panel.border = element_rect(fill = NA,size=1))+
  theme(legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=6),legend.key.size =unit(6,"mm"))+
  scale_fill_manual(values=refr_lus$colorCode,breaks=refr_lus$LUS) +
  guides(fill=guide_legend(override.aes = list(color="black",shape=21),ncol=3))



## calculation of a pland (percentage of forest in a 1km buffer around the point) factor
## how much pland to buffer 0.1°c of bioindicated temperature
whole_buffer<-mean(duets$cit_climplant_dif)/mean(duets$pland_1000m_dif)

## diverse inset tile

diverse_coord<-c(3853797 ,2875000)
inset_out_diverse<-get_inset(diverse_coord[1] ,diverse_coord[2],4620*2,diverse_tile,"cartolight",factor_pland = whole_buffer,"M.cropland-Forest") 

## forest poor tile

forest_poor_coord<-c(3439403,2687465)
inset_out_no_forest<-get_inset(  forest_poor_coord[1] , forest_poor_coord[2],4620*2,cropland_tile,"cartolight",factor_pland = whole_buffer,title="Cropland")


## forest rich tile

forest_rich_coord<-c(6203094,2285467)
inset_out_rich<-get_inset(forest_rich_coord[1] ,forest_rich_coord[2]  ,4620*2,forest_tile,"cartolight",factor_pland = whole_buffer,title="Forest")


coords_3_insets<-st_as_sf(as.data.frame(do.call(rbind,list(diverse_coord,forest_poor_coord,forest_rich_coord))),coords=c(1,2),crs=st_crs(crs(forest_tile)))
coords_3_insets<-st_transform(coords_3_insets,crs=st_crs((eucoreg_poly)))
coords_3_insets$label<-c("c","a","b")
coords_3_insets$label_num<-c("3","1","2")
coords_3_insets$x<-st_coordinates(coords_3_insets)[,1]
coords_3_insets$y<-st_coordinates(coords_3_insets)[,2]


## adding the coordinates of the insets to fig4.a
out_map_lus_coords<-out_map_lus+geom_sf(data=coords_3_insets,size=1.5)+geom_text(mapping=aes(x=x,y=y+c(0,0,100000),label=label_num),data=coords_3_insets,nudge_y = 80000,fontface="bold",size=5)

## adding the legend of fig 4.a
pal_d<-fread(file.path("Data","saved_calculation","bivariate_color.csv"))

legend<-ggplot(pal_d,aes(x=x_col,y=y_col,fill=fill))+geom_tile()+scale_fill_identity()+coord_equal()+theme_void()+
  theme(panel.border = element_rect(fill=NA,color="black",size=1))+
  labs(x="Maximal CIT \ncooling by forests",y="Minimal CIT \ncooling by forests") +
  theme(axis.title = element_text(size = 9),legend.title = element_text(size=8), axis.title.y = element_text(angle = 90,hjust=0),axis.title.x = element_text(hjust=0)) +
  scale_x_discrete(breaks=c(0,1),labels=c("0 C°","-0.37 C°"))+scale_y_discrete(breaks=c(0,1),labels=c("-0.37 C°","0 C°"))+
  theme(axis.text = element_text(size=8))


scale_leg<-1.2
out_map_lus_cow<-cowplot::ggdraw() +
  cowplot::draw_plot(out_map_lus_coords)+ 
  cowplot::draw_plot(legend, 0.7, 0.52, 0.22*scale_leg, 0.22*scale_leg)


## calculating the cooling effect for supplmentary barplots and figure 4.c
dt_lus_act[,buf_min:=Q10*whole_buffer]
dt_lus_act[,buf_max:=Q90*whole_buffer]
dt_lus_act[,buf_range:=buf_max-buf_min]



make_suppl_table<-function(x,round=T){
  
  if(round) res<-x[!is.na(Q90)|!is.na(Q10),.(Mean_forest_Q90=round(mean(Q90,na.rm=T),1),
                                             Mean_forest_Q10=round( mean(Q10,na.rm=T),1),
                                             Mean_buf_max= signif(mean(buf_max,na.rm=T),2),
                                             Mean_buf_min=signif(mean(buf_min,na.rm=T),2),
                                             mean_range_buf=signif(mean(buf_range,na.rm=T),2),
                                             N_tile=.N),by=LUS] else {
                                             res<-x[!is.na(Q90)|!is.na(Q10),.(Mean_forest_Q90=mean(Q90,na.rm=T),
                                                                              Mean_forest_Q10= mean(Q10,na.rm=T),
                                                                              Mean_buf_max= mean(buf_max,na.rm=T),
                                                                              Mean_buf_min=mean(buf_min,na.rm=T),
                                                                              mean_range_buf=mean(buf_range,na.rm=T),
                                                                              N_tile=.N),by=LUS] }

                                             res<-res[order(-Mean_forest_Q90),]
}

## we compute the cooling for every Land systems tiles, and the one where we don't extrapolate frome the dataset of 2012 pairs

summary_all_tile<-make_suppl_table(dt_lus_act,F)
summary_extrapo_ok_tile<-make_suppl_table(dt_lus_act[extrapo=="OK"],F)

## get the label of landsytems
summary_extrapo_ok_tile<-merge(summary_extrapo_ok_tile,refr_lus,by="LUS")
summary_all_tile<-merge(summary_all_tile,refr_lus,by="LUS")

lus_plot<-ggplot(summary_extrapo_ok_tile[LUS!=c("Bare","Bare with few livestock"),],aes(y=reorder(LUS,-order),yend=reorder(LUS,-order),x= Mean_forest_Q10,xend=Mean_forest_Q90,color=LUS))+geom_segment(show.legend = F,size=3)+scale_color_manual(values=refr_lus$colorCode,breaks=refr_lus$LUS)+theme_bw()+theme(axis.line.y=element_blank())+geom_text(mapping=aes(x=Mean_forest_Q90+6,label=N_tile),color="black")+scale_x_continuous(limits = c(0,105),sec.axis = sec_axis(~ . *whole_buffer, name = "Community Inferred Temperature cooled (C°)"))+xlab("Landscape forest cover in a 1km buffer (%)")+ylab("Land systems")
lus_plot_extrapo<-ggplot(summary_all_tile[LUS!=c("Bare","Bare with few livestock"),],aes(y=reorder(LUS,-order),yend=reorder(LUS,-order),x= Mean_forest_Q10,xend=Mean_forest_Q90,color=LUS))+geom_segment(show.legend = F,size=3)+scale_color_manual(values=refr_lus$colorCode,breaks=refr_lus$LUS)+theme_bw()+theme(axis.line.y=element_blank())+geom_text(mapping=aes(x=Mean_forest_Q90+6,label=N_tile),color="black")+scale_x_continuous(limits = c(0,105),sec.axis = sec_axis(~ . *whole_buffer, name = "Community Inferred Temperature cooled (C°)"))+xlab("Landscape forest cover in a 1km buffer (%)")+ylab("Land systems")


tiff(file.path("results_and_figures","fig_LUS_bars_suppl.tif"),width=850,height=400,res=100)

lus_plot

dev.off()

tiff(file.path("results_and_figures","fig_LUS_bars_suppl_extrapo.tif"),width=850,height=400,res=100)

lus_plot_extrapo

dev.off()


## aggregating of the information for the main categories of Land systems
summary_agg_LUS<-summary_extrapo_ok_tile[LUS!=c("Bare","Bare with few livestock"),.(Q90=weighted.mean(Mean_forest_Q90,N_tile),Q10=weighted.mean(Mean_forest_Q10,N_tile),N_tile=sum(N_tile)),by=landUse]
summary_agg_LUS$order<-1:8
summary_agg_LUS$color<-c("#fa9039","#70a800","#c7e371","#b393c2","#d3e5e9","#c9d7c2","#89cd66","#a80000")

## fig 4.c
aggregated_lus_plot<-ggplot(summary_agg_LUS[landUse!=c("M.Grassland-Bare"),],aes(y=reorder(landUse,Q90),yend=reorder(landUse,Q90),x= Q90 ,xend=Q10,color=landUse))+geom_segment(show.legend = F,size=3)+geom_text(mapping=aes(x=Q90+6,label=N_tile),color="black")+
  scale_color_manual(values=summary_agg_LUS$color,breaks=lus$landUse)+theme_bw()+theme(axis.line.y=element_blank())+xlim(c(0,-0.4))+scale_x_continuous(limits = c(0,105),sec.axis = sec_axis(~ . *whole_buffer, name = "Community Inferred Temperature cooled (C°)"))+xlab("Landscape forest cover in a 1km buffer (%)")+ylab("Land systems")


## arrange 4.a, the 3 insets of 4.b, and 4.c
out_fig_4_maps_lus<-ggarrange(out_map_lus_cow ,
                              ggarrange(inset_out_no_forest,inset_out_rich,inset_out_diverse,ncol=3,labels = c("1","2","3"),hjust = -5),
                              aggregated_lus_plot,nrow=3,heights = c(3.75,3,1.5),
                              labels=c("(a)","(b)","(c)"),
                              font.label = list(size=16 ,color = "black", face = "bold", family = NULL))



## write the figure
reso<-2
tiff(file.path("results_and_figures","fig_4_maps_lus.tif"),width=850*reso,height=1200*reso,res=100*reso)
out_fig_4_maps_lus
dev.off()
tiff(file.path("results_and_figures","fig_LUS_maps_suppl.tif"),width=850*reso,height=750*reso,res=100*reso)
suppl_map_LUS
dev.off()



## produce table 
summary_all<-make_suppl_table(dt_lus_act)
summary_extrapo_ok<-make_suppl_table(dt_lus_act[extrapo=="OK"])

write.table(summary_extrapo_ok,"table_suppl_3.csv",row.names=F,sep=";",dec = ",")
write.table(summary_all,"table_suppl_4.csv",row.names=F,sep = ";",dec=",")


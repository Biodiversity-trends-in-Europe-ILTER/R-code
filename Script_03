########################################################################################################################################################################
# Script 03 -  Figures
# Pilotto et al. Meta-analysis of multidecadal biodiversity trends in Europe, Nature Communications
#
# This is the code to create the figures. 
# The underlying data tables were computed and exported with the scripts #01 and #02.
#
########################################################################################################################################################################

require(rgdal)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)

# Figure 1 ------------------------------------------------------------------

# To run this code it is necessary to download the Biogeoregion shapefile from:
# https://www.eea.europa.eu/data-and-maps/data/biogeographical-regions-europe-3#tab-gis-data
# and the shapefiles of marine subregions from: 
# https://www.eea.europa.eu/data-and-maps/data/msfd-regions-and-subregions


# Import shapefiles:
MAR_sub_shape <- readOGR(dsn = ".", layer = "MarineSubregions")
shapeER <- readOGR(dsn = ".", layer = "BiogeoRegions2016")

# import site coordinates:
Sites <- read.table("Sites.csv", h=T, sep =";")

# Set projection: 
crs.laea <- CRS("+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs")  # Lambert Azimuthal Equal Area
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84

coordinates(Sites) <- c("Lon", "Lat") 
proj4string(Sites) <- crs.geo 
Sites.laea <- spTransform(Sites, crs.laea)  
proj4string(shapeER) <- crs.laea 
MAR_sub_shape.laea <- spTransform(MAR_sub_shape, crs.laea)


# Figure 1a, map: 

col.list.MAP=c(brewer.pal(12,"Set3")[2], #Alpine
               "grey","grey",
               brewer.pal(12,"Set3")[4], #Atlantic
               brewer.pal(12,"Set3")[6], #BlackSea
               brewer.pal(12,"Set3")[7], #Boreal
               brewer.pal(12,"Set3")[8], #Continental
               "grey", 
               brewer.pal(12,"Set3")[10], # Mediterranean
               "grey", 
               brewer.pal(12,"Set3")[12], #Pannonian
               "grey")

tiff(paste("Fig 1a.tif"),5,5,units="in",compression="lzw",res=600) 
par(las=1,mar=c(0,1,1,1))
plot(shapeER, col=col.list.MAP, xlim=c(700000, 6500000),ylim=c(500000, 5000000))
plot(MAR_sub_shape.laea[6,], col=brewer.pal(9,"Blues")[3],add=T) # North Sea
plot(MAR_sub_shape.laea[3,], col=brewer.pal(12,"Set3")[5],add=T) #Adriatic sea
ColLEgend=c(brewer.pal(12,"Set3")[c(2,4,6,7,8,10,12)],brewer.pal(9,"Blues")[3],brewer.pal(12,"Set3")[5], "grey","darkorchid")
legend("topleft", pch=c(rep(22,times=10), 21), col=1, pt.bg=ColLEgend, c("Alpine", "Atlantic", "Black Sea","Boreal","Continental", "Mediterranean","Pannonian","North Sea","Adriatic Sea","outside", "site"), cex=0.8, pt.cex=1.2)
points(Sites.laea, pch=21, bg="magenta")
dev.off()


# Figure 1, pie charts: 
# pie charts were placed on figure 1 in a second step and colors were changed using a graphic software 

# Pie charts for Figure 1a:

TaxGr_Biog_pivotTable <- read.table("TaxGr_Biog_pivotTable.csv", h=T, sep =";")
for (i in c(3:11)){
  tiff(paste("pie_TaxGr_", names(TaxGr_Biog_pivotTable)[i], ".tif", sep=""), 6, 8, units="in", compression="lzw", res=600) 
  par(las=1, mar=c(0,0,0,0))
  pie(TaxGr_Biog_pivotTable[,i], col=as.character(TaxGr_Biog_pivotTable$colB), label="", main=names(TaxGr_Biog_pivotTable)[i])
  dev.off()}
# Legend:
tiff(paste("pie_TaxGr_Legend.tif", sep=""), 15, 6, units="in", compression="lzw", res=600) 
par(las=1,mar=c(0,0,0,0))
pie(TaxGr_Biog_pivotTable[,i], col=as.character(TaxGr_Biog_pivotTable$colB), label="", main=names(TaxGr_Biog_pivotTable)[i])
legend("topleft", pch=20, col=as.character(TaxGr_Biog_pivotTable$colB), legend=c(  "Alg","Bir","Fish", "InvA", "InvT", "Mam","Pl","Veg"), pt.cex=3, h=T)
dev.off()


# Pie charts for Figure 1b (higher panel):

TaxGr_Real_pivotTable <- read.table("TaxGr_Real_pivotTable.csv", h=T, sep =";")
for (i in c(3:5)){
  tiff(paste("pie_TaxGr_", names(TaxGr_Real_pivotTable)[i], ".tif", sep=""), 6, 8, units="in", compression="lzw", res=600) 
  par(las=1,mar=c(0,0,0,0))
  pie(TaxGr_Real_pivotTable[,i], col=as.character(TaxGr_Real_pivotTable$colB), label="", main=names(TaxGr_Real_pivotTable)[i])
  dev.off()}


# Pie charts for Figure 1b (lower panel):

Biog_Real_pivotTable <- read.table("Biog_Real_pivotTable.csv", h=T, sep =";")
for (i in c(3:5)){
  tiff(paste("pie_Biog_", names(Biog_Real_pivotTable)[i], ".tif", sep=""), 6, 8, units="in", compression="lzw", res=600) 
  par(las=1,mar=c(0,0,0,0))
  pie(Biog_Real_pivotTable[,i], col=as.character(Biog_Real_pivotTable$colB), label="", main=names(Biog_Real_pivotTable)[i])
  dev.off()}


# Pie charts for Figure 1c:

Biog_TaxGr_pivotTable <- read.table("Biog_TaxGr_pivotTable.csv", h=T, sep =";")
for (i in c(3:10)){
  tiff(paste("pie_Biog_", names(Biog_TaxGr_pivotTable)[i], ".tif", sep=""), 6, 8, units="in", compression="lzw", res=600) 
  par(las=1, mar=c(0,0,0,0))
  pie(Biog_TaxGr_pivotTable[,i], col=as.character(Biog_TaxGr_pivotTable$colB), label="", main=names(Biog_TaxGr_pivotTable)[i])
  dev.off() }
# Legend:
tiff(paste("pie_TaxGrBiog_Legend.tif", sep=""), 15, 6, units="in", compression="lzw", res=600) 
par(las=1,mar=c(0,0,0,0))
pie(Biog_TaxGr_pivotTable[,i], col=as.character(Biog_TaxGr_pivotTable$colB), label="", main=names(Biog_TaxGr_pivotTable)[i])
legend("topleft", pch=20, col=as.character(Biog_TaxGr_pivotTable$colB), legend=unique(as.character(Biog_TaxGr_pivotTable$Biogeoregion)), pt.cex=3, h=T)
dev.off()

# Fig 2 a-d ----------------------------------

# Set colors of regions based on the results of the models: 
col.list.Abund=c("black", #Alpine
                 "grey","grey",
                 "#E69F00", #Atlantic
                 "black", #BlackSea
                 "black", #Boreal
                 "black", #Continental
                 "grey", 
                 "black", # Mediterranean
                 "grey", 
                 "black", #Pannonian
                 "grey")
col.list.NTaxa=c("black", #Alpine
                 "grey","grey",
                 "black", #Atlantic
                 "#009E73", #BlackSea
                 "#009E73", #Boreal
                 "black", #Continental
                 "grey", 
                 "black", # Mediterranean
                 "grey", 
                 "black", #Pannonian
                 "grey")
col.list.Simp=c("black", #Alpine
                "grey","grey",
                "black", #Atlantic
                "#009E73", #BlackSea
                "#009E73", #Boreal
                "black", #Continental
                "grey", 
                "black", # Mediterranean
                "grey", 
                "black", #Pannonian
                "grey")
col.list.Turnover=c("#009E73", #Alpine
                    "grey","grey",
                    "black", #Atlantic
                    "#E69F00", #BlackSea
                    "#009E73", #Boreal
                    "#009E73", #Continental
                    "black", 
                    "black", # Mediterranean
                    "grey", 
                    "black", #Pannonian
                    "grey")


# create and export the figure: 
tiff("Fig 2a-d.tif", 8, 8, units="in",compression="lzw", res=600) 
layout(matrix(1:4, 2, 2, byrow=T))
par(las=1, mar=c(0,0.5,2,0.5))
plot(shapeER, col=col.list.Abund, xlim=c(700000, 6500000),ylim=c(500000, 5000000), main="a. Abundance", cex.main=1.5)
plot(MAR_sub_shape.laea[6,], col="#009E73", add=T) # North Sea
plot(MAR_sub_shape.laea[3,], col="grey30", add=T) #Adriatic sea
legend("topleft", pch=22, col=1, pt.bg=c("#E69F00", "#009E73", "black", "grey"), c("Decrease", "Increase", "No changes", "Outside"), cex=1, pt.cex=2)

plot(shapeER, col=col.list.NTaxa, xlim=c(700000, 6500000),ylim=c(500000, 5000000), main="b. Richness", cex.main=1.5)
plot(MAR_sub_shape.laea[6,], col="#009E73", add=T) # North Sea
plot(MAR_sub_shape.laea[3,], col="grey30", add=T) #Adriatic sea 

plot(shapeER, col=col.list.Simp, xlim=c(700000, 6500000),ylim=c(500000, 5000000), main="c. Diversity", cex.main=1.5)
plot(MAR_sub_shape.laea[6,], col="#009E73", add=T) # North Sea
plot(MAR_sub_shape.laea[3,], col="grey30", add=T) #Adriatic sea 

plot(shapeER, col=col.list.Turnover, xlim=c(700000, 6500000),ylim=c(500000, 5000000), main="d. Turnover", cex.main=1.5)
plot(MAR_sub_shape.laea[6,], col="#E69F00", add=T) # North Sea
plot(MAR_sub_shape.laea[3,], col="grey30", add=T) #Adriatic sea 
dev.off()


# Panel "e" of figures 2, 3 and 4 -------------------------------------------

# Figure 2e:

# Read model results: 
RMA.Abund.Biogeoregion <- read.table("RMA.Abund.Biogeoregion.csv",h=T , sep =";")
RMA.NTaxa.Biogeoregion <- read.table("RMA.NTaxa.Biogeoregion.csv",h=T , sep =";")
RMA.Simp.Biogeoregion <- read.table("RMA.Simp.Biogeoregion.csv",h=T , sep =";")
RMA.Turn.Biogeoregion <- read.table("RMA.Turn.Biogeoregion.csv",h=T , sep =";")
# create and export figure 2e:
tiff(paste("Fig 2e.tif"),15,5,units="in",compression="lzw",res=300) 
par(mar=c(3,5,1,0.2))
plot(c(RMA.Abund.Biogeoregion$estimate,NA) , ylim=c(-100,150),xlim=c(1,9.5), axes=F, ylab="Trend (S-statistic)", xlab="", pch=c(21,21,20,21,21,21,21,20,0), col="#332288", main="", cex=2.5, cex.lab=1.5)
abline(h=0,lty=1)
arrows(1:9,y0=c(RMA.Abund.Biogeoregion$estimate,NA), y1=c(RMA.Abund.Biogeoregion$ci.lb,NA),length = 0.0, angle =90,lwd=5, col="#332288", lty=c(3,3,1,3,3,3,3,1,0)) #code=2,
arrows(1:9,y0=c(RMA.Abund.Biogeoregion$estimate,NA), y1=c(RMA.Abund.Biogeoregion$ci.ub,NA),length = 0.0, angle =90,lwd=5, col="#332288", lty=c(3,3,1,3,3,3,3,1,0))
axis(2,las=2, cex.axis=1.5)
axis(1, at=c(1:9)+0.2,c("Adr","Alp","Atl","BlS","Bor","Con","Med", "NoS","Pan"), cex.axis=2)
points(c(1:9)+0.1,RMA.NTaxa.Biogeoregion$estimate, pch=c(21,21,21,20,20,21,21,20,21), col="#AA4499", cex=2.5)
abline(h=0,lty=1)
arrows(1:length(RMA.NTaxa.Biogeoregion$estimate)+0.1,y0=RMA.NTaxa.Biogeoregion$estimate, y1=RMA.NTaxa.Biogeoregion$ci.lb,length = 0.0, angle =90,lwd=5, col="#AA4499", lty=c(3,3,3,1,1,3,3,1,3)) #code=2,
arrows(1:length(RMA.NTaxa.Biogeoregion$estimate)+0.1,y0=RMA.NTaxa.Biogeoregion$estimate, y1=RMA.NTaxa.Biogeoregion$ci.ub,length = 0.0, angle =90,lwd=5, col="#AA4499", lty=c(3,3,3,1,1,3,3,1,3))
points(c(1:9)+0.2,RMA.Simp.Biogeoregion$estimate,  pch=c(21,21,21,20,20,21,21,20,21), col="#999933", cex=2.5)
abline(h=0,lty=1)
arrows(1:9+0.2,y0=RMA.Simp.Biogeoregion$estimate, y1=RMA.Simp.Biogeoregion$ci.lb,length = 0.0, angle =90,lwd=5, col="#999933", lty=c(3,3,3,1,1,3,3,1,3)) #code=2,
arrows(1:9+0.2,y0=RMA.Simp.Biogeoregion$estimate, y1=RMA.Simp.Biogeoregion$ci.ub,length = 0.0, angle =90,lwd=5, col="#999933", lty=c(3,3,3,1,1,3,3,1,3))
points(c(1:9)+0.3,RMA.Turn.Biogeoregion$estimate , ylim=c(-100,50), pch=c(21,20,21,20,20,20,21,20,21), col="#6699CC", cex=2.5)
abline(h=0,lty=1)
arrows(1:length(RMA.Turn.Biogeoregion$estimate)+0.3,y0=RMA.Turn.Biogeoregion$estimate, y1=RMA.Turn.Biogeoregion$ci.lb,length = 0.0, angle =90,lwd=5, col="#6699CC", lty=c(3,1,3,1,1,1,3,1,3)) #code=2,
arrows(1:length(RMA.Turn.Biogeoregion$estimate)+0.3,y0=RMA.Turn.Biogeoregion$estimate, y1=RMA.Turn.Biogeoregion$ci.ub,length = 0.0, angle =90,lwd=5, col="#6699CC", lty=c(3,1,3,1,1,1,3,1,3))
legend("bottomright", horiz=T, pch=20, col=c("#332288","#AA4499","#999933","#6699CC"), legend=c("Abundance", "Richness", "Diversity","Turnover"), cex=1.3, pt.cex=3)
dev.off()


# Figure 3e:

# Read model results: 
RMA.Abund.Realm <- read.table("RMA.Abund.Realm.csv", sep =";")
RMA.NTaxa.Realm <- read.table("RMA.NTaxa.Realm.csv", sep =";")
RMA.Simp.Realm <- read.table("RMA.Simp.Realm.csv", sep =";")
RMA.Turn.Realm <- read.table("RMA.Turn.Realm.csv",h=T , sep =";")
# create and export figure 3e:
tiff(paste("Fig 3e.tif"),15,5,units="in",compression="lzw",res=300) 
par(mar=c(3,5,1,0.2))
plot(RMA.Abund.Realm$estimate , ylim=c(-40,80),xlim=c(1,3.5), axes=F, ylab="Trend (S-statistic)", xlab="", pch=c(21,21,21), col="#332288", main="", cex=2.5, cex.lab=1.5)
abline(h=0,lty=1)
arrows(1:3,y0=RMA.Abund.Realm$estimate, y1=RMA.Abund.Realm$ci.lb,length = 0.0, angle =90,lwd=5, col="#332288", lty=c(3,3,3)) #code=2,
arrows(1:3,y0=RMA.Abund.Realm$estimate, y1=RMA.Abund.Realm$ci.ub,length = 0.0, angle =90,lwd=5, col="#332288", lty=c(3,3,3))
axis(2,las=2, cex.axis=1.5)
axis(1, at=c(1:3)+0.15,c("FW","MA","TE"), cex.axis=2)
points(c(1:3)+0.1,RMA.NTaxa.Realm$estimate ,pch=c(20,20,21), col="#AA4499", cex=2.5)
abline(h=0,lty=1)
arrows(1:length(RMA.NTaxa.Realm$estimate)+0.1,y0=RMA.NTaxa.Realm$estimate, y1=RMA.NTaxa.Realm$ci.lb,length = 0.0, angle =90,lwd=5, col="#AA4499", lty=c(1,1,3)) #code=2,
arrows(1:length(RMA.NTaxa.Realm$estimate)+0.1,y0=RMA.NTaxa.Realm$estimate, y1=RMA.NTaxa.Realm$ci.ub,length = 0.0, angle =90,lwd=5, col="#AA4499", lty=c(1,1,3))
points(c(1:3)+0.2,RMA.Simp.Realm$estimate, ylim=c(-100,150), pch=c(21,20,21), col="#999933", cex=2.5)
abline(h=0,lty=1)
arrows(1:3+0.2,y0=RMA.Simp.Realm$estimate, y1=RMA.Simp.Realm$ci.lb,length = 0.0, angle =90,lwd=5, col="#999933", lty=c(3,1,3)) #code=2,
arrows(1:3+0.2,y0=RMA.Simp.Realm$estimate, y1=RMA.Simp.Realm$ci.ub,length = 0.0, angle =90,lwd=5, col="#999933", lty=c(3,1,3))
points(c(1:3)+0.3,RMA.Turn.Realm$estimate, ylim=c(-100,50), pch=c(21,21,20), col="#6699CC", cex=2.5)
abline(h=0,lty=1)
arrows(1:length(RMA.Turn.Realm$estimate)+0.3,y0=RMA.Turn.Realm$estimate, y1=RMA.Turn.Realm$ci.lb,length = 0.0, angle =90,lwd=5, col="#6699CC", lty=c(3,3,1)) #code=2,
arrows(1:length(RMA.Turn.Realm$estimate)+0.3,y0=RMA.Turn.Realm$estimate, y1=RMA.Turn.Realm$ci.ub,length = 0.0, angle =90,lwd=5, col="#6699CC", lty=c(3,3,1))
legend("bottomright", horiz=T, pch=20, col=c("#332288","#AA4499","#999933","#6699CC"), legend=c("Abundance", "Richness", "Diversity","Turnover"), cex=1.3, pt.cex=3)
dev.off()


# Figure 4e:

# Read model results: 
RMA.Abund.TaxonomicGroup <- read.table("RMA.Abund.TaxonomicGroup.csv", h=T, sep =";")
RMA.NTaxa.TaxonomicGroup <- read.table("RMA.NTaxa.TaxonomicGroup.csv", h=T, sep =";")
RMA.Simp.TaxonomicGroup <- read.table("RMA.Simp.TaxonomicGroup.csv", h=T, sep =";")
RMA.Turn.TaxonomicGroup <- read.table("RMA.Turn.TaxonomicGroup.csv", h=T, sep =";")
# create and export figure 4e:
tiff(paste("Fig 4e.tif"),15,5,units="in",compression="lzw",res=300) 
par(mar=c(3,5,1,0.2))
plot(RMA.Abund.TaxonomicGroup$estimate, ylim=c(-150,140), xlim=c(1,8.5), axes=F, ylab="Trend (S-statistic)", xlab="", pch=c(21,21,21,21,20,21,21,21), col="#332288", main="", cex=2.5, cex.lab=1.5)
abline(h=0, lty=1)
arrows(1:8, y0=RMA.Abund.TaxonomicGroup$estimate, y1=RMA.Abund.TaxonomicGroup$ci.lb, length=0.0, angle=90, lwd=5, col="#332288", lty=c(3,3,3,3,1,3,3,3)) 
arrows(1:8, y0=RMA.Abund.TaxonomicGroup$estimate, y1=RMA.Abund.TaxonomicGroup$ci.ub, length=0.0, angle=90, lwd=5, col="#332288", lty=c(3,3,3,3,1,3,3,3))
axis(2, las=2, cex.axis=1.5)
axis(1, at=c(1:8)+0.2, c("Alg", "Bir", "Fish", "InvA", "InvT", "Mam", "Pl", "Pla"), cex.axis=2)
points(c(1:8)+0.1, RMA.NTaxa.TaxonomicGroup$estimate, pch=c(21,20,21,20,21,21,21,21), col="#AA4499", cex=2.5)
abline(h=0, lty=1)
arrows(1:8+0.1, y0=RMA.NTaxa.TaxonomicGroup$estimate, y1=RMA.NTaxa.TaxonomicGroup$ci.lb, length=0.0, angle=90, lwd=5, col="#AA4499", lty=c(3,1,3,1,3,3,3,3)) 
arrows(1:8+0.1, y0=RMA.NTaxa.TaxonomicGroup$estimate, y1=RMA.NTaxa.TaxonomicGroup$ci.ub, length=0.0, angle=90, lwd=5, col="#AA4499", lty=c(3,1,3,1,3,3,3,3))
points(c(1:8)+0.2, RMA.Simp.TaxonomicGroup$estimate, ylim=c(-100,150), pch=c(20,20,21,20,21,21,21,21), col="#999933", cex=2.5)
abline(h=0, lty=1)
arrows(1:8+0.2, y0=RMA.Simp.TaxonomicGroup$estimate, y1=RMA.Simp.TaxonomicGroup$ci.lb,length=0.0, angle=90,lwd=5, col="#999933", lty=c(1,1,3,1,3,3,3,3)) 
arrows(1:8+0.2, y0=RMA.Simp.TaxonomicGroup$estimate, y1=RMA.Simp.TaxonomicGroup$ci.ub,length=0.0, angle=90,lwd=5, col="#999933", lty=c(1,1,3,1,3,3,3,3))
points(c(1:8)+0.3,RMA.Turn.TaxonomicGroup$estimate, ylim=c(-100,50), pch=c(21,21,21,21,21,21,21,20), col="#6699CC", cex=2.5)
abline(h=0, lty=1)
arrows(1:8+0.3,y0=RMA.Turn.TaxonomicGroup$estimate, y1=RMA.Turn.TaxonomicGroup$ci.lb,length=0.0, angle=90, lwd=5, col="#6699CC", lty=c(3,3,3,3,3,3,3,1)) 
arrows(1:8+0.3,y0=RMA.Turn.TaxonomicGroup$estimate, y1=RMA.Turn.TaxonomicGroup$ci.ub,length=0.0, angle=90, lwd=5, col="#6699CC", lty=c(3,3,3,3,3,3,3,1))
legend("bottomright", horiz=T, pch=20, col=c("#332288","#AA4499","#999933","#6699CC"), legend=c("Abundance", "Richness", "Diversity","Turnover"), cex=1.3, pt.cex=3)
dev.off()

# Supplementary Figure 1 --------------------------------------------------

# Abundance:

Table_interaction_plot_Ab <- read.table("Table_interaction_plot_Ab.csv", h=T, sep =";")
plot_TNat_Abund <- ggplot(Table_interaction_plot_Ab, aes(x = Temp, y = pred, color = Nat) ) +
  geom_ribbon( aes(ymin = ci.lb, ymax = ci.ub, fill = Nat, color = NULL), alpha = .2) +
  geom_line( aes(y = pred), size = 1)+
  ylab("Abundance trend (stand.)")+
  xlab("Temperature trend (stand.)")+
  theme_bw()+
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


# Richness:

Table_interaction_plot_NTaxa <- read.table("Table_interaction_plot_NTaxa.csv", h=T, sep =";")
plot_TNat_NTaxa <- ggplot(Table_interaction_plot_NTaxa, aes(x = Temp, y = pred, color = Nat) ) +
  geom_ribbon( aes(ymin = ci.lb, ymax = ci.ub, fill = Nat, color = NULL), alpha = .2) +
  geom_line( aes(y = pred), size = 1)+
  ylab("Richness trend (stand.)")+
  xlab("Temperature trend (stand.)")+
  theme_bw()+
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))

# Turnover:

Table_interaction_plot_Turn <- read.table("Table_interaction_plot_Turn.csv", h=T, sep =";")
plot_TNat_Turn <- ggplot(Table_interaction_plot_Turn, aes(x = Temp, y = pred, color = Nat) ) +
  geom_ribbon( aes(ymin = ci.lb, ymax = ci.ub, fill = Nat, color = NULL), alpha = .2) +
  geom_line( aes(y = pred), size = 1)+
  ylab("Turnover trend (stand.)")+
  xlab("Temperature trend (stand.)")+
  theme_bw()+
  scale_colour_manual(values=c("#999999", "#E69F00", "#56B4E9"))+ 
  scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))


# Export figure:
tiff("Supplementary Figure 1.tif",6,12,units="in",compression="lzw",res=300)  
par(las=1,mar=c(2,0,2,1))
ggarrange( plot_TNat_Abund, 
           plot_TNat_NTaxa,
           plot_TNat_Turn,
           ncol=1,
           nrow=3)
dev.off()

# Supplementary Figure 2 --------------------------------------------------

# Read data:
Start_end_year <- read.table("Start_end_year.csv", h=T, sep =";")

mycolors <- brewer.pal(8, name = 'Dark2')
mycolors2 <- data.frame(cbind(colB=mycolors, TaxonomicGroup= unique(as.character(Start_end_year$TaxonomicGroup))))
Start_end_year <- merge(Start_end_year, mycolors2, by="TaxonomicGroup")
Start_end_year <- Start_end_year[order(Start_end_year$startYear),]

tiff("Supplementary Figure 2.tif", 4, 5, units="in", compression="lzw", res=300) 
par(mar=c(4, 4, 1, 1)) 
plot(Start_end_year$startYear,c(1:length(Start_end_year$startYear)), xlim=c(1920, 2019), type="n", las=1, axes=F, xlab="Year", ylab= "Time series")
segments(x0=Start_end_year$startYear, x1=Start_end_year$endYear, y0=c(1:length(Start_end_year$startYear)),y1=c(1:length(Start_end_year$startYear)), col=as.character(Start_end_year$colB),lwd=1.7)
axis(1)
legend("topleft", lty=1,lwd=2, col=unique(as.character(Start_end_year$colB)), 
       legend=c("Pla", "Bir", "InvA", "InvT", "Pl", "Alg","Mam","Fish"))
dev.off()

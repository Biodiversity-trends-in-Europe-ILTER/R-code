########################################################################################################################################################################
# Script 02 -  Meta-analysis
# Pilotto et al. Meta-analysis of multidecadal biodiversity trends in Europe, Nature Communications
#
# This script includes the code for the following steps:
#    (1) combine and merge the results obtained with script 01 for the 161 time series,
#    (2) synthesize the results using meta-analytical models,
#    (3) export the source data for creating the figures (to be used in Script 03),
#    (4) run sensitivity analysis.
#
########################################################################################################################################################################

# Load libraries: 

library(vegan)
library(metafor)
library(rgeos)
library(sp)  
library(reshape2)
library(RColorBrewer)
library(glmulti)

# Read and merge all results computed with Script 01  ----------------------------------------
# All 161 .csv files exported with Script 01 are in a folder titled Results01

file_list <- list.files(path = "Results01")
DATA <- data.frame() #initiate a blank data frame
for (i in 1:length(file_list)){
  new.row <- read.table(paste("Results01/", file_list[i], sep=""),h=T, sep=";", quote="\"", fill=FALSE) # read each file in the specified folder
  DATA <- rbind(DATA, new.row) # bind the new data to the building dataset
}

# Transform variables (standardize continuous variables, and log- or sqrt- transform the ones that are not normally distributed)

DATA$StudyLength <- 1+(DATA$endYear-DATA$startYear)
DATA$Lat.s <- decostand(DATA$Lat, "standardize")
DATA$StudyLength.s <- decostand(DATA$StudyLength, "standardize")
DATA$Alt.Log.s <- decostand(log(DATA$Alt+2), "standardize")
DATA$TMean_S.s <- decostand(sqrt(DATA$TMean_S+154), "standardize")
DATA$Lon.s <- decostand(sqrt(DATA$Lon+9) , "standardize")
DATA$Naturalness.s <- decostand(DATA$Naturalness, "standardize")
DATA$PTot_S.s <- decostand(DATA$PTot_S, "standardize")

# Meta-analysis -----------------------------------------------------------


# Abundance:
res.rma.Abund.TaxonomicGroup <- rma(yi = Abund_S, vi = Abund_var, method="REML", mods= ~ TaxonomicGroup-1, data= DATA)
res.rma.Abund.Realm <- rma(yi = Abund_S, vi = Abund_var, method="REML", mods= ~ Realm-1, data= DATA)
res.rma.Abund.Biogeoregion <- rma(yi = Abund_S, vi = Abund_var, method="REML", mods= ~ Biogeoregion-1, data= DATA)
# export results
RMA.Abund.TaxonomicGroup <- data.frame(estimate = res.rma.Abund.TaxonomicGroup$beta,  pval = res.rma.Abund.TaxonomicGroup$pval, zval = res.rma.Abund.TaxonomicGroup$zval, ci.lb = res.rma.Abund.TaxonomicGroup$ci.lb, ci.ub = res.rma.Abund.TaxonomicGroup$ci.ub)
RMA.Abund.Realm <- data.frame(estimate = res.rma.Abund.Realm$beta,  pval = res.rma.Abund.Realm$pval, zval = res.rma.Abund.Realm$zval, ci.lb = res.rma.Abund.Realm$ci.lb, ci.ub = res.rma.Abund.Realm$ci.ub)
RMA.Abund.Biogeoregion <- data.frame(estimate = res.rma.Abund.Biogeoregion$beta,  pval = res.rma.Abund.Biogeoregion$pval, zval = res.rma.Abund.Biogeoregion$zval, ci.lb = res.rma.Abund.Biogeoregion$ci.lb, ci.ub = res.rma.Abund.Biogeoregion$ci.ub)
write.table(RMA.Abund.TaxonomicGroup, "RMA.Abund.TaxonomicGroup.csv", sep =";")
write.table(RMA.Abund.Realm, "RMA.Abund.Realm.csv", sep =";")
write.table(RMA.Abund.Biogeoregion, "RMA.Abund.Biogeoregion.csv", sep =";")

# Taxonomic richness:
res.rma.NTaxa.TaxonomicGroup <- rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods= ~ TaxonomicGroup-1, data= DATA)
res.rma.NTaxa.Realm <- rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods= ~ Realm-1, data= DATA)
res.rma.NTaxa.Biogeoregion <- rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods= ~ Biogeoregion-1, data= DATA)
# export results
RMA.NTaxa.TaxonomicGroup <- data.frame(estimate = res.rma.NTaxa.TaxonomicGroup$beta,  pval = res.rma.NTaxa.TaxonomicGroup$pval, zval = res.rma.NTaxa.TaxonomicGroup$zval, ci.lb = res.rma.NTaxa.TaxonomicGroup$ci.lb, ci.ub = res.rma.NTaxa.TaxonomicGroup$ci.ub)
RMA.NTaxa.Realm <- data.frame(estimate = res.rma.NTaxa.Realm$beta,  pval = res.rma.NTaxa.Realm$pval, zval = res.rma.NTaxa.Realm$zval, ci.lb = res.rma.NTaxa.Realm$ci.lb, ci.ub = res.rma.NTaxa.Realm$ci.ub)
RMA.NTaxa.Biogeoregion <- data.frame(estimate = res.rma.NTaxa.Biogeoregion$beta,  pval = res.rma.NTaxa.Biogeoregion$pval, zval = res.rma.NTaxa.Biogeoregion$zval, ci.lb = res.rma.NTaxa.Biogeoregion$ci.lb, ci.ub = res.rma.NTaxa.Biogeoregion$ci.ub)
write.table(RMA.NTaxa.TaxonomicGroup, "RMA.NTaxa.TaxonomicGroup.csv", sep =";")
write.table(RMA.NTaxa.Realm, "RMA.NTaxa.Realm.csv", sep =";")
write.table(RMA.NTaxa.Biogeoregion, "RMA.NTaxa.Biogeoregion.csv", sep =";")

# Simpson´s diversity:
res.rma.Simp.TaxonomicGroup <- rma(yi = Simp_S, vi = Simp_var, method="REML", mods= ~ TaxonomicGroup-1, data= DATA)
res.rma.Simp.Realm <- rma(yi = Simp_S, vi = Simp_var, method="REML", mods= ~ Realm-1, data= DATA)
res.rma.Simp.Biogeoregion <- rma(yi = Simp_S, vi = Simp_var, method="REML", mods= ~ Biogeoregion-1, data= DATA)
# export results
RMA.Simp.TaxonomicGroup <- data.frame(estimate = res.rma.Simp.TaxonomicGroup$beta,  pval = res.rma.Simp.TaxonomicGroup$pval, zval = res.rma.Simp.TaxonomicGroup$zval, ci.lb = res.rma.Simp.TaxonomicGroup$ci.lb, ci.ub = res.rma.Simp.TaxonomicGroup$ci.ub)
RMA.Simp.Realm <- data.frame(estimate = res.rma.Simp.Realm$beta,  pval = res.rma.Simp.Realm$pval, zval = res.rma.Simp.Realm$zval, ci.lb = res.rma.Simp.Realm$ci.lb, ci.ub = res.rma.Simp.Realm$ci.ub)
RMA.Simp.Biogeoregion <- data.frame(estimate = res.rma.Simp.Biogeoregion$beta,  pval = res.rma.Simp.Biogeoregion$pval, zval = res.rma.Simp.Biogeoregion$zval, ci.lb = res.rma.Simp.Biogeoregion$ci.lb, ci.ub = res.rma.Simp.Biogeoregion$ci.ub)
write.table(RMA.Simp.TaxonomicGroup, "RMA.Simp.TaxonomicGroup.csv", sep =";")
write.table(RMA.Simp.Realm, "RMA.Simp.Realm.csv", sep =";")
write.table(RMA.Simp.Biogeoregion, "RMA.Simp.Biogeoregion.csv", sep =";")

# Turnover:
res.rma.Turn.TaxonomicGroup <- rma(yi = Turn_S, vi = Turn_var, method="REML", mods= ~ TaxonomicGroup-1, data= DATA)
res.rma.Turn.Realm <- rma(yi = Turn_S, vi = Turn_var, method="REML", mods= ~ Realm-1, data= DATA)
res.rma.Turn.Biogeoregion <- rma(yi = Turn_S, vi = Turn_var, method="REML", mods= ~ Biogeoregion-1, data= DATA)
# export results
RMA.Turn.TaxonomicGroup <- data.frame(estimate = res.rma.Turn.TaxonomicGroup$beta,  pval = res.rma.Turn.TaxonomicGroup$pval, zval = res.rma.Turn.TaxonomicGroup$zval, ci.lb = res.rma.Turn.TaxonomicGroup$ci.lb, ci.ub = res.rma.Turn.TaxonomicGroup$ci.ub)
RMA.Turn.Realm <- data.frame(estimate = res.rma.Turn.Realm$beta,  pval = res.rma.Turn.Realm$pval, zval = res.rma.Turn.Realm$zval, ci.lb = res.rma.Turn.Realm$ci.lb, ci.ub = res.rma.Turn.Realm$ci.ub)
RMA.Turn.Biogeoregion <- data.frame(estimate = res.rma.Turn.Biogeoregion$beta,  pval = res.rma.Turn.Biogeoregion$pval, zval = res.rma.Turn.Biogeoregion$zval, ci.lb = res.rma.Turn.Biogeoregion$ci.lb, ci.ub = res.rma.Turn.Biogeoregion$ci.ub)
write.table(RMA.Turn.TaxonomicGroup, "RMA.Turn.TaxonomicGroup.csv", sep =";")
write.table(RMA.Turn.Realm, "RMA.Turn.Realm.csv", sep =";")
write.table(RMA.Turn.Biogeoregion, "RMA.Turn.Biogeoregion.csv", sep =";")


# Dataset descriptors -----------------------------------------------------
# This section allows to export data tables that will be used to create the figures: 

# Site coordinates:
Sites <- DATA[,c("Lat","Lon")]
write.table(Sites, "Sites.csv", sep =";")

# Study period:
Start_end_year = DATA[, c("TaxonomicGroup", "startYear", "endYear")]
write.table(Start_end_year, "Start_end_year.csv", sep =";")

# Data tables for pie charts of Figure 1: 
mycolors <- brewer.pal(8, name = 'Dark2')
mycolors2 <- data.frame(cbind(colB=mycolors, TaxonomicGroup= unique(as.character(DATA$TaxonomicGroup))))
TaxGr_Biog_pivotTable <- dcast(DATA, TaxonomicGroup ~ Biogeoregion, value.var="Realm", fun = length)
row.names(TaxGr_Biog_pivotTable) <- TaxGr_Biog_pivotTable$TaxonomicGroup
TaxGr_Biog_pivotTable <- merge(mycolors2, TaxGr_Biog_pivotTable ,by="TaxonomicGroup")
write.table(TaxGr_Biog_pivotTable, "TaxGr_Biog_pivotTable.csv", sep =";")

mycolors2 <-  data.frame(cbind(colB=mycolors, TaxonomicGroup=unique(as.character(DATA$TaxonomicGroup))))
TaxGr_Real_pivotTable <- dcast(DATA,TaxonomicGroup ~ Realm, value.var="Realm", fun = length)
row.names(TaxGr_Real_pivotTable) <- TaxGr_Real_pivotTable$TaxonomicGroup
TaxGr_Real_pivotTable <- merge(mycolors2, TaxGr_Real_pivotTable, by="TaxonomicGroup")
write.table(TaxGr_Real_pivotTable, "TaxGr_Real_pivotTable.csv", sep =";")

mycolors3 <-  data.frame(cbind(colB=mycolors, Biogeoregion=unique(as.character(DATA$Biogeoregion))))
Biog_Real_pivotTable <- dcast(DATA, Biogeoregion ~ Realm, value.var="Biogeoregion", fun = length)
row.names(Biog_Real_pivotTable) <- Biog_Real_pivotTable$Biogeoregion
Biog_Real_pivotTable <- merge(mycolors3, Biog_Real_pivotTable, by="Biogeoregion")
write.table(Biog_Real_pivotTable, "Biog_Real_pivotTable.csv", sep =";")

mycolors3= data.frame(cbind(colB=mycolors,Biogeoregion= unique(as.character(DATA$Biogeoregion))))
Biog_TaxGr_pivotTable=dcast(DATA,Biogeoregion~TaxonomicGroup, value.var="Biogeoregion", fun = length)
row.names(Biog_TaxGr_pivotTable)=Biog_TaxGr_pivotTable$Biogeoregion
Biog_TaxGr_pivotTable=merge(mycolors3,Biog_TaxGr_pivotTable,by="Biogeoregion")
write.table(Biog_TaxGr_pivotTable, "Biog_TaxGr_pivotTable.csv", sep =";")


# Information-theoretic approach for meta-analysis mixed model selection and multi-model inference ------------------------------------------------------------------

# This code is based on these two sources: https://github.com/nlkinlock/LDGmeta-analysis/blob/master/4-analysis.R
# and: http://www.metafor-project.org/doku.php/tips:model_selection_with_glmulti

# define function:
rma.glmulti <- function(formula, data, ...){
  rma(formula, vi, data=data, method="ML", ...)}
setOldClass("rma.uni")
setMethod('getfit', 'rma.uni', function(object, ...) {
  if (object$test=="z") {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=Inf)
  } else {
    cbind(estimate=coef(object), se=sqrt(diag(vcov(object))), df=object$k-object$p)
  }
})

# Check correlations among explanatory variables: 
cor(DATA[,c("TMean_S.s","PTot_S.s","Naturalness.s","Alt.Log.s","Lat.s","Lon.s")], use="pairwise.complete.obs")

#Remove timeseries with NAs in the target variables:
DATA.Abund.m <- DATA[!apply(DATA[,c("Abund_S","Abund_var","TMean_S.s","PTot_S.s","Naturalness.s","Lat.s","Lon.s","Alt.Log.s")], 1, anyNA),]
DATA.NTaxa.m <- DATA[!apply(DATA[,c("NTaxa_S","NTaxa_var","TMean_S.s","PTot_S.s","Naturalness.s","Lat.s","Lon.s","Alt.Log.s")], 1, anyNA),]
DATA_Simp.m <- DATA[!apply(DATA[,c("Simp_S","Simp_var","TMean_S.s","PTot_S.s","Naturalness.s","Lat.s","Lon.s","Alt.Log.s")], 1, anyNA),]
DATA.Turn.m <- DATA[!apply(DATA[,c("Turn_S","Turn_var","TMean_S.s","PTot_S.s","Naturalness.s","Lat.s","Lon.s","Alt.Log.s")], 1, anyNA),]

# Abundance:

vi <- DATA.Abund.m$Abund_var
modelsFULL.Abund <- glmulti(Abund_S ~ TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s,
                            data = DATA.Abund.m, level = 1,
                            fitfunction = rma.glmulti, crit = "aicc", confsetsize = 500)
rank.modelsFULL.Abund <- weightable(modelsFULL.Abund)
rank.modelsFULL.Abund <- rank.modelsFULL.Abund[rank.modelsFULL.Abund$aicc <= min(rank.modelsFULL.Abund$aicc) + 2, ]

# Model-averaged coefficients (estimate and 95% C.I.) and the relative importance (sum of Akaike weights, see main text for explanation) of each explanatory 
weights <- round(coef.glmulti(modelsFULL.Abund, select = nrow(rank.modelsFULL.Abund), icmethod = "Burnham"), 7)
weights  # Supplementary table 2

# average variation explained:
R2.Abund.sel <- numeric()
for (i in 1:nrow(rank.modelsFULL.Abund)) {
  R2.Abund.sel <- append(R2.Abund.sel,modelsFULL.Abund@objects[[i]]$R2)}
Averaged.R2.Abund <- sum(R2.Abund.sel*rank.modelsFULL.Abund[,3])/sum(rank.modelsFULL.Abund[,3])


# Richness:

vi <- DATA.NTaxa.m$NTaxa_var
modelsFULL.NTaxa <- glmulti(NTaxa_S ~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s,
                            data = DATA.NTaxa.m, level = 1,
                            fitfunction = rma.glmulti, crit = "aicc", confsetsize = 500)
rank.modelsFULL.NTaxa <- weightable(modelsFULL.NTaxa)
rank.modelsFULL.NTaxa <- rank.modelsFULL.NTaxa[rank.modelsFULL.NTaxa$aicc <= min(rank.modelsFULL.NTaxa$aicc) + 2, ]

# Model-averaged coefficients (estimate and 95% C.I.) and the relative importance (sum of Akaike weights, see main text for explanation) of each explanatory 
weights <- round(coef.glmulti(modelsFULL.NTaxa, select = nrow(rank.modelsFULL.NTaxa), icmethod = "Burnham"), 7)
weights # Supplementary table 2

# average variation explained:
R2.NTaxa.sel <- numeric()
for (i in 1:nrow(rank.modelsFULL.NTaxa)) {
  R2.NTaxa.sel <- append(R2.NTaxa.sel,modelsFULL.NTaxa@objects[[i]]$R2)}
Averaged.R2.NTaxa <- sum(R2.NTaxa.sel*rank.modelsFULL.NTaxa[,3])/sum(rank.modelsFULL.NTaxa[,3])


# Simpson´s diversity:

vi <- DATA_Simp.m$Simp_var
modelsFULL_Simp <- glmulti(Simp_S ~ TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s,
                           data = DATA_Simp.m, level = 1,
                           fitfunction = rma.glmulti, crit = "aicc", confsetsize = 500)
rank.modelsFULL_Simp <- weightable(modelsFULL_Simp)
rank.modelsFULL_Simp <- rank.modelsFULL_Simp[rank.modelsFULL_Simp$aicc <= min(rank.modelsFULL_Simp$aicc) + 2, ]
rank.modelsFULL_Simp

# Model-averaged coefficients (estimate and 95% C.I.) and the relative importance (sum of Akaike weights, see main text for explanation) of each explanatory 
weights <- round(coef.glmulti(modelsFULL_Simp, select = nrow(rank.modelsFULL_Simp), icmethod = "Burnham"), 7)
weights # Supplementary table 2

# average variation explained:
R2_Simp.sel <- numeric()
for (i in 1:nrow(rank.modelsFULL_Simp)) {
  R2_Simp.sel <- append(R2_Simp.sel,modelsFULL_Simp@objects[[i]]$R2)}
Averaged.R2_Simp <- sum(R2_Simp.sel*rank.modelsFULL_Simp[,3])/sum(rank.modelsFULL_Simp[,3])


# Turnover:

vi=DATA.Turn.m$Turn_var
modelsFULL.Turn <- glmulti(Turn_S ~ TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s,
                               data = DATA.Turn.m, level = 1,
                               fitfunction = rma.glmulti, crit = "aicc", confsetsize = 500)
rank.modelsFULL.Turn <- weightable(modelsFULL.Turn)
rank.modelsFULL.Turn <- rank.modelsFULL.Turn[rank.modelsFULL.Turn$aicc <= min(rank.modelsFULL.Turn$aicc) + 2, ]

# Model-averaged coefficients (estimate and 95% C.I.) and the relative importance (sum of Akaike weights, see main text for explanation) of each explanatory 
weights <- round(coef.glmulti(modelsFULL.Turn, select = nrow(rank.modelsFULL.Turn), icmethod = "Burnham"), 7)
weights  # Supplementary table 2

# average variation explained:
R2.Turn.sel <- numeric()
for (i in 1:nrow(rank.modelsFULL.Turn)) {
  R2.Turn.sel <- append(R2.Turn.sel,modelsFULL.Turn@objects[[i]]$R2)}
Averaged.R2.Turn <- sum(R2.Turn.sel*rank.modelsFULL.Turn[,3])/sum(rank.modelsFULL.Turn[,3])



# Test interactions between site Naturalness and Temperature and Precipitation--------------------

# Abundance:

# Fit the model as resulting from the model selection above
res.rma.Abund.Continuous0 <- rma(yi = Abund_S, vi = Abund_var, method="ML", 
                               mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+StudyLength.s, data= DATA.Abund.m)
# Add interaction between site naturalness and temperature and precipitation trends:
res.rma.Abund.Continuous_1 <- rma(yi = Abund_S, vi = Abund_var, method="ML", 
                                mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+StudyLength.s+TMean_S.s:Naturalness.s+PTot_S.s:Naturalness.s, data= DATA.Abund.m)
# Add interaction between site naturalness and temperature trend:
res.rma.Abund.Continuous_2 <- rma(yi = Abund_S, vi = Abund_var, method="ML", 
                                mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+StudyLength.s+TMean_S.s:Naturalness.s, data= DATA.Abund.m)
# Add interaction between site naturalness and precipitation trends:
res.rma.Abund.Continuous_3 <- rma(yi = Abund_S, vi = Abund_var, method="ML", 
                                mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+StudyLength.s+PTot_S.s:Naturalness.s, data= DATA.Abund.m)
# Compare the models and select the one with lowest AIC: 
AIC(res.rma.Abund.Continuous0, res.rma.Abund.Continuous_1, res.rma.Abund.Continuous_2, res.rma.Abund.Continuous_3)
res.rma.Abund.Continuous_selected <- rma(yi = Abund_S, vi = Abund_var, method="REML", 
                                       mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+StudyLength.s+TMean_S.s:Naturalness.s, data= DATA.Abund.m)


# Richness: 

res.rma.NTaxa.Continuous0<-rma(yi = NTaxa_S, vi = NTaxa_var, method="ML", 
                               mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s, data= DATA.NTaxa.m)
res.rma.NTaxa.Continuous_1<-rma(yi = NTaxa_S, vi = NTaxa_var, method="ML", 
                                mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+TMean_S.s:Naturalness.s+PTot_S.s:Naturalness.s, data= DATA.NTaxa.m)
res.rma.NTaxa.Continuous_2<-rma(yi = NTaxa_S, vi = NTaxa_var, method="ML", 
                                mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+TMean_S.s:Naturalness.s, data= DATA.NTaxa.m)
res.rma.NTaxa.Continuous_3<-rma(yi = NTaxa_S, vi = NTaxa_var, method="ML", 
                                mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+PTot_S.s:Naturalness.s, data= DATA.NTaxa.m)
# Compare the models and select the one with lowest AIC: 
AIC(res.rma.NTaxa.Continuous0, res.rma.NTaxa.Continuous_1, res.rma.NTaxa.Continuous_2, res.rma.NTaxa.Continuous_3)
res.rma.NTaxa.Continuous_selected<-rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", 
                                       mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+TMean_S.s:Naturalness.s, data= DATA.NTaxa.m)


# Simpson´s diversity:
# Naturalness is not in the selected model. 


# Turnover: 

res.rma.Turn.Continuous0<-rma(yi = Turn_S, vi = Turn_var, method="ML", 
                                  mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s, data= DATA.Turn.m)
res.rma.Turn.Continuous_1<-rma(yi = Turn_S, vi = Turn_var, method="ML", 
                                   mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+TMean_S.s:Naturalness.s+PTot_S.s:Naturalness.s, data= DATA.Turn.m)
res.rma.Turn.Continuous_2<-rma(yi = Turn_S, vi = Turn_var, method="ML", 
                                   mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+TMean_S.s:Naturalness.s, data= DATA.Turn.m)
res.rma.Turn.Continuous_3<-rma(yi = Turn_S, vi = Turn_var, method="ML", 
                                   mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+PTot_S.s:Naturalness.s, data= DATA.Turn.m)

AIC(res.rma.Turn.Continuous0, res.rma.Turn.Continuous_1, res.rma.Turn.Continuous_2, res.rma.Turn.Continuous_3)
res.rma.Turn.Continuous_selected<-rma(yi = Turn_S, vi = Turn_var, method="REML", 
                                          mods=~TMean_S.s+PTot_S.s+Naturalness.s+Lat.s+Lon.s+Alt.Log.s+StudyLength.s+TMean_S.s:Naturalness.s, data= DATA.Turn.m)



# Prepare and export data tables for Supplementary figure 1: 
# predict the values of the response variables using the selected models above, with minimum, mean and maximum level of naturalness, and the full range of temperature S-statistics, 
# the other variables are set to their mean values. Note that all variables are standardized.

# Abundance:

res.rma.Abund.Continuous_selected_minNat <- data.frame(predict(res.rma.Abund.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.Abund.m$PTot_S.s),-1.9755,  mean(DATA.Abund.m$Lat.s),
                                                                                                     mean(DATA.Abund.m$StudyLength.s), -3.0897:4.1233* -1.9755)))
res.rma.Abund.Continuous_selected_meanNat <- data.frame(predict(res.rma.Abund.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.Abund.m$PTot_S.s),0,  mean(DATA.Abund.m$Lat.s),
                                                                                                      mean(DATA.Abund.m$StudyLength.s), -3.0897:4.1233* 0)))
res.rma.Abund.Continuous_selected_maxNat <- data.frame(predict(res.rma.Abund.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.Abund.m$PTot_S.s),0.9877 ,  mean(DATA.Abund.m$Lat.s),
                                                                                                     mean(DATA.Abund.m$StudyLength.s), -3.0897:4.1233* 0.9877 )))
res.rma.Abund.Continuous_selected_minNat$Temp <- c(-3.0897:4.1233)
res.rma.Abund.Continuous_selected_minNat$Nat <- "min"
res.rma.Abund.Continuous_selected_meanNat$Temp <- c(-3.0897:4.1233)
res.rma.Abund.Continuous_selected_meanNat$Nat <-"mean"
res.rma.Abund.Continuous_selected_maxNat$Temp <- c(-3.0897:4.1233)
res.rma.Abund.Continuous_selected_maxNat$Nat <- "max"

Table_interaction_plot_Ab <- rbind(data.frame(res.rma.Abund.Continuous_selected_minNat), data.frame(res.rma.Abund.Continuous_selected_meanNat), data.frame(res.rma.Abund.Continuous_selected_maxNat))
write.table(Table_interaction_plot_Ab, "Table_interaction_plot_Ab.csv", sep =";")


# Richness:

res.rma.NTaxa.Continuous_selected_minNat <- data.frame(predict(res.rma.NTaxa.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.NTaxa.m$PTot_S.s),-1.9755,  mean(DATA.NTaxa.m$Lat.s),  mean(DATA.NTaxa.m$Lon.s),  mean(DATA.NTaxa.m$Alt.Log.s),
                                                                                                     mean(DATA.NTaxa.m$StudyLength.s), -3.0897:4.1233* -1.9755)))
res.rma.NTaxa.Continuous_selected_meanNat <- data.frame(predict(res.rma.NTaxa.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.NTaxa.m$PTot_S.s),0,  mean(DATA.NTaxa.m$Lat.s),  mean(DATA.NTaxa.m$Lon.s),  mean(DATA.NTaxa.m$Alt.Log.s),
                                                                                                      mean(DATA.NTaxa.m$StudyLength.s), -3.0897:4.1233* 0)))
res.rma.NTaxa.Continuous_selected_maxNat <- data.frame(predict(res.rma.NTaxa.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.NTaxa.m$PTot_S.s),0.9877 ,  mean(DATA.NTaxa.m$Lat.s),  mean(DATA.NTaxa.m$Lon.s),  mean(DATA.NTaxa.m$Alt.Log.s),
                                                                                                     mean(DATA.NTaxa.m$StudyLength.s), -3.0897:4.1233* 0.9877 )))
res.rma.NTaxa.Continuous_selected_minNat$Temp <- c(-3.0897:4.1233)
res.rma.NTaxa.Continuous_selected_minNat$Nat <- "min"
res.rma.NTaxa.Continuous_selected_meanNat$Temp <- c(-3.0897:4.1233)
res.rma.NTaxa.Continuous_selected_meanNat$Nat <-"mean"
res.rma.NTaxa.Continuous_selected_maxNat$Temp <- c(-3.0897:4.1233)
res.rma.NTaxa.Continuous_selected_maxNat$Nat <- "max"

Table_interaction_plot_NTaxa <- rbind(data.frame(res.rma.NTaxa.Continuous_selected_minNat), data.frame(res.rma.NTaxa.Continuous_selected_meanNat), data.frame(res.rma.NTaxa.Continuous_selected_maxNat))
write.table(Table_interaction_plot_NTaxa, "Table_interaction_plot_NTaxa.csv", sep =";")


# Turnover:

res.rma.Turn.Continuous_selected_minNat <- data.frame(predict(res.rma.Turn.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.Turn.m$PTot_S.s),-1.9755,  mean(DATA.Turn.m$Lat.s),  mean(DATA.Turn.m$Lon.s),  mean(DATA.Turn.m$Alt.Log.s),
                                                                                                           mean(DATA.Turn.m$StudyLength.s), -3.0897:4.1233* -1.9755)))
res.rma.Turn.Continuous_selected_meanNat <- data.frame(predict(res.rma.Turn.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.Turn.m$PTot_S.s),0,  mean(DATA.Turn.m$Lat.s),  mean(DATA.Turn.m$Lon.s),  mean(DATA.Turn.m$Alt.Log.s),
                                                                                                            mean(DATA.Turn.m$StudyLength.s), -3.0897:4.1233* 0)))
res.rma.Turn.Continuous_selected_maxNat <- data.frame(predict(res.rma.Turn.Continuous_selected, newmods=cbind(-3.0897:4.1233, mean(DATA.Turn.m$PTot_S.s),0.9877 ,  mean(DATA.Turn.m$Lat.s),  mean(DATA.Turn.m$Lon.s),  mean(DATA.Turn.m$Alt.Log.s),
                                                                                                           mean(DATA.Turn.m$StudyLength.s), -3.0897:4.1233* 0.9877 )))
res.rma.Turn.Continuous_selected_minNat$Temp <- c(-3.0897:4.1233)
res.rma.Turn.Continuous_selected_minNat$Nat <- "min"
res.rma.Turn.Continuous_selected_meanNat$Temp <- c(-3.0897:4.1233)
res.rma.Turn.Continuous_selected_meanNat$Nat <-"mean"
res.rma.Turn.Continuous_selected_maxNat$Temp <- c(-3.0897:4.1233)
res.rma.Turn.Continuous_selected_maxNat$Nat <- "max"

Table_interaction_plot_Turn <- rbind(data.frame(res.rma.Turn.Continuous_selected_minNat), data.frame(res.rma.Turn.Continuous_selected_meanNat), data.frame(res.rma.Turn.Continuous_selected_maxNat))
write.table(Table_interaction_plot_Turn, "Table_interaction_plot_Turn.csv", sep =";")


# Sensitivity analysis ----------------------------------------------------

# Randomly select datasets for overepresented groups and re-run models (repeat 5 times): 

#Subset only 8 Aquatic invertebrates from Boreal region:
set.seed(100)
DATA_Boreal_AqInv_Sel1 <- sample(DATA[DATA$Biogeoregion=="Boreal" & DATA$TaxonomicGroup == "InvertebratesA","TimeSeries"], size=8, replace=F)

set.seed(111)
DATA_Boreal_AqInv_Sel2 <- sample(DATA[DATA$Biogeoregion=="Boreal" & DATA$TaxonomicGroup == "InvertebratesA","TimeSeries"], size=8, replace=F)

set.seed(122)
DATA_Boreal_AqInv_Sel3 <- sample(DATA[DATA$Biogeoregion=="Boreal" & DATA$TaxonomicGroup == "InvertebratesA","TimeSeries"], size=8, replace=F)

set.seed(133)
DATA_Boreal_AqInv_Sel4 <- sample(DATA[DATA$Biogeoregion=="Boreal" & DATA$TaxonomicGroup == "InvertebratesA","TimeSeries"], size=8, replace=F)

set.seed(144)
DATA_Boreal_AqInv_Sel5 <- sample(DATA[DATA$Biogeoregion=="Boreal" & DATA$TaxonomicGroup == "InvertebratesA","TimeSeries"], size=8, replace=F)



#Subset only 14 Terrestrial invertebrates from Atlantic region:
set.seed(100)
DATA_Atlantic_TerrInv_Sel1 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=14, replace=F)

set.seed(111)
DATA_Atlantic_TerrInv_Sel2 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=14, replace=F)

set.seed(122)
DATA_Atlantic_TerrInv_Sel3 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=14, replace=F)

set.seed(133)
DATA_Atlantic_TerrInv_Sel4 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=14, replace=F)

set.seed(144)
DATA_Atlantic_TerrInv_Sel5 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=14, replace=F)


#Subset only 5 plants from Alpine region:
set.seed(100)
DATA_Alpine_Vegetation_Sel1 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=5, replace=F)

set.seed(111)
DATA_Alpine_Vegetation_Sel2 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=5, replace=F)

set.seed(122)
DATA_Alpine_Vegetation_Sel3 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=5, replace=F)

set.seed(133)
DATA_Alpine_Vegetation_Sel4 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=5, replace=F)

set.seed(144)
DATA_Alpine_Vegetation_Sel5 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=5, replace=F)


# Combine original DATA datasets  with only randomly  selected data for the over-represented groups:

DATA_Boreal_AqInv <- DATA[DATA$Biogeoregion=="Boreal" & DATA$TaxonomicGroup == "InvertebratesA",]
DATA_without_Boreal_AqInv <- DATA[!DATA$TimeSeries %in% DATA_Boreal_AqInv$TimeSeries,]
DATA_Boreal_AqInv_Sel_1 <- rbind(DATA_without_Boreal_AqInv,DATA[DATA$TimeSeries %in% DATA_Boreal_AqInv_Sel1,] )
DATA_Boreal_AqInv_Sel_2 <- rbind(DATA_without_Boreal_AqInv,DATA[DATA$TimeSeries %in% DATA_Boreal_AqInv_Sel2,] )
DATA_Boreal_AqInv_Sel_3 <- rbind(DATA_without_Boreal_AqInv,DATA[DATA$TimeSeries %in% DATA_Boreal_AqInv_Sel3,] )
DATA_Boreal_AqInv_Sel_4 <- rbind(DATA_without_Boreal_AqInv,DATA[DATA$TimeSeries %in% DATA_Boreal_AqInv_Sel4,] )
DATA_Boreal_AqInv_Sel_5 <- rbind(DATA_without_Boreal_AqInv,DATA[DATA$TimeSeries %in% DATA_Boreal_AqInv_Sel5,] )

DATA_Atlantic_TerrInv <- DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT",]
DATA_without_Atlantic_TerrInv <- DATA[!DATA$TimeSeries %in% DATA_Atlantic_TerrInv$TimeSeries,]
DATA_Atlantic_TerrInv_Sel_1 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel1,] )
DATA_Atlantic_TerrInv_Sel_2 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel2,] )
DATA_Atlantic_TerrInv_Sel_3 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel3,] )
DATA_Atlantic_TerrInv_Sel_4 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel4,] )
DATA_Atlantic_TerrInv_Sel_5 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel5,] )

DATA_Alpine_Vegetation <- DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation",]
DATA_without_Alpine_Vegetation <- DATA[!DATA$TimeSeries %in% DATA_Alpine_Vegetation$TimeSeries,]
DATA_Alpine_Vegetation_Sel_1 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel1,] )
DATA_Alpine_Vegetation_Sel_2 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel2,] )
DATA_Alpine_Vegetation_Sel_3 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel3,] )
DATA_Alpine_Vegetation_Sel_4 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel4,] )
DATA_Alpine_Vegetation_Sel_5 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel5,] )


# Re-run models:
# Run the models below 3 times, after selecting the data referring to (1) Boreal_AqInv, (2) Alpine_Vegetation and (3) Atlantic_TerrInv:
DATASET_SA1 <- DATA_Alpine_Vegetation_Sel_1  #  DATA_Atlantic_TerrInv_Sel_1 # DATA_Boreal_AqInv_Sel_1
DATASET_SA2 <- DATA_Alpine_Vegetation_Sel_2  #  DATA_Atlantic_TerrInv_Sel_2 # DATA_Boreal_AqInv_Sel_2
DATASET_SA3 <- DATA_Alpine_Vegetation_Sel_3  #  DATA_Atlantic_TerrInv_Sel_3 # DATA_Boreal_AqInv_Sel_3
DATASET_SA4 <- DATA_Alpine_Vegetation_Sel_4  #  DATA_Atlantic_TerrInv_Sel_4 # DATA_Boreal_AqInv_Sel_4
DATASET_SA5 <- DATA_Alpine_Vegetation_Sel_5  #  DATA_Atlantic_TerrInv_Sel_5 # DATA_Boreal_AqInv_Sel_5


# Abundance:
res.rma.Abund.Biota_SA1=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA1)
res.rma.Abund.Biota_SA2=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA2)
res.rma.Abund.Biota_SA3=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA3)
res.rma.Abund.Biota_SA4=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA4)
res.rma.Abund.Biota_SA5=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA5)

# Richness:
res.rma.NTaxa.Biota_SA1=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA1)
res.rma.NTaxa.Biota_SA2=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA2)
res.rma.NTaxa.Biota_SA3=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA3)
res.rma.NTaxa.Biota_SA4=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA4)
res.rma.NTaxa.Biota_SA5=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA5)

# Simpson´s diversity:
res.rma_Simp.Biota_SA1=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA1)
res.rma_Simp.Biota_SA2=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA2)
res.rma_Simp.Biota_SA3=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA3)
res.rma_Simp.Biota_SA4=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA4)
res.rma_Simp.Biota_SA5=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA5)

# Turnover:
res.rma.Turn.Biota_SA1=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA1)
res.rma.Turn.Biota_SA2=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA2)
res.rma.Turn.Biota_SA3=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA3)
res.rma.Turn.Biota_SA4=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA4)
res.rma.Turn.Biota_SA5=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~TaxonomicGroup-1, data= DATASET_SA5)



## Biogeoregions:
# sample 6 Alpine plants, and 8 Atlantic terrestrial invertebrates

#Subset only 8 Terrestrial invertebrates from Atlantic region:
set.seed(200)
DATA_Atlantic_TerrInv_Sel1 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=8, replace=F)

set.seed(211)
DATA_Atlantic_TerrInv_Sel2 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=8, replace=F)

set.seed(222)
DATA_Atlantic_TerrInv_Sel3 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=8, replace=F)

set.seed(233)
DATA_Atlantic_TerrInv_Sel4 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=8, replace=F)

set.seed(244)
DATA_Atlantic_TerrInv_Sel5 <- sample(DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT","TimeSeries"], size=8, replace=F)


#Subset only 6 plants from Alpine region:
set.seed(200)
DATA_Alpine_Vegetation_Sel1 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=6, replace=F)

set.seed(211)
DATA_Alpine_Vegetation_Sel2 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=6, replace=F)

set.seed(222)
DATA_Alpine_Vegetation_Sel3 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=6, replace=F)

set.seed(233)
DATA_Alpine_Vegetation_Sel4 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=6, replace=F)

set.seed(244)
DATA_Alpine_Vegetation_Sel5 <- sample(DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation","TimeSeries"], size=6, replace=F)


# Combine DATA datasets  with only selected sites for the over-represented groups:

DATA_Atlantic_TerrInv <- DATA[DATA$Biogeoregion=="Atlantic" & DATA$TaxonomicGroup == "InvertebratesT",]
DATA_without_Atlantic_TerrInv <- DATA[!DATA$TimeSeries %in% DATA_Atlantic_TerrInv$TimeSeries,]
DATA_Atlantic_TerrInv_Sel_1 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel1,] )
DATA_Atlantic_TerrInv_Sel_2 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel2,] )
DATA_Atlantic_TerrInv_Sel_3 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel3,] )
DATA_Atlantic_TerrInv_Sel_4 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel4,] )
DATA_Atlantic_TerrInv_Sel_5 <- rbind(DATA_without_Atlantic_TerrInv,DATA[DATA$TimeSeries %in% DATA_Atlantic_TerrInv_Sel5,] )

DATA_Alpine_Vegetation <- DATA[DATA$Biogeoregion=="Alpine" & DATA$TaxonomicGroup == "Vegetation",]
DATA_without_Alpine_Vegetation <- DATA[!DATA$TimeSeries %in% DATA_Alpine_Vegetation$TimeSeries,]
DATA_Alpine_Vegetation_Sel_1 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel1,] )
DATA_Alpine_Vegetation_Sel_2 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel2,] )
DATA_Alpine_Vegetation_Sel_3 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel3,] )
DATA_Alpine_Vegetation_Sel_4 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel4,] )
DATA_Alpine_Vegetation_Sel_5 <- rbind(DATA_without_Alpine_Vegetation,DATA[DATA$TimeSeries %in% DATA_Alpine_Vegetation_Sel5,] )


# Re-run models
# Run the models below 2 times, after selecting the data referring to (1) Atlantic_TerrInv and (2) Alpine_Vegetation:

DATASET_SA1 <- DATA_Atlantic_TerrInv_Sel_1 # DATA_Alpine_Vegetation_Sel_1
DATASET_SA2 <- DATA_Atlantic_TerrInv_Sel_2 # DATA_Alpine_Vegetation_Sel_2
DATASET_SA3 <- DATA_Atlantic_TerrInv_Sel_3 # DATA_Alpine_Vegetation_Sel_3
DATASET_SA4 <- DATA_Atlantic_TerrInv_Sel_4 # DATA_Alpine_Vegetation_Sel_4
DATASET_SA5 <- DATA_Atlantic_TerrInv_Sel_5 # DATA_Alpine_Vegetation_Sel_5

# Abundance:
res.rma.Abund.Biota_SA1=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA1)
res.rma.Abund.Biota_SA2=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA2)
res.rma.Abund.Biota_SA3=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA3)
res.rma.Abund.Biota_SA4=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA4)
res.rma.Abund.Biota_SA5=rma(yi = Abund_S, vi = Abund_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA5)

# Richness:
res.rma.NTaxa.Biota_SA1=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA1)
res.rma.NTaxa.Biota_SA2=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA2)
res.rma.NTaxa.Biota_SA3=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA3)
res.rma.NTaxa.Biota_SA4=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA4)
res.rma.NTaxa.Biota_SA5=rma(yi = NTaxa_S, vi = NTaxa_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA5)

# Simpson´s diversity:
res.rma_Simp.Biota_SA1=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA1)
res.rma_Simp.Biota_SA2=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA2)
res.rma_Simp.Biota_SA3=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA3)
res.rma_Simp.Biota_SA4=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA4)
res.rma_Simp.Biota_SA5=rma(yi = Simp_S, vi = Simp_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA5)

# Turnover:
res.rma.Turn.Biota_SA1=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA1)
res.rma.Turn.Biota_SA2=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA2)
res.rma.Turn.Biota_SA3=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA3)
res.rma.Turn.Biota_SA4=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA4)
res.rma.Turn.Biota_SA5=rma.uni(yi = Turn_S, vi = Turn_var, method="REML", mods=~Biogeoregion-1, data= DATASET_SA5)

#Title: Using landscape genomics to delineate seed and breeding zones for Lodgepole pine
#Code contributor: Yue Yu and Tongli Wang
#Contact info: yue.yu@ubc.ca

#This code includes the following five sections
# 1. Build GF model
#   1.1 Load genomic data
#   1.2 Load climate data
#   1.3 Merge data
#   1.4 GF model training
# 2. GF model outputs
#   2.1 R square (min,max,mean)
#   2.2 Predictor importance
#   2.3 Split importance
#   2.4 Cumulative importance
# 3. GF model prediction
#   3.1 Load climate data in raster format
#   3.2 Raster convert to dataframe
#   3.3 Transform
#   3.4 PCA
#   3.5 Map spatial pattern of predicted genomic composition
# 4. Determine number of zones
# 5. Seed and breeding zone delineation
#   5.1 Cluster into zones
#   5.2 GF-based seed and breeding zones
#   5.3 Convert map into raster format
#   5.4 Comparison with common garden-based zones in ArcGIS


#NOTES:
#File directories may vary, so make the changes in the code accordingly
#Typos might exist


#=========================
# Code starts here
#=========================

#=========================
# 1. Build GF model
#=========================
# Here we only use the full set as a coding example
# But the same code could be applied to two other sets of SNPs (Neutral and Candidate)

#-------------------------
# 1.1 Load genomic data
#-------------------------
# Load SNP data (Full/ Candidate/ Neutral SNP sets in genotype)
# Lodgepole pine SNPs data are archived at https://doi.org/10.5061/dryad.56j8vq8
# allSnp.csv is the Lodgepole pine SNPs found in https://doi.org/10.5061/dryad.56j8vq8
all <- read.csv('allSnp.csv');head(all);dim(all) 
#Keep SNPs with allele frequency > 5% and saved as variable c4

#-------------------------
# 1.2 Load climate data
#-------------------------
#climate data obtained from ClimateNA for 20 climate variables
getwd()
climate0 <- read.csv("originalData/output_climate1.csv");head(climate0);dim(climate0)


#-------------------------
# 1.3 Merge into one dataframe
#-------------------------
#merge genomic data and climate data based on sample location (ID) 
rawdatafinal_A <- merge(c4,climate0, all.x = T, by = c("ID"))

#-------------------------
# 1.4 GF model training
#-------------------------
#extract climate data as predictor variables, genomic data as response variables
predictor.vars <- colnames(rawdatafinal_A[,28955:28974]);predictor.vars 
#20 climate variables
response.vars <- colnames(rawdatafinal_A[,1:28954]);response.vars[1:10]
#28,954 SNPs 
imp.vars <- predictor.vars
#For GF model parameters, see Ellis et al., 2012
gf_all <- gradientForest(rawdatafinal_A,predictor.vars, response.vars,
                      ntree = 200, maxLevel = level,trace = T, corr.threshold = 0.5)
#view GF model
gf_all; importance(gf_all)
#With the same code, we could generate gf_EP for the selected candidate set and gf_control for the neutral set.
#For more detailed description on the GF model see Ellis et al., 2012



#=========================
# 2. GF model outputs
#=========================

#-------------------------
# 2.1 R square (min,max,mean)
#-------------------------
gf <- gf_all #or could assign gf_EP to see GF outputs for the candidate set.

sqar <- gf$imp.rsq
sqar[sqar == 0] <- NA
sqar[1:5,1:5]

max(sqar,na.rm = T)
min(sqar,na.rm = T)
mean(sqar,na.rm = T)

#-------------------------
# 2.2 Predictor importance
#-------------------------
names(gf)
plot(gf,plot.type = "O")
most_important <- names(importance(gf))[1:20] #20 environmental variables
par(mgp = c(2,0.75,0))

#-------------------------
# 2.3 Split importance
#-------------------------
plot(gf,plot.type = "S",imp.vars = most_important, leg.posn="topright", cex.legend=0.4, 
        cex.axis=0.6, cex.lab=0.7, line.ylab=0.9,
        par.args=list(mgp=c(1.5, 0.5, 0), mar=c(3.1,1.5,0.1,1)))

#-------------------------
# 2.4 Cumulative importance
#-------------------------
plot(gf, plot.type="C", imp.vars=most_important, show.species=F, common.scale=T,
        cex.axis=0.6, cex.lab=0.7, line.ylab=0.9, 
        par.args=list(mgp=c(1.5, 0.5, 0), mar=c(2.5,1.0,0.1,0.5), omi=c(0,0.3,0,0)))



#=========================
# 3. GF model prediction
#=========================

#-------------------------
# 3.1 Load climate data for areas to be predicted (in raster format)
#-------------------------
# For this research, it is the climate data (800mx800m resolution) for the distibution range of 
# Lodgepole pine in British Columbia and Alberta in Canada.

#Load packages
library(raster)
library(rgdal)
library(sf)


#Load outlines
#Lpine outline
range <- readOGR("C:/Users/hellooo/Documents/LpineRange/pc_range.shp")
plot(range)
range2 <- spTransform(range, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83
+units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
plot(range2) #distribution outline with projection
#British columbia and Alberta (ABBC) outline
outline <- readOGR("C:/Users/hellooo/Documents/bcab800/ABBC_outline.shp")
plot(outline) # ABBC outline with NO projection
outlineABBC <- spTransform(outline, CRS("+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83
+units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"))
plot(outlineABBC) #ABBC outline with projection


#Function to read tif files
toStk3 <- function(x,varList,rType='tif',vConvert=F){
   library(raster)
   for(var in varList){
      if(rType=='asc'){inF <- paste0(x,'/',var, '.asc')}
      r <- raster(inF)
      if(vConvert==T){
if(names(r)=='MAT'|names(r)=='MWMT'|names(r)=='MCMT'|names(r)=='TD'|names(r)=='EMT'|names(r)=='EXT'|names(r)=='AHM'|names(r)=='SHM'|names(r)=='MAR'){r <- r/10}
      }
      if(var==varList[1]){stk=r} else{stk=stack(stk,r)}
   }
   return(stk)
}

#Define list of 20 climate variables
varList <- c("MAT", "MWMT", "MCMT", "TD", "MAP", "MSP", "AHM", "SHM", "DD_0", "DD5", "NFFD", "FFP", 
            "bFFP", "eFFP","PAS", "EMT", "EXT", "Eref", "CMD", "RH") # For climate variables descriptions, see Table 1.

wd <- 'C:/Users/hellooo/Documents/800m_abbc/Normal_1961_1990Y'

#Load climate data using the function
stk_abbc2 <- toStk3(wd,varList,rType='asc',vConvert=T)
stk_abbc2
#assign crs
crs(stk_abbc2) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
stk_abbc2

#crop the ABBC area out with the NONE projected outline
stk_abbc3 <- crop(stk_abbc2, extent(range))
stk22 <- mask(stk_abbc3, range)
crs(stk22) <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
plot(stk22, 'MAT')

#add projection, final output (stk4 for ABBC)
stk44 <- projectRaster(stk22, crs ="+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83
+units=m +no_defs +ellps=GRS80 +towgs84=0,0,0")
stk44  #stk4 was used for ABBC historical 800m
plot(stk44, 'MAT')
plot(outlineABBC, add = TRUE,lwd = 0.005)


#-------------------------
# 3.2 Raster convert to dataframe (climate data to be predicted)
#-------------------------
csv <- rasterToPoints(stk44)
df <- as.data.frame(csv)
#tailor the data for input
names(df)[1] <- "Longitude"
names(df)[2] <- "Latitude"

#-------------------------
# 3.3 Transform (model predict)
#-------------------------
#vars: all climate variables used in the fitted gradient forest model
df_pca <- cbind(df[,c("Longitude","Latitude")],
                predict(gf,df[,vars]))


#-------------------------
# 3.4 PCA
#-------------------------
#Notes:
#vars: all climate variables used in the fitted gradient forest model
#df:climate information for all coordinates, ready to input for gf model prediction
#gf: the final fitted gradient forest model 
#df_pca: a dataframe containing all GF predicted outputs, ready for PCA

PC <- prcomp(df_pca[,vars])  
pcx <- PC$x

#Assign PC
pc1 <- pcx[,1]  
pc2 <- pcx[,2]
pc3 <- pcx[,3]

#define RGB color palette (your choice)
r <- pc2
g <- pc3+pc1-pc2
b <- pc3-pc2
r <- (r-min(r))/(max(r)-min(r))
g <- (g-min(g))/(max(g)-min(g))
b <- (b-min(b))/(max(b)-min(b))

summary(r)
summary(g)
summary(b)

#Biplot
plot(pcx[,1:2],pch = ".", cex = 1,col = rgb(r,g,b),asp = 1)

#plot arrows (climate variables) on the biplot
vec <- c("MCMT","TD","MSP","DD_0","PAS","EMT","Eref","CMD") #Here we chose some important climate variables
lv <- length(vec)
class(vind) 
vind <- rownames(PC$rotation) %in% vec

arrow1 <- PC$rotation[c(3,4,6,9,15,16,18,19),1]  
arrow2 <- PC$rotation[c(3,4,6,9,15,16,18,19),2]
arrow_scale <- 45  #set a scale for the length of the arrow
plot(pcx[,1:2],pch = ".", cex = 2,col = rgb(r,g,b),asp = 1)
arrows(rep(0,lv),rep(0,lv),arrow1/arrow_scale, arrow2/arrow_scale,length = 0.0625)
jit <- 0.0012 #distance between the text to the arrow
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec)


#-------------------------
# 3.5 Map spatial pattern of predicted genomic composition
#-------------------------
#Plot spatail map with the GF-predicted results
plot(df_pca[,c("Longitude","Latitude")],pch = ".",cex =0.1 ,asp = 1,col = rgb(r,g,b))
plot(outlineABBC,add = TRUE, lwd = 0.005)


#=========================
# 4. Determine number of zones
#=========================
# Run this section in the R 4.2.1 or higher version
#because some packages does not work in R 3.5.1 version
#Load packages
library("knitr")
library("ggplot2")
library("factoextra")
library("cluster")

#This function (fviz_nbclust) find the within cluster variation
#variation stands for the within cluster variation
f <- fviz_nbclust(pcx,clara, method = "wss", k.max = 16)
variation <- f$data$y;f2
#then calculate the reduction in within cluster variation (reduc_var)


#=========================
# 5. Seed and breeding zone delineation
#=========================

#-------------------------
# 5.1 Cluster points into zones
#-------------------------
ncl <- 6 #number of zones determined based on section 4 results
clPCs <- clara(pcx,ncl,sampsize=10000)

#set up the medoid color palette
medcolR <- r[clPCs$i.med]
medcolG <- g[clPCs$i.med]
medcolB <- b[clPCs$i.med]

summary(medcolR)
summary(medcolG)
summary(medcolB)

#PCA biplot into groups
plot(pcx[,1:2], pch = ".", cex = 1, 
     col=rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]),asp = 1)
arrows(rep(0,lv),rep(0,lv),arrow1/arrow_scale, arrow2/arrow_scale,length = 0.0625)
text(arrow1/arrow_scale+jit*sign(arrow1),arrow2/arrow_scale+jit*sign(arrow2), labels = vec)


#-------------------------
# 5.2 GF-based seed and breeding zones
#-------------------------
plot(df_pca[,c("Longitude","Latitude")],pch = ".",
     cex = 0.1,asp = 1,col = rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering]), 
    main = "", xlim = c(300000,2200000))
legend("bottomleft",as.character(seq(1,ncl)), pch=15, cex=1, col=rgb(medcolR,medcolG,medcolB))
points(df_pca[clPCs$i.med,c("Longitude","Latitude")], pch=as.character(seq(1,ncl)))
plot(outlineABBC, add = TRUE)


#-------------------------
# 5.3 Convert map into raster format
#-------------------------
# load color palette
hd(medcolR)
hd(medcolG)
hd(medcolB)
mix <- rgb(medcolR[clPCs$clustering], medcolG[clPCs$clustering], medcolB[clPCs$clustering])
#mix is the colors for each location

summary(mix)
str(mix)
unique(mix)

mix2 <- gsub("#8D3B4B", "1", mix) #extract color code
mix3 <- gsub("#80264B", "2", mix2)
mix4 <- gsub("#654C76", "3", mix3)
mix5 <- gsub("#643868", "4", mix4)
mix6 <- gsub("#4289A9", "5", mix5)
mix7 <- gsub("#58AC6E", "6", mix6)

mix9 <- as.numeric(mix7);mix9
unique(mix9)

mixcluster <- df_pca[,c("Longitude","Latitude")]
mixcluster$color <- mix9
hd(mixcluster)
#save in CSV format
fWrite(mixcluster,"rzone666_20210805.csv")

#Convert into raster
rzone6 <- rasterFromXYZ(mixcluster);rzone6 
crs(rzone6) <- "+proj=aea +lat_1=50 +lat_2=58.5 +lat_0=45 +lon_0=-126 +x_0=1000000 +y_0=0 +datum=NAD83 +
            units=m +no_defs +ellps=GRS80 +towgs84=0,0,0"
rzone6
#save it
writeRaster(rzone6, filename = "rzone6_20210805.grd") #grids
writeRaster(rzone6, filename = "rzone6_20210805.tiff")

#MOVE TO ARCGIS FROM HERE ON

#-------------------------
# 5.4 Comparison with common garden-based zones
#-------------------------
#This section is done in ArcGIS.
#See Liepe et al. 2016 and Ukrainetz et al. 2018 for shapefiles of the common garden-based zones.


#=========================
# CODE END
#=========================

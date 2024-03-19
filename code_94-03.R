#Code credit: Dr. Venki Uddameri
#Help taken from: ChatGPT and Copilot

#step1: Load Libraries
library(terra) 
library(gdata)
library(corrplot)
library(sf) 
library(pracma)
library(stats)
library(ggplot2)
library(cluster)
library(factoextra)
library(stats)

#step 2: set working directory
path <- 'C:/RS_Project/Merged bands/Merged 94-03'
setwd(path)

#step 3: Read all bands

bands <- c("aerosolx", "bluex", "greenx", "redx", "NIRx", "SWIR2x")
Nbands <- length(bands)

#Step 4 <- Read the CAL vector file
CAL <- vect("CAL_RP_94-03.gpkg") 
plot(CAL) 


#Step5 : Read the raster and crop the data
for( i in seq(1, Nbands, 1)){
  fname <- paste('Band', i, '.tif', sep = "" ) ##################Error due to 6 bands
  x <- rast(fname)
  x <- crop(x,CAL)
  mv('x', bands[i]) 
}

#step 6: Create a raster stack with all data
callsat <- c(aerosolx, bluex, greenx, redx, NIRx, SWIR2x)
writeRaster(callsat, "calallbands.TIF", filetype = "GTiff", overwrite = TRUE)

#Plot RGB Data bands
plotRGB(callsat, r= 4, g=3, b=2, stretch = 'hist') 
#false color plot
plotRGB(callsat, r= 5, g=4, b=3, stretch = 'hist')

#step 7: Create Boxplot
bandnames <- c("aerosolx", "bluex", "greenx", "redx", "NIRx", "SWIR2x")
f <- boxplot(callsat, axes = FALSE, outline = F, ylab = "Value", notch = FALSE) 
axis(1, at = 1:(Nbands), labels = bandnames, las=2) 
axis(2)
title(('Digital Numbers for Various bands in Calcasieu'))
grid()
box()

#Step 8: write the bands as a dataframe for aditional analysis
banddf <- as.data.frame(callsat, xy = T) 
head(banddf)
banddf <- na.omit(banddf) 
colnames(banddf) <- c('X', 'Y', bandnames)
head(banddf)

#Step 9: correlation
corbands <- cor(banddf[,3:8], method = 'pearson') 
corrplot(corbands, method = 'circle', type = 'lower', diag = F)

#Clamp values for contrast enhancement #histogram equalization
hist(banddf[, 3], main = "Histogram for column 3")
box()

#Step 11: Clamp values for contrast enhancement #histogram equalization
LL <- 0.05
UL <- 0.95
bandnamec <- c()
for (i in seq(1, Nbands, 1)){
  bname <- paste(bands[i], 'c', sep = "")
  xmin <- quantile(banddf[, (i+2)], LL, na.rm = T)    ########Error
  xmax <- quantile(banddf[, (i+2)], UL, na.rm = T)
  x <- clamp(callsat[[i]], lower = xmin, upper = xmax, values = T)
  mv('x', bname)
  bandnamec[i] <- bname
}

#plot RGB and FCC Plots
lsatcalc<- c(aerosolxc, bluexc, greenxc, redxc, NIRxc, SWIR2xc)

#Plot RGB Databands
plotRGB(lsatcalc, r= 4, g=3, b=2, stretch = 'lin') #before enhance
plotRGB(callsat, r= 4, g=3, b=2, stretch = 'lin') #after enhance

#histogram
hist(lsatcalc)

#Step 12:

#create a dataframe
bandcdf <- as.data.frame(lsatcalc, xy = TRUE, geom = 'WKT')
summary(bandcdf)
bandcdf <- na.omit(bandcdf)
colnames(bandcdf) <- c('X', 'Y', bandnames, 'WKT') 

#compute correlation between bands
corbands <- cor(bandcdf[,3:8], method = 'spearman')
corrplot(corbands, method = 'number', type = 'lower', diag = FALSE)

#Reproject data to lat-lon
crsaea <- crs(lsatcalc, proj = T) #Get crs
crs84 <- 4326 #EPSG code #Define the new crs

banddf.SP <- st_as_sf(bandcdf, coords = c('X', 'Y'), crs = crsaea)
banddf.SP$XAEA = st_coordinates(banddf.SP)[,1]
banddf.SP$YAEA = st_coordinates(banddf.SP)[,2]

bandcdf.SP <- st_transform(x = banddf.SP, crs = crs84) #reproject
bandcdf.SP$Lon = st_coordinates(bandcdf.SP)[,1]
bandcdf.SP$Lat = st_coordinates(bandcdf.SP)[,2]

#create a csv file
bandcdf.SP <- subset(bandcdf.SP, select = -c(WKT, geometry))
write.csv(bandcdf.SP, 'bandcdf.csv', row.names = F)

# Perform PCA with bands only
pca <- prcomp(bandcdf[,3:8], scale = TRUE)

# Get the first two principal components

scaled_data <- scale(bandcdf[, 3:8])
pca_result <- prcomp(scaled_data)
summary(pca_result)
pcadata <- pca$x
write.csv(pcadata[,1:3], 'pcdata.csv', row.names = F)

# Look at the first few rows of the rotation matrix
head(pca$rotation)
pc1 <- pca$x[ ,2]
pc2 <- pca$x[ ,3]

# Visualize the data points in the first two principal components
ggplot(bandcdf[, 3:8], aes(x = pc1, y = pc2, color = "darkblue")) +
  geom_point() +
  labs(title = " first two principal components",
       x = "PC1", y = "PC2")

#perform clustering
bpcaclust <- kmeans(pcadata[,1:3], centers=5, nstart = 10)
summary(bpcaclust)
fviz_cluster(bpcaclust, data = pcadata[,1:3],
             geom = "point",
             ellipse.type = 'convex')

bpcclust <- kmeans(pcadata[,1:3],5,10)
clusts <- bpcclust$cluster
pcadat.clust <- data.frame(bandcdf.SP,clusts)
write.csv(pcadat.clust,'bandcdfclus.csv',row.names=F)





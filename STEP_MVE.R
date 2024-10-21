# Packages needed ---------------------------------------------------------
library ("easypackages")
libraries (c("dismo", "maps", "maptools", "rgdal", "rgeos", "sp", "virtualspecies", "corrplot", 
             "raster", "usdm", "gtools", "tidyverse", "rJava", "ENMeval", "rgbif",
             "ENMTools", "geodata", "devtools", "RStoolbox", "glue", "rmaxent", "ENMeval",
             "sf", "xlsx", "rgl"))
library(ellipsenm)
library(RSQLite)

# path to save the information ----------------------------------------------------
path_wd<-"F:/Species_features/MVE"
path_species<-"D:/00013-Tesis_Doctoral/2- Base de datos/General_DB_Coreidae.xlsx"

# Open the climatic variables ---------------------------------------------------------
rasters<- worldclim_global("bio", res=5, 
                           path="C:/Users/amaru/OneDrive/Escritorio/worldclim")

#Open the databases needed for the analysis-------------------------------------------------
conn <- dbConnect(SQLite(), 
                  "D:/00013-Tesis_Doctoral/2- Base de datos/bd_chinches.db")
tablas <- dbListTables(conn)
df <- dbReadTable(conn, "Records_per_species") 

# remove species with leptoglossus gonagra and leptoglossus occidentalis------------------
df= df[!(df$Species=="Leptoglossus gonagra"),] # 1st species out of the american continent
df= df[!(df$Species=="Leptoglossus occidentalis"),] ### 2sd species out of the ameruican continent
df=df[!is.na(df["Longitude"]),] ## remove cells without any value in the databases

# Crop and Mask climatic variables -------------------------------------------------------
points=vect(df[c("Longitude", "Latitude")],geom=c("Longitude", "Latitude"))
range<-terra::convHull(points) # create the convexhull
crs(range)<-crs(rasters) ##3 change the Crs to a coordinate crs value
range<-buffer(x=range, width=200000 # create the buffer around each convhull
              , quadsegs=10, singlesided=FALSE)

# crop and mask
climatic_part<-crop(rasters, range)
climatic_part<-mask(climatic_part, range)

# Make the PCA of the archive --------------------------------------------------------
pca1<-raster.pca(climatic_part, 3)
writeRaster(pca1$rasters, glue("{path_wd}/MVE_PCA.tif"), 
            filetype = "GTiff",overwrite=T)

# calculations to store in the information
prop_varianza <- pca1$pca.object$sdev^2 / sum(pca1$pca.object$sdev^2) # make the calculation of the variance
informacion_pca<- data.frame("Standar Deviation"=pca1$pca.object$sdev,
                             "prop_variance"=prop_varianza,
                             "variance_acumulate"=cumsum(prop_varianza))
write.csv(informacion_pca, paste(path_wd,"/information_PCA_MVE.csv", sep=""))
save(pca1, file = paste(path_wd,"/Object_PCA_MVE.RData", sep="")) 



# calculated the number of individuals for each species----------------------------
# eliminate duplicates considering one point in each cell in the maxent model
table<-read.xlsx(path_species, "Sheet1")

volume=c()
archive=list()
length_species=c()
number=0
for ( e in table$Species_name){
  data_species=df[df$Species==e,]
  cells=cellFromXY(pca1$rasters,xy=data_species[c("Longitude","Latitude")])
  data_species$cells=cells
  data_species<-data_species[!duplicated(data_species$cells),] ### eliminate duplicates in each cells of the environmental layer
  length<-length(data_species$cells) # with this we calculate the length of separate points
  length_species<-c(length_species, length) ### add to the list neccessary with the number of ocurrences
  set.seed(100)
  if (length>=10){
    ellips <- ellipsoid_fit(data = data_species, 
                            longitude = "Longitude",
                            latitude = "Latitude", 
                            method = "mve1",
                            level = 95, 
                            raster_layers = stack(pca1$rasters))
    volume<-c(volume,ellips@niche_volume)
  }
  if (length<10){
    ellips<-"No calculations available"
    volume<-c(volume, 0) ### The species with a zero value are the one that do not have too much points
  }
  archive[e]<-ellips
  number=number+1
}


table$volume<-volume
table$length<-length_species

write.xlsx(table,path_species) # note all except leptoglossus gonagra and Leptoglossus occidentalis
save(archive, file = paste(path_wd,"/ellipsoids_per_species.RData", sep="")) 


#cor(table$volume, table$length) # we could notice that there are an positive correlation between the number of points and the volume of their niche
#rgl::plot3d(rgl::ellipse3d(ellips@covariance_matrix,centre = ellips@centroid,level = 0.95),alpha=0.9,col="black")
#niche1 <- overlap_object(data.frame(occurrences1)[c("Species","Longitud","Latitud")], species =  "Species", longitude = "Longitud", latitude = "Latitud", method = "mve1", level =95, variables = vars)
#overlap <- ellipsoid_overlap(niche1,niche2, niche3, overlap_type = "full")
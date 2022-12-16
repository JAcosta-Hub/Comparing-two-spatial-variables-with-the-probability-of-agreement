##########################################
##--------------------------------------##
## Librerias

 library(GeoModels)
 library(OpenImageR)
 library(jpeg)

##--------------------------------------##
## Upload Original Images 
# Specify correctly the directory where the images are located.

im2008=readJPEG("Data Images Selected/Original/harvard_2008_10_15_133138.jpg")
im2009=readJPEG("Data Images Selected/Original/harvard_2009_10_12_120136.jpg")
im2010=readJPEG("Data Images Selected/Original/harvard_2010_10_13_120935.jpg")
im2011=readJPEG("Data Images Selected/Original/harvard_2011_10_15_130138.jpg")
im2012=readJPEG("Data Images Selected/Original/harvard_2012_10_13_120136.jpg")
im2013=readJPEG("Data Images Selected/Original/harvard_2013_10_14_132032.jpg")
im2014=readJPEG("Data Images Selected/Original/harvard_2014_10_12_120138.jpg")
im2015=readJPEG("Data Images Selected/Original/harvard_2015_10_15_120138.jpg")
im2016=readJPEG("Data Images Selected/Original/harvard_2016_10_14_120138.jpg")
im2017=readJPEG("Data Images Selected/Original/harvard_2017_10_13_120138.jpg")
im2018=readJPEG("Data Images Selected/Original/harvard_2018_10_14_133137.jpg")
im2019=readJPEG("Data Images Selected/Original/harvard_2019_10_15_120137.jpg")
im2020=readJPEG("Data Images Selected/Original/harvard_2020_10_15_122906.jpg")
im2021=readJPEG("Data Images Selected/Original/harvard_2021_10_13_122907.jpg")
im2022=readJPEG("Data Images Selected/Original/harvard_2022_10_15_135906.jpg")

imagenes=list(im2008=im2008, im2009=im2009, im2010=im2010,
              im2011=im2011, im2012=im2012, im2013=im2013,
              im2014=im2014, im2015=im2015, im2016=im2016,
              im2017=im2017, im2018=im2018, im2019=im2019,
              im2020=im2020, im2021=im2021, im2022=im2022)

##---------------------------------------##
## Image cropping


 dim(im2008) # dimentions
 
 # Choice of pixels
 dim1=300:960 
 dim2=250:820

 # Example 
 plot(NA, xlim=c(0,1), ylim=c(0,1), frame.plot=FALSE,  xaxt="n", yaxt="n", xlab="", ylab="")
 OpenImageR::imageShow(imagenes[[7]])
 plot(NA, xlim=c(0,1), ylim=c(0,1), frame.plot=FALSE,  xaxt="n", yaxt="n", xlab="", ylab="")
 OpenImageR::imageShow(imagenes[[7]][dim1, dim2,])
 
 # Name of the cutout images
 new_name<-paste0("Data Images Selected/Recortadas/",
                  paste("new","harvard",2008:2022,10,sep="_"),".jpg")
 new_name2 <- paste0("new.im",2008:2022)

 new_images<-list()
 for(i in 1:length(imagenes)){
  new_images[[new_name2[i]]] <- imagenes[[i]][dim1, dim2,]
  jpeg(filename=new_name[i],
      width = 820-250, height = 960-300, quality = 100)
  OpenImageR::imageShow(new_images[[i]])
  dev.off()
 }

 plot(NA, xlim=c(0,1), ylim=c(0,1), frame.plot=FALSE,  xaxt="n", yaxt="n", xlab="", ylab="")
 OpenImageR::imageShow(new_images[[7]])
 
 ##-------------------------------------##
 ## Rasterización
 
 # Name of raster images
 raster_name<-paste0("Data Images Selected/Rasterizadas/",
                     paste("raster","harvard",2008:2022,10,sep="_"),".jpg")
 
 raster15_images <- rasterization(new_images, px=15, raster_name=raster_name)

 # Examples
 k=16 # (choose between 1-15)
 plot(NA, xlim=c(0,1), ylim=c(0,1), frame.plot=FALSE,  xaxt="n", yaxt="n", xlab="", ylab="")
 OpenImageR::imageShow(new_images[[k]])
 plot(NA, xlim=c(0,1), ylim=c(0,1), frame.plot=FALSE,  xaxt="n", yaxt="n", xlab="", ylab="")
 OpenImageR::imageShow(raster15_images[[k]]) 
 
 
 ##-------------------------------------##
 ## Save cropped and raster images in RData formt
 save(new_images, raster15_images, file="Imagenes.RData" ) 
 
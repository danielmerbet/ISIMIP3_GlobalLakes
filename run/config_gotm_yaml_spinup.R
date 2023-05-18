save_coordinates <- read.csv("/scratch/brussel/106/vsc10623/run/coord_area_depth.csv")
pixel <- 10
year_ini <- 1901
year_end <- 2019

gotmyaml <- readLines("gotm.yaml")
gotmyaml[3] <- paste0("   name: ", pixel)
gotmyaml[4] <- paste0("   latitude: ", save_coordinates$lat[pixel])
gotmyaml[5] <- paste0("   longitude: ", save_coordinates$lon[pixel])
gotmyaml[6] <- paste0("   depth: ", round(save_coordinates$depth[pixel],2))

gotmyaml[10] <- paste0("   start: ", year_ini-20, "-01-01 00:00:00") #adding -20 for spinup
gotmyaml[11] <- paste0("   stop: ", year_end, "-12-31 00:00:00")


if (save_coordinates$depth[pixel]<20){
  gotmyaml[14] <- paste0("   nlev: ", round(save_coordinates$depth[pixel]*2))
}else{
  gotmyaml[14] <- paste0("   nlev: ", round(save_coordinates$depth[pixel]))
}

if (save_coordinates$depth[pixel]<1){
  gotmyaml[6] <- "   depth: 1"
  gotmyaml[14] <- "   nlev: 4"
}

meteo <- read.table("meteo_file.dat")
meteotas <- mean(meteo[,6]) #air temperature column

if (meteotas<0){
  gotmyaml[26] <- "      t_1: 0"
  gotmyaml[28] <- "      t_2: 0"
}else{
  gotmyaml[26] <- paste0("      t_1: ", round(meteotas,2))
  gotmyaml[28] <- paste0("      t_2: ", round(meteotas,2))
}

writeLines(gotmyaml, con="gotm.yaml")
date_spinup <- as.character(seq(as.Date(paste0((year_ini-20),"-01-01")), as.Date(paste0((year_ini-1),"-12-31")), by=1))
meteo_spinup <-meteo[1:length(date_spinup),]
meteo_spinup[,1] <- date_spinup
meteo <- rbind(meteo_spinup, meteo)

#to fix model with date yyyymmdd 12:00:00
if(meteo[1,2]=="12:00:00"){
 meteo[,2] <- rep("00:00:00", nrow(meteo)) 
}

write.table(meteo, file="meteo_file.dat", quote=F, row.names=F, col.names=F, sep="\t") 
 

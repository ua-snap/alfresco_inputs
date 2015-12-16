require(sp)
require(maptools)

t <-read.table("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/extracted_files/grid_10min_reh.dat")

colnames(t)<-c("lat","lon","reh_01","reh_02","reh_03","reh_04","reh_05","reh_06","reh_07","reh_08","reh_09","reh_10","reh_11","reh_12")
head(t)
count=0
for(i in 3:ncol(t)){
	count=count+1
	# do the change the data thing here.
	#pc.lon <- t[,2]+180
	xyz <- data.frame(t[,2],t[,1],t[,i])
	pts.cru <- SpatialPoints(xyz[,1:2])
	xyz.df <- as.data.frame(xyz)
	colnames(xyz.df) <- c("lon","lat","rhum")
	pts.cruSPDF <- SpatialPointsDataFrame(pts.cru, xyz.df)
	names(pts.cruSPDF) <- c("lon","lat","rhum")

	if(nchar(count)<2){ month=paste("0",count,sep="")}else{month=paste(count,sep="")}

	writeSpatialShape(pts.cruSPDF, paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/hur_climatology/v2/","hur_cru_10min_",month,"_1961_1990", sep=""))

# 	r <- rasterize(pcll.cruSPDF, cru.current, data.frame(pcll.cruSPDF[,3]), fun=mean)
# 	writeRaster(r, filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/hur_climatology/v2/","hur_cru_10min_",month,"_1961_1990.tif"), overwrite=T)
}
library(raster)

# some needed functions for calculation
source("/big_storage/malindgren/AIEM/functions/CalcAnom.r")
source("/big_storage/malindgren/AIEM/functions/ListFiles.r")

# BEGIN TEMPERATURE ANOMALY CREATION --> CCCMA

historical <- brick("/Data/Base_Data/Climate/World/GCM_raw/IPCC_AR4_Data/historical/tas/tas_Amon_cccma-cgcm3-1-t47_historical_run1_185001-200012.nc")
future <- brick("/Data/Base_Data/Climate/World/GCM_raw/IPCC_AR4_Data/sresa1b/tas/tas_Amon_cccma-cgcm3-1-t47_sresa1b_run1_200101-210012.nc")
fullSeries <- stack(historical,future)

# the new fullSeries is from 1850-2100 at a monthly timestep
#((1961-1850)*12)+1 = 1333  <<<<  this is the begin point on the subset for 1961-1990
#((1990-1850)*12)+12 = 1692

fullSeries.anom <- calcAnom(fullSeries,1333,1692,absolute=T)

# write out the anomalies here or write in this step into the full code
monthList <- c("01","02","03","04","05","06","07","08","09","10","11","12")
count=0
for(y in 1850:2100){
	for(m in monthList){
	count=count+1
	print(count)
	writeRaster(subset(fullSeries.anom,count,T), filename=paste("/big_storage/malindgren/AIEM/RHUM/anomalies/tas/cccma/","tas_mean_cccma_cgcm3_1_absAnom_",m,"_",y,".tif",sep=""), overwrite=T, options="COMPRESS=LZW") # tas_mean_cccma_cgcm3_1_anom_01_2001
	
	}
}

# BEGIN RELATIVE HUMIDITY ANOMALY CREATION --> CCCMA

historical <- brick("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.hur_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")
future <- brick("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.hur_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc")
fullSeries <- stack(historical,future)

# the new fullSeries is from 1850-2100 at a monthly timestep
#((1961-1850)*12)+1 = 1333  <<<<  this is the begin point on the subset for 1961-1990
#((1990-1850)*12)+12 = 1692

fullSeries.anom <- calcAnom(fullSeries,1333,1692,absolute=F)

# write out the anomalies here or write in this step into the full code
monthList <- c("01","02","03","04","05","06","07","08","09","10","11","12")
count=0
for(y in 1850:2100){
	for(m in monthList){
	count=count+1
	print(count)
	writeRaster(subset(fullSeries.anom,count,T), filename=paste("/big_storage/malindgren/AIEM/RHUM/anomalies/hur/cccma/","hur_mean_cccma_cgcm3_1_prpAnom_",m,"_",y,".tif",sep=""), overwrite=T, options="COMPRESS=LZW") # tas_mean_cccma_cgcm3_1_anom_01_2001
	
	}
}

# BEGIN TEMPERATURE ANOMALY CREATION --> ECHAM

historical <- stack("/Data/Base_Data/Climate/World/GCM_raw/IPCC_AR4_Data/historical/tas/tas_Amon_mpi-echam5_historical_run1_186001-210012.nc", bands=1:1692)
future <- stack("/Data/Base_Data/Climate/World/GCM_raw/IPCC_AR4_Data/sresa1b/tas/tas_Amon_mpi-echam5_sresa1b_run1_200101-220012.nc", bands=1:1200)
fullSeries <- stack(historical,future)

# the new fullSeries is from 1850-2100 at a monthly timestep
#((1961-1860)*12)+1 = 1213  <<<<  this is the begin point on the subset for 1961-1990
#((1990-1860)*12)+12 = 1572

fullSeries.anom <- calcAnom(fullSeries,1213,1572,absolute=F)

# write out the anomalies here or write in this step into the full code
monthList <- c("01","02","03","04","05","06","07","08","09","10","11","12")
count=0
for(y in 1860:2100){
	for(m in monthList){
	count=count+1
	print(count)
	writeRaster(subset(fullSeries.anom,count,T), filename=paste("/big_storage/malindgren/AIEM/RHUM/anomalies/tas/echam5/","tas_mean_mpi_echam5_absAnom_",m,"_",y,".tif",sep=""), overwrite=T, options="COMPRESS=LZW") # tas_mean_cccma_cgcm3_1_anom_01_2001
	
	}
}

# BEGIN RELATIVE HUMIDITY ANOMALY CREATION --> ECHAM

historical <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.mpi_echam5.20c3m.run1.monthly.hur_A1_1960-2009.nc", bands=1:492)
future1 <- brick("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.mpi_echam5.sresa1b.run1.monthly.hur_A1_2001-2050.nc")
future2 <- brick("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.mpi_echam5.sresa1b.run1.monthly.hur_A1_2051-2100.nc")
fullSeries <- stack(historical,future1,future2)

# the new fullSeries is from 1850-2100 at a monthly timestep
#((1961-1960)*12)+1 = 13  <<<<  this is the begin point on the subset for 1961-1990
#((1990-1860)*12)+12 = 372

fullSeries.anom <- calcAnom(fullSeries,13,1572,absolute=F)

# write out the anomalies here or write in this step into the full code
monthList <- c("01","02","03","04","05","06","07","08","09","10","11","12")
count=0
for(y in 1960:2100){
	for(m in monthList){
	count=count+1
	print(count)
	writeRaster(subset(fullSeries.anom,count,T), filename=paste("/big_storage/malindgren/AIEM/RHUM/anomalies/hur/echam5/","hur_mean_mpi_echam5_prpAnom_",m,"_",y,".tif",sep=""), overwrite=T, options="COMPRESS=LZW") # tas_mean_cccma_cgcm3_1_anom_01_2001
	
	}
}




library(raster)

# some needed functions for calculation
source("/big_storage/malindgren/AIEM/functions/CalcAnom.r")
source("/big_storage/malindgren/AIEM/functions/ListFiles.r")

# BEGIN CLOUD ANOMALY CREATION --> CCCMA

historical <- brick("/Data/Base_Data/Climate/World/GCM_raw/clt/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.clt_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")
future <- brick("/Data/Base_Data/Climate/World/GCM_raw/clt/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.clt_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc")
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
	writeRaster(subset(fullSeries.anom,count,T), filename=paste("/big_storage/malindgren/AIEM/RSDS/anomalies/clt/cccma/","clt_mean_cccma_cgcm3_1_prpAnom_",m,"_",y,".tif",sep=""), overwrite=T, options="COMPRESS=LZW") # tas_mean_cccma_cgcm3_1_anom_01_2001
	
	}
}


# BEGIN CLOUD ANOMALY CREATION --> ECHAM5

library(raster)

# some needed functions for calculation
source("/big_storage/malindgren/AIEM/functions/CalcAnom.r")
source("/big_storage/malindgren/AIEM/functions/ListFiles.r")

# BEGIN CLOUD ANOMALY CREATION --> CCCMA

historical <- brick("/Data/Base_Data/Climate/World/GCM_raw/clt/pcmdi.ipcc4.mpi_echam5.20c3m.run1.monthly.clt_A1.nc")
historical <- subset(stack(historical),1:1692,T)
future <- brick("/Data/Base_Data/Climate/World/GCM_raw/clt/pcmdi.ipcc4.mpi_echam5.sresa1b.run1.monthly.clt_A1.nc")
future <- subset(stack(future), 1:1200, T)
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
	writeRaster(subset(fullSeries.anom,count,T), filename=paste("/big_storage/malindgren/AIEM/RSDS/anomalies/cld/echam5/10min/","clt_mean_mpi_echam5_prpAnom_",m,"_",y,".tif",sep=""), overwrite=T, options="COMPRESS=LZW") # tas_mean_cccma_cgcm3_1_anom_01_2001
	
	}
}



# BEGIN CLOUD ANOMALY CREATION --> TS31
library(raster)

# some needed functions for calculation
source("/big_storage/malindgren/AIEM/functions/CalcAnom.r")
source("/big_storage/malindgren/AIEM/functions/ListFiles.r")

# BEGIN CLOUD ANOMALY CREATION --> CCCMA
fullSeries <- brick("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/RSDS_working/September2012_finalRuns/10min_sunshine/original_data/TS31/cru_ts_3_10.1901.2009.cld.dat.nc")
fullSeries <- stack(fullSeries)
# the new fullSeries is from 1850-2100 at a monthly timestep
#((1961-1901)*12)+1 = 721  <<<<  this is the begin point on the subset for 1961-1990
#((1990-1901)*12)+12 = 1080

fullSeries.anom <- calcAnom(fullSeries,721,1080,absolute=F)

# write out the anomalies here or write in this step into the full code
monthList <- c("01","02","03","04","05","06","07","08","09","10","11","12")
count=0
for(y in 1901:2009){
	for(m in monthList){
		count=count+1
		print(count)
		writeRaster(subset(fullSeries.anom,count,T), filename=paste("/big_storage/malindgren/AIEM/RSDS/anomalies/cld/ts31/10min/","clt_mean_cru_ts31_prpAnom_",m,"_",y,".tif",sep=""), overwrite=T, options="COMPRESS=LZW") # tas_mean_cccma_cgcm3_1_anom_01_2001
	
	}
}






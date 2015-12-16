# TEM actually uses vapor pressure to drive the model. So basically the rel hum I sent to you was converted from vapor pressure. The conversion is conducted as following:
# (1) calculate saturated vapor pressure (esa, hPa or mbar) according to Meteorological Instruments and Methods of Observation (CIMO Guide)
#      (WMO, 2008):
#         esa = 6.112 exp(17.62*Tair/(243.12 + Tair)),  where Tair in °C and >0; 
#        OR, 
#         esa = 6.112 exp(22.46 *Tair/(272.62 + Tair)), where Tair  in °C and <=0

# (2) relative humility RH:
#        RH = ea/esa * 100, where ea is the vapor pressure in hPa or mbar.


# so we need to set up a couple of variables for the testing procedure

tas.2001 <- list.files("/Data/Base_Data/Climate/AK_CAN_1km_from2km/cccma_cgcm3_1/sresa1b/tas", pattern="*_2001.tif$", full.names=T)
rhum.2001 <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/downscaled/AKalbers", pattern="*_2001.tif$", full.names=T)


tas.year1 <- stack(tas.2001)

rhum.year1 <- stack(rhum.2001)

#here is the conversion procedure that is used to convert from hte input units of 
#	RHUM to the Needed units of Vapor Presure 
for(i in 1:nlayers(tas.year1)){

	#here we subset the stack to the layer of interest at this stage in the conversion
	tas.year1.sub <- subset(tas.year1, i, drop=T)
	rhum.year1.sub <- subset(rhum.year1, i, drop=T)

	# Conversion Procedure from Fengming vetted by Steph McAfee
	# source: WMO_CIMO_Guide-7th_Edition-2008
	ESA = 6.112 ^(22.46*tas.year1.sub / (272.62 + tas.year1.sub))
	EA = (ESA*rhum.year1.sub)/100

	# write out the new raster file
	writeRaster(EA, filename=paste(outPath,"",".tif", sep=""))
}




#### TESTING AREA #####

test <- raster(ncstack, layer=1)
test.xyz <- cbind(coordinates(test), getValues(test))
test.xyz
class(test.xyz)
as.data.frame(test.xyz)
test.xyz <- as.data.frame(test.xyz)

colnames(test.xyz) <- c("x","y","z")

m <- matrix(test.xyz$z, nrow=nrow(test), ncol=ncol(test))

if(any(coords[,1]>180)) coords[,1][coords[,1]>180] <- coords[,1][coords[,1]>180]-360

#-----------------------------------------------------------------------------
# here is how we will create the new flipped raster file from an XYZ table

# create an xyz table from 0-360 raster
xyz <- data.frame(coordinates(test),getValues(test))

# set the colnames to something simple
colnames(xyz) <- c("x","y","z")

# turn that 0-360 XYZ table into a -180-180 table
if(any(xyz[,1]>180)) xyz[,1][xyz[,1]>180] <- xyz[,1][xyz[,1]>180]-360

# turn that table into a raster
xyz.r <- rasterFromXYZ(xyz,res=res(test),crs=NA)

#-----------------------------------------------------------------------------

ncstack <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.hur_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc", xmn=0,xmx=360,ymn=-90,ymx=90)

if(any(c.f[,1]>180)) c.f[,1][c.f[,1]>180] <- c.f[,1][c.f[,1]>180]-360


testing <- projectRaster(xxx, newProj)


writeRaster(testing, filename="/workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/test.tif", overwrite=TRUE)


system("gdalwarp -t_srs '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs' /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test.tif /workspace/UA/malindgren/projects/iem/PHASE2_DATA/testing_dir/monthFlip_test_gdal.tif")

gdalCall <- paste("gdalwarp -t_srs '+proj=aea +lat_1=55 +lat_2=65 +lat_0=50 +lon_0=-154 +x_0=0 +y_0=0 +ellps=GRS80 +nadgrids=@alaska +datum=NAD83 +units=m +no_defs' ", <in_file>," " ,<out_file>, sep="")



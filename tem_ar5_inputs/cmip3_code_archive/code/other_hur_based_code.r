library(raster)

l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/cru_files/v2_rasterized", pattern=".tif$", full.names=T)

s <- stack(l)

for(i in 1:nlayers(s)){
	r <- raster(s, layer=i)
	name.r <- paste(substr(l[i], 1, nchar(l[i])-3),"rst", sep="")
	writeRaster(r, filename=name.r, overwrite=T)
}



l <- list.files("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/cru_files", pattern="_PCLL.tif", full.names=T)

s <- stack(l)


#r2 <- raster(ncstack,layer=1,xmn=xmin(cru.current), xmx=xmax(cru.current), ymn=ymin(cru.current), ymx=ymax(cru.current))

#OLD ncstack extent
class       : Extent
xmin        : -1.875
xmax        : 358.125
ymin        : -89.01354
ymax        : 89.01354

# give ncstack the bounding coordinates that it should have lon(0,360)
extent(ncstack) <- c(0,360,-89.01354,89.01354)

# create a new extent() object that is a union of the two existing inputs extents
u <- union(extent(ncstack), extent(s))

#expand the smaller one to the new union extent # in this case it is the cru_files
s.expand <- expand(s,u)



for(i in 1:nlayers(s.expand)){
	r.s.expand <- raster(s.expand, layer=1)
	r.v <- getValues(r.s.expand)
	r.s.expand2 <- r.s.expand
	ind <- which(r.v < 0); values(r.s.expand2)[ind] <- 0
	values(s.expand)

	#plot(r.s.expand2)
}	


l.names <- layerNames(s.expand2)
for(i in 1:nlayers(s.expand2)){

	writeRaster(subset(s.expand2,i), filename=paste("/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/cru_files/v2_fixed_oob/",substr(l.names[i],1,nchar(l.names[i])),".tif", sep=""), overwrite=T)
}

# this little code area is to get the ncstack and climstack into the right lon(min, max)
e <- extent(c(0,360,-89.01354,89.01354))

# read in the entire series of data from a netCDF using stack()
ncstack <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.sresa1b.run1.monthly.hur_a1_sresa1b_1_cgcm3.1_t47_2001_2100.nc") 
# read in the 20c3m stack
climstack <- stack("/Data/Base_Data/Climate/World/GCM_raw/hur/pcmdi.ipcc4.cccma_cgcm3_1.20c3m.run1.monthly.hur_a1_20c3m_1_cgcm3.1_t47_1850_2000.nc")


#ncstack.expand <- expand(ncstack,e)
#climstack.expand <- expand(climstack,e)
ncstack2 <- ncstack
climstack2 <- climstack

extent(ncstack2) <- e
extent(climstack2) <- e


writeRaster(subset(ncstack2,1), filename="/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/testing_extent/ncstack_1.tif", overwrite=T)
writeRaster(subset(climstack2,1), filename="/workspace/UA/malindgren/projects/iem/PHASE2_DATA/CRU_TS20/testing_extent/climstack_1.tif", overwrite=T)


if(any(test[,1]>180)) test[,1][test[,1]>180] <- test[,1][test[,1]>180]-360
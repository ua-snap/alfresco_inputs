library(raster)
library(sp)
library(akima)

# here is the iterator for poormans parallel
mon=12

# read in the .dat file from the cru group
xyz <- read.table("/Data/Base_Data/Climate/World/CRU_grids/CRU_TS20/grid_10min_tmp.dat")

# re-order the columns so that 1:2 are lon/lat not lat/lon
xyz <- cbind(xyz[,2], xyz[,1], xyz[,3:ncol(xyz)])

if(any(xyz[,1] < 0))  xyz[xyz[,1]<0,1] <- xyz[xyz[,1]<0,1]+360

# HARDWIRED
# make a template raster from a known extent of the ts20 maps
template <- raster(xmn=0, xmx=360, ymn=-59.08312, ymx=83.58371)
# change the res to match what we know is correct
res(template) <- 0.1666667
# END HARDWIRED

# this just gets the code running since I am tired of this
xyz.hur <- xyz

# input xy
in_xy <- xyz[,1:2]

# output xy
out_xy <- coordinates(template)

# now lets do the interp and rasterize it
print(mon)
z_in <- xyz[,mon+2]

hur.spline	<- interp(x=in_xy[,1],y=in_xy[,2],z=z_in,xo=seq(min(out_xy[,1]),max(out_xy[,1]),l=ncol(template)), yo=seq(min(out_xy[,2]),max(out_xy[,2]),l=nrow(template),linear=F))

mat <- matrix(hur.spline$z,nrow=nrow(template),ncol=ncol(template),byrow=T)[nrow(template):1,]

hur.interp.test <- raster(mat, xmn=xmin(template), xmx=xmax(template), ymn=ymin(template), ymx=ymax(template), crs=projection(template))

if(nchar(mon) == 1) { mon=paste("0",mon,sep="") }else{ mon=mon}

writeRaster(hur.interp.test, filename=paste("/big_storage/malindgren/AIEM/RHUM/tas_cru_10min_1961_1990/pcll_spline_final//","tmp_cru_ts20_",mon,"_1961_1990_spline.tif", sep=""), options="COMPRESS=LZW", overwrite=T)



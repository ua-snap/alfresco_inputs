esa.stack = 6.112^((22.46*t1) / (272.62 + t1))

hur.stack = (vap.stack/esa.stack)*100

# pull out the values of the stack

x <- getValues(subset(tas.stack,1,T))



test<-apply(tas.stack.v,1,FUN=function(x) if(x <= 0){	esa = 6.112 ^(22.46 *x/(272.62 + x))}else{	esa = 6.112 ^(17.62*x/(243.12 + x))})

outStack <- stack()

for(i in 1:nlayers(tas.stack)){
	print(i)
	cur <- subset(tas.stack,i,T)
	new.cur <- numeric()
	cur.v <- getValues(cur)
	for(v in cur.v){
		if(is.na(v) == FALSE){ 
			if(v > 0){
				esa <- 6.112 ^(17.62*v/(243.12 + v))
				new.cur <- append(new.cur,esa, after=length(new.cur))
			}else{
				esa <- 6.112 ^(22.46 *v/(272.62 + v))
				new.cur <- append(new.cur,esa, after=length(new.cur))
			}
		}else{
			new.cur <- append(new.cur,v, after=length(new.cur))

		}
	}
	values(cur) <- new.cur

	outStack<-addLayer(outStack, cur)
}


# now look at outStack and do some conversion with that thing.
# importatnt to note that the outStack used here is the same as the ESA calculated above.

vap.stack2 <- subset(vap.stack,1:17,T)

hur.stack = (vap.stack2 / outStack)*100


esa = 6.112^(17.62*t1/(243.12+t1))


testing = 6.112 exp(((17.62*t1)/(243.12+t1)))




6.112^(17.62*test[,1]/(243.12+test[,1]))



# your input raster
x=vap1

# get the xyz from it
xyz <- cbind(coordinates(x), getValues(x))

# flip the coordinates to pacific centered latlong
if(any(xyz[,1] < 0))  xyz[xyz[,1]<0,1] <- xyz[xyz[,1]<0,1]+360

# make a new raster object
x1 <- rasterFromXYZ(xyz, res=res(x), projection(x), digits=0)




# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

library(raster)

#projection
prj <- "+proj=longlat +ellps=WGS84 +datum=WGS84"

# list that
l <- list.files("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/clipped_buffer_IDW_shoreline", pattern=".tif",full.names=T)

# stack that 
ts.test <- stack(l)

# convert the stack to a data frame with xy
ts.table <- data.frame(coordinates(ts.test), getValues(ts.test))

# flip it into a pacific centered latlong
# flip the coordinates to pacific centered latlong
if(any(ts.table[,1] < 0))  ts.table[ts.table[,1]<0,1] <- ts.table[ts.table[,1]<0,1]+360


# rasterize it to a stack
# make a stack placeholder
s <- stack()

# now loop through the data in the table and rasterize that 
for(i in 3:ncol(ts.table)){

	xyz <- data.frame(ts.table[,1:2], ts.table[,i])

	r <- rasterFromXYZ(xyz, res=res(ts.test)[1], crs=prj, digits=0)
	s <- addLayer(s, r)
}


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

esa.stack <- stack()

for(i in 1:nlayers(tas.stack)){
	tas <- subset(tas.stack,i)
	vap <- subset(vap.stack,i)
	tmat <- as.matrix(tas)
	vmat <- as.matrix(vap)

	esa = 6.112*exp(((22.46*tmat) / (272.62 + tmat)))
	hur = (vmat/esa)*100

	rst <- raster(hur, extent(tas), crs=projection(tas))

	esa.stack <- addLayer(esa.stack, rst)



}

esa = 6.112*exp((17.62*tas1)/(243.12+tas1))

esa = 6.112*exp(22.46*tas1 / (272.62 + tas1))
hur = (vap/esa)*100


fun1<-function(x,y,filename=''){

	out <- brick(x, layers=nlayers(x))
	
	out <- writeStart(out, filename, overwrite=TRUE)
	
	bs <- blockSize(x)
	pb <- pbCreate(bs$n, progress)
	
	for (i in 1:bs$n) {
		v <- getValuesBlock(x, row=bs$row[i], nrows=bs$nrows[i])
		v2 <- getValuesBlock(y, row=bs$row[i], nrows=bs$nrows[i])

		v <-  6.112*exp(((22.46*v) / (272.62 + v)))
		v <- (v2/v)*100

		setValues(out, v, bs$row[i])
	
		pbStep(pb, i)
	}
	pbClose(pb)
	return(out)
	out <- writeStop(out)
}







dbf <- read.dbf("/big_storage/malindgren/AIEM/RHUM/hur_cru_10min_1961_1990/testing_area/hur_cru_10min_01_1961_1990.dbf")

lon <- dbf[,1]

# look for all of the data that are less than 0
ind0 <- which(dbf[,1] <= 0)
lon[ind0] <- dbf[ind0,1] + 360



# now the stuff on the other end of the scale
ind1 <- which(dbf[,1] > 0)
lon[ind1] <- abs(dbf[ind1,1] +180)


vap1.test <- cbind(coordinates(vap1), getValues(vap1))


vap1.test[,1] <- vap1.test[,1] + 180

vap1.flip <- rasterFromXYZ(vap1.test, res=res(vap1),digits=0)

rot <- rotate(vap1.flip)




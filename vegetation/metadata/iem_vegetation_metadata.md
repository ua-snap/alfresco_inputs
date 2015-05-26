Data record on SNAP Data portal:
http://ckan.snap.uaf.edu/dataset/land-cover-v0-4  

This 1km land cover dataset represent highly modified output originating from the Alaska portion of the North American Land Change Monitoring System (NALCMS) 2005 dataset as well as the National Land Cover Dataset 2001.  This model input dataset was developed solely for use in the ALFRESCO, TEM, GIPL and the combined Integrated Ecosystem Model landscape scale modeling studies and is not representative of any ground based observations.

**Original landcover data, including legends:**  
[NALCMS data](http://www.cec.org/Page.asp?PageID=924&ContentID=2819&AA_SiteLanguageID=1)  
[NLCD 2001](http://www.mrlc.gov/nlcd01_data.php)

------------
###Final Legend:

| **value** | **class name** |
|---|---|
0 | Not Modeled  
1 | Black Spruce Forest  
2 | White Spruce Forest  
3 | Deciduous Forest  
4 | Shrub Tundra  
5 | Graminoid Tundra  
6 | Wetland Tundra  
7 | Barren lichen-moss  
8 | Heath  
9 | Maritime Upland Forest  
10 | Maritime Forested Wetland  
11 | Maritime Fen  
12 | Maritime Alder Shrubland**

###Methods of production:
Due to specific models' land cover input requirements, including the fact that each model is primarily focused on different descriptive aspects of land cover (i.e. ALFRESCO considers land cover in respect to how it burns, TEM considers land cover in respect to how it cycles carbon through the system, and GIPL considers land cover with respect to its influence on the insulative qualities of the soil), 4 distinct subregions were delineated with unique sets of reclassification rules. Those regions are generally named below and also list the Unified Ecoregions of Alaska: 2001 (Nowacki et al 2001) Commoner names that make up each subregion.

1) **Kodiak Island**  
2) **Mountainous South Central Alaska and Canada** - Chugach-St Elias Mountains  
3) **Maritime influenced Southeast and South Central Alaska and Canada** - The Alexander Archipelago, Gulf of Alaska Coast, Boundary Range



**code for 1,2,3:** [landcarbon_vegetation_maritime.py](https://github.com/EarthScientist/alfresco_inputs/blob/master/vegetation/landcarbon_vegetation_maritime.py)

4) **Mainland Alaska and Canada** - all remaining ecoregions excluding the Aleutian Islands  

**code for 4:** [alfresco_vegetation.py](https://github.com/EarthScientist/alfresco_inputs/blob/master/vegetation/alfresco_vegetation.py)

**Note**, *legend values in these code bases do not represent final values in legend above.  
Input data were processed in their native resolutions and only resampled to 1km in the final step where they were mosiacked together.*

### Kodiak Island reclassification
This reclassification used the National Land Cover Database (NLCD) 2001 for Alaska as its starting point and was a straight 1 to 1 or many to 1 reclassification approach.The Kodiak Island Existing Vegetation Type, circa 2000 map for the Kodiak Archipelago produced by Michael Fleming (http://akevt.gina.alaska.edu/) was used to help inform the NLCD reclassification into the IEM land cover types.

**NLCD Classes** |**NLCD descriptions** | **Final IEM class value**  
--- | --- | ---  
11, 12, 22, 23, 31 | open water, perennial snow/ice, developed low, developed med, barren | 0 - Not Modeled  
41 | deciduous forest | 3 - Deciduous  Forest
42 | evergreen forest | 9 - Maritime Upland Forest  
51, 90 | dwarf shrub, woody wetland | 4 - Shrub Tundra  
52 | shrub/scrub | 12 - Maritime Alder Shrubland  
71 | grassland/herbaceous | 5 - Graminoid Tundra  
72,95 | sedge/herbaceous, emergent herbaceous wetland | 6 - Wetland Tundra  

### Mountainous South Central Alaska and Canada
This reclassification used the NLCD 2001 dataset as its starting point and was a straight 1 to 1 or many to 1 reclassification approach. .  Existing Vegetation Type, circa 2000 map for the region around Chitina produced by Michael Fleming (http://akevt.gina.alaska.edu/) was used to help develop the NLCD reclassification.

**NLCD Classes** |**NLCD descriptions** | **Final IEM class value**  
--- | --- | ---  
0, 11, 12, 21, 22, 23, 24, 31 | nodata, open water, perennial snow/ice, developed open space, developed low, developed med,developed high, barren | 0 - Not Modeled
41 | deciduous forest | 12 - Maritime Alder Shrubland  
42 | evergreen forest | 2 - White Spruce Forest  
43 | mixed forest | 3 - Deciduous  Forest
51,52,90 | dwarf shrub, shrub/scrub, woody wetland | 4 - Shrub Tundra  
71,72 | grassland/herbaceous,sedge/herbaceous | 8 - Heath  
95 | emergent herbaceous wetland | 6 - Wetland Tundra

### Maritime influenced Southeast and South Central Alaska and Canada
This reclassification used the NLCD 2001 dataset as its starting point for a 1 to 1 or many to 1 reclassification, but also utilized ancillary data including NLCD percent canopy cover (above and below 20% cover) and historical timber harvest regions for those areas within the Tongas National Forest. The percent canopy cover was used to distinguish between maritime *upland forest* vs maritime *fen* as well as between maritime *forested wetland* vs maritime *fen*. The timber harvest dataset covers the region between Yakataga and Dixon Entrance.  Because young-growth forests are not easily discernable in the NLCD and often misclassified, timber harvest extents were used to identify areas of logged maritme upland forest.

The West Copper River Delta Landscape Assessment (2003, updated 2007) was used to help develop the NLCD reclassification for the Maritime South Central region. The USDA Forest Service Tongass National Forest’s cover type map was used to help develop the NLCD reclassification for the Maritime Southeast region.

Details of the reclassification methods can be seen in the supplied code linked above. It is generally described below.

**NLCD Classes** |**NLCD descriptions** | **Final IEM class value**  
--- | --- | ---  
0,11,12,21,22,23,24,31 | no data, open water, perennial snow/ice, developed open space, developed low, developed med,developed high, barren | 0 - Not Modeled
42 | evergreen forest | 9 - Maritime Upland Forest or 11 - Maritime Fen  
41, 43, 51, 52, 71, 72, 90, 95 | deciduous,mixed forest,dwarf shrub,shrub/scrub,grassland/herbaceous,sedge/herbaceous,woody wetland,emergent herbaceous wetland | 10 - Maritime Forested Wetland or 11 - Maritime Fen
81,82 | pasture/hay, cultivated crops  | 12 - Maritime Alder Shrubland

**Then parse by canopy cover:**

 - if originally 42 and now 9 or 11 and canopy cover >20%, then = 9 - Maritime Upland Forest
 - if originally 42 and now 9 or 11 and canopy cover <=20%, then = 11 - Maritime Fen  
 - if originally 41, 43, 51, 52, 71, 72, 90, 95 and now 10 or 11 and canopy cover >20%, then = 10 - Maritime Forested Wetland  
 - if originally 41, 43, 51, 52, 71, 72, 90, 95 and now 10 or 11 and canopy cover <=20%, then = 11 - Maritime Fen  

**Then correct misclassified areas that have been logged in the past:**
if the pixels fall within the harvested masked area, then = 9 - Maritime Upland Forest


### Mainland Alaska and Canada
___
This reclassification used the NALCMS 2005 dataset as its starting point for a 1 to 1 or many to 1 reclassification, but also utilized ancillary data including:

**north south** - this layer is derived from an aspect calculation (in degrees) of the PRISM 2km Digital Elevation Model (DEM). It constitutes a reclassification into 3 classes: North, South, Water.  The classification scheme used to reclassify the aspect map is as follows: south = greater than 90 degrees and less than 301 degrees (value=1). north is all else (value=2), water (value=999). It is important to note here that when calculating aspect there are situations that arise where there is no slope and therefore no aspect.  In these situations the flat areas were reclassified as NORTH (2) unless they are a body of water.  Water bodies were extracted from the NALCMS Landcover Data and those areas that were both flat and waterbodies were reclassed as water (999).

**growing season temperature** - This layer  constitutes the average of the months of May, June, July, August, from the 1961-1990 PRISM climatology, resampled to 1km using a bilinear interpolation.

**coastal interior** - this layer was derived through combining features of the Nowacki Unified Ecoregions in Alaska (Nowacki, Gregory; Spencer, Page; Fleming, Michael; Brock, Terry; and Jorgenson, Torre. Ecoregions of Alaska: 2001.), and the Ecozones and Ecoregions of Canada (A National Ecological Framework for Canada: Attribute Data. Marshall, I.B., Schut, P.H., and Ballard, M. 1999. Agriculture and Agri-Food Canada, Research Branch, Centre for Land and Biological Resources Research, and Environment Canada, State of the Environment Directorate, Ecozone Analysis Branch, Ottawa/Hull.).

**Unified Ecoregions of Alaska** features include 2 classes of the "LEVEL_2" ecoregions in that dataset where Intermontane Boreal was classed as 2 and all other ecoregions (all coastal) was classed as 1.  This gives a differentiation between the wetlands classes that fall into each of these two reclassified regions, allowing for the creation of an Spruce Bog class from the "wetland" class in the NALCMS 2005 Land Cover Classification.  

**Ecozones/Ecoregions of Canada**  
access this ecozone/ecoregion shapefile data: [Canada Ecozones/Ecoregions](http://sis.agr.gc.ca/cansis/nsdb/ecostrat/gis_data.html)

With the area of interest including sites in Western Canada, it was important to differentiate the coastal and inland boreal ecozones/regions in order to mimick what was done on the Alaska side.  To do this the Canada Ecozone/ecoregion data involved the following:

i. In the areas around south Saskatchewan there are large areas of prairie, which also coincides with the bread basket of Canada since their green revolution.  Since this area is not boreal forest and is a heavily human-will dictated environment (since it is mainly farms) it is not a good candidate to include in this classification.  **removed the Canada Ecozone "Prairie" and Canada Ecoregion "Interlake Plain" & "Boreal Transition"; both of which exist to the North of the excluded "Prairie" Canda Ecozone.  This was determined to be removed by examining the input NALCMS Land Cover map and inspecting visually that there were little to no trees in these areas.

ii. To extend the "coastal" region beyond the southern extent of southeast Alaska, it was determined that the Canada Ecozone "Pacific Maritime" should be classed as coast to differentiate between coastal and non-coastal wetlands. Therefore this was added to the non-"Intermontane Boreal" classes from the Unified Ecoregions of Alaska map.

iii. The southern-most extent of the new coastal vs spruce bog layer is the international border between Canada and the U.S.

iv. The areas classified as "Boreal" on the Canada side include Ecozones of: Montane Cordillera, Boreal Cordillera, Taiga Cordillera, Taiga Plain, Boreal Plain, Taiga Shield, Boreal Shield, Hudson Plain, Mixed Wood Plain, Atlantique Maritime, Hudson Plain, Arctic Cordillera.

v. Included Areas to the North of the new "boreal" class on the Canada side include the Canada Ecozones of: Southern Arctic AND the Canada Ecoregions of: Wager Bay Plateau, Boothia Peninsula Plateau, Meta Incognita Peninsula, Central Ungava Peninsula, Foxe Basin Plain, Melville Peninsula Plateau, Baffin Island Uplands. *** Areas North of these locations are not considered for this layer.  

**treeline** - This layer was created by rasterizing the Circumpolar Arctic Vegetation Map [CAVM](http://www.geobotany.uaf.edu/cavm/) and defining the treeline using the extent of this data.  The result is a boolean map where 0=no trees and 1=trees.

**North Pacific Maritime** -  This layer was created by rasterizing the Unified Ecoregions of Alaska data for the Northern Pacific Rainforest region and creating a layer which reclassifies this region as the North Pacific Maritime Region.


The steps of reclassification of the NALCMS Land Cover Map are:

1. Once the area of interest (AOI) was determined the input map was clipped to that extent.  The resulting classes in the AOI are:

| **value** | **NALCMS type** | 
|---|---|
0 | out of bounds 
1 | Temperate or sub-polar needleleaf forest
2 | Sub-polar taiga needleleaf forest
5 | Temperate or sub-polar broadleaf deciduous
6 | Mixed Forest
8 | Temperate or sub-polar shrubland
10 | Temperate or sub-polar grassland
11 | Sub-polar or polar shrubland-lichen-moss
12 | Sub-polar or polar grassland-lichen-moss 
13 | Sub-polar or polar barren-lichen-moss
14 | Wetland
15 | Cropland
16 | Barren Lands
17 | Urban and Built-up
18 | Water
19 | Snow and Ice

    These initial classes were reclassed as follows:
    * classes [ 15, 17, 18, 19 ] were reclassed as no vegetation.
    * classes [ 1, 2 ] were reclassed as a spruce class.
    * classes [ 5, 6 ] were reclassed as deciduous forest
    * class [ 11 ] was reclassed as shrub tundra
    * The remaining initial classes were left as-is for further classification using ancillary layers.

2. The wetland class was reclassified using the coastal_interior layer into coastal wetlands and interior spruce bogs (spruce class).

3. The newly derived coastal wetland layer was further reclassified into wetland tundra or no veg by using the growing_season_temperature layer.  Any coastal wetland pixel with a growing season temperature <6.5 C is classified as wetland tundra and the remainder of the data are classified as no veg.

4.  Next the class 8 (Temperate or sub-polar shrubland) is reclassified into deciduous forest or shrub tundra using the growing_season_temperature layer.  All pixels with a growing season temperature <6.5 C is classified as shrub tundra and the remainder are classified as deciduous forest.

5. The classes 12 "Sub-polar or polar grassland-lichen-moss" and 10 "Temperate or sub-polar grassland" are reclassified into graminoid tundra or grassland based on the growing_season_temperature layer.  Where pixels of class 10 or 12 with growing season temperature values of <6.5 are classified as graminoid tundra and the values >6.5 are classified as grassland.

6. Then, we reclassify the spruce class we created in step 1 into black or white spruce using the north_south layer.  If a spruce pixel is north facing it is classified as black spruce, and if the pixel is south facing then it is classified as white spruce.

7. Due to some deficiencies in the NALCMS map there are spruce trees on the North Slope of Alaska, which is known to not be factual.  Instead of simply reclassifying these data into a single class, we found that it would be better to run a focal analysis on the data that would convert the suspect cells into most common class in its surrounding 16 neighbors.  This allowed for a more realistic reclassification of these suspect pixels with the land cover classes that are in the surrounding area.  

8. Using the north_pacific_maritime layer, the pixels that are within that extent are reclassified as temperate rainforest in the output map.

9. Class 13 (Sub-polar or polar barren-lichen-moss) was reclassified into class 7 "Barren lichen-moss"

###end of Mainland Alaska and Canada process

**All 4 regions were then mosaicked together at 1km spatial resolution into a single final map.**

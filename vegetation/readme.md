# Vegetation

Set of scripts for generating the needed inputs for alfresco and the iem 
project as a whole.  

- **alfresco_vegetation.py** - reclassifies NALCMS 2005 LandCover over Alaska and Western Canada using ancillary layers to best classify the landscape for a model input to the ALFRESCO Fire Model.

- **landcarbon_vegetation_maritime.py** - reclassifies NLCD 2001 LandCover over Alaska using ancillary layers to best classify the landscape for use as model input to the LandCarbon Project.  This layer has also subsequently been added to the above alfresco vegetation map to fill-in the maritime regions that were not previously classified.

- **landcarbon_resample_mosaic.py** - resamples and mosaic's the two above generated maps ( alfresco_vegetation & landcarbon_vegetation_maritime ) into a single *somewhat* harmonized raster.

[ Metadata ]:https://github.com/EarthScientist/alfresco_inputs/blob/master/vegetation/metadata/iem_vegetation_metadata.md

[ Metadata ]

![alt text](https://github.com/EarthScientist/alfresco_inputs/blob/master/vegetation/metadata/iem_vegetation_model_input_v0_4_map.png "IEM Vegetation Model Input Map Composition")

		[ More Documentation to follow ]
# download and unzip the fire history perimeter file GDB from the AICC
import os, glob

# setup args
output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/raw_downloaded'
link = 'https://fire.ak.blm.gov/content/maps/aicc/Data%20Downloads/Data%20(zipped%20filegeodatabases)/FireHistory_Perimeters_1940_2017.zip'

# download the data...
os.chdir( output_dir )
os.system( 'wget "{}"'.format( link ) )

fn, = glob.glob( os.path.join( output_dir, '*.zip' ) )
os.system( 'unzip {}'.format( fn ) )
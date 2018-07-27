# download and unzip the fire history perimeter file GDB from the AICC
import os, glob

# setup args
output_dir = '/workspace/Shared/Tech_Projects/ALFRESCO_Inputs/project_data/Fire/FireHistory/raw_downloaded'
link = 'http://cwfis.cfs.nrcan.gc.ca/downloads/nfdb/fire_poly/current_version/NFDB_poly.zip'

# download the data...
os.chdir( output_dir )
os.system( 'wget "{}"'.format( link ) )
# unzip
os.system( 'unzip {}'.format( os.path.basename(link) ) )
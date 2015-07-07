# pre-processing area
def group_input_filenames( prefix, root_dir ):
	import fnmatch, functools, itertools, os, glob
	import pandas as pd
	''' function that wraps some ugliness regarding returning the files we want to process '''
	def find_files( dir_path, patterns ):
		"""
		Returns a generator yielding files matching the given patterns
		:type dir_path: str
		:type patterns: [str]
		:rtype : [str]
		:param dir_path: Directory to search for files/directories under. Defaults to current dir.
		:param patterns: Patterns of files to search for. Defaults to ["*"]. Example: ["*.json", "*.xml"]
		"""
		import itertools, functools
		path = dir_path
		if not patterns:
			path_patterns = [ "*" ]
		else:
			path_patterns = patterns

		for root_dir, dir_names, file_names in os.walk( path ):
			filter_partial = functools.partial(fnmatch.filter, file_names)

			for file_name in itertools.chain( *map( filter_partial, path_patterns ) ):
				yield os.path.join( root_dir, file_name )
	def version_grouper( x ):
		''' groupby function for grouping by filenames '''
		dir_path = os.path.dirname( x )
		fn, _ = os.path.splitext( os.path.basename( x ) )
		# remove dates from filename -- they have a hyphen
		fn_base = '_'.join([ i for i in fn.split( '_' ) if '-' not in i ])
		# return the path element that startswith 'v' this is the version attribute
		version = [ x for x in dir_path.split( os.path.sep ) if x.startswith( 'v' ) ]
		return '_'.join([ fn_base, version[0] ])
	def drop_old_versions( df ):
		rows,cols = df.shape
		if rows > 1 & rows < 3:
			version_nums = df[ df.columns[-1] ].apply( lambda x : int( x.replace( 'v', '' ) ) )
			# max( version_nums )
			return df.drop( df[ df[ df.columns[-1] ] != 'v' + str( max( version_nums ) )].index )
		elif rows > 3:
			# potentially unnecessary
			None
		else:
			return df

	# get all matches with prefix
	matches = pd.Series([ match for match in find_files( root_dir, [ prefix ] ) ])
	# group by version
	grouped = dict([ group for group in matches.groupby( matches.apply( version_grouper ) )])
	# group keys to DataFrame
	keys_df = pd.DataFrame({ key:key.split( '_' ) for key in grouped.keys() }).T
	# parse the keys / values and keep only latest versions
	keys_df_grouped = pd.concat([ drop_old_versions(i[1]) for i in keys_df.groupby( keys_df.columns[-3] ) ])
	# make a new dictionary holding the filenames grouped the way we want
	final_out = { k:v for k,v in grouped.iteritems() if k in keys_df_grouped.index.tolist() }
	return final_out
def get_file_years( filename ):
	path, fn = os.path.split( filename )
	fn, ext = os.path.splitext( fn )
	split = fn.split( '_' )
	dates = split[ len( split ) - 1 ] # grab last element
	begin, end = dates.split( '-' )
	return [begin, end]
def get_modelname( filename ):
	path, fn = os.path.split( filename )
	return [ i for i in path.split( '/' ) if i in models ][0]
def concat_to_nc( filelist, output_filename, dim='time', begin_time=None, end_time=None, nc_format='NETCDF4', **kwargs ):
	'''
	take list of consecutive netcdf files (made for CMIP5 data) and stack them into a 
	single larger netcdf file.  This was necessary to overcome some bugginess in how 
	MFDataset is dealing with different calendar units on different files.  This is 
	technically valid CF-Compliant metadata, but is tricky to work with.  This hack allows
	us to get around some of this unpredictable behavior.

	PARAMETERS:
	-----------
	filelist = [list] list of string file paths to the sorted netcdf files to stack together
	output_filename = [str] path to and name of the output file to be generated (.nc extension)
	dim = [str] dimension to stack on -- default is 'time'
	begin_time = [str] PANDAS style datetime string syntax -- used in xray
	end_time = [str] PANDAS style datetime string syntax -- used in xray
	format = [str] output NetCDF format desired. valid strings are:
					'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_64BIT', 'NETCDF3_CLASSIC'
					default is 'NETCDF4'
	**kwargs -- potential future arguments or overloaded args to pass through (none implemented)

	RETURNS:
	--------

	output_filename as string, with the important side-effect of writing data to disk

	'''
	import xray
	with xray.concat([ xray.open_dataset( i ).load() for i in filelist ], dim ) as ds:
		# time slicer condition
		if begin_time != None and end_time != None:
			ds = ds.loc[ { dim:slice( begin_time, end_time ) } ]
		if os.path.exists( output_filename ):
			os.remove( output_filename )
		ds.to_netcdf( output_filename, mode='w', format=nc_format )
	return output_filename

# # hand processing of a file that exceeds PANDAS DateRange functionality > 2300-07-06 00:00:00, which 
# # I am going to put a bug report in about
# fn = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data/cmip5/output1/IPSL/IPSL-CM5A-LR/rcp26/mon/atmos/Amon/r1i1p1/v20120114/tas/tas_Amon_IPSL-CM5A-LR_rcp26_r1i1p1_200601-230012.nc'
# ds = xray.open_dataset( fn )
# fn_split = os.path.splitext( os.path.basename( fn ) )[0].split( '_' )
# begin_time, end_time = fn_split[ len( fn_split ) - 1 ].split( '-' )

# begin_idx = 0
# end_idx = np.repeat( range( 2006, 2300+1), 12 ).tolist()[(int(end_year[:4]) - int(end_time[:4]))]

# variable = [str] abbreviated name of the variable being converted. i.e 'tas'/'pr'/'hur'
if __name__ == '__main__':
	import os, sys, re, glob, xray
	import numpy as np
	import pandas as pd

	# models = [ 'GISS-E2-R','IPSL-CM5A-LR', 'MRI-CGCM3', 'CCSM4', 'GFDL-CM3' ]
	models = [ 'CCSM4' ]
	variables = [ 'tas', 'hur' ]
	base_path = '/workspace/Shared/Tech_Projects/ESGF_Data_Access/project_data/data'

	# problem_files
	output_base_path = os.path.join( base_path, 'prepped' )
	if not os.path.exists( output_base_path ):
		os.makedirs( output_base_path )

	problem_files_log = open( os.path.join( output_base_path, 'error_files.txt' ), mode='w' )

	for model in models:
		for variable in variables:
			fn_prefix_filter = variable + '_*' + model + '*'
			file_groups = group_input_filenames( fn_prefix_filter, os.path.join( base_path, 'cmip5' ) )
			for files in file_groups.values():
				try:
					files = sorted( files.tolist() )
					output_path = os.path.join( output_base_path, model, variable )

					if not os.path.exists( os.path.dirname( output_path ) ):
						os.makedirs( os.path.dirname( output_path ) )

					fn = files[ 0 ]
					model = get_modelname( fn )
					begin_year ='190001'
					end_year = '210012'

					# a handler for the historical (1900-2005) and the modeled (2006-2100) filebnameing
					if '_historical_' in os.path.basename( fn ):
						begin_year_fnout = '190001'
						end_year_fnout = '200512'
					else:
						begin_year_fnout = '200601'
						end_year_fnout = '210012'

					# this is a hacky sort of thing... but we need a way to output the [proper] naming convention
					output_filename = os.path.join( output_path, '_'.join([ '_'.join( os.path.splitext( os.path.basename( fn ) )[0].split( '_' )[:-1] ), begin_year_fnout, end_year_fnout + '.nc' )

					# this logic can be fine tuned to subset the data down to only the files we need
					# for this project it is 1900-2100.
					df = pd.DataFrame([ get_file_years(fn) for fn in files ])

					# this is the way to interrogate that dataframe for the values we want
					df = df.astype( int )
					begin_idx = (np.abs(df[0].astype( int ) - int( begin_year ) ) ).argmin()
					end_idx = (np.abs(df[1].astype( int ) - int( end_year ) ) ).argmin()

					# return the files between the desired date ranges
					if begin_idx == end_idx:
						files = [files[ begin_idx ]]
					else:
						files = files[ begin_idx:end_idx + 1 ]

					print files
					print '\n'
					
					# run the concatenation and the output to a new netcdf file
					concat_to_nc( files, output_filename, dim='time', begin_time=begin_year[:4], end_time=end_year[:4] )

				except:
					print '\n--> ERROR !!!\n\n%s\n\n' % files
					problem_files_log.writelines( files )
					pass

	problem_files_log.close()


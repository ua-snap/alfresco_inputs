def replace_erroneous_treeline( lc2d, tl2d ):
	''' replace values based on neighbors of offenders and a where condition 
		arguments:
			lc2d = 2d landcover array with erroneous values within0
			tl2d = 2d treeline boolean array to subset the landcover array
		returns:
			2d numpy array with the erroneous values removed
	'''
	out = []
	ind = np.where( (((lc2d == 1) | (lc2d == 2)) & (tl2d == 1)) ) # this is hardwired...
	ind_zip = zip( *ind )
	print len( ind_zip )
	if len( ind_zip ) == 0:
		return lc2d
	else:
		index_groups = [ [( i-2, j-2 ),
						( i-2, j-1 ),
						( i-2, j+0 ),
						( i-2, j+1 ),
						( i-2, j+2 ),
						( i-1, j-2 ),
						( i-1, j-1 ),
						( i-1, j+0 ),
						( i-1, j+1 ),
						( i-1, j+2 ),
						( i+0, j-2 ),
						( i+0, j-1 ),
						( i+0, j+1 ),
						( i+0, j+2 ),
						( i+1, j-2 ),
						( i+1, j-1 ),
						( i+1, j+0 ),
						( i+1, j+1 ),
						( i+1, j+2 ),
						( i+2, j-2 ),
						( i+2, j-1 ),
						( i+2, j+0 ),
						( i+2, j+1 ),
						( i+2, j+2 ) ]
						for i,j in ind_zip ]
		
		for count, group in enumerate( index_groups ):
			cols = np.array( [j for i,j in group] )
			rows = np.array( [i for i,j in group] )
			vals = lc2d[ ( rows, cols ) ]
			vals = vals[ (vals > 0) & (vals!=1) & (vals!=2) ]
			out.append( vals )
	return out




def replace_erroneous_treeline( lc2d, tl2d ):
	''' replace values based on neighbors of offenders and a where condition 
		arguments:
			lc2d = 2d landcover array with erroneous values within
			tl2d = 2d treeline boolean array to subset the landcover array
		returns:
			2d numpy array with the erroneous values removed

	'''
	ind = np.where( (((lc2d == 1) | (lc2d == 2)) & (tl2d == 1)) ) # this is hardwired...
	ind_zip = zip( *ind )
	print len( ind_zip )
	if len( ind_zip ) == 0:
		return lc2d
	else:
		index_groups = [ [(i-1,j-1), 
						  (i-1, j+0), 
						  (i-1,j+1), 
						  (i+0,j-1), 
						  (i+0,j+1), 
						  (i+1,j-1), 
						  (i+1,j+0), 
						  (i+1,j+1)]
						for i,j in ind_zip ]
		
		for count, group in enumerate( index_groups ):
			cols = np.array( [j for i,j in group] )
			rows = np.array( [i for i,j in group] )
			vals = lc2d[ ( rows, cols ) ]
			vals = vals[ (vals!=1) & (vals!=2) ]
			uniques, counts = np.unique( vals, return_counts = True )
			new_val = uniques[ np.argmax( counts ) ]
			
			if new_val == 0:
				try:
					vals = vals[ (vals > 0) ]
					uniques, counts = np.unique( vals, return_counts = True )
					new_val = uniques[ np.argmax( counts ) ]
				except:
					new_val = 0	
			lc2d[ ind_zip[ count ] ] = new_val
		return replace_erroneous_treeline( lc2d, tl2d )



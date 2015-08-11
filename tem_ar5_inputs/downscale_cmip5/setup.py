from setuptools import setup

dependencies_list = [ 'xray','rasterio','pandas','numpy','rasterio','pathos' ]

scripts_list = ['bin/downscaling_launcher.py',
				 'bin/hur_ar5_model_data_downscaling.py',
				 'bin/hur_ar5_model_data_preprocess.py']

classifiers = [
		# How mature is this project? Common values are
		#   3 - Alpha
		#   4 - Beta
		#   5 - Production/Stable
		'Development Status :: 3 - Alpha',

		'Intended Audience :: Users',
		'Topic :: Data Downscaling :: CMIP5 - Integrated Ecosystem Model Inputs',

		# Pick your license as you wish (should match "license" above)
		 'License :: OSI Approved :: MIT License',

		# Specify the Python versions you support here. In particular, ensure
		# that you indicate whether you support Python 2, Python 3 or both.
		'Programming Language :: Python :: 2.7'
		]

setup(	name='downscale_cmip5',
		version='2.0',
		description='tool to downscale CMIP5 model outputs for use in the Integrated Ecosystem Model Project',
		url='https://github.com/ua-snap/alfresco-calibration',
		author='Michael Lindgren',
		author_email='malindgren@alaska.edu',
		license='MIT',
		packages=['downscale_cmip5'],
		install_requires=dependencies_list,
		zip_safe=False,
		include_package_data=True,
		dependency_links=['https://github.com/uqfoundation/pathos'],
		scripts=scripts_list,
		classifiers=classifiers
	)


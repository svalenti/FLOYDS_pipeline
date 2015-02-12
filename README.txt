###########################################################################
# 
#                                   Floyds Pipeline
#
#                              INSTALLATION
###########################################################################


FLOYDS is written in python and requires the following package:

- Python 2.5 or Python 2.6 or Python 2.7
   these modules have to be installed:
            - numpy
	    - pyraf
	    - matplotlib
	    - pyfits

- Iraf
- xhtml2pdf 
	(pip install  --allow-external pyPdf --allow-unverified pyPdf \
	xhtml2pdf was needed to make it install with recent pip versions; your
	mileage may vary)
	    
##############################################################################
extract the files from the tarball
> tar -xvf floyds-version.tar

> cd floyds-version
> python setup.py install  (--record files.txt)  (--prefix=<install location>)

##########################################################################
To uninstall a previus version 

- delete the floyds directory in your site-package path
- delete the floyds****.egg-info from the same directory
- delete the floyd executable: floydsspec 

or if during installation  you used the option: --record files.txt
you can run the following command in theterminal:
> cat files.txt | xargs sudo rm -rf

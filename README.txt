###########################################################################
# 
#                                   Floyds Pipeline
#
#                              INSTALLATION
###########################################################################

FURTHER INFORMATION

Go to this document for more instructions on how to run the pipeline:

https://www.authorea.com/users/598/articles/6566/_show_article

And this website gives an overview of the pipeline steps:

https://lco.global/documentation/floyds-pipeline/

##########################################################################

FLOYDS is written in python and requires the following packages:

- Python 2.5 or Python 2.6 or Python 2.7
   these modules have to be installed:
            - numpy
	    - pyraf
	    - matplotlib
	    - astropy

- Iraf
- xhtml2pdf 
	(pip install  --allow-external pyPdf --allow-unverified pyPdf \
	xhtml2pdf was needed to make it install with recent pip versions; your
	mileage may vary)
	    
#########################################################################

There are two options for installation.  The first is a straight-up
installation given that you have the required packages described
above.  The second is the setup of a conda environment for a 'clean'
installation.  We describe both here.

#########################################################################

STRAIGHT-UP INSTALLATION

extract the files from the tarball  [*****right now there is *no* tarball, right?  The user must get the FLOYDS code from github*****]
> tar -xvf floyds-version.tar

> cd floyds-version
> python setup.py install  (--record files.txt)  (--prefix=<install location>)

At this point, you should be ready to go.  Go to this document for
more instructions on how to run the pipeline:

https://www.authorea.com/users/598/articles/6566/_show_article



****To uninstall a previous version 

- delete the floyds directory in your site-package path
- delete the floyds****.egg-info from the same directory
- delete the floyd executable: floydsspec 

or if during installation  you used the option: --record files.txt
you can run the following command in the terminal:
> cat files.txt | xargs sudo rm -rf


############################################################################

Another option for installation is the creation of a 'clean' conda.
Note that several commands below are for the bash environment.

Preferably in the directory where you have placed the FLOYDS pipeline.


> conda create -n floyds python=2 numpy astropy ipython matplotlib
> source activate floyds
> conda install -c pkgw pyraf iraf
> mkiraf
> export IRAFARCH=linux64
> export iraf=/home/dsand/.conda/envs/floyds/lib/iraf/
> python setup.py install

The FLOYDS pipeline should now be setup.  Anytime you want to run it,
you must activate the floyds conda environment, and run the two export
lines above (if you are using bash).

> source activate floyds
> export IRAFARCH=linux64
> export iraf=/home/dsand/.conda/envs/floyds/lib/iraf/




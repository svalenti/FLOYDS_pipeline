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
	    

############################################################################

We recommend a installation in something like anaconda python, in a
'clean' conda environment. Note that several commands in the sample
below are for the bash environment.  If you have all of the
dependencies mentioned above, you can begin from the 'git clone'
command below.


> conda create -n floyds python=2 numpy astropy ipython matplotlib
> source activate floyds
> conda install -c pkgw pyraf iraf
> mkiraf
> export IRAFARCH=linux64
> export iraf=/home/dsand/.conda/envs/floyds/lib/iraf/
> git clone https://github.com/lcogt/FLOYDS_pipeline.git
> cd FLOYDS_pipeline
> python setup.py install

The FLOYDS pipeline should now be setup.  Anytime you want to run it,
you must activate the floyds conda environment, and run the two export
lines above (if you are using bash).

> source activate floyds
> export IRAFARCH=linux64
> export iraf=/home/dsand/.conda/envs/floyds/lib/iraf/

If you want to update your version, 

> git pull origin master
> python setup.py install


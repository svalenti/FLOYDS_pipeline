from distutils.core import setup
from distutils.command.install import INSTALL_SCHEMES
from os import sys, path
import os,shutil,re
from glob import glob
for scheme in INSTALL_SCHEMES.values():
    scheme['data'] = scheme['purelib']

#  test
from imp import find_module
try: find_module('numpy')
except: sys.exit('### Error: python module numpy not found')

try:
    find_module('astropy')
except:
    sys.exit('### Error: python module astropy not found')

try:
    find_module('pyraf')
except:
    sys.exit('### Error: python module pyraf not found')

try:
    find_module('matplotlib')
except:
    sys.exit('### Error: python module matplotlib not found')

setup(
    name='floyds',
    version='2.6.0',
    author='S. Valenti',
    author_email='svalenti@lcogt.net',
    scripts=['bin/floydsspec','bin/floydsfixheaders','bin/floydsauto'],
    url='lcogt.net',
    license='LICENSE.txt', 
    description='floyds is a package to reduce floyds spectra',
    long_description=open('README.txt').read(),
    install_requires=['numpy==1.16.6','astropy==5.3.3','pyraf==2.1.15','matplotlib==2.2.4', 'xhtml2pdf==0.2.4', 'pathlib2==2.3.5', 'requests==2.22.0', 'stsci.tools==3.6.0'],
    packages=['floyds'],
    package_dir={'':'src','doc':'doc'},
    package_data = {'floyds' : ["standard/MAB/*","standard/ident/*","standard/cat/*","standard/extinction/*",\
                                "standard/ident/en*/*",
                                "standard/fits/*","standard/sex/*","standard/stdlist/*","standard/flux/*",\
                                "archive/ftn/*/arc/blu/*/*fits","archive/ftn/*/arc/blu/*/*/id*",\
                                "archive/ftn/*/arc/red/*/*fits","archive/ftn/*/arc/red/*/*/id*",\
                                "archive/ftn/*/atmo/red/*fits","archive/fts/*/atmo/red/*fits",\
                                "archive/ftn/flat/*/*fits","archive/ftn/*/sens/*/*fits","archive/ftn/*/arc/*/*fits",\
                                "archive/fts/*/arc/blu/*/*fits","archive/fts/*/arc/blu/*/*/id*",\
                                "archive/fts/*/arc/red/*/*fits","archive/fts/*/arc/red/*/*/id*",\
                                "doc/*pdf",\
                                "archive/fts/flat/*/*fits","archive/fts/*/sens/*/*fits","archive/fts/*/arc/*/*fits"]}
)

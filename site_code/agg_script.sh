date
export LD_LIBRARY_PATH=":/icc/bin/lib"
#export PYTHONPATH=$PYTHONPATH:/home/eng:/opt/epd/lib/python2.7/site-packages/pyraf-1.11-py2.7-linux-x86_64.egg:/opt/epd/lib/python2.7/site-packages/lcogt_dateutils-0.1.0-py2.7.egg:/opt/epd/lib/python27.zip:/opt/epd/lib/python2.7:/opt/epd/lib/python2.7/plat-linux2:/opt/epd/lib/python2.7/lib-tk:/opt/epd/lib/python2.7/lib-old:/opt/epd/lib/python2.7/lib-dynload:/opt/epd/lib/python2.7/site-packages:/opt/epd/lib/python2.7/site-packages/PIL
export TERM=xterm
cd /home/eng/
THISHOST=$(hostname)
./agg_floyds.py -s $THISHOST

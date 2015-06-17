#!/opt/epd/bin/python

import pyfits
import os
#import matplotlib
#matplotlib.use('Agg')
import scipy
from numpy import array, where, nan
from xml.dom import minidom
from plot_guideinfo import plot_guideflux, plot_guidepos, plot_xposwtime, plot_yposwtime, plot_fwhmwtime, plot_guidestate

def mk_guideinfo_plots(rootname, utstart, utstop, logname, debug=True):
    '''Generate plots from FITS guide images'''

    if debug: print 'mk_guide_plots: Got:', rootname, utstart, utstop, logname
    
    (timeimg_list, fileimg_list) = readlogfile(logname)
    timeimg = scipy.array(timeimg_list) 
    fileimg = scipy.array(fileimg_list) 

#pull out all lines from the guide file that happened between utstart and utstop

    time_cull = timeimg * (timeimg > float(utstart)) * (timeimg < float(utstop))
    timeimg_gd = timeimg[where(time_cull != 0)]
    fileimg_gd = fileimg[where(time_cull != 0)]

    if fileimg_gd.size <= 0:
        print 'No guide frames taken during exposure!'
    else:
        guidecheck = 1
        mjd_list = read_mjd_from_guide_fits(fileimg_gd)
        mjd_arr = array(mjd_list)
        mjd_arr = mjd_arr - (min(mjd_arr)-0.1/60.0/24.0)
        mjd_arr_sec = mjd_arr*24.0*60.0*60.0
        if debug >2: print mjd_list

        guidestate_list = read_guidestate_from_guide_fits(fileimg_gd)

        plot_guidestate(mjd_arr_sec,guidestate_list)
        guidestateplotname = rootname+'_guidestate.png'
        retval = os.system('mv yo.png '+guidestateplotname)
        if debug: print "mv from yo.png to",guidestateplotname,"had status",retval

#        xmlfiles_gd = [f.replace('.fits', '.fits.inst.guide.xml') for f in fileimg_gd]
        xmlfiles_gd = [f.replace('.fits', '.fits.guide.xml').replace('flash/', 'cat/') for f in fileimg_gd]
        if debug: print xmlfiles_gd
        (totcnts_gd, xcen_gd, ycen_gd, fwhm_gd) = read_stats_from_xml_files(xmlfiles_gd)

        plot_guideflux(mjd_arr_sec,totcnts_gd)
        guidecountplotname = rootname+'_guidecounts.png'
        retval = os.system('mv yo.png '+guidecountplotname)
        if debug: print "mv from yo.png to",guidecountplotname,"had status",retval

        plot_guidepos(xcen_gd,ycen_gd)
        guidexypos_plotname = rootname+'_guidexy.png'
        retval = os.system('mv yo_xy.png '+guidexypos_plotname)
        if debug: print "mv from yo_xy.png to",guidexypos_plotname,"had status",retval
        
        plot_xposwtime(mjd_arr_sec,xcen_gd)
        guidext_plotname = rootname+'_guidext.png'
        retval = os.system('mv yo_xt.png '+guidext_plotname)
        if debug: print "mv from yo_xt.png to",guidext_plotname,"had status",retval
        
        plot_yposwtime(mjd_arr_sec,ycen_gd)
        guideyt_plotname = rootname+'_guideyt.png'
        retval = os.system('mv yo_yt.png '+guideyt_plotname)
        if debug: print "mv from yo_yt.png to",guideyt_plotname,"had status",retval
        
        plot_fwhmwtime(mjd_arr_sec,fwhm_gd)
        guidefwhmt_plotname = rootname+'_guidefwhmt.png'
        retval = os.system('mv yo_fwhm.png '+guidefwhmt_plotname)
        if debug: print "mv from yo_fwhm.png to",guidefwhmt_plotname,"had status",retval

    return

def readlogfile(logfile):
    '''Read in the new-style guide log which consists of UTC datetime of guide 
    frame start (e.g. 20140411050644) and the path and filename of the guide
    image.
    The HHMMSS part of the datetime and the filenames are returned in 2 lists'''

    log_fh = open(logfile, 'r')
    timeimg = []
    fileimg = []
    for line in log_fh.readlines():
        chunks = line.split()
        if len(chunks) == 2:
            timestr = chunks[0]
            timestr = timestr[8:]
            timeimg.append(float(timestr))
            fileimg.append(chunks[1])
    log_fh.close()
    
    return [timeimg, fileimg]

def read_mjd_from_guide_fits(fits_images):
    '''Loops through the <fits_images> list of FITS guide frames, reading in 
    the MJD and returning this as a list of floats'''

    mjd_list = []    
    for _fileimg in fits_images:
        headlist = pyfits.open(_fileimg)
        headnow = headlist[0].header

        mjdnow = headnow['MJD-OBS']
        try:
            mjd_list.append(float(mjdnow))
        except NameError:
            mjd_list = [float(mjdnow)]

        headlist.close()

    return mjd_list


def read_guidestate_from_guide_fits(fits_images):
    '''Loops through the <fits_images> list of FITS guide frames, reading in 
    the AGSTATE and returning this as a list of strings'''

    guidestate_list = []    
    for _fileimg in fits_images:
        headlist = pyfits.open(_fileimg)
        headnow = headlist[0].header

        guidestatenow = headnow['AGSTATE']
        try:
            guidestate_list.append(guidestatenow)
        except NameError:
            guidestate_list = [guidestatenow]

        headlist.close()

    return guidestate_list

def read_stats_from_xml_files(xmlfiles):

    totcnts_gd = []
    xcen_gd = []
    ycen_gd = [] 
    fwhm_gd = []
    for x in xmlfiles:
        (counts, fwhm, x_center, y_center) = read_xml_guide_file(x)
        if None not in  [counts, fwhm, x_center, y_center]:
            totcnts_gd.append(counts)
            xcen_gd.append(x_center)
            ycen_gd.append(y_center)
            fwhm_gd.append(fwhm)
        else:
            print "No values found for", x

    return (totcnts_gd, xcen_gd, ycen_gd, fwhm_gd)

def getXMLvalues(element, tag, debug=False):
    '''Extract the value from the searched-for tag in the passed element and
    return the value.
    Note: the above are semirandom words which probably don't have meaning or
    understanding...'''
    
    thing = element.getElementsByTagName(tag)[0]
    thing_value = thing.childNodes[0].data
    if debug: print thing_value
    
    return thing_value
    

def read_xml_guide_file(guidefile, debug=False):

    guide_doc = minidom.parse(guidefile)
    centroids = guide_doc.getElementsByTagName("centroids")
    
    flux_value = nan
    fwhm_value = nan
    x_pixel = nan
    y_pixel = nan
    for centroid in centroids:
        guidestar = centroid.getElementsByTagName("guideStar")[0]
        gs_value = guidestar.childNodes[0]
        if debug: print gs_value
        if gs_value.data == 'true':
            if debug: print "Found guidestar", gs_value.nodeType, gs_value.nodeValue
            flux_value = getXMLvalues(centroid, "totalFlux")
            fwhm_value = getXMLvalues(centroid, "fwhm")

            guidepixel = centroid.getElementsByTagName("pixel")[0]
            x_pixel = getXMLvalues(guidepixel, "x")
            y_pixel = getXMLvalues(guidepixel, "y")

    return (flux_value, fwhm_value, x_pixel, y_pixel)

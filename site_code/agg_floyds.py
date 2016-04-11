#!/usr/bin/env python

import os
import matplotlib
matplotlib.use('Agg')

import pyfits
import numpy
import scipy
import pylab
import sys
import math
import logging
import glob
import string
# New style plotting
import plot_xmlguideinfo as plot_guideinfo
from pyfits import open as popen
from numpy import array, abs, where
from datetime import datetime, timedelta


def determine_camnames(site):
    '''Method to use the site name (only do a partial match to each 2m site) to
    determine what our acquisition/guide and spectrograph camera codes are. Halts
    if the passed site isn't a recognized 2m site.'''

    if 'ogg' in site:
        acq_cam = 'kb42'
        spec_cam = 'en06'
    elif 'coj' in site:
        acq_cam = 'kb37'
        spec_cam = 'en05'
    else:
        print 'Unrecognized site', site
        sys.exit(-1)

    return acq_cam, spec_cam


def agg_floyds(nightlist, site='floyds.coj.lco.gtn', tmp_dir="./", debug=False):
    """Set up basic logging config"""
    logfile = 'agglog.txt'
    path, logfile = os.path.split(logfile)
    logging.basicConfig(filename=tmp_dir + logfile, filemode="w", level=logging.DEBUG,
                        format='%(asctime)s %(levelname)s:%(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info('Log for agg_floyds -- mostly for debug purposes')

    #os.system('export LD_LIBRARY_PATH=":/icc/bin/lib"')

# split the nights up, comma delimited
    night_list = string.split(nightlist, sep=',')

    acq_cam, spec_cam = determine_camnames(site)
    if debug:
        print "Cameras are:", acq_cam, spec_cam

    acq_root = os.path.join('/mnt/data/daydirs', acq_cam + os.sep)
    spec_root = os.path.join('/mnt/data/daydirs', spec_cam + os.sep)

# aggregate data for each night
    i = 0
    for _night in night_list:

        # grab last four characters of _night
        _night_four = _night[4:]
        _night_year = _night[:4]
        dirnow = '/var/www/html/night_summary/' + _night
        dir_for_guideimgs = '/var/www/html/images/' + _night
        if not os.path.exists(dirnow):
            os.makedirs(dirnow)
        if not os.path.exists(dir_for_guideimgs):
            os.makedirs(dir_for_guideimgs)

#        acq_images = glob.glob('/icc2/tmp/[s,r]*'+_night+'*.fits')
        acq_glob = acq_root + _night + '/flash/*' + acq_cam + '*' + _night + '*01.fits'
#        print "acq_glob=",acq_glob
        acq_images = glob.glob(acq_glob)
        acq_images.sort()
#        spec_images = glob.glob('/icc/tmp/[g,f]*'+_night+'*.fits')
        spec_glob = spec_root + _night + '/flash/*' + \
            spec_cam + '*' + _night + '*02.fits'
        if debug > 2:
            print "spec_glob=", spec_glob
        spec_images = glob.glob(spec_glob)
        spec_images.sort()
        if (acq_images == []) & (spec_images == []):
            print 'No images take on this night'
            sys.exit(0)
# Make a guider log from the XML files

#        gxml_glob = acq_root+'/*'+acq_cam+'*'+_night+'*01.fits.inst.guide.xml'
        gxml_glob = acq_root + _night + '/cat/*' + \
            acq_cam + '*' + _night + '*01.fits.guide.xml'
        gxml_images = glob.glob(gxml_glob)
        gxml_images.sort()
        if len(gxml_images) > 0:
            if debug:
                print "No. of guide XML files=", len(gxml_images)
            guidetag = 1
        else:
            print 'No guide images found!!!'
            guidetag = 0


# which acquisition images have region files?  acquisition Log files?
# figure out which ones do and: 1) link them

        os.system('rm -f ' + dir_for_guideimgs + '/*.reg')
        os.system('rm -f ' + dir_for_guideimgs + '/*.log')
        os.system('rm -f ' + dir_for_guideimgs + '/*.cat')
        os.system('rm -f ' + dir_for_guideimgs + '/*.fits')
        os.system('rm -f ' + dir_for_guideimgs + '/*.jpg')

        # print acq_images
        home = os.path.expanduser('~')
        guidelog = os.path.join(home, _night + '_guide.log')
        guide_fh = open(guidelog, 'w')
        previous_blkuid = ''
        new_acq_images = []
        for _acqimage in acq_images:
            # print os.path.getsize(_acqimage)
            if os.path.getsize(_acqimage) < 100000:
                badname = _acqimage
                acq_images.remove(badname)

            hdulist = pyfits.open(_acqimage)
            # print _acqimage
            prihdr = hdulist[0].header

            blkuid_acq = str(prihdr['BLKUID'])
# Remove manually taken images (which have BLKUID='N/A') and ones where BLKUID
# has the same value as the previous frame (which indicates its a guide frame
# and not an acquisition image
#
            if blkuid_acq == 'N/A' or blkuid_acq == previous_blkuid:
                if debug:
                    print "Rejecting frame", _acqimage, blkuid_acq, previous_blkuid
#                acq_images.remove(_acqimage)
            else:
                if debug:
                    print "Keeping   frame", _acqimage, blkuid_acq, previous_blkuid
                new_acq_images.append(_acqimage)
            previous_blkuid = blkuid_acq

# Extract date/time of observation from guide frame header, strip out unwanted
# characters and write this into our new guidelog
# Only do this if we have a <foo>-fits.inst.guide.xml file

#            guide_xml_file = _acqimage.replace('.fits','.fits.inst.guide.xml')
            guide_xml_file = _acqimage.replace('.fits', '.fits.guide.xml')
            guide_xml_file = guide_xml_file.replace('flash/', 'cat/')
            guide_dateobs = prihdr['DATE-OBS']
            if guide_dateobs != '' and guide_dateobs != 'N/A' and os.path.exists(guide_xml_file):
                for char in ['T', '-', ':']:
                    if char in guide_dateobs:
                        guide_dateobs = guide_dateobs.replace(char, '')
                if '.' in guide_dateobs:
                    guide_dateobs = guide_dateobs.split('.')[0]
                print >> guide_fh, "%s %s" % (guide_dateobs, _acqimage)
            hdulist.close()

        guide_fh.close()
        acq_images = new_acq_images
        if debug:
            print "Trimmed acq image list=",  len(acq_images), "images"
        for _acqimage in acq_images:
            if debug:
                print _acqimage

        for _acqimage in acq_images:
            # for each acquisition images, link the file & the jpg and any region files, acquisition Log files,
            # and catalog files that might exist and link them too

            acqfile = _acqimage
            acqfile_dest = os.path.basename(acqfile)
            acqfile_destpath = os.path.join(dir_for_guideimgs, acqfile_dest)

            yo = _acqimage  # .replace(acq_root,'/lco/floyds/tmp/')
            jpgfile = yo.replace('.fits', '.jpg').replace(
                'flash/', 'flash/jpg/')
            jpgfile_dest = os.path.basename(yo)
            jpgfile_dest = jpgfile_dest.replace(
                '.fits', '.jpg').replace('flash/', 'flash/jpg/')
            jpgfile_destpath = os.path.join(dir_for_guideimgs, jpgfile_dest)

            yo = _acqimage  # .replace(acq_root,'/lco/floyds/tmp/')
            regfile = yo.replace('.fits', '.fits.regs')
            regfile_dest = os.path.basename(yo)
            regfile_dest = regfile_dest.replace('.fits', '.reg')
            regfile_destpath = os.path.join(dir_for_guideimgs, regfile_dest)

            logfile = yo.replace('.fits', '.log')

            catfile = yo.replace('.fits', '.fits.sex')
            catfile_dest = os.path.basename(yo)
            catfile_dest = catfile_dest.replace('.fits', '.cat')
            catfile_destpath = os.path.join(dir_for_guideimgs, catfile_dest)

            if os.path.exists(acqfile):
                link_cmd = 'cp ' + acqfile + ' ' + acqfile_destpath
                if debug >= 2:
                    print "acq link_cmd=", link_cmd
                os.system(link_cmd)
            if os.path.exists(jpgfile):
                link_cmd = 'cp ' + jpgfile + ' ' + jpgfile_destpath
                if debug >= 2:
                    print "jpg link_cmd=", link_cmd
                os.system(link_cmd)
            else:
                print "Acquisition jpg missing, creating"
                status = make_jpg_image(jpgfile, jpgfile_destpath)
                print "JPG creation status =", status
            if os.path.exists(regfile):
                link_cmd = 'ln -s ' + regfile + ' ' + regfile_destpath
                if debug >= 2:
                    print "reg link_cmd=", link_cmd
                os.system(link_cmd)
            if os.path.exists(logfile):
                os.system('ln -s ' + logfile + ' ' + dir_for_guideimgs + '/')
            if os.path.exists(catfile):
                link_cmd = 'ln -s ' + catfile + ' ' + catfile_destpath
                if debug >= 2:
                    print "cat link_cmd=", link_cmd
                os.system(link_cmd)

# read in each acquisition header
            hdulist = pyfits.open(_acqimage)
            # print _acqimage
            prihdr = hdulist[0].header

            propid_acq = prihdr['PROPID']
            try:
                propid_acqlist.append(propid_acq)
            except NameError:
                propid_acqlist = [propid_acq]

            MJD_acq = prihdr['MJD-OBS']
            try:
                MJD_acqlist.append(MJD_acq)
            except NameError:
                MJD_acqlist = [MJD_acq]

            grpuid_acq = prihdr['BLKUID']
            if grpuid_acq == 'N/A':
                grpuid_acq = -1
            grpuid_acq = int(grpuid_acq)
            try:
                grpuid_acqlist.append(grpuid_acq)
            except NameError:
                grpuid_acqlist = [grpuid_acq]

            grpnumob_acq = int(prihdr['MOLFRNUM'])
            try:
                grpnumob_acqlist.append(grpnumob_acq)
            except NameError:
                grpnumob_acqlist = [grpnumob_acq]

            groupid_acq = prihdr['OBJECT']

# Turn any brackets and other odd chars into underscores
            for badchar in ['(', ')']:
                groupid_acq = groupid_acq.replace(badchar, '')
            for badchar in [' ', "'", "/"]:
                groupid_acq = groupid_acq.replace(badchar, '_')

            try:
                grpid_acqlist.append(groupid_acq)
            except NameError:
                grpid_acqlist = [groupid_acq]

            acq_utstop = prihdr['UTSTOP']
            try:
                acq_utstoplist.append(acq_utstop)
                acq_utstop_nocolon.append(float(acq_utstop.replace(':', '')))
            except NameError:
                acq_utstoplist = [acq_utstop]
                acq_utstop_nocolon = [float(acq_utstop.replace(':', ''))]

            acq_start = prihdr['UTSTART']
            try:
                acq_utstartlist.append(acq_start)
                acq_utstart_nocolon.append(float(acq_start.replace(':', '')))
            except NameError:
                acq_utstartlist = [acq_start]
                acq_utstart_nocolon = [float(acq_start.replace(':', ''))]

            hdulist.close()

        if guidetag != 0:
            first_guideimage = find_first_guide(
                _night, acq_images, acq_utstart_nocolon, acq_utstop_nocolon, MJD_acqlist, site)
        else:
            first_guideimage = 'Null'
#        print "First guideimage=",first_guideimage
        print dir_for_guideimgs

        if first_guideimage != 'Null':
            for _first in first_guideimage:
                if _first != 'Null':
                    print 'ln -s ' + _first + ' ' + dir_for_guideimgs + '/'
                    os.system('ln -sf ' + _first + ' ' +
                              dir_for_guideimgs + '/')
                    _first_jpg = _first.replace('.fits', '.jpg')
                    if os.path.exists(_first_jpg):
                        print "Linking", _first_jpg
                        os.system('ln -sf ' + _first_jpg +
                                  ' ' + dir_for_guideimgs + '/')
                    else:
                        print "Making", _first_jpg
                        os.system('/usr/local/bin/gpp ' + _first + ' ' +
                                  _first_jpg + ' 2 /usr/local/bin/gpp_cfg.txt')
                        os.system('mv ' + _first_jpg + ' ' + dir_for_guideimgs)

        for _specimage in spec_images:

            # Copy into same directory as guideimages
            #
            if debug > 2:
                print "Specimage=", _specimage
            copy_cmd = 'cp ' + _specimage + ' ' + dir_for_guideimgs + '/'
            if debug > 2:
                print copy_cmd
            os.system(copy_cmd)
            _specimage_jpg = _specimage.replace(
                '02.fits', '01.jpg').replace('flash/', 'flash/jpg/')
            if debug > 2:
                print "Specjpg  =", _specimage_jpg
            _specimage_jpg_dest = os.path.basename(_specimage_jpg)
            _specimage_jpg_dest = _specimage_jpg_dest.replace(
                '01.jpg', '02.jpg')
            if debug > 2:
                print "Specjpgds=", _specimage_jpg_dest
            copy_cmd = 'cp ' + _specimage_jpg + ' ' + \
                dir_for_guideimgs + '/' + _specimage_jpg_dest
            if debug > 2:
                print copy_cmd
            os.system(copy_cmd)

# read in and collect header info
            hdulist = pyfits.open(_specimage)
            prihdr = hdulist[0].header
            propid_spec = prihdr['PROPID']
            try:
                propid_speclist.append(propid_spec)
            except NameError:
                propid_speclist = [propid_spec]

            grpuid_spec = prihdr['BLKUID']
            if grpuid_spec == 'N/A':
                grpuid_spec = -1
            grpuid_spec = int(grpuid_spec)
            try:
                grpuid_speclist.append(grpuid_spec)
            except NameError:
                grpuid_speclist = [grpuid_spec]

            obstype_spec = prihdr['OBSTYPE']
            try:
                obstype_speclist.append(obstype_spec)
            except NameError:
                obstype_speclist = [obstype_spec]

            utstart_spec = prihdr['UTSTART']
            try:
                utstart_speclist.append(utstart_spec)
                utstart_speclist_nocolon.append(
                    float(utstart_spec.replace(':', '')))
            except NameError:
                utstart_speclist = [utstart_spec]
                utstart_speclist_nocolon = [
                    float(utstart_spec.replace(':', ''))]

            exptime_spec = prihdr['EXPTIME']
            try:
                exptime_speclist.append(float(exptime_spec))

            except NameError:
                exptime_speclist = [float(exptime_spec)]


# Fudge time to take off end of expsoure time to work around UTSTOP!=UTSTART+EXPTIME
# FLOYDS bug
            end_fudge = 22

            utstop_spec = prihdr['UTSTOP']
            try:
                print "UTSTOP before=", utstop_spec
                utstop_dt = datetime.strptime(utstop_spec, '%H:%M:%S.%f')
                utstop_dt = utstop_dt - timedelta(seconds=end_fudge)
                utstop_spec = utstop_dt.strftime('%H:%M:%S.%f')
                print "UTSTOP after=", utstop_spec
                utstop_speclist.append(utstop_spec)
                utstop_speclist_nocolon.append(
                    float(utstop_spec.replace(':', '')))
            except NameError:
                utstop_speclist = [utstop_spec]
                utstop_speclist_nocolon = [float(utstop_spec.replace(':', ''))]
            hdulist.close()

        try:
            grpuid_speclist_arr = array(grpuid_speclist)
            spec_images_arr = array(spec_images)
            utstart_speclist_arr = array(utstart_speclist_nocolon)
            utstop_speclist_arr = array(utstop_speclist_nocolon)
            obstype_speclist_arr = array(obstype_speclist)
            exptime_speclist_arr = array(exptime_speclist)
        except UnboundLocalError:
            grpuid_speclist_arr = numpy.zeros(1)
            spec_images_arr = numpy.zeros(1)
            utstart_speclist_arr = numpy.zeros(1)
            utstop_speclist_arr = numpy.zeros(1)
            obstype_speclist_arr = numpy.zeros(1)
            exptime_speclist_arr = numpy.zeros(1)

        # print obstype_speclist_arr
        # print utstart_speclist_arr
        # print utstop_speclist_arr

    # sys.exit(0)
        i = 0
        for _acq in acq_images:
            # find those spectra with the same grpuid
            grpuid_curr = grpuid_acqlist[i]
            print _acq, grpuid_curr  # , type(grpuid_curr)
#	    print grpuid_speclist
#	    print grpuid_speclist_arr
            specs_forthis_acq = spec_images_arr[
                where(grpuid_speclist_arr == grpuid_curr)]
            specs_forthis_acqlist = specs_forthis_acq.tolist()
            print "specs_forthis_acqlist=", specs_forthis_acqlist

            science_cull = (grpuid_speclist_arr == grpuid_curr) * \
                (obstype_speclist_arr == 'SPECTRUM')

            science_index = where(science_cull == True)

            sciencespecs_forthis_acq = spec_images_arr[science_index]

            sciencespecs_forthis_acqlist = sciencespecs_forthis_acq.tolist()

            utstart_science = utstart_speclist_arr[science_index]
            utstop_science = utstop_speclist_arr[science_index]
            exptime_science = exptime_speclist_arr[science_index]
            j = 0
            for _utstart in utstart_science:

                sciencespecnow = sciencespecs_forthis_acq[j]
                junk1 = os.path.basename(sciencespecnow)
                guideplot_rootname = os.path.splitext(junk1)[0]
                # print guidetag
                if guidetag > 0:
                    # print guideplot_rootname
                    # print _utstart
                    # print utstop_science[j]
                    guidecheck = plot_guideinfo.mk_guideinfo_plots(
                        guideplot_rootname, _utstart, utstop_science[j], guidelog)
                else:
                    guidecheck = 0
                    pass
                j = j + 1

            mk_obs_website(acq_images[i], grpid_acqlist[i], propid_acqlist[i],
                           acq_utstart_nocolon[i], _night, first_guideimage[
                               i], dir_for_guideimgs,
                           specs_forthis_acqlist, sciencespecs_forthis_acq)
            i = i + 1


# read in each spectroscopy header
#    for _specimage in spec_images:
#        hdulist = pyfits.open(_specimage)
#        prihdr = hdulist[0].header

# for each acquisition image, look for acquisition log

def find_first_guide(night, acqimages, acq_utstart, acq_utstop, acq_mjd, site):

    from numpy import array, abs, where

# Determine acquisition camera name and therefore path
    acq_cam, spec_cam = determine_camnames(site)
    imgdir = os.path.join('/mnt/data/daydirs/', acq_cam,
                          night, 'flash' + os.sep)
#    print acq_cam, imgdir


# find full frame guide images from the night in question

    # print imgdir+nightnext+'*f.fits'
    # +glob.glob(imgdir+nightnext+'*f.fits') # no longer needed, same root
    full_guides_all = glob.glob(imgdir + '*' + night + '*g01.fits')
    print "No. of full guide images=", len(full_guides_all)
    acq_utstop_arr = array(acq_utstop)
    acq_mjd_arr = array(acq_mjd)

# get the UT start time for all of the full frame guide images
    for _guides in full_guides_all:
        hdulist = pyfits.open(_guides)
        prihdr = hdulist[0].header
        utstart_guides = prihdr['UTSTART']
        MJD_guides = prihdr['MJD-OBS']
        try:
            utstart_guidelist.append(float(utstart_guides.replace(':', '')))
        except NameError:
            utstart_guidelist = [float(utstart_guides.replace(':', ''))]
        try:
            MJD_guideslist.append(MJD_guides)
        except NameError:
            MJD_guideslist = [MJD_guides]
        hdulist.close()

    print "UTSTART, MJD lists=", len(utstart_guidelist), len(MJD_guideslist)
    full_guides_arr = array(full_guides_all)
    i = 0
    utstart_guide_arr = array(utstart_guidelist)
    MJD_guides_arr = array(MJD_guideslist)
    for _acq in acqimages:

        timediff = utstart_guide_arr - acq_utstop_arr[i]
        mjddiff = MJD_guides_arr - acq_mjd_arr[i]

        timediff_cull = timediff * \
            (timediff < 150.0) * (timediff > 0.0) * (mjddiff > 0.0)
        timediff_gd = timediff[where(timediff_cull != 0)]
        guides_gd = full_guides_arr[where(timediff_cull != 0)]

        try:
            if not guides_gd.size:
                guide_list.append('Null')
            else:
                guide_list.append(guides_gd[0])

        except:
            if not guides_gd.size:
                guide_list = ['Null']
            else:
                guide_list = [guides_gd[0]]

        i = i + 1

    return guide_list

# def
# mk_obs_website(acqimage,grpid,propid,UTstartnocolon,night,guideimage,data_dir,specslist,guidecountplots):


def mk_obs_website(acqimage, grpid, propid, UTstartnocolon, night, guideimage, data_dir, specslist, sciencespecs):

    new_dir = '/var/www/html/night_summary/' + night + '/' + \
        grpid + '_' + str(UTstartnocolon) + '_' + propid
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)

#   acq_root = '/icc2/tmp/'
#    spec_root='/icc/tmp/'
    acq_root = os.path.dirname(acqimage)
    print sciencespecs
    spec_root = '/icc/tmp/'
    if len(sciencespecs) > 0:
        spec_root = os.path.dirname(sciencespecs[0])
    print acq_root, spec_root

    acqimage = acqimage.replace(acq_root, data_dir + '/')
    acq_jpg = acqimage.replace('.fits', '.jpg')
    acq_reg = acq_jpg.replace('.jpg', '.reg')
    acq_cat = acq_jpg.replace('.jpg', '.cat')
    acq_log = acq_jpg.replace('.jpg', '.log')
    guideimage = guideimage.replace(acq_root, data_dir + '/')
    guide_jpg = guideimage.replace('.fits', '.jpg')
    print "  acqimage,jpg=", acqimage, acq_jpg
    print "guideimage,jpg=", guideimage, guide_jpg

    yoindex = acqimage.rfind('/')
    acq_in_curr = acqimage[yoindex + 1:]
    acqjpg_in_curr = acq_in_curr.replace('.fits', '.jpg')
    acqreg_in_curr = acq_in_curr.replace('.fits', '.reg')
    acqlog_in_curr = acq_in_curr.replace('.fits', '.log')
    acqcat_in_curr = acq_in_curr.replace('.fits', '.cat')

    yoyoindex = guideimage.rfind('/')
    guide_in_curr = guideimage[yoyoindex + 1:]
    guidejpg_in_curr = guide_in_curr.replace('.fits', '.jpg')

    try:
        os.remove(new_dir + '/' + acqjpg_in_curr)
    except:
        pass
    try:
        os.remove(new_dir + '/' + acq_in_curr)
    except:
        pass
    try:
        os.remove(new_dir + '/' + acqreg_in_curr)
    except:
        pass
    try:
        os.remove(new_dir + '/' + acqcat_in_curr)
    except:
        pass
    try:
        os.remove(new_dir + '/' + acqlog_in_curr)
    except:
        pass
    try:
        os.remove(new_dir + '/' + guide_in_curr)
    except:
        pass
    try:
        os.remove(new_dir + '/' + guidejpg_in_curr)
    except:
        pass

    link_cmd = 'ln -s ' + acqimage + ' ' + new_dir + '/'
    print "link_cmd1=", link_cmd
    os.system(link_cmd)
    link_cmd = 'ln -s ' + acq_jpg + ' ' + new_dir + '/'
    print "link_cmd2=", link_cmd
    os.system(link_cmd)

    if guideimage != 'Null':
        link_cmd = 'ln -s ' + guideimage + ' ' + new_dir + '/'
    else:
        link_cmd = 'ln -sf /var/www/html/images/no_guide.png ' + new_dir + '/'
    os.system(link_cmd)
    print "link_cmd3=", link_cmd

    if guide_jpg != 'Null':
        link_cmd = 'ln -s ' + guide_jpg + ' ' + new_dir + '/'
        os.system(link_cmd)
    else:
        link_cmd = "No guide image"
    print "link_cmd4=", link_cmd

    link_cmd = 'ln -s ' + acq_reg + ' ' + new_dir + '/'
    print "link_cmd5=", link_cmd
    os.system(link_cmd)
    link_cmd = 'ln -s ' + acq_cat + ' ' + new_dir + '/'
    print "link_cmd6=", link_cmd
    os.system(link_cmd)
    link_cmd = 'ln -s ' + acq_log + ' ' + new_dir + '/'
    print "link_cmd7=", link_cmd
    os.system(link_cmd)

    tarname = new_dir + '/' + grpid + '_' + propid + '.tar'
#    os.system('rm -f *guide*.png')
#    os.system('cp '+new_dir+'/*guide*.png .')
    os.system('cp ' + acqimage + ' ./')
    tarfiles = acq_in_curr
    print tarfiles
    if os.path.exists(guideimage):
        os.system('cp ' + guideimage + ' ./')
        tarfiles = tarfiles + ' ' + guide_in_curr
        print tarfiles
    home = os.path.expanduser('~')

    for _sciencespec in sciencespecs:
        science_name = os.path.basename(_sciencespec)
        guidecounts_pl = science_name.replace('.fits', '_guidecounts.png')
        print "guidecounts_pl =", guidecounts_pl
        guidecounts_pl = os.path.join(home, guidecounts_pl)
        print "guidecounts_pl =", guidecounts_pl
        guidecounts_curr = os.path.basename(guidecounts_pl)
        guidefwhmt_pl = guidecounts_pl.replace(
            'guidecounts.png', 'guidefwhmt.png')
        guidestate_pl = guidecounts_pl.replace(
            'guidecounts.png', 'guidestate.png')
        guidext_pl = guidecounts_pl.replace('guidecounts.png', 'guidext.png')
        guideyt_pl = guidecounts_pl.replace('guidecounts.png', 'guideyt.png')
        guidexy_pl = guidecounts_pl.replace('guidecounts.png', 'guidexy.png')

        guidefwhmt_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidefwhmt.png')
        guidestate_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidestate.png')
        guidext_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidext.png')
        guideyt_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guideyt.png')
        guidexy_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidexy.png')

#        science_name = os.path.basename(_sciencespec)
#        guidecounts_pl = _sciencespec.replace('.fits','_guidecounts.png')
#        guidecounts_pl = guidecounts_pl.replace(spec_root + 'f_','/home/eng/f_')
#        guidecounts_pl = guidecounts_pl.replace(spec_root + 'g_','/home/eng/g_')
#        guidecounts_curr = guidecounts_pl.replace('/home/eng/','')
#        guidefwhmt_pl = guidecounts_pl.replace('guidecounts.png','guidefwhmt.png')
#        guidext_pl = guidecounts_pl.replace('guidecounts.png','guidext.png')
#        guideyt_pl = guidecounts_pl.replace('guidecounts.png','guideyt.png')
#        guidexy_pl = guidecounts_pl.replace('guidecounts.png','guidexy.png')

#        guidefwhmt_curr = guidecounts_curr.replace('guidecounts.png','guidefwhmt.png')
#        guidext_curr = guidecounts_curr.replace('guidecounts.png','guidext.png')
#        guideyt_curr = guidecounts_curr.replace('guidecounts.png','guideyt.png')
#        guidexy_curr = guidecounts_curr.replace('guidecounts.png','guidexy.png')

        if os.path.exists(guidecounts_pl):
            guidecounts_tar = os.path.basename(guidecounts_pl)
            tarfiles = tarfiles + ' ' + guidecounts_tar
        if os.path.exists(guidefwhmt_curr):
            guidefwhm_tar = os.path.basename(guidefwhmt_curr)
            tarfiles = tarfiles + ' ' + guidefwhm_tar
        if os.path.exists(guidext_curr):
            guidext_tar = os.path.basename(guidext_curr)
            tarfiles = tarfiles + ' ' + guidext_tar
        if os.path.exists(guideyt_curr):
            guideyt_tar = os.path.basename(guideyt_curr)
            tarfiles = tarfiles + ' ' + guideyt_tar
        if os.path.exists(guidexy_curr):
            guidexy_tar = os.path.basename(guidexy_curr)
            tarfiles = tarfiles + ' ' + guidexy_tar
        if os.path.exists(guidestate_curr):
            guidestate_tar = os.path.basename(guidestate_curr)
            tarfiles = tarfiles + ' ' + guidestate_tar

    os.system('rm -f ' + tarname)
    os.system('tar -cvf ' + tarname + ' ' + tarfiles)

    os.system('rm -f ' + acq_in_curr)
    if os.path.exists(guide_in_curr):
        os.system('rm -f ' + guide_in_curr)

# file translations and path additions

    guide_jpg = guideimage.replace('.fits', '.jpg')
    guide_wpath = data_dir + '/' + guideimage
    guidejpg_wpath = data_dir + '/' + guide_jpg

    html_file = new_dir + '/' + grpid + '_' + propid + '.html'
    if os.path.exists(html_file):
        os.remove(html_file)
    outfile = open(html_file, 'w')
    name = grpid + '_' + propid

    no_image_html = '<IMG src="no_guide.png" height="400">\n'
    outfile.write('<html>\n')
    outfile.write('<head>\n')
    outfile.write('<title>' + name + '</title>\n')
    outfile.write('</head>\n')
    outfile.write('<center>\n')
    outfile.write('<h1><font COLOR=red> Quick summary <br> ' +
                  name + '</h1>\n')
    outfile.write('</center>\n')
    outfile.write('<br><br><br>\n')
    outfile.write('<table>\n')
    outfile.write('<tr ALIGN="center">\n')
    outfile.write('<td COLSPAN=1><b> Acquisition image </b></td>\n')
    outfile.write('<td COLSPAN=1><b> First Guide Image</b></td>\n')
    outfile.write('</tr>\n')
    outfile.write('<tr>\n')
    outfile.write('<th ROWSPAN=2>\n')
    outfile.write('<IMG src="' + acqjpg_in_curr + '" height="400">\n')
    outfile.write('</th>\n')
    outfile.write('<th ROWSPAN=2>\n')
    if guidejpg_in_curr != 'Null':
        outfile.write('<IMG src="' + guidejpg_in_curr + '" height="400">\n')
    else:
        outfile.write(no_image_html)
    outfile.write('</th>\n')
    outfile.write('</tr>\n')
    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    outfile.write('<td COLSPAN=1> <a href=./' +
                  acq_in_curr + '> ' + acq_in_curr + '</td>\n')
    if guide_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./' +
                      guide_in_curr + '> ' + guide_in_curr + '</td>')
    else:
        outfile.write('<td COLSPAN=1></td>')
    outfile.write('</tr>')
    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    if acqreg_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./' +
                      acqlog_in_curr + '> ' + acqlog_in_curr + '</td>\n')
    outfile.write('</tr>')
    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    if acqreg_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./' +
                      acqreg_in_curr + '> ' + acqreg_in_curr + '</td>\n')
    outfile.write('</tr>')

    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    if acqcat_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./' +
                      acqcat_in_curr + '> ' + acqcat_in_curr + '</td>\n')
    outfile.write('</tr>')

    outfile.write('</table>')

# search for guide plots corresponding to each science spec
    home = os.path.expanduser('~')
    for _sciencespec in sciencespecs:
        science_name = os.path.basename(_sciencespec)
        guidecounts_pl = science_name.replace('.fits', '_guidecounts.png')
        print "guidecounts_pl =", guidecounts_pl
        guidecounts_pl = os.path.join(home, guidecounts_pl)
        print "guidecounts_pl =", guidecounts_pl
        guidecounts_curr = os.path.basename(guidecounts_pl)
        guidefwhmt_pl = guidecounts_pl.replace(
            'guidecounts.png', 'guidefwhmt.png')
        guidestate_pl = guidecounts_pl.replace(
            'guidecounts.png', 'guidestate.png')
        guidext_pl = guidecounts_pl.replace('guidecounts.png', 'guidext.png')
        guideyt_pl = guidecounts_pl.replace('guidecounts.png', 'guideyt.png')
        guidexy_pl = guidecounts_pl.replace('guidecounts.png', 'guidexy.png')

        guidefwhmt_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidefwhmt.png')
        guidestate_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidestate.png')
        guidext_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidext.png')
        guideyt_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guideyt.png')
        guidexy_curr = guidecounts_curr.replace(
            'guidecounts.png', 'guidexy.png')

        os.system('mv ' + guidecounts_pl + ' ' + new_dir)
        os.system('mv ' + guidefwhmt_pl + ' ' + new_dir)
        os.system('mv ' + guidestate_pl + ' ' + new_dir)
        os.system('mv ' + guidext_pl + ' ' + new_dir)
        os.system('mv ' + guideyt_pl + ' ' + new_dir)
        os.system('mv ' + guidexy_pl + ' ' + new_dir)

        outfile.write('<table> \n <tr ALIGN="center"> \n <td COLSPAN=1><b> Guide Plots For ' +
                      science_name + '  </b></td> \n </tr> \n')
        outfile.write('<tr>')
        outfile.write('<th ROWSPAN=2>')
# if guide
        if guidejpg_in_curr != 'Null':
            outfile.write('<IMG src="' + guidecounts_curr + '" height="300">')
            outfile.write('</th>')

            outfile.write('<th ROWSPAN=2>')
            outfile.write('<IMG src="' + guidefwhmt_curr + '" height="300">')
            outfile.write('</th>')

            outfile.write('<th ROWSPAN=2>')
            outfile.write('<IMG src="' + guidestate_curr + '" height="300">')
            outfile.write('</th>')
            outfile.write('</tr>')

            outfile.write('<tr>')
            outfile.write('</tr>')
            outfile.write('</table>')

            outfile.write('<table>')
            outfile.write('<tr>')
            outfile.write('<th ROWSPAN=2>')
            outfile.write('<IMG src="' + guidexy_curr + '" height="300">')
            outfile.write('</th>')
            outfile.write('<th ROWSPAN=2>')

            outfile.write('<IMG src="' + guidext_curr + '" height="300">')
            outfile.write('</th>')
            outfile.write('<th ROWSPAN=2>')

            outfile.write('<IMG src="' + guideyt_curr + '" height="300">')
            outfile.write('</th>')
            outfile.write('</tr>')
        else:
            outfile.write('<td COLSPAN=2>' + no_image_html + '</td>')
            outfile.write('</th>')
            outfile.write('</tr>')

        outfile.write('</table>')

    outfile.write(
        '<table> \n <tr ALIGN="center"> \n <td COLSPAN=1><b> Spectra </b></td> \n </tr> \n')

    for _specs in specslist:
        if debug > 2:
            print "specs before=", _specs
        _specs = _specs.replace(spec_root, data_dir + '/')
        if '/var/www/html' not in _specs:
            _specsjpg = _specs.replace(
                '02.fits', '01.jpg').replace('flash/', 'flash/jpg/')
        else:
            _specsjpg = _specs.replace('.fits', '.jpg')

        if debug > 2:
            print "specs  after=", _specs
        if debug > 2:
            print "specsjpg=", _specsjpg
        specindex = _specs.rfind('/')
        spec_in_curr = _specs[specindex + 1:]
        specjpg_in_curr = spec_in_curr.replace('.fits', '.jpg')
        if debug > 2:
            print "specjpg_in_curr=", specjpg_in_curr
        try:
            os.remove(new_dir + '/' + spec_in_curr)
        except:
            pass
        try:
            os.remove(new_dir + '/' + specjpg_in_curr)
        except:
            pass

        link_spec = 'ln -s ' + _specs + ' ' + new_dir + '/'
        retval = os.system(link_spec)
        print "Spec link cmd,status=", link_spec, retval
        link_spec = 'ln -s ' + _specsjpg + ' ' + new_dir + '/' + specjpg_in_curr
        retval = os.system(link_spec)
        print "Spec link cmd,status=", link_spec, retval

        outfile.write('<tr>')
        outfile.write('<th ROWSPAN=2>')
        outfile.write('<IMG src="' + specjpg_in_curr + '" height="200">')
        outfile.write('</th>')
        outfile.write('</tr>')
        outfile.write('<tr></tr>')
        outfile.write('<tr></tr>')
        outfile.write('<tr ALIGN="center">')
        outfile.write('<td COLSPAN=1> <a href=./' +
                      spec_in_curr + '> ' + spec_in_curr + '</td>')
        outfile.write('</tr>')

    try:
        os.remove(new_dir + '/Null')
    except:
        pass


def which(program):
    import os

    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


def make_jpg_image(jpgfile, jpgfile_destname, debug=False):
    '''Create a JPG in <jpgfile_destpath> in the same way as the reduction agent should have done:
    /usr/local/astrometry/bin/image2pnm.py -i $1 -o $1.tmp1.pnm
    pnmtojpeg --greyscale --quality=85 $1.tmp1.pnm > $1.tmp1.jpg
    jpegtran -copy none -optimize -flip vertical -outfile $1.tmp2.jpg $1.tmp1.jpg
    jpegtran -copy none -progressive -outfile $2 $1.tmp2.jpg
    rm $1.tmp*
    '''

    if debug:
        print jpgfile, jpgfile_destname

    needed_programs = [
        '/usr/local/astrometry/bin/image2pnm.py', 'pnmtojpeg', 'jpegtran']
    status = 0
    for prog in needed_programs:
        if debug:
            print "Checking for", prog
        if which(prog) == None:
            if debug:
                print "Not found"
            status = -1

    if status != -1:
        fitsfilepath = jpgfile.replace('.jpg', '.fits')
        fitsfile = os.path.basename(fitsfilepath)
        jpgfile_destpath = os.path.dirname(jpgfile_destname) + os.path.sep

# Convert FITS to PNM
        infile = fitsfilepath
        outfile = os.path.join(jpgfile_destpath, fitsfile + '.tmp1.pnm')
        cmd = '/usr/local/astrometry/bin/image2pnm.py -i ' + fitsfilepath + ' -o ' + outfile
        if debug:
            print cmd
        else:
            status = os.system(cmd)

        if status == 0:
            # Convert PNM to greyscale JPG
            infile = outfile
            outfile = os.path.join(jpgfile_destpath, fitsfile + '.tmp1.jpg')
            cmd = 'pnmtojpeg --greyscale --quality=85 ' + infile + ' > ' + outfile
            if debug:
                print cmd
            else:
                status = os.system(cmd)

# Optimize and flip JPG
            infile = outfile
            outfile = outfile.replace('.tmp1.jpg', '.tmp2.jpg')
            cmd = 'jpegtran -copy none -optimize -flip vertical -outfile ' + outfile + ' ' + infile
            if debug:
                print cmd
            else:
                status = os.system(cmd)

# Convert JPG to progressive format
            infile = outfile
            outfile = jpgfile_destname
            cmd = 'jpegtran -copy none -progressive -outfile ' + outfile + ' ' + infile
            if debug:
                print cmd
            else:
                status = os.system(cmd)
# Cleanup
            cmd = 'rm ' + jpgfile_destpath + fitsfile + '.tmp*'
            if debug:
                print cmd
            else:
                status = os.system(cmd)
    return status


def readlogfile(logname):
    import csv
    """ read in the guide log file """
    f = open(logname)
    headerlines = 1
    # There is one header line in these guider log files.
    nlines = sum(1 for line in f) - headerlines
    f.close()
    logf = open(logname, "rb")
    for line in range(headerlines):
        logf.next()  # skip header
    logreadin = csv.reader(logf, delimiter=" ", skipinitialspace=True)
    """ initialize arrays """
    time = nlines * [None]

    nrow = 0
    for row in logreadin:
        time[nrow] = float(row[0])
        xguide[nrow] = float(row[5])
        yguide[nrow] = float(row[6])
        nrow = nrow + 1

if __name__ == '__main__':
    from optparse import OptionParser
    from datetime import datetime
    parser = OptionParser()
    parser.add_option("-n", "--nightlist", dest="nightlist",
                      help="Nights to aggregate data", default='now')
    parser.add_option('-t', '--tmp', dest="tmp_dir",
                      help="Temporary directory", default="./")
    parser.add_option('-s', '--site', dest="site",
                      help="site name, ogg or coj", default="floyds.coj.lco.gtn")
    parser.add_option('-d', dest="debug", help="debug output enabled",
                      action="store_true", default=False)

    # Add others
    (options, args) = parser.parse_args()
    site = options.site
    nightlist = options.nightlist
    debug = options.debug

    if nightlist == 'now':
        yo = str(datetime.date(datetime.now()))
        yo = yo.replace('-', '')
        nightlist = yo
        if site == 'floyds.ogg.lco.gtn':
            nightlist_four = nightlist[4:]
            nightlist_year = nightlist[:4]
            # print nightlist_four
            if nightlist_four == '0101':
                nightlist = nightlist_year + '1231'
            elif nightlist_four == '0201':
                nightlist = nightlist_year + '0131'
            elif nightlist_four == '0301':
                nightlist = nightlist_year + '0228'
            elif nightlist_four == '0401':
                nightlist = nightlist_year + '0331'
            elif nightlist_four == '0501':
                nightlist = nightlist_year + '0430'
            elif nightlist_four == '0601':
                nightlist = nightlist_year + '0531'
            elif nightlist_four == '0701':
                nightlist = nightlist_year + '0630'
            elif nightlist_four == '0801':
                nightlist = nightlist_year + '0731'
            elif nightlist_four == '0901':
                nightlist = nightlist_year + '0831'
            elif nightlist_four == '1001':
                nightlist = nightlist_year + '0930'
            elif nightlist_four == '1101':
                nightlist = nightlist_year + '1031'
            elif nightlist_four == '1201':
                nightlist = nightlist_year + '1130'
            else:
                nightlist = str(int(nightlist) - 1)
    if nightlist == 'yesterday':
        yo = str(datetime.date(datetime.now()))

    tmp_dir = options.tmp_dir

    print "Running for nightlist", nightlist
    agg_floyds(nightlist, site=site, tmp_dir=tmp_dir, debug=debug)

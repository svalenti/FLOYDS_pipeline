#!/opt/epd/bin/python

import os
import matplotlib
matplotlib.use('Agg')

import pyfits
import numpy, scipy, pylab, matplotlib, sys, math, signal,subprocess,shlex,time,logging,textwrap,glob
import string
import plot_guideinfo
from pyfits import open as popen
import numpy
from numpy import array,abs,where


def agg_floyds(nightlist,site='floyds.coj.lco.gtn',tmp_dir="./"):

    """Set up basic logging config"""
    logfile = 'agglog.txt'
    path,logfile = os.path.split(logfile)
    logging.basicConfig(filename=tmp_dir+logfile,filemode="w",level=logging.DEBUG,format='%(asctime)s %(levelname)s:%(message)s',datefmt='%m/%d/%Y %I:%M:%S %p')
    logging.info('Log for agg_floyds -- mostly for debug purposes')

    #os.system('export LD_LIBRARY_PATH=":/icc/bin/lib"')

#split the nights up, comma delimited
    night_list = string.split(nightlist,sep=',')


#aggregate data for each night
    i=0
    for _night in night_list:

#grab last four characters of _night
        _night_four = _night[4:]
        _night_year = _night[:4]
        dirnow = '/var/www/html/night_summary/'+_night
        dir_for_guideimgs = '/var/www/html/images/'+_night
        if not os.path.exists(dirnow):
            os.makedirs(dirnow)
        acq_images = glob.glob('/icc2/tmp/[s,r]*'+_night+'*.fits')
        spec_images = glob.glob('/icc/tmp/[g,f]*'+_night+'*.fits')
        if (acq_images == []) & (spec_images == []):
            print 'No images take on this night'
            sys.exit(0)
#Get the guider log
        if site == 'floyds.ogg.lco.gtn':
            if _night_four == '0131':
                guidelog = '/icc2/tmp/'+_night_year+'0201.log'
            elif _night_four == '0228':
                guidelog = '/icc2/tmp/'+_night_year+'0301.log'
            elif _night_four == '0331':
                guidelog = '/icc2/tmp/'+_night_year+'0401.log'
            elif _night_four == '0430':
                guidelog = '/icc2/tmp/'+_night_year+'0501.log'
            elif _night_four == '0531':
                guidelog = '/icc2/tmp/'+_night_year+'0601.log'
            elif _night_four == '0630':
                guidelog = '/icc2/tmp/'+_night_year+'0701.log'
            elif _night_four == '0731':
                guidelog = '/icc2/tmp/'+_night_year+'0801.log'
            elif _night_four == '0831':
                guidelog = '/icc2/tmp/'+_night_year+'0901.log'
            elif _night_four == '0930':
                guidelog = '/icc2/tmp/'+_night_year+'1001.log'
            elif _night_four == '1031':
                guidelog = '/icc2/tmp/'+_night_year+'1101.log'
            elif _night_four == '1130':
                guidelog = '/icc2/tmp/'+_night_year+'1201.log'
            elif _night_four == '1231':
                guidelog = '/icc2/tmp/'+_night_year+'0101.log'
            else:
                guidelog = '/icc2/tmp/'+str(int(_night)+1)+'.log'


        else:
            guidelog = '/icc2/tmp/'+str(int(_night))+'.log'
        if os.path.isfile(guidelog):
            #print guidelog
#            pass
            guidetag = 1
        else:
            print 'No guide file found!!!'
            guidetag = 0

#which acquisition images have region files?  acquisition Log files?
#figure out which ones do and: 1) link them

        os.system('rm -f '+dir_for_guideimgs+'/*.reg')
        os.system('rm -f '+dir_for_guideimgs+'/*.log')
        os.system('rm -f '+dir_for_guideimgs+'/*.cat')
        #print acq_images
        for _acqimage in acq_images:
            #print os.path.getsize(_acqimage)
            if os.path.getsize(_acqimage) < 100000:
                badname = _acqimage
                acq_images.remove(badname)


            #try:
            #    testit = pyfits.info(_acqimage,ignore_missing_end=False)
            #except:
            #    badname = _acqimage
            #    print 'did I make it here'
            #    acq_images.remove[badname]
        #print acq_images

        for _acqimage in acq_images:
#which acquisition images have region files?  acquisition Log files?
#figure out which ones do and: 1) link them

            yo = _acqimage.replace('/icc2/tmp/','/lco/floyds/tmp/')
            regfile = yo.replace('.fits','.reg')
            logfile = yo.replace('.fits','.log')
            catfile = yo.replace('.fits','.cat')

            if os.path.exists(regfile):
                os.system('ln -s '+regfile+' '+ dir_for_guideimgs+'/')
            if os.path.exists(logfile):
                os.system('ln -s '+logfile+' '+dir_for_guideimgs+'/')
            if os.path.exists(catfile):
                os.system('ln -s '+catfile+' '+dir_for_guideimgs+'/')

#read in each acquisition header
            #print _acqimage
            hdulist = pyfits.open(_acqimage)
            #print _acqimage
            prihdr = hdulist[0].header

            propid_acq = prihdr['PROPID']
            try:
                propid_acqlist.append(propid_acq)
            except NameError:
                propid_acqlist = [propid_acq]

            MJD_acq = prihdr['MJD']
            try:
                MJD_acqlist.append(MJD_acq)
            except NameError:
                MJD_acqlist = [MJD_acq]

            grpuid_acq = int(prihdr['GRPUID'])
            try:
                grpuid_acqlist.append(grpuid_acq)
            except NameError:
                grpuid_acqlist = [grpuid_acq]

            grpnumob_acq = int(prihdr['GRPNUMOB'])
            try:
                grpnumob_acqlist.append(grpnumob_acq)
            except NameError:
                grpnumob_acqlist = [grpnumob_acq]

            groupid_acq = prihdr['GROUPID']
            try:
                grpid_acqlist.append(groupid_acq)
            except NameError:
                grpid_acqlist = [groupid_acq]

            acq_utstop = prihdr['UTSTOP']
            try:
                acq_utstoplist.append(acq_utstop)
                acq_utstop_nocolon.append(float(acq_utstop.replace(':','')))
            except NameError:
                acq_utstoplist = [acq_utstop]
                acq_utstop_nocolon = [float(acq_utstop.replace(':',''))]

            acq_start = prihdr['UTSTART']
            try:
                acq_utstartlist.append(acq_start)
                acq_utstart_nocolon.append(float(acq_start.replace(':','')))
            except NameError:
                acq_utstartlist = [acq_start]
                acq_utstart_nocolon = [float(acq_start.replace(':',''))]



        if guidetag != 0:
            first_guideimage = find_first_guide(_night,acq_images,acq_utstart_nocolon,acq_utstop_nocolon,MJD_acqlist)
        else:
            first_guideimage = 'Null'
        print first_guideimage
        print dir_for_guideimgs
        os.system('rm -f '+dir_for_guideimgs+'/*f.fits')
        os.system('rm -f '+dir_for_guideimgs+'/*f.jpg')
        for _first in first_guideimage:
            if _first != 'Null':
                print 'ln -s '+_first+' '+dir_for_guideimgs+'/'
                os.system('ln -s '+_first+' '+dir_for_guideimgs+'/')
                _first_jpg = _first.replace('.fits','.jpg')
                os.system('/usr/local/bin/gpp '+_first+' '+_first_jpg+' 2 /usr/local/bin/gpp_cfg.txt')
                os.system('mv '+_first_jpg+' '+dir_for_guideimgs)

        for _specimage in spec_images:
#read in and collect header info
            hdulist = pyfits.open(_specimage)
            prihdr = hdulist[0].header
            propid_spec = prihdr['PROPID']
            try:
                propid_speclist.append(propid_spec)
            except NameError:
                propid_speclist = [propid_spec]

            grpuid_spec = int(prihdr['GRPUID'])
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
                utstart_speclist_nocolon.append(float(utstart_spec.replace(':','')))
            except NameError:
                utstart_speclist = [utstart_spec]
                utstart_speclist_nocolon = [float(utstart_spec.replace(':',''))]

            exptime_spec = prihdr['EXPTIME']
            try:
                exptime_speclist.append(float(exptime_spec))

            except NameError:
                exptime_speclist = [float(exptime_spec)]


            utstop_spec = prihdr['UTSTOP']
            try:
                utstop_speclist.append(utstop_spec)
                utstop_speclist_nocolon.append(float(utstop_spec.replace(':','')))
            except NameError:
                utstop_speclist = [utstop_spec]
                utstop_speclist_nocolon = [float(utstop_spec.replace(':',''))]


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


        #print obstype_speclist_arr
        #print utstart_speclist_arr
        #print utstop_speclist_arr

    #sys.exit(0)
        i=0
        for _acq in acq_images:
            # find those spectra with the same grpuid
            grpuid_curr = grpuid_acqlist[i]
            specs_forthis_acq = spec_images_arr[where(grpuid_speclist_arr == grpuid_curr)]
            specs_forthis_acqlist = specs_forthis_acq.tolist()

            science_cull = (grpuid_speclist_arr == grpuid_curr) * (obstype_speclist_arr == 'SPECTRUM')

            science_index = where(science_cull == True)

            sciencespecs_forthis_acq = spec_images_arr[science_index]

            sciencespecs_forthis_acqlist = sciencespecs_forthis_acq.tolist()

            utstart_science = utstart_speclist_arr[science_index]
            utstop_science = utstop_speclist_arr[science_index]
            exptime_science = exptime_speclist_arr[science_index]
            j=0
            for _utstart in utstart_science:

                sciencespecnow = sciencespecs_forthis_acq[j]
                junk1 = sciencespecnow.replace('.fits','')
                guideplot_rootname = junk1.replace('/icc/tmp/','')
                #print guidetag
                if guidetag > 0:
                    #print guideplot_rootname
                    #print _utstart
                    #print utstop_science[j]
                    guidecheck = plot_guideinfo.mk_guideinfo_plots(guideplot_rootname,_utstart,utstop_science[j],guidelog)
                else:
                    guidecheck = 0
                    pass
                j=j+1

            mk_obs_website(acq_images[i],grpid_acqlist[i],propid_acqlist[i],acq_utstart_nocolon[i],_night,first_guideimage[i],dir_for_guideimgs,specs_forthis_acqlist,sciencespecs_forthis_acq)
            i=i+1



#read in each spectroscopy header
#    for _specimage in spec_images:
#        hdulist = pyfits.open(_specimage)
#        prihdr = hdulist[0].header

# for each acquisition image, look for acquisition log 

def find_first_guide(night,acqimages,acq_utstart,acq_utstop,acq_mjd):

    from numpy import array,abs,where

    imgdir = '/icc2/tmp/'

#I need a kludge in here to include nights one day before and after.  there is some weird timing issue at FTN.

    night_year = night[:4]
    night_four = night[4:]

    if night_four == '0131':
        nightnext = night_year+'0201'
    elif night_four == '0228':
        nightnext = night_year+'0301'
    elif night_four == '0331':
        nightnext = night_year+'0401'
    elif night_four == '0430':
        nightnext = night_year+'0501'
    elif night_four == '0531':
        nightnext = night_year+'0601'
    elif night_four == '0630':
        nightnext = night_year+'0701'
    elif night_four == '0731':
        nightnext = night_year+'0801'
    elif night_four == '0831':
        nightnext = night_year+'0901'
    elif night_four == '0930':
        nightnext = night_year+'1001'
    elif night_four == '1031':
        nightnext = night_year+'1101'
    elif night_four == '1130':
        nightnext = night_year+'1201'
    elif night_four == '1231':
        nightnext = night_year+'0101'
    else:
        nightnext = str(int(night)+1)
    

#find full frame guide images from the night in question


    #print imgdir+nightnext+'*f.fits'
    full_guides_all = glob.glob(imgdir+night+'*f.fits')+glob.glob(imgdir+nightnext+'*f.fits')
    #print full_guides_all
    acq_utstop_arr = array(acq_utstop)
    acq_mjd_arr = array(acq_mjd)

#get the UT start time for all of the full frame guide images
    for _guides in full_guides_all:
        hdulist = pyfits.open(_guides)
        prihdr = hdulist[0].header
        utstart_guides = prihdr['UTSTART']
        MJD_guides = prihdr['MJD']
        try:
            utstart_guidelist.append(float(utstart_guides.replace(':','')))
        except NameError:
            utstart_guidelist = [float(utstart_guides.replace(':',''))]
        try:
            MJD_guideslist.append(MJD_guides)
        except NameError:
            MJD_guideslist = [MJD_guides]


    full_guides_arr = array(full_guides_all)
    i=0
    utstart_guide_arr = array(utstart_guidelist)
    MJD_guides_arr = array(MJD_guideslist)
    for _acq in acqimages:

        timediff = utstart_guide_arr - acq_utstop_arr[i]
        mjddiff = MJD_guides_arr - acq_mjd_arr[i]

        timediff_cull = timediff * (timediff < 150.0) * (timediff > 0.0) * (mjddiff > 0.0)
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



        i=i+1

    return guide_list

#def mk_obs_website(acqimage,grpid,propid,UTstartnocolon,night,guideimage,data_dir,specslist,guidecountplots):
def mk_obs_website(acqimage,grpid,propid,UTstartnocolon,night,guideimage,data_dir,specslist,sciencespecs):

    new_dir = '/var/www/html/night_summary/'+night+'/'+grpid+'_'+str(UTstartnocolon)+'_'+propid
    if not os.path.exists(new_dir):
        os.makedirs(new_dir)
    



    acqimage = acqimage.replace('/icc2/tmp/',data_dir+'/')
    acq_jpg = acqimage.replace('.fits','.jpg')
    acq_reg = acq_jpg.replace('.jpg','.reg')
    acq_cat = acq_jpg.replace('.jpg','.cat')
    acq_log = acq_jpg.replace('.jpg','.log')
    guideimage = guideimage.replace('/icc2/tmp/',data_dir+'/')
    guide_jpg = guideimage.replace('.fits','.jpg')
    print "guideimage,guide_jpg",guideimage,guide_jpg


    yoindex = acqimage.rfind('/')
    acq_in_curr = acqimage[yoindex+1:]
    acqjpg_in_curr = acq_in_curr.replace('.fits','.jpg')
    acqreg_in_curr = acq_in_curr.replace('.fits','.reg')
    acqlog_in_curr = acq_in_curr.replace('.fits','.log')
    acqcat_in_curr = acq_in_curr.replace('.fits','.cat')


    yoyoindex = guideimage.rfind('/')
    guide_in_curr = guideimage[yoyoindex+1:]
    guidejpg_in_curr = guide_in_curr.replace('.fits','.jpg')

    print "guide_in_curr,guidejpg_in_curr",guide_in_curr,guidejpg_in_curr
    try:
        os.remove(new_dir+'/'+acqjpg_in_curr)
    except:
        pass
    try:
        os.remove(new_dir+'/'+acq_in_curr)
    except:
        pass
    try:
        os.remove(new_dir+'/'+acqreg_in_curr)
    except:
        pass
    try:
        os.remove(new_dir+'/'+acqcat_in_curr)
    except:
        pass
    try:
        os.remove(new_dir+'/'+acqlog_in_curr)
    except:
        pass
    try:
        os.remove(new_dir+'/'+guide_in_curr)
    except:
        pass
    try:
        os.remove(new_dir+'/'+guidejpg_in_curr)
    except:
        pass

    os.system('ln -s '+acqimage+' '+new_dir+'/')
    os.system('ln -s '+acq_jpg+' '+new_dir+'/')
    os.system('ln -s '+guideimage+' '+new_dir+'/')
    os.system('ln -s '+guide_jpg+' '+new_dir+'/')
    os.system('ln -s '+acq_reg+' '+new_dir+'/')
    os.system('ln -s '+acq_cat+' '+new_dir+'/')
    os.system('ln -s '+acq_log+' '+new_dir+'/')

    tarname = new_dir+'/'+grpid+'_'+propid+'.tar'
#    os.system('rm -f *guide*.png')
#    os.system('cp '+new_dir+'/*guide*.png .')
    os.system('cp '+acqimage+' ./')
    tarfiles = acq_in_curr
    print tarfiles
    if os.path.exists(guideimage):
        os.system('cp '+guideimage+' ./')
        tarfiles = tarfiles+' '+guide_in_curr
        print tarfiles

    for _sciencespec in sciencespecs:
        science_name = _sciencespec.replace('/icc/tmp/','')
        guidecounts_pl = _sciencespec.replace('.fits','_guidecounts.png')
        guidecounts_pl = guidecounts_pl.replace('/icc/tmp/f_','/home/eng/f_')
        guidecounts_pl = guidecounts_pl.replace('/icc/tmp/g_','/home/eng/g_')
        guidecounts_curr = guidecounts_pl.replace('/home/eng/','')
        guidefwhmt_pl = guidecounts_pl.replace('guidecounts.png','guidefwhmt.png')
        guidext_pl = guidecounts_pl.replace('guidecounts.png','guidext.png')
        guideyt_pl = guidecounts_pl.replace('guidecounts.png','guideyt.png')
        guidexy_pl = guidecounts_pl.replace('guidecounts.png','guidexy.png')

        guidefwhmt_curr = guidecounts_curr.replace('guidecounts.png','guidefwhmt.png')
        guidext_curr = guidecounts_curr.replace('guidecounts.png','guidext.png')
        guideyt_curr = guidecounts_curr.replace('guidecounts.png','guideyt.png')
        guidexy_curr = guidecounts_curr.replace('guidecounts.png','guidexy.png')

        if os.path.exists(guidecounts_pl):
            guidecounts_tar = guidecounts_pl.replace('/home/eng/','')
            tarfiles = tarfiles+' '+guidecounts_tar
        if os.path.exists(guidefwhmt_curr):
            guidefwhm_tar = guidefwhmt_curr.replace('/home/eng/','')
            tarfiles = tarfiles+' '+guidefwhm_tar
        if os.path.exists(guidext_curr):
            guidext_tar = guidext_curr.replace('/home/eng/','')
            tarfiles = tarfiles+' '+guidext_tar
        if os.path.exists(guideyt_curr):
            guideyt_tar = guideyt_curr.replace('/home/eng/','')
            tarfiles = tarfiles+' '+guideyt_tar
        if os.path.exists(guidexy_curr):
            guidexy_tar = guidexy_curr.replace('/home/eng/','')
            tarfiles = tarfiles+' '+guidexy_tar


    os.system('rm -f '+tarname)
    os.system('tar -cvf '+tarname+' '+tarfiles)

    os.system('rm -f '+acq_in_curr)
    if os.path.exists(guide_in_curr):
        os.system('rm -f '+guide_in_curr)

#file translations and path additions

    guide_jpg = guideimage.replace('.fits','.jpg')
    guide_wpath = data_dir+'/'+guideimage
    guidejpg_wpath = data_dir+'/'+guide_jpg


    if os.path.exists(new_dir+'/'+grpid+'_'+propid+'.html'):
        os.remove(new_dir+'/'+grpid+'_'+propid+'.html')
    outfile = open(new_dir+'/'+grpid+'_'+propid+'.html','w')
    name = grpid+'_'+propid    



    outfile.write('<html>\n')
    outfile.write('<head>\n')
    outfile.write('<title>'+name+'</title>\n')
    outfile.write('</head>\n')
    outfile.write('<center>\n')
    outfile.write('<h1><font COLOR=red> Quick summary <br> '+name+'</h1>\n')
    outfile.write('</center>\n')
    outfile.write('<br><br><br>\n')
    outfile.write('<table>\n')
    outfile.write('<tr ALIGN="center">\n')
    outfile.write('<td COLSPAN=1><b> Acquisition image </b></td>\n')
    outfile.write('<td COLSPAN=1><b> First Guide Image</b></td>\n')
    outfile.write('</tr>\n')
    outfile.write('<tr>\n')
    outfile.write('<th ROWSPAN=2>\n')
    outfile.write('<IMG src="'+acqjpg_in_curr+'" height="400">\n')
    outfile.write('</th>\n')
    outfile.write('<th ROWSPAN=2>\n')
    if guidejpg_in_curr != 'Null':
        outfile.write('<IMG src="'+guidejpg_in_curr+'" height="400">\n')
    outfile.write('</th>\n')
    outfile.write('</tr>\n')
    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    outfile.write('<td COLSPAN=1> <a href=./'+acq_in_curr+'> '+acq_in_curr+'</td>\n')
    if guide_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./'+guide_in_curr+'> '+guide_in_curr+'</td>')
    outfile.write('</tr>')
    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    if acqreg_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./'+acqlog_in_curr+'> '+acqlog_in_curr+'</td>\n')
    outfile.write('</tr>')
    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    if acqreg_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./'+acqreg_in_curr+'> '+acqreg_in_curr+'</td>\n')
    outfile.write('</tr>')

    outfile.write('<tr></tr><tr></tr>\n')
    outfile.write('<tr ALIGN="center">\n')
    if acqcat_in_curr != 'Null':
        outfile.write('<td COLSPAN=1> <a href=./'+acqcat_in_curr+'> '+acqcat_in_curr+'</td>\n')
    outfile.write('</tr>')

    outfile.write('</table>')

#search for guide plots corresponding to each science spec
    for _sciencespec in sciencespecs:
        science_name = _sciencespec.replace('/icc/tmp/','')
        guidecounts_pl = _sciencespec.replace('.fits','_guidecounts.png')
        guidecounts_pl = guidecounts_pl.replace('/icc/tmp/f_','/home/eng/f_')
        guidecounts_pl = guidecounts_pl.replace('/icc/tmp/g_','/home/eng/g_')
        guidecounts_curr = guidecounts_pl.replace('/home/eng/','')
        guidefwhmt_pl = guidecounts_pl.replace('guidecounts.png','guidefwhmt.png')
        guidext_pl = guidecounts_pl.replace('guidecounts.png','guidext.png')
        guideyt_pl = guidecounts_pl.replace('guidecounts.png','guideyt.png')
        guidexy_pl = guidecounts_pl.replace('guidecounts.png','guidexy.png')

        guidefwhmt_curr = guidecounts_curr.replace('guidecounts.png','guidefwhmt.png')
        guidext_curr = guidecounts_curr.replace('guidecounts.png','guidext.png')
        guideyt_curr = guidecounts_curr.replace('guidecounts.png','guideyt.png')
        guidexy_curr = guidecounts_curr.replace('guidecounts.png','guidexy.png')

        

        os.system('mv '+guidecounts_pl+' '+new_dir)
        os.system('mv '+guidefwhmt_pl+' '+new_dir)
        os.system('mv '+guidext_pl+' '+new_dir)
        os.system('mv '+guideyt_pl+' '+new_dir)
        os.system('mv '+guidexy_pl+' '+new_dir)

        

        outfile.write('<table> \n <tr ALIGN="center"> \n <td COLSPAN=1><b> Guide Plots For '+science_name+'  </b></td> \n </tr> \n')
        outfile.write('<tr>')
        outfile.write('<th ROWSPAN=2>')
#if guide
        outfile.write('<IMG src="'+guidecounts_curr+'" height="400">')
        outfile.write('</th>')
        outfile.write('<th ROWSPAN=2>')

        outfile.write('<IMG src="'+guidefwhmt_curr+'" height="400">')
        outfile.write('</th>')
        outfile.write('</tr>')
        outfile.write('<tr>')
        outfile.write('</tr>')
        outfile.write('</table>')

        outfile.write('<table>')
        outfile.write('<tr>')
        outfile.write('<th ROWSPAN=2>')
        outfile.write('<IMG src="'+guidexy_curr+'" height="300">')
        outfile.write('</th>')
        outfile.write('<th ROWSPAN=2>')

        outfile.write('<IMG src="'+guidext_curr+'" height="300">')
        outfile.write('</th>')
        outfile.write('<th ROWSPAN=2>')

        outfile.write('<IMG src="'+guideyt_curr+'" height="300">')
        outfile.write('</th>')
        outfile.write('</tr>')


        outfile.write('</table>')

    outfile.write('<table> \n <tr ALIGN="center"> \n <td COLSPAN=1><b> Spectra </b></td> \n </tr> \n')


    for _specs in specslist:
        _specs = _specs.replace('/icc/tmp/',data_dir+'/')
        _specsjpg = _specs.replace('.fits','.jpg')
        specindex = _specs.rfind('/')
        spec_in_curr = _specs[specindex+1:]
        specjpg_in_curr = spec_in_curr.replace('.fits','.jpg')
        try:
            os.remove(new_dir+'/'+spec_in_curr)
        except:
            pass
        try:
            os.remove(new_dir+'/'+specjpg_in_curr)
        except:
            pass

        os.system('ln -s '+_specs+' '+new_dir+'/')
        os.system('ln -s '+_specsjpg+' '+new_dir+'/')

        outfile.write('<tr>')
        outfile.write('<th ROWSPAN=2>')
        outfile.write('<IMG src="'+specjpg_in_curr+'" height="200">')
        outfile.write('</th>')
        outfile.write('</tr>')
        outfile.write('<tr></tr>')
        outfile.write('<tr></tr>')
        outfile.write('<tr ALIGN="center">')
        outfile.write('<td COLSPAN=1> <a href=./'+spec_in_curr+'> '+spec_in_curr+'</td>')
        outfile.write('</tr>')



    try:
        os.remove(new_dir+'/Null')
    except:
        pass



def readlogfile(logname):
    import csv
    """ read in the guide log file """
    f = open(logname)
    headerlines = 1
    nlines = sum(1 for line in f) - headerlines   #There is one header line in these guider log files.
    f.close()
    logf = open(logname,"rb")
    for line in range(headerlines): logf.next()  #skip header
    logreadin = csv.reader(logf,delimiter=" ",skipinitialspace=True)
    """ initialize arrays """
    time = nlines*[None]

    nrow=0
    for row in logreadin:
        time[nrow] = float(row[0])
        xguide[nrow] = float(row[5])
        yguide[nrow] = float(row[6])
        nrow = nrow +1

if __name__ == '__main__':
    from optparse import OptionParser
    from datetime import datetime
    parser = OptionParser()
    parser.add_option("-n", "--nightlist", dest="nightlist", help="Nights to aggregate data",default='now')
    parser.add_option('-t', '--tmp', dest="tmp_dir", help="Temporary directory", default="./")
    parser.add_option('-s', '--site', dest="site", help="site name, ogg or coj", default="floyds.coj.lco.gtn")


    # Add others
    (options, args) = parser.parse_args()
    site = options.site
    nightlist = options.nightlist
    if nightlist == 'now':
        yo = str(datetime.date(datetime.now()))
        yo = yo.replace('-','')
        nightlist = yo
        if site == 'floyds.ogg.lco.gtn':
            nightlist_four = nightlist[4:]
            nightlist_year = nightlist[:4]
            #print nightlist_four
            if nightlist_four == '0101':
                nightlist = nightlist_year+'1231'
            elif nightlist_four == '0201':
                nightlist = nightlist_year+'0131'
            elif nightlist_four == '0301':
                nightlist = nightlist_year+'0228'
            elif nightlist_four == '0401':
                nightlist = nightlist_year+'0331'
            elif nightlist_four == '0501':
                nightlist = nightlist_year+'0430'
            elif nightlist_four == '0601':
                nightlist = nightlist_year+'0531'
            elif nightlist_four == '0701':
                nightlist = nightlist_year+'0630'
            elif nightlist_four == '0801':
                nightlist = nightlist_year+'0731'
            elif nightlist_four == '0901':
                nightlist = nightlist_year+'0831'
            elif nightlist_four == '1001':
                nightlist = nightlist_year+'0930'
            elif nightlist_four == '1101':
                nightlist = nightlist_year+'1031'
            elif nightlist_four == '1201':
                nightlist = nightlist_year+'1130'
            else:
                nightlist = str(int(nightlist)-1)
    if nightlist == 'yesterday':
        yo = str(datetime.date(datetime.now()))

    tmp_dir = options.tmp_dir

    agg_floyds(nightlist,site=site,tmp_dir=tmp_dir)

#!/opt/epd/bin/python

import pyfits,os
#matplotlib.use('Agg')
from pylab import *

def mk_guideinfo_plots(rootname,utstart,utstop,logname):
    import scipy,pylab

#    matplotlib.use('Agg')
#read in the guider log file
    guideinfo = readlogfile(logname)


#    timeimg = scipy.array(guideinfo[0]); fileimg = scipy.array(guideinfo[1]); totcnts = scipy.array(guideinfo[9]); xguide = scipy.array(guideinfo[2]); yguide = scipy.array(guideinfo[3]); guideval = scipy.array(guideinfo[4]); xcen = scipy.array(guideinfo[5]); ycen = scipy.array(guideinfo[6]); fwhm = scipy.array(guideinfo[7]); status = scipy.array(guideinfo[8])
    timeimg = scipy.array(guideinfo[0]); fileimg = scipy.array(guideinfo[1]); totcnts = scipy.array(guideinfo[7]); guideval = scipy.array(guideinfo[2]); xcen = scipy.array(guideinfo[3]); ycen = scipy.array(guideinfo[4]); fwhm = scipy.array(guideinfo[5]); status = scipy.array(guideinfo[6])


#pull out all lines from the guide file that happened between utstart and utstop

    #print utstart
    #print 'hey hey hey'
    #print utstop
    time_cull = timeimg * (timeimg > float(utstart)) * (timeimg < float(utstop))
    timeimg_gd = timeimg[where(time_cull != 0)]
    fileimg_gd = fileimg[where(time_cull != 0)]
    totcnts_gd = totcnts[where(time_cull != 0)]
    #xguide_gd = xguide[where(time_cull != 0)]
    #yguide_gd = yguide[where(time_cull != 0)]
    xcen_gd = xcen[where(time_cull != 0)]          
    ycen_gd = ycen[where(time_cull != 0)]
    fwhm_gd = fwhm[where(time_cull != 0)]
    #print totcnts_gd

#read in all the guide images that are within the appropriate time range.  Read their MJD's
    #print fileimg_gd
    #print timeimg_gd
    #print fileimg_gd.size
    if fileimg_gd.size > 0:
        pass
    else: 
        print 'No guide frames taken during exposure!'


    if fileimg_gd.size > 0:
        guidecheck = 1
        for _fileimg in fileimg_gd:
            _fileimg = _fileimg+'.fits'
            headlist = pyfits.open(_fileimg)
            headnow = headlist[0].header
            
            mjdnow = headnow['MJD']
            try:
                mjd_list.append(float(mjdnow))
            except NameError:
                mjd_list = [float(mjdnow)]
            
            
        mjd_arr = array(mjd_list)
        mjd_arr = mjd_arr - (min(mjd_arr)-0.1/60.0/24.0)
        mjd_arr_sec = mjd_arr*24.0*60.0*60.0
    #print mjd_list
        plot_guideflux(mjd_arr_sec,totcnts_gd)
        guidecountplotname = rootname+'_guidecounts.png'
        os.system('mv yo.png '+guidecountplotname)

        plot_guidepos(xcen_gd,ycen_gd)
        guidexypos_plotname = rootname+'_guidexy.png'
        os.system('mv yo_xy.png '+guidexypos_plotname)
        
        plot_xposwtime(mjd_arr_sec,xcen_gd)
        guidext_plotname = rootname+'_guidext.png'
        os.system('mv yo_xt.png '+guidext_plotname)
        
        plot_yposwtime(mjd_arr_sec,ycen_gd)
        guideyt_plotname = rootname+'_guideyt.png'
        os.system('mv yo_yt.png '+guideyt_plotname)
        
        plot_fwhmwtime(mjd_arr_sec,fwhm_gd)
        guidefwhmt_plotname = rootname+'_guidefwhmt.png'
        os.system('mv yo_fwhm.png '+guidefwhmt_plotname)


#remake all plots, put on one pdf page
#        plcolor=True
#        ioff()
#mp = subplot(mjd

    else:
        guidecheck = 0

def readlogfile(logname):
    import csv
    """ read in the guide log file """
    f = open(logname)
    headerlines = 1
    nlines = sum(1 for line in f) - headerlines   #There is one header line in these guider log files.
    f.close()
    logf = open(logname,"rb")
    for line in range(headerlines): logf.next()  #skip header
    logreadin = csv.reader(logf,delimiter="\t",skipinitialspace=True)
    """ initialize arrays """
    timeimg = nlines*[None]
    fileimg = nlines*[None]
    #xguide = nlines*[None]
    #yguide = nlines*[None]
    guideval = nlines*[None]
    xcen = nlines*[None]
    ycen = nlines*[None]
    fwhm = nlines*[None]
    status = nlines*[None]
    totcounts = nlines*[None]

    nrow=0
    for row in logreadin:
        timestr = str(row[0])
        timestr = timestr[8:]
        timeimg[nrow] = float(timestr)
        fileimg[nrow] = str(row[1])
        try:
            #if isinstance(row[5],str):
                #print row[5]
                #print row[6]
                #print row[7]
            guideval[nrow] = str(row[5])
            xcen[nrow] = float(row[6])
            ycen[nrow] = float(row[7])
            fwhm[nrow] = float(row[8])
            status[nrow] = str(row[9])
            totcounts[nrow] = float(row[10])
            
        except ValueError:
            guideval[nrow] = str(row[7])
            xcen[nrow] = float(row[8])
            ycen[nrow] = float(row[9])
            fwhm[nrow] = float(row[10])
            status[nrow] = str(row[11])
            totcounts[nrow] = float(row[12])
        nrow = nrow +1

#    return [timeimg,fileimg,xguide,yguide,guideval,xcen,ycen,fwhm,status,totcounts]
    return [timeimg,fileimg,guideval,xcen,ycen,fwhm,status,totcounts]


def plot_guideflux(mjd,totcounts):

#    print 'yo'
    plcolor=True
    ioff()
    p = plot(mjd,totcounts,'bo')
    xlabel('Seconds')
    ylabel('Total Counts')
    #show()
    #SetPlot()
    #raw_input('press enter to close')
    os.system('rm -f yo.png')
    savefig('yo.png')
    close()

def plot_xposwtime(mjd,guidecoordx):
    plcolor=True
    ioff()
    pl_xy = plot(mjd,guidecoordx,'bo')
    xlabel('Seconds')
    ylabel('xpos (pix)')
    #show()
    #SetPlot()
    #raw_input('press enter to close')
    os.system('rm -f yo_xt.png')
    savefig('yo_xt.png')

    close()

def plot_yposwtime(mjd,guidecoordy):
    plcolor=True
    ioff()
    pl_xy = plot(mjd,guidecoordy,'bo')
    xlabel('Seconds')
    ylabel('ypos (pix)')
    #show()
    #SetPlot()
    #raw_input('press enter to close')
    os.system('rm -f yo_yt.png')
    savefig('yo_yt.png')

    close()

def plot_fwhmwtime(mjd,fwhm):
    plcolor=True
    ioff()
    pl_fwhm = plot(mjd,fwhm,'bo')
    xlabel('Seconds')
    ylabel('FWHM (pix)')
    #show()
    #SetPlot()
    #raw_input('press enter to close')
    os.system('rm -f yo_fwhm.png')
    savefig('yo_fwhm.png')

    close()





def plot_guidepos(xguide,yguide):
    plcolor=True
    ioff()
    pl_xy = plot(xguide,yguide,'bo')
    xlabel('xpos (pix)')
    ylabel('ypos (pix)')
    #show()
    #SetPlot()
    #raw_input('press enter to close')
    os.system('rm -f yo_xy.png')
    savefig('yo_xy.png')
 
    close()


def SetPlot():           ############################ set plot defaults
    yticklabels = getp(gca(),'yticklabels')
    xticklabels = getp(gca(),'xticklabels')
    setp(xticklabels,fontsize='20')   
    setp(yticklabels,fontsize='20')   
    
    legend(numpoints=1,markerscale=1.5)
    leg = gca().get_legend()
    ltext = leg.get_texts()
    setp(ltext,fontsize=15)
    
def SY(plcolor):  ###################################    define plot symbol
    sym = 'osDd1234hHpx+<>^v' # simboli
    col = 'bgrcmk'            # colori
    _sy = []
    if plcolor: colrange = col
    else:     colrange = 'k'
    for s in sym:
        for c in colrange:     
            _sy.append(c+s)
    return _sy


if __name__ == '__main__':
    from optparse import OptionParser
    from datetime import datetime
    parser = OptionParser()

    parser.add_option("-l", "--logname", dest="logname",help="Guider log file",default="test.log")
    parser.add_option("-b", "--utstart", dest="utstart",help="UT Start time of interest",default="10000.")
    parser.add_option("-e", "--utstop", dest="utstop",help="UT Stop time of interest",default="10001.")

    (options, args) = parser.parse_args()
    logname = options.logname
    utstart = options.utstart
    utstop = options.utstop

    mk_guideinfo_plots(utstart,utstop,logname)

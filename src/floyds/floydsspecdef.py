from pathlib2 import Path

def rebin(y, s):
    from numpy import array
    # s=raw_input('which value for rebining ? [1]')
    if not s:
        s = 1
    else:
        s = int(s)
    ytemp = array(y)
    yreb = []
    for jj in range(len(ytemp)):
        yreb.append(sum(ytemp[jj: jj + s]) / len(ytemp[jj: jj + s]))
    return yreb


def correctfringing_auto(img0, imgN):
    '''correct for fringing automatically 
       giving as input the 1D spectra: object and flat '''
    import floyds
    from astropy.io import fits
    import string
    from numpy import interp as ninterp
    from numpy import trapz, argmin, array, float32
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.continuum']
    for t in toforget: iraf.unlearn(t)
    floyds.util.delete('sky1.fits')
    floyds.util.delete('tsky1.fits')
    floyds.util.delete('tsky11.fits')
    floyds.util.delete('sky2.fits')
    iraf.scopy(img0, 'sky1.fits', w1=7000, w2=10000)
    iraf.scopy(imgN, 'sky2.fits', w1=7000, w2=10000)
    # 
    iraf.specred.continuum('sky1.fits', output='tsky1.fits', line='*', type='difference',
                           interact='yes', function='spline3', niterat=100, low_rej=4, high_re=5,
                           sample='7000:7550,7730:10000', order=10, ask='no')
    iraf.specred.continuum('sky1.fits', output='tsky11.fits', type='fit', interact='yes', function='spline3',
                           niterat=100, low_rej=4, high_re=5, sample='7000:7550,7730:10000', order=10, ask='no')
    #   sky1   standard
    #   tsky1  fringing on the standard 
    #   tsky11 continuum  on the standard 
    #   sky2   fringing
    #
    xx1, yy1 = floyds.util.readspectrum('tsky1.fits')  #  fringing on the standrad
    xx11, yy11 = floyds.util.readspectrum('sky1.fits')  #  standard
    xx111, yy111 = floyds.util.readspectrum('tsky11.fits')  #  continuum
    xx2, yy2 = floyds.util.readspectrum('sky2.fits')  #  fringing form flat
    # shift fringing on std and fringing from flat
    shift = floyds.floydsspecdef.checkwavelength_arc(xx1, yy1 - min(yy1), xx2, yy2, '', '', False) * (-1)
    if abs(shift) > 15:
        shift = 0
    xx3 = array(xx2, float)[:]
    yy3 = yy2[:]
    xx4 = xx3 + shift
    yy3 = ninterp(xx1, xx4, yy3)
    _scaleo2 = []
    integral_o2 = []
    for i in range(1, 81):
        j = 0.6 + i * 0.04
        ll = abs((yy11 / ((yy3 * j) - (j - 1))) - yy111)  #  standard/(fringing*J)
        integraleo2 = trapz(ll, xx11)
        integral_o2.append(integraleo2)
        _scaleo2.append(j)
    so2 = _scaleo2[argmin(integral_o2)]
    floyds.util.delete('sky1.fits')
    floyds.util.delete('tsky1.fits')
    floyds.util.delete('tsky11.fits')
    floyds.util.delete('sky2.fits')
    data1, hdr = fits.getdata(img0, 0, header=True)
    if 'GRISM' in hdr:
        arm = hdr['GRISM']
    else:
        arm = ''
    x1, y1 = floyds.util.readspectrum(img0)
    xn, yn = floyds.util.readspectrum(imgN)
    xx4 = xn + shift
    yy3 = ninterp(x1, xx4, yn)
    yynew = y1 / ((yy3 * so2) - (so2 - 1))
    data1[0] = array(yynew)
    floyds.util.delete(img0)
    fits.writeto(img0, float32(data1), hdr)
    floyds.util.updateheader(img0, 0, {'FRSCALE': [so2, 'fringing scale factor ']})
    floyds.util.updateheader(img0, 0, {'FRSHIFT': [shift, 'fringing shift factor ']})
    floyds.util.updateheader(img0, 0, {'FLAT' + arm: [string.split(imgN, '/')[-1], 'flat field file ']})
    return img0, so2, shift


def normflat(img):
    import floyds
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.continuum']
    for t in toforget:
        iraf.unlearn(t)
    floyds.util.delete('sky1.fits')
    floyds.util.delete('sky2.fits')
    floyds.util.delete('tsky1.fits')
    floyds.util.delete('tsky2.fits')
    floyds.util.delete('N' + img)
    floyds.util.delete('M' + img)
    _order1 = 50
    _order2 = 30
    iraf.scopy(img, 'sky1.fits', w1='INDEF', w2=6750)
    iraf.scopy(img, 'sky2.fits', w1=6700, w2='INDEF')
    iraf.specred.continuum('sky1.fits', output='tsky1.fits', type='fit', interact='no', function='spline3',
                           niterat=100, low_rej=0, high_re=0, sample='*', order=_order1, ask='YES')
    iraf.specred.continuum('sky2.fits', output='tsky2.fits', type='fit', interact='no', function='chebyshev',
                           niterat=100, low_rej=0, high_re=0, sample='*', order=_order2, ask='YES')
    sss = iraf.specred.scombine('tsky1.fits,tsky2.fits', output='M' + img, combine='average', reject='none',
                                scale='none', weight='none', Stdout=1)
    iraf.specred.sarith(img, op='/', input2='M' + img, output='N' + img, format='multispec')
    floyds.util.delete('sky1.fits')
    floyds.util.delete('sky2.fits')
    floyds.util.delete('tsky1.fits')
    floyds.util.delete('tsky2.fits')
    floyds.util.delete('M' + img)
    return 'N' + img


def correctfringing(imgex, fringingmask):
    import floyds
    from pyraf import iraf
    from astropy.io import fits
    from numpy import float32

    iraf.images(_doprint=0)
    xx1, yy1 = floyds.util.readspectrum(imgex)
    xx2, yy2 = floyds.util.readspectrum(fringingmask)
    from pylab import plot, ion, clf, legend

    ion()
    clf()
    _shift = 0
    _imge = '_tmpimg.fits'
    plot(xx1, yy1, color='red', label='spectrum')
    legend(numpoints=1, markerscale=1.5)
    from numpy import array

    xx3 = array(xx2, float)[:]
    yy3 = yy2[:]
    yynew = yy1 / yy3
    from numpy import interp as ninterp

    while _shift not in ['n', 'no']:
        lines = plot(xx1, yynew, color='blue', label='fringing corrected cor.')
        legend(numpoints=1, markerscale=1.5)
        _shift = raw_input('\n### shift the correction (1,1.2,1.8,2.0,....)   [[n]/num] ? ')
        if not _shift:
            _shift = 'n'
        try:
            _shift = float(_shift)
            if _shift > 0.0:
                xx4 = xx3 + _shift
                yy3 = ninterp(xx2, xx4, yy3)
            else:
                xx4 = xx3 + _shift
                yy3 = ninterp(xx2, xx4, yy3)
            yynew = yy1 / yy3
        except:
            yynew = yy1 / yy3
        lines.pop(0).remove()
    _scale = 1
    while _scale not in ['n', 'no']:
        lines = plot(xx1, yynew, color='green', label='fringing corrected cor.')
        legend(numpoints=1, markerscale=1.5)
        _scale = raw_input('\n### scale the correction (1,1.2,1.8,2.0,....)   [[n]/num] ? ')
        if not _scale:
            _scale = 'n'
        try:
            _scale = float(_scale)
            yynew = yy1 / ((yy3 * _scale) - (_scale - 1))
        except:
            yynew = yy1 / yy3
        lines.pop(0).remove()
    data1, hdr = fits.getdata(imgex, 0, header=True)
    data1[0] = array(yynew)
    floyds.util.delete(imgex)
    fits.writeto(imgex, float32(data1), hdr)
    return imgex


def sensfunction(standardfile, _outputsens, _function, _order, _interactive, sample='*', fts_contaminated=False):
    import re
    import os
    import sys
    import string
    import floyds
    import datetime
    # import matplotlib
    MJDtoday = floyds.util.mjdtoday()
    from floyds.util import delete, readhdr, readkey3, updateheader, name_duplicate
    from astropy.io import fits
    from numpy import float32
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.scopy', 'specred.sensfunc', 'specred.standard']
    for t in toforget:
        iraf.unlearn(t)
    iraf.specred.scopy.format = 'multispec'
    iraf.specred.verbose = 'no'
    iraf.set(direc=floyds.__path__[0] + '/')
    _caldir = 'direc$standard/MAB/'
    _extinctdir = 'direc$standard/extinction/'
    print standardfile, _outputsens
    _outputsens = name_duplicate(standardfile, _outputsens, '')
    if os.path.isfile(_outputsens):
        if _interactive.lower() != 'yes':
            floyds.util.delete(_outputsens)
        else:
            answ = raw_input('sensitivity function already computed, do you want to do it again [[y]/n] ? ')
            if not answ:
                answ = 'y'
            if answ.lower() in ['y', 'yes']:
                floyds.util.delete(_outputsens)

    if not os.path.isfile(_outputsens):
        hdrs = readhdr(standardfile)
        _airmass = readkey3(hdrs, 'airmass')
        _exptime = readkey3(hdrs, 'exptime')
        _tel = readkey3(hdrs, 'TELID')
        if _tel not in ['ftn', 'fts']:
            _tel = readkey3(hdrs, 'SITEID')
        if _tel in ['fts', 'coj']:  # to be checked
            _extinction = 'ssoextinct.dat'
            _observatory = 'sso'
        elif _tel in ['ftn', 'ogg']:
            _extinction = 'maua.dat'
            _observatory = 'cfht'
        else:
            sys.exit('ERROR: observatory not recognised')

        refstar = 'm' + re.sub('.dat', '', fits.open(standardfile)[0].header['stdname'])
        if sample == '*':
            _outputstd = 'std_' + str(readkey3(hdrs, 'grism')) + '_' + str(readkey3(hdrs, 'filter')) + '.fits'
            floyds.util.delete(_outputstd)
            floyds.util.delete(_outputsens)
            iraf.specred.standard(input=standardfile, output=_outputstd, extinct=_extinctdir + _extinction,
                                  caldir=_caldir, observa=_observatory, star_nam=refstar, airmass=_airmass,
                                  exptime=_exptime, interac=_interactive)

            if str(readkey3(hdrs, 'grism')) == 'blu':
                floyds.floydsspecdef.cutstd(_outputstd, start=5980, end=6012, out=False)

            iraf.specred.sensfunc(standard=_outputstd, sensitiv=_outputsens, extinct=_extinctdir + _extinction,
                                  ignorea='yes', observa=_observatory, graphs='sri', functio=_function, order=_order,
                                  interac=_interactive)
        else:
            _outputstd = 'std_' + str(readkey3(hdrs, 'grism')) + '_' + str(readkey3(hdrs, 'filter')) + '.fits'
            floyds.util.delete(_outputstd)
            iraf.specred.standard(input=standardfile, output=_outputstd, extinct=_extinctdir + _extinction,
                                  caldir=_caldir, observa=_observatory, star_nam=refstar, airmass=_airmass,
                                  exptime=_exptime, interac=_interactive)
            #     sens full range
            if str(readkey3(hdrs, 'grism')) == 'blu' and _tel in ['fts', 'coj'] and fts_contaminated:
                print 'split blue sens'
                sss=''
                for jj in range(0, len(string.split(sample, ','))):
                    aa, bb = string.split(string.split(sample, ',')[jj], ':')
                    std0 = '_stdcut' + str(jj) + '.fits'
                    std1 = '_stdcut' + str(jj)
                    sens1 = '_senscut' + str(jj) + '.fits'
                    floyds.util.delete(std1 + ',' + sens1 + ',' + std0)
                    iraf.scopy(standardfile, std0, w1=aa, w2=bb)
                    iraf.specred.standard(input=std0, output=std1, extinct=_extinctdir + _extinction, caldir=_caldir,
                                      observa=_observatory, star_nam=refstar, airmass=_airmass,
                                      exptime=_exptime, interac=_interactive)
                    iraf.specred.sensfunc(standard=std1, sensitiv=sens1, extinct=_extinctdir + _extinction, ignorea='yes',
                                      observa=_observatory, graphs='sri', functio=_function, order=100, interac=_interactive)
                    if sss:
                        sss = sss + ',' + sens1
                    else:
                        sss = sens1

                    floyds.util.delete(std1 + ',' + std0)
                _outputsens = combineredsens(sss, _outputsens)
                for i in string.split(sss, ','):
                    floyds.util.delete(i)
            elif str(readkey3(hdrs, 'grism')) == 'blu':
                floyds.floydsspecdef.cutstd(_outputstd, start=3600, end=3900, out=False)
                floyds.floydsspecdef.cutstd(_outputstd, start=4300, end=4600, out=False)
                floyds.floydsspecdef.cutstd(_outputstd, start=5990, end=6012, out=False)

                floyds.util.delete('sens0.fits')
                iraf.specred.sensfunc(standard=_outputstd, sensitiv='sens0.fits', extinct=_extinctdir + _extinction,
                                  ignorea='yes', observa=_observatory, graphs='sri', functio=_function, order=_order,
                                  interac=_interactive)
                sss = 'sens0.fits'
                for jj in range(0, len(string.split(sample, ','))):
                    aa, bb = string.split(string.split(sample, ',')[jj], ':')
                    std0 = '_stdcut' + str(jj) + '.fits'
                    std1 = '_stdcut' + str(jj)
                    sens1 = '_senscut' + str(jj) + '.fits'
                    floyds.util.delete(std1 + ',' + sens1 + ',' + std0)
                    iraf.scopy(standardfile, std0, w1=aa, w2=bb)
                    iraf.specred.standard(input=std0, output=std1, extinct=_extinctdir + _extinction, caldir=_caldir,
                                      observa=_observatory, star_nam=refstar, airmass=_airmass,
                                      exptime=_exptime, interac=_interactive)
                    iraf.specred.sensfunc(standard=std1, sensitiv=sens1, extinct=_extinctdir + _extinction, ignorea='yes',
                                      observa=_observatory, graphs='sri', functio=_function, order=30, interac=_interactive)
                    sss = sss + ',' + sens1
                    floyds.util.delete(std1 + ',' + std0)
                _outputsens = combineblusens(sss, _outputsens)
                print _outputsens
                for i in string.split(sss, ','):
                    floyds.util.delete(i)
            else: # red arm
                print 'split red sens'
                sss=''
                for jj in range(0, len(string.split(sample, ','))):
                    aa, bb = string.split(string.split(sample, ',')[jj], ':')
                    std0 = '_stdcut' + str(jj) + '.fits'
                    std1 = '_stdcut' + str(jj)
                    sens1 = '_senscut' + str(jj) + '.fits'
                    floyds.util.delete(std1 + ',' + sens1 + ',' + std0)
                    iraf.scopy(standardfile, std0, w1=aa, w2=bb)
                    iraf.specred.standard(input=std0, output=std1, extinct=_extinctdir + _extinction, caldir=_caldir,
                                      observa=_observatory, star_nam=refstar, airmass=_airmass,
                                      exptime=_exptime, interac=_interactive)
                    iraf.specred.sensfunc(standard=std1, sensitiv=sens1, extinct=_extinctdir + _extinction, ignorea='yes',
                                      observa=_observatory, graphs='sri', functio=_function, order=30, interac=_interactive)
                    if sss:
                        sss = sss + ',' + sens1
                    else:
                        sss = sens1

                    floyds.util.delete(std1 + ',' + std0)
                _outputsens = combineredsens(sss, _outputsens)
                for i in string.split(sss, ','):
                    floyds.util.delete(i)

        hdr = fits.getheader(standardfile)
        data1, hdr1 = fits.getdata(_outputsens, 0, header=True)
        for key in ['ctype1', 'crval1', 'crpix1', 'cdelt1', 'cd1_1']:
            hdr[key] = hdr1[key] # keep wavelength calibration from sensfuncs
        floyds.util.delete(_outputsens)
        floyds.util.delete(_outputstd)
        fits.writeto(_outputsens, float32(data1), hdr)
    return _outputsens


def telluric_atmo(imgstd):
    from numpy import interp as ninterp
    from numpy import trapz, compress, argmin, array, float32
    from pyraf import iraf
    from astropy.io import fits

    iraf.images(_doprint=0)
    iraf.noao(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.onedspec(_doprint=0)
    toforget = ['imfilter.gauss', 'specred.apall', 'longslit.identify', 'longslit.reidentify', 'specred.standard',
                'onedspec.wspectext']
    for t in toforget:
        iraf.unlearn(t)
    import floyds
    from floyds.util import readhdr, readkey3, readspectrum

    _grism = readkey3(readhdr(imgstd), 'grism')
    _tel = readkey3(readhdr(imgstd), 'TELID')
    # if _tel not in ['ftn','fts']:    _tel=readkey3(readhdr(imgstd),'SIDEID')
    imgout = 'invers_atmo_' + imgstd
    floyds.util.delete(imgout)
    iraf.set(direc=floyds.__path__[0] + '/')
    _cursor = 'direc$standard/ident/cursor_sky_0'
    #   iraf.noao.onedspec.bplot(imgstd, cursor=_cursor, spec2= imgstd, new_ima=imgout, overwri='yes')
    imgout = floyds.floydsspecdef.atmofile(imgstd, imgout)
    xxstd, ffstd = readspectrum(imgout)
    if _grism in ['red']:
        llo2 = compress((array(xxstd) >= 7550) & (array(xxstd) <= 7750), array(xxstd))
        llh2o = compress((array(xxstd) >= 7100) & (array(xxstd) <= 7500), array(xxstd))
        ffo2 = compress((array(xxstd) >= 7550) & (array(xxstd) <= 7750), array(ffstd))
        ffh2o = compress((array(xxstd) >= 7100) & (array(xxstd) <= 7500), array(ffstd))
    if _grism in ['red']:
        _skyfileh2o = 'direc$standard/ident/ATLAS_H2O.fits'
        _skyfileo2 = 'direc$standard/ident/ATLAS_O2.fits'
        atlas_smooto2 = '_atlas_smoot_o2.fits'
        atlas_smooth2o = '_atlas_smoot_h2o.fits'
        _sigma = 200
        floyds.util.delete(atlas_smooto2)
        floyds.util.delete(atlas_smooth2o)
        iraf.imfilter.gauss(_skyfileh2o, output=atlas_smooth2o, sigma=_sigma)
        iraf.imfilter.gauss(_skyfileo2, output=atlas_smooto2, sigma=_sigma)
        llskyh2o, ffskyh2o = readspectrum(atlas_smooth2o)
        llskyo2, ffskyo2 = readspectrum(atlas_smooto2)
        ffskyo2cut = ninterp(llo2, llskyo2, ffskyo2)
        ffskyh2ocut = ninterp(llh2o, llskyh2o, ffskyh2o)
        _scaleh2o = []
        integral_h2o = []
        for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyh2ocut = list((array(ffskyh2ocut) * j) + 1 - j)
            diff_h2o = abs(_ffskyh2ocut - ffh2o)
            integraleh2o = trapz(diff_h2o, llh2o)
            integral_h2o.append(integraleh2o)
            _scaleh2o.append(j)
        _scaleo2 = []
        integral_o2 = []
        for i in range(1, 21):
            j = 0.6 + i * 0.04
            _ffskyo2cut = list((array(ffskyo2cut) * j) + 1 - j)
            diff_o2 = abs(_ffskyo2cut - ffo2)
            integraleo2 = trapz(diff_o2, llo2)
            integral_o2.append(integraleo2)
            _scaleo2.append(j)
        sh2o = _scaleh2o[argmin(integral_h2o)]
        so2 = _scaleo2[argmin(integral_o2)]
        telluric_features = ((array(ffskyh2o) * sh2o) + 1 - sh2o) + ((array(ffskyo2) * so2) + 1 - so2) - 1
        telluric_features = array([1] + list(telluric_features) + [1])
        llskyo2 = array([1000] + list(llskyo2) + [15000])
        telluric_features_cut = ninterp(xxstd, llskyo2, telluric_features)

        _imgout = 'atmo_' + _tel + '_' + imgstd

        data1, hdr = fits.getdata(imgstd, 0, header=True)
        data1[0] = array(telluric_features_cut)
        data1[1] = data1[1] / data1[1]
        data1[2] = data1[2] / data1[2]
        data1[3] = data1[3] / data1[3]
        floyds.util.delete(_imgout)
        fits.writeto(_imgout, float32(data1), hdr)
        floyds.util.delete(atlas_smooto2)
        floyds.util.delete(atlas_smooth2o)
        floyds.util.delete(imgout)
    else:
        _imgout = ''
        print '### telluric correction with model not possible '
    return _imgout


def checkwavestd(imgex, _interactive, _type=1):
    import floyds
    from astropy.io import fits
    from numpy import arange, array

    print '\n### Warning: check in wavelength with sky lines not performed\n'
    if _interactive.upper() in ['YES', 'Y']:
        answ = raw_input('\n### Do you want to check the wavelength calibration with telluric lines [[y]/n]? ')
        if not answ: answ = 'y'
    else:
        answ = 'y'
    if answ in ['y', 'yes']:
        print '\n### check wavelength calibration with telluric lines \n'
        _skyfile = floyds.__path__[0] + '/standard/ident/sky_new_0.fits'
        skydata, skyhdr = fits.getdata(_skyfile, header=True)
        skyff = 1 - skydata
        crval1 = skyhdr['CRVAL1']
        cd1 = skyhdr['CD1_1']
        skyxx = arange(len(skyff))
        skyaa = crval1 + skyxx * cd1
        _tel = fits.getheader(imgex)['TELID']
        atmofile = floyds.floydsspecdef.atmofile(imgex, 'atmo2_' + _tel + '_' + imgex)
        atmodata, atmohdr = fits.getdata(atmofile, header=True)
        if _type == 1:
            atmoff = 1 - atmodata[0][0]
        else:
            atmoff = 1 - atmodata
        crval1 = atmohdr['CRVAL1']
        cd1 = atmohdr['CD1_1']
        atmoxx = arange(len(atmoff))
        atmoaa = crval1 + atmoxx * cd1
        shift = floyds.floydsspecdef.checkwavelength_arc(atmoaa, atmoff, skyaa, skyff, 6800, 7800, _interactive)
        floyds.util.delete('atmo2_' + _tel + '_' + imgex)
    else:
        shift = 0
    return shift

def checkwavelength_obj(fitsfile, skyfile, _interactive='yes', usethirdlayer=True):
    import floyds
    from astropy.io import fits
    import numpy as np
    from pyraf import iraf

    if _interactive.lower() in ['yes', 'y']:
        do_shift = raw_input('### Do you want to check the wavelength calibration with telluric lines? [[y]/n] ')
    else:
        print '### Checking wavelength calibration with telluric lines'
        do_shift = ''
    if do_shift != 'n':
        if usethirdlayer:
            iraf.scopy(fitsfile + '[*,1,3]', 'skylayer.fits') # iraf.continuum doesn't allow slices
            subtracted = floyds.floydsspecdef.continumsub('skylayer.fits', 6, 1)
            floyds.util.delete('skylayer.fits')
        else:
            subtracted = floyds.floydsspecdef.continumsub(fitsfile, 6, 1)
        sky_spec = fits.open(subtracted)[0]
        y1 = sky_spec.data
        crval1 = sky_spec.header['CRVAL1']
        x1 = crval1 + np.arange(len(y1)) * sky_spec.header['CD1_1']
        sky_arch = fits.open(skyfile)[0]
        y2 = sky_arch.data
        x2 = sky_arch.header['CRVAL1'] + np.arange(len(y2)) * sky_arch.header['CD1_1']
        shift = checkwavelength_arc(x1, y1, x2, y2, 5500, 6500, _interactive)
        if _interactive.lower() in ['yes', 'y']:
            answ = raw_input('By how much do you want to shift the wavelength calibration? [{}] '.format(shift))
            if answ:
                shift = float(answ)
        arm = sky_spec.header['GRISM'].upper()
        floyds.util.updateheader(fitsfile, 0, {'CRVAL1': (crval1 + shift, ''), 'SHIFT'+arm: (shift, '')})
        floyds.util.delete(subtracted)
    else:
        shift = 0
    return shift

############################################

def checkwavelength_arc(xx1, yy1, xx2, yy2, xmin, xmax, _interactive='yes'):
    from numpy import array, trapz, compress
    from numpy import interp as ninterp

    minimo = max(min(xx1), min(xx2)) + 60
    massimo = min(max(xx1), max(xx2)) - 60
    yy1 = [0 if e < 0 else e for e in array(yy1)]
    yy2 = [0 if e < 0 else e for e in array(yy2)]
    _shift, integral = [], []
    for shift in range(-600, 600, 1):
        xxnew = xx1 + shift / 10.
        yy2interp = ninterp(xxnew, xx2, yy2)
        yy2timesyy = yy2interp * yy1
        xxcut = compress((array(xxnew) >= minimo) & (array(xxnew) <= massimo), array(xxnew))
        yycut = compress((array(xxnew) >= minimo) & (array(xxnew) <= massimo), array(yy2timesyy))
        integrale = trapz(yycut, xxcut)
        integral.append(integrale)
        _shift.append(shift / 10.)
    result = _shift[integral.index(max(integral))]
    if _interactive in ['YES', 'y', 'Y', 'yes', 'Yes', True]:
        #   import matplotlib as mpl  
        #   mpl.use("TKAgg")  
        from pylab import plot, show, ion, clf, legend, xlim, ylim

        ion()
        clf()
        ratio = trapz(yy1, xx1) / trapz(yy2, xx2)
        yy3 = array(yy2) * float(ratio)
        xx4 = xx1 + result
        plot(xx1, yy1, label='spectrum')
        plot(xx2, yy3, label='reference sky')
        plot(xx4, yy1, label='shifted spectrum')
        legend(numpoints=1, markerscale=1.5)
        if xmin != '' and xmax != '':
            xlim(xmin, xmax)
    return result


def atmofile(imgstd, imgout=''):
    import os
    from astropy.io import fits
    import numpy as np

    data, hdr = fits.getdata(imgstd, 0, header=True)
    if len(data) < 10:
        yy = data[0][0]
    else:
        yy = data
    #yy=data
    crvals = hdr['CRVAL1']
    cds = hdr['CD1_1']
    xx = np.arange(len(yy))
    aasens = crvals + (xx) * cds
    aasens1 = np.compress(((aasens < 6829) | (aasens > 7090)) & ((aasens < 7140) | (aasens > 7420)) & \
                          ((aasens < 7570) | (aasens > 7750)) & ((aasens < 7820) | (aasens > 8450)) & \
                          ((aasens < 8910) | (aasens > 9225)) & ((aasens < 9267) | (aasens > 9890)), aasens)
    yy1 = np.compress(((aasens < 6829) | (aasens > 7090)) & ((aasens < 7140) | (aasens > 7420)) & \
                      ((aasens < 7570) | (aasens > 7750)) & ((aasens < 7820) | (aasens > 8450)) & \
                      ((aasens < 8910) | (aasens > 9225)) & ((aasens < 9267) | (aasens > 9890)), yy)
    y11 = np.interp(aasens, aasens1, yy1)
    if not imgout:   imgout = 'atmo_' + imgstd
    os.system('rm -rf ' + imgout)
    #    if len(data)<10:   data=np.float32(yy/y11)
    #    else:              data=np.float32(yy/y11)
    mask = y11 == 0
    y11[mask] = 1
    data = np.float32(yy / y11)
    #    mask = np.isnan(data)     # masking   nan values 
    #    data[mask]=1
    fits.writeto(imgout, data, hdr)
    return imgout


def extractspectrum(img, dv, _ext_trace, _dispersionline, _interactive, _type, automaticex=False):
    import glob, os, string, sys, re
    from floyds.util import delete, readhdr, readkey3, dvex, updateheader
    import floyds
    import datetime
    from numpy import arange, float32
    from astropy.io import fits
    MJDtoday = floyds.util.mjdtoday()
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.apall', 'specred.transform']
    for t in toforget: iraf.unlearn(t)

    dv = dvex()
    hdr = readhdr(img)
    _gain = readkey3(hdr, 'gain')
    _rdnoise = readkey3(hdr, 'ron')
    _grism = readkey3(hdr, 'grism')
    iraf.specred.dispaxi = 1
    imgex = re.sub('.fits', '_ex.fits', img)
    imgfast = re.sub(str(MJDtoday) + '_', '', img)
    if _type=='agn' or (not os.path.isfile(imgex)
        and not os.path.isfile('database/ap' + re.sub('.fits', '', img))
        and not os.path.isfile('database/ap' + re.sub('.fits', '', imgfast))):
        _new = 'yes'
        _extract = 'yes'
    else:
        if automaticex:
            if _interactive.upper() in ['YES', 'Y']:
                answ = 'x'
                while answ not in ['o', 'n', 's']:
                    answ = raw_input('\n### New extraction [n], extraction with old parameters [o], skip extraction [s] ? [o] ')
                    if not answ: answ = 'o'
                if answ == 'o':
                    _new, _extract = 'no', 'yes'
                elif answ == 'n':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'no', 'yes'
        else:
            if _interactive.upper() in ['YES', 'Y']:
                answ = 'x'
                while answ not in ['y', 'n']:
                    answ = raw_input('\n### do you want to extract again [[y]/n] ? ')
                    if not answ: answ = 'y'
                if answ == 'y':
                    _new, _extract = 'yes', 'yes'
                else:
                    _new, _extract = 'yes', 'no'
            else:
                _new, _extract = 'yes', 'yes'
    if _extract == 'yes':
        floyds.util.delete(imgex)
        if _dispersionline and _new in ['Yes', 'yes', 'YES', 'y', 'Y']:
            question = 'yes'
            while question == 'yes':
                _z1, _z2, goon = floyds.util.display_image(img, 1, '', '', False)
                dist = raw_input(
                    '\n### At which line do you want to extract the spectrum [' + str(dv['line'][_grism]) + '] ? ')
                if not dist: dist = dv['line'][_grism]
                try:
                    dist = int(dist)
                    question = 'no'
                except:
                    print '\n### input not valid, try again:'
        else:
            dist = dv['line'][_grism]
        if _ext_trace in ['yes', 'Yes', 'YES', True]:
            lista = glob.glob('*ex.fits')
            if lista:
                for ii in lista: print ii
                _reference = raw_input('\### which object do you want to use for the trace [' + str(lista[0]) + '] ? ')
                if not _reference: _reference = lista[0]
                _reference = re.sub('_ex', '', _reference)
                _fittrac = 'no'
                _trace = 'no'
            else:
                sys.exit('\n### error: no extracted spectra in the directory')
        else:
            _reference = ''
            _fittrac = 'yes'
            _trace = 'yes'
        if _new == 'no':
            if not os.path.isfile('database/ap' + re.sub('.fits', '', img)):
                floyds.util.repstringinfile('database/ap' + re.sub('.fits', '', imgfast),
                                            'database/ap' + re.sub('.fits', '', img), re.sub('.fits', '', imgfast),
                                            re.sub('.fits', '', img))
            _find = 'no'
            _recenter = 'no'
            _edit = 'no'
            _trace = 'no'
            _fittrac = 'no'
            _mode = 'h'
            _resize = 'no'
            _review = 'no'
            iraf.specred.mode = 'h'
            _interactive = 'no'
        else:
            iraf.specred.mode = 'q'
            _mode = 'q'
            _find = 'yes'
            _recenter = 'yes'
            _edit = 'yes'
            _review = 'yes'
            _resize = dv[_type]['_resize']
        #   extraction for agn, we always reset the parameters 
        if _type=='agn':
            if os.path.isfile('database/ap' + re.sub('.fits', '', img)):
                os.system('rm database/ap' + re.sub('.fits', '', img))

        iraf.specred.apall(img, output=imgex, referen=_reference, trace=_trace, fittrac=_fittrac, find=_find,
                           recenter=_recenter, edit=_edit,
                           nfind=1, extract='yes', backgro='fit', gain=_gain, readnoi=_rdnoise, lsigma=4, usigma=4,
                           format='multispec',
                           b_function='legendre', b_sample=dv[_type]['_b_sample'], clean='yes', pfit='fit1d',
                           b_naver=dv[_type]['_b_naver'],
                           lower=dv[_type]['_lower'], upper=dv[_type]['_upper'], t_niter=dv[_type]['_t_niter'],
                           width=dv[_type]['_width'],
                           radius=dv[_type]['_radius'], line=dist, nsum=dv[_type]['_nsum'], t_step=dv[_type]['_t_step'],
                           t_nsum=dv[_type]['_t_nsum'],
                           t_nlost=dv[_type]['_t_nlost'], t_sample=dv[_type]['_t_sample'], resize=_resize,
                           t_order=dv[_type]['_t_order'],
                           weights=dv[_type]['_weights'], interactive=_interactive, review=_review, mode=_mode)
        floyds.util.repstringinfile('database/ap' + re.sub('.fits', '', img),
                                    'database/ap' + re.sub('.fits', '', imgfast), re.sub('.fits', '', img),
                                    re.sub('.fits', '', imgfast))

        data, hdr = fits.getdata(imgex, 0, header=True)
        xxex = arange(len(data[0][0]))
        aaex = readkey3(hdr, 'CRVAL1') + (xxex) * readkey3(hdr, 'CD1_1')
        floyds.util.updateheader(imgex, 0, {'XMIN': [aaex[0], 'min wavelength [Angstrom]'],
                                            'XMAX': [aaex[-1], 'max wavelength [Angstrom]']})
    else:
        print '\n### skipping new extraction'
    return imgex


################################################

def choseflat(obj, listflat, setup, _JD0, _interactive):
    from floyds.util import sortbyJD
    from floyds.util import readhdr, readkey3, display_image
    from numpy import abs, argsort, array, compress, argmin

    differences = []
    obidflat = []
    flatgood = []
    _OBID = readkey3(readhdr(obj), 'esoid')
    listflat = sortbyJD(listflat)
    for flat in listflat:
        hdrf = readhdr(flat)
        _JDf = readkey3(hdrf, 'JD')
        differences.append(abs(float(_JD0) - float(_JDf)))
        obidflat.append(readkey3(hdrf, 'esoid'))
    if _interactive:
        for flat in listflat:
            _JDf = readkey3(readhdr(flat), 'JD')
            display_image(flat, 1, '', '', False)
            answ = raw_input('### good/bad/stop(enough files, go on) [[g],b,s]')
            if not answ: answ = 'g'
            if answ in ['G', 'g', 'good', 'Good']:
                flatgood.append(flat)
            elif answ in ['stop', 'S', 'Stop', 's']:
                break
    else:
        if obidflat.count(_OBID) >= 3:
            flatgood = list(compress(array(obidflat) == _OBID, listflat))
            print '### Flat field in the same OB !!'
        else:
            print '### ', str(_OBID)
            inds = array(differences).argsort()
            for i in range(0, 3):
                flatgood.append(listflat[inds[i]])
    return flatgood


############################################################
def continumsub(imagefile, _order1, _order2):
    from floyds.util import delete
    from pyraf import iraf

    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['specred.continuum']
    for t in toforget: iraf.unlearn(t)
    delete('tsky.fits')
    output = 'subtracted.fits'
    iraf.specred.continuum(imagefile, output='tsky.fits', type='difference', interact='no', function='legendre',
                           niterat=300, low_rej=3, high_re=2, sample='*', order=_order1, ask='YES')
    iraf.specred.continuum('tsky.fits', output=output, type='difference', interact='no', function='spline1',
                   overrid='yes', niterat=10, low_rej=3, high_re=1, sample='*', order=_order2, ask='YES')
    delete('tsky.fits')
    return output

##########################################
def imreplace_region(img):
    from floyds.util import readhdr, readkey3

    _grism = readkey3(readhdr(img), 'grism')
    from pyraf import iraf

    iraf.imutil(_doprint=0)
    iraf.unlearn('imutil.imreplace')
    if _grism == 'blue':
        iraf.imutil.imreplace(img + '[*,1:200]', value=1, lower='INDEF', upper='INDEF')
        print '### replace pixel 1:200 with 1 (y axes)'
    elif _grism == 'red':
        iraf.imutil.imreplace(img + '[*,1:200]', value=1, lower='INDEF', upper='INDEF')
        print '### replace pixel 1:200 with 1 (y axes)'
    else:
        print '### no replace '


#####################################################
def floydsspecreduction(files, _interactive, _dobias, _doflat, _listflat, _listbias, _listarc, _cosmic, _ext_trace,
                        _dispersionline, liststandard, listatmo, _automaticex, _classify=False, _verbose=False,
                        smooth=1, fringing=1, _typefromuser='obj', fts_contaminated=False):
    import floyds
    from floyds.util import readhdr, readkey3
    import string, re, os, sys, glob
    from numpy import arange, pi, arccos, sin, cos, argmin, sqrt
    from astropy.io import fits
    from pyraf import iraf
    import datetime

    os.environ["PYRAF_BETA_STATUS"] = "1"
    iraf.set(direc=floyds.__path__[0] + '/')
    _extinctdir = 'direc$standard/extinction/'
    header = readhdr(re.sub('\n', '', files[0]))
    _tel = readkey3(header, 'TELID')
    camera = readkey3(header, 'INSTRUME')
    if _tel in ['fts', 'coj']:
        _extinction = 'ssoextinct.dat'
        _observatory = 'sso'
    elif _tel in ['ftn', 'ogg']:
        _extinction = 'maua.dat'
        _observatory = 'cfht'
    else:
        sys.exit('ERROR: observatory not recognised')
    dv = floyds.util.dvex()
    scal = pi / 180.
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['ccdred.flatcombine', 'ccdred.zerocombine', 'ccdproc', 'specred.apall', 'longslit.identify',
                'longslit.reidentify',
                'specred.standard', 'longslit.fitcoords', 'specred.transform', 'specred.response']
    for t in toforget: iraf.unlearn(t)
    iraf.longslit.dispaxi = 2
    iraf.longslit.mode = 'h'
    iraf.identify.fwidth = 7
    iraf.identify.order = 2
    iraf.specred.dispaxi = 2
    iraf.specred.mode = 'h'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.trim = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.overscan = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.biassec = ''
    iraf.ccdproc.ccdtype = ''
    iraf.ccdred.instrument = "/dev/null"
    if _verbose:
        iraf.ccdred.verbose = 'yes'
        iraf.specred.verbose = 'yes'
    else:
        iraf.specred.verbose = 'no'
        iraf.ccdred.verbose = 'no'
    now = datetime.datetime.now()
    datenow = now.strftime('20%y%m%d%H%M')
    MJDtoday = floyds.util.mjdtoday()
    outputlist = []
    hdr0 = floyds.util.readhdr(re.sub('\n', '', files[0]))
    _gain = readkey3(hdr0, 'gain')
    _rdnoise = readkey3(hdr0, 'ron')
    std, rastd, decstd, magstd = floyds.util.readstandard('standard_floyds_mab.txt')
    _naxis2 = hdr0.get('NAXIS2')
    _naxis1 = hdr0.get('NAXIS1')
    if not _naxis1: _naxis1 = 2079
    if not _naxis2:
        if not hdr0.get('HDRVER'):
            _naxis1 = 511
        else:
            _naxis1 = 512

    _overscan = '[2049:' + str(_naxis1) + ',1:' + str(_naxis2) + ']'
    _biassecblu = '[379:2047,325:511]'
    _biassecred = '[1:1801,1:351]'
    lista = {}
    objectlist = {}
    biaslist = {}
    flatlist = {}
    flatlistd = {}
    arclist = {}
    max_length = 22 # max length of FITS header values is 68, and filenames must fit in header
    for img in files:
        hdr0 = readhdr(img)
        if readkey3(hdr0, 'naxis2') >= 500:
            if 'blu' not in lista: lista['blu'] = []
            if 'red' not in lista: lista['red'] = []
            _object0 = floyds.util.readkey3(hdr0, 'object')
            # remove all the following characters from the filename:
            _object0=floyds.util.readkey3(hdr0,'object')
            _object0 = re.sub(':', '', _object0) # colon
            _object0 = re.sub('/', '', _object0) # slash
            _object0 = re.sub('\s', '', _object0) # any whitespace
            _object0 = re.sub('\(', '', _object0) # open parenthesis
            _object0 = re.sub('\[', '', _object0) # open square bracket
            _object0 = re.sub('\)', '', _object0) # close parenthesis
            _object0 = re.sub('\]', '', _object0) # close square bracket
            if len(_object0) > max_length:
                _object0 = _object0[:max_length]
            _date0 = readkey3(hdr0, 'date-night')
            _tel = readkey3(hdr0, 'TELID')
            _type = readkey3(hdr0, 'OBSTYPE')
            if not _type:    _type = readkey3(hdr0, 'imagetyp')
            _slit = readkey3(hdr0, 'slit')
            if _type:
                if _type.lower() in ['expose', 'sky', 'spectrum']:
                    nameoutb = str(_object0) + '_' + _tel + '_' + str(_date0) + '_blue_' + str(_slit) + '_' + str(
                        MJDtoday)
                    nameoutr = str(_object0) + '_' + _tel + '_' + str(_date0) + '_red_' + str(_slit) + '_' + str(
                        MJDtoday)
                elif _type.lower() in ['lamp', 'arc', 'l']:
                    nameoutb = 'arc_' + str(_object0) + '_' + _tel + '_' + str(_date0) + '_blue_' + str(
                        _slit) + '_' + str(MJDtoday)
                    nameoutr = 'arc_' + str(_object0) + '_' + _tel + '_' + str(_date0) + '_red_' + str(
                        _slit) + '_' + str(MJDtoday)
                elif _type.lower() in ['flat', 'f', 'lampflat', 'lamp-flat']:
                    nameoutb = 'flat_' + str(_object0) + '_' + _tel + '_' + str(_date0) + '_blue_' + str(
                        _slit) + '_' + str(MJDtoday)
                    nameoutr = 'flat_' + str(_object0) + '_' + _tel + '_' + str(_date0) + '_red_' + str(
                        _slit) + '_' + str(MJDtoday)
                else:
                    nameoutb = str(_type.lower()) + '_' + str(_object0) + '_' + _tel + '_' + str(
                        _date0) + '_blue_' + str(_slit) + '_' + str(MJDtoday)
                    nameoutr = str(_type.lower()) + '_' + str(_object0) + '_' + _tel + '_' + str(
                        _date0) + '_red_' + str(_slit) + '_' + str(MJDtoday)

                bimg = floyds.util.name_duplicate(img, nameoutb, '')
                rimg = floyds.util.name_duplicate(img, nameoutr, '')
                ####
                floyds.util.delete(bimg)
                floyds.util.delete(rimg)

                iraf.imcopy(img, bimg, verbose='no')
                iraf.imcopy(img, rimg, verbose='no')

                if _type.lower() in ['expose', 'sky', 'spectrum']:
                    #                    iraf.noao.observatory.observa='ftn'
                    #                    iraf.noao.observatory.name='ftn'
                    #                    iraf.noao.observatory.latitude='20:42:29.88'
                    #                    iraf.noao.observatory.longitude='156:15:25.56'
                    #                    iraf.noao.observatory.timezone='10'
                    #                    iraf.noao.observatory.altitude='3050'
                    _observatory = {'ftn': 'cfht', 'ogg': 'cfht', 'fts': 'sso', 'coj': 'sso'}
                    #                    _observatory={'ftn':'obspars','ogg':'obspars','fts':'sso','coj':'sso'}
                    iraf.specred.setjd(bimg, date='DATE-OBS', time='UTSTART', exposure='EXPTIME', ra='ra', dec='dec',
                                       epoch='', observa=_observatory[_tel])
                    iraf.specred.setjd(img, date='DATE-OBS', time='UTSTART', exposure='EXPTIME', ra='ra', dec='dec',
                                       epoch='', observa=_observatory[_tel])

                aaa = iraf.hedit(bimg, 'CCDSEC', delete='yes', update='yes', verify='no', Stdout=1)
                aaa = iraf.hedit(bimg, 'TRIMSEC', delete='yes', update='yes', verify='no', Stdout=1)
                aaa = iraf.hedit(rimg, 'CCDSEC', delete='yes', update='yes', verify='no', Stdout=1)
                aaa = iraf.hedit(rimg, 'TRIMSEC', delete='yes', update='yes', verify='no', Stdout=1)

                iraf.ccdproc(bimg, output='', overscan="yes", trim="yes", zerocor='no', flatcor='no', zero='',
                             ccdtype='',
                             fixpix='no', trimsec=_biassecblu, biassec=_overscan, readaxi='line', Stdout=1)
                iraf.ccdproc(rimg, output='', overscan="yes", trim="yes", zerocor='no', flatcor='no', zero='',
                             ccdtype='',
                             fixpix='no', trimsec=_biassecred, biassec=_overscan, readaxi='line', Stdout=1)
                floyds.util.updateheader(bimg, 0, {'GRISM': ['blu', ' blue order']})
                floyds.util.updateheader(rimg, 0, {'GRISM': ['red', ' blue order']})
                floyds.util.updateheader(bimg, 0, {'arcfile': [img, 'file name in the archive']})
                floyds.util.updateheader(rimg, 0, {'arcfile': [img, 'file name in the archive']})
                lista['blu'].append(bimg)
                lista['red'].append(rimg)
            else:
                'print warning type not defined'
    for arm in lista.keys():
        for img in lista[arm]:
            hdr = readhdr(img)
            _type = readkey3(hdr, 'obstype')
            if not _type:    _type = readkey3(hdr, 'imagetyp')
            _slit = readkey3(hdr, 'slit')
            _grpid = readkey3(hdr, 'grpid')
            if _type.lower() in ['flat', 'f', 'lamp-flat', 'lampflat']:
                if (arm, _slit) not in flatlist:  flatlist[arm, _slit] = {}
                if _grpid not in flatlist[arm, _slit]:
                    flatlist[arm, _slit][_grpid] = [img]
                else:
                    flatlist[arm, _slit][_grpid].append(img)
            elif _type.lower() in ['lamp', 'l', 'arc']:
                if (arm, _slit) not in arclist:  arclist[arm, _slit] = {}
                if _grpid not in arclist[arm, _slit]:
                    arclist[arm, _slit][_grpid] = [img]
                else:
                    arclist[arm, _slit][_grpid].append(img)
            elif _type in ['bias', 'b']:
                if arm not in biaslist: biaslist[arm] = []
                biaslist[arm].append(img)
            elif _type.lower() in ['sky', 's', 'spectrum', 'expose']:
                try:
                    _ra = float(readkey3(hdr, 'RA'))
                    _dec = float(readkey3(hdr, 'DEC'))
                except:
                    ra00 = string.split(readkey3(hdr, 'RA'), ':')
                    ra0, ra1, ra2 = float(ra00[0]), float(ra00[1]), float(ra00[2])
                    _ra = ((ra2 / 60. + ra1) / 60. + ra0) * 15.
                    dec00 = string.split(readkey3(hdr, 'DEC'), ':')
                    dec0, dec1, dec2 = float(dec00[0]), float(dec00[1]), float(dec00[2])
                    if '-' in str(dec0):
                        _dec = (-1) * ((dec2 / 60. + dec1) / 60. + ((-1) * dec0))
                    else:
                        _dec = (dec2 / 60. + dec1) / 60. + dec0
                dd = arccos(sin(_dec * scal) * sin(decstd * scal) + cos(_dec * scal) * cos(decstd * scal) * cos(
                    (_ra - rastd) * scal)) * ((180 / pi) * 3600)
                if _verbose:
                    print _ra, _dec
                    print std[argmin(dd)], min(dd)

                if min(dd) < 1200:
                    _typeobj = 'std'
                else:
                    if not _typefromuser:
                        _typeobj = 'obj'
                    else:
                        _typeobj = _typefromuser

                if min(dd) < 1200:
                    floyds.util.updateheader(img, 0, {'stdname': [std[argmin(dd)], '']})
                    floyds.util.updateheader(img, 0, {'magstd': [float(magstd[argmin(dd)]), '']})
                if _typeobj not in objectlist:      objectlist[_typeobj] = {}

                if (arm, _slit) not in objectlist[_typeobj]:
                    objectlist[_typeobj][arm, _slit] = [img]
                else:
                    objectlist[_typeobj][arm, _slit].append(img)
    if _verbose:
        print 'object'
        print objectlist
        print 'flat'
        print flatlist
        print 'bias'
        print biaslist
        print 'arc'
        print arclist

    if liststandard and 'std' in objectlist.keys():
        print 'external standard, raw standard not used'
        del objectlist['std']

    sens = []
    outputfile = []
    atmo = {}
    wavecalib = {}
    for tpe in objectlist:
        for setup in objectlist[tpe]:
#            if setup not in sens:   
#                sens[setup] = []
            print '\n### setup= ', setup, '\n### objects= ', objectlist[tpe][setup], '\n'
            for img in objectlist[tpe][setup]:
                print '\n\n### next object= ', img, ' ', floyds.util.readkey3(readhdr(img), 'object'), '\n'
                hdr = floyds.util.readhdr(img)
                _gain = floyds.readkey3(hdr, 'gain')
                _rdnoise = floyds.readkey3(hdr, 'ron')
                _grism = floyds.readkey3(hdr, 'grism')
                _grpid = floyds.readkey3(hdr, 'grpid')
                #####################      flat   ###############
                if _listflat:
                    flatgood = _listflat  # flat list from reducer
                elif setup in flatlist:
                    if _grpid in flatlist[setup]:
                        print '\n###FLAT WITH SAME GRPID'
                        flatgood = flatlist[setup][_grpid]  # flat in the  raw data
                    else:
                        flatgood = []
                        for _grpid0 in flatlist[setup].keys():
                            for ii in flatlist[setup][_grpid0]:
                                flatgood.append(ii)
                else:
                    flatgood = []
                if len(flatgood) != 0:
                    if len(flatgood) > 1:
                        f = open('_oflatlist', 'w')
                        for fimg in flatgood:
                            f.write(fimg + '\n')
                        f.close()
                        floyds.util.delete('flat' + img)
                        iraf.ccdred.flatcombine('"@_oflatlist"', output='flat' + img, combine='average', reject='none',
                                                ccdtype=' ', rdnoise=_rdnoise, gain=_gain, process='no', Stdout=1)
                        flatfile = 'flat' + img
                        floyds.util.delete('_oflatlist')
                    elif len(flatgood) == 1:
                        os.system('cp ' + flatgood[0] + ' flat' + img)
                        flatfile = 'flat' + img
                else:
                    flatfile = ''
                ##########################   find arcfile            #######################
                arcfile = ''
                if _listarc:
                    arcfile = [floyds.util.searcharc(img, _listarc)[0]][0]  # take arc from list
                if not arcfile and setup in arclist.keys():
                    if _grpid in arclist[setup]:
                        print '\n###ARC WITH SAME GRPID'
                        arcfile = arclist[setup][_grpid]  # arc in the  raw data
                    else:
                        arcfile = []
                        for _grpid0 in arclist[setup].keys():
                            for ii in arclist[setup][_grpid0]:
                                arcfile.append(ii)
                if arcfile:
                    if len(arcfile) > 1:  # more than one arc available
                        #                     _arcclose=floyds.util.searcharc(imgex,arcfile)[0]   # take the closest in time 
                        _arcclose = floyds.util.sortbyJD(arcfile)[-1]  #  take the last arc of the sequence
                        if _interactive.upper() in ['YES', 'Y']:
                            for ii in floyds.util.sortbyJD(arcfile):
                                print '\n### ', ii
                            arcfile = raw_input(
                                '\n### more than one arcfile available, which one to use [' + str(_arcclose) + '] ? ')
                            if not arcfile: arcfile = _arcclose
                        else:
                            arcfile = _arcclose
                    else:
                        arcfile = arcfile[0]
                else:
                    print '\n### Warning: no arc found'

                ###################################################################
                if setup[0] == 'red':
                    fcfile = floyds.__path__[0] + '/standard/ident/' + camera + '/fcrectify_' + _tel + '_red'
                    fcfile1 = floyds.__path__[0] + '/standard/ident/' + camera + '/fcrectify1_' + _tel + '_red'
                    fcfile_untilt = floyds.__path__[0] + '/standard/ident/' + camera + '/fcuntilt_' + _tel + '_red'
                    print fcfile
                else:
                    fcfile = floyds.__path__[0] + '/standard/ident/' + camera + '/fcrectify_' + _tel + '_blue'
                    fcfile1 = floyds.__path__[0] + '/standard/ident/' + camera + '/fcrectify1_' + _tel + '_blue'
                    fcfile_untilt = floyds.__path__[0] + '/standard/ident/' + camera + '/fcuntilt_' + _tel + '_blue'

                    print fcfile

                if not img:  # or not arcfile:
                    print '\n### no calibration for this setup available'
                else:
                    print img, arcfile, flatfile, fcfile, fcfile1, _cosmic
                    img, arcfile, flatfile = floyds.floydsspecdef.rectifyspectrum(img, arcfile, flatfile, fcfile,
                                                                                  fcfile1, fcfile_untilt, 'no', _cosmic)

                ###############     flat correction  method 1 or 3  ###############
                if img[:2] != 'tt':
                    img = 'tt' + img
                print img, arcfile, flatfile, setup
                imgn = 'n' + img
                if fringing in [1, 3]:
                    if flatfile and setup[0] == 'red':
                        print flatfile, img
                        if fringing == 1:
                            print 'fringing 1 classical'
                            imgn = floyds.fringing_classicmethod(flatfile, img, 'no', '*', 15, setup[0])
                        elif fringing == 3:
                            print 'fringing 3 apsum'
                            imgn = floyds.fringing_classicmethod2(flatfile, img, 'no', '*', 15, setup[0])
                    else:
                        floyds.util.delete(imgn)
                        iraf.imcopy(img, imgn, verbose='yes')
                else:
                    floyds.util.delete(imgn)
                    iraf.imcopy(img, imgn, verbose='no')
                ################    extraction         ####################################
                imgex = floyds.floydsspecdef.extractspectrum(imgn, dv, _ext_trace, _dispersionline, _interactive, tpe,
                                                             automaticex=_automaticex)
                ####################################################
                if setup[0] == 'blu' and smooth != 1:
                    print 'smooth blu part of ' + str(smooth) + ' pixels'
                    from numpy import float32

                    x1, y1 = floyds.util.readspectrum(imgex)
                    y2 = rebin(y1, smooth)
                    data1, hdr = fits.getdata(imgex, 0, header=True)
                    data1[0][0] = y2
                    floyds.util.delete(imgex)
                    fits.writeto(imgex, float32(data1), hdr)
                    floyds.util.updateheader(imgex, 0, {'REBIN': [smooth, ' spectrum rebinned']})
                ####################################################################
                if not os.path.isfile(imgex): sys.exit('### error, extraction not computed')
                if imgex not in outputfile: outputfile.append(imgex)
                ###############     flat correction methid 2  ###############
                if fringing == 2:
                    if setup[0] == 'red' and flatfile:
                        floyds.util.delete('flat' + imgex)
                        iraf.specred.apsum(flatfile, output='flat' + imgex, referen=imgn, interac='no', find='no',
                                           recente='no', resize='no', edit='no', trace='no', fittrac='no',
                                           extract='yes', extras='no', review='no', backgro='none')
                        fringingmask = normflat('flat' + imgex)
                        print '\n### fringing correction'
                        imgex, so2, shift = correctfringing_auto(imgex, fringingmask)  #  automatic correction
                        #imgex=correctfringing(imgex,fringingmask)      #  manual correction
                    #############################################
                if arcfile:
                    floyds.util.delete('arc_' + imgex)
                    iraf.specred.apsum(arcfile, output='arc_' + imgex, referen=imgn, interac='no', find='no',
                                       recente='no', resize='no', edit='no', trace='no', fittrac='no',
                                       extract='yes', extras='no', review='no', backgro='none')
                    arcfile = 'arc_' + imgex
                if not arcfile:  #   search arc in the archive
                    arcfile = floyds.util.searcharc(imgex, '')[0]
                    if arcfile:
                        os.system('cp ' + arcfile + ' arc_' + imgex)
                        arcfile = ' arc_' + imgex
                else:
                    print 'arcfile = ' + arcfile
                ############################################
                if not arcfile:
                    print '\n### warning no arcfile \n exit '
                    imgl = ''
                else:
                    imgl = imgex.replace('_ex.fits', '_l.fits')
                    oldfiletoday = os.path.isfile(imgl)
                    if not oldfiletoday:
                        oldfiles = [f.replace('_l', '_ex') for f in glob.glob(re.sub(str(floyds.util.mjdtoday()), '*', imgl))]
                    if _automaticex:
                        default = 'o'
                    else:
                        default = 'n'
                    ans3 = ''
                    if _interactive.lower() in ['yes', 'y']:
                        if oldfiletoday:
                            while ans3 not in ['n', 'o', 's']:
                                ans3 = raw_input('\n### [n]ew wavelength calibration, redo calibration with [o]ld parameters, '
                                                 '[s]kip calibration (use previous file) [' + default + '] ')
                                if not ans3:
                                    ans3 = default
                        elif oldfiles:
                            while ans3 not in ['n', 'o']:
                                ans3 = raw_input('\n### [n]ew wavelength calibration, redo calibration with [o]ld parameters, [' + default + '] ')
                                if not ans3:
                                    ans3 = default
                        else:
                            ans3 = 'n'
                    else:
                        ans3 = default

#                    if _interactive.upper() in ['YES', 'Y']:
#                        if os.path.isfile(imgl):
#                            answ = raw_input('\n### wavelength calibrated file already there, '
#                                             'do you want to calibrate again [[y]/n] ? ')
#                            if not answ: answ = 'y'
#                        else:
#                            answ = 'y'
#                    else:
#                        answ = 'y'

                    if ans3 == 'n':
                        arcref = floyds.util.searcharc(imgex, '')[0]
                        wlident = True
                        wlcal = True
                    elif ans3 == 'o' and not oldfiletoday:
                        arcref = 'arc_' + max(oldfiles)
                        wlident = True
                        wlcal = True
                    elif ans3 == 'o':
                        wlident = False
                        wlcal = True
                    else:
                        wlident = False
                        wlcal = False

#                    if answ in ['y', 'Y', 'YES', 'yes', 'Yes']:
                    if wlident:
                        print imgex, setup
                        if setup[0] == 'blu':
                            _order = 3
                        else:
                            _order = 6
                        if setup[1] in ['6.0']:
                            _fwidth = 30
                            _specres = 40
                            _cradius = 30
                        else:
                            _fwidth = 7
                            _specres = 25
                            _cradius = 10

#                        arcref = floyds.util.searcharc(imgex, '')[0]
                        os.system('cp ' + floyds.__path__[0] + '/standard/ident/FLOYDS_lines.txt .')

                        if not arcref:
                            identific = iraf.specred.identify(images=arcfile, section='middle line',
                                                              coordli='FLOYDS_lines.txt', nsum=10, fwidth=_fwidth,
                                                              cradius=_cradius, function='legendre', order=_order,
                                                              mode='h', Stdout=1)
                        else:
                            if arcref[0] == '/':
                                os.system('cp ' + arcref + ' .')
                                arcref = string.split(arcref, '/')[-1]
                            Path('database/').mkdir(parents=True, exist_ok=True)
                            if os.path.isfile(floyds.util.searcharc(imgex, '')[1] + '/database/id' + re.sub('.fits', '',
                                                                                                            arcref)):
                                os.system('cp ' + floyds.util.searcharc(imgex, '')[1] +
                                          '/database/id' + re.sub('.fits', '', arcref) + ' database/')
                            identific = iraf.specred.reidentify(referenc=arcref, images=arcfile, interac=_interactive,
                                                                section='middle line', coordli='FLOYDS_lines.txt',
                                                                overrid='yes', cradius=_cradius, step=0, newaps='no',
                                                                nsum=5, nlost=2, mode='h', verbose='yes', Stdout=1)
                            if _interactive.upper() in ['YES', 'Y']:
                                answ = raw_input('### do you like the identification [[y]/n]')
                                if not answ: answ = 'y'
                            else:
                                import time

                                time.sleep(1)
                                answ = 'y'
                            if answ in ['n', 'N', 'no', 'NO', 'No']:
                                yy1 = fits.open(arcref)[0].data
                                xx1 = arange(len(yy1))
                                yy2 = fits.open(arcfile)[0].data
                                xx2 = arange(len(yy2))
                                _shift = floyds.floydsspecdef.checkwavelength_arc(xx1, yy1, xx2, yy2, '', '',
                                                                                  _interactive)  #*(-1)
                                print '\n### computed shift= ' + str(_shift)
                                identific = iraf.specred.reidentify(referenc=arcref, images=arcfile,
                                                                    interac=_interactive, section='middle line',
                                                                    shift=_shift, coordli='FLOYDS_lines.txt',
                                                                    overrid='yes', cradius=_cradius, step=0,
                                                                    newaps='no', nsum=5, nlost=2, mode='h',
                                                                    verbose='yes', Stdout=1)
                                answ = raw_input('### is it ok now [[y]/n]')
                                if not answ: answ = 'y'
                                if answ in ['n', 'N', 'no', 'NO', 'No']:
                                    if setup[0] == 'blu':
                                        _order = 3
                                    else:
                                        _order = 6
                                    identific = iraf.specred.identify(images=arcfile, section='middle line',
                                                                      coordli='FLOYDS_lines.txt',
                                                                      nsum=10, fwidth=_fwidth, cradius=_cradius,
                                                                      function='legendre', order=_order, mode='h',
                                                                      Stdout=1)

                        floyds.util.delete('FLOYDS_lines.txt')
                        imgl = re.sub('_ex.fits', '_l.fits', imgex)
                        print arcfile
                        try:
                            specred = floyds.util.spectraresolution3(arcfile, _specres)
                        except:
                            specred = 0
                        print identific
                        if identific:
                            try:
                                _rms = float(identific[-1].split()[-1])
                                _num = float(identific[-1].split()[4].split('/')[0])
                            except:
                                _rms = 9999
                                _num = 9999
                            hedvec = {'LAMRMS_' + setup[0][0]: [_rms * .1, 'residual RMS [nm]'],
                                      'LAMNLIN' + setup[0][0]: [_num, 'Nb of arc lines used in the fit of '
                                                                      'the wavel. solution'],
                                      'SPE_ER_' + setup[0][0]: [(_rms * .1) / sqrt(float(_num)),
                                                                'statistical uncertainty'],
                                      'REFSPEC1': [re.sub('.fits', '', arcfile), ' reference arc'],
                                      'arc' + setup[0]: [re.sub('.fits', '', arcfile), ' reference arc']}
                        else:
                            hedvec = {'REFSPEC1': [re.sub('.fits', '', arcfile), ' reference arc'],
                                      'arc' + setup[0]: [re.sub('.fits', '', arcfile), ' reference arc']}

                        if specred:         
                            hedvec['SPERES_' + setup[0][0]] = [specred, 'Spectral resolving power']

                        floyds.util.updateheader(imgex, 0, hedvec)

                    if wlcal:
                        floyds.util.delete(imgl)

                        iraf.specred.dispcor(imgex, output=imgl, flux='yes')

                        if imgl not in outputfile: 
                            outputfile.append(imgl)
                        if tpe == 'std' or floyds.util.readkey3(floyds.util.readhdr(imgex), 'exptime') < 300:
                            if setup[0] == 'red':
                                print '\n### check standard wave calib'
                                print imgl
                                shiftred = floyds.floydsspecdef.checkwavestd(imgl, _interactive, 2)
                                zro = fits.open(imgl)[0].header['CRVAL1']
                                if _interactive.upper() in ['YES', 'Y']:
                                    print 'Figure out the shift of H-alpha (should be near 6563).'
                                    iraf.onedspec.splot(imgl + '[*,1,1]')
                                    num = raw_input('By how much do you want to shift the wavelength calibration? [' + str(shiftred) + '] ')
                                    if not num:
                                        num = shiftred
                                    floyds.util.updateheader(imgl, 0, {'CRVAL1': [zro + float(num), '']})
                                    floyds.util.updateheader(imgl, 0, {'shift' + setup[0][0]: [float(num), '']})
                                else:
                                    print 'shifting wavelength calibration by', shiftred
                                    floyds.util.updateheader(imgl, 0, {'CRVAL1': [zro + float(shiftred), '']})
                                    floyds.util.updateheader(imgl, 0, {'shift' + setup[0][0]: [float(shiftred), '']})
                            else:
                                print re.sub('blue','red',imgl)
                                if os.path.isfile(re.sub('blue','red',imgl)):
                                    hdrred = floyds.util.readhdr(re.sub('blue','red',imgl))
                                    shiftred = floyds.readkey3(hdrred,'SHIFTR')
                                    if shiftred:
                                        zro = fits.open(imgl)[0].header['CRVAL1']
                                        if _interactive.upper() in ['YES', 'Y']:
                                            print 'Figure out the shift of H-beta (should be near 4861).'
                                            iraf.onedspec.splot(imgl + '[*,1,1]')
                                            num = raw_input('By how much do you want to shift the wavelength calibration? [' + str(shiftred) + '] ')
                                            if not num:
                                                num = shiftred
                                            floyds.util.updateheader(imgl, 0, {'CRVAL1': [zro + float(num), '']})
                                            floyds.util.updateheader(imgl, 0, {'shift' + setup[0][0]: [float(num), '']})
                                        else:
                                            print 'shifting wavelength calibration by', shiftred
                                            floyds.util.updateheader(imgl, 0, {'CRVAL1': [zro + float(shiftred), '']})
                                            floyds.util.updateheader(imgl, 0, {'shift' + setup[0][0]: [float(shiftred), '']})
                                    else:
                                        print shiftred
                                        print 'no shift'
                                else:
                                    print '\n### warning check in wavelength not possible for short exposure in the blu range'
                        else:
                            print '\n### check object wave calib'
                            _skyfile = floyds.__path__[0] + '/standard/ident/sky_' + setup[0] + '.fits'
                            floyds.floydsspecdef.checkwavelength_obj(imgl, _skyfile, _interactive)
                if imgl:
                    if tpe not in wavecalib: 
                        wavecalib[tpe] = {}
                    if setup not in wavecalib[tpe]:
                        wavecalib[tpe][setup] = [imgl]
                    else:
                        wavecalib[tpe][setup].append(imgl)

                ###################################
    print wavecalib
    #sens = {}
    atmo = {}
    if 'std' in wavecalib.keys():
        for setup in wavecalib['std']:
#            if setup not in sens:   
#                sens[setup] = []
            print '\n### setup= ', setup, '\n### objects= ', wavecalib['std'][setup], '\n'
            for imgl in wavecalib['std'][setup]:
                ######################################################
                hdrs = readhdr(imgl)
                _tel = readkey3(hdrs, 'TELID')
                datenight = readkey3(hdrs, 'date-night')
                if _tel not in ['fts', 'ftn']:  _tel = readkey3(hdrs, 'SITEID')
                try:
                    _outputsens2 = 'sens_' + _tel + '_' + datenight + '_' + str(readkey3(hdrs, 'grism')) + \
                                   '_' + re.sub('.dat', '', readkey3(hdrs, 'stdname')) + '_' + str(MJDtoday)
                except:
                    sys.exit('Error: missing header -stdname- in standard ' + str(standardfile) + '  ')

                print '\n### compute sensitivity function and atmofile'
                if fts_contaminated is not None: fts_contam = fts_contaminated         # if an option is given on the command line, use that
                else: fts_contam = (datenight > '20140915' and datenight < '20150616') # by default, use the period that FTS FLOYDS had glycol on the mirrors
                if setup[0] == 'red':
                    atmofile = floyds.floydsspecdef.telluric_atmo(imgl)
                    print atmofile
                    stdusedclean = re.sub('_l.fits', '_clean.fits', imgl)
                    floyds.util.delete(stdusedclean)
                    _function = 'spline3'
                    iraf.specred.sarith(input1=imgl, op='/', input2=atmofile, output=stdusedclean, format='multispec')
                    if _tel in ['fts','coj'] and fts_contam:
                        _outputsens2 = floyds.floydsspecdef.sensfunction(stdusedclean, _outputsens2, _function, 8, _interactive, '4600:6730,6720:10000', fts_contam)
                    else:
                        _outputsens2 = floyds.floydsspecdef.sensfunction(stdusedclean, _outputsens2, _function, 8, _interactive)

                    if setup not in atmo:
                        atmo[setup] = [atmofile]
                    else:
                        atmo[setup].append(atmofile)
                else: # blue arm
                    _function = 'spline3'
                    if _tel in ['fts','coj'] and fts_contam:
                        _outputsens2 = floyds.floydsspecdef.sensfunction(imgl, _outputsens2, _function, 8, _interactive, '3200:4700,4600:5900', fts_contam)
                    else:
                        _outputsens2 = floyds.floydsspecdef.sensfunction(imgl, _outputsens2, _function, 12, _interactive, '3400:4700')  #,3600:4300')

                if _outputsens2 not in sens:
                    print _outputsens2
                    sens.append(_outputsens2)
#                    if _verbose:
                    #    calibrate the standard using the sensitivity just obtained  
#                    _airmass = readkey3(hdrs, 'airmass')
#                    _exptime = readkey3(hdrs, 'exptime')
#                    imgf = re.sub('_l.fits', '_f.fits', imgl)
#                    floyds.util.delete(imgf)
#                    qqq = iraf.specred.calibrate(input=imgl, output=imgf, sensiti=_outputsens2, extinct='yes', flux='yes',
#                                                 extinction=_extinctdir + _extinction, observatory=_observatory,
#                                                 airmass=_airmass, ignorea='yes', exptime=_exptime, fnu='no')


    if _verbose:
        print wavecalib
        print sens
        print atmo

    for tpe in ['std', 'obj', 'agn']:
        if tpe in wavecalib.keys():
            for setup in wavecalib[tpe].keys():
                for img in wavecalib[tpe][setup]:
                    hdr = floyds.util.readhdr(img)
                    _sens = ''
                    if liststandard:  
                        _sens = floyds.util.searchsens(img, liststandard)[0]  # search in the list from reducer
                    if not _sens:
                        try:
#                            _sens = floyds.util.searchsens(img, sens[setup])[0]  # search in the reduced data
                            _sens = floyds.util.searchsens(img, sens)[0]  # search in the reduced data
                        except:
                            print setup
                            print 'no standard'
                            _sens = floyds.util.searchsens(img, '')[0]  # search in the archive

                    _atmo = ''
                    if listatmo:  _atmo = floyds.util.searchatmo(img, listatmo)[
                        0]  # search atmo in the list from reducer 
                    if not _atmo:
                        try:
                            _atmo = floyds.util.searchatmo(img, atmo[setup])[0]  # search in the reduced data
                        except:
                            _atmo = ''
                    if not _atmo:
                        try:
                            _atmo = floyds.util.searchatmo(img, '')[0]  # search in the archive
                        except:
                            _atmo = ''

                    imgf = ''
                    if _sens:
                        _airmass = readkey3(hdr, 'airmass')
                        _exptime = readkey3(hdr, 'exptime')
                        imgf = re.sub('_l.fits', '_f.fits', img)
                        floyds.util.delete(imgf)
                        qqq = iraf.specred.calibrate(input=img, output=imgf, sensiti=_sens, extinct='yes', flux='yes',
                                                     extinction=_extinctdir + _extinction, observatory=_observatory,
                                                     airmass=_airmass, ignorea='yes', exptime=_exptime, fnu='no')
                        floyds.util.updateheader(imgf, 0, {
                        'sensfun' + setup[0][0]: [string.split(_sens, '/')[-1], 'sensitivity curve']})
                        hdr = floyds.readhdr(imgf)
                        fileident = str(floyds.readkey3(hdr, 'JD'))[:9] + ' ' + floyds.readkey3(hdr,
                                                                                                'object') + ' ' + str(
                            floyds.readkey3(hdr, 'date-obs'))[:4] + '-' + \
                                    str(floyds.readkey3(hdr, 'date-obs'))[4:6] + '-' + str(
                            floyds.readkey3(hdr, 'date-obs'))[6:8] + ' ' + str(
                            floyds.readkey3(hdr, 'grism')) + ' floyds.' + floyds.__version__
                        floyds.util.updateheader(imgf, 0, {'IDENT': [fileident, 'file identification']})
                        if imgf not in outputfile: outputfile.append(imgf)
                        if _sens not in outputfile: outputfile.append(_sens)
                        if _atmo:
                            imge = re.sub('_f.fits', '_e.fits', imgf)
                            floyds.util.delete(imge)
                            iraf.specred.sarith(input1=imgf, op='/', input2=_atmo, output=imge, w1='INDEF', w2='INDEF',
                                                format='multispec')
                            try:
                                iraf.imutil.imcopy(input=imgf + '[*,1,2]', output=imge + '[*,1,2]', verbose='no')
                            except:
                                pass
                            try:
                                iraf.imutil.imcopy(input=imgf + '[*,1,3]', output=imge + '[*,1,3]', verbose='no')
                            except:
                                pass
                            try:
                                iraf.imutil.imcopy(input=imgf + '[*,1,4]', output=imge + '[*,1,4]', verbose='no')
                            except:
                                pass
                            if imge not in outputfile: 
                                outputfile.append(imge)
                            floyds.util.updateheader(imge, 0,
                                                     {'ATMO' + setup[0][0]: [string.split(_atmo, '/')[-1], '']})
                            imgin = imge
                            if _atmo not in outputfile: 
                                outputfile.append(_atmo)
                        else:
                            imgin = imgf
                        imgasci = re.sub('.fits', '.asci', imgin)
                        floyds.util.delete(imgasci)
                        iraf.onedspec(_doprint=0)
                        iraf.onedspec.wspectext(imgin + '[*,1,1]', imgasci, header='no')
                        if imgasci not in outputfile: 
                            outputfile.append(imgasci)
    coppie = {}
    for obj in outputfile:
        if obj[-4:] == 'fits':
            hdr = floyds.util.readhdr(obj)
            if '_e.fits' in obj:
                MJD = floyds.util.readkey3(hdr, 'MJD')
                if MJD not in coppie:
                    coppie[MJD] = [obj]
                else:
                    coppie[MJD].append(obj)
            elif '_f.fits' in obj:
                if re.sub('_f.fits', '_e.fits', obj) not in outputfile:
                    MJD = floyds.util.readkey3(hdr, 'MJD')
                    if MJD not in coppie:
                        coppie[MJD] = [obj]
                    else:
                        coppie[MJD].append(obj)

    print coppie
    _i = _interactive.lower() in ['y', 'yes']
    for mjd in coppie.keys():
        lista = coppie[mjd]
        _output = re.sub('_red_', '_merge_', lista[0])
        _output = re.sub('_blue_', '_merge_', _output)
        _output = combspec2(lista[0], lista[1], _output, scale=True, num=None)
        print _output
        if '_e.fits' in _output:
            _output_f = re.sub('_e.fits', '_f.fits', _output)
            if '_e.fits' in lista[0]: lista[0] = re.sub('_e.fits', '_f.fits', lista[0])
            if '_e.fits' in lista[1]: lista[1] = re.sub('_e.fits', '_f.fits', lista[1])
            _output_f = combspec2(lista[0], lista[1], _output_f, scale=True, num=None)
            if _classify: aa, bb, cc = floyds.util.classifyfast(_output, program='snid')
        if _i: _trim = raw_input('Do you want to trim the edges of the spectrum? [[y]/n] ')
        if not _i or _trim != 'n':  # if not interactive, spectrum is trimmed at default boundaries (3200,10000)
            trimmed = 'trim_' + _output
            floyds.util.delete(trimmed)
            while True:
                if _i:
                    print "Find the lower and upper wavelength limits. Then press 'q' to continue."
                    iraf.onedspec.splot(_output+'[*,1,1]')
                    low = raw_input_num('Lower limit (\xc3\x85)',3200)
                    up  = raw_input_num('Upper limit (\xc3\x85)',10000)
                else:
                    low = 3200
                    up = 10000
                iraf.scopy(_output,trimmed,w1=low,w2=up,rebin='no')
########################################################
#                    iraf.onedspec.splot(_output + '[*,1,1]')
#                    low = raw_input('Lower limit (in angstroms) [3200]: ')
#                    up = raw_input('Upper limit (in angstroms) [10000]: ')
#                if not _i or not low: low = 3200
#                if not _i or not up:  up = 10000
#                iraf.scopy(_output, trimmed, w1=low, w2=up, rebin='no')
#########################################################
                if _i:
                    print "See if this looks better. Then press 'q' to continue."
                    iraf.onedspec.splot(trimmed + '[*,1,1]')
                    okay = raw_input('Is this okay? [[y]/n] ')
                if not _i or okay != 'n':
                    break
                else:
                    floyds.util.delete(trimmed)
                    again = raw_input('Do you want to try again? [[y]/n] ')
                    if again == 'n':
                        trimmed = _output  # so it says the right output file below
                        break
            print 'Output saved as', trimmed
        else:
            print 'Output saved as', _output  # only not trimmed if you specifically say no
    return outputfile


def raw_input_num(base_prompt,default=None):
    if default is None:
        prompt = base_prompt+': '
    else:
        prompt = base_prompt+' ['+str(default)+']: '
    while True:
        num = raw_input(prompt)
        try: # if there was a response, cast to float and continue
            if num: num = float(num)
            else:   num = default
            break
        except ValueError: # if response can't be cast to float, reprompt
            print 'Not a valid number.'
    return num

##############################################################################

def cutstd(stdfile, start=1, end=1e10, out=True):
    import os, string, re, sys
    import floyds

    f = open(stdfile, 'r')
    ss = f.readlines()
    f.close()
    floyds.util.delete(stdfile)
    f = open(stdfile, 'w')
    if out:  # use only the values inside start and end
        f.write(re.sub(string.split(ss[0])[-3], str(start), re.sub(string.split(ss[0])[-2], str(end), ss[0])))
    else:
        f.write(ss[0])
    for line in ss[1:]:
        if out:  # use only the values inside start and end
            if float(string.split(line)[0]) >= float(start) and float(string.split(line)[0]) <= float(end):
                f.write(line)
        else:  # don't use the value inside start and end
            if float(string.split(line)[0]) <= float(start) or float(string.split(line)[0]) >= float(end):
                f.write(line)
    f.close()


######## currently using combspec2 instead of combspec ########
#def combspec(_img0, _img1, _output, scale=True, num=None):
#    import numpy as np
#    #    from numpy import compress, array, trapz, argsort
#    #    from numpy import interp as ninterp
#    import re, string, os
#    import floyds
#    from astropy.io import fits
#    _x0, _y0 = floyds.util.readspectrum(_img0)
#    _x1, _y1 = floyds.util.readspectrum(_img1)

#    hdr0 = fits.getheader(_img0)
#    hdr1 = fits.getheader(_img1)
#    if 'NAXIS3' in hdr0 and 'NAXIS3' in hdr1:
#        if hdr0['NAXIS3'] == hdr1['NAXIS3']:
#            dimension = hdr1['NAXIS3']
#        else:
#            dimension = 1
#    else:
#        dimension = 1

#    if min(_x0) < min(_x1):
#        x0, y0 = _x0, _y0
#        x1, y1 = _x1, _y1
#        img0, img1 = _img0, _img1
#    else:
#        x0, y0 = _x1, _y1
#        x1, y1 = _x0, _y0
#        img0, img1 = _img1, _img0

#    limmin = max(min(x0), min(x1))
#    limup = min(max(x0), max(x1))
#    x01 = np.compress((np.array(x0) > limmin) & (np.array(x0) < limup), x0)
#    y01 = np.compress((np.array(x0) > limmin) & (np.array(x0) < limup), y0)
#    x11 = x01
#    if not num:   
#        num = int(len(x01) / 7)
#    y11 = np.interp(x01, x1, y1)
#    if scale:
#        integral0 = np.trapz(y01, x01)
#        integral1 = np.trapz(y11, x11)
#        if integral0 < integral1:
#            A0, A1 = integral1 / integral0, 1
#        else:
#            A0, A1 = 1, integral0 / integral1
#    else:
#        A0, A1 = 1, 1
#    a = abs(y11 * A1 - y01 * A0)
#    b = (y11 * A1 - y01 * A0)
#    a1 = a[0:num]
#    c1 = x01[0:num]
#    a2 = a[-num:]
#    c2 = x01[-num:]
#    from pyraf import iraf

#    floyds.util.delete('s1.fits,s2.fits,s11.fits,s22.fits')
#    floyds.util.delete(_output)

#    if dimension == 1:
#        iraf.scopy(img0, 's1.fits', w1='INDEF', w2=c2[np.argsort(a2)[0]], rebin='no')
#        iraf.scopy(img1, 's2.fits', w1=c1[np.argsort(a1)[0]], w2='INDEF', rebin='no')
#        iraf.sarith(input1='s1.fits', op='*', input2=A0, output='s11.fits')
#        iraf.sarith(input1='s2.fits', op='*', input2=A1, output='s22.fits')
#        iraf.specred.scombine(input='s11.fits,s22.fits', w1='INDEF', w2='INDEF', output=_output)
#    else:
#        print 'more than one dimension'
#        datavec = {}
#        hdrvec = {}
#        for ii in range(dimension):
#            outputn = re.sub('.fits', '', _output) + '_' + str(ii + 1) + '.fits'
#            floyds.util.delete(outputn)
#            floyds.util.delete('s1.fits,s2.fits,s11.fits,s22.fits')
#            iraf.scopy(img0 + '[*,1,' + str(ii + 1) + ']', 's1.fits', w1='INDEF', w2=c2[np.argsort(a2)[0]], rebin='no')
#            iraf.scopy(img1 + '[*,1,' + str(ii + 1) + ']', 's2.fits', w1=c1[np.argsort(a1)[0]], w2='INDEF', rebin='no')
#            iraf.sarith(input1='s1.fits', op='*', input2=A0, output='s11.fits')
#            iraf.sarith(input1='s2.fits', op='*', input2=A1, output='s22.fits')
#            iraf.specred.scombine(input='s11.fits,s22.fits', w1='INDEF', w2='INDEF', output=outputn)
#            datavec[ii], hdrvec[ii] = fits.getdata(outputn, 0, header=True)

#        datat = np.array([[datavec[0]], [datavec[1]], [datavec[2]], [datavec[3]]])
#        floyds.util.delete(_output)
#        try:
#            hdr0.update('NAXIS1', hdrvec[1]['NAXIS1'], hdrvec[1].comments['NAXIS1'])
#            hdr0.update('CRVAL1', hdrvec[1]['CRVAL1'], hdrvec[1].comments['CRVAL1'])
#            hdr0.update('CD1_1', hdrvec[1]['CD1_1'], hdrvec[1].comments['CD1_1'])
#            hdr0.update('CRPIX1', hdrvec[1]['CRPIX1'], hdrvec[1].comments['CRPIX1'])
#        except:
#            hdr0.update('NAXIS1', hdrvec[1]['NAXIS1'], 'Width of image data')
#            hdr0.update('CRVAL1', hdrvec[1]['CRVAL1'], 'wavelength ref.')
#            hdr0.update('CD1_1', hdrvec[1]['CD1_1'], '')
#            hdr0.update('CRPIX1', hdrvec[1]['CRPIX1'], 'Pixel ref.')
#            #       hdr0.update('NAXIS1',hdrvec[1]['NAXIS1'],hdrvec[1].comments['NAXIS1'])
#        fits.writeto(_output, datat, hdr0)

#    floyds.util.delete('s1.fits,s2.fits,s11.fits,s22.fits')
#    for ii in range(dimension):
#        outputn = re.sub('.fits', '', _output) + '_' + str(ii + 1) + '.fits'
#        floyds.util.delete(outputn)
#    dicto = {}
#    listahed = ['ATMOR', 'ATMOB', 'SENSFUNB', 'SENSFUNR', 'ARCBLU', 'ARCRED', 'FLATRED', 'FLATBLUE',
#                'SHIFTR', 'SHIFTB', 'LAMRMS_R', 'LAMNLINR', 'SPE_ER_R', 'LAMRMS_B', 'LAMNLINB', 'SPE_ER_B',
#                'SPERES_R', 'SPERES_B']

#    for hed in listahed:
#        for hh in [hdr0, hdr1]:
#            if hed in hh:
#                try:
#                    dicto[hed] = [hh[hed], hh.comments[hed]]
#                except:
#                    dicto[hed] = [hh[hed], '']

#    if 'XMIN' in hdr0 and 'XMIN' in hdr1: _xmin = min(hdr0['XMIN'], hdr1['XMIN'])
#    if 'XMAX' in hdr0 and 'XMAX' in hdr1: _xmax = max(hdr0['XMAX'], hdr1['XMAX'])
#    dicto['XMIN'] = [_xmin, '']
#    dicto['XMAX'] = [_xmax, '']
#    dicto['GRISM'] = ['red/blu', 'full range spectrum']
#    floyds.util.updateheader(_output, 0, dicto)

#    #    iraf.scopy(img0,'s1.fits',w1='INDEF',w2=c2[argsort(a2)[0]],rebin='no')
#    #    iraf.scopy(img1,'s2.fits',w1=c1[argsort(a1)[0]],w2='INDEF',rebin='no')
#    #    iraf.sarith(input1='s1.fits',op='*',input2=A0,output='s11.fits')
#    #    iraf.sarith(input1='s2.fits',op='*',input2=A1,output='s22.fits')
#    #    iraf.specred.scombine(input='s11.fits,s22.fits',w1='INDEF',w2='INDEF',output=_output)
#    import time

#    time.sleep(1)
#    return _output


######## currently using combspec2 instead of combspec ########
def combspec2(_img0, _img1, _output, scale=True, num=None):
    import numpy as np
    import re, string, os, floyds, time
    from astropy.io import fits
    from pyraf import iraf

    # read in spectra and headers, determine third dimension of images
    _x0, _ = floyds.util.readspectrum(_img0)
    _x1, _ = floyds.util.readspectrum(_img1)
    hdr0 = fits.getheader(_img0)
    hdr1 = fits.getheader(_img1)
    if 'NAXIS3' in hdr0 and 'NAXIS3' in hdr1 and hdr0['NAXIS3'] == hdr1['NAXIS3']:
        dimension = hdr1['NAXIS3']  # will be 4 for Floyds spectra
    else:
        dimension = 1

    # figure out which is blue and which is red
    if min(_x0) < min(_x1):
        xblue = np.array(_x0)
        xred  = np.array(_x1)
        imgblue, imgred = _img0, _img1
    else:
        xblue = np.array(_x1)
        xred  = np.array(_x0)
        imgblue, imgred = _img1, _img0

    # num is the index in xred corresponding to 1st decile of the overlap
    if not num: num = np.searchsorted(xred, max(xblue)) / 10

    # rescale red & blue to match
    if scale: scomb_scale = 'median'
    else:     scomb_scale = 'none'

    floyds.util.delete('s1.fits,s2.fits')
    floyds.util.delete(_output)
    if dimension == 1:
        iraf.scopy(imgblue, 's1.fits', w1='INDEF', w2=xred[2*num], rebin='no') # take 1st-2nd deciles of blue part in overlap
        iraf.scopy(imgred, 's2.fits', w1=xred[num], w2='INDEF', rebin='no')    # take all but 1st decile of red part in overlap
        iraf.specred.scombine(input='s1.fits,s2.fits', w1='INDEF', w2='INDEF', output=_output,
                              scale=scomb_scale, sample=str(xred[num])+':'+str(xred[2*num])) # combine by averaging 2nd decile of overlap
    else:
        print 'more than one dimension'
        datavecs = []
        hdrvec = []
        for layer in np.arange(dimension)+1: # IRAF is 1-indexed
            outputn = re.sub('.fits', '', _output) + '_' + str(layer) + '.fits'
            floyds.util.delete(outputn)
            floyds.util.delete('s1.fits,s2.fits')
            iraf.scopy(imgblue + '[*,1,' + str(layer) + ']', 's1.fits', w1='INDEF', w2=xred[2*num], rebin='no')
            iraf.scopy(imgred + '[*,1,' + str(layer) + ']', 's2.fits', w1=xred[num], w2='INDEF', rebin='no')
            iraf.specred.scombine(input='s1.fits,s2.fits', w1='INDEF', w2='INDEF', output=outputn,
                                  scale=scomb_scale, sample=str(xred[num])+':'+str(xred[2*num]))
            datavec, head = fits.getdata(outputn, header=True) # these have shape (4440,)
            datavecs.append(datavec)
            hdrvec.append(head)
            floyds.util.delete(outputn)
        floyds.util.delete(_output)
        for key in ['NAXIS1', 'CRVAL1', 'CD1_1', 'CRPIX1']:
            hdr0[key] = hdrvec[0][key]
#        try:
#            hdr0.update('NAXIS1', hdrvec[0]['NAXIS1'], hdrvec[0].comments['NAXIS1'])
#            hdr0.update('CRVAL1', hdrvec[0]['CRVAL1'], hdrvec[0].comments['CRVAL1'])
#            hdr0.update('CD1_1', hdrvec[0]['CD1_1'], hdrvec[0].comments['CD1_1'])
#            hdr0.update('CRPIX1', hdrvec[0]['CRPIX1'], hdrvec[0].comments['CRPIX1'])
#        except:
#            hdr0.update('NAXIS1', hdrvec[0]['NAXIS1'], 'Width of image data')
#            hdr0.update('CRVAL1', hdrvec[0]['CRVAL1'], 'wavelength ref.')
#            hdr0.update('CD1_1', hdrvec[0]['CD1_1'], '')
#            hdr0.update('CRPIX1', hdrvec[0]['CRPIX1'], 'Pixel ref.')
        data3d = np.rollaxis(np.dstack(datavecs), 2)
        fits.writeto(_output, data3d, hdr0) # this must have shape (4, 1, 4440)
    floyds.util.delete('s1.fits,s2.fits')

    header = {}
    keywords = ['ATMOR', 'ATMOB', 'SENSFUNB', 'SENSFUNR', 'ARCBLU', 'ARCRED', 'FLATRED', 'FLATBLUE',
                'SHIFTR', 'SHIFTB', 'LAMRMS_R', 'LAMNLINR', 'SPE_ER_R', 'LAMRMS_B', 'LAMNLINB', 'SPE_ER_B',
                'SPERES_R', 'SPERES_B']
    for key in keywords:
        for hdr in [hdr0, hdr1]:
            if key in hdr:
                try: comment = hdr.comments[key]
                except: comment = ''
                header[key] = (hdr[key], comment)

    if 'XMIN' in hdr0 and 'XMIN' in hdr1: _xmin = min(hdr0['XMIN'], hdr1['XMIN'])
    if 'XMAX' in hdr0 and 'XMAX' in hdr1: _xmax = max(hdr0['XMAX'], hdr1['XMAX'])
    header['XMIN'] = [_xmin, '']
    header['XMAX'] = [_xmax, '']
    header['GRISM'] = ['red/blu', 'full range spectrum']
    floyds.util.updateheader(_output, 0, header)
    time.sleep(1)  # needed for iraf.specred.scombine to work reliably
    return _output


#############################################################
def combineblusens(imglist, imgout='pippo.fits'):
    import string
    from astropy.io import fits
    from numpy import compress, where, array, arange, float32
    from numpy import interp as ninterp
    import floyds

    img1, img2 = string.split(imglist, ',')
    xx1, yy1 = floyds.util.readspectrum(img1)
    xx2, yy2 = floyds.util.readspectrum(img2)
    yy22 = ninterp(xx2, xx1, yy1)
    xxm = compress((xx2 > 3800) & (xx2 < 4400), xx2)
    yym2 = compress((xx2 > 3800) & (xx2 < 4400), yy2)
    yym1 = compress((xx2 > 3800) & (xx2 < 4400), yy22)
    if len(where(yym2 - yym1 > 0)) > 0:
        xfix1 = xxm[where(yym2 - yym1 > 0)][0]
        xfix2 = xxm[where(yym2 - yym1 > 0)][-1]
    else:
        xfix1 = xxm[argmin(abs(yym2 - yym1))]
        xfix2 = xxm[argmin(abs(yym2 - yym1))]
    xxr = compress((xx2 > 4450), xx2)
    yyr2 = compress((xx2 > 4450), yy2)
    yyr1 = compress((xx2 > 4450), yy22)
    if len(where(yyr2 - yyr1 > 0)):
        xfix3 = xxr[where(yyr2 - yyr1 > 0)][0]
    else:
        xfix3 = xxr[argmin(abs(yyr2 - yyr1))]
    xxb = compress((xx2 < 3800), xx2)
    yyb2 = compress((xx2 < 3800), yy2)
    yyb1 = compress((xx2 < 3800), yy22)
    if len(where(yyb2 - yyb1 > 0)):
        xfix0 = xxb[where(yyb2 - yyb1 > 0)][-1]
    else:
        xfix0 = xxb[argmin(abs(yyb2 - yyb1))]
    xx30 = compress(xx2 <= xfix0, xx2)
    yy30 = compress(xx2 <= xfix0, yy22)
    xx31 = compress((xx2 > xfix0) & (xx2 <= xfix1), xx2)
    yy31 = compress((xx2 > xfix0) & (xx2 <= xfix1), yy2)
    xx32 = compress((xx2 > xfix1) & (xx2 <= xfix2), xx2)
    yy32 = compress((xx2 > xfix1) & (xx2 <= xfix2), yy22)
    xx33 = compress((xx2 > xfix2) & (xx2 <= xfix3), xx2)
    yy33 = compress((xx2 > xfix2) & (xx2 <= xfix3), yy2)
    xx34 = compress(xx2 > xfix3, xx2)
    yy34 = compress(xx2 > xfix3, yy22)
    yyfin = array(list(yy30) + list(yy31) + list(yy32) + list(yy33) + list(yy34))
    yyfin2 = ninterp(xx1, xx2, yyfin)
    xx11 = compress(xx1 <= xx2[0], xx1)
    yy12 = compress(xx1 <= xx2[0], yy1)
    xx12 = compress((xx1 > xx2[0]) & (xx1 <= xx2[-1]), xx1)
    yy13 = compress((xx1 > xx2[0]) & (xx1 <= xx2[-1]), yyfin2)
    xx13 = compress(xx1 > xx2[-1], xx1)
    yy14 = compress(xx1 > xx2[-1], yy1)
    yyfin3 = array(list(yy12) + list(yy13) + list(yy14))
    data1, hdr = fits.getdata(img1, 0, header=True)
    floyds.util.delete(imgout)
    data1 = array(yyfin3)
    fits.writeto(imgout, float32(data1), hdr)
    return imgout
    ##################################################

def combineredsens(imglist, imgout='pippo.fits'): # used for both arms of FTS
    from pyraf import iraf
    import floyds

    floyds.util.delete(imgout)
    sss = iraf.specred.scombine(imglist, output=imgout, combine='average', reject='none',
                                scale='none', weight='none', Stdout=1)
    return imgout

    ##################################################


def rectify_single_image(img, imgrect, imgrect1, fcuntilt_file, xa, xb, ya, yb, lambda1, lambda2, y2, _cosmic=False):
    import floyds
    from numpy import float32
    from astropy.io import fits
    from pyraf import iraf
    import re
    import os

    first_rectified_image = 't' + img
    floyds.util.delete(first_rectified_image)
    iraf.specred.transform(input=img, output=first_rectified_image, minput='', fitnames=re.sub('.fits', '', imgrect),
                           databas='database', x1='INDEF', x2='INDEF', dx=1, y1='INDEF', y2=y2, dy=1,
                           flux='yes', blank=0, logfile='logfile')

    data, hdr = fits.getdata(first_rectified_image, 0, header=True)

    if not hdr.get('HDRVER'):
        fits.writeto(first_rectified_image, float32(data[0][ya:yb, xa:xb]), hdr, overwrite=True)
    else:
        fits.writeto(first_rectified_image, float32(data[ya:yb, xa:xb]), hdr, overwrite=True)

    iraf.hedit(first_rectified_image, 'CCDSEC', delete='yes', update='yes', verify='no', Stdout=1)
    iraf.hedit(first_rectified_image, 'TRIM', delete='yes', update='yes', verify='no', Stdout=1)

    if _cosmic:
        _gain = floyds.util.readkey3(hdr, 'gain')
        _rdnoise = floyds.util.readkey3(hdr, 'ron')
        floyds.cosmics.lacos(first_rectified_image, output='', gain=_gain, readn=_rdnoise, xorder=9, yorder=9,
                             sigclip=4.5, sigfrac=0.5, objlim=1, verbose=True, interactive=True)
        floyds.util.updateheader(first_rectified_image, 0, {'LACOSMIC': [True, 'Laplacian cosmic ray rejection']})
        print '\n### cosmic rays rejections ........ done '
    else:
        floyds.util.updateheader(first_rectified_image, 0, {'LACOSMIC': [False, 'Laplacian cosmic ray rejection']})

    untilted_image = 'u' + first_rectified_image
    iraf.specred.transform(input=first_rectified_image, output=untilted_image, fitnames=os.path.basename(fcuntilt_file)[2:],
                           databas='database', x1='INDEF', x2='INDEF', dx=1,
                           y1='INDEF', y2='INDEF', dy=1, flux='yes', blank=0,
                           logfile='logfile')
    wavelength_rectified_image = 'tt' + img
    floyds.util.delete(wavelength_rectified_image)
    iraf.specred.transform(input=untilted_image, output=wavelength_rectified_image, minput='',
                           fitnames=re.sub('.fits', '', imgrect1),
                           databas='database', x1=lambda1, x2=lambda2, dx='INDEF',
                           y1='INDEF', y2='INDEF', dy=1, flux='yes', blank=0,
                           logfile='logfile')  # , mode='h')
    floyds.util.updateheader(wavelength_rectified_image, 0, {'DISPAXIS': [1, 'dispersion axis'],
                                                             'CUNIT1': ['Angstrom', 'Units of dispersion axis']})
    floyds.util.delete(first_rectified_image)
    floyds.util.delete(untilted_image)


def rectifyspectrum(img, arcfile, flatfile, fcfile, fcfile1, fcfile_untilt, _interactive=True, _cosmic=True):
    import floyds
    from floyds.util import delete, updateheader, readhdr, readkey3, display_image
    import string, re, os, glob, sys
    from numpy import array, arange, argmin, float32
    from astropy.io import fits
    from pyraf import iraf
    import datetime

    os.environ["PYRAF_BETA_STATUS"] = "1"
    iraf.set(direc=floyds.__path__[0] + '/')
    dv = floyds.util.dvex()
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['ccdred.flatcombine', 'ccdred.zerocombine', 'ccdproc', 'specred.apall', 'longslit.identify',
                'longslit.reidentify',
                'specred.standard', 'longslit.fitcoords', 'specred.transform', 'specred.response']
    for t in toforget: iraf.unlearn(t)
    iraf.longslit.dispaxi = 2
    iraf.longslit.mode = 'h'
    iraf.identify.fwidth = 7
    iraf.identify.order = 2
    iraf.specred.dispaxi = 2
    iraf.specred.mode = 'h'
    iraf.ccdproc.darkcor = 'no'
    iraf.ccdproc.fixpix = 'no'
    iraf.ccdproc.trim = 'no'
    iraf.ccdproc.flatcor = 'no'
    iraf.ccdproc.overscan = 'no'
    iraf.ccdproc.zerocor = 'no'
    iraf.ccdproc.biassec = ''
    iraf.ccdproc.ccdtype = ''
    iraf.ccdred.verbose = 'yes'
    iraf.specred.verbose = 'yes'

    imgrect = string.split(fcfile, '/fc')[-1] + '.fits'
    imgrect1 = string.split(fcfile1, '/fc')[-1] + '.fits'

    Path('database/').mkdir(parents=True, exist_ok=True)
    os.system('cp ' + fcfile + ' database/')
    os.system('cp ' + fcfile1 + ' database/')
    os.system('cp ' + fcfile_untilt + ' database/')

    ###################################    transform img, arc, flat
    data, hdr = fits.getdata(img, 0, header=True)
    _arm = floyds.util.readkey3(hdr, 'GRISM')
    _slit = floyds.util.readkey3(hdr, 'slit')
    _arm = floyds.util.readkey3(hdr, 'GRISM')
    _tel = floyds.util.readkey3(hdr, 'TELID')
    camera = floyds.util.readkey3(hdr, 'INSTRUME')

    if _tel not in ['ftn', 'fts']:     _tel = floyds.util.readkey3(hdr, 'SITEID')
    if _arm == 'red':
        lambda1, lambda2 = 4850.0, 10180.0
        if _tel in ['ftn', 'ogg']:
            if camera == 'en06':
                xa, xb = 0, 1792
                ya, yb = 231, 322
                y2 = 'INDEF'
            else:
                raise ValueError('Camera not supported by pipeline')
        else:
            if camera == 'en05':
                xa, xb = 0, 1800
                ya, yb = 186, 285
            elif camera == 'en12':
                xa, xb = 0, 1792
                ya, yb = 221, 312
                y2 = 'INDEF'
            else:
                raise ValueError('Camera not supported by pipeline')
    else:
        lambda1, lambda2 = 3300.0, 5700.0
        if _tel in ['ftn', 'ogg']:
            if camera == 'en06':
                xa, xb = 104, 1515
                ya, yb = 131, 224
                y2 = 100
            else:
                raise ValueError('Camera not supported by pipeline')
        else:
            if camera == 'en05':
                xa, xb = 0, 1669
                ya, yb = 125, 226
            elif camera == 'en12':
                xa, xb = 206, 1587
                ya, yb = 182, 273
                y2 = 100
            else:
                raise ValueError('Camera not supported by pipeline')

    rectify_single_image(img, imgrect, imgrect1, fcfile_untilt, xa, xb, ya, yb, lambda1, lambda2, y2, _cosmic=_cosmic)

    if arcfile:
        rectify_single_image(arcfile, imgrect, imgrect1, fcfile_untilt, xa, xb, ya, yb, lambda1, lambda2, y2, _cosmic=_cosmic)
        output_arcfile = 'tt' + arcfile
    else:
        output_arcfile = ''
    if flatfile:
        rectify_single_image(flatfile, imgrect, imgrect1, fcfile_untilt, xa, xb, ya, yb, lambda1, lambda2, y2, _cosmic=False)
        output_flatfile = 'tt' + flatfile
    else:
        output_flatfile = ''
    return 'tt' + img, output_arcfile, output_flatfile


#############################################################33333
def fringing_classicmethod(flatfile, img, _inter, _sample, _order, arm):
    from pyraf import iraf
    import floyds
    import re, sys, string, os
    import numpy as np
    from astropy.io import fits

    floyds.util.delete('n' + flatfile)
    floyds.util.delete('norm.fits')
    floyds.util.delete('n' + img)
    datax, hdrx = fits.getdata(flatfile, 0, header=True)
    xdim = hdrx['NAXIS1']
    ydim = hdrx['NAXIS2']
    if arm == 'red':
        floyds.util.delete(re.sub('.fits', 'c.fits', flatfile))
        iraf.imcopy(flatfile + '[350:' + str(xdim) + ',*]', re.sub('.fits', 'c.fits', flatfile), verbose='no')
        flatfile = re.sub('.fits', 'c.fits', flatfile)
        floyds.util.delete('n' + flatfile)
        iraf.specred.response(flatfile, normaliz=flatfile + '[*,1:' + str(int(ydim) - 1) + ']',
                              response='n' + flatfile, interac=_inter, thresho='INDEF', sample=_sample, naverage=2,
                              function='spline3', low_rej=3, high_rej=3, order=_order, niterat=20, grow=0,
                              graphic='stdgraph')
        iraf.imarith(img, '/', img, 'norm.fits', verbose='no')
        iraf.imcopy('n' + flatfile, 'norm.fits[350:' + str(xdim) + ',*]', verbose='no')
    else:
        iraf.specred.response(flatfile, normaliz=flatfile + '[*,1:' + str(int(ydim) - 1) + ']',
                              response='norm.fits', interac=_inter, thresho='INDEF', sample=_sample, naverage=2,
                              function='spline3', low_rej=3, high_rej=3, order=40, niterat=20, grow=0,
                              graphic='stdgraph')
    floyds.util.delete('n' + flatfile)
    floyds.util.delete('n' + img)
    iraf.imrename('norm.fits', 'n' + flatfile, verbose='no')
    print img, 'n' + flatfile, 'n' + img
    floyds.util.delete('n' + img)
    data, hdr = fits.getdata(img, 0, header=True)
    datan, hdrn = fits.getdata('n' + flatfile, 0, header=True)
    #########
    mask = datan == 0
    datan[mask] = 1
    data2 = np.float32(data / datan)
    #########
    _grism = hdr['grism']
    fits.writeto('n' + img, data2, hdr)
    floyds.util.updateheader('n' + img, 0, {'FLAT' + str(_grism)[0]: ['n' + flatfile, 'flat file']})
    return 'n' + img


####################################################################
def applyflat(img, flatimg, output='', scale='', shift=''):
    from astropy.io import fits
    import floyds
    import string
    import numpy as np
    from pyraf import iraf
    from iraf import specred

    data, hdr = fits.getdata(img, 0, header=True)
    datan, hdrn = fits.getdata(flatimg, 0, header=True)
    if 'GRISM' in hdr:
        arm = hdr['GRISM']
    else:
        arm = ''
    if scale == '' and shift == '':
        y = data.mean(1)
        x = np.arange(len(y))
        if np.argmax(y) < 80 and np.argmax(y) > 15:
            y2 = data[np.argmax(y) - 3:np.argmax(y) + 3].mean(0)
            yn = datan[np.argmax(y) - 3:np.argmax(y) + 3].mean(0)
            yy2 = data[np.argmax(y) - 9:np.argmax(y) - 3].mean(0)
            x2 = np.arange(len(y2))
            floyds.delete('_sky1.fits')
            floyds.delete('_sky2.fits')
            floyds.delete('sky1.fits')
            floyds.delete('sky2.fits')
            floyds.delete('tsky1.fits')
            floyds.delete('tsky11.fits')
            fits.writeto('_sky1.fits', np.float32(y2 - yy2), hdr)
            fits.writeto('_sky2.fits', np.float32(yn), hdr)
            iraf.scopy('_sky1.fits', 'sky1.fits', w1=7000, w2=10000)
            iraf.scopy('_sky2.fits', 'sky2.fits', w1=7000, w2=10000)
            iraf.specred.continuum('sky1.fits', output='tsky1.fits', line='*', type='difference',
                                   interact='yes', function='spline3', niterat=100, low_rej=4, high_re=5,
                                   sample='7000:7550,7730:10000', order=10, ask='no')
            iraf.specred.continuum('sky1.fits', output='tsky11.fits', type='fit',
                                   interact='yes', function='spline3', niterat=100, low_rej=4, high_re=5,
                                   sample='7000:7550,7730:10000', order=10, ask='no')

            data1, hdr1 = fits.getdata('tsky1.fits', 0, header=True)  #  fringing on the standrad
            yy1 = data1
            xx1 = np.arange(len(yy1))
            data11, hdr11 = fits.getdata('sky1.fits', 0, header=True)  #  standard
            yy11 = data11
            xx11 = np.arange(len(yy11))
            data111, hdr111 = fits.getdata('tsky11.fits', 0, header=True)  #  continuum
            yy111 = data111
            xx111 = np.arange(len(yy111))
            data2, hdr2 = fits.getdata('sky2.fits', 0, header=True)  #  fringing form flat
            yy2 = data2
            xx2 = np.arange(len(yy2))
            shift = floyds.floydsspecdef.checkwavelength_arc(xx1, yy1, xx2, yy2, '', '', False) * (
            -1.)  #   shift fringing
            if abs(shift) > 10: shift = 0
            xx3 = np.array(xx2, float)[:]
            yy3 = yy2[:]
            xx4 = xx3 + shift
            yy3 = np.interp(xx1, xx4, yy3)
            _scaleo2 = []
            integral_o2 = []
            for i in range(1, 381):
                j = -0.6 + i * 0.04
                ll = abs((yy11 / ((yy3 * j) - (j - 1))) - yy111)
                integraleo2 = np.trapz(ll, xx11)
                integral_o2.append(integraleo2)
                _scaleo2.append(j)
            scale = _scaleo2[np.argmin(integral_o2)]
            if abs(scale) >= 2 or abs(scale) <= 0.5: scale = 1
        else:
            scale = 1
            shift = 0
    else:
        print 'use ' + str(scale) + ' ' + str(shift)
    ####################################
    datanew = datan * scale + (1 - scale)  # scale the normalizedflat
    if int(shift) >= 0:  # shift the flat
        datanew[:, int(shift):len(datanew[0])] = datanew[:, 0:len(datanew[0]) - int(shift)]
    else:
        datanew[:, 0:len(datanew[0]) - int(abs(shift))] = datanew[:, int(abs(shift)):len(datanew[0])]
    datanew2 = data / datanew
    if not output: output = 'NN' + img
    floyds.delete(output)
    fits.writeto(output, np.float32(datanew2), hdr)
    floyds.util.updateheader(output, 0, {'FRSCALE': [scale, 'fringing scale factor ']})
    floyds.util.updateheader(output, 0, {'FRSHIFT': [shift, 'fringing shift factor ']})
    floyds.util.updateheader(output, 0, {'FLAT' + arm: [string.split(flatimg, '/')[-1], 'flat field file ']})
    floyds.delete('_sky1.fits')
    floyds.delete('_sky2.fits')
    floyds.delete('sky1.fits')
    floyds.delete('sky2.fits')
    floyds.delete('tsky1.fits')
    floyds.delete('tsky11.fits')
    return output


#####################################################
#############################################################33333
def fringing_classicmethod2(flatfile, img, _inter, _sample, _order, arm):
    from pyraf import iraf

    iraf.specred(_doprint=0)
    import floyds
    import re, sys, string, os
    from numpy import float32
    from astropy.io import fits

    datax, hdrx = fits.getdata(flatfile, 0, header=True)
    xdim = hdrx['NAXIS1']
    ydim = hdrx['NAXIS2']
    floyds.util.delete('n' + flatfile)
    floyds.util.delete('norm.fits')
    floyds.util.delete('n' + img)
    iraf.specred.apedit.nsum = 15
    iraf.specred.apedit.width = 100.
    iraf.specred.apedit.line = 1024
    iraf.specred.apfind.minsep = 20.
    iraf.specred.apfind.maxsep = 1000.
    iraf.specred.apresize.bkg = 'no'
    iraf.specred.apresize.ylevel = 0.5
    iraf.specred.aptrace.nsum = 10
    iraf.specred.aptrace.step = 10
    iraf.specred.aptrace.nlost = 10
    if arm == 'red':
        floyds.util.delete(re.sub('.fits', 'c.fits', flatfile))
        iraf.imcopy(flatfile + '[350:' + str(xdim) + ',*]', re.sub('.fits', 'c.fits', flatfile), verbose='no')
        iraf.imarith(flatfile, '/', flatfile, 'norm.fits', verbose='no')
        flatfile = re.sub('.fits', 'c.fits', flatfile)
        floyds.util.delete('n' + flatfile)
        iraf.unlearn(iraf.specred.apflatten)
        floyds.floydsspecdef.aperture(flatfile)
        iraf.specred.apflatten(flatfile, output='n' + flatfile, interac=_inter, find='no', recenter='no', resize='no',
                               edit='no', trace='no',
                               fittrac='no', fitspec='no', flatten='yes', aperture='',
                               pfit='fit2d', clean='no', function='legendre', order=_order, sample='*', mode='ql')
        iraf.imcopy('n' + flatfile, 'norm.fits[350:' + str(xdim) + ',*]', verbose='no')
    else:
        iraf.specred.response(flatfile, normaliz=flatfile + '[*,1:' + str(int(ydim) - 1) + ']',
                              response='norm.fits', interac=_inter, thresho='INDEF', sample=_sample, naverage=2,
                              function='spline3', low_rej=3, high_rej=3, order=40, niterat=20, grow=0,
                              graphic='stdgraph')
    floyds.util.delete('n' + flatfile)
    floyds.util.delete('n' + img)
    iraf.imrename('norm.fits', 'n' + flatfile, verbose='no')
    print img, 'n' + flatfile, 'n' + img
    floyds.util.delete('n' + img)
    data, hdr = fits.getdata(img, 0, header=True)
    datan, hdrn = fits.getdata('n' + flatfile, 0, header=True)
    fits.writeto('n' + img, float32(data / datan), hdr)
    floyds.util.updateheader('n' + img, 0, {'FLAT' + arm: ['n' + flatfile, 'flat file']})
    return 'n' + img


####################################################################
def aperture(img):
    from astropy.io import fits
    import re
    import time

    hdr = fits.open(img)[0].header
    xmax = hdr['NAXIS2']
    xmin = -39
    img2 = re.sub('.fits', '', img)
    line = "# Sun 13:10:40 16-Jun-2013\nbegin	aperture " + img2 + " 1 650. 40.0\n" + \
           "	 image	" + img2 + "\n	aperture	1\n	beam	1\n	center	650. 40.0\n" + \
           "	 low	-950. " + str(xmin) + "\n	high	950. " + str(xmax) + "\n" \
                                                                                  "	 background\n	 xmin -110.\n" + \
           "	 xmax 110.\n	 function chebyshev\n		order 1\n		sample *\n" + \
           "	 naverage -3\n	 niterate 0\n		low_reject 3.\n		high_reject 3.\n" + \
           "	 grow 0.\n	 axis	2\n	 curve	5\n		2.\n		1.\n" + \
           "	 10.\n		840.\n		-7.418439\n"
    f = open('database/ap' + img2, 'w')
    f.write(line)
    f.close()
    time.sleep(1)

##########################################################################

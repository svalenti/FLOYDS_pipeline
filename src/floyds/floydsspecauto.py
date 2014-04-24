def fluxcalib2d(img2d,sensfun):  # flux calibrate 2d images
    import pyfits
    import re,string
    from pyfits import open as popen
    from numpy import arange, array, float32
    from numpy import interp as ninterp
    import floyds
    from floyds.util import readhdr,readkey3,delete

    _tel=floyds.util.readkey3(floyds.util.readhdr(re.sub('\n','',img2d)),'TELID')
    if _tel=='fts':
        _extinction='ssoextinct.dat'
        _observatory='sso'
    elif _tel=='ftn':
        _extinction='maua.dat' 
        _observatory='cfht'
    else: sys.exit('ERROR: observatory not recognised')

    data2d, hdr2d = pyfits.getdata(img2d, 0, header=True)
    if 'GRISM' in hdr2d:  _arm=hdr2d['GRISM']
    else:                 _arm=''
    xxd=arange(len(data2d[0,:]))
    crvald=popen(img2d)[0].header.get('CRVAL1')
    cdd=popen(img2d)[0].header.get('CD1_1')
    _exptime,_airmass=readkey3(readhdr(img2d),'exptime'),readkey3(readhdr(img2d),'airmass')
    #  read sensfunction and interpole pixel of 2D image
    yys=popen(sensfun)[0].data
    crvals=popen(sensfun)[0].header.get('CRVAL1')
    cds=popen(sensfun)[0].header.get('CD1_1')    
    yys=(10**(yys/2.5))*cds                  # from sens iraf in sens flux
    xxs=arange(len(yys))
    aasens=crvals+(xxs)*cds
    xxs2=(aasens-crvald)/cdd
    aasens2=ninterp(xxd,xxs2,yys)
    #  read atmosferic function and interpole pixel of 2D image
    aae,yye=floyds.util.ReadAscii2(floyds.__path__[0]+'/standard/extinction/'+_extinction)
    aae,yye=array(aae,float),array(yye,float)
    xxe=(aae-crvald)/cdd
    atm_xx=ninterp(xxd,xxe,yye)
    aircorr=10**(0.4*array(atm_xx)*_airmass)
    img2df=re.sub('.fits','_2df.fits',img2d)
    for i in range(len(data2d[:,0])):  data2d[i,:]=((array(data2d[i,:]/_exptime)*array(aircorr))/aasens2)*1e20
    floyds.util.delete(img2df)
    pyfits.writeto(img2df, float32(data2d), hdr2d)
    floyds.util.updateheader(img2df,0,{'SENSFUN'+_arm[0]:[string.split(sensfun,'/')[-1],'']})
    floyds.util.updateheader(img2df,0,{'BUNIT':['erg/cm2/s/A  10^20','Physical unit of array values']})
    return img2df

def gettar(img):
    import pyfits
    import floyds
    from urllib import urlopen
    import re,string,os
    from datetime import datetime, timedelta
    data, hdr = pyfits.getdata(img, 0, header=True)

    #
    #
    #   need to be change depending what we are doing downloding from floyds machine
    #
    #
    # print "DEBUG: img=", img
    if '_' in img:
        # Old-style filenames
        imgg=re.sub(string.split(img,'_')[3],re.sub('0','',string.split(img,'_')[3]),img)
        img1=re.sub(string.split(imgg,'_')[4],re.sub('0','',string.split(imgg,'_')[4]),imgg)
    else:
        imgg = ''
        img1 = img.replace('e00', 'e02')
    # print "DEBUG: imgg,img1=", imgg,img1

    # FTN filenames will be the day prior
    _tel=hdr['TELID']
    if '2m0a' in _tel:
        _tel = hdr['SITEID']
    _tel = _tel.lower()
    if _tel == 'fts' or _tel == 'coj': delta=0
    else:          delta=1

    # Extract date from header. We only want the whole seconds part (not sure
    # why, could be read by .%f in strptime)
    _date=floyds.readkey3(hdr,'DATE-OBS')
    _date_wholesecs = _date
    # Check if there is a '.' (Previous version will break if the header already
    # has whole seconds and there is no period)
    if '.' in _date:
      _date_wholesecs = _date.split('.')[0]
    # Parse the string, subtract off any delta and reformat to YYYYMMDD (without hyphens)
    _date_dt = datetime.strptime(_date_wholesecs,"%Y-%m-%dT%H:%M:%S")
    _date_dt = _date_dt - timedelta(days=delta)
    p=_date_dt.strftime('%Y%m%d')

    obj=hdr['object']
    propid=hdr['PROPID']
    if 'fts' in _tel or 'coj' in _tel: 
        i=r'http://floyds.coj.lco.gtn/night_summary/%s/' % (p)
    else:
        i=r'http://floyds.ogg.lco.gtn/night_summary/%s/' % (p)

    try:
        webpage=urlopen(i).read()
        for data_dir in webpage.split('href='):
            if propid in data_dir:
                dir_name = data_dir.split('>')[0]
                dir_name = dir_name.replace('"','')
                dir_url = i+dir_name
                cc=urlopen(dir_url).read()
                if img1 in cc:
                    for jj in string.split(cc,'href='):
                        if '.tar' in jj:
                            # Store url of directory where we found tarfile so
                            # we can get to HTML file later
                            dir_url_for_tar = dir_url
                            tar_name=jj.split('>')[0].replace('"','')
                            if os.path.isfile(tar_name): floyds.util.delete(tar_name)
                            com='wget -nv %s%s ' % (dir_url, tar_name)
        try:
            os.system(com)
        except:
            tar_name=''

        # print "DEBUG: ", tar_name, dir_url, dir_url_for_tar
        htmlfile_url = dir_url_for_tar + tar_name.replace('.tar', '.html')
        pdffile = tar_name.replace('.tar', '.pdf')
        # Fetch html file from remote URL and convert to PDF
        line='xhtml2pdf  ' + htmlfile_url + ' ' + pdffile
        # print "DEBUG: ", line
        try:    
            os.system(line)
        except:      pdffile=''
    except:  tar_name,i,dir_url,line,pdffile='','','','',''
    return tar_name,pdffile


def floydsautoredu(files,_interactive,_dobias,_doflat,_listflat,_listbias,_listarc,_cosmic,_ext_trace,_dispersionline,liststandard,listatmo,_automaticex,_classify=False,_verbose=False,smooth=1,fringing=1):
    import floyds
    import string,re,os,glob,sys,pickle
    from numpy import array, arange, mean,pi,arccos,sin,cos,argmin
    import pyfits
    from pyraf import iraf
    import datetime
    os.environ["PYRAF_BETA_STATUS"] = "1"
    iraf.set(direc=floyds.__path__[0]+'/')
    _extinctdir='direc$standard/extinction/'
    _tel=floyds.util.readkey3(floyds.util.readhdr(re.sub('\n','',files[0])),'TELID')
    if _tel=='fts':
        _extinction='ssoextinct.dat'
        _observatory='sso'
    elif _tel=='ftn':
        _extinction='maua.dat' 
        _observatory='cfht'   
    else: sys.exit('ERROR: observatory not recognised')
    dv=floyds.util.dvex()
    scal=pi/180.
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    iraf.twodspec(_doprint=0)
    iraf.longslit(_doprint=0)
    iraf.specred(_doprint=0)
    toforget = ['ccdred.flatcombine','ccdred.zerocombine','ccdproc','specred.apall','longslit.identify','longslit.reidentify',\
                    'specred.standard','longslit.fitcoords','specred.transform','specred.response']
    for t in toforget: iraf.unlearn(t)
    iraf.longslit.dispaxi=2
    iraf.longslit.mode='h'
    iraf.identify.fwidth=7 
    iraf.identify.order=2 
    iraf.specred.dispaxi=2
    iraf.specred.mode='h'
    iraf.ccdproc.darkcor='no'
    iraf.ccdproc.fixpix='no'
    iraf.ccdproc.trim='no'
    iraf.ccdproc.flatcor='no'
    iraf.ccdproc.overscan='no'
    iraf.ccdproc.zerocor='no'
    iraf.ccdproc.biassec=''
    iraf.ccdproc.ccdtype=''
    iraf.ccdred.instrument = "/dev/null"
    if _verbose: 
        iraf.ccdred.verbose='yes'
        iraf.specred.verbose='yes'
    else: 
        iraf.specred.verbose='no'
        iraf.ccdred.verbose='no'
    now=datetime.datetime.now()
    datenow=now.strftime('20%y%m%d%H%M')
    MJDtoday=55928+(datetime.date.today()-datetime.date(2012, 01, 01)).days
    outputlist=[]
    hdra=floyds.util.readhdr(re.sub('\n','',files[0]))
    _gain=floyds.util.readkey3(hdra,'gain')
    _rdnoise=floyds.util.readkey3(hdra,'ron')
    std,rastd,decstd,magstd=floyds.util.readstandard('standard_floyds_mab.txt')
    _naxis2=hdra.get('NAXIS2')
    _naxis1=hdra.get('NAXIS1')
    if not _naxis1: _naxis1=2079
    if not _naxis2: 
        if not hdr0.get('HDRVER'):   _naxis1=511
        else:                        _naxis1=512
    _overscan='[2049:'+str(_naxis1)+',1:'+str(_naxis2)+']'
    _biassecblu='[380:2048,325:'+str(_naxis2)+']'    
    _biassecred='[1:1800,1:350]'    
    lista={}
    objectlist={}
    biaslist={}
    flatlist={}
    flatlistd={}
    arclist={}
    for img in files:
        hdr0=floyds.util.readhdr(img)
        if  floyds.util.readkey3(hdr0,'naxis2')>=500:
            if 'blu' not in lista: lista['blu']=[]
            if 'red' not in lista: lista['red']=[]
            _object0=floyds.util.readkey3(hdr0,'object')
            _object0=re.sub(' ' ,'',_object0)
            _object0=re.sub('/' ,'_',_object0)
            _object0=re.sub('\(' ,'',_object0)
            _object0=re.sub('\[' ,'',_object0)
            _object0=re.sub('\)' ,'',_object0)
            _object0=re.sub('\]' ,'',_object0)
            _date0=floyds.util.readkey3(hdr0,'date-night')
            _tel=floyds.util.readkey3(hdr0,'TELID')
            _type=floyds.util.readkey3(hdr0,'OBSTYPE')
            if not _type:    _type=floyds.util.readkey3(hdr0,'imagetyp')
            _slit=floyds.util.readkey3(hdr0,'slit')
            if _type:
                _type = _type.lower()
                if _type in ['sky','spectrum','expose']:
                    nameoutb=str(_object0)+'_'+_tel+'_'+str(_date0)+'_blue_'+str(_slit)+'_'+str(MJDtoday)
                    nameoutr=str(_object0)+'_'+_tel+'_'+str(_date0)+'_red_'+str(_slit)+'_'+str(MJDtoday)
                elif _type in ['lamp','arc','l']:
                    nameoutb='arc_'+str(_object0)+'_'+_tel+'_'+str(_date0)+'_blue_'+str(_slit)+'_'+str(MJDtoday)
                    nameoutr='arc_'+str(_object0)+'_'+_tel+'_'+str(_date0)+'_red_'+str(_slit)+'_'+str(MJDtoday)
                elif _type in ['flat','f','lampflat','lamp-flat']:
                    nameoutb='flat_'+str(_object0)+'_'+_tel+'_'+str(_date0)+'_blue_'+str(_slit)+'_'+str(MJDtoday)
                    nameoutr='flat_'+str(_object0)+'_'+_tel+'_'+str(_date0)+'_red_'+str(_slit)+'_'+str(MJDtoday)
                else:
                    nameoutb=str(_type.lower())+'_'+str(_object0)+'_'+_tel+'_'+str(_date0)+'_blue_'+str(_slit)+'_'+str(MJDtoday)
                    nameoutr=str(_type.lower())+'_'+str(_object0)+'_'+_tel+'_'+str(_date0)+'_red_'+str(_slit)+'_'+str(MJDtoday)

                bimg=floyds.util.name_duplicate(img,nameoutb,'')
                rimg=floyds.util.name_duplicate(img,nameoutr,'')
####
                floyds.util.delete(bimg)
                floyds.util.delete(rimg)
                iraf.imcopy(img,bimg,verbose='no')
                iraf.imcopy(img,rimg,verbose='no')

                aaa=iraf.hedit(bimg,'CCDSEC',delete='yes',update='yes',verify='no',Stdout=1)
                aaa=iraf.hedit(bimg,'TRIMSEC',delete='yes',update='yes',verify='no',Stdout=1)
                aaa=iraf.hedit(rimg,'CCDSEC',delete='yes',update='yes',verify='no',Stdout=1)
                aaa=iraf.hedit(rimg,'TRIMSEC',delete='yes',update='yes',verify='no',Stdout=1)

                iraf.ccdproc(bimg,output='', overscan="yes", trim="yes", zerocor='no', flatcor='no', zero='', ccdtype='',\
                                 fixpix='no', trimsec=_biassecblu, biassec=_overscan, readaxi='line', Stdout=1)
                iraf.ccdproc(rimg,output='', overscan="yes", trim="yes", zerocor='no', flatcor='no', zero='', ccdtype='',\
                                 fixpix='no', trimsec=_biassecred, biassec=_overscan, readaxi='line', Stdout=1)
                floyds.util.updateheader(bimg,0,{'GRISM':['blu',' blue order']})
                floyds.util.updateheader(rimg,0,{'GRISM':['red',' blue order']})
                floyds.util.updateheader(bimg,0,{'arcfile':[img,'file name in the archive']})
                floyds.util.updateheader(rimg,0,{'arcfile':[img,'file name in the archive']})
                lista['blu'].append(bimg)
                lista['red'].append(rimg)
            else: 
                print 'warning type not defined'
    for arm in lista.keys():
        for img in lista[arm]:
            print img
            hdr=floyds.util.readhdr(img)
            _type=floyds.util.readkey3(hdr,'OBSTYPE')
            if _type=='EXPOSE':  
                      _type=floyds.util.readkey3(hdr,'imagetyp')
                      if not _type: _type='EXPOSE'

            if _type=='EXPOSE':  
                print 'warning obstype still EXSPOSE, are this old data ?  run manually floydsfixheader'

            _slit=floyds.util.readkey3(hdr,'slit')
            _grpid=floyds.util.readkey3(hdr,'grpid')
            if _type.lower() in ['flat','f','lamp-flat'] :
                if (arm,_slit) not in flatlist:  flatlist[arm,_slit]={}
                if _grpid not in flatlist[arm,_slit]: flatlist[arm,_slit][_grpid]=[img]
                else: flatlist[arm,_slit][_grpid].append(img)
            elif _type.lower() in ['lamp','l','arc']:
                if (arm,_slit) not in arclist:  arclist[arm,_slit]={}
                if _grpid not in arclist[arm,_slit]: arclist[arm,_slit][_grpid]=[img]
                else: arclist[arm,_slit][_grpid].append(img)
            elif _type in ['bias','b']:
                if arm not in biaslist: biaslist[arm]=[]
                biaslist[arm].append(img)
            elif _type.lower() in ['sky','s','spectrum']:
                try:
                    _ra=float(floyds.util.readkey3(hdr,'RA'))
                    _dec=float(floyds.util.readkey3(hdr,'DEC'))
                except:
                    ra00=string.split(floyds.util.readkey3(hdr,'RA'),':')
                    ra0,ra1,ra2=float(ra00[0]),float(ra00[1]),float(ra00[2])
                    _ra=((ra2/60.+ra1)/60.+ra0)*15.
                    dec00=string.split(floyds.util.readkey3(hdr,'DEC'),':')
                    dec0,dec1,dec2=float(dec00[0]),float(dec00[1]),float(dec00[2])
                    if '-' in str(dec0):       _dec=(-1)*((dec2/60.+dec1)/60.+((-1)*dec0))
                    else:                      _dec=(dec2/60.+dec1)/60.+dec0
                dd=arccos(sin(_dec*scal)*sin(decstd*scal)+cos(_dec*scal)*cos(decstd*scal)*cos((_ra-rastd)*scal))*((180/pi)*3600)
                if _verbose:
                    print _ra,_dec
                    print std[argmin(dd)],min(dd)
                if min(dd)<5200: _typeobj='std'
                else: _typeobj='obj'
                if min(dd)<5200:
                    floyds.util.updateheader(img,0,{'stdname':[std[argmin(dd)],'']})
                    floyds.util.updateheader(img,0,{'magstd':[float(magstd[argmin(dd)]),'']})
                if _typeobj not in objectlist:      objectlist[_typeobj]={}

                if (arm,_slit) not in objectlist[_typeobj]:     objectlist[_typeobj][arm,_slit]=[img]
                else: objectlist[_typeobj][arm,_slit].append(img)
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

    sens={}
    outputfile={}
    atmo={}
    for tpe in objectlist:
      if tpe not in outputfile:  outputfile[tpe]={}
      for setup in objectlist[tpe]:
        if setup not in sens:   sens[setup]=[]
        print '\n### setup= ',setup,'\n### objects= ',objectlist[tpe][setup],'\n'
        for img in objectlist[tpe][setup]:
              print '\n\n### next object= ',img,' ',floyds.util.readkey3(floyds.util.readhdr(img),'object'),'\n'
              hdr=floyds.util.readhdr(img)
              archfile=floyds.util.readkey3(hdr,'arcfile')
              _gain=floyds.util.readkey3(hdr,'gain')
              _rdnoise=floyds.util.readkey3(hdr,'ron')
              _grism=floyds.util.readkey3(hdr,'grism')
              _grpid=floyds.util.readkey3(hdr,'grpid')
              if archfile not in outputfile[tpe]: outputfile[tpe][archfile]=[]
#####################      flat   ###############
              if _listflat:   flatgood=_listflat    # flat list from reducer
              elif setup in flatlist:  
                  if _grpid in flatlist[setup]:
                      print '\n###FLAT WITH SAME GRPID'
                      flatgood= flatlist[setup][_grpid]     # flat in the  raw data
                  else:  
                      flatgood=[]
                      for _grpid0 in flatlist[setup].keys():
                          for ii in flatlist[setup][_grpid0]:
                              flatgood.append(ii)
              else: flatgood=[]
              if len(flatgood)!=0:
                  if len(flatgood)>1:
                      f=open('_oflatlist','w')
                      for fimg in flatgood:
                          print fimg
                          f.write(fimg+'\n')
                      f.close()
                      floyds.util.delete('flat'+img)
                      iraf.ccdred.flatcombine('"@_oflatlist"',output='flat'+img,combine='average',reject='none',ccdtype=' ',rdnoise=_rdnoise,gain=_gain, process='no', Stdout=1)
                      floyds.util.delete('_oflatlist')
                      flatfile='flat'+img
                  elif len(flatgood)==1:
                      os.system('cp '+flatgood[0]+' flat'+img)
                      flatfile='flat'+img
              else: flatfile=''
##########################   find arcfile            #######################
              arcfile=''
              if _listarc:       arcfile= [floyds.util.searcharc(img,_listarc)[0]][0]   # take arc from list 
              if not arcfile and setup in arclist.keys():
                    if _grpid in arclist[setup]:  
                        print '\n###ARC WITH SAME GRPID'
                        arcfile= arclist[setup][_grpid]     # flat in the  raw data
                    else:  
                        arcfile=[]
                        for _grpid0 in arclist[setup].keys():
                            for ii in arclist[setup][_grpid0]:
                                arcfile.append(ii)                   
              if arcfile:
                  if len(arcfile)>1:                           # more than one arc available
                      print arcfile
#                     _arcclose=floyds.util.searcharc(imgex,arcfile)[0]   # take the closest in time 
                      _arcclose=floyds.sortbyJD(arcfile)[-1]               #  take the last arc of the sequence
                      if _interactive.upper() in ['YES','Y']:
                              for ii in floyds.floydsspecdef.sortbyJD(arcfile):
                                  print '\n### ',ii 
                              arcfile=raw_input('\n### more than one arcfile available, which one to use ['+str(_arcclose)+'] ? ')
                              if not arcfile: arcfile=_arcclose
                      else: arcfile=_arcclose
                  else: arcfile=arcfile[0]
              else:   print '\n### Warning: no arc found'

###################################################################   rectify 
              if setup[0]=='red':
                  fcfile=floyds.__path__[0]+'/standard/ident/fcrectify_'+_tel+'_red'
                  fcfile1=floyds.__path__[0]+'/standard/ident/fcrectify1_'+_tel+'_red'
                  print fcfile
              else:
                  fcfile=floyds.__path__[0]+'/standard/ident/fcrectify_'+_tel+'_blue'
                  fcfile1=floyds.__path__[0]+'/standard/ident/fcrectify1_'+_tel+'_blue'
                  print fcfile
              print img,arcfile,flatfile
              img0=img
              if img      and not img in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(img)
              if arcfile  and arcfile not in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(arcfile)
              if flatfile and flatfile not in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(flatfile)

              img,arcfile,flatfile=floyds.floydsspecdef.rectifyspectrum(img,arcfile,flatfile,fcfile,fcfile1,'no',_cosmic)
              print outputfile
              if img      and not img in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(img)
              if arcfile  and arcfile not in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(arcfile)
              if flatfile and flatfile not in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(flatfile)

###################################################################         check wavecalib  
              if tpe=='std' or floyds.util.readkey3(floyds.util.readhdr(img),'exptime') < 300:
                  if setup[0]=='red':
                      print '\n### check standard wave calib'
                      #print img
                      data, hdr = pyfits.getdata(img, 0, header=True) 
                      y=data.mean(1)
                      import numpy as np
                      if np.argmax(y) < 80 and np.argmax(y) > 15:                      
                          y2=data[np.argmax(y)-3:np.argmax(y)+3].mean(0)
                          yy2=data[np.argmax(y)-9:np.argmax(y)-3].mean(0)
                          floyds.util.delete('_std.fits')
                          pyfits.writeto('_std.fits', np.float32(y2-yy2), hdr)
                          #print '_std.fits',_interactive
                          shift=floyds.floydsspecdef.checkwavestd('_std.fits',_interactive,2)
                          zro=hdr['CRVAL1']
                          floyds.util.updateheader(img,0,{'CRVAL1':[zro+int(shift),'']})
                          floyds.util.updateheader(img,0,{'shift':[float(shift),'']})
                          floyds.util.delete('_std.fits')
                      else:
                          print 'object not found'
                  else: 
                      print '\n### warning check in wavelength not possible for short exposure in the blu range '
              else:
                    print '\n### check object wave calib'
                    _skyfile=floyds.__path__[0]+'/standard/ident/sky_'+setup[0]+'.fits'
                    data, hdr = pyfits.getdata(img, 0, header=True) 
                    y=data.mean(1)
                    import numpy as np
                    if np.argmax(y) < 80 and np.argmax(y) > 15:
                        yy1=data[10:np.argmax(y)-9].mean(0)
                        yy2=data[np.argmax(y)+9:-10].mean(0)
                        floyds.util.delete('_sky.fits')
                        pyfits.writeto('_sky.fits', np.float32(yy1+yy2), hdr)
                        shift=floyds.floydsspecdef.checkwavelength_obj('_sky.fits',_skyfile,_interactive,2)
                        floyds.util.delete('_sky.fits')
                        zro=hdr['CRVAL1']
                        floyds.util.updateheader(img,0,{'CRVAL1':[zro+int(shift),'']})
                        floyds.util.updateheader(img,0,{'shift':[float(shift),'']})
                    else:  print 'object not found'
####################################################     flat field
              if img and flatfile and setup[0]=='red':
                      imgn='n'+img
                      hdr1 = floyds.readhdr(img)
                      hdr2 = floyds.readhdr(flatfile)
                      _grpid1=floyds.util.readkey3(hdr1,'grpid')
                      _grpid2=floyds.util.readkey3(hdr2,'grpid')
                      if _grpid1==_grpid2:
                          imgn=floyds.fringing_classicmethod2(flatfile,img,'no','*',15,setup[0])
                      else:
                          print 'Warning flat not the same OB'
                          imgex=floyds.floydsspecdef.extractspectrum(img,dv,_ext_trace,_dispersionline,_interactive,tpe,automaticex=_automaticex)
                          floyds.delete('flat'+imgex)
                          iraf.specred.apsum(flatfile,output='flat'+imgex,referen=img,interac='no',find='no',recente='no',resize='no',\
                                             edit='no',trace='no',fittrac='no',extract='yes',extras='no',review='no',backgro='none')
                          fringingmask=floyds.normflat('flat'+imgex)
                          print '\n### fringing correction'
                          print imgex,fringingmask
                          imgex,scale,shift=floyds.correctfringing_auto(imgex,fringingmask)  #  automatic correction
                          shift=int(.5+float(shift)/3.5)        # shift from correctfringing_auto in Angstrom
                          print '\n##### flat scaling: ',str(scale),str(shift)
########################################################
                          datax, hdrx = pyfits.getdata(flatfile, 0, header=True)
                          xdim=hdrx['NAXIS1']
                          ydim=hdrx['NAXIS2']
                          iraf.specred.apedit.nsum=15 
                          iraf.specred.apedit.width=100.  
                          iraf.specred.apedit.line=1024 
                          iraf.specred.apfind.minsep=20.  
                          iraf.specred.apfind.maxsep=1000.  
                          iraf.specred.apresize.bkg='no' 
                          iraf.specred.apresize.ylevel=0.5 
                          iraf.specred.aptrace.nsum=10
                          iraf.specred.aptrace.step=10
                          iraf.specred.aptrace.nlost=10
                          floyds.util.delete('n'+flatfile)
                          floyds.util.delete('norm.fits')
                          floyds.util.delete('n'+img)
                          floyds.util.delete(re.sub('.fits','cut.fits',flatfile))
                          iraf.imcopy(flatfile+'[500:'+str(xdim)+',*]',re.sub('.fits','cut.fits',flatfile),verbose='no')
                          iraf.imarith(flatfile,'/',flatfile,'norm.fits',verbose='no')
                          flatfile=re.sub('.fits','cut.fits',flatfile)
                          floyds.util.delete('n'+flatfile)
                          iraf.unlearn(iraf.specred.apflatten)
                          floyds.floydsspecdef.aperture(flatfile)
                          iraf.specred.apflatten(flatfile,output='n'+flatfile,interac=_interactive,find='no',recenter='no', resize='no',edit='no',trace='no',\
                                                 fittrac='no',fitspec='no', flatten='yes', aperture='',\
                                                 pfit='fit2d',clean='no',function='legendre',order=15,sample = '*', mode='ql')
                          iraf.imcopy('n'+flatfile,'norm.fits[500:'+str(xdim)+',*]',verbose='no')
                          floyds.util.delete('n'+flatfile)
                          floyds.util.delete('n'+img)
                          iraf.imrename('norm.fits','n'+flatfile,verbose='no')
                          imgn=floyds.floydsspecdef.applyflat(img,'n'+flatfile,'n'+img,scale,shift)
              else:                  imgn=''

              if imgn    and not imgn in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(imgn)

###################################################      2D flux calib
              print '####### '+imgn
              hdr=floyds.util.readhdr(img)
              _sens=''
              if liststandard:  _sens=floyds.util.searchsens(img,liststandard)[0]   # search in the list from reducer
              if not _sens:
                  try:      _sens=floyds.util.searchsens(img,sens[setup])[0]        # search in the reduced data
                  except:   _sens=floyds.util.searchsens(img,'')[0]              # search in tha archive
              if _sens:
                  if _sens[0]=='/': 
                      os.system('cp '+_sens+' .')
                      _sens=string.split(_sens,'/')[-1]
                  imgd=fluxcalib2d(img,_sens)
                  if imgn:     imgdn=fluxcalib2d(imgn,_sens)
                  else: imgdn=''
                  if _sens not in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(_sens)
                  else:        imgdn='' 
                  print '\n### do 2D calibration'
              else:
                  imgd=''
                  imgdn=''
################    extraction         ####################################
              if imgdn:
                  try:
                      imgdnex=floyds.floydsspecdef.extractspectrum(imgdn,dv,_ext_trace,_dispersionline,_interactive,tpe,automaticex=_automaticex)
                  except:
                      imgdnex=''
              else:       
                  imgdnex=''
              if imgd:
                  try:
                      imgdex=floyds.floydsspecdef.extractspectrum(imgd,dv,_ext_trace,_dispersionline,_interactive,tpe,automaticex=_automaticex)  
                  except:
                      imgdex=''
              else:
                  imgdex=''
              if imgd    and not imgd in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(imgd)
              if imgdn   and not imgdn in outputfile[tpe][archfile]: outputfile[tpe][archfile].append(imgdn)
              if imgdnex and imgdnex not in outputfile[tpe][archfile]:   outputfile[tpe][archfile].append(imgdnex)
              if imgdex  and imgdex not in outputfile[tpe][archfile]:    outputfile[tpe][archfile].append(imgdex)
              if tpe=='std':
                  if imgn:
                      try:
                          imgnex=floyds.floydsspecdef.extractspectrum(imgn,dv,_ext_trace,_dispersionline,_interactive,tpe,automaticex=_automaticex)  
                      except:
                          imgnex=''
                  elif img:
                      try:
                          imgnex=floyds.floydsspecdef.extractspectrum(img,dv,_ext_trace,_dispersionline,_interactive,tpe,automaticex=_automaticex)  
                      except:
                          imgnex=''
                  if imgnex:
                    hdrs=floyds.util.readhdr(imgnex)
                    _tel=floyds.util.readkey3(hdrs,'TELID')
                    try:
                      _outputsens2='sens_'+_tel+'_'+str(floyds.util.readkey3(hdrs,'date-night'))+'_'+str(floyds.util.readkey3(hdrs,'grism'))+\
                          '_'+re.sub('.dat','',floyds.util.readkey3(hdrs,'stdname'))+'_'+str(MJDtoday)
                    except:  sys.exit('Error: missing header -stdname- in standard '+str(standardfile)+'  ')                          
                    print '\n### compute sensitivity function and atmofile'
                    if setup[0]=='red':
                          atmofile=floyds.floydsspecdef.telluric_atmo(imgnex)
                          #print atmofile
                          if atmofile  and atmofile not in outputfile[tpe][archfile]:    outputfile[tpe][archfile].append(atmofile)
                          stdusedclean=re.sub('_ex.fits','_clean.fits',imgnex)
                          floyds.util.delete(stdusedclean)
                          _function='spline3'
                          iraf.specred.sarith(input1=imgnex,op='/',input2=atmofile,output=stdusedclean, format='multispec')
                          _outputsens2=floyds.floydsspecdef.sensfunction(stdusedclean,_outputsens2,_function,8,_interactive)
                          if setup not in atmo: atmo[setup]=[atmofile]
                          else: atmo[setup].append(atmofile)
                    else:
                          _function='spline3'
                          _outputsens2=floyds.floydsspecdef.sensfunction(imgnex,_outputsens2,_function,12,_interactive,'3400:4700')#,3600:4300')
                    if _outputsens2  and _outputsens2 not in outputfile[tpe][archfile]:    outputfile[tpe][archfile].append(_outputsens2)
    ###################################################
    print outputfile
    if 'obj' in outputfile:
      for imm in outputfile['obj']:
        lista= outputfile['obj'][imm]
        lista1=[]
        for i in lista: 
            if '_ex.fits' in i:  
                lista1.append(i)
                lista.pop(lista.index(i))
        done=[]
        for img in lista1:
            _output=''
            if img not in done:
                if '_red_' in img:
                    redfile=img
                    if re.sub('_red_','_blue_',redfile)[1:] in lista1:
                        bluefile=lista1[lista1.index(re.sub('_red_','_blue_',redfile)[1:])]
                        _output=re.sub('_red_','_merge_',redfile)
                        try:
                            _output=floyds.floydsspecdef.combspec(bluefile,redfile,_output,scale=True,num=None)
                            lista.append(_output)
                            done.append(redfile)
                            done.append(bluefile)
                            floyds.util.delete(bluefile)
                            floyds.util.delete(redfile)
                            floyds.util.delete(redfile[1:])
                        except:
                            done.append(redfile)
                            lista.append(redfile)
                    else:
                        done.append(redfile)
                        lista.append(redfile)
                elif '_blue_' in img:
                    bluefile=img
                    if 'n'+re.sub('_blue_','_red_',bluefile) in lista1:
                        redfile=lista1[lista1.index('n'+re.sub('_blue_','_red_',bluefile))]
                        _output=re.sub('_red_','_merge_',redfile)
                        try:
                            _output=floyds.floydsspecdef.combspec(bluefile,redfile,_output,scale=True,num=None)
                            lista.append(_output)
                            done.append(redfile)
                            done.append(bluefile)
                            floyds.util.delete(bluefile)
                            floyds.util.delete(redfile)
                            floyds.util.delete(redfile[1:])
                        except:
                            done.append(bluefile)
                            lista.append(bluefile)
                    else:
                        done.append(bluefile)
                        lista.append(bluefile)
        outputfile['obj'][imm]= lista
    readme=floyds.floydsspecauto.writereadme()
    return outputfile,readme

######################################################
def badimage(img,_type):
    from pyraf import iraf
    import string
    import floyds
    hdr=floyds.util.readhdr(img)
    _tel=floyds.util.readkey3(hdr,'telescop')
    rr1=float(string.split(iraf.imstat(img+'[1250:1270,140:200]',Stdout=1)[1])[2])
    rr2=float(string.split(iraf.imstat(img+'[1570:1620,170:210]',Stdout=1)[1])[2])
    bb1=float(string.split(iraf.imstat(img+'[1090:1110,370:400]',Stdout=1)[1])[2])
    bb2=float(string.split(iraf.imstat(img+'[1560:1700,340:360]',Stdout=1)[1])[2])
    bg1=float(string.split(iraf.imstat(img+'[1500:1540,290:310]',Stdout=1)[1])[2])
    #print rr1,rr2,bb1,bb2,bg1
    good=1
    if _type in ['flat','arc','lamp']:
        if rr1/rr2<=1.25 and bb1/bb2<=1.25:   
            good=0
    if good==1:
        if _type in ['flat']:
            if rr1/bg1<20 and bb1/bg1>2.5 and 'South' in _tel:
                good=0
            if bb1/bg1<1.1:
                good=0
            if rr1<1000:
                good=0
    if good==0:  print '\n### warning:  image '+img+' rejected'
    return good

#########################################
def writereadme():
    readme='#######################\n#\n#\n#  floyds automatic pipeline \n#\n#######################\n\n\n'+\
        'The files included in this tar have been reduced automatically\n'+\
        'The main steps of the automatic reduction include:\n'+\
        '- rectification of the frames: science,flat,arc along y axes\n'+\
        '- rectification of the frames: science,flat,arc along x axes\n'+\
        '- fringing correction (only for the red part of the spectrum)\n'+\
        '- wavelength check using sky line or telluric lines\n'+\
        '- flux calibration using sensitivity average function\n'+\
        '- fast extraction\n'+\
        '- if the science file is a standard star, sensitivity function and atmo file are also computed\n\n#####################\n'+\
        'The user should use the 2D wavelength and flux calibrated files (extension _2df.fits) and optimize the extraction for his/her science\n\n'+\
        'the file Ntt*****red***.2df.fits has been corrected for fringing using the normalized flat field\n'+\
        'fringing file, shift applied to the original flatfield, scale applied to the original flatfield are specified in the header \n'+\
        'the file tt*****red***.2df.fits has NOT been corrected for fringing\n'+\
        'The fast extracted spectra have extension _2df_ex.fits\n'+\
        'files starting with tt have been rectified along x and y axes\n\n##############\n\n        GO FLOYDS !!!!!!!\n'
    f=open('README','w')
    f.write(readme)
    f.close()
    return 'README'

#####################################################
def archivespectrum(img,_force=True):
    import pyfits
    import datetime
    import floyds
    import string,re,os,glob
    hdr=pyfits.open(img)[0].header
    _tel=hdr['telescop']
    if 'South' in _tel: _tel='fts'
    if 'North' in _tel: _tel='ftn'
    _object=hdr['OBJECT']
    _date=hdr['DATE-OBS']
    _grism=hdr['grism']
    user=os.environ['USER']
    try:        _UT=hdr['UTSTART']
    except:     _UT=hdr['UT']
    a=(datetime.datetime.strptime(string.split(_date[2:],'.')[0],"%y-%m-%dT%H:%M:%S")-datetime.timedelta(.0)).isoformat()
    a=re.sub('-','',string.split(a,'T')[0])
    directory='/science/'+str(user)+'/data/WEB/floyds/'+a+'_'+_tel
    try:
        if os.path.isdir(directory): print 'directory there'
        else:                        os.system('mkdir '+directory)
        imglist=glob.glob(directory+'/*_'+_tel+'_*2df*fits')
        filethere=0
        for imgold in imglist:
            hdr1=pyfits.open(imgold)[0].header
            _date1=hdr1['DATE-OBS']
            _grism1=hdr1['grism']
            if _date1==_date and _grism1==_grism:
                filethere=1
                break
        if filethere:
            if _force:
                os.system('rm '+imgold)      
                os.system('cp '+img+' '+directory)
        else:  os.system('cp '+img+' '+directory)
    except:  pass
#########################################################

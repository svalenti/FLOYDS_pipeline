#!/bin/tcsh

set prog = 'Process_Night_FLOYDS'
if ( $#argv != 2 && $#argv != 3 ) then
  echo "$prog FATAL: Wrong number of command line arguments"
  echo "Usage: $prog <date> <camera> [0|1]"
  exit(-1)
endif
set UTnight = $argv[1]
set camera = $argv[2]
set pushtoarchive = 1
if ( $#argv == 3 ) then
  if ( $argv[3] == 0 ) then
      set pushtoarchive = 0
  endif
endif  
echo "$prog DBG: Setting processing night,camera to $UTnight,$camera with archive push=$pushtoarchive from cmdline."

set bodge = 0
set oldstyle = 0 # Set to 1 for old-style (RCS-based) data and data locations
set finalroot = '/data/archive/data'
set pipepath =  '/home/data/FLOYDS_pipeline'
set tel = ''
set camprefix = $camera
if ( $camera == 'en06' ) then
  set site = 'ogg'
  set camprefix = $site
  if ( $oldstyle == 1 ) then
    set tel = 'ft1'
    set camprefix = 'f_'
  endif
else if ( $camera == 'en05' ) then
  set site = 'coj'
  set camprefix = $site
  if ( $oldstyle == 1 ) then
    set tel = 'ft2'
    set camprefix = 'g_'
  endif
else
  echo "$prog FATAL: Unrecognized camera (${camera})"
  exit (-1)
endif

set rawroot = /data/archive/data
set outpath = /data/archive/floyds/$site/proc
# Testing
set plumber = 'tlister@lcogt.net'
# Production
#set plumber = 'tlister@lcogt.net,svalenti@lcogt.net'

# Setup
if ( $oldstyle == 0 ) then
  set rawroot = /mfs-sba/engineering/
  setenv FLOYDS_DATA_IN $rawroot/$site/$camera/$UTnight/raw
else
  setenv FLOYDS_DATA_IN $rawroot/$tel/Raw/$UTnight
endif
setenv FLOYDS_DATA_OUT $outpath/$UTnight
setenv FINAL_DATA_OUT /data/archive/webfiles/new/
set LOGDIR = ~/Pipeline_Logs
set FITS_GET_HEADER = "/usr/local/bin/wcs/imhead"

# Check raw data directory exists
if ( ! -rd  $FLOYDS_DATA_IN ) then
  echo "Specified night of raw data (${UTnight}) does not exist."
  echo "(Looked in $FLOYDS_DATA_IN)"
  exit (-1)
endif

echo "FLOYDS_DATA_IN= $FLOYDS_DATA_IN"
echo "FLOYDS_DATA_OUT= $FLOYDS_DATA_OUT"
echo "FINAL_DATA_OUT= $FINAL_DATA_OUT"
#exit(-42)

# Make directories
if ( ! -rd $FLOYDS_DATA_OUT ) then
  mkdir -p $FLOYDS_DATA_OUT
  set pipestatus = $status
  if ( $status != 0 ) then
    echo "$prog FATAL: Unable to make directory $FLOYDS_DATA_OUT"
    exit (-3)
  endif
endif
if ( ! -rd $FINAL_DATA_OUT ) then
  mkdir -p $FINAL_DATA_OUT
  set pipestatus = $status
  if ( $status != 0 ) then
    echo "$prog FATAL: Unable to make directory $FINAL_DATA_OUT"
    exit (-3)
  endif
endif
if ( ! -rd $LOGDIR ) then
  mkdir -p $LOGDIR
  set pipestatus = $status
  if ( $status != 0 ) then
    echo "$prog FATAL: Unable to make directory $LOGDIR"
    exit (-3)
  endif
endif

# Copy over data
set rsync_options = '-tOrL --timeout=90'
set rsync_excludes = ''
rsync ${rsync_options} $FLOYDS_DATA_IN/${camprefix}*  $FLOYDS_DATA_OUT/ ${rsync_excludes}
#exit(-42)

cd $FLOYDS_DATA_OUT
# Process data
if ( ! -r ${UTnight}.processed ) then
  if ( `ls -1 $FLOYDS_DATA_OUT/ | grep '.fits' | wc -l | awk '{print $1}'` > 0 ) then
    if ( `ls -1 $FLOYDS_DATA_OUT/ | grep '.fits.gz' | wc -l | awk '{print $1}'` > 0 ) then
      gunzip -f $FLOYDS_DATA_OUT/${camprefix}*.fits.gz
    endif
    $pipepath/bin/floydsauto -X >>& $LOGDIR/Pipeline_${UTnight}_${camera}.log
    set floydsstatus = $status
    echo "Pipeline Status was=$floydsstatus"
    if ( $floydsstatus != 0 ) then
      tail -n 20 $LOGDIR/Pipeline_${UTnight}_${camera}.log  >! /tmp/mail.message
      mail -s "Error running FLOYDS pipeline for $UTnight $camera" $plumber < /tmp/mail.message
      \rm -f /tmp/mail.message
      exit(-4)
    endif
    # move _1.fits and tarballs to sci-archive
    if ( `ls -1 $FLOYDS_DATA_OUT/ | grep '90.fits' | wc -l | awk '{print $1}'` > 0 ) then
      foreach donefile ( *90.fits )
	set donefilename = `basename $donefile`
	set donefileprefix = `echo $donefilename | cut -d. -f1 | cut -d- -f-6`
        set propid = `fits_get_keyword_value_static $donefile PROPID STRING`
        if ( ! -rd $FINAL_DATA_OUT/$propid ) then
          mkdir -p $FINAL_DATA_OUT/$propid
        endif
        set pubpriv = `fits_get_keyword_value_static $donefile L1PUBPRV STRING`
        if ( $pubpriv != 'public' ) then
	    echo "AuthUserFile /data/archive/webfiles/.htpasswd" > ${FINAL_DATA_OUT}/"${propid}"/.htaccess
	    echo "AuthGroupFile /dev/null" >> ${FINAL_DATA_OUT}/"${propid}"/.htaccess
	    echo "AuthName ${propid}" >> ${FINAL_DATA_OUT}/"${propid}"/.htaccess
	    echo "AuthType Basic" >> ${FINAL_DATA_OUT}/"${propid}"/.htaccess
	    echo "<limit GET>" >> ${FINAL_DATA_OUT}/"${propid}"/.htaccess
	    echo "require user oraceng ${propid}" >> ${FINAL_DATA_OUT}/"${propid}"/.htaccess
	    echo "</Limit>" >> ${FINAL_DATA_OUT}/"${propid}"/.htaccess
        else
            echo "Public data, not creating .htaccess file."
        endif
        if ( ! -rd $FINAL_DATA_OUT/$propid/$UTnight ) then
          mkdir -p $FINAL_DATA_OUT/$propid/$UTnight
        endif
        cp -v $donefile $FINAL_DATA_OUT/$propid/$UTnight/
        set tarfile = `fits_get_keyword_value_static $donefile TARFILE STRING`
        mv -v $tarfile $FINAL_DATA_OUT/$propid/$UTnight/
        \rm -f $donefile
        
	set fits_in = "${FINAL_DATA_OUT}/${propid}/${UTnight}/${donefilename}"
	set txt_out = "${FINAL_DATA_OUT}/${propid}/${UTnight}/${donefileprefix}.head"
	if ( ! -e $fits_in ) then
	    echo "     $fits_in cannot be found for some reason; skipping this one"
	    continue
	else
	    $FITS_GET_HEADER $fits_in > $txt_out  
	    if ( $status != 0 ) then
	        echo "     fits_get_header did not work. Continuing anyway."
	    endif
	endif
        
      end
      /usr/local/bin/recent_webindex
    else
      echo "No processed frames produced"
    endif
  else
    echo "No actual object frames to do."
  endif
else
  echo "$UTnight already processed ('${UTnight}.processed' flagfile exists)"
  set floydsstatus = 0
endif

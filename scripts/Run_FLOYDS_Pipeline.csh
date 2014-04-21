#!/bin/tcsh

set oprog = 'Run_FLOYDS_Pipeline'
set pushtoarchive = 1
if ( $#argv == 1 ) then
  set camera = $argv[1]
  set UTnight = `date -u --date="-1 days" +%Y%m%d`
  if ( $camera == 'en05' ) then
    set UTnight = `date -u +%Y%m%d`
  endif
else if ( $#argv == 2 ) then
  # Test if the first two characters of $argv[1] are digits or not.
  # If they are digits, then we got a date and a camera otherwise we got a
  # camera and a push instruction
  echo $argv[1] | cut -c1-2 | egrep '^[0-9][0-9]*$' >&! /dev/null
  if ( $status != 0 ) then
    echo "First two characters not digits, assuming camera+push instruction"
    set UTnight = `date -u --date="-1 days" +%Y%m%d`
    set camera = $argv[1]
    if ( $camera == 'en05' ) then
      set UTnight = `date -u +%Y%m%d`
    endif
    if ( $argv[2] == 0 ) then
	set pushtoarchive = 0
    endif
  else
    set UTnight = $argv[1]
    set camera = $argv[2]
  endif
else if ( $#argv == 3 ) then
  set UTnight = $argv[1]
  set camera = $argv[2]
  if ( $argv[3] == 0 ) then
      set pushtoarchive = 0
  endif
else  
  echo "$oprog FATAL: Wrong number of command line arguments"
  echo "Usage: $oprog [<date>] <camera> [push to archive 0|1]"
  echo "(if [<date>] is omitted, the current processing date is used)"
  exit(-1)
endif

set logfile = ~/Pipeline_Logs/Pipeline_${UTnight}_${camera}.log
touch $logfile
echo "$oprog INFO: Calling Process_Night_FLOYDS.csh to process to $UTnight,$camera with archive push=$pushtoarchive from cmdline." >> $logfile
source ~/Process_Night_FLOYDS.csh $UTnight $camera $pushtoarchive>> $logfile
set pipestate = $status

echo "$oprog INFO: Pipeline status was $pipestate" >> $logfile

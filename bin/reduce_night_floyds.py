"""
Reduce a night's worth of FLOYDS data

Author
    Curtis McCully (cmccully@lco.global)

January 2016
"""
import argparse
import datetime
import os
import shutil
from glob import glob

import smtplib
from email import mime
from astropy.io import fits


def main():
    """
    Main driver script to reduce the night of FLOYDS data
    :return: None
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--dayobs', default=None, dest='dayobs',
                        help='Dayobs of data to reduce. Default is last night. Example: 20160101')
    parser.add_argument('--site', dest='site', choices=['ogg', 'coj'],
                        help='Telescope site to reduce.')
    parser.add_argument('--maintainer-email', dest='maintainer',
                        help='Email of pipeline maintainer.')
    parser.add_argument('--sender-email', dest='sender',
                        help='Email account to send crash reports as. Must be a gmail account')
    parser.add_argument('--sender-password', dest='sender_password',
                        help='Password for sender email.')
    args = parser.parse_args()

    # Get the UT date from either the command line. Default to last night
    if args.dayobs is None:
        last_night = datetime.datetime.now()
        if args.site == 'ogg':
            last_night -= datetime.timedelta(days=1)
        args.dayobs = last_night.strftime('%Y%m%d')

    # Figure out which camera we are reducing
    if args.site == 'coj':
        args.camera = 'en05'
    else:
        args.camera = 'en06'

    # Copy the raw files into a temporary directory here
    working_path = '{site}/{dayobs}'.format(args)
    if not os.path.exists(working_path):
        os.makedirs(working_path)

    os.chdir(working_path)

    raw_files = glob('/archive/engineering/{site}/{camera}/{dayobs}/raw/*.fits.fz'.format(args))

    if len(raw_files) > 0:
        # Copy and unpack the raw data
        for f in raw_files:
            shutil.copy(f, '.')
            uncompressed_file = os.path.splitext(os.path.basename(f))
            os.system('funpack -O {output} {input}'.format(input=f, output=uncompressed_file))

        gather_guider_frames()

        # reduce the data
        success = os.system('floydsauto -X')
        # If there is an issue, send an email alert
        if success > 0:
            send_error_report(args)
        else:
            # Copy the tarballs to the output directory
            archive_tar_files()


def send_error_report(args):
    msg = mime.multipart()
    msg['From'] = args.sender
    msg['To'] = args.maintainer
    msg['Subject'] = "Error Running FLOYDS Pipeline for {dayobs} at {site}".format(args)

    server = smtplib.SMTP('smtp.gmail.com', 587)
    server.starttls()
    server.login(args.sender, args.sender_password)

    msg = "The FLOYDS pipeline exited with an exit code > 0"
    server.sendmail(args.sender, args.maintainer, msg.as_string())
    server.quit()


def archive_tar_files():
    # Find all of the tar files and put them in the specproc directory and post them to the archive
    
    pass


def gather_guider_frames():
    # Find all of the corresponding guider frames to put into a multi-extension fits file.
    # Get the unique block ids for each of the science frames (e00)
    raw_science_files = glob('*e00.fits')
    block_meta_data = [{'blockid': fits.getval(f, 'BLKUID'), 'agcamera': fits.getval(f, 'AGCAM'),
                         'site': fits.getval(f, 'SITEID'), 'dayobs': fits.getval(f, 'DAY-OBS')}
                        for f in raw_science_files]

    # Filter out the duplicates
    block_meta_data = list({d['blockid']:d for d in block_meta_data}.values())

    for observation_block in block_meta_data:
        guider_frames = glob('/archive/engineering/{site}/{camera}/{dayobs}/raw/*.fz'
                             .format(site=observation_block['site'],
                                     camera=observation_block['agcamera'],
                                     dayobs=observation_block['dayobs']))
        for frame in guider_frames:
            uncompressed_file = os.path.splitext(os.path.basename(f))
            os.system('funpack -O {output} {input}'.format(input=f, output=uncompressed_file))
            # Run SEP to get the FWHM of the frame
            # Use astrometry.net to solve for the WCS
        # Put all of the frames into a single multi-extension fits file in the specproc directory
        # Post it to the archive
        # Make the summary plot

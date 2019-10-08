#!/usr/bin/python
# -*- coding: utf-8 -*-

from dateutil.parser import parse
import os

from pytplot import cdf_to_tplot
from .orbit_time import orbit_time
from .download_files_utilities import *
from .read_kp import *
from .kp_to_tplot import *


def maven_filenames(instruments=None,
                    start_date='2014-01-01',
                    end_date='2020-01-01',
                    update_prefs=False,
                    only_update_prefs=False,
                    local_dir=None,
                    public=True):
    """
    This function identifies which MAVEN data to download.
    """

    # Check for orbit num rather than time string
    if isinstance(start_date, int) and isinstance(end_date, int):
        start_date, end_date = orbit_time(start_date, end_date)
        start_date = parse(start_date)
        end_date = parse(end_date)
        start_date = start_date.replace(hour=0, minute=0, second=0)
        end_date = end_date.replace(day=end_date.day+1, hour=0, minute=0, second=0)
        start_date = start_date.strftime('%Y-%m-%d')
        end_date = end_date.strftime('%Y-%m-%d')
        
    if update_prefs or only_update_prefs:
        set_new_data_root_dir()
        if only_update_prefs:
            return

    if not public:
        get_uname_and_password()

    # No instruments indicated --> error
    if instruments is None:
        print("You must specify at least one instrument from which you want data downloaded.")
        return

    # Default data level to download for each instrument is L2, but if KP data is downloaded,
    # it will be reset for that "instrument" below
    level = ['l2']*len(instruments)

    # If downloading KP dada, figure out if you're downloading insitu or iuvs (but not both)
    if all(x in instruments for x in ['kp-insitu', 'kp-iuvs']):
        print("Can't request both INSITU and IUVS in one query.")
        return
    elif 'kp-insitu' in instruments and 'kp-iuvs' not in instruments:
        kp_index = instruments.index('kp-insitu')
        instruments[kp_index] = 'kp'
        level[kp_index] = 'insitu'
    elif 'kp-iuvs' in instruments and 'kp-insitu' not in instruments:
        kp_index = instruments.index('kp-iuvs')
        instruments[kp_index] = 'kp'
        level[kp_index] = 'iuvs'

    # Set data download location
    if local_dir is None:
        mvn_root_data_dir = get_root_data_dir()
    else:
        mvn_root_data_dir = local_dir

    # Keep track of files to download
    maven_files = {}
    for i, instrument in enumerate(instruments):
        # Build the query to the website
        query_args = []
        query_args.append("instrument=" + instrument)
        query_args.append("level=" + level[i])
        query_args.append("start_date=" + start_date)
        query_args.append("end_date=" + end_date)
        if level == 'iuvs':
            query_args.append("file_extension=tab")

        data_dir = os.path.join(mvn_root_data_dir, 'maven', 'data', 'sci', instrument, level[i])
        
        query = '&'.join(query_args)
        
        s = get_filenames(query, public)

        if not s:
            print("No files found for {}.".format(instrument))
            maven_files[instrument] = []
            continue

        s = s.split(',')

        maven_files[instrument] = [s, data_dir, public, level[i]]

    return maven_files


def load_data(instruments=None,
              kp_instruments=None,
              start_date='2014-01-01',
              end_date='2020-01-01',
              update_prefs=False,
              only_update_prefs=False,
              local_dir=None,
              list_files=False,
              new_files=False,
              exclude_orbit_file=False,
              download_only=False,
              varformat=None,
              prefix='',
              suffix='',
              get_support_data=False,
              public=True):
    """
    This function downloads MAVEN data loads it into tplot variables, if applicable.
    """

    # 1. Download files

    maven_files = maven_filenames(instruments, start_date, end_date, update_prefs, only_update_prefs, local_dir, public)

    # Keep track of what files are downloaded
    downloaded_files = []

    # Track if KP files were downloaded or not
    kp_downloaded = False

    # Actually download files here
    for instr in maven_files.keys():
        if maven_files[instr]:
            s = maven_files[instr][0]
            data_dir = maven_files[instr][1]
            public = maven_files[instr][2]
            level = maven_files[instr][3]
            if list_files:
                for f in s:
                    print(f)
                return

            print("Your request will download a total of: "+str(len(s))+" files for instrument "+str(instr))
            print('Would you like to proceed with the download? ')
            valid_response = False
            cancel = False
            while not valid_response:
                response = (input('(y/n) >  '))
                if response == 'y' or response == 'Y':
                    valid_response = True
                    cancel = False
                elif response == 'n' or response == 'N':
                    print('Cancelled download. Returning...')
                    valid_response = True
                    cancel = True
                else:
                    print('Invalid input.  Please answer with y or n.')

            if cancel:
                continue

            if not exclude_orbit_file:
                print("Before downloading data files, checking for updated orbit # file from naif.jpl.nasa.gov")
                print("")
                get_orbit_files()

            # If we want to only download new files (i.e., files not already on the user's machine),
            # make of a note of what requested downloads already are on the user's machine.
            if new_files:
                files_already_on_hd = get_new_files(s, data_dir, instr, level)

            i = 0
            display_progress(i, len(s))
            for f in s:
                i = i+1
                full_path = create_dir_if_needed(f, data_dir, level)
                if new_files:
                    get_file_from_site(f, public, full_path, files_on_hd=files_already_on_hd)
                else:
                    get_file_from_site(f, public, full_path)
                display_progress(i, len(s))

                downloaded_files.append(os.path.join(full_path, f))

            if instr == 'kp':
                kp_downloaded = True

    # 2. Load files into tplot

    if downloaded_files:
        # Flatten out downloaded files from list of lists of filenames
        if isinstance(downloaded_files[0], list):
            downloaded_files = [item for sublist in downloaded_files for item in sublist]

        # If CDF files downloaded, grab them to load into tplot
        cdf_files = [f for f in downloaded_files if '.cdf' in f]

        # If KP files downloaded, grab them to load into tplot
        kp_files = [f for f in downloaded_files if '.tab' in f]

        if not download_only:
            stored_vars = []
            # Create CDF tplot variables
            stored_vars.extend(cdf_to_tplot(cdf_files, varformat=varformat, get_support_data=get_support_data,
                                            prefix=prefix, suffix=suffix, merge=True))
            # Create KP in situ tplot vars
            # Currently, no functionality to put KP iuvs data into tplot vars,
            # but could use the read_kp() function elsewhere to load KP iuvs data
            # into a list object
            if kp_downloaded:
                kp_read = read(kp_files, kp_instruments)
                if 'insitu' in kp_files[0]:
                    stored_kp_vars = kp_to_tplot(kp_read)
                    stored_vars.extend(stored_kp_vars)

            return stored_vars


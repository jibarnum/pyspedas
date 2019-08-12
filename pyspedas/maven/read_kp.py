# -*- coding: utf-8 -*-
"""
File:
    read_kp.py

Description:
    Functions used by maven_load to read KP files and transform them into a python dict of dicts/list of dicts.
"""

import calendar
import numpy as np
from datetime import datetime

from .kp_utilities import param_dict
from .kp_utilities import remove_inst_tag
from .kp_utilities import read_iuvs_file
from .kp_utilities import get_header_info
from _collections import OrderedDict
import builtins


def read(filenames=None,
         instruments=None):
    '''
    Read MAVEN Key Parameter (KP) data into a dictionary object.

    Input:
        filenames: KP files to be read in
        instruments: Optional keyword listing the instruments to include
        in the returned dictionary/structure.

    Output:
        A dictionary containing up to all of the columns
        included in a MAVEN KP in situ data file. Or, a list of dictionaries
        in the case of a MAVEN KP iuvs data file.
    '''
    import pandas as pd
    import re

    if instruments is not None:
        if not isinstance(instruments, builtins.list):
            instruments = [instruments]

    if 'insitu' in filenames[0]:
        # Loads in KP in situ data

        # Get column names from first file
        names, inst = get_header_info(filenames[0])

        # Strip off the first name for now (Time), and use that as the dataframe index.
        names = names[1:len(names)]
        inst = inst[1:len(inst)]

        # Break up dictionary into instrument groups
        lpw_group, euv_group, swe_group, swi_group, sta_group, sep_group, mag_group, ngi_group, app_group, sc_group = \
            [], [], [], [], [], [], [], [], [], []

        for i, j in zip(inst, names):
            if re.match('^LPW$', i.strip()):
                lpw_group.append(j)
            elif re.match('^LPW-EUV$', i.strip()):
                euv_group.append(j)
            elif re.match('^SWEA$', i.strip()):
                swe_group.append(j)
            elif re.match('^SWIA$', i.strip()):
                swi_group.append(j)
            elif re.match('^STATIC$', i.strip()):
                sta_group.append(j)
            elif re.match('^SEP$', i.strip()):
                sep_group.append(j)
            elif re.match('^MAG$', i.strip()):
                mag_group.append(j)
            elif re.match('^NGIMS$', i.strip()):
                ngi_group.append(j)
            elif re.match('^SPICE$', i.strip()):
                # NB Need to split into APP and SPACECRAFT
                if re.match('(.+)APP(.+)', j):
                    app_group.append(j)
                else:  # Everything not APP is SC in SPICE
                    # But do not include Orbit Num, or IO Flag
                    # Could probably stand to clean this line up a bit
                    if not re.match('(.+)(Orbit Number|Inbound Outbound Flag)', j):
                        sc_group.append(j)
            else:
                pass

        delete_groups = []
        if instruments is not None:
            if 'LPW' not in instruments and 'lpw' not in instruments:
                delete_groups += lpw_group
            if 'MAG' not in instruments and 'mag' not in instruments:
                delete_groups += mag_group
            if 'EUV' not in instruments and 'euv' not in instruments:
                delete_groups += euv_group
            if 'SWI' not in instruments and 'swi' not in instruments:
                delete_groups += swi_group
            if 'SWE' not in instruments and 'swe' not in instruments:
                delete_groups += swe_group
            if 'NGI' not in instruments and 'ngi' not in instruments:
                delete_groups += ngi_group
            if 'SEP' not in instruments and 'sep' not in instruments:
                delete_groups += sep_group
            if 'STA' not in instruments and 'sta' not in instruments:
                delete_groups += sta_group

        # Read in all relavent data into a pandas dataframe
        temp_data = []
        for filename in filenames:
            # Determine number of header lines
            nheader = 0
            with open(filename) as f:
                for line in f:
                    if line.startswith('#'):
                        nheader += 1

                temp_data.append(pd.read_fwf(filename, skiprows=nheader, index_col=0,
                                             widths=[19] + len(names) * [16], names=names))
                for i in delete_groups:
                    del temp_data[-1][i]

        temp_unconverted = pd.concat(temp_data)

        # Need to convert columns
        if 'SWEA.Electron Spectrum Shape' in temp_unconverted and 'NGIMS.Density NO' in temp_unconverted:
            temp = temp_unconverted.astype(dtype={'SWEA.Electron Spectrum Shape': np.float64,
                                                  'NGIMS.Density NO': np.float64})
        elif 'SWEA.Electron Spectrum Shape' in temp_unconverted and 'NGIMS.Density NO' not in temp_unconverted:
            temp = temp_unconverted.astype(dtype={'SWEA.Electron Spectrum Shape': np.float64})
        elif 'SWEA.Electron Spectrum Shape' not in temp_unconverted and 'NGIMS.Density NO' in temp_unconverted:
            temp = temp_unconverted.astype(dtype={'NGIMS.Density NO': np.float64})
        else:
            temp = temp_unconverted

        # Assign the first-level only tags
        time_unix = [calendar.timegm(datetime.strptime(i, '%Y-%m-%dT%H:%M:%S').timetuple()) for i in temp.index]
        time = temp.index
        time_unix = pd.Series(time_unix)  # convert into Series for consistency
        time_unix.index = temp.index
        orbit = temp['SPICE.Orbit Number']
        io_flag = temp['SPICE.Inbound Outbound Flag']

        # Build the sub-level DataFrames for the larger dictionary/structure
        app = temp[app_group]
        spacecraft = temp[sc_group]
        if instruments is not None:
            if 'LPW' in instruments or 'lpw' in instruments:
                lpw = temp[lpw_group]
            else:
                lpw = None
            if 'MAG' in instruments or 'mag' in instruments:
                mag = temp[mag_group]
            else:
                mag = None
            if 'EUV' in instruments or 'euv' in instruments:
                euv = temp[euv_group]
            else:
                euv = None
            if 'SWE' in instruments or 'swe' in instruments:
                swea = temp[swe_group]
            else:
                swea = None
            if 'SWI' in instruments or 'swi' in instruments:
                swia = temp[swi_group]
            else:
                swia = None
            if 'NGI' in instruments or 'ngi' in instruments:
                ngims = temp[ngi_group]
            else:
                ngims = None
            if 'SEP' in instruments or 'sep' in instruments:
                sep = temp[sep_group]
            else:
                sep = None
            if 'STA' in instruments or 'sta' in instruments:
                static = temp[sta_group]
            else:
                static = None
        else:
            lpw = temp[lpw_group]
            euv = temp[euv_group]
            swea = temp[swe_group]
            swia = temp[swi_group]
            static = temp[sta_group]
            sep = temp[sep_group]
            mag = temp[mag_group]
            ngims = temp[ngi_group]

        # Strip out the duplicated instrument part of the column names
        for i in [lpw, euv, swea, swia, sep, static, ngims, mag, app, spacecraft]:
            if i is not None:
                i.columns = remove_inst_tag(i)

        if lpw is not None:
            lpw = lpw.rename(index=str, columns=param_dict)
        if euv is not None:
            euv = euv.rename(index=str, columns=param_dict)
        if swea is not None:
            swea = swea.rename(index=str, columns=param_dict)
        if swia is not None:
            swia = swia.rename(index=str, columns=param_dict)
        if sep is not None:
            sep = sep.rename(index=str, columns=param_dict)
        if static is not None:
            static = static.rename(index=str, columns=param_dict)
        if ngims is not None:
            ngims = ngims.rename(index=str, columns=param_dict)
        if mag is not None:
            mag = mag.rename(index=str, columns=param_dict)
        if app is not None:
            app = app.rename(index=str, columns=param_dict)
        if spacecraft is not None:
            spacecraft = spacecraft.rename(index=str, columns=param_dict)

        # Do not forget to save units
        # Define the list of first level tag names
        tag_names = ['TimeString', 'Time', 'Orbit', 'IOflag',
                     'LPW', 'EUV', 'SWEA', 'SWIA', 'STATIC',
                     'SEP', 'MAG', 'NGIMS', 'APP', 'SPACECRAFT']
        # Define list of first level data structures
        data_tags = [time, time_unix, orbit, io_flag,
                     lpw, euv, swea, swia, static,
                     sep, mag, ngims, app, spacecraft]
        kp_insitu = (OrderedDict(zip(tag_names, data_tags)))

        return kp_insitu

    else:
        # Now for IUVS
        kp_iuvs = []
        for file in filenames:
            kp_iuvs.append(read_iuvs_file(file))
        return kp_iuvs

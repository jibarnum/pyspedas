# -*- coding: utf-8 -*-
"""
File:
    utilities.py

Description:
    Functions used by read_kp.
"""

import re
import collections


def remove_inst_tag(df):
    '''
    Remove the leading part of the column name that includes the instrument
    identifier for use in creating the parameter names for the toolkit.

    Input:
        A DataFrame produced from the insitu KP data

    Output:
        A new set of column names
    '''

    newcol = []
    for i in df.columns:
        if len(i.split('.')) >= 2:
            j = i.split('.')
            newcol.append('.'.join(j[1:]))

    return newcol


def get_header_info(filename):
    """
    Retrieve column names for a given KP in situ file

    Input:
        A KP in situ filename

    Output:
        Column names and instrument names
    """
    # Determine number of header lines
    nheader = 0
    with open(filename) as f:
        for line in f:
            if line.startswith('#'):
                nheader += 1

    # Parse the header (still needs special case work)
    read_param_list = False
    index_list = []
    with open(filename) as fin:
        icol = -2  # Counting header lines detailing column names
        iname = 1  # for counting seven lines with name info
        ncol = -1  # Dummy value to allow reading of early headerlines?
        col_regex = '#\s(.{16}){%3d}' % ncol  # needed for column names
        for iline in range(nheader):
            line = fin.readline()
            if re.search('Number of parameter columns', line):
                ncol = int(re.split("\s{3}", line)[1])
                col_regex = '#\s(.{16}){%3d}' % ncol  # needed for column names
            elif re.search('Line on which data begins', line):
                nhead_test = int(re.split("\s{3}", line)[1]) - 1
            elif re.search('Number of lines', line):
                ndata = int(re.split("\s{3}", line)[1])
            elif re.search('PARAMETER', line):
                read_param_list = True
                param_head = iline
            elif read_param_list:
                icol += 1
                if icol > ncol:
                    read_param_list = False
            elif re.match(col_regex, line):
                # OK, verified match now get the values
                temp = re.findall('(.{16})', line[3:])
                if iname == 1:
                    index = temp
                elif iname == 2:
                    obs1 = temp
                elif iname == 3:
                    obs2 = temp
                elif iname == 4:
                    obs3 = temp
                elif iname == 5:
                    inst = temp
                elif iname == 6:
                    unit = temp
                elif iname == 7:
                    format_code = temp
                else:
                    print('More lines in data descriptor than expected.')
                    print('Line %d' % iline)
                iname += 1
            else:
                pass

        # Generate the names list.
        # NB, there are special case redundancies in there
        # (e.g., LPW: Electron Density Quality (min and max))
        # ****SWEA FLUX electron QUALITY *****
        first = True
        parallel = None
        names = []
        for h, i, j, k in zip(inst, obs1, obs2, obs3):
            combo_name = (' '.join([i.strip(), j.strip(), k.strip()])).strip()
            if re.match('^LPW$', h.strip()):
                # Max and min error bars use same name in column
                # SIS says first entry is min and second is max
                if re.match('(Electron|Spacecraft)(.+)Quality', combo_name):
                    if first:
                        combo_name = combo_name + ' Min'
                        first = False
                    else:
                        combo_name = combo_name + ' Max'
                        first = True
            elif re.match('^SWEA$', h.strip()):
                # electron flux qual flags do not indicate whether parallel or anti
                # From context it is clear; but we need to specify in name
                if re.match('.+Parallel.+', combo_name):
                    parallel = True
                elif re.match('.+Anti-par', combo_name):
                    parallel = False
                else:
                    pass
                if re.match('Flux, e-(.+)Quality', combo_name):
                    if parallel:
                        p = re.compile('Flux, e- ')
                        combo_name = p.sub('Flux, e- Parallel ', combo_name)
                    else:
                        p = re.compile('Flux, e- ')
                        combo_name = p.sub('Flux, e- Anti-par ', combo_name)
                if re.match('Electron eflux (.+)Quality', combo_name):
                    if parallel:
                        p = re.compile('Electron eflux ')
                        combo_name = p.sub('Electron eflux  Parallel ', combo_name)
                    else:
                        p = re.compile('Electron eflux ')
                        combo_name = p.sub('Electron eflux  Anti-par ', combo_name)
            # Add inst to names to avoid ambiguity
            # Will need to remove these after splitting
            names.append('.'.join([h.strip(), combo_name]))
            names[0] = 'Time'

    return names, inst


def read_iuvs_file(file):
    """
    Read in a KP IUVS file

    Input:
        KP IUVS file

    Output:
        Dictionary with IUVS data
    """
    iuvs_dict = {}
    periapse_num = 0
    occ_num = 0
    with open(file) as f:
        line = f.readline()
        while line is not '':
            if line.startswith('*'):
                # Read the header
                line = f.readline()
                obs_mode = line[19:len(line) - 1].strip()
                
                header = {}
                f.readline()
                line = f.readline()
                header['time_start'] = line[19:len(line) - 1].strip()
                line = f.readline()
                header['time_stop'] = line[19:len(line) - 1].strip()
                line = f.readline()
                if obs_mode == "OCCULTATION":
                    header['target_name'] = line[19:len(line) - 1].strip()
                    line = f.readline()
                header['sza'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['local_time'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['lat'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['lon'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['lat_mso'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['lon_mso'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['orbit_number'] = int(line[19:len(line) - 1].strip())
                line = f.readline()
                header['mars_season'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_geo_x'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_geo_y'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_geo_z'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_mso_x'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_mso_y'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_mso_z'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sun_geo_x'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sun_geo_y'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sun_geo_z'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sun_geo_lat'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sun_geo_lon'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sun_mso_lat'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sun_mso_lon'] = float(line[19:len(line) - 1].strip())
                line = f.readline()     
                header['subsol_geo_lon'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['subsol_geo_lat'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_sza'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_local_time'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['sc_alt'] = float(line[19:len(line) - 1].strip())
                line = f.readline()
                header['mars_sun_dist'] = float(line[19:len(line) - 1].strip())
                
                if obs_mode == "PERIAPSE":
                    periapse_num += 1
                    line = f.readline()
                    n_alt_bins = int(line[19:len(line) - 1].strip())
                    header['n_alt_bins'] = float(n_alt_bins)
                    line = f.readline()
                    n_alt_den_bins = int(line[19:len(line) - 1].strip())
                    header['n_alt_den_bins'] = float(n_alt_den_bins)
                    
                    iuvs_dict['periapse' + str(periapse_num)] = {}
                    iuvs_dict['periapse' + str(periapse_num)].update(header)
                    
                    # Empty space
                    f.readline()
                    
                    # Read the Temperature
                    line = f.readline()
                    temp_labels = line[19:len(line) - 1].strip().split()
                    temperature = collections.OrderedDict((x, []) for x in temp_labels)
                    temperature_unc = collections.OrderedDict((x, []) for x in temp_labels)
                    line = f.readline()
                    vals = line[20:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        temperature[list(temperature.keys())[index]].append(val)
                        index += 1
                    line = f.readline()
                    vals = line[20:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        temperature_unc[list(temperature_unc.keys())[index]].append(val)
                        index += 1
                    iuvs_dict['periapse' + str(periapse_num)]['temperature'] = temperature
                    iuvs_dict['periapse' + str(periapse_num)]['temperature_unc'] = temperature_unc
                    
                    # Empty space
                    f.readline()
                    
                    # Read the Scale Heights
                    line = f.readline()
                    scale_height_labels = line[19:len(line) - 1].strip().split()
                    scale_height = collections.OrderedDict((x, []) for x in scale_height_labels)
                    scale_height_unc = collections.OrderedDict((x, []) for x in scale_height_labels)
                    line = f.readline()
                    vals = line[20:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        scale_height[list(scale_height.keys())[index]].append(val)
                        index += 1
                    line = f.readline()
                    vals = line[20:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        scale_height_unc[list(scale_height_unc.keys())[index]].append(val)
                        index += 1
                    
                    iuvs_dict['periapse' + str(periapse_num)]['scale_height'] = scale_height
                    iuvs_dict['periapse' + str(periapse_num)]['scale_height_unc'] = scale_height_unc
                    
                    # Empty space
                    f.readline()
                    f.readline()
                    
                    # Read in the density
                    line = f.readline()
                    density_labels = line.strip().split()
                    density = collections.OrderedDict((x, []) for x in density_labels)
                    for i in range(0, n_alt_den_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            density[list(density.keys())[index]].append(val)
                            index += 1
                    iuvs_dict['periapse' + str(periapse_num)]['density'] = density
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the density systematic uncertainty
                    density_sys_unc = collections.OrderedDict((x, []) for x in density_labels)
                    line = f.readline()
                    vals = line.strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        density_sys_unc[list(density.keys())[index + 1]].append(val)
                        index += 1
                        
                    iuvs_dict['periapse' + str(periapse_num)]['density_sys_unc'] = density_sys_unc
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the density uncertainty
                    density_unc = collections.OrderedDict((x, []) for x in density_labels)
                    for i in range(0, n_alt_den_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            density_unc[list(density.keys())[index]].append(val)
                            index += 1
                    iuvs_dict['periapse' + str(periapse_num)]['density_sys_unc'] = density_sys_unc
                    
                    f.readline()
                    f.readline()
                    
                    line = f.readline()
                    radiance_labels = line.strip().split()
                    if "Cameron" in radiance_labels:
                        radiance_labels.remove('Cameron')
                    radiance = collections.OrderedDict((x, []) for x in radiance_labels)
                    for i in range(0, n_alt_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            radiance[list(radiance.keys())[index]].append(val)
                            index += 1
                            
                    iuvs_dict['periapse' + str(periapse_num)]['radiance'] = radiance
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the radiance systematic uncertainty
                    radiance_sys_unc = collections.OrderedDict((x, []) for x in radiance_labels)
                    line = f.readline()
                    vals = line.strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        radiance_sys_unc[list(radiance.keys())[index + 1]].append(val)
                        index += 1
                    
                    iuvs_dict['periapse' + str(periapse_num)]['radiance_sys_unc'] = radiance_sys_unc
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the radiance uncertainty
                    radiance_unc = collections.OrderedDict((x, []) for x in radiance_labels)
                    for i in range(0, n_alt_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            radiance_unc[list(radiance.keys())[index]].append(val)
                            index += 1
                            
                    iuvs_dict['periapse'+str(periapse_num)]['radiance_unc'] = radiance_unc
                            
                elif obs_mode == "OCCULTATION":
                    occ_num += 1
                    line = f.readline()
                    n_alt_den_bins = int(line[19:len(line) - 1].strip())
                    header['n_alt_den_bins'] = float(n_alt_den_bins)
                    
                    iuvs_dict['occultation' + str(occ_num)] = {}
                    iuvs_dict['occultation' + str(occ_num)].update(header)
                    
                    # Empty space
                    f.readline()
                    
                    # Read the Scale Heights
                    line = f.readline()
                    scale_height_labels = line[19:len(line) - 1].strip().split()
                    scale_height = collections.OrderedDict((x, []) for x in scale_height_labels)
                    scale_height_unc = collections.OrderedDict((x, []) for x in scale_height_labels)
                    line = f.readline()
                    vals = line[20:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        scale_height[list(scale_height.keys())[index]].append(val)
                        index += 1
                    line = f.readline()
                    vals = line[20:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        scale_height_unc[list(scale_height_unc.keys())[index]].append(val)
                        index += 1
                    
                    iuvs_dict['occultation' + str(occ_num)]['scale_height'] = scale_height
                    iuvs_dict['occultation' + str(occ_num)]['scale_height_unc'] = scale_height_unc
                    
                    # Empty space
                    f.readline()
                    f.readline()
                    
                    # Read in the retrieval
                    line = f.readline()
                    retrieval_labels = line.strip().split()
                    retrieval = collections.OrderedDict((x, []) for x in retrieval_labels)
                    for i in range(0, n_alt_den_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            retrieval[list(retrieval.keys())[index]].append(val)
                            index += 1
                    iuvs_dict['occultation' + str(occ_num)]['retrieval'] = retrieval
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the retrieval systematic uncertainty
                    retrieval_sys_unc = collections.OrderedDict((x, []) for x in retrieval_labels)
                    line = f.readline()
                    vals = line.strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        retrieval_sys_unc[list(retrieval.keys())[index + 1]].append(val)
                        index += 1
                        
                    iuvs_dict['occultation' + str(occ_num)]['retrieval_sys_unc'] = retrieval_sys_unc
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the retrieval uncertainty
                    retrieval_unc = collections.OrderedDict((x, []) for x in retrieval_labels)
                    for i in range(0, n_alt_den_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            retrieval_unc[list(retrieval.keys())[index]].append(val)
                            index += 1
                    iuvs_dict['occultation' + str(occ_num)]['retrieval_sys_unc'] = retrieval_sys_unc

                elif obs_mode == "CORONA_LORES_HIGH":
                    line = f.readline()
                    n_alt_bins = int(line[19:len(line) - 1].strip())
                    header['n_alt_bins'] = float(n_alt_bins)
                    
                    iuvs_dict['corona_lores_high'] = {}
                    iuvs_dict['corona_lores_high'].update(header)

                    f.readline()
                    
                    # Read the Half int
                    line = f.readline()
                    half_int_dist_labels = line[19:len(line) - 1].strip().split()
                    half_int_dist = collections.OrderedDict((x, []) for x in half_int_dist_labels)
                    half_int_dist_unc = collections.OrderedDict((x, []) for x in half_int_dist_labels)
                    line = f.readline()
                    vals = line[26:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        half_int_dist[list(half_int_dist.keys())[index]].append(val)
                        index += 1
                    line = f.readline()
                    vals = line[26:len(line) - 1].strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        half_int_dist_unc[list(half_int_dist_unc.keys())[index]].append(val)
                        index += 1
                        
                    iuvs_dict['corona_lores_high']['half_int_dist'] = half_int_dist
                    iuvs_dict['corona_lores_high']['half_int_dist_unc'] = half_int_dist_unc
                    
                    # Blank space
                    f.readline()
                    f.readline()
                    
                    # Read in the density
                    line = f.readline()
                    density_labels = line.strip().split()
                    density = collections.OrderedDict((x, []) for x in density_labels)
                    for i in range(0, n_alt_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            density[list(density.keys())[index]].append(val)
                            index += 1
                    
                    iuvs_dict['corona_lores_high']['density'] = density
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the density systematic uncertainty
                    density_sys_unc = collections.OrderedDict((x, []) for x in density_labels)
                    line = f.readline()
                    vals = line.strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        density_sys_unc[list(density.keys())[index + 1]].append(val)
                        index += 1
                        
                    iuvs_dict['corona_lores_high']['density_sys_unc'] = density_sys_unc
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()

                    # Read in the density uncertainty
                    density_unc = collections.OrderedDict((x, []) for x in density_labels)
                    for i in range(0, n_alt_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            density_unc[list(density.keys())[index]].append(val)
                            index += 1
                            
                    iuvs_dict['corona_lores_high']['density_unc'] = density_unc
                    
                    f.readline()
                    f.readline()
                    
                    line = f.readline()
                    radiance_labels = line.strip().split()
                    if "Cameron" in radiance_labels:
                        radiance_labels.remove('Cameron')
                    radiance = collections.OrderedDict((x, []) for x in radiance_labels)
                    for i in range(0, n_alt_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            radiance[list(radiance.keys())[index]].append(val)
                            index += 1
                            
                    iuvs_dict['corona_lores_high']['radiance'] = radiance
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the radiance systematic uncertainty
                    radiance_sys_unc = collections.OrderedDict((x, []) for x in radiance_labels)
                    line = f.readline()
                    vals = line.strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        radiance_sys_unc[list(radiance.keys())[index + 1]].append(val)
                        index += 1
                        
                    iuvs_dict['corona_lores_high']['radiance_sys_unc'] = radiance_sys_unc
                    
                    # Not needed lines
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the radiance uncertainty
                    radiance_unc = collections.OrderedDict((x, []) for x in radiance_labels)
                    for i in range(0, n_alt_bins):
                        line = f.readline()
                        vals = line.strip().split()
                        index = 0
                        for val in vals:
                            if val == '-9.9999990E+09':
                                val = float('nan')
                            else:
                                val = float(val)
                            radiance_unc[list(radiance.keys())[index]].append(val)
                            index += 1
                            
                    iuvs_dict['corona_lores_high']['radiance_unc'] = radiance_unc
                    
                elif obs_mode == 'APOAPSE':

                    f.readline()
                    maps = {}
                    for j in range(0, 17):
                        var = f.readline().strip()
                        line = f.readline()
                        lons = line.strip().split()
                        lons = [float(x) for x in lons]
                        lats = []
                        data = []
                        for k in range(0, 45):
                            line = f.readline().strip().split()
                            lats.append(float(line[0]))
                            line_data = line[1:]
                            line_data = [float(x) if x != '-9.9999990E+09' else float('nan') for x in line_data]
                            data.append(line_data)
                            
                        maps[var] = data
                        f.readline()
    
                    maps['latitude'] = lats
                    maps['longitude'] = lons
                    
                    iuvs_dict['apoapse'] = {}
                    iuvs_dict['apoapse'].update(header)
                    iuvs_dict['apoapse'].update(maps)
                    
                    f.readline()
                    f.readline()
                    f.readline()
                    
                    # Read in the radiance systematic uncertainty
                    line = f.readline()
                    radiance_labels = line.strip().split()
                    radiance_sys_unc = collections.OrderedDict((x, []) for x in radiance_labels)
                    line = f.readline()
                    vals = line.strip().split()
                    index = 0
                    for val in vals:
                        if val == '-9.9999990E+09':
                            val = float('nan')
                        else:
                            val = float(val)
                        radiance_sys_unc[list(radiance.keys())[index + 1]].append(val)
                        index += 1
                    
                    iuvs_dict['apoapse']['radiance_sys_unc'] = radiance_sys_unc
                    
            line = f.readline()    
    
    return iuvs_dict


param_dict = {'Electron Density': 'ELECTRON_DENSITY',
              'Electron Density Quality Min': 'ELECTRON_DENSITY_QUAL_MIN',
              'Electron Density Quality Max': 'ELECTRON_DENSITY_QUAL_MAX',
              'Electron Temperature': 'ELECTRON_TEMPERATURE',
              'Electron Temperature Quality Min': 'ELECTRON_TEMPERATURE_QUAL_MIN',
              'Electron Temperature Quality Max': 'ELECTRON_TEMPERATURE_QUAL_MAX',
              'Spacecraft Potential': 'SPACECRAFT_POTENTIAL',
              'Spacecraft Potential Quality Min':  'SPACECRAFT_POTENTIAL_QUAL_MIN',
              'Spacecraft Potential Quality Max':  'SPACECRAFT_POTENTIAL_QUAL_MAX',
              'E-field Power 2-100 Hz':  'EWAVE_LOW_FREQ',
              'E-field 2-100 Hz Quality':  'EWAVE_LOW_FREQ_QUAL_QUAL',
              'E-field Power 100-800 Hz':  'EWAVE_MID_FREQ',
              'E-field 100-800 Hz Quality':  'EWAVE_MID_FREQ_QUAL_QUAL',
              'E-field Power 0.8-1.0 Mhz':  'EWAVE_HIGH_FREQ',
              'E-field 0.8-1.0 Mhz Quality':  'EWAVE_HIGH_FREQ_QUAL_QUAL',
              'EUV Irradiance 0.1-7.0 nm':  'IRRADIANCE_LOW',
              'Irradiance 0.1-7.0 nm Quality':  'IRRADIANCE_LOW_QUAL',
              'EUV Irradiance 17-22 nm':  'IRRADIANCE_MID',
              'Irradiance 17-22 nm Quality':  'IRRADIANCE_MID_QUAL',
              'EUV Irradiance Lyman-alpha':  'IRRADIANCE_LYMAN',
              'Irradiance Lyman-alpha Quality':  'IRRADIANCE_LYMAN_QUAL',
              'Solar Wind Electron Density':  'SOLAR_WIND_ELECTRON_DENSITY',
              'Solar Wind E- Density Quality':  'SOLAR_WIND_ELECTRON_DENSITY_QUAL',
              'Solar Wind Electron Temperature':  'SOLAR_WIND_ELECTRON_TEMPERATURE',
              'Solar Wind E- Temperature Quality':  'SOLAR_WIND_ELECTRON_TEMPERATURE_QUAL',
              'Flux, e- Parallel (5-100 ev)':  'ELECTRON_PARALLEL_FLUX_LOW',
              'Flux, e- Parallel (5-100 ev) Quality':  'ELECTRON_PARALLEL_FLUX_LOW_QUAL',
              'Flux, e- Parallel (100-500 ev)':  'ELECTRON_PARALLEL_FLUX_MID',
              'Flux, e- Parallel (100-500 ev) Quality':  'ELECTRON_PARALLEL_FLUX_MID_QUAL',
              'Flux, e- Parallel (500-1000 ev)':  'ELECTRON_PARALLEL_FLUX_HIGH',
              'Flux, e- Parallel (500-1000 ev) Quality':  'ELECTRON_PARALLEL_FLUX_HIGH_QUAL',
              'Flux, e- Anti-par (5-100 ev)':  'ELECTRON_ANTI_PARALLEL_FLUX_LOW',
              'Flux, e- Anti-par (5-100 ev) Quality':  'ELECTRON_ANTI_PARALLEL_FLUX_LOW_QUAL',
              'Flux, e- Anti-par (100-500 ev)':  'ELECTRON_ANTI_PARALLEL_FLUX_MID',
              'Flux, e- Anti-par (100-500 ev) Quality':  'ELECTRON_ANTI_PARALLEL_FLUX_MID_QUAL',
              'Flux, e- Anti-par (500-1000 ev)':  'ELECTRON_ANTI_PARALLEL_FLUX_HIGH',
              'Flux, e- Anti-par (500-1000 ev) Quality':  'ELECTRON_ANTI_PARALLEL_FLUX_HIGH_QUAL',
              'Electron eflux Parallel (5-100 ev)':  'ELECTRON_PARALLEL_FLUX_LOW',
              'Electron eflux Parallel (5-100 ev) Quality':  'ELECTRON_PARALLEL_FLUX_LOW_QUAL',
              'Electron eflux Parallel (100-500 ev)':  'ELECTRON_PARALLEL_FLUX_MID',
              'Electron eflux Parallel (100-500 ev) Quality':  'ELECTRON_PARALLEL_FLUX_MID_QUAL',
              'Electron eflux Parallel (500-1000 ev)':  'ELECTRON_PARALLEL_FLUX_HIGH',
              'Electron eflux Parallel (500-1000 ev) Quality':  'ELECTRON_PARALLEL_FLUX_HIGH_QUAL',
              'Electron eflux Anti-par (5-100 ev)':  'ELECTRON_ANTI_PARALLEL_FLUX_LOW',
              'Electron eflux Anti-par (5-100 ev) Quality':  'ELECTRON_ANTI_PARALLEL_FLUX_LOW_QUAL',
              'Electron eflux Anti-par (100-500 ev)':  'ELECTRON_ANTI_PARALLEL_FLUX_MID',
              'Electron eflux Anti-par (100-500 ev) Quality':  'ELECTRON_ANTI_PARALLEL_FLUX_MID_QUAL',
              'Electron eflux Anti-par (500-1000 ev)':  'ELECTRON_ANTI_PARALLEL_FLUX_HIGH',
              'Electron eflux Anti-par (500-1000 ev) Quality':  'ELECTRON_ANTI_PARALLEL_FLUX_HIGH_QUAL',
              'Electron Spectrum Shape':  'ELECTRON_SPECTRUM_SHAPE_PARAMETER',
              'Spectrum Shape Quality':  'ELECTRON_SPECTRUM_SHAPE_PARAMETER_QUAL',
              'H+ Density':  'HPLUS_DENSITY',
              'H+ Density Quality':  'HPLUS_DENSITY_QUAL',
              'H+ Flow Velocity MSO X':  'HPLUS_FLOW_VELOCITY_MSO_X',
              'H+ Flow MSO X Quality':  'HPLUS_FLOW_VELOCITY_MSO_X_QUAL',
              'H+ Flow Velocity MSO Y':  'HPLUS_FLOW_VELOCITY_MSO_Y',
              'H+ Flow MSO Y Quality':  'HPLUS_FLOW_VELOCITY_MSO_Y_QUAL',
              'H+ Flow Velocity MSO Z':  'HPLUS_FLOW_VELOCITY_MSO_Z',
              'H+ Flow MSO Z Quality':  'HPLUS_FLOW_VELOCITY_MSO_Z_QUAL',
              'Solar Wind Dynamic Pressure':  'SOLAR_WIND_DYNAMIC_PRESSURE',
              'Solar Wind Pressure Quality':  'SOLAR_WIND_DYNAMIC_PRESSURE_QUAL',
              'STATIC Quality Flag':  'STATIC_QUALITY_FLAG',
              'O+ Density':  'OPLUS_DENSITY',
              'O+ Density Quality':  'OPLUS_DENSITY_QUAL',
              'O2+ Density':  'O2PLUS_DENSITY',
              'O2+ Density Quality':  'O2PLUS_DENSITY_QUAL',
              'H+ Temperature':  'HPLUS_TEMPERATURE',
              'H+ Temperature Quality':  'HPLUS_TEMPERATURE_QUAL',
              'O+ Temperature':  'OPLUS_TEMPERATURE',
              'O+ Temperature Quality':  'OPLUS_TEMPERATURE_QUAL',
              'O2+ Temperature':  'O2PLUS_TEMPERATURE',
              'O2+ Temperature Quality':  'O2PLUS_TEMPERATURE_QUAL',
              'O2+ Flow Velocity MAVEN_APP X':  'O2PLUS_FLOW_VELOCITY_MAVEN_APP_X',
              'O2+ Flow MAVEN_APP X Quality':  'O2PLUS_FLOW_VELOCITY_MAVEN_APP_X_QUAL',
              'O2+ Flow Velocity MAVEN_APP Y':  'O2PLUS_FLOW_VELOCITY_MAVEN_APP_Y',
              'O2+ Flow MAVEN_APP Y Quality':  'O2PLUS_FLOW_VELOCITY_MAVEN_APP_Y_QUAL',
              'O2+ Flow Velocity MAVEN_APP Z':  'O2PLUS_FLOW_VELOCITY_MAVEN_APP_Z',
              'O2+ Flow MAVEN_APP Z Quality':  'O2PLUS_FLOW_VELOCITY_MAVEN_APP_Z_QUAL',
              'O2+ Flow Velocity MSO X':  'O2PLUS_FLOW_VELOCITY_MSO_X',
              'O2+ Flow MSO X Quality':  'O2PLUS_FLOW_VELOCITY_MSO_X_QUAL',
              'O2+ Flow Velocity MSO Y':  'O2PLUS_FLOW_VELOCITY_MSO_Y',
              'O2+ Flow MSO Y Quality':  'O2PLUS_FLOW_VELOCITY_MSO_Y_QUAL',
              'O2+ Flow Velocity MSO Z':  'O2PLUS_FLOW_VELOCITY_MSO_Z',
              'O2+ Flow MSO Z Quality':  'O2PLUS_FLOW_VELOCITY_MSO_Z_QUAL',
              'H+ Omni Flux':  'HPLUS_OMNI_DIRECTIONAL_FLUX',
              'H+ Energy':  'HPLUS_CHARACTERISTIC_ENERGY',
              'H+ Energy Quality':  'HPLUS_CHARACTERISTIC_ENERGY_QUAL',
              'He++ Omni Flux':  'HEPLUS_OMNI_DIRECTIONAL_FLUX',
              'He++ Energy':  'HEPLUS_CHARACTERISTIC_ENERGY',
              'He++ Energy Quality':  'HEPLUS_CHARACTERISTIC_ENERGY_QUAL',
              'O+ Omni Flux':  'OPLUS_OMNI_DIRECTIONAL_FLUX',
              'O+ Energy':  'OPLUS_CHARACTERISTIC_ENERGY',
              'O+ Energy Quality':  'OPLUS_CHARACTERISTIC_ENERGY_QUAL',
              'O2+ Omni Flux':  'O2PLUS_OMNI_DIRECTIONAL_FLUX',
              'O2+ Energy':  'O2PLUS_CHARACTERISTIC_ENERGY',
              'O2+ Energy Quality':  'O2PLUS_CHARACTERISTIC_ENERGY_QUAL',
              'H+ Direction MSO X':  'HPLUS_CHARACTERISTIC_DIRECTION_MSO_X',
              'H+ Direction MSO Y':  'HPLUS_CHARACTERISTIC_DIRECTION_MSO_Y',
              'H+ Direction MSO Z':  'HPLUS_CHARACTERISTIC_DIRECTION_MSO_Z',
              'H+ Angular Width':  'HPLUS_CHARACTERISTIC_ANGULAR_WIDTH',
              'H+ Width Quality':  'HPLUS_CHARACTERISTIC_ANGULAR_WIDTH_QUAL',
              'Pickup Ion Direction MSO X':  'DOMINANT_PICKUP_ION_CHARACTERISTIC_DIRECTION_MSO_X',
              'Pickup Ion Direction MSO Y':  'DOMINANT_PICKUP_ION_CHARACTERISTIC_DIRECTION_MSO_Y',
              'Pickup Ion Direction MSO Z':  'DOMINANT_PICKUP_ION_CHARACTERISTIC_DIRECTION_MSO_Z',
              'Pickup Ion Angular Width':  'DOMINANT_PICKUP_ION_CHARACTERISTIC_ANGULAR_WIDTH',
              'Pickup Ion Width Quality':  'DOMINANT_PICKUP_ION_CHARACTERISTIC_ANGULAR_WIDTH_QUAL',
              'Ion Flux FOV 1 F':  'ION_ENERGY_FLUX__FOV_1_F',
              'Ion Flux FOV 1F Quality':  'ION_ENERGY_FLUX__FOV_1_F_QUAL',
              'Ion Flux FOV 1 R':  'ION_ENERGY_FLUX__FOV_1_R',
              'Ion Flux FOV 1R Quality':  'ION_ENERGY_FLUX__FOV_1_R_QUAL',
              'Ion Flux FOV 2 F':  'ION_ENERGY_FLUX__FOV_2_F',
              'Ion Flux FOV 2F Quality':  'ION_ENERGY_FLUX__FOV_2_F_QUAL',
              'Ion Flux FOV 2 R':  'ION_ENERGY_FLUX__FOV_2_R',
              'Ion Flux FOV 2R Quality':  'ION_ENERGY_FLUX__FOV_2_R_QUAL',
              'Electron Flux FOV 1 F':  'ELECTRON_ENERGY_FLUX___FOV_1_F',
              'Electron Flux FOV 1F Quality':  'ELECTRON_ENERGY_FLUX___FOV_1_F_QUAL',
              'Electron Flux FOV 1 R':  'ELECTRON_ENERGY_FLUX___FOV_1_R',
              'Electron Flux FOV 1R Quality':  'ELECTRON_ENERGY_FLUX___FOV_1_R_QUAL',
              'Electron Flux FOV 2 F':  'ELECTRON_ENERGY_FLUX___FOV_2_F',
              'Electron Flux FOV 2F Quality':  'ELECTRON_ENERGY_FLUX___FOV_2_F_QUAL',
              'Electron Flux FOV 2 R':  'ELECTRON_ENERGY_FLUX___FOV_2_R',
              'Electron Flux FOV 2R Quality':  'ELECTRON_ENERGY_FLUX___FOV_2_R_QUAL',
              'Look Direction 1-F MSO X':  'LOOK_DIRECTION_1_F_MSO_X',
              'Look Direction 1-F MSO Y':  'LOOK_DIRECTION_1_F_MSO_Y',
              'Look Direction 1-F MSO Z':  'LOOK_DIRECTION_1_F_MSO_Z',
              'Look Direction 1-R MSO X':  'LOOK_DIRECTION_1_R_MSO_X',
              'Look Direction 1-R MSO Y':  'LOOK_DIRECTION_1_R_MSO_Y',
              'Look Direction 1-R MSO Z':  'LOOK_DIRECTION_1_R_MSO_Z',
              'Look Direction 2-F MSO X':  'LOOK_DIRECTION_2_F_MSO_X',
              'Look Direction 2-F MSO Y':  'LOOK_DIRECTION_2_F_MSO_Y',
              'Look Direction 2-F MSO Z':  'LOOK_DIRECTION_2_F_MSO_Z',
              'Look Direction 2-R MSO X':  'LOOK_DIRECTION_2_R_MSO_X',
              'Look Direction 2-R MSO Y':  'LOOK_DIRECTION_2_R_MSO_Y',
              'Look Direction 2-R MSO Z':  'LOOK_DIRECTION_2_R_MSO_Z',
              'Magnetic Field MSO X':  'MSO_X',
              'Magnetic MSO X Quality':  'MSO_X_QUAL',
              'Magnetic Field MSO Y':  'MSO_Y',
              'Magnetic MSO Y Quality':  'MSO_Y_QUAL',
              'Magnetic Field MSO Z':  'MSO_Z',
              'Magnetic MSO Z Quality':  'MSO_Z_QUAL',
              'Magnetic Field GEO X':  'GEO_X',
              'Magnetic GEO X Quality':  'GEO_X_QUAL',
              'Magnetic Field GEO Y':  'GEO_Y',
              'Magnetic GEO Y Quality':  'GEO_Y_QUAL',
              'Magnetic Field GEO Z':  'GEO_Z',
              'Magnetic GEO Z Quality':  'GEO_Z_QUAL',
              'Magnetic Field RMS Dev':  'RMS_DEVIATION',
              'Magnetic RMS Quality':  'RMS_DEVIATION_QUAL',
              'Density He':  'HE_DENSITY',
              'Density He Precision':  'HE_DENSITY_PRECISION',
              'Density He Quality':  'HE_DENSITY_QUAL',
              'Density O':  'O_DENSITY',
              'Density O Precision':  'O_DENSITY_PRECISION',
              'Density O Quality':  'O_DENSITY_QUAL',
              'Density CO':  'CO_DENSITY',
              'Density CO Precision':  'CO_DENSITY_PRECISION',
              'Density CO Quality':  'CO_DENSITY_QUAL',
              'Density N2':  'N2_DENSITY',
              'Density N2 Precision':  'N2_DENSITY_PRECISION',
              'Density N2 Quality':  'N2_DENSITY_QUAL',
              'Density NO':  'NO_DENSITY',
              'Density NO Precision':  'NO_DENSITY_PRECISION',
              'Density NO Quality':  'NO_DENSITY_QUAL',
              'Density Ar':  'AR_DENSITY',
              'Density Ar Precision':  'AR_DENSITY_PRECISION',
              'Density Ar Quality':  'AR_DENSITY_QUAL',
              'Density CO2':  'CO2_DENSITY',
              'Density CO2 Precision':  'CO2_DENSITY_PRECISION',
              'Density CO2 Quality':  'CO2_DENSITY_QUAL',
              'Density 32+':  'O2PLUS_DENSITY',
              'Density 32+ Precision':  'O2PLUS_DENSITY_PRECISION',
              'Density 32+ Quality':  'O2PLUS_DENSITY_QUAL',
              'Density 44+':  'CO2PLUS_DENSITY',
              'Density 44+ Precision':  'CO2PLUS_DENSITY_PRECISION',
              'Density 44+ Quality':  'CO2PLUS_DENSITY_QUAL',
              'Density 30+':  'NOPLUS_DENSITY',
              'Density 30+ Precision':  'NOPLUS_DENSITY_PRECISION',
              'Density 30+ Quality':  'NOPLUS_DENSITY_QUAL',
              'Density 16+':  'OPLUS_DENSITY',
              'Density 16+ Precision':  'OPLUS_DENSITY_PRECISION',
              'Density 16+ Quality':  'OPLUS_DENSITY_QUAL',
              'Density 28+':  'CO2PLUS_N2PLUS_DENSITY',
              'Density 28+ Precision':  'CO2PLUS_N2PLUS_DENSITY_PRECISION',
              'Density 28+ Quality':  'CO2PLUS_N2PLUS_DENSITY_QUAL',
              'Density 12+':  'CPLUS_DENSITY',
              'Density 12+ Precision':  'CPLUS_DENSITY_PRECISION',
              'Density 12+ Quality':  'CPLUS_DENSITY_QUAL',
              'Density 17+':  'OHPLUS_DENSITY',
              'Density 17+ Precision':  'OHPLUS_DENSITY_PRECISION',
              'Density 17+ Quality':  'OHPLUS_DENSITY_QUAL',
              'Density 14+':  'NPLUS_DENSITY',
              'Density 14+ Precision':  'NPLUS_DENSITY_PRECISION',
              'Density 14+ Quality':  'NPLUS_DENSITY_QUAL',
              'APP Attitude GEO X':  'ATTITUDE_GEO_X',
              'APP Attitude GEO Y':  'ATTITUDE_GEO_Y',
              'APP Attitude GEO Z':  'ATTITUDE_GEO_Z',
              'APP Attitude MSO X':  'ATTITUDE_MSO_X',
              'APP Attitude MSO Y':  'ATTITUDE_MSO_Y',
              'APP Attitude MSO Z':  'ATTITUDE_MSO_Z',
              'Spacecraft GEO X':  'GEO_X',
              'Spacecraft GEO Y':  'GEO_Y',
              'Spacecraft GEO Z':  'GEO_Z',
              'Spacecraft MSO X':  'MSO_X',
              'Spacecraft MSO Y':  'MSO_Y',
              'Spacecraft MSO Z':  'MSO_Z',
              'Spacecraft GEO Longitude':  'SUB_SC_LONGITUDE',
              'Spacecraft GEO Latitude':  'SUB_SC_LATITUDE',
              'Spacecraft Solar Zenith Angle':  'SZA',
              'Spacecraft Local Time':  'LOCAL_TIME',
              'Spacecraft Altitude Aeroid':  'ALTITUDE',
              'Spacecraft Attitude GEO X':  'ATTITUDE_GEO_X',
              'Spacecraft Attitude GEO Y':  'ATTITUDE_GEO_Y',
              'Spacecraft Attitude GEO Z':  'ATTITUDE_GEO_Z',
              'Spacecraft Attitude MSO X':  'ATTITUDE_MSO_X',
              'Spacecraft Attitude MSO Y':  'ATTITUDE_MSO_Y',
              'Spacecraft Attitude MSO Z':  'ATTITUDE_MSO_Z',
              'Mars Season (Ls)':  'MARS_SEASON',
              'Mars-Sun Distance':  'MARS_SUN_DISTANCE',
              'Subsolar Point GEO Longitude':  'SUBSOLAR_POINT_GEO_LONGITUDE',
              'Subsolar Point GEO Latitude':  'SUBSOLAR_POINT_GEO_LATITUDE',
              'Sub-Mars Point on the Sun Longitude':  'SUBMARS_POINT_SOLAR_LONGITUDE',
              'Sub-Mars Point on the Sun Latitude':  'SUBMARS_POINT_SOLAR_LATITUDE',
              'Rot matrix MARS -> MSO Row 1, Col 1':  'T11',
              'Rot matrix MARS -> MSO Row 1, Col 2':  'T12',
              'Rot matrix MARS -> MSO Row 1, Col 3':  'T13',
              'Rot matrix MARS -> MSO Row 2, Col 1':  'T21',
              'Rot matrix MARS -> MSO Row 2, Col 2':  'T22',
              'Rot matrix MARS -> MSO Row 2, Col 3':  'T23',
              'Rot matrix MARS -> MSO Row 3, Col 1':  'T31',
              'Rot matrix MARS -> MSO Row 3, Col 2':  'T32',
              'Rot matrix MARS -> MSO Row 3, Col 3':  'T33',
              'Rot matrix SPCCRFT -> MSO Row 1, Col 1':  'SPACECRAFT_T11',
              'Rot matrix SPCCRFT -> MSO Row 1, Col 2':  'SPACECRAFT_T12',
              'Rot matrix SPCCRFT -> MSO Row 1, Col 3':  'SPACECRAFT_T13',
              'Rot matrix SPCCRFT -> MSO Row 2, Col 1':  'SPACECRAFT_T21',
              'Rot matrix SPCCRFT -> MSO Row 2, Col 2':  'SPACECRAFT_T22',
              'Rot matrix SPCCRFT -> MSO Row 2, Col 3':  'SPACECRAFT_T23',
              'Rot matrix SPCCRFT -> MSO Row 3, Col 1':  'SPACECRAFT_T31',
              'Rot matrix SPCCRFT -> MSO Row 3, Col 2':  'SPACECRAFT_T32',
              'Rot matrix SPCCRFT -> MSO Row 3, Col 3':  'SPACECRAFT_T33'}

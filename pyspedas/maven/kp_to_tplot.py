# -*- coding: utf-8 -*-
"""
File:
    kp_to_tplot.py

Description:
    Functions used by maven_load to transform KP in situ python dicts into tplot variables.

Parameters:
    dates: str/list of str ['yyyy-mm-dd']
        List of dates to be downloaded (eg. ['2015-12-31']).
"""

import pytplot


def kp_to_tplot(insitu):
    """Creates tplot variables from the insitu variable
    """
    # Keep track of stored KP variables
    stored_variables = []

    # initialize each instrument
    inst_list = ["EUV", "LPW", "STATIC", "SWEA", "SWIA", "MAG", "SEP", "NGIMS"]

    for instrument in inst_list:
        # for each observation for each instrument
        if insitu[instrument] is not None:
            for obs in insitu[instrument]:
                # create variable name
                obs_specific = "mvn_kp::" + instrument.lower() + "::" + obs.lower()
                # if NaN or string, continue
                if insitu[instrument][obs].isnull().all() or insitu[instrument][obs].dtype == 'O':
                    continue
                # store data in tplot variable
                pytplot.store_data(obs_specific, data={'x': insitu['Time'], 'y': insitu[instrument][obs]})
                if obs_specific not in stored_variables:
                    stored_variables.append(obs_specific)

    return stored_variables

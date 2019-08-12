"""
This module contains routines for loading MAVEN data.
"""

from .maven_load import load_data


def maven_load(instruments=None,
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
    Main function for downloading MAVEN data and loading it into tplot variables (if applicable).

    Parameters:
        instruments: str/list of str
            Instruments from which you want to download data. This is where you can indicate that you want to
            download KP data (via 'kp-insitu' for KP in situ data or 'kp-iuvs' for KP iuvs data).
        kp_instruments: str/list of str
            Instruments from which you want to grab KP in situ data. Only needed if you're downloading KP in situ data,
            and is optional (if you don't specify, and you chose to download KP in situ data,
            KP in situ data will be downloaded for all instruments with in situ data).
        list_files: bool (True/False0
            If true, lists the files instead of downloading them.
        new_files: bool (True/False)
            Checks downloaded files and only downloads those that haven't already been downloaded.
        start_date: str
            String that is the start date for downloading data (YYYY-MM-DD)
        end_date: str
            String that is the end date for downloading data (YYYY-MM-DD)
        update_prefs: bool (True/False)
            If true, updates where you want to store data locally
        only_update_prefs: bool (True/False)
            If true, *only* updates where to store dat alocally, doesn't download files.
        exclude_orbit_file: bool (True/False)
            If true, won't download the latest orbit tables.
        local_dir: str
            If indicated, specifies where to download files for a specific implementation of this function.
        download_only: bool (True/False)
            If True then files are downloaded only,
            if False then CDF files are also loaded into pytplot using cdf_to_tplot.
        varformat : str
            The file variable formats to load into tplot.  Wildcard character
            "*" is accepted.  By default, all variables are loaded in.
        prefix: str
            The tplot variable names will be given this prefix.
            By default, no prefix is added.
        suffix: str
            The tplot variable names will be given this suffix.
            By default, no suffix is added.
        get_support_data: bool
            Data with an attribute "VAR_TYPE" with a value of "support_data"
            will be loaded into tplot.  By default, only loads in data with a
            "VAR_TYPE" attribute of "data".
        public: bool
            If True, downloads data from the MAVEN public website.
            If False, downloads data from the MAVEN private website (will ask for username/password).
    """
    tvars = load_data(instruments=instruments, kp_instruments=kp_instruments, start_date=start_date, end_date=end_date,
                      update_prefs=update_prefs, only_update_prefs=only_update_prefs, local_dir=local_dir,
                      list_files=list_files, new_files=new_files, exclude_orbit_file=exclude_orbit_file,
                      download_only=download_only, varformat=varformat, prefix=prefix, suffix=suffix,
                      get_support_data=get_support_data, public=public)
    return tvars

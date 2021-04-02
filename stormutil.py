# Last modified: Time-stamp: <2020-12-09 14:02:49 haines>
""" Storm utilities (stormutil)

"""

import os
import re
import time
import datetime

import netCDF4

import numpy as np
import metpy
import metpy.calc
from metpy.units import units


def scanf_datetime(ts, fmt='%Y-%m-%dT%H:%M:%S'):
    """Convert string representing date and time to datetime object"""
    # default string format follows convention YYYY-MM-DDThh:mm:ss
    
    try:
        t = time.strptime(ts, fmt)
        # the '*' operator unpacks the tuple, producing the argument list.
        dt = datetime.datetime(*t[0:6])
    except ValueError as e:
        # value error if something not valid for datetime
        # e.g. month 1...12, something parsed wrong
        dt = None
        
    return dt

def find_months(year, month=1):
    """Find prev, this, and next month to process

    :Parameters:
        year : int value or str 'yyyy_mm'
        month : int value
    :Returns:
        which_months : list of  datetime objects
             [this_month, next_month]
    Examples
    --------
    >>> find_months(2007, 2)
    >>> find_months('2007_02')
    
    """
    if type(year) == int and type(month) == int :
        dt = datetime.datetime(year, month, day=1)
        this_month = dt
    elif type(year) == str :
        dt = scanf_datetime(year, fmt='%Y_%m')
        this_month = dt
    #
    if dt.month == 1: # if January
        # prev_month = datetime.datetime(dt.year-1, month=12, day=1) # Dec
        next_month = datetime.datetime(dt.year, dt.month+1, day=1) # Feb
    elif dt.month == 12: # if December
        # prev_month = datetime.datetime(dt.year, dt.month-1, day=1) # Nov
        next_month = datetime.datetime(dt.year+1, month=1, day=1)  # Jan
    else:
        # prev_month = datetime.datetime(dt.year, dt.month-1, day=1)
        next_month = datetime.datetime(dt.year, dt.month+1, day=1)
    #
    # return (prev_month, this_month, next_month)
    return [this_month, next_month]

def update_storms(storms):
    """
    Add key for computed datetime (dt) of each position of storm to use in find_storms()
    """
    for ed in range(len(storms)):
        age = storms[ed]['age']
        dt = np.array([], dtype=np.double)
        for idx in range(age):
            yyyy=int(storms[ed]['year'][idx])
            mm=int(storms[ed]['month'][idx])
            dd=int(storms[ed]['day'][idx])
            HH=int(storms[ed]['hour'][idx])
            dt = np.append(dt, datetime.datetime(yyyy,mm,dd,HH,0,0))
        storms[ed]['dt'] = dt
    
    return storms

def find_storms(storms, dt, dtidx=0):
    """
    Find storms and tracks at the given time of dtidx
    """
    # dt[dtidx] datetime of current time
    # storms[ed]['dt'] all datetimes of the track calc'd in update_storms()

    track_lons = np.array([], dtype=np.float)
    track_lats = np.array([], dtype=np.float)
    storm_lons = np.array([], dtype=np.float)
    storm_lats = np.array([], dtype=np.float)
    
    for ed in range(len(storms)):
        # if these are not found or skipping may have to add range within small tolerance
        which = (dt[dtidx] == storms[ed]['dt']) 
        if which.any() * (storms[ed]['type'] == 'cyclonic'):
            track_lons = np.append(track_lons, storms[ed]['lon'])
            track_lons = np.append(track_lons, np.nan) # append NaN to pickup pen for line
            #
            track_lats = np.append(track_lats, storms[ed]['lat'])
            track_lats = np.append(track_lats, np.nan) # append NaN to pickup pen for line 
            #
            idx, = which.nonzero()[0] 
            storm_lons = np.append(storm_lons, storms[ed]['lon'][idx]) 
            storm_lats = np.append(storm_lats, storms[ed]['lat'][idx])
           

    # container for current storm positions and storm tracks
    ss = dict()
    ss['track_lons'] = track_lons
    ss['track_lats'] = track_lats
    ss['storm_lons'] = storm_lons
    ss['storm_lats'] = storm_lats
    
    return ss


def get_data(indir, BB):
    """ Read in 4d-var ERA5 data

    Parameter
    ---------
    indir : string
       The input directory path. Data loaded by param and by year.
       All param and year files (param.YYYY.nc) in indir must be
       of the same space and time (time, lvl, lat, lon).
    BB : dictionary 
       Requires 4 keys (lat,lon,lvl,dt)
       Each key has value [min, max]

    Returns
    -------
    d : dict of ndarrays and computed quantities

    """

    # key words are used in filename but values are names with netcdf file
    # let's only load hgt and msl for stormviz
    # uwnd, and vwnd have data for 100 to 500 hPa, but hgt has levels 100 to 1000 hPa
    # removed from params
    #           'uwnd': 'u_component_of_wind',
    #           'vwnd': 'v_component_of_wind',
    params = {'hgt' : 'geopotential',
              'msl' : 'mean_sea_level_pressure'
              }
    
    sfc_params = ['msl']
    press_params = ['hgt', 'uwnd', 'vwnd']
    dt1 = BB['dt'][0]
    dt2 = BB['dt'][1]

    #
    print('Reading ERA5 data from: %s' % indir)

    for param in list(params.keys()):
        fn = '%s.%04d.nc' % (param, dt1.year) # each file year has one param
        # ifn = os.path.join(indir, fn)
        ifn = '/'.join([indir, fn])
        nc = netCDF4.Dataset(ifn)
        varnames = list(nc.variables.keys())
        print(varnames)
        t = nc.variables['time']
        dt = netCDF4.num2date(t[:], units=t.units, calendar=t.calendar)
        lat = nc.variables['latitude'][:].data
        lon = nc.variables['longitude'][:].data
        
        # nonzero returns a tuple of idx per dimension
        # we're unpacking the tuple for each of these idx-vars
        (dtidx,) = np.logical_and(dt >= BB['dt'][0], dt < BB['dt'][1]).nonzero()
        (latidx,) = np.logical_and(lat >= BB['lat'][0], lat <= BB['lat'][1]).nonzero()
        (lonidx,) = np.logical_and(lon >= BB['lon'][0], lon <= BB['lon'][1]).nonzero()
       
        if param in press_params:
            level = nc.variables['level'][:].data
            (levidx,) =  np.logical_and(level >= BB['lvl'][0], level <= BB['lvl'][1]).nonzero()
            level_units = nc.variables['level'].units
        # get subset of data from file
        if param=='uwnd':
            uwnd = nc.variables['u'][dtidx, levidx, latidx, lonidx].data * units(nc.variables['u'].units)
        elif param=='vwnd':
            vwnd = nc.variables['v'][dtidx, levidx, latidx, lonidx].data * units(nc.variables['v'].units)
        elif param=='hgt':
            geopot = nc.variables['z'][dtidx, levidx, latidx, lonidx].data * units(nc.variables['z'].units)
        elif param=='msl':
            msl = nc.variables['msl'][dtidx, latidx, lonidx].data * units(nc.variables['msl'].units).to('hPa')
        # close the param datafile
        nc.close()

    # -------------------------------
    # do a calculation using metpy functions -- wspd(dt,level,lat,lon)
    # wspd = metpy.calc.wind_speed(uwnd, vwnd)
    # compute geopotential height -- hgt(dt,level,lat,lon)
    hgt = metpy.calc.geopotential_to_height(geopot)
    
    # compute height (1d) based on standard pressure -- ht(level)
    ht_std = metpy.calc.pressure_to_height_std(level * units(level_units))
    # we will be adjusting ht from ht_std by difference in msl from
    # std pressure (1013.25 * units.hPa) at sea level
    p0 = 1013.25 * units.hPa
    # difference in pressure to add for each msl(dt,lat,lon)
    pdiff = msl-p0

    # collection of data for plots
    d = dict()
    # dimensions within BB -- nparrays
    d['dt'] = dt[dtidx]
    d['lat']= lat[latidx]
    d['lon']= lon[lonidx]
    d['level']=level[levidx]
    # data already subset to BB and units attached -- Quantity object (nparray * units)
    d['ht_std'] = ht_std # ht(level)
    d['msl'] = msl       # msl(dt,lat,lon)
    d['hgt'] = hgt       # hgt(dt,level,lat,lon)
    # d['uwnd']= uwnd      # uwnd(dt,level,lat,lon)
    # d['vwnd']= vwnd      # vwnd(dt,level,lat,lon)
    # d['wspd']= wspd      # wspd(dt,level,lat,lon)
    d['pdiff']=pdiff     # pdiff(dt,lat,lon)

    return d

def get_coastlines():
    # --------------------------
    # Global Self-consistent, Hierarchical, High-resolution Shoreline Database (gshhs)
    # http://opendap.deltares.nl/thredds/catalog/opendap/noaa/gshhs/catalog.html
    # --- 5 resolutions ----- netCDF versions available on deltars.nl
    # gshhs_f = finest
    # gshhs_h = high
    # gshhs_i = intermediate
    # gshhs_l = low
    # gshhs_c = coarse
    # lineurl  = 'http://opendap.deltares.nl/thredds/dodsC/opendap/noaa/gshhs/gshhs_i.nc';
    # lineurl  = 'http://whewell.marine.unc.edu/dods/gshhs/gshhs_i.nc'
    lineurl  = 'http://whewell.marine.unc.edu/dods/gshhs/gshhs_l.nc'

    # Get coatline line data: 1D vectors are small, so we can get all data
    # opendap(url_line) # when netCDF4 was not compiled with OPeNDAP
    linedata = netCDF4.Dataset(lineurl)
    
    lines = dict(
        lon=linedata.variables['lon'][:],
        lat=linedata.variables['lat'][:]
        )
    linedata.close()
    # -----------------------------
    return lines

def generate_columns(types_str='JSDT JSLVL JSLAT JSLON JSHT WSPD UWND VWND HGT'):
    # use dict to store column label and it's column number
    #c = col.defaultdict(int)
    c = {}
    column_labels = types_str.strip().split(' ')
    m = re.findall(r'\w{2,}', types_str) # 2 or more char per word
    for label in column_labels:
        c[label]=m.index(label) # c['VFLG']=4
    return c

def write_jet_data(ofn, header, js):
    """Write header and data. """
    f = open(ofn, 'w')
    if header[-1] == '\n':
        f.write(header)
    else:
        f.write(header+'\n')
    # if there is any data, save to the file)
    if js.size > 0:
        np.savetxt(f, js, fmt='%s')
    f.close()

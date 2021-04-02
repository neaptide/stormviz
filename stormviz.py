# Last modified: Time-stamp: <2021-03-30 08:31:20 haines>
r""" Storm vizualization (stormviz) tool using ECMWF Reanalysis v5 (ERA5) data
and storm tracks from Eric Oliver's stormTracking code.

Plots:

(1) For a given time, show map (pmap) of mean sea level pressure (msl).
Show complete track path of any storm that has any point for the given time. (line with small solid marker)
Show current storm position (red dot) fot the given time. 

GUI:
   Time slider with prev and next buttons

Usage:
Using IPython console, use magic to run code as if at unix prompt and
provide year and month to view e.g.
%run stormviz.py [yyyy_mm]

Start ipython in era5 python environment
(era5) C:\Users\haines>ipython

In[]: cd Dropbox/peach/era5/stormviz
In[]: %run stormviz.py 2018_01
In[]: plt.show()

"""

import sys
from stormutil import *

import matplotlib.pyplot as plt
import matplotlib.gridspec as gs
from matplotlib.widgets import Slider, Button, TextBox, CheckButtons

# suppress warnings
import warnings
warnings.filterwarnings("ignore")

# Define default data bounds for analysis
BB = dict( lon=[-140,   0],
           lat=[   0,  80],
           lvl=[ 100, 1000],
           dt = [datetime.datetime(2017,1,1), datetime.datetime(2017,2,1)]
           )

# Define default bounding box for the figure
BB_fig = dict( lon=[-140,   0],
               lat=[  10,  80],
               lvl=[ 100, 1000],
               dt = [datetime.datetime(2017,1,1), datetime.datetime(2017,2,1)]
               )

# grab the coastline dataset
lines = get_coastlines()

# empty array for storm data
storm = np.array([])

# setup figure layout 
fig = plt.figure(figsize=(10, 7.5))
axs = [fig.add_axes((.1,.1,.8,.7)),0,0]

# main map
t1 = axs[0].set_title('YYYY_MM_DD_HHMM', loc='left')
title_str = 'hgt (%d hPa), msl pressure (hPa)' % 500
t2 = axs[0].set_title(title_str, loc='right')
# set aspect to simply mimic equidistant projection
axs[0].set_aspect(1/np.cos(np.pi*np.mean(BB_fig['lat'])/180.)) 
axs[0].set_xlim(BB_fig['lon'][0],BB_fig['lon'][1])
axs[0].set_ylim(BB_fig['lat'][0],BB_fig['lat'][1])
axs[0].set_xlabel('Longitude (deg)')
axs[0].set_ylabel('Latitude (deg)')
# plot coastline/lakes
axs[0].plot(lines['lon'],lines['lat'],'k',linewidth=0.5)

# wpsd color bar
# get_positions returns Bbox, we want Bbox.bounds
l,b,w,h = axs[0].get_position().bounds
axs[1] = fig.add_axes([l,b-0.075,w,0.02])

# do this now so that we can get adjusted ax get_position
plt.draw()
plt.pause(0.01)

# blank data for initiating contours
blank = np.array(np.ones((2,2), dtype=float))

# contour lines for hmap(lon,lat)
toggle_hgt_on = True
cslines1 = np.arange(0, 21000, 100)
cs11 = axs[0].contour(blank, blank, blank, cslines1, colors='k', linewidths=1.0, linestyles='solid')

# plot filled contour for pmap(lon,lat) 
toggle_msl_on = True
cflines = np.arange(900, 1110, 10)
# cmap = plt.cm.get_cmap('BuPu')
cmap = plt.cm.get_cmap('PuOr')
cf1 = axs[0].contourf(blank, blank, blank, cflines, cmap=cmap)
# contour lines for pmap(lon,lat)
cs12 = axs[0].contour(blank, blank, blank, np.arange(900,1015,5),colors='gray', linewidths=1.0, linestyles='dashed')
cs13 = axs[0].contour(blank, blank, blank, np.arange(1015,1110,5), colors='gray', linewidths=1.0, linestyles='solid')

# plot track and storm locations on map
trackmap, = axs[0].plot([],[], 'bo-', markersize=4)
stormmap, = axs[0].plot([],[], 'ro', markersize=6)


def update_plot(val):
    # when dt slider changes
    global ss,stormmap,trackmap,cf1,cs11,cs12,cs13
    dtidx = int(sdt.val)
    lvlidx = int(slvl.val)

    # find storms each time step
    ss = find_storms(storms, d['dt'], dtidx)
    # print(ss)
    stormmap.set_ydata(ss['storm_lats'])
    stormmap.set_xdata(ss['storm_lons'])
    
    trackmap.set_ydata(ss['track_lats'])
    trackmap.set_xdata(ss['track_lons'])
    
    dt_str = d['dt'][dtidx].strftime("%Y_%m_%d_%H%M")
    t1.set_text(dt_str)
    title_str = 'hgt (%d hPa), msl pressure (hPa)' % d['level'][lvlidx]
    t2.set_text(title_str)
    lons, lats = np.meshgrid(d['lon'], d['lat'])

    # pick level (hPa) of hgt
    hmap = d['hgt'][dtidx,lvlidx,:,:].squeeze()
    
    pmap = d['msl'][dtidx,:,:].squeeze()
    
    if toggle_hgt_on:
        # remove all previous contours and labels
        remove_hgt_layer()
    
        # contour lines for hmap(lon,lat)
        cs11 = axs[0].contour(lons, lats, hmap, cslines1, colors='k', linewidths=1.0, linestyles='solid')
        cs11_lab = axs[0].clabel(cs11, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
                                 rightside_up=True, use_clabeltext=True)
    
    if toggle_msl_on:
        # remove all previous contours and labels
        remove_msl_layer()
        # plot filled contour for pmap(lon,lat) # using clabel on cf1 fill contour does weird stuff
        cf1 = axs[0].contourf(lons, lats, pmap, cflines, cmap=cmap)
        
        # contour lines for pmap(lon,lat) (low pmap<1013 dashed) (high pmap>1013 solid)
        # ever recorded lowest 870 hPa (typhoon), highest 1085 hPa
        cs12 = axs[0].contour(lons, lats, pmap, np.arange(900,1015,5),colors='gray', linewidths=1.0, linestyles='dashed')
        cs12_lab = axs[0].clabel(cs12, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
                                 rightside_up=True, use_clabeltext=True)
        cs13 = axs[0].contour(lons, lats, pmap, np.arange(1015,1110,5), colors='gray', linewidths=1.0, linestyles='solid')
        cs13_lab = axs[0].clabel(cs13, fontsize=8, inline=1, inline_spacing=10, fmt='%i',
                                 rightside_up=True, use_clabeltext=True)
    plt.draw()


def remove_hgt_layer():
    global cs11
    c = cs11.collections
    for tp in c:
        try:
            tp.remove()
        except ValueError:
            pass
    l = cs11.clabel()
    for lb in l:
        try:
            lb.remove()
        except ValueError:
            pass

def remove_msl_layer():
    global cf1, cs12, cs13
    c = cf1.collections
    c.extend(cs12.collections)
    c.extend(cs13.collections)
    for tp in c:
        try:
            tp.remove()
        except ValueError:
            pass
    l = cf1.clabel()
    l.extend(cs12.clabel())
    l.extend(cs13.clabel())
    for lb in l:
        try:
            lb.remove()
        except ValueError:
            pass


def prev_lvl(val):
    lvlidx = int(slvl.val)
    slvl.set_val(lvlidx-1)

def next_lvl(val):
    lvlidx = int(slvl.val)
    slvl.set_val(lvlidx+1)

def prev_dt(val):
    dtidx = int(sdt.val)
    sdt.set_val(dtidx-1)

def next_dt(val):
    dtidx = int(sdt.val)
    sdt.set_val(dtidx+1)

def toggle_storm_tracks(val):
    global stormmap, trackmap
    stormmap.set_visible(not stormmap.get_visible())
    trackmap.set_visible(not trackmap.get_visible())
    if stormmap.get_visible():
        cjs1.label.set_text('Storm Tracks ON')
        cjs1.ax.set_facecolor('green')
    else:
        cjs1.label.set_text('Storm Tracks OFF')
        cjs1.ax.set_facecolor('red')
    plt.draw()

def toggle_msl(val):
    global toggle_msl_on
    toggle_msl_on = (not toggle_msl_on)
    if toggle_msl_on:
        cjs2.label.set_text('MSL Pressure -- ON')
        cjs2.ax.set_facecolor('green')
    else:
        remove_msl_layer()
        cjs2.label.set_text('MSL Pressure -- OFF')
        cjs2.ax.set_facecolor('red')
    update_plot(val)

def toggle_hgt(val):
    global toggle_hgt_on
    toggle_hgt_on = (not toggle_hgt_on)
    if toggle_hgt_on:
        cjs3.label.set_text('HGT -- ON')
        cjs3.ax.set_facecolor('green')
    else:
        remove_hgt_layer()
        cjs3.label.set_text('HGT -- OFF')
        cjs3.ax.set_facecolor('red')
    update_plot(val)

# outer grid to frame inner grid of gui, 
# use the object handle of figure (fig) and method add_gridspec
ogs = fig.add_gridspec(5,4, left=0.05, right=0.95,  top=0.95, bottom=0.05)
# otherwise this direct call in jupyter-notebooks put the grid and widgets in new figure
# ogs = gs.GridSpec(4,3, left=0.1, right=0.95,  top=0.95, bottom=0.05)

# use top row of ogs for inner grids
# ogs[0,0] for storm track on/off
# ogs[0,1] for sliders
# ogs[0,2] fpr next and prev button sets


igs = gs.GridSpecFromSubplotSpec(4,2,subplot_spec=ogs[0,0], hspace=0.1)
# storm track hide/show
cjs1 = Button(fig.add_subplot(igs[0,0:]), 
              label='Storm Tracks -- ON', color='green', hovercolor='green')
cjs1.on_clicked(toggle_storm_tracks)

# MSL hide/show 
cjs2 = Button(fig.add_subplot(igs[1,0:]), 
              label='MSL Pressure -- ON', color='green', hovercolor='green')
cjs2.on_clicked(toggle_msl)

# HGT hide/show  
cjs3 = Button(fig.add_subplot(igs[2,0:]), 
              label='HGT -- ON', color='green', hovercolor='green')
cjs3.on_clicked(toggle_hgt)

igs = gs.GridSpecFromSubplotSpec(4,1,subplot_spec=ogs[0,1], hspace=0.2)
# Level slider
# use the object handle of figure (fig) and method add_subplot to add
axlvl = fig.add_subplot(igs[0])
slvl = Slider(axlvl, 'Level', 0, 22, valinit=0, valfmt='%d')
slvl.on_changed(update_plot)

# Date slider 
axdt = fig.add_subplot(igs[1])
sdt = Slider(axdt, 'Date', 0, 31*4, valinit=0, valfmt='%d')
sdt.on_changed(update_plot)

igs = gs.GridSpecFromSubplotSpec(4,4,subplot_spec=ogs[0,2], hspace=0.2)
# Level prev button
axlvlprev = fig.add_subplot(igs[0,0])
blvlprev = Button(axlvlprev, '<')
blvlprev.on_clicked(prev_lvl)
# Level next button
axlvlnext = fig.add_subplot(igs[0,1])
blvlnext = Button(axlvlnext, '>')
blvlnext.on_clicked(next_lvl)
# Date prev button
axdtprev = fig.add_subplot(igs[1,0])
bdtprev = Button(axdtprev, '<')
bdtprev.on_clicked(prev_dt)
# Date next button
axdtnext = fig.add_subplot(igs[1,1])
bdtnext = Button(axdtnext, '>')
bdtnext.on_clicked(next_dt)


def init_plot():
    """ initialize plots, finish setting up, and set slider limits
    """
    global ss,stormmap,trackmap,cf1,cs12,cs13,sdt,slvl
    dtidx = 0

    dt_str = d['dt'][dtidx].strftime("%Y_%m_%d_%H%M")
    t1.set_text(dt_str)

    sdt.valinit = dtidx
    sdt.valmin = 0
    sdt.valmax = len(d['dt'])-1

    level = 500
    (lvlidx,) = (d['level']==level).nonzero()
    slvl.valinit = lvlidx
    slvl.valmin = 0
    slvl.valmax = len(d['level'])-1
    slvl.valstep = 1

    update_plot(0)

    # plot map and set up colorbar
    # draw colorbar
    cb = fig.colorbar(cf1, cax=axs[1], orientation='horizontal') 
    cb.set_label('MSL Pressure (hPa)')


if len(sys.argv)==2:
    yyyy_mm = sys.argv[1]
else:
    yyyy_mm = '2018_01'

BB['dt'] = find_months(yyyy_mm)
BB_fig['dt'] = find_months(yyyy_mm)

# input path of netcdf files
# local data
# indir = os.path.join('/data', 'era5', 'test')
# d = get_data(indir, BB)

# use data on dap server
# dapdir = 'http://whewell.marine.unc.edu/dods/era5/test' # 10/60 N
dapdir = 'http://whewell.marine.unc.edu/dods/era5' # 0/80 N
d = get_data(dapdir, BB)

# Load storm data (local file in current working directory
# py3 needs allow_pickle=True
# data = np.load('storm_track_slp.npz',allow_pickle=True)
data = np.load('../storm_tracking_ERA5/npz_data_1979-2019/storm_track_ERA5_msl.npz',allow_pickle=True)
# data = np.load('../storm_tracking_ERA5/storm_track_ERA5_msl.npz',allow_pickle=True)
# data = np.load('../storm_tracking_ERA5/storm_track_ERA5_msl_peach.npz',allow_pickle=True)
# data = np.load('../storm_tracking_ERA5/storm_track_ERA5_hgt_500.npz',allow_pickle=True)
# data = np.load('../storm_tracking_ERA5/storm_track_ERA5_hgt_850.npz',allow_pickle=True)
# data = np.load('../storm_tracking_ERA5/storm_track_ERA5_hgt_1000.npz',allow_pickle=True)

# get storms and update with dt key field
storms = update_storms(data['storms'])

init_plot()
plt.draw()

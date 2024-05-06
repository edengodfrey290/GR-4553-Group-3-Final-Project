


"""

TO RUN IN SPYDER USE: runfile('pathtofile/MetPy_skewT_winds_hodograph.py', args='STATIONID YYMMDDHH')
where 
-pathtofile: the folder where the file is (you still have to specify this even if you are running in the current folder)
-STATIONID: 3 letter station ID, e.g. JAN
-YYMMDDHH: the date/time of interest, eg. 2023033112

Example script usage to get a sounding WITHOUT a parcel path:
    python MetPy_skewT_winds_hodograph.py JAN 2023033112

Example script usage to get a sounding WITH the default parcel path (100-hPa mixed layer parcel starting at surface):
    python MetPy_skewT_winds_hodograph.py JAN 2023033112 true

Example script usage to get a sounding WITH a different parcel path (provided mixed layer parcel starting at surface):
(uses lowest 50 hPa as the mixed layer instead of the lowest 100 -- example can be changed to 90 or other values)
    python MetPy_skewT_winds_hodograph.py JAN 2023033112 50

Example script usage to get a sounding with a SURFACE parcel path:
    python MetPy_skewT_winds_hodograph.py JAN 2023033112 sfc

Example script usage to get a sounding with a SURFACE parcel path **and** a mixed-layer parcel path:
    python MetPy_skewT_winds_hodograph.py JAN 2023033112 sfc 100

Package versions:
    python 3.11.0, matplotlib 3.6.2, metpy 1.3.1, numpy 1.23.4, pandas 1.5.1, scipy 1.9.3, siphon 0.9
        (can use 'skewTenv.yml' with conda to create an environment:  conda env create -f skewTenv.yml)
"""
from datetime import datetime
from metpy.plots import SkewT, Hodograph
from metpy.units import pandas_dataframe_to_unit_arrays, units
from scipy.interpolate import interp1d
from siphon.simplewebservice.wyoming import WyomingUpperAir
import matplotlib.pyplot as plt
import metpy.calc as mpcalc
import numpy as np
import pandas as pd
import sys

station = sys.argv[1].upper()                   # examples: JAN, TUS
dt = datetime.strptime(sys.argv[2],'%Y%m%d%H')  # example: 2023033112 (12 UTC 31 March 2023)
plotParcelPath = False
mixedParcel = False
surfaceParcel = False
if len(sys.argv) > 3:
    plotParcelPath = True
    arguments = sys.argv[3:]
    ## check for 'sfc' request
    try:
        ii = arguments.index('sfc')
        surfaceParcel = True
        del arguments[ii]  # remove 'sfc' from the 'arguments' list
    except ValueError:
        mixedParcel = True
    if len(arguments) == 1 or mixedParcel is True:
        mixedParcel = True
        ## define settings for a mixed-layer parcel
        mixedLayerDepth = 100
        if len(arguments) > 0:
            try:
                mixedLayerDepth = int(arguments[0])  # use provided value to re-define the mixed layer depth in hPa
            except ValueError:
                pass  # make no changes

# ----------------------------------------------------------------------------------------------------------------------
# DO NOT CHANGE ANYTHING BELOW THIS COMMENT (unless you're comfortable with experimenting!)
# ----------------------------------------------------------------------------------------------------------------------
df = WyomingUpperAir.request_data(dt, station)
data = pandas_dataframe_to_unit_arrays(df)
p = data['pressure']
end = max(np.where(p.m>95.)[0]) + 1  # limit data to pressure levels >= 95 hPa
p = data['pressure'][0:end]
T = data['temperature'][0:end]
Td = data['dewpoint'][0:end]
u = data['u_wind'][0:end]
v = data['v_wind'][0:end]
hgt = data['height'][0:end]

#### Do a data quality check for repeated height and/or pressure values (which can break the script's plotting)
check = np.ediff1d(p.m)
if check.min() <= 0.:
    bad = np.where(check==0)[0]
    if len(bad) > 0:
        p_new = np.delete(p.m, bad)
        p = p_new * data['pressure'].units  # assign proper units to the corrected 'pressure' array
        T_new = np.delete(T.m, bad)
        T = T_new * data['temperature'].units  # assign proper units to the corrected 'temperature' array
        Td_new = np.delete(Td.m, bad)
        Td = Td_new * data['dewpoint'].units  # assign proper units to the corrected 'dewpoint' array
        u_new = np.delete(u.m, bad)
        u = u_new * data['u_wind'].units  # assign proper units to the corrected 'u_wind' array
        v_new = np.delete(v.m, bad)
        v = v_new * data['v_wind'].units  # assign proper units to the corrected 'v_wind' array
        hgt_new = np.delete(hgt.m, bad)
        hgt = hgt_new * data['height'].units  # assign proper units to the corrected 'height' array

#### Define display ranges for temperature and pressure
Tmin, Tmax = [-50, 50]
Pmax, Pmin = [1020, 100]

#### Figure out where to plot temperature labels on isotherms for the given T and P display ranges
TT = np.arange(Tmax-10,-121,-10)
LL = np.logspace(np.log10(850), np.log10(Pmin+20), num=TT.size, endpoint=True, base=10.)

#### Compute dewpoint values for desired saturation mixing ratio lines (to plot from 1000 to 500 hPa)
pressure = np.arange(1000, 499, -50) * units('hPa')
mixing_ratio = np.array([1, 2, 4, 6, 8, 10, 13, 16, 20, 25, 30, 36, 42]).reshape(-1, 1) * units('g/kg')
plotTd = mpcalc.dewpoint(mpcalc.vapor_pressure(pressure, mixing_ratio))

#### Compute various parameters (currently set up for 100-hPa mixed layer)
parcel_prof = None
if plotParcelPath is True:
    if surfaceParcel is True:
        parcel_prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
        lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td, parcel_prof, Td[0])
        mixT, mixTd = [T[0], Td[0]]  # for later use
        LI = mpcalc.lifted_index(p, T, parcel_prof)
        CAPE, CIN = mpcalc.cape_cin(p, T, Td, parcel_prof, which_lfc='bottom', which_el='top')
        if mixedParcel is True:
            ## Save above parameters with '_sfc' variable names
            ## (not the cleanest method, but it does the job)
            parcel_prof_sfc = parcel_prof.copy()
            lcl_pressure_sfc, lcl_temperature_sfc = [lcl_pressure.copy(), lcl_temperature.copy()]
            lfc_pressure_sfc, lfc_temperature_sfc = [lfc_pressure.copy(), lfc_temperature.copy()]
            mixT_sfc, mixTd_sfc = [mixT.copy(), mixTd.copy()]
            LI_sfc, CAPE_sfc, CIN_sfc = [LI, CAPE, CIN]
    if mixedParcel is True:
        ## Compute mixed-layer parcel characteristics, then bring down to surface pressure
        ref_T, ref_Td = mpcalc.mixed_layer(p, T, Td, depth=mixedLayerDepth*units('hPa'))
        ref_p = p.m[0] - (float(mixedLayerDepth)/2.)  # midpoint of mixed layer in hPa
        mixT = (ref_T.m * (p.m[0]/ref_p)**(287./1004.)) * units('degC')  # follow dry adiabat to surface P
        q = mpcalc.specific_humidity_from_dewpoint(ref_p*units('hPa'), ref_Td)
        mixTd = mpcalc.dewpoint_from_specific_humidity(p[0], mixT, q)  # compute Td using surface P and above q
        ## Now compute mixed-layer parcel profile and related parameters
        parcel_prof = mpcalc.parcel_profile(p, mixT, mixTd).to('degC')
        lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], mixT, mixTd)
        lfc_pressure, lfc_temperature = mpcalc.lfc(p, T, Td, parcel_prof, mixTd)
        LI = mpcalc.lifted_index(p, T, parcel_prof)
        CAPE, CIN = mpcalc.cape_cin(p, T, Td, parcel_prof, which_lfc='bottom', which_el='top')

# ----------------------------------------------------------------------------------------------------------------------
# EVERYTHING BELOW THIS COMMENT IS FOR PLOTTING PURPOSES
# ----------------------------------------------------------------------------------------------------------------------
#### Create the figure and skew-T diagram (sized for letter paper to be able to print)
fig = plt.figure(figsize=(8.5, 11))
skew = SkewT(fig, rotation=45)
skew.ax.set_xlim(Tmin, Tmax)
skew.ax.set_ylim(Pmax, Pmin)
if plotParcelPath is True:
    if surfaceParcel is True and mixedParcel is True:
        skew.plot(p, parcel_prof, 'k', linewidth=0.75)                      # mixed parcel path is solid
        skew.plot(p, parcel_prof_sfc, 'k', linewidth=0.75, linestyle='--')  # surface parcel path is dashed
    elif surfaceParcel is True and mixedParcel is False:
        skew.plot(p, parcel_prof, 'k', linewidth=0.75, linestyle='--')
    elif surfaceParcel is False and mixedParcel is True:
        skew.plot(p, parcel_prof, 'k', linewidth=0.75)

#### Add data to skew-T plot
skew.plot(p, T, 'r')
skew.plot(p, Td, 'g')
end = np.where(p.m>=100.)[0][-1]  # grab the LAST index value from 'p' with values at or above 100 hPa
skew.plot_barbs(p[:end:4], u[:end:4], v[:end:4], y_clip_radius=0.03)  # plot every 4th wind barb
skew.plot_dry_adiabats(t0=np.arange(233, 533, 10) * units.K, alpha=0.25, colors='orangered', linewidths=1)
skew.plot_moist_adiabats(t0=np.arange(233, 400, 5) * units.K, alpha=0.25, colors='tab:green', linewidths=1)
skew.plot_mixing_lines(pressure=pressure, mixing_ratio=mixing_ratio, linestyles='dotted', colors='dodgerblue', 
                       linewidths=1)
skew.ax.grid(axis='y', c='k', lw=0.5)
for l in range(150,951,100):
    ## This for loop adds dashed horizontal lines for pressure levels at 950 hPa, 850 hPa, etc.
    pp = l * np.ones(10)
    tt = np.linspace(-120,50,10)
    skew.plot(pp, tt, c='k', ls='--', lw=0.5)

#### Add some additional plot features (change tick label sizes, add x- and y-axis labels, add a title)
plt.xticks(fontsize=12)
plt.xlabel('temperature [\N{DEGREE SIGN}C]', fontsize=13)
plt.yticks(fontsize=12)
plt.ylabel('pressure [hPa]', fontsize=13)
plt.title(station+'\n'+dt.strftime('%H UTC %d %b %Y'), loc='left', size=14)

#### Consistently position 'powered by MetPy' label
bbox = skew.ax.get_position()
_ = skew.ax.set_position([bbox.x0, bbox.y0-0.06, bbox.width, bbox.height])  # shift skew-T downward for hodograph
bbox = skew.ax.get_position()  # update 'bbox' values based on previous line
x = bbox.x1
y = bbox.y0 - 0.04
_ = skew.ax.text(x, y, 'powered by MetPy', size=9, weight='bold', ha='right', va='center', transform=fig.transFigure)

#### Add isotherm labels
for i in range(TT.size):
    if TT[i] > 0.:
        color = 'tab:red'
    elif TT[i] == 0.:
        color = '#777777'
    else:
        color = 'tab:blue'
    ## When you see '_ = ' it's a way to hide output to the terminal when working within the Python interpreter
    _ = skew.ax.text(TT[i], LL[i], '%d'%TT[i], ha='center', va='center', weight='bold', color=color, size=8, 
                     rotation=45)

#### Add saturation mixing ratio labels
ll0 = max(skew.ax.get_ylim())
for i in range(mixing_ratio.m.size):
    if plotTd.m[i,0] > Tmax or plotTd.m[i,0] < Tmin:
        continue
    _ = skew.ax.text(plotTd.m[i,0], ll0, '%d'%mixing_ratio.m[i,0], ha='center', va='bottom', weight='bold', size=8, 
                     color='dodgerblue')

#### Identify the location of desired reference altitudes for the input data (1 km, 3 km, 6 km, 12 km, 15 km)
ref_levels = [1., 3., 6., 9., 12., 15.]
h_km = np.array(hgt.m, dtype=float) / 1000.  # convert to float numbers and divide by 1000 to get km
check = np.ediff1d(h_km)
if check.min() <= 0.:
    ## Handle soundings that have bad height values
    rows = np.where(check<=0.)[0]
    h_km = h_km[0:rows[0]]
    p_for_interp = p.m[0:rows[0]]
else:
    p_for_interp = p.m

#### Add height labels found above
for z in ref_levels:
    try:
        ref_p = interp1d(h_km, p_for_interp, kind='linear')(z)
        _ = skew.ax.text(1.06, ref_p, '%d km'%z, transform=skew.ax.get_yaxis_transform(which='tick2'))
    except ValueError:
        continue

#### Print parcel parameter values above skew-T diagram
if plotParcelPath is True:
    if (surfaceParcel is True and mixedParcel is False) or (surfaceParcel is False and mixedParcel is True):
        if surfaceParcel is True:
            title = 'surface-based parcel'
        else:
            title = '%d-hPa mixed layer parcel' % mixedLayerDepth
        #### Add parameter text
        string = '\nT: %.1f\N{DEGREE SIGN}C    Td: %.1f\N{DEGREE SIGN}C'%(mixT.m,mixTd.m)
        if str(np.isnan(lcl_pressure.m)) == 'True':
            string += '    LCL: N/A'
        else:
            string += '    LCL: %.1f hPa' % lcl_pressure.m
        if str(np.isnan(lfc_pressure.m)) == 'True':
            string += '    LFC: N/A'
        else:
            string += '    LFC: %.1f hPa' % lfc_pressure.m
        string += '\nCAPE: %.1f J/kg    CIN: %.1f J/kg    LI: %.1f' % (CAPE.m,CIN.m,LI.m)
        x = bbox.x0
        y = bbox.y1 - 0.1
        _ = skew.ax.text(x, y, title, ha='left', va='top', size=9, weight='bold', 
                         transform=fig.transFigure)
        _ = skew.ax.text(x, y, string, ha='left', va='top', size=9, transform=fig.transFigure)
    elif surfaceParcel is True and mixedParcel is True:
        ## surface parameters first
        sfcString = '\nT: %.1f\N{DEGREE SIGN}C    Td: %.1f\N{DEGREE SIGN}C'%(mixT_sfc.m,mixTd_sfc.m)
        if str(np.isnan(lcl_pressure_sfc.m)) == 'True':
            sfcString += '    LCL: N/A'
        else:
            sfcString += '    LCL: %.1f hPa' % lcl_pressure_sfc.m
        if str(np.isnan(lfc_pressure_sfc.m)) == 'True':
            sfcString += '    LFC: N/A'
        else:
            sfcString += '    LFC: %.1f hPa' % lfc_pressure_sfc.m
        sfcString += '\nCAPE: %.1f J/kg    CIN: %.1f J/kg    LI: %.1f' % (CAPE_sfc.m,CIN_sfc.m,LI_sfc.m)
        print('surface-based parcel' + sfcString)
        x = bbox.x0
        y = bbox.y1 + 0.19
        _ = skew.ax.text(x, y, 'surface-based parcel', ha='left', va='top', size=9, weight='bold', 
                         transform=fig.transFigure)
        _ = skew.ax.text(x, y, sfcString, ha='left', va='top', size=9, transform=fig.transFigure)
        ## mixed-layer parameters next
        title = '%d-hPa mixed layer parcel' % mixedLayerDepth
        string = '\nT: %.1f\N{DEGREE SIGN}C    Td: %.1f\N{DEGREE SIGN}C'%(mixT.m,mixTd.m)
        if str(np.isnan(lcl_pressure.m)) == 'True':
            string += '    LCL: N/A'
        else:
            string += '    LCL: %.1f hPa' % lcl_pressure.m
        if str(np.isnan(lfc_pressure.m)) == 'True':
            string += '    LFC: N/A'
        else:
            string += '    LFC: %.1f hPa' % lfc_pressure.m
        string += '\nCAPE: %.1f J/kg    CIN: %.1f J/kg    LI: %.1f' % (CAPE.m,CIN.m,LI.m)
        print(title + string)
        x = bbox.x0
        y = bbox.y1 + 0.1 + 0.045
        _ = skew.ax.text(x, y, title, ha='left', va='top', size=9, weight='bold', 
                         transform=fig.transFigure)
        _ = skew.ax.text(x, y, string, ha='left', va='top', size=9, transform=fig.transFigure)

#### Add hodograph panel
stop = max(np.where(hgt.m<=10000.)[0]) + 1  # display hodograph below 10 km
ax2 = fig.add_axes([bbox.x1-0.22, bbox.y1, 0.22, 0.22])
h = Hodograph(ax2, component_range=80.)
h.add_grid(increment=20)
h.ax.set_xlim(-80,80)
h.ax.set_xticks(np.arange(-75,76,25))
h.ax.set_ylim(-80,80)
h.ax.set_yticks(np.arange(-75,76,25))
h.ax.tick_params(labelsize=8.5)
h.plot_colormapped(u[0:stop], v[0:stop], hgt[0:stop])

#### Save and display produced skew-T diagram
#plt.savefig(station+dt.strftime('_%Y%m%d%H.png'), bbox_inches='tight', pad_inches=0.05)
#plt.savefig(station+dt.strftime('_%Y%m%d%H.pdf'))
plt.show()  # add a pound sign before this line if you don't want Python to display the plot
plt.close()




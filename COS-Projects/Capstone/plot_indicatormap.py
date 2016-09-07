import numpy as np
import csv as csv
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

#from geonamescache import GeonamesCache
#from helpers import slug
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from mpl_toolkits.basemap import Basemap
from readcsv import *
from data_features import *

datadir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/PERM'
plotdir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Figures/'
inddir = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/Indicators/Country/'
shapefile = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/GEO/ne_10m_admin_0_countries'
isofile = '/Users/sudhirraskutti/Desktop/COS424/Homework/Project/Data/GEO/wikipedia-iso-country-codes.csv'

year = '2008'
years, countries, indicators, inddata = read_wb_indicator(inddir, 'WB_Development')
countries = np.asarray([c.upper() for c in countries])
isodata = read_any(isofile)
selindc = 'GINI index (World Bank estimate)'
selind = np.squeeze(np.where(indicators == selindc))
yind = np.squeeze(np.where(years == year))
indvals = np.squeeze(inddata[yind,:,selind])
indcountries = countries

countries = isodata['English short name lower case']
countries = np.asarray([c.upper() for c in countries])
isos = np.asarray(isodata['Alpha-3 code'])
ckey = 'COUNTRY_OF_CITZENSHIP'
alldata = read_cfs(datadir, ['PERM_FY' + year])
cleandata = clean_data(alldata, 1)
countryaccs = acc_rate(cleandata, ckey, 40)
countrykey = countryaccs.keys()

countryaccs = {}
for i in xrange(0, len(indcountries)):
    k = indcountries[i]
    v = indvals[i]
    countryaccs[k] = v
cvals = countryaccs.values()
altcountrykey = countryaccs.keys()
countrykey = list(set(countrykey) & set(altcountrykey))

num_colors = 12
cm = plt.get_cmap('Reds')
scheme = [cm(float(i) / num_colors) for i in range(num_colors)]
bins = np.linspace(-1.0, 1.0, num_colors)
print bins

if (True):
    #mpl.style.use('map')

    ysize = 0.25 * 11.69
    xsize = 0.4 * 8.27
    fontsize = '10'
    fig = plt.figure(figsize = (22,12))
    fig.suptitle('{1:s} in {0:s}'.format(year, selindc), fontsize=25, y=.95)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif', size=fontsize)

    ax = plt.subplot(1,1,1, axisbg='w', frame_on=False)
    m = Basemap(lon_0=0, projection='robin')
    m.drawmapboundary(color='w')
    
    m.readshapefile(shapefile, 'units', color='#444444', linewidth=.2)
    for info, shape in zip(m.units_info, m.units):
        iso3 = info['ADM0_A3']
        niso = np.where(isos == iso3)
        country = countries[niso]
        nlc = len(country)
        if nlc > 0:
            country = country[0]
        if country in countrykey and nlc > 0:
            ar = countryaccs[country]
            arbin = np.digitize(ar, bins) - 1
            color = scheme[arbin]
        else:
            color = 'w'
        patches = [Polygon(np.array(shape), True)]
        pc = PatchCollection(patches)
        pc.set_facecolor(color)
        ax.add_collection(pc)

# Cover up Antarctica so legend can be placed over it.
    ax.axhspan(0, 1000 * 1800, facecolor='w', edgecolor='w', zorder=2)

# Draw color legend.
    ax_legend = fig.add_axes([0.35, 0.14, 0.3, 0.03], zorder=3)
    cmap = mpl.colors.ListedColormap(scheme)
    cb = mpl.colorbar.ColorbarBase(ax_legend, cmap=cmap, ticks=bins, boundaries=bins, orientation='horizontal')
    cb.ax.set_xticklabels([' ' for i in bins])

    plt.savefig(plotdir + '/WorldMap' + year + selindc + '.png', bbox_inches='tight', pad_inches=.2)
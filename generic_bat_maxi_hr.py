import numpy as np
from astropy.io import fits, ascii
from astropy.time import Time

import matplotlib.pyplot as plt
from matplotlib import rc
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)
rc('font', family='serif')
colorset = ['#000000','#00270C','#00443C','#005083','#034BCA','#483CFC','#9C2BFF','#EB24F4','#FF2DC2','#FF4986','#FF7356','#FFA443','#EBD155','#D3F187','#D7FFC8','#FFFFFF']

def lc_binning(t,y,dy,width):
    """
    lc_binning(t,y,dy,width)

    A simple function to bin lightcurves with uncertainties.
    
    Parameters
    ----------

    t : array_like
        Array containing time values
    
    y : array_like
        Array containing flux/rate/magnitude values
    
    dy : array_like
         Array containing uncertainties values (not weights) for y
    
    width : number
            Width of each bin in the same units as t
    
    Returns
    -------
    binned_t : array
               binned time array
    
    binned_y : array
               binned flux/rate/magnitude array
    
    binned_dy : array
                binned flux/rate/magnitude uncertainty array
    
    """
    width = float(width)
    nbins = int(np.ceil((t[-1]-t[0])/width))+1
    binned_y = np.zeros(nbins)
    binned_y_weights = np.zeros(nbins)
    
    for i in range(len(y)):
        bin_no = int((t[i]-t[0])/width)
        weight = dy[i]**-2.0
        binned_y[bin_no] += y[i] * weight
        binned_y_weights[bin_no] += weight
        
    binned_y /= binned_y_weights
    binned_dy = np.sqrt(1.0/(binned_y_weights)) 
    binned_t = np.arange(start=t[0]+(width/2.0),stop=t[0]+nbins*width,step=width)
    return binned_t, binned_y, binned_dy


# Inputs:
src_name = 'NGC 6624'
bat_data_path = 'https://swift.gsfc.nasa.gov/results/transients/H1820-303.lc.fits'
maxi_data_path = 'http://maxi.riken.jp/star_data/J1823-303/J1823-303_g_lc_1day_all.dat'
#comments = 'Persistent source'
# if you dont want the above line, just comment it out

# LC properties
x_init = 57950
x_end = Time.now().mjd+10
binwidth = 3.0

# Reading and  binning the light curves:
bat_data = fits.open(bat_data_path,cache=False)[1].data
maxi_data = ascii.read(maxi_data_path,names=('MJD','RATE_2_20','RATE_2_20_ERR','RATE_2_4','RATE_2_4_ERR','RATE_4_10','RATE_4_10_ERR','RATE_10_20','RATE_10_20_ERR'))
bat_binned_t, bat_binned_y, bat_binned_dy = lc_binning(bat_data['TIME'],bat_data['RATE'],bat_data['ERROR'],width=binwidth)
maxi_binned_t, maxi_binned_y, maxi_binned_dy = lc_binning(maxi_data['MJD'],maxi_data['RATE_4_10'],maxi_data['RATE_4_10_ERR'], width=binwidth)

# Estimating hardness ratio:
hardness_t_maxi = []
hardness_t_bat = []
hardness = []
hardness_err = []
for i in range(len(maxi_binned_t)):
    for j in range(len(bat_binned_t)):
        if abs(maxi_binned_t[i]+0.5-bat_binned_t[j]) < binwidth/2:
            hardness_t_maxi.append(maxi_binned_t[i])
            hardness_t_bat.append(bat_binned_t[j])
            hardness.append(((bat_binned_y[j]/0.220)-(maxi_binned_y[i]/1.15))/((bat_binned_y[j]/0.220)+(maxi_binned_y[i]/1.15)))
            
            xdif = (bat_binned_y[j]/0.220)-(maxi_binned_y[i]/1.15)
            xsum = (bat_binned_y[j]/0.220)+(maxi_binned_y[i]/1.15)
            xdiff_er = np.sqrt( (bat_binned_dy[j]/0.220)**2 + (maxi_binned_dy[i]/1.15)**2 )
            xsum_er = np.sqrt( (bat_binned_dy[j]/0.220)**2 + (maxi_binned_dy[i]/1.15)**2 )
            hardness_err.append(abs(hardness[-1])*np.sqrt( (xdiff_er/xdif)**2 + (xsum_er/xsum)**2 ))

"""
Note on converting rates to crab flux:

BAT rates are converted to Crab assuming 1 mCrab = 0.000220 ct/cm2/sec 
Reference for conversion: https://swift.gsfc.nasa.gov/results/transients/

MAXI rates are converted to Crab assuming 1 Crab = 1.15 ct/s/cm2 (in the 4-10 keV band)
Reference for conversion: http://maxi.riken.jp/top/readme.html

"""

# Plotting
plt.figure(figsize=(10,7))
plt.subplot2grid((3,1),(0,0))
plt.errorbar(bat_binned_t, bat_binned_y/0.220, bat_binned_dy/0.220,fmt='o',capsize=3,ms=2,alpha=0.8, markeredgecolor=colorset[0],color=colorset[4],markeredgewidth=0.5)
plt.ylabel('BAT Rate\n15-50 keV\n(Crab)',fontsize=12)
plt.xlim(x_init,x_end)
plt.ylim(0)
#plt.yscale('log')
plt.gca().axes.set_xticklabels([])
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='major', length=9)
plt.tick_params(axis='both', which='minor', length=4.5)
plt.tick_params(axis='both', which='both',direction='in',right='on',top='on')
plt2 = plt.twiny()
plt2.set_xticks(np.arange(0,1,1.0/5.001))
plt2.set_xticklabels(Time(Time((plt2.get_xticks()*(x_end-x_init))+x_init,format='mjd').iso,out_subfmt='date').value,rotation=45,ha='left')
plt2.tick_params(axis='both', which='major', labelsize=13)
plt2.tick_params(axis='both', which='major', length=9)
plt2.tick_params(axis='both', which='minor', length=4.5)
plt2.tick_params(axis='x', which='both',direction='in')
plt2.minorticks_on()

plt.subplot2grid((3,1),(1,0))
plt.errorbar(maxi_binned_t, maxi_binned_y/1.15, maxi_binned_dy/1.15,fmt='o',capsize=3,ms=2,alpha=0.8, markeredgecolor=colorset[0],color=colorset[6],markeredgewidth=0.5)
plt.ylabel('MAXI Rate\n4-10 keV\n(Crab)',fontsize=12)
plt.xlim(x_init,x_end)
plt.ylim(0)
#plt.yscale('log')
plt.gca().axes.set_xticklabels([])
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='major', length=9)
plt.tick_params(axis='both', which='minor', length=4.5)
plt.tick_params(axis='both', which='both',direction='in',right='on',top='on')

plt.subplot2grid((3,1),(2,0))
plt.errorbar(hardness_t_maxi,hardness,hardness_err,fmt='o',capsize=3,ms=2,alpha=0.8, markeredgecolor=colorset[0],color=colorset[9],markeredgewidth=0.5)
plt.ylabel('Hardness Ratio\n'+r'($\frac{BAT-MAXI}{BAT+MAXI}$)',fontsize=12)
plt.xlim(x_init,x_end)
plt.ylim(-0.99,0.99)
plt.minorticks_on()
plt.tick_params(axis='both', which='major', labelsize=12)
plt.tick_params(axis='both', which='major', length=9)
plt.tick_params(axis='both', which='minor', length=4.5)
plt.tick_params(axis='both', which='both',direction='in',right='on',top='on')

plt.xlabel('MJD',fontsize=14)
try:
    plt.text(x_init,-2.0,src_name+'\nLight curves binned by '+str(binwidth)+' days\nLast Update: '+Time.now().isot+' UTC\n'+comments,ha='left')
except:
    plt.text(x_init,-2.0,src_name+'\nLight curves binned by '+str(binwidth)+' days\nLast Update: '+Time.now().isot+' UTC',ha='left')
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.08, hspace=0.0)
plt.savefig(src_name.replace(' ','')+'_bat_maxi_lc.pdf',bbox_inches='tight')

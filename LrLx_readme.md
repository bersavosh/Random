# XRB-LrLx
A python script for plotting radio vs. X-ray luminosity for X-ray binaries with quick and easy modification. I use python dictionaries to classify input and make it easy to separate them based on source class, GC association, etc. Additionally this script produces a table of sources currently present in its database. 

The sources currently available in the database are listed in [source_catalog.html](http://htmlpreview.github.io/?https://github.com/bersavosh/Random/blob/master/source_catalog.html).

**Python libraries used:** [Numpy](http://www.numpy.org/), [Astropy](http://www.astropy.org/), [Matplotlib](http://matplotlib.org/)

## Table of Contents:

[Using the script](https://github.com/bersavosh/XRB-LrLx/blob/master/README.md#using-the-script)

[Adding new data to the plot](https://github.com/bersavosh/XRB-LrLx/blob/master/README.md#adding-new-data-to-the-plot)

[Type of object](https://github.com/bersavosh/XRB-LrLx/blob/master/README.md#type-of-object)

[Handling upper limits](https://github.com/bersavosh/XRB-LrLx/blob/master/README.md#handling-upper-limits)

[Adding a new type of object](https://github.com/bersavosh/XRB-LrLx/blob/master/README.md#adding-a-new-type-of-object)

[Source Catalog](https://github.com/bersavosh/XRB-LrLx/blob/master/README.md#source-catalog)

[Help and accessibility](https://github.com/bersavosh/XRB-LrLx/blob/master/README.md#help-and-accessibility)

## Using the script:
This is a simple Python script. You need to download the `Lfile` directory which contains the luminosity database to be able to run it. By default the script saves the final plot as a pdf file in the same directory as the script.

## Adding new data to the plot:
Depending on how you have stored the luminosity data, there are two ways to add a new source to the plot: You can either add a source with luminosity values directly in the command (easy for sources with just one data point or checking the location of a new source) or direct the script to read luminosity values from an ascii file (appropriate for sources with multiple data points, and/or data points with asymmetric errorbars).

For either of these options, you can add one of the following lines to [Part II of `lrlx_plot.py`](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L99):

### 1: Reading luminosity values from a file:
If reading luminosities from an ascii file, the file should be stored in the "Lfiles" folder and should contain at least two columns with labels `Lx` and `Lr`. If there are columns for uncertainties in the data, they should be labeled `Lx_er` and `Lr_er`. It is also possible to add asymmetric uncertainties. If so, columns should be named `Lx_ler` and `Lx_uer` (`Lr_ler` and `Lr_uer` for radio). You can check the current files in the "Lfiles" directory for examples.

```
src_list.append(src_dict(<source name>, <source class>, <globular cluster>,\
                Lfile=<filename.txt>, uplim=<the axis with upper limit values>, ref=<reference>))
```

### 2: Reading luminosity values in-line:
If reading luminosity values in-line, then uncertainties should follow a format recognizable by [`pyplot.errorbar`](http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.errorbar).
```
src_list.append(src_dict(<source name>, <source class>, <globular cluster>,\
                Lx=[lx1,lx2,...], Lr=[lr1,lr2,...],\
                uplim=<the axis with upper limit values>, ref=<reference>))
```

## Type of object:
Currently, the following classes are defined and their members from the source list in the script are automatically plotted:

* `BH`: Quiescent and hard-state black hole X-ray binaries (both confirmed and candidates) from the literature.

* `LrLx_BH`: BH candidates identified based on LR-Lx correlation.

* `tMSP`: transitional millisecond pulsars.

* `NS`: Hard-state neutron star X-ray binaries.

* `AMXP`: Accreting millisecond X-ray pulsars.

* `CV`: Cataclysmic variables.

* `UI`: Unidentified.

## Handling upper limits:
This script is able to handle plotting upper limits (either on Lx or LR) as well. You can check the lines in the code for [M62](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L131), [M22](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L132) and [VLA J2130+12](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L134) as examples for radio sources with no X-ray detection and the line for [EXO 1745-248](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L143) for an X-ray source with no radio detection. If all the data points for your source are just upperlimits (similar to M62 and M22), you only need one entry for `src_list`. However, if some data points are detected values and some are upper limits, then you need two different entries; one for the detections and one for upper limits. See the case for [EXO 1745-248](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L143) as an example.


## Adding a new type of object:

1- Write the new class while adding the source to the catalog

2- Add the new class to the [plot handler](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L196). 
  
Example:
  
```
if i['Class'] == 'newclass':
    newclass,=plt.loglog(i['Lx'],i['Lr'],'o',ms=4, c='k',mec='k',zorder=2,label='New class')

```
  
3- Add the artist object `newclass` defined above to the [legend handler](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L223).
  
Example:
  
```
plt.legend(handles=[BHs,LrLxBHs,NSs,AMXPs,tMSPs,CVs,newclass],loc=2,numpoints=1,fontsize=12)
```


The script [`lrlx_plot.py`](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py) contains description and example of adding new data and modifying previous data as well. 

## Source Catalog:
Besides producing a plot, the script also creates a searchable, sortable html table of sources currently present in the database. This table is stored in [source_catalog.html](http://htmlpreview.github.io/?https://github.com/bersavosh/Random/blob/master/source_catalog.html) and is updated everytime the script is executed. Note that this is performed via [`astropy.table`](http://docs.astropy.org/en/stable/table/) implemented in [Part II of the script](https://github.com/bersavosh/XRB-LrLx/blob/master/lrlx_plot.py#L158). Thus you can change the table format to your desired format (e.g., latex, ascii).

## Help and accessibility: 

This is a private repository and is only accessible by invited collaborators. I've prepared this script to ease plotting and cataloging sources studied by the MAVERICS collaboration. 

I'm still occasionally improving this script and the database. If you have questions/comments/suggestions/updates, please don't hesitate to contact me (either through available methods on Github or by email). 

<img src="https://github.com/bersavosh/Random/blob/master/lrlx_plot.jpg" width="600">

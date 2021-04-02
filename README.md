# stormviz
Storm Visualization

### Overview
The Storm Visualization Tool (stormviz) is a graphical tool to help see the variable paths of cyclonic storm tracks relative to atmospheric pressure fields and how these tracks evolve with time and space.  Combined with data layers of geopotential height and pressure at mean sea level, the tool helps provide understanding of the basic atmospheric dynamics and structure associated with storm generation, duration, and termination. 

The tool provides a map to view tracks with surface pressure highs and lows.  [ECMWF Reanalysis (ERA5)](https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era5) [1] data of height of the geopotential surface (hgt), and mean sea-level pressure (msl) fields are displayed.  Currently, the extents of the map and data focus on North America and hence the storms in the northern hemisphere.  These extents can be adjusted in code and data downloaded from ERA5 for other regions of the world.

### The Graphical Tool

The main map shows the current storm positions (red dots) for the current date and time and their associated "detected" tracks (discussed below) in blue lines with blue dots for each time step beyhond the present. Also plotted are fields of msl with filled contours corresponding to the scale below the map and hgt (100-1000 hPa pressure level) for North America. By using the graphical interface, times can be selected to show (or "play") the 3-dimensional (3D), dynamic nature of the storm tracks along with changes in msl. Furthermore, the hgt levels can be selected to view the upper atmospheric structure and how far up in the atmosphere the storm may affect this structure. 

Here is a screen image of the map and interface. 

![Image of stormviz map](https://github.com/neaptide/stormviz/blob/main/images/stormviz_gui_map.png)

### Storm Detection and Tracking

Storm tracks are determined and stored using the [stormTracking](https://github.com/ecjoliver/stormTracking) [2] algorithm by Eric Oliver.  This code automatically detects and tracks atmospheric storms (cyclones) and high-pressure centers (anticyclones) using a series of surface pressure maps. While the initial code used NCEP data, it is easily adapted to use ERA5 data and resolution.  A sample storm track dataset (`storm_track_slp.npz`) is provided for testing `stormviz` following the steps to run the tool. 

### Get the code 

You can clone the [stormviz Github Repository](https://github.com/neaptide/stormviz) or download the zipped code from Github: 

```bash
git clone https://github.com/neaptide/stormviz.git
```

or download and unzip the code:

```bash
wget https://github.com/neaptide/stormviz/archive/main.zip
```
### Build environment 

Once the code is downloaded locally, you will be using the `environment.yml` file in the code folder to get all the dependencies to build an environment with default name `test_stormviz`.  You can use a different enviroment name by editing the `environment.yml`. 

For this example the code is unzipped or cloned into `C:\your\home\stormviz`
 
Open a terminal that activates Anaconda base environment. 
Change the working directory to where you unpacked the code and create the `test_stormviz` environment  

```
(base) C:\your\home> cd stormviz
(base) C:\your\home\stormviz>conda env create -f environment.yml
```

Activate the `test_stormviz` environment
 
```
(base) C:\your\home\stormviz>conda activate test_stormviz
(test_stormviz) C:\your\home\stormviz>
```
 
You can now run `stormviz.py` in either IPython or in a Jupyter Notebook. 

### Running `stormviz.py` in Ipython

Open an Ipython terminal.
```
(test_stormviz) C:\your\home\stormviz>ipython
```

In Ipython, use the magic command `%run` to run the `stormviz.py` code for the specified year and month (YYYY_MM). 

``` 
[1] %run stormviz.py 2018_01
```

This will initialize the interactive graph, previously described, to the first day and hour of that month on the map and you can begin using the interface. For example, if `2018_01` is used, the map will show the data for 2018-01-01 at 00:00 (UTC). The hgt data will initialize to the the lowest level (1000 hPa).  

### Using the interactive display

- Select another hgt level by moving the "Level" slider or pressing left- (<) and right-arrow (>) associated with it.  
- Select a different time and date by moving the "Date" slider or pressing left- (<) and right-arrow (>) associated with it.
  - Press "Storm Tracks ON/OFF" button to toggle display of storm tracks.
  - Press "MSL ON/OFF" button to toggle display of mean sea level pressure.
  - Press "HGT ON/OFF" button to toggle display of geopotential contours.


### Dependencies

`stormviz` has dependency on the following Python modules:

  - numpy
  - matplotlib
  - metpy
  - netCDF4
 
 ### Acknowledgements

The `stormviz` tool is based upon work supported by National Science Foundation under grant OCE-1558920. 

### References

[1] Copernicus Climate Change Service (C3S) (2017): ERA5: Fifth generation of ECMWF atmospheric reanalyses of the global climate . Copernicus Climate Change Service Climate Data Store (CDS), date of access. https://cds.climate.copernicus.eu/cdsapp#!/home

[2] Oliver, E.: stormTracking, https://github.com/ecjoliver/stormTracking, @ecjoliver

# Godunov-type finite volume scheme for 2d, free-surface, inviscid shallow water equations on a Cartesian grid

Date: Jan. 31, 2018

Dr.-Ing. Roberto Tatis-Muvdi

Institute of Hydraulic Engineering and Technical Hydromechanics

Technische Universität  Dresden

## General information:

This code implements a combination of the approaches proposed by Liang and Borthwick (2009) and Liang (2010) for solving the inviscid 2d shallow water equations with wetting and drying on a Cartesian grid. Face fluxes are computed using an HLLC approximate Riemann solver. The model uses ASCII raster files as inputs for bottom elevation (z), Manning roughness (n) and initial conditions (water level (eta), x-velocity (u), y-velocity (v) (last two optional)), and simple text files for the control file and boundary conditions. Output files are also in ASCII raster format.

The code makes use of the [Diffpack library](http://diffpack.de/), mostly concerning array definitions and indexing.

- Coordinate system: Cartesian, x increases to the right (E), y increases upwards (N). Same for velocities.
- Unit system: SI (meter, second, kg)
- Comments are not yet allowed in any of the input files.
- All input files must be in the same folder. All file names given by the user must have extension.
- The only fixed file name is that of the control file “control.txt”
- The model can be started through a batch file with the following content:
  swRiemann.exe
  pause

## Demo - Dam break simulation over irregular terrain for testing wettin-drying scheme
![](3hump.gif)

## Variables used in the model:
- eta: Water level \[m\]
- Q: discharge \[m<sup>3</sup>/s\]
- u: y-velocity component	\[m/s\]
- v: x-velocity component	\[m/s\]
- qx: x-flux \[m<sup>2</sup>/s\]
- qy: y-flux \[m<sup>2</sup>/s\]

## Input files:

### Raster files
These are ASCII text files containing the input arrays for the model. File names and extensions are user-defined. The model requires:
  - Terrain
  - Manning
  - Initial water level (for dry start, terrain file can be used)
  - Initial x-velocity component (optional)
  - Initial y-velocity component (optional)

![File format](https://github.com/robetatis/sweRiemann/blob/master/fileFormat.png)

  File header:
    - NCOLS: number of columns
    - NROWS: number of rows
    - XLLCORNER: x coordinate of lower-left corner
    - YLLCORNER: y coordinate of lower-left corner
    - CELLSIZE: pixel size
    - NODATA_VALUE: flag for pixels with non-existent data

*Note*: All the above files must have the same cell size, number of columns and number of rows. In its current version, the model does not check file formats, for instance, whether the information in the header corresponds to the actual size of the data matrix, or if the terrain, initial conditions and manning rasters have different number of columns, rows, cell size, etc.

### Text files

These are ASCII text files describing the boundary conditions and the configuration of the model (control file, see below):

- Boundary conditions: the model works with *open* (transmissive), *closed* (reflective), *discharge* or *water level* boundary conditions. 
Since the model is raster-based, there are four domain boundaries: north (N), south (S), east (E) and west (W). Along each domain boundary there can be as many boundary segments of different types as required. For instance:

  ![Boundary segments](https://github.com/robetatis/sweRiemann/blob/master/bc.png). 
  Blue: water level boundary segments, green: discharge boundary segments, black: closed boundary segments.

  By default, all model boundaries are closed and only those segments indicated in the boundary conditions file will be either open, have a discharge or a water level. Information must be provided for open, discharge or water level boundaries. If nothing is entered, the corresponding model boundary is treated as closed. The keywords needed for this are (last column of boundary conditions file):
  - Q: flow boundary
  - eta: water level boundary
  - open: open (transmissive) boundary

  The boundary conditions text file must have the following format:

  ![Boundary condition file format](https://github.com/robetatis/sweRiemann/blob/master/bcFile.png)

  Each line in this file corresponds to a boundary segment. The first line must contain column names. This is simply a header for identfying columns. Columns are read and interpreted based on position, not name!

  Columns:

  - edge: can be N (north), S (south), E (east) or W (west). Indicates the location of the boundary segment.
   start: the starting coordinate of the boundary segment. If this point does not coincide exactly with the position of a cell edge, the model automatically shifts it upwards or to the right (depending on whether the boundary segment is along a vertical or a horizontal domain edge).
  - end: the end coordinate of the desired boundary. If this point does not coincide exactly with the position of a cell edge, the model automatically shifts it upwards or to the right (depending on whether the boundary segment is along a vertical or a horizontal domain edge).
  - btype: boundary type. This simply indicates what is the variable contained in the boundary condition. As indicated in the keywords above, can be only “open”, “Q” or “eta”. These keywords are case sensitive! 
  - bvaluefile: this is the name of the file containing the discharge or water level time series for the corresponding boundary segment, i.e., the time series file. The model only works with boundary conditions in the form of time series. For steady inflow/outflow or water level, simply provide a flat time series. These file names and extensions are user-defined. In the example above, line 5 of the boundary conditions file:

  S		45.104 		100.23		eta		eta.txt     

  points to file “eta.txt”. This file must also be in the same folder. Make sure to name the file exactly as indicated in the boundary conditions file. 

  The format for the time series file is very simple:

  ![Time series file format](https://github.com/robetatis/sweRiemann/blob/master/bcTimeSeries.png)

The 1st line is simply a header, not used for identification. The 1st column contains the time in seconds. The 2nd column contains the value of the variable. The time series is interpolated to the model’s computational time steps automatically, such that only inflection points in the curve must be provided. There must be one time series file per boundary segment.

- Control file: this file contains all run control parameters. The ordering of the lines in this file is fixed. Keywords in this file are only a guide for the user, as they are not really used for the assignment of parameter values.
  - tstop: simulation duration, in seconds
  - dt: time step, in seconds
  - Ctolup: maximum tolerable Courant number
  - Ctollow: minimum tolerable Courant number
  - outprf: ouput frequency of results files, in seconds
  - outpc: output frequency for console, in seconds
  - demfile: name of terrain ASCII file
  - bcfile: name of boundary conditions file
  - etaicfile: name of initial water level ASCII file. It can be the same as the terrain for a dry start. Cannot be “no”.
  - uicfile: name of ASCII file containing initial x-velocity component. This is optional, i.e., if not provided, all initial x-velocities are set to zero. In this case, “no” must be entered (case sensitive)
  - vicfile: name of ASCII file containing y-velocity component. This is optional, i.e., if not provided, all initial u-velocities are set to zero. In this case, “no” must be entered (case sensitive)
  - manfile: name of Manning ASCII file
  - rfpre: results file prefix. This is simply a root for the name of result files.
  - outdir: name of results directory. All output files will be written here
  - hout: “yes” or “no” for controlling output of water depth
  - etaout: “yes” or “no” for controlling output of water level
  - vout: “yes” or “no” for controlling output of v-velocity component
  - uout: “yes” or “no” for controlling output of u-velocity component
  - vmagout: “yes” or “no” for controlling output of velocity magnitude, √(u2 + v2)
  - flowdout: “yes” or “no” for controlling output of flow direction in degrees

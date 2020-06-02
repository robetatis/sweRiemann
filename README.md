# Godunov-type finite volume scheme for two-dimensional inviscid shallow water equations on a Cartesian grid

Date: Jan. 31, 2018

Dr.-Ing. Roberto Tatis-Muvdi

Institute of Hydraulic Engineering and Technical Hydromechanics

Technische Universität  Dresden

## General information:

This code implements a combination of the approaches proposed by Liang and Borthwick (2009) and Liang (2010) for solving the inviscid 2d shallow water equations with wetting and drying on a Cartesian grid. It uses raster files as input for bottom elevation (z) and initial conditions, and simple text files for the control file and boundary conditions. Outputs are also in raster format.

The code makes use of the [Diffpack library](http://diffpack.de/), mostly concerning array definitions and indexing.

- Coordinate system: Cartesian, x increases to the right (E), y increases upwards (N). Same for velocities.
- Unit system: SI (meter, second, kg)
- Comments are not yet allowed in any of the input files.
- All input files must be in the same folder. All file names given by the user must have extension.
- The only fixed file name is that of the control file  “control.txt”
- The model can be started through a batch file with the following content:
  swRiemann.exe
  pause

## Variables used in the model:
- eta: Water level \[m\]
- Q: discharge \[m<sup>3</sup>/s\]
- u: y-velocity component	\[m/s\]
- v: x-velocity component	\[m/s\]

## Input files:

### ASCII text files. These file names and extensions are user-defined.
  - Terrain
  - Manning
  - Initial water level (for dry start, terrain file can be used)
  - Initial x-velocity component (optional)
  - Initial y-velocity component (optional)

![File format]:(https://github.com/robetatis/sweRiemann/blob/master/fileFormat.png)

  File header:
    - NCOLS: number of columns
    - NROWS: number of rows
    - XLLCORNER: x coordinate of lower-left corner
    - YLLCORNER: y coordinate of lower-left corner
    - CELLSIZE: pixel size
    - NODATA_VALUE: flag for pixels with non-existent data

*Note*: All the above files must have the same cell size, number of columns and number of rows. In its current version, the model does not check file formats, for instance, whether the information in the header corresponds to the actual size of the data matrix, or if the terrain, initial conditions and manning rasters have different number of columns, rows, cell size, etc.

### Text files
- Boundary conditions: the model works with *open* (transmissive), *closed* (reflective), *discharge* or *water level* boundary conditions. 
Since the model is raster-based, there are four domain boundaries: north (N), south (S), east (E) and west (W). Along each domain boundary there can be as many boundary segments of different types as required. For instance:

![Boundary segments]:(https://github.com/robetatis/sweRiemann/blob/master/bc.png)
Blue: water level boundary segments, green: discharge boundary segments, black: closed boundary segments.

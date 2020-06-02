# Godunov-type finite volume scheme for two-dimensional inviscid shallow water equations on a Cartesian grid

#### Date: Jan. 31, 2018
#### Dr.-Ing. Roberto Tatis-Muvdi
#### Institute of Hydraulic Engineering and Technical Hydromechanics
#### Technische Universität  Dresden

This code implements a combination of the approaches proposed by Liang and Borthwick (2009) and Liang (2010) for solving the inviscid 2d shallow water equations with wetting and drying on a Cartesian grid. It uses raster files as input for bottom elevation (z) and initial conditions, and simple text files for the control file and boundary conditions. Outputs are also in raster format.

The code makes heavy use of the http://diffpack.de/[Diffpack library]

## General information:
- Coordinate system: Cartesian, x increases to the right (E), y increases upwards (N). Same for velocities.
- Unit system: SI (meter, second, kg)
- Comments are not yet allowed in any of the input files.
- All input files must be in the same folder. All file names given by the user must have extension.
- The only fixed file name is that of the control file  “control.txt”
- The model can be started through a batch file with the following content:
  swRiemann.exe
  pause

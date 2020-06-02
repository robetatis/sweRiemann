# Godunov-type finite volume scheme for two-dimensional inviscid shallow water equations on a Cartesian grid

## Date: 31.1.2018
## Dr.-Ing. Roberto Tatis-Muvdi
## Institut für Wasserbau und THM
## Technische Universität  Dresden

This code implements a combination of the approaches proposed by Liang and Borthwick (2009) and Liang (2010) for solving the inviscid 2d shallow water equations with wetting and drying on a Cartesian grid. It uses classical raster file inputs for bottom elevation (z) and initial conditions, and simple text files for the control file and boundary conditions.

## General information:
- Coordinate system: Cartesian, x increases to the right (E), y increases upwards (N). Same for velocities.
- Unit system: SI (meter, second, kg)
- Comments are not yet allowed in any of the input files.
- All input files must be in the same folder. All file names given by the user must have extension.
- The only fixed file name is that of the control file  “control.txt”
- The model can be started through a batch file with the following content:
  swRiemann.exe
  pause

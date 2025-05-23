# WFTEM3D

A quick, flexible, and extensible 3D TEM modeling open-source software.

**WFTEM3D1.0:**
Available in Matlab/Fortran/EXE versions. 
Supports loop source, and ground/airborne/marine/tunnel/borehole scenarios.

**WFTEM3D2.0:**
Available in Matlab/Fortran/EXE versions.
Supports wire/loop sources, and ground/airborne/marine/tunnel/borehole scenarios.

WFTEM3D2.0 is simpler and more applicable than WFTEM3D1.0.

## Method

### WFTEM3D1.0:
First the scheme calculates the approximate initial field excited by
whole-space magnetic dipole sources at an initial time after the current 
is switched off. Then the scheme steps Maxwell’s equations in time 
using a staggered grid and a modified DuFort-Frankel method. See the 
paper for detail: ***Fei Li and Jiulong Cheng, (2023), 3D finite-difference
transient electromagnetic modeling with a whole-space initial field, 
GEOPHYSICS 88: F15-F27. https://doi.org/10.1190/geo2021-0828.1.***

### WFTEM3D2.0:
The scheme steps Maxwell’s equations in time using a staggered grid 
and a modified DuFort-Frankel method: First, wires are modeled as volume 
currents, and the scheme calculates the primary field based on a whole-space 
homogeneous model; then it calculates the secondary field using the true model. 
The relevant paper will be provided in subsequent updates.

## Installation

**The Matlab version:** first the **Matlab** should be installed, 
then the code can run in Matlab without installation.

**The Fortran version:** first the **Visual Studio** and the **Intel Fortran Compiler** 
should be installed, then the code can run in Visual Studio without installation. 
You can also compile and run in VS Code or other development environments (no Makefile 
provided as the project is simple, you can easily create one based on the .f90 files).

**The exe version:** requires no installation—just run it directly.

## Usage

### 1. Edit the input file according to your model.

For the input file format, refer to the "README.md" in each version's folder.
Note: WFTEM3D1.0 and WFTEM3D2.0 use different input file formats.

### 2. Run the program.

**The Matlab version:** first open the **Main.m** file with Matlab, then click on 
the "Run" button in the toolbar or use the shortcut key F5 to run the program.

**The Fortran version:** first open the **WFTEM3D.vfproj** file with Visual Studio, 
then click on the "Run" button in the toolbar or use the shortcut key F5 to run the program.

**The exe version:** double-click **WFTEM3D2.0.exe** to run the program.

### 3. View calculation results.

Results (dBz/dt at every time instants) will be output automatically when the calculation is finished. 
There are three output files:

* Result_time.txt

The time instants are saved in this file. 

* Result_dBz.txt

The dBz/dt at receivers are saved in this file. 

* Run_time.txt

The run time is saved in this file.

## Contributing
Pull requests are welcome. We expect contributions via email with the 
corresponding author (email: figo1@163.com).
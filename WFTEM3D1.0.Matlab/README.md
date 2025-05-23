# WFTEM3D1.0

WFTEM3D1.0 is a quick, flexible and extensible 3D finite-difference 
transient electromagnetic modeling open-source software. It is 
applicable to the tunnel model, the non-flat topography model, 
the multi-scale model, etc. WFTEM3D has two language versions: 
a Matlab version (in the **WFTEM3D.Matlab** folder) and a Fortran version (in 
the **WFTEM3D.Fortran** folder). A simplified Python version is also provided,
but with slow calculation speed.

## Method

First the scheme calculates the approximate initial field excited by
whole-space magnetic dipole sources at an initial time after the current 
is switched off. Then the scheme steps Maxwell’s equations in time 
using a staggered grid and a modified DuFort-Frankel method. See the 
paper ***LiFei and JiulongCheng, 2023, 3D finite-difference transient 
electromagnetic modeling with a whole-space initial field, Geophysics, online: 
https://doi.org/10.1190/geo2021-0828.1*** for detail.

## Installation

The Matlab version: first the Matlab software should be installed, 
then the code can run in Matlab without installation.

The Fortran version: first the Visual Studio software and the Intel 
Visual Fortran software should be installed, then the code can run 
in Visual Studio without installation.

## Usage

We use the example of ***a Conductive brick in a half-space*** 
to show how to use WFTEM3D. Detailed description of this model can be
seen in the paper *3D Finite-difference Transient Electromagnetic 
Modeling with Whole-space Initial Field*.

### 1. Edit the input file according to your model.

***Example1_Conductive brick in a half-space.dat*** is the input file 
of the model ***a Conductive brick in a half-space***. It is in 
the **Data** folder and it is a standard input file：
``` 
100
17,18,18
1280,640,320,160,80,40,20,20,20,20,20,40,80,160,320,640,1280
1280,640,320,160,80,40,20,10,10,15,20,20,40,80,160,320,640,1280
1920,960,480,240,120,60,30,15,15,15,15,15,15,30,90,270,810,2430
9,8,9
8,8
2400
64
0.7
3
1,17,1,18,1,18,0.1
1,17,1,18,1,9,0.0003
7,11,11,12,12,13,2
```
Description of the input file:

``` 
Line 1: Length of Tx loop (m), current is 1 A by default:
        100
Line 2: Number of cells in the x-, y-, and z-directions:
        17,18,18
Line 3: Grid size in the x direction (m):
        1280,640,320,160,80,40,20,20,20,20,20,40,80,160,320,640,1280
Line 4: Grid size in the y direction (m):
        1280,640,320,160,80,40,20,10,10,15,20,20,40,80,160,320,640,1280
Line 5: Grid size in the z direction (m):
        1920,960,480,240,120,60,30,15,15,15,15,15,15,30,90,270,810,2430
Line 6: Tx loop center (i0,j0,k0, at the bottom face center of this cell):
        9,8,9
Line 7: Receivers (from i0,Rx_st,k0 to i0,Rx_end,k0, at the bottom face center of these cells):
        8,8
Line 8: Maximum number of iterations:
        2400
Line 9: Number of subloops:
        64
Line 10: Coefficient used to calculate time step size (typically ranges from 0.4 to 0.8)：
        0.7
Line 11: Number of model subareas:
        3
Line 12: Subarea 1. Cells range (x=1–17,y=1–18,z=1–18) and conductivity (0.1 S/m) of the background: 
        1,17,1,18,1,18,0.1
Line 13: Subarea 2, Cells range (x=1–17,y=1–18,z=1–9) and conductivity (0.0003 S/m) of the air:
        1,17,1,18,1,9,0.0003
Line 14: Subarea 3, Cells range (x=7–11,y=11–12,z=12–13) and conductivity (2 S/m) of the conductive brick:
        7,11,11,12,12,13,2
```
### 2. Run the program.

The Matlab version: first open the **Main.m** file with Matlab, then click on 
the "Run" button in the toolbar or use the shortcut key F5 to run the program.

The Fortran version: first open the **WFTEM3D.vfproj** file with Visual Studio, 
then click on the "Run" button in the toolbar or use the shortcut key F5 to run the program.

### 3. View calculation results.

Results (dBz/dt at every time instants) will be output automatically in 
the **Data** folder when the calculation is finished. 
There are three output files:

* Result_time.txt

The time instants are saved in this file. The example is shown as follows:
```
1.087338257278733E-007
1.361303417188709E-007
...
1.000056346847610E-002
1.000888668992117E-002
```
The unit of time is s.

* Result_dBz.txt

The dBz/dt at receivers are saved in this file. The example is shown as follows:
```
-0.2916895094E-004
-0.2574270181E-002
... 
-0.5137402423E-009
-0.5129357698E-009
```
The unit of dBz/dt is V/Am^2.

* Run_time.txt

The run time is saved in this file. The example is shown as follows:
```
 Computation finished. Run-time is 0.40241540000000 s.
```

## Contributing
Pull requests are welcome. We expect contributions via email with the 
corresponding author (email: figo1@163.com).
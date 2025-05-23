# WFTEM3D1.0

A quick, flexible, and extensible 3D TEM modeling open-source software. 
Supports loop source, and ground/airborne/marine/tunnel/borehole scenarios.

## Method

First the scheme calculates the approximate initial field excited by
whole-space magnetic dipole sources at an initial time after the current 
is switched off. Then the scheme steps Maxwell’s equations in time 
using a staggered grid and a modified DuFort-Frankel method. See the 
paper for detail: ***Fei Li and Jiulong Cheng, (2023), 3D finite-difference
transient electromagnetic modeling with a whole-space initial field, 
GEOPHYSICS 88: F15-F27.https://doi.org/10.1190/geo2021-0828.1.***

## Installation

This software requires no installation—just run it directly.

## Usage

We use the example of ***a conductive brick in a half-space*** 
to show how to use WFTEM3D1.0. 

### 1. Edit the input file according to your model.

***Example_conductive_brick_in_a_half-space.txt*** is the input file 
of the model ***a conductive brick in a half-space***：
``` 
##### Comment lines start with #, data is separated by spaces, ',', or ';'
##### Model description: A 0.1 S/m half-space contains a 2 S/m conductive brick (100m*40m*30m) 
##### at 30 m depth. Central-loop configuration, loop size: 100m*100m. 

##### Length of Tx loop (m), current is 1 A by default:
100
##### Number of cells in the x-, y-, and z-directions:
17,18,18
##### Grid size in the x direction (m):
1280,640,320,160,80,40,20,20,20,20,20,40,80,160,320,640,1280
##### Grid size in the y direction (m):
1280,640,320,160,80,40,20,10,10,15,20,20,40,80,160,320,640,1280
##### Grid size in the z direction (m):
1920,960,480,240,120,60,30,15,15,15,15,15,15,30,90,270,810,2430
##### Tx loop center (i0,j0,k0, at the bottom face center of this cell):
9,8,9
##### Receivers (from i0,Rx_st,k0 to i0,Rx_end,k0, at the bottom face center of these cells):
8,8
##### Maximum number of iterations:
2400
##### Number of subloops:
64
##### Coefficient used to calculate time step size (typically ranges from 0.4 to 0.8)：
0.7
##### Number of model subdomains:
3
##### Subdomain 1. Cells range (x=1–17,y=1–18,z=1–18) and conductivity (0.1 S/m) of the background: 
1,17,1,18,1,18,0.1
##### Subdomain 2, Cells range (x=1–17,y=1–18,z=1–9) and conductivity (0.0003 S/m) of the air:
1,17,1,18,1,9,0.0003
##### Subdomain 3, Cells range (x=7–11,y=11–12,z=12–13) and conductivity (2 S/m) of the conductive brick:
7,11,11,12,12,13,2
```
### 2. Run the program.

First, double-click ***WFTEM3D1.0.exe***, then enter the input filename, and press Enter to start the calculation.

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
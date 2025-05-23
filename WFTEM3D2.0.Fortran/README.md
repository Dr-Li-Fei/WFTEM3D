# WFTEM3D2.0

A quick, flexible, and extensible 3D TEM modeling open-source software. 
Supports wire/loop sources, half/whole-space, and ground/airborne/marine/tunnel/borehole scenarios.

## Method

The scheme steps Maxwell’s equations in time using a staggered grid 
and a modified DuFort-Frankel method: First, wires are modeled as volume 
currents, and the scheme calculates the primary field based on a whole-space 
homogeneous model; then it calculates the secondary field using the true model. 
The relevant paper will be provided in subsequent updates.

## Installation

first the **Visual Studio** and the **Intel Fortran Compiler** 
should be installed, then the code can run in Visual Studio without installation. 
You can also compile and run in VS Code or other development environments (no Makefile 
provided as the project is simple, you can easily create one based on the .f90 files).

## Usage

We use the example of ***a conductive brick in a half-space*** 
to show how to use WFTEM3D2.0. 

### 1. Edit the input file according to your model.

***Example_conductive_brick_in_a_half-space.txt*** is the input file 
of the model ***a conductive brick in a half-space***:
``` 
##### Comment lines start with #, data is separated by spaces, ',', or ';'.
##### Model description: A 0.1 S/m half-space contains a 2 S/m conductive brick (100 m × 40 m × 30 m) 
##### at 30 m depth. Central-loop configuration: 100 m × 100 m. 

##### Number of cells in the x-, y-, and z-directions.
27 25 25
##### Grid size in the x direction (m).
2560 1280 640 320 160 80 40 20 10 10 15 10 10 10 15 10 10 10 10 20 40 80 160 320 640 1280 2560
##### Grid size in the y direction (m).
2560 1280 640 320 160 80 40 20 10 10 15 10 10 10 15 10 10 20 40 80 160 320 640 1280 2560
##### Grid size in the z direction (m).
5120 2560 1280 640 320 160 80 40 20 10 5 5 10 15 15 15 20 40 80 160 320 640 1280 2560 5120

##### Tx position:
9 9 11 # (x1, y1, z1).
17 17 11 # (x2, y2, z2).
##### Loop source: Aligns with the outer edge of the cells between (x1, y1, z1) and (x2, y2, z2).
##### Wire source: Aligns with the -x or -y edge of the cells between (x1, y1, z1) and (x2, y2, z2).
##### Transmitter depth: at the bottom plane of these cells. Current=1 A.
# To display following graphic correctly, please view in monospaced font (e.g., Courier New/Lucida Console).
#    +=======+=======+=======+          +-------+-------+-------+          +-------+-------+-------+
#   || x1,y1 |       |       ||         |      || x1,y1 |       |          |       |       |       |
#    +-------+-------+-------+          +-------+-------+-------+          +=======+=======+=======+
#   ||       |       |       ||         |      ||       |       |          | x1,y1 |       | x2,y2 |  
#    +-------+-------+-------+          +-------+-------+-------+          +-------+-------+-------+ 
#   ||       |       | x2,y2 ||         |      || x2,y2 |       |          |       |       |       |
#    +=======+=======+=======+          +-------+-------+-------+          +-------+-------+-------+
#               Loop                         Wire (if x1=x2)                    Wire (if y1=y2)

##### Rx position:
13 13 1 # Rx along x-direction: start/end/interval.
13 13 1 # Rx along y-direction: start/end/interval.
11 11 1 # Rx along z-direction: start/end/interval.
##### Receivers are located at the centers of the bottom surfaces of these cells.

##### Iteration numbers of primary field and secondary field.
600 6400
##### Coefficients for calculating time step size of primary field and secondary field.
0.8 0.8

##### Number of model subdomains.
3
##### Cells range (x=?–?, y=?–?, z=?–?) and conductivity (?) of each subdomain.
1	27	1	25	1	25	0.1 #Cells range (x=1–27, y=1–25, z=1–25) and conductivity (0.1 S/m) of the background.
1	27	1	25	1	11	0.0003 #Cells range (x=1–27, y=1–25, z=1–11) and conductivity (0.0003 S/m) of the air.
16	19	9	17	15	16	2 #Cells range (x=16–19, y=9–17, z=15–16) and conductivity (2 S/m) of the conductive brick.
```
### 2. Run the program.

first open the **WFTEM3D2.0.vfproj** file with Visual Studio, 
then click on the "Run" button in the toolbar or use the shortcut key F5 to run the program.

### 3. View calculation results.

Results (dBz/dt at every time instants) will be output automatically when the calculation is finished. 
There are three output files:

* Result_time.txt

The time instants are saved in this file. The example is shown as follows:
```
1.2286042e-08
1.5792647e-08
...
1.0299405e-02
1.0302622e-02
```
The unit of time is s.

* Result_dBz.txt

The dBz/dt at receivers are saved in this file. The example is shown as follows:
```
-1.0170416e-04
-1.0553606e-04
... 
-4.7702894e-10
-4.7667189e-10
```
The unit of dBz/dt is V/Am^2.

For multiple receivers, each as a column, the columns are arranged from left to right in the following order:
the x-direction varies first, followed by the y-direction, and then the z-direction, e.g., 
(1,1,1), (2,1,1), (1,2,1), (2,2,1), (1,1,2), (2,1,2), (1,2,2), (2,2,2).

* Run_time.txt

The run time is saved in this file. The example is shown as follows:
```
 Computation finished. Run-time is 1.09375000000000 s.
```

## Contributing
Pull requests are welcome. We expect contributions via email with the 
corresponding author (email: figo1@163.com).
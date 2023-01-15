# MMA-CompoundMatrixMethod
Solving boundary-value eigenvalue problems in Mathematica using the Compound Matrix Method to construct the Evans function.

Mathematical details may be found <a href=http://www.maths.gla.ac.uk/~xl/FYB-background.pdf>here </a> and at this <a href=https://mathematica.stackexchange.com/questions/155079/finding-eigenvalues-for-a-boundary-value-problem>stack exchange question</a>. I'm working on a publication with the details, please feel free to email me (address below) to see my current version of it.


## How to Download the latest release
 The initial version is available as a  `.paclet` file. Download and install it using the `PacletInstall` function in Mathematica:
 
     Needs["PacletManager`"]
     PacletInstall["CompoundMatrixMethod", "Site" -> "http://raw.githubusercontent.com/paclets/Repository/master"]
        
Alternatively, download the paclet locally and install using `PacletInstall` on the appropriate directory. For example, if the file was downloaded to the directory `~/Downloads`, evaluate  `PacletInstall["~/Downloads/CompoundMatrixMethod-0.9.paclet"]`. If all else fails, you can copy the files from github to a folder called CompoundMatrixMethod in  the AddOns/Applications folder on your computer.

The package can then be loaded by calling 

        Needs["CompoundMatrixMethod`"]

## Usage

The Compound Matrix Method is a package for finding eigenvalues of boundary-value ordinary differential equations.

This includes problems with a single interface, decaying conditions at one or both ends, in an upto 10th order differential system.

First we need to transform the boundary-value problem (BVP) into a set of first order matrix equations. The function ToMatrixSystem will do this, linearising the equations if necessary (with a warning if it does). 

        sys=ToMatrixSystem[y''[x] + k^2 y[x] == 0, {y[0] == 0, y[1] == 0}, y, {x, 0, 1}, k]

This will store the system into the variable `sys`. The syntax is similar to that of ParametricNDSolve, with the differential equations, boundary conditions, dependent variables, independent variable and eigenvalue.

We can then evaluate the Evans function at a given guess of the eigenvalue `k` (here k=1):

        Evans[1,sys]
        -0.841471
    
Zeros of this function correspond to eigenvalues of the original BVP: 

    FindRoot[Evans[k, sys], {k, 1}]
    {k -> 3.14159}
    
The function is smooth and can be plotted by the built-in routines:
    
    Plot[Evans[k, sys], {k, 0, 15}]
   
A number of further examples are shown in the file `CMMExamples.nb`, available from this respository. This includes examples with boundary conditions at infinity, higher order equations (up to 10th order), split domains with interface conditions and when higher precision is required. Also check out where I have <a href="https://mathematica.stackexchange.com/search?q=compoundmatrixmethod" > answered questions on stackexchange using my package</a>.


## Citations

I used this method to solve an eigenvalue problem in <a href=https://doi.org//10.1093/imamat/hxq026>my 2010 paper </a>, and the package itself in both <a href=https://link.springer.com/article/10.1007/s11538-018-0505-4>a tenth-order ODE </a> as well as <a href=https://journals.aps.org/pre/abstract/10.1103/PhysRevE.98.033003>an example with an interface</a>. 

I have a half-written expository paper (expansion of the <a href=http://www.maths.gla.ac.uk/~xl/FYB-background.pdf>note</a> written by Professor Yibin Fu) to give more details on how the method works and introduce the package, please email me for a copy of that. 

## Contact

Feel free to contact me if you have any questions, suggestions or issues, I'm also interested in collaborations involving this work, but please note my time is severely limited as I am no longer active in this field. My email address is simon (dot) pearce (at) cruk (dot) manchester (dot) ac (dot) uk. 

## Funding Acknowledgement

This code was initially written while I held an Early Career Fellowship from the <a href=https://www.leverhulme.ac.uk>Leverhulme Trust</a>. I'm now based at the Cancer Research UK Manchester Institute, part of the University of Manchester, UK, and work on bioinformatics.

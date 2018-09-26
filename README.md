# MMA-CompoundMatrixMethod
Solving boundary-value eigenvalue problems in Mathematica using the Compound Matrix Method to construct the Evans function.

Mathematical details may be found <a href=http://www.maths.gla.ac.uk/~xl/FYB-background.pdf>here </a> and at this <a href=https://mathematica.stackexchange.com/questions/155079/finding-eigenvalues-for-a-boundary-value-problem>stack exchange question</a>. 


## How to Download the latest release
 The initial version is available as a  `.paclet` file. Download and install it using the `PacletInstall` function in Mathematica:
 
     Needs["PacletManager`"]
     PacletInstall["CompoundMatrixMethod", "Site" -> "http://raw.githubusercontent.com/paclets/PacletServer/master"]
        
 Alternatively, download the paclet locally and install using `PacletInstall` on the appropriate directory. For example, if the file was downloaded to the directory `~/Downloads`, evaluate  `PacletInstall["~/Downloads/CompoundMatrixMethod-0.9.paclet"]`

The package can then be loaded by calling 

        Needs["CompoundMatrixMethod`"]

## Usage

The Compound Matrix Method is a package for finding eigenvalues of boundary-value ordinary differential equations.

First we need to transform the boundary-value problem (BVP) into a set of first order matrix equations. The function ToLinearMatrixForm will do this, linearising the equations if necessary (with a warning if it does). 

        sys=ToLinearMatrixForm[y''[x] + k^2 y[x] == 0, {y[0] == 0, y[1] == 0}, y, {x, 0, 1}, k]

This will store the system into the variable `sys`. We can then evaluate the Evans function at a given guess of the eigenvalue `k` (here k=1):

        Evans[1,sys]
        -0.841471
    
Zeros of this function correspond to eigenvalues of the original BVP: 

    FindRoot[Evans[k, sys], {k, 1}]
    {k -> 3.14159}
    
The function is smooth and can be plotted by the built-in routines:
    
    Plot[Evans[k, sys], {k, 0, 15}]
   
A further set of examples are shown in the file `CMMExamples.nb`.

## Citations

I used this method to solve an eigenvalue problem in <a href=https://doi.org//10.1093/imamat/hxq026>my 2010 paper </a> (also available on the arXiv), as well as <a href=https://journals.aps.org/pre/abstract/10.1103/PhysRevE.98.033003>an example with an interface</a>. I'm currently working on an expository paper to detail how the method works and introduce the package.

## Contact

Feel free to contact me if you have any questions, suggestions or issues, I'm also interested in collaborations involving this work.
My email address is simon (dot) pearce (at) manchester (dot) ac (dot) uk. 

## Funding Acknowledgement

This code was initially written while I held an Early Career Fellowship from the <a href=https://www.leverhulme.ac.uk>Leverhulme Trust</a>.

# MMA-CompoundMatrixMethod
Solving boundary-value eigenvalue problems in Mathematica using the Compound Matrix Method

 - [Download the latest release](https://github.com/krazug/CompoundMatrixMethod/releases), 
 The initial version  is available as a  `.paclet` file. Download and install it using the `PacletInstall` function in Mathematica.  For example, assuming that the file `CompoundMatrixMethod-0.1.paclet` was downloaded into the directory `~/Downloads`, evaluate

        Needs["PacletManager`"]
        PacletInstall["~/Downloads/CompoundMatrixMethod-0.1.paclet"]

The package can then be loaded by calling Needs["CompoundMatrixMethod`"].

## Usage

First we need to transform the BVP into a set of first order matrix equations. The function ToLinearMatrixForm will do this, linearising the equations if necessary. 

        sys=ToLinearMatrixForm[y''[x] + k^2 y[x] == 0, {y[0] == 0, y[1] == 0}, y, {x, 0, L}]

This will return the matrices A, B and C that are required. We can then use the Compound Matrix Method to evaluate the Evans function at a given guess of the eigenvalue `k` (here k=1):

        CompoundMatrixMethod[{k,1},sys]
        -0.841471
    
Zeros of this function correspond to eigenvalues of the original BVP: 

    FindRoot[CompoundMatrixMethod[{k, kk}, sys], {kk, 1}]
    {kk -> 3.14159}
    
The function is smooth and can be plotted by the built-in routines:
    
    Plot[CompoundMatrixMethod[{k, kk}, sys], {kk, 0, 15}]
   

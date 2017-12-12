# MMA-CompoundMatrixMethod
Solving boundary-value eigenvalue problems in Mathematica using the Compound Matrix Method
Mathematical details on the Compound Matrix Method may be found <a href=http://www.maths.gla.ac.uk/~xl/FYB-background.pdf>here </a> and at this <a href=https://mathematica.stackexchange.com/questions/155079/finding-eigenvalues-for-a-boundary-value-problem>stack exchange question</a>. 


## How to Download the latest release
 The initial version is available as a  `.paclet` file. Download and install it using the `PacletInstall` function in Mathematica:  
 
        Needs["PacletManager`"]
        PacletInstall["https://github.com/krazug/CompoundMatrixMethod/releases/download/v0.3/CompoundMatrixMethod-0.3.paclet"]
        
 Alternatively, download the paclet locally and install using `PacletInstall` on the appropriate directory. For example, if the file was downloaded to the directory `~/Downloads`, evaluate  `PacletInstall["~/Downloads/CompoundMatrixMethod-0.3.paclet"]`

The package can then be loaded by calling 

        Needs["CompoundMatrixMethod`"]

## Usage

The Compound Matrix Method is a package for finding eigenvalues to boundary-value differential equations with an unknown parameter.

First we need to transform the BVP into a set of first order matrix equations. The function ToLinearMatrixForm will do this, linearising the equations if necessary. 

        sys=ToLinearMatrixForm[y''[x] + k^2 y[x] == 0, {y[0] == 0, y[1] == 0}, y, {x, 0, 1}]

This will return the matrices A, B and C that are required. We can then use the Compound Matrix Method to evaluate the Evans function at a given guess of the eigenvalue `k` (here k=1):

        CompoundMatrixMethod[{k,1},sys]
        -0.841471
    
Zeros of this function correspond to eigenvalues of the original BVP: 

    FindRoot[CompoundMatrixMethod[{k, k0}, sys], {k0, 1}]
    {k0 -> 3.14159}
    
The function is smooth and can be plotted by the built-in routines:
    
    Plot[CompoundMatrixMethod[{k, k0}, sys], {k0, 0, 15}]
   
A further six examples are shown in the file `CMMExamples.nb`.

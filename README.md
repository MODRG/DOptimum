# About
This repository contains a Python wrapper for Algorithm 295 published in the Journal of the Royal Statistical Society in 1994.  Algorithm 295 provides a FORTRAN implementation of a Fedorov Exchange algorithm for D-Optimal design.  The original source code, a description of the algorithm and a test case and example is provided.

## Algorithm obtained from
-  **Algorithm AS 295:** A Fedorov Exchange Algorithm for D-Optimal Design
-  **Author(s):** Alan J. Miller and Nam-Ky Nguyen
-  **Source:** Journal of the Royal Statistical Society. Series C (Applied Statistics), Vol. 43, No. 4, pp. 669-677, 1994
-  **Stable URL:** http://www.jstor.org/stable/2986264

## Source code from
-  http://ftp.uni-bayreuth.de/math/statlib/apstat/

# Information
This project was setup to work for both Windows and Linux environments although it was only tested in a Linux environment.  Minor tweaks may thus be required to make this work on Windows.  The [f2py](https://docs.scipy.org/doc/numpy/f2py/) tool from Python is used to compile the supplied FORTRAN code into a shared library.  Typically this is a two step process as follows:

1. Create a 295.pyf signature file that corresponds to the 295.f FORTRAN code
```
f2py3 -h 295.pyf -m dopt 295.f
```
2. Compile the shared library using the created signature file
```
f2py3 -c 295.pyf 295.f
```

For the case presented here, the 295.f signature file was heavily modified to specify arguments that are already known.  This simplifies the Python interface as the use have to specify fewer arguments.  When compiling the shared library, **please do not** recreate the signature file (step 1 above) as this will overwrite the these changes, only do the compilation (step 2 above).

Note that the block size has been hard coded in the signature file to be 0 and that the user must provide the full X matrix (including constant term if required) to the algorithm.

# Requirements
The test case and example problem is provided as [Jupyter](https://jupyter.org/) notebooks.  Please make sure that you have a [Jupyter](https://jupyter.org/) environment to run these notebooks.
- [NumPy](http://www.numpy.org/) Numerical Python library
- [Pandas](https://pandas.pydata.org/) Python data analysis library
- [StatsModels](https://www.statsmodels.org/stable/index.html) Statistical package used for linear regression

# File description
| File        | Description  |
| ----------- |------------- |
| [295.f]()  | The Algorithm 295 FORTRAN source code |
| [295.pyf]() | The corresponding signature file, created with f2py and manually modified |
| [295.pdf]() | A description of the algorithm |
| [test_295.ipynb]() | A basic test case to test the implementation (Jupyter notebook) |
| [example.ipynb]() | A two variable, quadratic model example problem from Meyers and Montgommery (Jupyter notebook)|
| [MyersExample.xlsx]() | The data required for example.ipynb from Meyers and Montgommery |


# Getting started
1. Clone this library (command line command for Linux are provided below)
```
git clone https://github.com/MODRG/DOptimum.git
```
2. Change to the DOptimum folder
```
cd DOptimum
```
3. Compile the shared library for your Python environment:
```
f2py3 -c 295.pyf 295.f
```
4. Run the example
```
jupyter notebook example.ipynb
```

# Authors
- Gerhard Venter - First publication 3 September 2019

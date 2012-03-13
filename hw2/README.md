## Homework 2 -- Comparing Coordinate Descent and Interior Point Methods for sparse regression

### First things first
* Before running any code, open up `MATLAB`, cd to the directory that contains this file and type **`mex tdma.c`** on the command line to build the necessary MEX file.
* Most of the functions in the tarball have halfway decent documentation, so look there for more details, .e.g, type `doc descent` at the `MATLAB` command line to bring up the documentation for `descent.m`.

### Files
* `README.html` -- View in your favorite web browser.
* `src/bumpgen.m` -- function that generates some fake data in the form of weighted Gaussian basis functions.
* `src/center.m` -- goofy little function to center and scale a matrix and vector
* `src/tdma.c` -- code to solve a tridiagonal system of equations
* `src/coordesc.m` -- implementation of the coordinate descent algorithm for sparse regression
* `src/test_intpoint.m` -- code to run an example of the interior point method for sparse regression
* `src/intpoint.m` -- interior point solver
* `src/test_coordesc.m` -- code to run an example of the coordinate descent method for sparse regression
* `src/descent.m` -- the main loop of the coordinate descent algorithm; see the file for more detailed documentation
* `src/pdfsave.m` -- extremely simplistic function to generate a PDF of the current figure
* `pdf/svm.pdf` -- The first part of the homework assignment, which was to write down the dual of a support vector machine quadratic program.
* `pdf/hw2_interior_point.pdf` -- Plots of interior point runs.
* `pdf/hw2_coordinate_descent.pdf` -- Plots of coordinate descent runs.
* `pdf/boyd-l1-regression.pdf` -- Paper on the interior point method that I used

### Results
* Coordinate descent seems to be **much** faster than the interior point method, but doesn't give as good of a fit.
* Coordinate descent would probably benefit from some C code speed up.
* The interior point method can fit the data perfectly, but is much slower than the coordinate descent method when there are a large number of functions. Indeed, the running time complexity scales with O(*d*) (as expected) where *d* is the number of functions.

README.md rev. 02 January 2014 by Stuart Ambler.
Copyright (c) 2014 Stuart Ambler.
Distributed under the Boost License in the accompanying file LICENSE.

# Automatic Multiscale-based Peak Detection (AMPD)

ampd.cc implements automatic multiscale-based peak detection (AMPD) algorithm
as in An Efficient Algorithm for Automatic Peak Detection in Noisy Periodic and
Quasi-Periodic Signals, by Felix Scholkmann, Jens Boss and Martin Wolf,
Algorithms 2012, 5, 588-603.

ampd-arm.cc is the same except it uses Mat, Col, etc. from the Armadillo 18.5
C++ library rather than using std::valarray etc.

The code uses some C++ 11 syntax and was compiled with gcc 4.8.1.  It includes a
main function with a number of parameters for generation of test data and a few
parameters for operation of algorithm.  The test outputs command and data files
for gnuplot; tested with version 4.6 patchlevel 3.  Environment lubuntu 13.10,
intel processor, 8 GB RAM.  Requires installation of argtable2; tested with
version 12-1.  Build script usage:

buildcc11 ampd

Since the program writes a number of output files with long filenames including
all the command line parameters or defaults used, it may be best to create a
directory for each test run with a certain set of parameters.  In that
directory, after executing the program, for example, to use defaults, with

path-to-executable/ampd

replacing path-to-executable with the actual path, gnuplot can be run with

gnuplot *cmd.txt

and the resulting plots can be viewed with gimp (test with version 2.8.6) with

gimp *png &

There are comments indicating that certain code, useful for the plots, could be
removed for simple use of the algorithm.  There is a long comment at the start
describing the implementation, including a minor correction and a minor change
made to the published algorithm.
# DAS
This folder holds the 2011 version of DAS3 used by the Eric Schearer lab group

## Building the model

Preparation:
* A "make" tool should be installed
* Install the CLANG C compiler. This is used to compile das3_al.c and is much faster than other compilers (compile time and run time).
* Use "mex -setup" in Matlab and select Microsoft Visual C++ 2019 for C compilation.  This is used to build the final MEX binary.
Currently, MinGW gives an error message "undefined reference to `__chkstk'".

Edit the script make.m so it can find "make.exe" or equivalent on your computer.  Then build as follows:

	>> make

We assume that the kinematic structure of the model is never changed.  If you do that, das3_al_raw.c should be rebuilt which requires Autolev.  The Makefile will automatically do this, but only on Ton's computer which has an Autolev license.

When the model is initialized, the model parameters are loaded from das3.bio.  If you ever change the model parameters (mass properties,
muscle paths, muscle properties) in the das3.osim file, you would go into the tools folder and rebuild the das3.bio file.  These tools
have not been tested recently.  Some parameters can be easily edited manually in das3.bio.

## Using the model

Information about can be found in das3implicit.docx.

Before using the model, it is a good idea to run some of the tests in das3test.m.  Some of the tests are currently not
working correctly (please post issues on github as needed!).

When using the fixed-step simulation method (das3step.m or the corresponding C function),
the simulation can go unstable when the time step is too large.  A time step of 3 ms seems to work
well, as long as the control input u does not change too quickly.  A change in u of 0.5 in one time
step seems to be OK with a 3 ms time step.

The C version of das3step (which is coded in das3mex.c) executes in about 1 ms, so time steps of 
2 or 3 ms can be comfortably done in real time.




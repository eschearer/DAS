# DAS
This folder holds the 2011 version of DAS3 used by the Eric Schearer lab group

## Building the model

Tools needed:
* A "make" tool should be installed and in the OS search path
* A C compiler that is suitable for building MEX functions (follow Matlab help instructions)
* The CLANG C compiler for compiling the lengthy das3_al.c file.  It may be possible to use the other C compiler for this, but it could be very slow or run out of memory.

Edit the script make.m so it can find "make" on your computer.  Then build as follows:

	>> make

We assume that the kinematic structure of the model is never changed.  If you do that, das3_al_raw.c should be rebuilt which requires Autolev.  The Makefile will automatically do this, but only on Ton's computer which has an Autolev license.

When the model is initialized, the model parameters are loaded from das3.bio.  If you ever change the model parameters (mass properties,
muscle paths, muscle properties) in the das3.osim file, you would go into the tools folder and rebuild the das3.bio file.  These tools
have not been tested recently.  Some parameters can be easily edited manually in das3.bio.

## Using the model

Information about can be found in das3implicit.docx.

Before using the model, it is a good idea to run some of the tests in das3test.m.  Some of the tests are currently not
working correctly (please post issues on github as needed!).




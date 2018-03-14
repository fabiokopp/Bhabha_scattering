ErrorBarLogPlots

Original code Copyright 1988-2007 Wolfram Research, Inc.
Modifications to code Copyright 2007, California Institute of Technology
Modifications written by Frank Rice

The ErrorBarLogPlots package is written for Mathematica 6.0. It includes the
functions:
	ErrorListLogPlot
	ErrorListLogLinearPlot
	ErrorListLogLogPlot

which add log-scale plotting to suplement Mathematica's ErrorListPlot function.

Also included is the function ErrorBarScale which can be used as a setting of
the ErrorBarFunction option to ErrorListPlot.

Unzip the files into a convenient directory. To install the package automatically
to the proper Mathematica Applications directory, launch the notebook:
	ErrorBarLogPlots_Installer.nb

It will execute when loaded by the Mathematica front-end (you may be prompted to
run the initialization cells - you should allow it to do so). The notebook will
copy the file ErrorBarLogPlots.m into the appropriate directory for Mathematica
to be able to find it; you may need Administrator or Root privileges to install
the package for all users.

Once installed, you may load the package during a Mathematica session using
	<<ErrorBarLogPlots`
or
	Needs["ErrorBarLogPlots`"]

The installer notebook contains some very brief examples of uses of the functions.


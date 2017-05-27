# vsop2013

## A Python implementation of the VSOP2013 Ephemeris functionality.

This module provides basic support for the VSOP2013 ephemerides which can be used to compute the
heliocentric positions and velocities of the 8 planets Mercury, Venus, the Earth-Moon barycenter,
Mars, Jupiter, Saturn, Uranus, Neptune and the dwarf planet Pluto (only in the period from +1500
to +3000).

The ephemerides are provided in 6 sequential text files (ASCII) and can be downloaded at:

    ftp://ftp.imcce.fr/pub/ephem/planets/vsop2013/ephemerides
    
The following provides the name of each file and the period that it covers:

* `VSOP2013.m4000: -4500 to -3000`
* `VSOP2013.m2000: -3000 to -1500`
* `VSOP2013.m1000: -1500 to 0`
* `VSOP2013.p1000: 0 to +1500`
* `VSOP2013.p2000: +1500 to +3000`
* `VSOP2013.p4000: +3000 to +4500`

While these text files could be used directly, this would have a serious impact on the performance
of computing any significant number of planetary positions. Instead it is suggested that these
files are parsed and stored into binary format for the platform on which they will be used.
**NOTE THAT THE BINARY FILE GENERATED FOR ONE PLATFORM MAY NOT BE COMPATIBLE FOR USE ON A DIFFERENT
PLATFORM.** Once converted to binary form the text files may be deleted to save storage space where
necessary. The binary form of the files provides a much more efficient mechanism for performing
large numbers of position calculations.

For reference the authors have provided the programs (written in Fortran) that they used for binary
file generation and basic calculations.

These are:

VSOP2013_binfile.f : This program converts the sequential files into direct access files.

VSOP2013_compute.f : This program computes planetary coordinates from a direct access file.

The file VSOP2013_ctl.txt contains planetary coordinates computed by the program VSOP2013_compute.f
and given as control values for the users.

The file README.pdf also provides this information as well as other info that may be of interest.

This module has been developed to provide a degree of convenience for the Python user, not just
a reference implementation.

TXT_FILES_ROOT is the path to the directory containing the VSOP2013 source text files. You may
change this to reflect where your source text files are actually located. The default is to place
the source text files into a subdirectory called "ephemerides" at the location of the vsop2013.py
script. If the vsop2013.py functionality is being used by other python scripts then this variable
can be changed once the module is loaded.

TXT_FILES is an array of the filenames of the VSOP2013 source text files.

BIN_FILES_ROOT is the path to the directory containing the VSOP2013 binary files. You may
change this to reflect where your binary files are actually located. This variable is similar to
the TXT_FILES_ROOT variable discussed previously.

BIN_FILES is an array of the filenames of the VSOP2013 binary files.

PLANET NAMES is an array of the names of the supported planets and is used ONLY for the display
of the results by the print_results function.

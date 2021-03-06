"README"

(I have no affiliation with Agilent, I've just put this repository up to have
a place to collect improvements to the open-source *Xrecon* that's distributed
with VnmrJ.  I'm happy to incorporate any improvements/ports/updates!  Just
submit a pull request over github.)

Requirements
- GSL
- FFTW3
- libTIFF
- (cmake)

Build
=====
```sh
cmake .
make
```

or

```
make -f Makefile.orig
```

Usage
=====
- as a drop-in replacement for the existing
- for stand-alone command-line image reconstruction

```
Usage:
    Xrecon -h <...>
        print this help message

    Xrecon -v <paths to VnmrJ acquisition>...
        called from VnmrJ on unsaved acquisition data, outputs
        reconstruction to (TODO...?)

    Xrecon <paths to saved data>...
        command line invocation on saved acquisition data, outputs
        reconstruction to (.../data.img?)
```

--- snip: original README follows: ---
Xrecon - External Reconstruction
================================

This is Xrecon, a collection of magnetic resonance image reconstruction
routines tuned to read data acquired on Varian MRI Systems.
It has been developed to work externally and independently of the Varian data
acquisition software VnmrJ.

The routines provided here link against the GNU Scientific Library (GSL) and
the Fastest Fourier Transform in the West (FFTW) library of subroutines.
The GSL is free software published under the GNU General Public License (GPL).
The FFTW is free software published under the GNU General Public License (GPL).
As such it is a requirement that the routines provided here are also published
under the GPL.

To publish under the GPL two elements must be added to each source file:
1. A copyright notice with the names of all contributors.
2. A statement of copying permission saying the program is distributed under
   the terms of the GNU General Public License.

A copy of the license itself should also be included.

Xrecon is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

Xrecon is distributed in the hope that it will be useful, but WITHOUT ANY
WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE. See the GNU General Public License for more details.


Installation
============

Please consult the INSTALL file in this distribution for instructions.


Files
=====

- COPYING         The GPL
- INSTALL         Installation instructions
- README          This file
- Makefile        Xrecon makefile
- data.h          Stripped down Varian data file handler
- Xrecon.c        Main External Reconstruction
- Xrecon.h        Global includes and defines

recon1D
-------
- recon1D.c       1D recon
- default1D.c     Default 1D recon

recon2D
-------
- recon2D.c       2D recon
- default2D.c     Default 2D recon
- profile2D.c     Profile of 2D experiment

reconEPI
--------
- reconEPI.c      EPI recon
- defaultEPI.c    Default EPI recon
- dprocEPI.c      EPI data processing
- prescanEPI.c    EPI prescan recon

recon3D
-------
- recon3D.c       3D recon
- default3D.c     Default 3D recon

common1D
--------
- dproc1D.c       1D data processing
- dread1D.c       1D data read
- dutils1D.c      1D data utilities
- write1D.c       1D phasefile and datafaile writing

common2D
--------
- dmask2D.c       2D data masking
- dproc2D.c       2D data processing
- dread2D.c       2D data read
- dutils2D.c      2D data utilities
- noise2D.c       2D routines using noise measurements
- fdfwrite2D.c    2D fdf writing
- rawwrite2D.c    2D raw binary writing
- tifwrite2D.c    2D TIFF writing

common3D
--------
- dproc3D.c       3D data processing
- fdfwrite3D.c    3D fdf writing
- rawIO3D.c       3D raw binary writing

common
------
- dhead.c         Data header routines
- dutils.c        Data utilities
- options.c       Input options
- pars.c          Parameter read routines
- utils.c         General utilities

nifti
-----
- niftiwrite.c    NIFTI-1/Analyze7.5 writing
- nifti1.h        NIFTI-1 header


Contacts
========

~~Paul Kinchesh,
Agilent Technologies (formerly Varian Ltd), 10 Mead Road, Oxford Industrial Park, Oxford OX5 1QU, UK
paul.kinchesh@agilent.com~~


Bug Reports
===========

~~Paul Kinchesh, paul.kinchesh@agilent.com~~
Paul wrote XRecon, but AFAIK, he is no longer actively working on it.  Please add bug reports to this github project.

Contributors
============

Paul Kinchesh, Magnetic Resonance Systems, Agilent Technologies (formerly Varian Ltd).
Alexandr Khrapichev, Department of Physiology, Oxford University.
Martyn Klassen, Robarts Research Institute, University of Western Ontario.

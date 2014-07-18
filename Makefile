#
# Makefile #
#*---------------------------------------------------------------------------*#
#*                                                                           *#
#* Makefile: Makefile for Xrecon                                             *#
#*                                                                           *#
#* This file is part of Xrecon.                                              *#
#*                                                                           *#
#* Xrecon is free software: you can redistribute it and/or modify            *#
#* it under the terms of the GNU General Public License as published by      *#
#* the Free Software Foundation, either version 3 of the License, or         *#
#* (at your option) any later version.                                       *#
#*                                                                           *#
#* Xrecon is distributed in the hope that it will be useful,                 *#
#* but WITHOUT ANY WARRANTY; without even the implied warranty of            *#
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              *#
#* GNU General Public License for more details.                              *#
#*                                                                           *#
#* You should have received a copy of the GNU General Public License         *#
#* along with Xrecon. If not, see <http://www.gnu.org/licenses/>.            *#
#*                                                                           *#
#*---------------------------------------------------------------------------*#
#**#

# Installation Directory
BINDIR = $(HOME)/bin
# Windoze
WBINDIR = ~/bin

# Source
SRC = \
Xrecon.c \
recon1D/recon1D.c \
recon1D/default1D.c \
recon2D/recon2D.c \
recon2D/default2D.c \
recon2D/profile2D.c \
recon3D/recon3D.c \
recon3D/default3D.c \
reconEPI/reconEPI.c \
reconEPI/defaultEPI.c \
reconEPI/dprocEPI.c \
reconEPI/prescanEPI.c \
nifti/niftiwrite.c \
common1D/dread1D.c \
common1D/dproc1D.c \
common1D/dutils1D.c \
common1D/write1D.c \
common2D/dread2D.c \
common2D/dproc2D.c \
common2D/noise2D.c \
common2D/dmask2D.c \
common2D/dutils2D.c \
common2D/fdfwrite2D.c \
common2D/rawwrite2D.c \
common2D/tifwrite2D.c \
common3D/fdfwrite3D.c \
common3D/rawIO3D.c \
common3D/dproc3D.c \
common/dproc.c \
common/dhead.c \
common/dutils.c \
common/options.c \
common/pars.c \
common/utils.c

# Objects
OBJS  = $(SRC:.c=.o)

# Includes, Library Defines
INCL  =
LIBS  = -lgsl -lgslcblas -lfftw3 -ltiff -lm

# Executable
PROG  = Xrecon
# Windoze
WPROG = $(PROG).exe

# Compiler, Linker Defines
CC      = gcc
CFLAGS  = -pedantic -Wall -O4 -std=gnu9x
LIBPATH = -L /usr/local/lib
LDFLAGS = -o $(PROG) $(LIBPATH) $(LIBS)
DEBUG  = -pedantic -Wall -O4 -g -std=gnu9x -DDEBUG

# Compile and Assemble C Source Files into Object Files
%.o: %.c
	$(CC) -c $(CFLAGS) -o $*.o $*.c

# Link all Object Files with external Libraries into Binaries
$(PROG): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS)

# Objects depend on these Libraries
$(OBJS): $(INCL)

# DEBUG flags turned on
debug:
	$(CC) $(DEBUG) $(SRC) $(LDFLAGS)

install:
	if [ ! -d $(BINDIR) ]; then mkdir $(BINDIR) ; fi ; mv $(PROG) $(BINDIR)

# Windoze
winstall:
	if [ ! -d $(WBINDIR) ]; then mkdir $(WBINDIR) ; fi ; if [ -f $(WPROG) ]; then mv $(WPROG) $(PROG) ; fi ; mv $(PROG) $(WBINDIR)

clean:
	rm -f $(OBJS) $(PROG) $(WPROG) core

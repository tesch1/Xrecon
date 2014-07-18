/* default2D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* default2D.c: default 2D recon                                             */
/*                                                                           */
/* Copyright (C) 2009 Paul Kinchesh                                          */
/*                                                                           */
/* This file is part of Xrecon.                                              */
/*                                                                           */
/* Xrecon is free software: you can redistribute it and/or modify            */
/* it under the terms of the GNU General Public License as published by      */
/* the Free Software Foundation, either version 3 of the License, or         */
/* (at your option) any later version.                                       */
/*                                                                           */
/* Xrecon is distributed in the hope that it will be useful,                 */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the              */
/* GNU General Public License for more details.                              */
/*                                                                           */
/* You should have received a copy of the GNU General Public License         */
/* along with Xrecon. If not, see <http://www.gnu.org/licenses/>.            */
/*                                                                           */
/*---------------------------------------------------------------------------*/
/**/

#include "../Xrecon.h"

#define SOURCEFILE "recon2D/default2D.c"

void default2D(struct data *d)
{

#ifdef DEBUG
  char function[20];
  strcpy(function,"default2D"); /* Set function name */
#endif

  dimorder2D(d);                /* Sort ascending slice and phase order */

  setnvols2D(d);                /* Set the number of data volumes */

  /* Loop over requested volumes in data file */
  for (d->vol=d->startvol;d->vol<d->endvol;d->vol++) {

    /* Loop over data blocks */
    for (d->block=0;d->block<d->nblocks;d->block++) {

      getblock2D(d,d->vol,NDCC);    /* Get block without applying dbh.lvl and dbh.tlt */

      w2Dfdfs(d,VJ,FLT32,d->vol);   /* Write 2D fdf raw data from volume */

      getmax2D(d);                  /* Get coordinates of maximum */

      shiftdata2D(d,OPT);           /* Shift FID data for fft */

      equalizenoise2D(d,STD);       /* Scale for equal noise in all receivers */

      phaseramp2D(d,PHASE);         /* Phase ramp the data to correct for phase encode offset ppe */

      weightdata2D(d,STD);          /* Weight data using standard VnmrJ parameters */

      zerofill2D(d,STD);            /* Zero fill data using standard VnmrJ parameters */

      fft2D(d);                     /* 2D fft */

      phasedata2D(d,VJ);            /* Phase data if required */

      shiftdata2D(d,STD);           /* Shift data to get images */

      w2Dfdfs(d,VJ,FLT32,d->vol);   /* Write 2D fdf image data from volume */

      wnifti(d,VJ,FLT32,d->vol);    /* Write NIFTI-1/Analyze7.5 data */

      clear2Ddata(d);               /* Clear data volume from memory */

    }

  }

  clear2Dall(d);                    /* Clear everything from memory */

}

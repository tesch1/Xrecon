/* default3D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* default3D.c: default 3D recon                                             */
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

#define SOURCEFILE "recon3D/default3D.c"

void default3D(struct data *d)
{
  int MRCV,MBLK;

#ifdef DEBUG
  char function[20];
  strcpy(function,"default3D"); /* Set function name */
#endif

  MRCV = (d->nr > 1);           /* Multiple receivers */
  MBLK = (d->nblocks > 1);      /* Process in multiple data blocks */

  dimorder3D(d);                /* Sort ascending slice and phase order */

  setnvols3D(d);                /* Set the number of data volumes */

  if (MBLK) { /* Process in multiple data blocks */

    /* Loop over requested volumes in data file */
    for (d->vol=d->startvol;d->vol<d->endvol;d->vol++) {

      zeromax(d);                   /* Zero max structure & coordinates of maximum */

      zeronoise(d);                 /* Zero values in noise structure */

      /* Loop over data blocks */
      for (d->block=0;d->block<d->nblocks;d->block++) {

        getblock3D(d,d->vol,NDCC);  /* Get block without applying dbh.lvl and dbh.tlt */

        getmax(d);                  /* Get coordinates of maximum */

        if (MRCV) getnoise(d,STD);  /* Sample noise data if multiple receivers */

        w3Dfdfs(d,VJ,FLT32,d->vol); /* Write 3D fdf raw data from processing block */

        clear2Ddata(d);             /* Clear data processing block from memory */

      }

      /* Loop over data blocks */
      for (d->block=0;d->block<d->nblocks;d->block++) {

        getblock3D(d,d->vol,NDCC);  /* Get block without applying dbh.lvl and dbh.tlt */

        shiftdata2D(d,OPT);         /* Shift FID data for 2D dim1*dim2 fft */

        if (MRCV) equalizenoise2D(d,STD); /* Scale for equal noise in all receivers */

        phaseramp2D(d,PHASE);       /* Phase ramp the data in dim1*dim2 to correct for phase encode offset ppe */

        weightdata2D(d,STD);        /* Weight data in dim1*dim2 using standard VnmrJ parameters */

        zerofill2D(d,STD);          /* Zero fill data in dim1*dim2 using standard VnmrJ parameters */

        fft2D(d);                   /* 2D dim1*dim2 fft */

        phasedata2D(d,VJ);          /* Phase data in dim1*dim2 if required */

        shiftdata2D(d,STD);         /* Shift data in dim1*dim2 to get 2D images */

        wrawbin3D(d,D12,CX,DBL64);  /* Write complex raw data from "D12" processing block */

        clear2Ddata(d);             /* Clear data processing block from memory */

      }

      /* Loop over data blocks */
      for (d->block=0;d->block<d->nblocks;d->block++) {

        rrawbin3D(d,D12,CX,DBL64);  /* Read complex "D12" processing block into "D3" processing block */

        shiftdata1D(d,OPT,D3);      /* Shift FID data for 1D dim3 fft */

        phaseramp1D(d,PHASE2);      /* Phase ramp the data to correct for phase encode 2 offset ppe2 */

        weightdata1D(d,STD,D3);     /* Weight data in dim3 using standard VnmrJ parameters */

        zerofill1D(d,STD,D3);       /* Zero fill data in dim3 using standard VnmrJ parameters */

        fft1D(d,D3);                /* 1D dim3 fft */

        phasedata1D(d,VJ,D3);       /* Phase data in dim3 if required */

        shiftdata1D(d,STD,D3);      /* Shift data in dim3 to get image */

        wrawbin3D(d,D3,CX,DBL64);   /* Write complex raw data from "D3" processing block */

        cleardim3data(d);           /* Clear data processing block from memory */

      }

      delrawbin3D(d,D12,CX);        /* Delete the raw data we have finished with */

      /* Loop over data blocks */
      for (d->block=0;d->block<d->nblocks;d->block++) {

        rrawbin3D(d,D3,CX,DBL64);   /* Read complex "D3" processing block into "D12" processing block */

        w3Dfdfs(d,VJ,FLT32,d->vol); /* Write 3D fdf image data from volume */

        wnifti(d,VJ,FLT32,d->vol);  /* Write NIFTI-1/Analyze7.5 data */

        clear2Ddata(d);             /* Clear data volume from memory */

      }

      delrawbin3D(d,D3,CX);         /* Delete the raw data we have finished with */

    } /* end of volume for loop */

  } else { /* Process in a single data block */

    /* Loop over requested volumes in data file */
    for (d->vol=d->startvol;d->vol<d->endvol;d->vol++) {

      zeromax(d);                   /* Zero max structure & coordinates of maximum */

      zeronoise(d);                 /* Zero values in noise structure */

      getblock3D(d,d->vol,NDCC);    /* Get block without applying dbh.lvl and dbh.tlt */

      getmax(d);                    /* Get coordinates of maximum */

      if (MRCV) getnoise(d,STD);    /* Sample noise data if multiple receivers */

      w3Dfdfs(d,VJ,FLT32,d->vol);   /* Write 3D fdf raw data from processing block */

      shiftdata2D(d,OPT);           /* Shift FID data for 2D dim1*dim2 fft */

      if (MRCV) equalizenoise2D(d,STD); /* Scale for equal noise in all receivers */

      phaseramp2D(d,PHASE);         /* Phase ramp the data to correct for phase encode offset ppe */

      weightdata2D(d,STD);          /* Weight data in dim1*dim2 using standard VnmrJ parameters */

      zerofill2D(d,STD);            /* Zero fill data in dim1*dim2 using standard VnmrJ parameters */

      fft2D(d);                     /* 2D dim1*dim2 fft */

      phasedata2D(d,VJ);            /* Phase data in dim1*dim2 if required */

      shiftdata2D(d,STD);           /* Shift data in dim1*dim2 to get 2D images */

      shiftdatadim3(d,D12,OPT);     /* Shift FID data for 1D dim3 fft */

      phaserampdim3(d,PHASE2);      /* Phase ramp the data to correct for phase encode 2 offset ppe2 */

      weightdatadim3(d,STD);        /* Weight data in dim3 using standard VnmrJ parameters */

      zerofilldim3(d,STD);          /* Zero fill data in dim3 using standard VnmrJ parameters */

      fftdim3(d);                   /* 1D dim3 fft */

      phasedatadim3(d,VJ);          /* Phase data in dim3 if required */

      shiftdatadim3(d,D12,STD);     /* Shift data in dim3 to get image */

      w3Dfdfs(d,VJ,FLT32,d->vol);   /* Write 3D fdf image data from volume */

      wnifti(d,VJ,FLT32,d->vol);    /* Write NIFTI-1/Analyze7.5 data */

      clear2Ddata(d);               /* Clear data volume from memory */

    } /* end of volume for loop */

  } /* end process in a single data block */

  clear2Dall(d);                    /* Clear everything from memory */

}

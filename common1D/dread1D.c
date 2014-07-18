/* dread1D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dread1D.c: 1D Data reading routines                                       */
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

#define SOURCEFILE "common1D/dread1D.c"

void getvol1D(struct data *d,int volindex,int DCCflag)
{
  int dim1,dim2,dim3,nr;
  int i,j;
  int datatype;
  char function[20];
  strcpy(function,"getvol1D"); /* Set function name */

  /* Set datatype */
  datatype=0;
  if (d->fh.status & S_FLOAT) /* 32-bit float */
    datatype=FLT32;
  else if (d->fh.status & S_32) /* 32-bit int */
    datatype=INT32;
  else /* 16-bit int */
    datatype=INT16;  
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  datatype = %d\n",datatype);
  fflush(stdout);
#endif

  /* If a previous block has been zero filled matrix size will be incorrect */
  /* Refresh from procpar values */
  if (d->zerofill) setdatapars(d);

  /* Set data dimensions */
  dim1=d->np/2; dim2=1; dim3=d->fh.ntraces; nr=d->nr;

  /* Allocate memory according to nr for all apptypes */
  if ((d->data = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();

  /* Allocate memory according to nr and number of traces */
  if ((d->data = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();
  for (i=0;i<nr;i++) {
    if ((d->data[i] = (fftw_complex **)fftw_malloc(dim3*sizeof(fftw_complex *))) == NULL) nomem();
    for (j=0;j<dim3;j++) 
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(dim1*sizeof(fftw_complex))) == NULL) nomem();
  }
  /* Read according to d->type */
  switch(datatype) {
    case FLT32: /* 32-bit float */
      read1Df(d,volindex,DCCflag);
      break;
    case INT32: /* 32-bit int */
      read1Dl(d,volindex,DCCflag);
      break;
    case INT16: /* 16-bit int */
      read1Ds(d,volindex,DCCflag);
      break;
    default:
      break;
  } /* end datatype switch */
  /* Set datamode for fid data */
  d->datamode=FID;
  /* Update shifted flag */
  d->shift=FALSE; /* not shifted */
  /* Update zerofill flag */
  d->zerofill=FALSE; /* not zerofilled */
  /* Zero max measurement storage */
  zeromax(d);
  /* Zero noise measurement storage */
  zeronoise(d);

}

int read1Df(struct data *d,int volindex,int DCCflag)
{
  float *fdata;        /* Pointer for 32-bit floating point data */
  int dim1,dim2,dim3,nr;
  int i,j,k;

#ifdef DEBUG
  char function[20];
  strcpy(function,"read1Df"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 32-bit floating point data\n");
  fflush(stdout);
#endif    

  /* Data dimensions */
  dim1=d->np/2; dim2=1; dim3=d->fh.ntraces; nr=d->nr;

  /* Allocate memory */
  if ((fdata = (float *)malloc(d->fh.np*d->fh.ebytes)) == NULL) nomem();

#ifdef DEBUG
  if (DCCflag) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Volume %d:\n",volindex);
    fprintf(stdout,"  DC correction using dbh.lvl and dbh.tlt\n");
    fflush(stdout);
  }
#endif

 /* Read nr consecutive blocks for the 'volume' */
  for (i=0;i<nr;i++) { /* loop over receiver blocks */
    if (DCCflag) {
      getdbh(d,nr*volindex+i); /* Get block header */
    }
    for (j=0;j<dim3;j++) { /* loop over traces */
      /* Set data pointer */
      if (set1Doffset(d,volindex,i,j)) {
        free(fdata);
        return(1);
      }
      fread(fdata,d->fh.ebytes,d->fh.np,d->fp);
      if (reverse_byte_order)
        reverse4ByteOrder(d->fh.np,(char *)fdata);
      for (k=0;k<dim1;k++) {
        d->data[i][j][k][0]=(double)fdata[2*k];
        d->data[i][j][k][1]=(double)fdata[2*k+1];
      }
      if (DCCflag) {
        for (k=0;k<dim1;k++) {
          d->data[i][j][k][0]-=(double)d->bh.lvl;
          d->data[i][j][k][1]-=(double)d->bh.tlt;
        }
      }
    } /* end trace loop */
  } /* end receiver block loop */

  /* Free memory */
  free(fdata);
  return(0);
}

int read1Dl(struct data *d,int volindex,int DCCflag)
{
  long int *ldata;     /* Pointer for 32-bit integer data */
  int dim1,dim2,dim3,nr;
  int i,j,k;

#ifdef DEBUG
  char function[20];
  strcpy(function,"read1Dl"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 32-bit integer data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=1; dim3=d->fh.ntraces; nr=d->nr;

  /* Allocate memory */
  if ((ldata = (long *)malloc(d->fh.np*d->fh.ebytes)) == NULL) nomem();

#ifdef DEBUG
  if (DCCflag) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Volume %d:\n",volindex);
    fprintf(stdout,"  DC correction using dbh.lvl and dbh.tlt\n");
    fflush(stdout);
  }
#endif

 /* Read nr consecutive blocks for the 'volume' */
  for (i=0;i<nr;i++) { /* loop over receiver blocks */
    if (DCCflag) {
      getdbh(d,nr*volindex+i); /* Get block header */
    }
    for (j=0;j<dim3;j++) { /* loop over traces */
      /* Set data pointer */
      if (set1Doffset(d,volindex,i,j)) {
        free(ldata);
        return(1);
      }
      fread(ldata,d->fh.ebytes,d->fh.np,d->fp);
      if (reverse_byte_order)
        reverse4ByteOrder(d->fh.np,(char *)ldata);
      for (k=0;k<dim1;k++) {
        d->data[i][j][k][0]=(double)ldata[2*k];
        d->data[i][j][k][1]=(double)ldata[2*k+1];
      }
      if (DCCflag) {
        for (k=0;k<dim1;k++) {
          d->data[i][j][k][0]-=(double)d->bh.lvl;
          d->data[i][j][k][1]-=(double)d->bh.tlt;
        }
      }
    } /* end trace loop */
  } /* end receiver block loop */

  /* Free memory */
  free(ldata);
  return(0);
}

int read1Ds(struct data *d,int volindex,int DCCflag)
{
  short int *sdata;    /* Pointer for 16-bit integer data */
  int dim1,dim2,dim3,nr;
  int i,j,k;

#ifdef DEBUG
  char function[20];
  strcpy(function,"read1Ds"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 16-bit integer data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=1; dim3=d->fh.ntraces; nr=d->nr;

  /* Allocate memory */
  if ((sdata = (short *)malloc(d->fh.np*d->fh.ebytes)) == NULL) nomem();

#ifdef DEBUG
  if (DCCflag) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Volume %d:\n",volindex);
    fprintf(stdout,"  DC correction using dbh.lvl and dbh.tlt\n");
    fflush(stdout);
  }
#endif

 /* Read nr consecutive blocks for the 'volume' */
  for (i=0;i<nr;i++) { /* loop over receiver blocks */
    if (DCCflag) {
      getdbh(d,nr*volindex+i); /* Get block header */
    }
    for (j=0;j<dim3;j++) { /* loop over traces */
      /* Set data pointer */
      if (set1Doffset(d,volindex,i,j)) {
        free(sdata);
        return(1);
      }
      fread(sdata,d->fh.ebytes,d->fh.np,d->fp);
      if (reverse_byte_order)
        reverse2ByteOrder(d->fh.np,(char *)sdata);
      for (k=0;k<dim1;k++) {
        d->data[i][j][k][0]=(double)sdata[2*k];
        d->data[i][j][k][1]=(double)sdata[2*k+1];
      }
      if (DCCflag) {
        for (k=0;k<dim1;k++) {
          d->data[i][j][k][0]-=(double)d->bh.lvl;
          d->data[i][j][k][1]-=(double)d->bh.tlt;
        }
      }
    } /* end trace loop */
  } /* end receiver block loop */

  /* Free memory */
  free(sdata);
  return(0);
}

int set1Doffset(struct data *d,int volindex,int receiver,int trace)
{
  int dim1,dim2,dim3,nr;
  long offset=0;
  int blockindex;

#ifdef DEBUG
  char function[20];
  strcpy(function,"set1Doffset"); /* Set function name */
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=1; dim3=d->fh.ntraces; nr=d->nr;

  blockindex=volindex*nr+receiver;

  /* Calculate offset to the required trace */
  offset=sizeof(d->fh)+(long)(blockindex)*d->fh.bbytes+sizeof(d->bh)+d->fh.tbytes*(long)(trace);

  /* Set offset */
  if (d->buf.st_size >= offset+d->fh.tbytes) fseek(d->fp,offset,SEEK_SET);
  else {
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Unable to set offset for volume %d receiver %d trace %d\n",
    volindex,receiver,trace);
  fflush(stdout);
#endif
    return(1);
  }
  return(0);
}

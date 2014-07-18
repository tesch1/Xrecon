/* dread2D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dread2D.c: 2D Data read routines                                          */
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

#define SOURCEFILE "common2D/dread2D.c"

static int currentblockhead=-1;

void getblock2D(struct data *d,int volindex,int DCCflag)
{
  int dim1,dim2,dim3,nr;
  int i,j;
  int datatype;
  char function[20];
  strcpy(function,"getblock2D"); /* Set function name */

  /* Check d->nvols has been set */
  if (d->nvols == 0) setnvols2D(d); /* Set the number of data volumes */

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

  /* Set start and end slice of block */
  setblock2D(d);

  /* Set data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Allocate memory according to nr for all apptypes */
  if ((d->data = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();

  /* Data from each slice of the 'volume' requires 2D fft */
  /* Allocate memory according to number of slices to be processed */
  for (i=0;i<nr;i++) { /* loop over receiver blocks */
    if ((d->data[i] = (fftw_complex **)fftw_malloc(dim3*sizeof(fftw_complex *))) == NULL) nomem();
    for (j=0;j<dim3;j++) /* loop over slices */
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();
  }

  /* Read according to d->type */
  switch(datatype) {
    case FLT32: /* 32-bit float */
      read2Df(d,volindex,DCCflag);
      break;
    case INT32: /* 32-bit int */
      read2Dl(d,volindex,DCCflag);
      break;
    case INT16: /* 16-bit int */
      read2Ds(d,volindex,DCCflag);
      break;
    default:
      break;
  } /* end datatype switch */
  /* Reset blockhead read for next volume */
  currentblockhead=-1;
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

void getvol2D(struct data *d,int volindex,int DCCflag)
{
  int dim1,dim2,ns,nr;
  int i,j;
  int datatype;
  char function[20];
  strcpy(function,"getvol2D"); /* Set function name */

  /* Check d->nvols has been set */
  if (d->nvols == 0) setnvols2D(d); /* Set the number of data volumes */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; ns=d->ns; nr=d->nr;

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

  /* If a previous volume has been zero filled matrix size will be incorrect */
  /* Refresh from procpar values */
  if (d->zerofill) setdatapars2D(d);

  /* Allocate memory according to nr for all apptypes */
  if ((d->data = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();

  /* Data from each slice of the 'volume' requires 2D fft */
  /* Allocate memory according to number of slices */
  for (i=0;i<nr;i++) { /* loop over receiver blocks */
    if ((d->data[i] = (fftw_complex **)fftw_malloc(ns*sizeof(fftw_complex *))) == NULL) nomem();
    for (j=0;j<ns;j++) /* loop over slices */
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();
  }
  /* Read according to d->type */
  switch(datatype) {
    case FLT32: /* 32-bit float */
      read2Df(d,volindex,DCCflag);
      break;
    case INT32: /* 32-bit int */
      read2Dl(d,volindex,DCCflag);
      break;
    case INT16: /* 16-bit int */
      read2Ds(d,volindex,DCCflag);
      break;
    default:
      break;
  } /* end datatype switch */
  /* Reset blockhead read for next volume */
  currentblockhead=-1;
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

int read2Df(struct data *d,int volindex,int DCCflag)
{
  float *fdata;        /* Pointer for 32-bit floating point data */
  int dim1,dim2,dim3,nr;
  int startslice,slice,scale;
  int i,j,k,l;

#ifdef DEBUG
  char function[20];
  strcpy(function,"read2Df"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 32-bit floating point data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Start slice */
  startslice=d->startpos;

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
    for (j=0;j<dim3;j++) { /* loop over slices */
      slice=startslice+j;
      for (k=0;k<dim2;k++) { /* loop over phase */
        /* Set data pointer */
        if (set2Doffset(d,volindex,i,slice,k)) {
          free(fdata);
          return(1);
        }
        fread(fdata,d->fh.ebytes,d->fh.np,d->fp);
        if (reverse_byte_order)
          reverse4ByteOrder(d->fh.np,(char *)fdata);
        scale=d->bh.ctcount;
        if (scale<1) scale=1; /* don't want to divide by 0 if there's no data */
        for (l=0;l<dim1;l++) {
          d->data[i][j][k*dim1+l][0]=(double)fdata[2*l]/scale;
          d->data[i][j][k*dim1+l][1]=(double)fdata[2*l+1]/scale;
        }
        if (DCCflag) {
          for (l=0;l<dim1;l++) {
            d->data[i][j][k*dim1+l][0]-=(double)d->bh.lvl/scale;
            d->data[i][j][k*dim1+l][1]-=(double)d->bh.tlt/scale;
          }
        }
      } /* end phase loop */
    } /* end slice loop */
  } /* end receiver block loop */

  /* Free memory */
  free(fdata);
  return(0);
}

int read2Dl(struct data *d,int volindex,int DCCflag)
{
  long int *ldata;     /* Pointer for 32-bit integer data */
  int dim1,dim2,dim3,nr;
  int startslice,slice,scale;
  int i,j,k,l;

#ifdef DEBUG
  char function[20];
  strcpy(function,"read2Dl"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 32-bit integer data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Start slice */
  startslice=d->startpos;

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
    for (j=0;j<dim3;j++) { /* loop over slices */
      slice=startslice+j;
      for (k=0;k<dim2;k++) { /* loop over phase */
        /* Set data pointer */
        if (set2Doffset(d,volindex,i,slice,k)) {
          free(ldata);
          return(1);
        }
        fread(ldata,d->fh.ebytes,d->fh.np,d->fp);
        if (reverse_byte_order) 
          reverse4ByteOrder(d->fh.np,(char *)ldata);
        scale=d->bh.ctcount;
        if (scale<1) scale=1; /* don't want to divide by 0 if there's no data */
        for (l=0;l<dim1;l++) {
          d->data[i][j][k*dim1+l][0]=(double)ldata[2*l]/scale;
          d->data[i][j][k*dim1+l][1]=(double)ldata[2*l+1]/scale;
        }
        if (DCCflag) {
          for (l=0;l<d->fh.np/2;l++) {
            d->data[i][j][k*dim1+l][0]-=(double)d->bh.lvl/scale;
            d->data[i][j][k*dim1+l][1]-=(double)d->bh.tlt/scale;
          }
        }
      } /* end phase loop */
    } /* end slice loop */
  } /* end receiver block loop */

  /* Free memory */
  free(ldata);
  return(0);
}

int read2Ds(struct data *d,int volindex,int DCCflag)
{
  short int *sdata;    /* Pointer for 16-bit integer data */
  int dim1,dim2,dim3,nr;
  int startslice,slice,scale;
  int i,j,k,l;

#ifdef DEBUG
  char function[20];
  strcpy(function,"read2Ds"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 16-bit integer data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Start slice */
  startslice=d->startpos;

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
    for (j=0;j<dim3;j++) { /* loop over slices */
      slice=startslice+j;
      for (k=0;k<dim1;k++) { /* loop over phase */
        /* Set data pointer */
        if (set2Doffset(d,volindex,i,slice,k)) {
          free(sdata);
          return(1);
        }
        fread(sdata,d->fh.ebytes,d->fh.np,d->fp);
        if (reverse_byte_order) 
          reverse2ByteOrder(d->fh.np,(char *)sdata);
        scale=d->bh.ctcount;
        if (scale<1) scale=1; /* don't want to divide by 0 if there's no data */
        for (l=0;l<dim1;l++) {
          d->data[i][j][k*dim1+l][0]=(double)sdata[2*l]/scale;
          d->data[i][j][k*dim1+l][1]=(double)sdata[2*l+1]/scale;
        }
        if (DCCflag) {
          for (l=0;l<dim1;l++) {
            d->data[i][j][k*dim1+l][0]-=(double)d->bh.lvl/scale;
            d->data[i][j][k*dim1+l][1]-=(double)d->bh.tlt/scale;
          }
        }
      } /* end phase loop */
    } /* end slice loop */
  } /* end receiver block loop */

  /* Free memory */
  free(sdata);
  return(0);
}

int set2Doffset(struct data *d,int volindex,int receiver,int slice,int phase)
{
  int dim1,dim2,ns,nr;
  int index;
  int arraydim,ne,etl,psscycle,seqmode;
  long offset=0;
  int blockindex=0;
  char seqcon[6];

#ifdef DEBUG
  char function[20];
  strcpy(function,"set2Doffset"); /* Set function name */
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; ns=d->ns; nr=d->nr;

  /* Number of echoes */
  ne=(int)*val("ne",&d->p);
  if (ne < 1) ne=1; /* Set ne to 1 if 'ne' does not exist */

  /* Allow for compressed multi-echo loop */
  index=volindex/ne;

  /* Set slice according to its actual index */
  slice=sliceindex(d,slice);

  /* Set phase according to its actual index */
  phase=phaseindex(d,phase);

  /* Set sequence mode */
  seqmode=0;
  /* Standard cases */
  strcpy(seqcon,*sval("seqcon",&d->p));
  if ((seqcon[1] == 'c') && (seqcon[2] == 'c')) seqmode=1;
  if ((seqcon[1] == 'c') && (seqcon[2] == 's')) seqmode=2;
  if ((seqcon[1] == 's') && (seqcon[2] == 'c')) seqmode=3;
  if ((seqcon[1] == 's') && (seqcon[2] == 's')) seqmode=4;
  /* Special cases */
  if (spar(d,"apptype","im2Dfse")) seqmode=5;
  if (spar(d,"apptype","im2Depi")) seqmode=6;

  /* Calculate offset according to how the expt is run */
  switch(seqmode) {
    case 1: /* seqcon = *cc** */
      blockindex=index*nr+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(phase*ns*ne+slice*ne+volindex%ne);
      break;
    case 2: /* seqcon = *cs** */
      /* The standard phase encode loop cycles more slowly than any
         array elements */
      arraydim=(int)*val("arraydim",&d->p);
      blockindex=index*nr+phase*arraydim/dim2+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(slice*ne+volindex%ne);
      break;
    case 3: /* seqcon = *sc** */
      /* Somewhat trickier as we must use 'psscycle' to account for where 
         and how 'pss' appears amongst the array elements */
      if (ns > 1) { /* there is more than one slice */
        psscycle=getcycle("pss",&d->a);
        blockindex=(index/psscycle)*psscycle*ns + index%psscycle;
        blockindex=blockindex*nr+slice*psscycle*nr+receiver;
      } else { /* there's a single slice */
        blockindex=index;
      }
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(phase*ne+volindex%ne);
      break;
    case 4: /* seqcon = *ss** */
      /* Even trickier. The standard phase encode loop cycles more slowly 
         than any array elements and we must use 'psscycle' to account for 
         where and how 'pss' appears amongst the array elements */
      arraydim=(int)*val("arraydim",&d->p);
      if (ns > 1) { /* there is more than one slice */
        psscycle=getcycle("pss",&d->a);
        blockindex=(index/psscycle)*psscycle*ns + index%psscycle;
        blockindex=blockindex*nr+slice*psscycle*nr+phase*arraydim/dim2+receiver;
      } else { /* there's a single slice */
        blockindex=index*nr+phase*arraydim/dim2+receiver;
      }
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(volindex%ne);
      break;
    case 5: /* im2Dfse */
      arraydim=(int)*val("arraydim",&d->p);
      etl=(int)*val("etl",&d->p);
      if (seqcon[2] == 'c') { /* seqcon = *cc** */
        blockindex=index*nr+receiver;
        /* Calculate offset to end of the required block */
        offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
        /* Adjust for required trace */
        offset-=d->fh.ntraces*d->fh.tbytes;
        offset+=d->fh.tbytes*(long)((phase/etl)*ns*etl+slice*etl+phase%etl);
      } else { /* seqcon = *cs** */
        blockindex=index*nr+(phase/etl)*(arraydim*etl)/dim2+receiver;
        /* Calculate offset to end of the required block */
        offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
        /* Adjust for required trace */
        offset-=d->fh.ntraces*d->fh.tbytes;
        offset+=d->fh.tbytes*(long)(slice*etl+phase%etl);
      }
      break;
    case 6: /* imEPI */
      blockindex=index*nr+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(ns*phase+slice);
      break;
    default:
      break;
  } /* end seqmode switch */

  /* Get block header if required */
  if (blockindex != currentblockhead) {
    if (blockindex<d->fh.nblocks) {
      getdbh(d,blockindex);
      currentblockhead=blockindex;
    } else {
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Unable to read block header %d\n",blockindex);
  fflush(stdout);
#endif
    }
  }

  /* Set offset */
  if (d->buf.st_size >= offset+d->fh.tbytes) fseek(d->fp,offset,SEEK_SET);
  else {
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Unable to set offset for volume %d receiver %d slice %d phase %d\n",
    volindex,receiver,slice,phase);
  fflush(stdout);
#endif
    return(1);
  }
  return(0);
}

int sliceindex(struct data *d,int slice)
{
  int i;

  if (d->pssorder[0] == -1) return(slice);

  for (i=0;i<d->ns;i++)
    if (d->pssorder[i] == slice) return(i);

  return(-1);
}

int phaseindex(struct data *d,int phase)
{
  int i;

  if (d->dim2order[0] == -1) return(phase);

  for (i=0;i<d->nv;i++)
    if (d->dim2order[i] == phase) return(i);

  return(-1);
}

int phase2index(struct data *d,int phase2)
{
  int i;

  if (d->dim3order[0] == -1) return(phase2);

  for (i=0;i<d->nv2;i++)
    if (d->dim3order[i] == phase2) return(i);

  return(-1);
}

void setnvols2D(struct data *d)
{
  int ns,nr;
  long nvols;
  long volbytes,databytes;
  int ne,etl;
  char seqcon[6];
  char function[20];
  strcpy(function,"setnvols2D"); /* Set function name */

  /* To figure the number of volumes insist we must have the following */
  if ((ptype("np",&d->p) < 0) || (ptype("nv",&d->p) < 0) 
    || (ptype("pss",&d->p) < 0) || (ptype("rcvrs",&d->p) < 0)
    || (ptype("seqcon",&d->p) < 0) || (ptype("arraydim",&d->p) < 0)
    || (ptype("nf",&d->p) < 0)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  The following values are required from a 'procpar':\n");
    fprintf(stdout,"  np nv pss rcvrs seqcon arraydim nf\n\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }

  /* Some additional parameters may be set */
  ne=(int)*val("ne",&d->p);
  if (ne < 1) ne=1; /* Set ne to 1 if 'ne' does not exist */
  etl=(int)*val("etl",&d->p);
  if (etl < 1) etl=1; /* Set etl to 1 if 'etl' does not exist */

  /* Set data parameters from procpar values */
  setdatapars2D(d);

  /* Set the sequence mode from seqcon and apptype parameters */
  setseqmode(d);

  /* Number of slices and receivers */
  ns=d->ns; nr=d->nr;

  /* Sanity check */
  if ((d->np != d->fh.np) || ((int)*val("nf",&d->p) != d->fh.ntraces)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Check data file and procpar file are compatible\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }

  /* Use value of arraydim to calculate the total number of volumes */
  nvols=(long)*val("arraydim",&d->p);
  /* Multiple echoes generate 'ne' volumes that are not arrayed */
  nvols*=ne;
  /* Multiple receivers generate 'nr' volumes included in arraydim */
  if (nvols%nr != 0) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Problem calculating number of volumes\n");
    fprintf(stdout,"  Check data file and procpar file are compatible\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }
  nvols/=nr;

  /* Check for 'standard' slice and phase loops as these
     generate volumes included in arraydim */
  strcpy(seqcon,*sval("seqcon",&d->p));
  if (seqcon[1] == 's') /* 'standard' slice loop has pss arrayed */
    nvols/=ns;
  if (seqcon[2] == 's') { /* 'standard' phase loop */
    nvols*=etl;
    nvols/=d->nv;
  }
  /* We now have the total number of volumes (nvols) according to
     the input parameter set */

  /* So lets check with a byte count ... */

  /* Calculate the number of bytes in a volume */
  volbytes=(long)d->np*d->fh.ebytes; /* a single trace */
  if (d->nv > 0) volbytes*=d->nv;    /* a single slice */
  volbytes*=ns;                      /* a single volume */
  volbytes*=nr;                      /* Include data from all receivers */

  /* Calculate the total number of data bytes in the file */
  databytes=(long)d->fh.nblocks*d->fh.tbytes*d->fh.ntraces;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Number of volumes from arraydim = %ld\n",nvols);
  fprintf(stdout,"  Calculated volume bytes = %ld\n",volbytes);
  fprintf(stdout,"  Actual data bytes = %ld\n",databytes);
  fprintf(stdout,"  Number of volumes according to byte count = %f\n",
    (float)databytes/volbytes);
  fflush(stdout);
#endif

  if ((databytes%volbytes != 0) 
    || (nvols != databytes/volbytes)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Hmmm - something is not quite right\n");
    fflush(stdout);
  }

  /* Set nvols */
  d->nvols=nvols;

  /* Set start volume and end volume */
  if (!strcmp(*sval("allvolumes",&d->p),"n")) {  /* Don't process all volumes */
    d->startvol=(int)*val("startvol",&d->p)-1;
    d->endvol=(int)*val("endvol",&d->p);
    if (d->startvol > d->endvol) {
      /* Swap them */
      d->vol=d->endvol; d->endvol=d->startvol; d->startvol=d->vol;
    }
    if (d->startvol < 0) d->startvol=0;
    if (d->endvol > d->nvols) d->endvol=d->nvols;
  } else { /* Process all volumes */
    d->startvol=0;
    d->endvol=d->nvols;
  }

  /* Initialise volume counter */
  d->vol=d->startvol;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  d->nvols = %d\n",d->nvols);
  fprintf(stdout,"  d->startvol = %d\n",d->startvol);
  fprintf(stdout,"  d->endvol = %d\n",d->endvol);
  fflush(stdout);
#endif
}

void setblock2D(struct data *d)
{
  int blockns;

#ifdef DEBUG
  char function[20];
  strcpy(function,"setblock2D"); /* Set function name */
#endif

  /* Figure how many slices are required per block */
  blockns=d->ns/d->nblocks;
  if (d->ns%d->nblocks) blockns++;

  /* Set start slice and end slice */
  d->startpos=d->block*blockns;
  d->endpos=(d->block+1)*blockns;
  if (d->endpos>d->ns) d->endpos=d->ns;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Block %d (of %d)\n",d->block+1,d->nblocks);
  fprintf(stdout,"  Number of slices in block = %d (slice %d to slice %d)\n",d->endpos-d->startpos,d->startpos+1,d->endpos);
  fflush(stdout);
#endif
}

int process2Dblock(struct data *d,char *startpar,char *endpar)
{
  int blockns;
  int startpos,endpos;
  int start,end;

#ifdef DEBUG
  char function[20];
  strcpy(function,"process2Dblock"); /* Set function name */
#endif

  /* Number of pss values = slices */
  d->ns=nvals("pss",&d->p);

  /* There must be at least one block per volume */
  d->nblocks=(int)*val("nblocks",&d->p);
  if (!d->nblocks) d->nblocks++;

  /* There can not be more blocks than slices */
  if (d->nblocks>d->ns) d->nblocks=d->ns;

  /* Figure how many slices are required per block */
  blockns=d->ns/d->nblocks;
  if (d->ns%d->nblocks) blockns++;

  /* Figure start slice and end slice of the block */
  startpos=d->block*blockns;
  endpos=(d->block+1)*blockns;
  if (endpos>d->ns) endpos=d->ns;

  /* Get start slice and end slice to be processed */
  start=(int)*val(startpar,&d->p);
  if ((start<1) || (start>d->ns)) start=1;

  end=(int)*val(endpar,&d->p);
  if ((end<1) || (end>d->ns)) end=d->ns;

  /* Return whether to process or not */
  if ((start>endpos) || (end<=startpos)) {
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Skipping processing of block %d (of %d)\n",d->block+1,d->nblocks);
#endif
    return(FALSE);
  }

  return(TRUE);

}

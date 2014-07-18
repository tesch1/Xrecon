/* dproc.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dproc.c: Data processing routines                                         */
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

#define SOURCEFILE "common/dproc.c"

static int currentblockhead=-1;

/*-----------------------*/
/*---- Some defaults ----*/
/*-----------------------*/
/* Knoisefraction: Fraction of k-space FOV to use to sample noise */
static double Knoisefraction=0.05;

/* IMnoisefraction: Fraction of image space FOV to use to sample noise */
static double IMnoisefraction=0.05;

int *sliceorder(struct data *d,int dim,char *par)
{
  int *dimorder;
  double *pss;
  int index;
  int standardorder;
  int i,j;


#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"sliceorder"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Set dimorder=-1 if par is not found */
  standardorder=TRUE;
  index=parindex(par,&d->p);
  if (index < 0) {
    if ((dimorder = (int *)malloc(sizeof(int))) == NULL) nomem();
    dimorder[0]=-1;
  } else {
    /* Take a copy of the starting slice order */
    if ((pss = (double *)malloc(dim*sizeof(double))) == NULL) nomem();
    for (i=0;i<dim;i++) pss[i]=d->p.d[index][i];
    /* Sort slice positions into ascending order */
    qsort(pss,dim,sizeof(pss[0]),doublecmp);
    /* Fill dimorder with an index of where successive slices are */
    if ((dimorder = (int *)malloc(dim*sizeof(int))) == NULL) nomem();
    for (i=0;i<dim;i++) {
      for (j=0;j<dim;j++) {
        if (d->p.d[index][i] == pss[j]) {
          dimorder[i]=j;
          break;
        }
      }
    }
    /* Set dimorder=-1 if we simply have standard order */
    for (i=0;i<dim;i++) {
      if (dimorder[i] != i) {
        standardorder=FALSE;
        break;
      }
    }
    if (standardorder) {
      free(dimorder);
      if ((dimorder = (int *)malloc(sizeof(int))) == NULL) nomem();
      dimorder[0]=-1;
    }

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Slice ordering took %f secs\n",t2-t1);
  if (!standardorder) {
    fprintf(stdout,"            pss\n");
    fprintf(stdout,"  original\tnew\t    dimorder\n");
    for (i=0;i<dim;i++)
      fprintf(stdout,"  %f\t%f\t%d\n",d->p.d[index][i],pss[i],dimorder[i]);
  }
  fflush(stdout);
#endif

    free(pss);
  }

  return(dimorder);

}

int *phaseorder(struct data *d,int dim,char *par)
{
  int *dimorder;
  int index,min;
  int standardorder;
  int i,j;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"phaseorder"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Set dimorder=-1 if par is not found or has wrong size */
  standardorder=TRUE;
  index=parindex(par,&d->p);
  if ((index < 0) || (nvals(par,&d->p) != dim)) {
    if ((dimorder = (int *)malloc(sizeof(int))) == NULL) nomem();
    dimorder[0]=-1;
  } else {
    /* Fill dimorder with an index of where successive phase encodes are */
    if ((dimorder = (int *)malloc(dim*sizeof(int))) == NULL) nomem();
    /* Find the minimum value (phase encode multipliers can run from -nv/2 to +nv/2-1 or -nv/2+1 to +nv/2) */
    min=0;
    for (i=0;i<dim;i++) if (d->p.d[index][i]<min) min=d->p.d[index][i];
    for (i=0;i<dim;i++) {
      for (j=0;j<dim;j++) {
        if (d->p.d[index][i] == j+min) {
          dimorder[i]=j;
          break;
        }
      }
    }
    /* Set dimorder=-1 if we simply have standard order */
    for (i=0;i<dim;i++) {
      if (dimorder[i] != i) {
        standardorder=FALSE;
        break;
      }
    }
    if (standardorder) {
      free(dimorder);
      if ((dimorder = (int *)malloc(sizeof(int))) == NULL) nomem();
      dimorder[0]=-1;
    }
  }

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Phase ordering took %f secs\n",t2-t1);
  if (!standardorder) {
  fprintf(stdout,"  dimorder:\n");
    for (i=0;i<dim;i++) {
      fprintf(stdout,"%5.1d",dimorder[i]);
      if ((i+1)%16 == 0) fprintf(stdout,"\n");
    }
  }
  fflush(stdout);
#endif

  return(dimorder);

}

void setnvols(struct data *d)
{
  int nvols;
  int ne,etl;
  char seqcon[6];
  char function[20];
  strcpy(function,"setnvols"); /* Set function name */

  /* To figure the number of volumes insist we must have the following */
  if ((ptype("np",&d->p) < 0) || (ptype("nv",&d->p) < 0)
    || (ptype("nv",&d->p) < 0) || (ptype("pss",&d->p) < 0) 
    || (ptype("rcvrs",&d->p) < 0) || (ptype("seqcon",&d->p) < 0) 
    || (ptype("arraydim",&d->p) < 0) || (ptype("nf",&d->p) < 0)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  The following values are required from a 'procpar':\n");
    fprintf(stdout,"  np nv nv2 pss rcvrs seqcon arraydim nf\n\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }

  /* Some additional parameters may be set */
  ne=(int)*val("ne",&d->p);
  if (ne < 1) ne=1; /* Set ne to 1 if 'ne' does not exist */
  etl=(int)*val("etl",&d->p);
  if (etl < 1) etl=1; /* Set etl to 1 if 'etl' does not exist */

  /* Set data parameters from procpar values */
  setdatapars(d);

  /* Set the sequence mode from seqcon and apptype parameters */
  setseqmode(d);

  /* Sanity check */
  if ((d->np != d->fh.np) || ((int)*val("nf",&d->p) != d->fh.ntraces)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Check data file and procpar file are compatible\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }

  /* Use value of arraydim to calculate the total number of volumes */
  nvols=(int)*val("arraydim",&d->p);
  /* Multiple echoes generate 'ne' volumes that are not arrayed */
  nvols*=ne;
  /* Multiple receivers generate 'nr' volumes included in arraydim */
  if (nvols%d->nr != 0) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Problem calculating number of volumes\n");
    fprintf(stdout,"  Check data file and procpar file are compatible\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }
  nvols/=d->nr;

  /* Check for 'standard' slice and phase loops as these
     generate volumes included in arraydim */
  strcpy(seqcon,*sval("seqcon",&d->p));
  if (d->seqmode<IM3D) {  /* 2D */
    /* 'standard' slice loop has pss arrayed */
    if (seqcon[1] == 's') nvols/=d->ns;
  } else { /* 3D */
    /* 'compressed' slice loop generates 'ns' extra volumes  */
    if (seqcon[1] == 'c') nvols*=d->ns;
  }
  if (seqcon[2] == 's') { /* 'standard' phase loop */
    nvols*=etl;
    nvols/=d->nv;
  }
  if (seqcon[3] == 's') { /* 'standard' phase loop */
    nvols/=d->nv2;
  }
  /* We now have the total number of volumes (nvols) according to
     the input parameter set */

  /* Set nvols */
  d->nvols=nvols;

  /* Set start volume and end volume */
  if (spar(d,"allvolumes","n")) {  /* Don't process all volumes */
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

  /* Initialise block counter */
  d->block=0;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Total number of volumes (d->nvols) = %d\n",d->nvols);
  fprintf(stdout,"  Start volume (d->startvol+1) = %d\n",d->startvol+1);
  fprintf(stdout,"  End volume (d->endvol) = %d\n",d->endvol);
  fflush(stdout);
#endif
}

void setseqmode(struct data *d)
{
  char seqcon[6];
  char function[20];
  strcpy(function,"setseqmode"); /* Set function name */

  /* To figure the sequence mode insist we must have the following */
  if ((ptype("seqcon",&d->p) < 0) || (ptype("apptype",&d->p) < 0)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  The 'procpar' must contain 'seqcon' and 'apptype' parameters:\n\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }

  /* Set sequence mode */
  d->seqmode=0;
  /* Standard cases */
  strcpy(seqcon,*sval("seqcon",&d->p));
  if ((seqcon[3] == 'n') && (seqcon[4] == 'n')) { /* standard 2D multislice */
    if ((seqcon[1] == 'c') && (seqcon[2] == 'c')) d->seqmode=IM2DCC;
    if ((seqcon[1] == 'c') && (seqcon[2] == 's')) d->seqmode=IM2DCS;
    if ((seqcon[1] == 's') && (seqcon[2] == 'c')) d->seqmode=IM2DSC;
    if ((seqcon[1] == 's') && (seqcon[2] == 's')) d->seqmode=IM2DSS;
  } else { /* standard 3D */
    if ((seqcon[2] == 'c') && (seqcon[3] == 'c')) d->seqmode=IM3DCC;
    if ((seqcon[2] == 'c') && (seqcon[3] == 's')) d->seqmode=IM3DCS;
    if ((seqcon[2] == 's') && (seqcon[3] == 'c')) d->seqmode=IM3DSC;
    if ((seqcon[2] == 's') && (seqcon[3] == 's')) d->seqmode=IM3DSS;
  }
  /* Special cases */
  if (spar(d,"apptype","im2Dfse")) {
    if (seqcon[2] == 'c') d->seqmode=IM2DCFSE;
    if (seqcon[2] == 's') d->seqmode=IM2DSFSE;
  }
  if (spar(d,"apptype","im2Depi")) d->seqmode=IM2DEPI;
  if (spar(d,"apptype","im3Dfse")) {
    if (seqcon[3] == 'c') d->seqmode=IM3DCFSE;
    if (seqcon[3] == 's') d->seqmode=IM3DSFSE;
  }

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Sequence mode set to ");
  switch (d->seqmode) {
    case IM2DCC:   fprintf(stdout,"2D with seqcon='*ccnn'\n"); break;
    case IM2DCS:   fprintf(stdout,"2D with seqcon='*csnn'\n"); break;
    case IM2DSC:   fprintf(stdout,"2D with seqcon='*scnn'\n"); break;
    case IM2DSS:   fprintf(stdout,"2D with seqcon='*ssnn'\n"); break;
    case IM2DCFSE: fprintf(stdout,"2D with appmode='im2Dfse' and seqcon='nccnn'\n"); break;
    case IM2DSFSE: fprintf(stdout,"2D with appmode='im2Dfse' and seqcon='ncsnn'\n"); break;
    case IM2DEPI:  fprintf(stdout,"2D with appmode='im2Depi'\n"); break;
    case IM3DCC:   fprintf(stdout,"3D with seqcon='**ccn'\n"); break;
    case IM3DCS:   fprintf(stdout,"3D with seqcon='**csn'\n"); break;
    case IM3DSC:   fprintf(stdout,"3D with seqcon='**scn'\n"); break;
    case IM3DSS:   fprintf(stdout,"3D with seqcon='**ssn'\n"); break;
    case IM3DCFSE: fprintf(stdout,"3D with appmode='im3Dfse' and seqcon='ncccn'\n"); break;
    case IM3DSFSE: fprintf(stdout,"3D with appmode='im3Dfse' and seqcon='nccsn'\n"); break;
  }
  fflush(stdout);
#endif
}

void getblock(struct data *d,int volindex,int DCCflag)
{
  int dim1,dim2,dim3,nr;
  int i,j;
  int datatype;
  char function[20];
  strcpy(function,"getblock"); /* Set function name */

  /* Check d->nvols has been set */
  if (d->nvols == 0) return;

  /* Set datatype */
  datatype=0;
  if (d->fh.status & S_FLOAT) /* 32-bit float */
    datatype=FLT32;
  else if (d->fh.status & S_32) /* 32-bit int */
    datatype=INT32;
  else /* 16-bit int */
    datatype=INT16;

  /* If a previous block has been zero filled matrix size will be incorrect */
  /* Refresh from procpar values */
  if (d->zerofill) setdatapars(d);

  /* Set data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Allocate memory according to nr */
  if ((d->data = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();

  /* Allocate memory according to block size */
  for (i=0;i<nr;i++) { /* loop over receivers */
    if ((d->data[i] = (fftw_complex **)fftw_malloc(dim3*sizeof(fftw_complex *))) == NULL) nomem();
    for (j=0;j<dim3;j++) /* loop over dim3 (slices or 2nd phase encodes) */
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();
  }

  /* Read according to d->type */
  switch(datatype) {
    case FLT32: /* 32-bit float */
      readfvol(d,volindex,DCCflag);
      break;
    case INT32: /* 32-bit int */
      readlvol(d,volindex,DCCflag);
      break;
    case INT16: /* 16-bit int */
      readsvol(d,volindex,DCCflag);
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
/*  zeromax(d);*/
  /* Zero noise measurement storage */
/*  zeronoise(d);*/

}

int readfvol(struct data *d,int volindex,int DCCflag)
{
  float *fdata;        /* Pointer for 32-bit floating point data */
  int dim1,dim2,dim3,nr;
  int startpos,dim3index,scale;
  int i,j,k,l;

#ifdef DEBUG
  char function[20];
  strcpy(function,"readfvol"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex+1);
  fprintf(stdout,"  Reading 32-bit floating point data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Start position */
  startpos=d->startpos;

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
  for (i=0;i<nr;i++) { /* loop over receivers */
    for (j=0;j<dim3;j++) { /* loop over dim3 (slices or 2nd phase encodes) */
      dim3index=startpos+j;
      for (k=0;k<dim2;k++) { /* loop over phase */
        /* Set data pointer */
        if (setoffset(d,volindex,i,dim3index,k)) {
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
    } /* end dim3 loop (slice or 2nd phase encode) */
  } /* end receiver loop */

  /* Free memory */
  free(fdata);
  return(0);
}

int readlvol(struct data *d,int volindex,int DCCflag)
{
  long int *ldata;     /* Pointer for 32-bit integer data */
  int dim1,dim2,dim3,nr;
  int startpos,dim3index,scale;
  int i,j,k,l;

#ifdef DEBUG
  char function[20];
  strcpy(function,"readlvol"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 32-bit integer data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Start position */
  startpos=d->startpos;

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
  for (i=0;i<nr;i++) { /* loop over receivers */
    for (j=0;j<dim3;j++) { /* loop over dim3 (slices or 2nd phase encodes) */
      dim3index=startpos+j;
      for (k=0;k<dim2;k++) { /* loop over phase */
        /* Set data pointer */
        if (setoffset(d,volindex,i,dim3index,k)) {
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
    } /* end dim3 loop (slice or 2nd phase encode) */
  } /* end receiver loop */

  /* Free memory */
  free(ldata);
  return(0);
}

int readsvol(struct data *d,int volindex,int DCCflag)
{
  short int *sdata;    /* Pointer for 16-bit integer data */
  int dim1,dim3,nr;
  int startpos,dim3index,scale;
  int i,j,k,l;

#ifdef DEBUG
  char function[20];
  strcpy(function,"readsvol"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Volume %d:\n",volindex);
  fprintf(stdout,"  Reading 16-bit integer data\n");
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2;
  dim3=d->endpos-d->startpos;
  nr=d->nr;

  /* Start slice */
  startpos=d->startpos;

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
  for (i=0;i<nr;i++) { /* loop over receivers */
    for (j=0;j<dim3;j++) { /* loop over dim3 (slices or 2nd phase encodes) */
      dim3index=startpos+j;
      for (k=0;k<dim1;k++) { /* loop over phase */
        /* Set data pointer */
        if (setoffset(d,volindex,i,dim3index,k)) {
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
    } /* end dim3 loop (slice or 2nd phase encode) */
  } /* end receiver loop */

  /* Free memory */
  free(sdata);
  return(0);
}

int setoffset(struct data *d,int volindex,int receiver,int dim3index,int dim2index)
{
  int dim2,ns,nv2,nr;
  int index,phase,slice=0,phase2=0;
  int arraydim,ne,etl,psscycle;
  long offset=0;
  int blockindex=0;

#ifdef DEBUG
  char function[20];
  strcpy(function,"setoffset"); /* Set function name */
#endif

  /* Data dimensions */
  dim2=d->nv;
  ns=d->ns;
  nv2=d->nv2;
  nr=d->nr;

  /* Number of echoes */
  ne=(int)*val("ne",&d->p);
  if (ne < 1) ne=1; /* Set ne to 1 if 'ne' does not exist */

  /* Allow for compressed multi-echo loop */
  index=volindex/ne;

  /* For 3D scans multiple slices always give different volumes */
  if (d->seqmode > IM3D) {
    ns=(int)*val("ns",&d->p); /* We can use the value of ns */
    index/=ns;                /* Correct the index accordingly */
  }

  /* Set phase according to its actual index */
  phase=phaseindex(d,dim2index);

  /* Set slice/phase2 according to its actual index */
  if (d->seqmode > IM3D) phase2=phase2index(d,dim3index);
  else if (d->seqmode > IM2D) slice=sliceindex(d,dim3index);

  /* Calculate offset according to how the expt is run */
  switch(d->seqmode) {
    case IM2DCC: /* seqcon = *ccnn */
      blockindex=index*nr+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(phase*ns*ne+slice*ne+volindex%ne);
      break;
    case IM2DCS: /* seqcon = *csnn */
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
    case IM2DSC: /* seqcon = *scnn */
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
    case IM2DSS: /* seqcon = *ssnn */
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
    case IM2DCFSE: /* im2Dfse with seqcon = nccnn */
      arraydim=(int)*val("arraydim",&d->p);
      etl=(int)*val("etl",&d->p);
      blockindex=index*nr+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)((phase/etl)*ns*etl+slice*etl+phase%etl);
      break;
    case IM2DSFSE: /* im2Dfse with seqcon = ncsnn */
      arraydim=(int)*val("arraydim",&d->p);
      etl=(int)*val("etl",&d->p);
      blockindex=index*nr+(phase/etl)*(arraydim*etl)/dim2+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(slice*etl+phase%etl);
      break;
    case IM2DEPI: /* im2Depi */
      blockindex=index*nr+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(ns*phase+slice);
      break;
    case IM3DCC: /* seqcon = **cc* */
      blockindex=index*nr+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(phase2*dim2*ns*ne+phase*ns*ne+volindex%(ne*ns));
      break;
    case IM3DCS: /* seqcon = **cs* */
      arraydim=(int)*val("arraydim",&d->p);
      blockindex=index*nr+phase2*arraydim/nv2+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(phase*ns*ne+volindex%(ne*ns));
      break;
    case IM3DSC: /* seqcon = **sc* */
      arraydim=(int)*val("arraydim",&d->p);
      blockindex=index*nr+phase*arraydim/dim2+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(phase2*ns*ne+volindex%(ne*ns));
      break;
    case IM3DSS: /* seqcon = **ss* */
      arraydim=(int)*val("arraydim",&d->p);
      blockindex=index*nr+phase2*arraydim/nv2+phase*arraydim/(nv2*dim2)+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(volindex%(ne*ns));
      break;
    case IM3DCFSE: /* im3Dfse with seqcon = ncccn */
      arraydim=(int)*val("arraydim",&d->p);
      etl=(int)*val("etl",&d->p);
      blockindex=index*nr+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)(phase2*dim2*ns*ne+(phase/etl)*ns*etl+(volindex%ns)*etl+phase%etl);
      break;
    case IM3DSFSE: /* im3Dfse with seqcon = nccsn */
      arraydim=(int)*val("arraydim",&d->p);
      etl=(int)*val("etl",&d->p);
      blockindex=index*nr+phase2*arraydim/nv2+receiver;
      /* Calculate offset to end of the required block */
      offset=sizeof(d->fh)+(long)(blockindex+1)*d->fh.bbytes;
      /* Adjust for required trace */
      offset-=d->fh.ntraces*d->fh.tbytes;
      offset+=d->fh.tbytes*(long)((phase/etl)*ns*etl+(volindex%ns)*etl+phase%etl);
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

void setblock(struct data *d,int dim)
{
  int blockdim;

#ifdef DEBUG
  char function[20];
  strcpy(function,"setblock"); /* Set function name */
#endif

  /* Figure the dim required per block */
  blockdim=dim/d->nblocks;
  if (dim%d->nblocks) blockdim++;

  /* Correct the number of blocks */
  d->nblocks=dim/blockdim;
  if (dim%blockdim) d->nblocks++;

  /* Set start position and end position */
  d->startpos=d->block*blockdim;
  d->endpos=(d->block+1)*blockdim;
  if (d->endpos>dim) d->endpos=dim;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Block %d (of %d)\n",d->block+1,d->nblocks);
  fprintf(stdout,"  Dim = %d\n",dim);
  fprintf(stdout,"  Block dim = %d (from %d to %d)\n",d->endpos-d->startpos,d->startpos,d->endpos-1);
  fflush(stdout);
#endif
}

void getmax(struct data *d)
{
  int dim1,dim2,dim3,nr;
  double re,im,M;
  int i,j,k,l;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"getmax"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Make sure Mval holds M^2 */
  d->max.Mval*=d->max.Mval;
  /* Get maximum and coordinates */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for(k=0;k<dim2;k++) {
        for (l=0;l<dim1;l++) {
          re=fabs(d->data[i][j][k*dim1+l][0]);
          im=fabs(d->data[i][j][k*dim1+l][1]);
          M=re*re+im*im;
          if (M > d->max.Mval) {
            d->max.Mval=M;
            d->max.np=l;
            d->max.nv=k;
            d->max.nv2=d->startpos+j;
          }
          if (re > d->max.Rval) d->max.Rval=re;
          if (im > d->max.Ival) d->max.Ival=im;
        }
      }
    }
  }
  d->max.Mval=sqrt(d->max.Mval);

  /* Set d->max.Rval = d->max.Ival */
  if (d->max.Rval>d->max.Ival) d->max.Ival=d->max.Rval;
  else d->max.Rval=d->max.Ival;

  /* Set data flag */
  d->max.data=TRUE;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs:\n",t2-t1);
  fprintf(stdout,"  d->max.Mval = %f\n",d->max.Mval);
  fprintf(stdout,"  d->max.Rval = %f\n",d->max.Rval);
  fprintf(stdout,"  d->max.Ival = %f\n",d->max.Ival);
  fprintf(stdout,"  d->max.np = %d\n",d->max.np);
  fprintf(stdout,"  d->max.nv = %d\n",d->max.nv);
  fprintf(stdout,"  d->max.nv2 = %d\n",d->max.nv2);
  fflush(stdout);
#endif

}

void getnoise(struct data *d,int mode)
{
  int dim1,dim2,dim3,nr;
  int start1,end1,start2,end2;
  double noisefrac;
  int i,j,k,l;
  double re,im,M2;
  int factor=0;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"getnoise"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  noisefrac=0.0;
  switch(mode) {
    case MK: /* Mask parameters */
      switch(d->datamode) {
        case FID:
          noisefrac=*val("masknoisefrac",&d->p);
          if (!noisefrac) noisefrac=Knoisefraction;
          break;
        default:
          noisefrac=*val("masklvlnoisefrac",&d->p);
          if (!noisefrac) noisefrac=IMnoisefraction;
      } /* end datamode switch */
      break;
    case SM: /* Sensitivity map parameters */
      switch(d->datamode) {
        case FID:
          noisefrac=*val("smapnoisefrac",&d->p);
          if (!noisefrac) noisefrac=Knoisefraction;
          break;
        default:
          noisefrac=*val("masklvlnoisefrac",&d->p);
          if (!noisefrac) noisefrac=IMnoisefraction;
      } /* end datamode switch */
      break;
    default: /* Default parameters */
      switch(d->datamode) {
        case FID:
          noisefrac=Knoisefraction;
          break;
        default:
          noisefrac=IMnoisefraction;
      } /* end datamode switch */
  } /* end mode switch */

  /* For 3D data block check to see if we should sample noise */
  if (d->seqmode >= IM3D) { /* 3D data */
    /* Default start1, end1, start2, end2 for a single processing block */
    start1=0; end1=d->nv2*noisefrac/2+1;
    start2=d->nv2-d->nv2*noisefrac/2; end2=d->nv2;
    /* See if we should skip the processing block */
    if ((d->startpos>=end1) && (d->endpos<start2)) {
#ifdef DEBUG
  fprintf(stdout,"  Skipping data block %d\n",d->block+1);
  fflush(stdout);
#endif
      return;
    }
    /* If not, figure the start1, end1, start2, end2 for the processing block */
    if (d->startpos<end1) {
      if (d->endpos<end1) {
        end1=d->endpos-d->startpos;
        end2=0;
      } else {
        end1-=d->startpos;
      }
    } else {
      end1=0;
    }
    if (d->endpos>=start2) {
      start2-=d->startpos;
      if (start2<0) start2=0;
      end2=d->endpos-d->startpos;
    } else {
      end2=0;
    }
    if (!end2) start2=0;
#ifdef DEBUG
  fprintf(stdout,"  Sampling data block %d (start1 = %d, end1 = %d, start2 = %d, end2 = %d)\n",d->block+1,start1,end1,start2,end2);
  fflush(stdout);
#endif
  } else { /* 2D data */
    start1=0; end1=dim3;
    start2=0; end2=0;    /* Just skip the loop from start2 to end2 */
  }

  factor=d->noise.samples/nr;
  for (i=0;i<nr;i++) {
    d->noise.M[i]  *=factor;
    d->noise.M2[i] *=factor;
    d->noise.Re[i]  *=factor;
    d->noise.Im[i]  *=factor;
    /* For not shifted FID and shifted IMAGE the noise is at edges */
    for (j=start1;j<end1;j++) {
      for(k=0;k<dim2*noisefrac/2;k++) {
        for (l=0;l<dim1*noisefrac/2;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
        for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
      }
      for(k=dim2-dim2*noisefrac/2;k<dim2;k++) {
        for (l=0;l<dim1*noisefrac/2;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
        for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
      }
    }
    for (j=start2;j<end2;j++) {
      for(k=0;k<dim2*noisefrac/2;k++) {
        for (l=0;l<dim1*noisefrac/2;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
        for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
      }
      for(k=dim2-dim2*noisefrac/2;k<dim2;k++) {
        for (l=0;l<dim1*noisefrac/2;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
        for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
          re=d->data[i][j][k*dim1+l][0];
          im=d->data[i][j][k*dim1+l][1];
          M2=re*re+im*im;
          d->noise.M[i]+=sqrt(M2);
          d->noise.M2[i]+=M2;
          d->noise.Re[i]+=fabs(re);
          d->noise.Im[i]+=fabs(im);
          d->noise.samples++;
        }
      }
    }
  }
  /* calculate the mean */
  factor=d->noise.samples/nr;
  for (i=0;i<nr;i++) {
    d->noise.M[i] /=factor;
    d->noise.M2[i] /=factor;
    /* For Real and Imaginary we must consider console type.
       The DDR in VNMRS produces equal noise Re and Im channels - no quad images */
    if (spar(d,"console","vnmrs")) { /* VNMRS */
      d->noise.Re[i] += d->noise.Im[i];
      d->noise.Re[i] /=2.0;
      d->noise.Im[i] = d->noise.Re[i];
    }
    d->noise.Re[i] /=factor;
    d->noise.Im[i] /=factor;
  }

  /* Now average over all receivers */
  d->noise.avM =0;
  d->noise.avM2 =0;
  d->noise.avRe =0;
  d->noise.avIm =0;
  for (i=0;i<nr;i++) {
    d->noise.avM += d->noise.M[i];
    d->noise.avM2 += d->noise.M2[i];
    d->noise.avRe += d->noise.Re[i];
    d->noise.avIm += d->noise.Im[i];
  }
  d->noise.avM /=nr;
  d->noise.avM2 /=nr;
  d->noise.avRe /=nr;
  d->noise.avIm /=nr;

  /* Set data flag */
  d->noise.data=TRUE;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Noise data averaged over %d points: took %f secs\n",d->noise.samples/nr,t2-t1);
  for (i=0;i<nr;i++) {
  fprintf(stdout,"  Receiver %d: M = %.3f, M2 = %.3f, Re = %.3f, Im = %.3f\n",
    i,d->noise.M[i],d->noise.M2[i],d->noise.Re[i],d->noise.Im[i]);
  }
  fprintf(stdout,"  Average:    M = %.3f, M2 = %.3f, Re = %.3f, Im = %.3f\n",
    d->noise.avM,d->noise.avM2,d->noise.avRe,d->noise.avIm);
  fflush(stdout);
#endif

}

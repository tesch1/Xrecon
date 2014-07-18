/* dproc1D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dproc1D.c: 1D Data processing routines                                    */
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

#define SOURCEFILE "common1D/dproc1D.c"

void fft1D(struct data *d,int dataorder)
{
  fftw_complex *data;
  fftw_plan p;
  int dim1=0,dim3=0,nr=0;
  int i,j,k;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"fft1D"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

/*
  fftw_import_wisdom_from_file(FILE *input_file);
*/

  /* Set data dimensions */
  switch (dataorder) {
    case D1: dim1=d->np/2; dim3=d->fh.ntraces; nr=d->nr; break;
    case D3: dim1=d->nv2; dim3=(d->endpos-d->startpos)*d->np/2; nr=d->nr; break;
  }

  /* Allocate memory for a trace */
  if ((data = (fftw_complex *)fftw_malloc(dim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Measured plan for ft */
  p=fftw_plan_dft_1d(dim1,data,data,FFTW_FORWARD,FFTW_MEASURE);

  /* ft  profiles ... */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim1;k++) { data[k][0]=d->data[i][j][k][0]; data[k][1]=d->data[i][j][k][1]; }
      /* Estimated plan */
/*
      p=fftw_plan_dft_1d(dim1,d->data[i][j],d->data[i][j],FFTW_FORWARD,FFTW_ESTIMATE);
*/
      /* ft */
      fftw_execute(p);
      for (k=0;k<dim1;k++) { d->data[i][j][k][0]=data[k][0]; d->data[i][j][k][1]=data[k][1]; }
    }
  }
  /* ... and tidy up */
  fftw_destroy_plan(p);
  fftw_free(data);

  /* Set datamode for image data */
  d->datamode=IMAGE;

  /* Update shifted flag */
  d->shift=FALSE; /* not shifted */

  /* Reset max pars to default values */
/*  zeromax(d); */

  /* Zero noise */
/*  zeronoise(d); */

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  ft (%d): %d trace(s), %d receiver(s): took %f secs\n",
    dim1,dim3,nr,t2-t1);
  fflush(stdout);
#endif

}

int shiftdata1D(struct data *d,int mode,int dataorder)
{
  int dim1=0;
  int shft=0;
  char function[20];
  strcpy(function,"shiftdata1D"); /* Set function name */

  /* Set data dimensions */
  switch (dataorder) {
    case D1: dim1=d->np/2; break;
    case D3: dim1=d->nv2;  break;
  }

  /* Set shift points */
  switch(mode) {
    case STD:
      shft=dim1/2;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  switch(dataorder) {
    case D1: fprintf(stdout,"  Standard shift: npshft = %d\n",shft); break;
    case D3: fprintf(stdout,"  Standard shift: nv2shft = %d\n",shft); break;
  }
  fflush(stdout);
#endif
      break;
    case OPT:
      /* getmax1D sets maximum and its coordinates */
      /* if (!d->max.data) getmax1D(d); */
      switch(dataorder) {
        case D1: shft=d->max.np;  break;
        case D3: shft=d->max.nv2; break;
      }
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Optimized shift: shft = %d\n",shft);
  fflush(stdout);
#endif
      break;
    default:
      /* Invalid mode */
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Invalid 2nd argument %s(*,'type',*)\n",function);
      return(1);
      break;
  } /* end mode switch */

  /* Shift the data */
  shift1Ddata(d,shft,dataorder);

  /* Update shift flag */
  if (d->shift) d->shift=FALSE;  /* if 'shifted' set to 'not shifted' */
  else d->shift=TRUE;            /* else set flag to 'shifted' */

  return(0);
}

int shift1Ddata(struct data *d,int shft,int dataorder)
{
  fftw_complex *data=0;
  int dim1=0,dim3=0,nr=0;
  int i,j,k;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"shift1Ddata"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  shft = %d: ",shft);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Set data dimensions */
  switch (dataorder) {
    case D1: dim1=d->np/2; dim3=d->fh.ntraces; nr=d->nr; break;
    case D3: dim1=d->nv2; dim3=(d->endpos-d->startpos)*d->np/2; nr=d->nr; break;
  }

  /* Allocate memory for data */
  if ((data = (fftw_complex *)fftw_malloc(dim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Loop over receivers and slices */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      /* Fill data with shifted data */
      for (k=0;k<shft;k++) {
        data[k+dim1-shft][0]=d->data[i][j][k][0];
        data[k+dim1-shft][1]=d->data[i][j][k][1];
      }
      for (k=shft;k<dim1;k++) {
        data[k-shft][0]=d->data[i][j][k][0];
        data[k-shft][1]=d->data[i][j][k][1];
      }
      /* Copy shifted data back to d->data */
      for (k=0;k<dim1;k++) {
        d->data[i][j][k][0]=data[k][0];
        d->data[i][j][k][1]=data[k][1];
      }
    }
  }

  fftw_free(data);

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"took %f secs\n",t2-t1);
  fflush(stdout);
#endif

  return(0);
}

int shifttrace1D(struct data *d,int trace,int mode)
{
  int dim1;
  int npshft=0;
  char function[20];
  strcpy(function,"shifttrace1D"); /* Set function name */

  /* Data dimensions */
  dim1=d->np/2;

  /* Set shift points */
  switch(mode) {
    case STD:
      npshft=dim1/2;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Standard shift for trace %d: npshft = %d\n",trace,npshft);
  fflush(stdout);
#endif
      break;
    case OPT:
      /* getmax2D sets maximum and its coordinates */
      if (!d->max.data) getmaxtrace1D(d,trace);
      npshft=d->max.np;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Optimized shift for trace %d: npshft = %d\n",trace,npshft);
  fflush(stdout);
#endif
      break;
    default:
      /* Invalid mode */
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Invalid 3rd argument %s(*,*,'type')\n",function);
      return(1);
      break;
  } /* end mode switch */

  /* Shift the data */
  shift1Dtrace(d,trace,npshft);

  return(0);
}

int shift1Dtrace(struct data *d,int trace,int npshft)
{
  fftw_complex *data;
  int dim1,nr;
  int i,j;

  /* Data dimensions */
  dim1=d->np/2;
  nr=d->nr;

  /* Allocate memory for slice data */
  if ((data = (fftw_complex *)fftw_malloc(dim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Loop over receivers */
  for (i=0;i<nr;i++) {
    /* Fill trace with shifted data */
    for (j=0;j<npshft;j++) {
      data[j+dim1-npshft][0]=d->data[i][trace][j][0];
      data[j+dim1-npshft][1]=d->data[i][trace][j][1];
    }
    for (j=npshft;j<dim1;j++) {
      data[j-npshft][0]=d->data[i][trace][j][0];
      data[j-npshft][1]=d->data[i][trace][j][1];
    }
    /* Copy shifted data back to d->data */
    for (j=0;j<dim1;j++) {
      d->data[i][trace][j][0]=data[j][0];
      d->data[i][trace][j][1]=data[j][1];
    }
  }

  fftw_free(data);

  return(0);
}

void weightdata1D(struct data *d,int mode,int dataorder)
{
  int i,j,k;
  double f;
  int dim1=0,dim3=0,nr=0;
  double lb=0,gf=0,sb=0,at;
  double *weight;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"weightdata1D"); /* Set function name */
#endif

  switch(mode) {
    case STD: /* Standard parameters */
      switch (dataorder) {
        case D1: lb=*val("lb",&d->p); gf=*val("gf",&d->p); sb=*val("sb",&d->p); break;
        case D3: lb=*val("lb2",&d->p); gf=*val("gf2",&d->p); sb=*val("sb2",&d->p); break;
      }
      break;
    default:
      break;
  } /* end mode switch */

  /* Return if no weighting is active */
  if (!lb && !gf && !sb) return;

  /* Return if not FID data and not shifted */
/*  if ((d->datamode != FID) || !d->shift) return; */

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  switch (dataorder) {
    case D1: dim1=d->np/2; dim3=d->fh.ntraces; nr=d->nr; break;
    case D3: dim1=d->nv2; dim3=(d->endpos-d->startpos)*d->np/2; nr=d->nr; break;
  }

  at=*val("at",&d->p);

  /* Calculate weighting */
  if ((weight = (double *)malloc(dim1*sizeof(double))) == NULL) nomem();
  for (k=0;k<dim1;k++) weight[k]=1.0;

  /* Lorentzian line broadening */
  if (lb) { /* lb is active */
    for (k=0;k<dim1/2;k++) {
      /* Lorentzian broadening as fraction of FOV */
/* f=exp(-k*M_PI*lb); */
      /* Lorentzian broadening in # pixels */
      /* f=exp(-k*M_PI*lb/dim1); */
      /* Lorentzian broadening in Hz (as VnmrJ) */
f=exp(-k*at*M_PI*lb/dim1);
      weight[k] *=f;
      weight[dim1-1-k] *=f;
    }
  }

  /* Gaussian line broadening */
  if (gf) { /* gf is active */
    for (k=0;k<dim1/2;k++) {
      /* Gaussian broadening as fraction of FOV */
      f=k*M_PI*gf;
      f=exp(-f*f);
      weight[k] *=f;
      weight[dim1-1-k] *=f;
    }
  }

  /* Sinebell line broadening */
  if (sb) { /* sb is active */
    for (k=0;k<dim1/2;k++) {
      /* Sinebell broadening as fraction of FOV */
      if (k*M_PI*sb < M_PI/2)
        f=cos(k*M_PI*sb);
      else 
        f=0.0;
      weight[k] *=f;
      weight[dim1-1-k] *=f;
    }
  }

  /* Weight the data */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim1;k++) {
        d->data[i][j][k][0] *=weight[k];
        d->data[i][j][k][1] *=weight[k];
      }
    }
  }

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  switch (dataorder) {
    case D1:
      if (lb) fprintf(stdout,"  lb weighting = %f\n",lb);
      if (gf) fprintf(stdout,"  gf weighting = %f\n",gf);
      if (sb) fprintf(stdout,"  sb weighting = %f\n",sb);
      break;
    case D3:
      if (lb) fprintf(stdout,"  lb2 weighting = %f\n",lb);
      if (gf) fprintf(stdout,"  gf2 weighting = %f\n",gf);
      if (sb) fprintf(stdout,"  sb2 weighting = %f\n",sb);
      break;
  }
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fflush(stdout);
#endif
}

void zerofill1D(struct data *d,int mode,int dataorder)
{
  int i,j,k;
  int ix1;
  int dim1=0,dim3=0,nr=0;
  int fn=0;
  int ndim1;
  fftw_complex *data;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"zerofill1D"); /* Set function name */
#endif

  switch(mode) {
    case MK: /* Mask parameters */
      fn=*val("maskfn",&d->p);
      break;
    default: /* Internally set */
      switch (dataorder) {
        case D1: fn=d->fn;  break;
        case D3: fn=d->fn2; break;
      }
      break;
  } /* end mode switch */

  /* Check that fn is active */
  if (!fn) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Set data dimensions */
  switch (dataorder) {
    case D1: dim1=d->np/2; dim3=d->fh.ntraces; nr=d->nr; break;
    case D3: dim1=d->nv2; dim3=(d->endpos-d->startpos)*d->np/2; nr=d->nr; break;
  }

  /* Make sure fn is exactly divisible by 4 */
  fn /=4; fn *=4;

  /* New data dimension */
  ndim1=fn/2;

  /* Make data the new size and initialise to zero */
  if ((data = (fftw_complex *)fftw_malloc(ndim1*sizeof(fftw_complex))) == NULL) nomem();
  for (k=0;k<ndim1;k++) {
    data[k][0]=0.0;
    data[k][1]=0.0;
  }

  /* Copy from d->data[i][j] to data, zero filling as we go */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      if (ndim1 < dim1) {
        for (k=0;k<ndim1/2;k++) {
          data[k][0]=d->data[i][j][k][0];
          data[k][1]=d->data[i][j][k][1];
        }
        for (k=dim1-ndim1/2;k<dim1;k++) {
          ix1=k-(dim1-ndim1);
          data[ix1][0]=d->data[i][j][k][0];
          data[ix1][1]=d->data[i][j][k][1];
        }
      } else {
        for (k=0;k<dim1/2;k++) {
          data[k][0]=d->data[i][j][k][0];
          data[k][1]=d->data[i][j][k][1];
        }
        for (k=dim1/2;k<dim1;k++) {
          ix1=k-(dim1-ndim1);
          data[ix1][0]=d->data[i][j][k][0];
          data[ix1][1]=d->data[i][j][k][1];
        }
      }

      /* free d->data[i][j], reallocate and copy data back in its place */
      fftw_free(d->data[i][j]);
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(ndim1*sizeof(fftw_complex))) == NULL) nomem();
      for (k=0;k<ndim1;k++) {
        d->data[i][j][k][0]=data[k][0];
        d->data[i][j][k][1]=data[k][1];
      }
    }
  }

  /* Update parameters to reflect new data size */
  switch (dataorder) {
    case D1: d->np=fn;    break;
    case D3: d->nv2=fn/2; break;
  }

  /* Set zerofill flag */
  d->zerofill=TRUE;

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Zero filling to give %d data: took %f secs\n",ndim1,t2-t1);
  fflush(stdout);
#endif

}

void phaseramp1D(struct data *d,int mode)
{
  int dim1,dim3,nr;
  int i,j,k;
  double re,im,M,theta,factor;
  double offset,fov;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"phaseramp1D"); /* Set function name */
#endif

  /* Data dimensions */
  dim1=d->nv2; dim3=(d->endpos-d->startpos)*d->np/2; nr=d->nr;

  switch (mode) {
    case PHASE2:
      /* Return if there is no phase ramp to apply */
      offset=*val("ppe2",&d->p);
      if (offset == 0.0) return;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif
      fov=*val("lpe2",&d->p);
      /* Set phase ramp factor to correct for offset */
      factor=-2*M_PI*offset/fov;
      /* Apply phase ramp to generate frequency shift */
      for (i=0;i<nr;i++) {
        for (j=0;j<dim3;j++) {
          for (k=0;k<dim1/2;k++) {
            re=d->data[i][j][k][0];
            im=d->data[i][j][k][1];
            M=sqrt(re*re+im*im);
            theta = atan2(im,re) + factor*(k);
            d->data[i][j][k][0]=M*cos(theta);
            d->data[i][j][k][1]=M*sin(theta);
          }
          for (k=dim1/2;k<dim1;k++) {
            re=d->data[i][j][k][0];
            im=d->data[i][j][k][1];
            M=sqrt(re*re+im*im);
            theta = atan2(im,re) + factor*(k-dim1);
            d->data[i][j][k][0]=M*cos(theta);
            d->data[i][j][k][1]=M*sin(theta);
          }
        }
      }
#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Phase encode phase ramp (%d traces): took %f secs\n",dim3,t2-t1);
  fflush(stdout);
#endif
      break;
  }

}

void phasedata1D(struct data *d,int mode,int dataorder)
{
  int dim1=0,dim3=0,nr=0;
  double rp=0,lp=0;
  double re,im,M,theta;
  double factor;
  int i,j,k;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"phasedata2D"); /* Set function name */
#endif

  switch(mode) {
    case VJ: /* Standard VnmrJ parameters */
      if ((strcmp(*sval("imRE",&d->p),"y"))
        && (strcmp(*sval("imIM",&d->p),"y")))
        return;
      break;
  } /* end mode switch */

  /* Get phasing parameters */
  switch (dataorder) {
    case D1: rp=*val("rp",&d->p); lp=*val("lp",&d->p); break;
    case D3: rp=0; lp=*val("lp2",&d->p); break;
  }

  /* Return if no phasing is required */
  if (!rp && !lp) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Set data dimensions */
  switch (dataorder) {
    case D1: dim1=d->np/2; dim3=d->fh.ntraces; nr=d->nr; break;
    case D3: dim1=d->nv2; dim3=(d->endpos-d->startpos)*d->np/2; nr=d->nr; break;
  }

  /* Phase */
  /* Set factor so that lp adjusts the phase by the specified
     phase in degrees accross the profile */
  factor=DEG2RAD/dim1;
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim1/2;k++) {
        re=d->data[i][j][k][0];
        im=d->data[i][j][k][1];
        M=sqrt(re*re+im*im);
        theta = atan2(im,re) + rp*DEG2RAD + lp*factor*k;
        d->data[i][j][k][0]=M*cos(theta);
        d->data[i][j][k][1]=M*sin(theta);
      }
      for (k=dim1/2;k<dim1;k++) {
        re=d->data[i][j][k][0];
        im=d->data[i][j][k][1];
        M=sqrt(re*re+im*im);
        theta = atan2(im,re) + rp*DEG2RAD + lp*factor*(k-dim1);
        d->data[i][j][k][0]=M*cos(theta);
        d->data[i][j][k][1]=M*sin(theta);
      }
    }
  }

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs:\n",t2-t1);
  fflush(stdout);
#endif

}

void getmaxtrace1D(struct data *d,int trace)
{
  int dim1,nr;
  double re,im,M;
  int i,j;

#ifdef DEBUG
  struct timeval tp;
  double t1=0,t2=0;
  char function[20];
  strcpy(function,"getmaxtrace1D"); /* Set function name */
  if (trace == 0) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    gettimeofday(&tp, NULL);
    t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    fflush(stdout);
  }
#endif

  /* Data dimensions */
  dim1=d->np/2;
  nr=d->nr;

  /* Set some defaults */
  zeromax(d);

  /* Now get maximum and coordinates */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim1;j++) {
      re=fabs(d->data[i][trace][j][0]);
      im=fabs(d->data[i][trace][j][1]);
      M=re*re+im*im;
      if (M > d->max.Mval) {
        d->max.Mval=M;
        d->max.np=j;
      }
      if (re > d->max.Rval) d->max.Rval=re;
      if (im > d->max.Ival) d->max.Ival=im;
    }
  }
  d->max.Mval=sqrt(d->max.Mval);

  /* Set d->max.Rval = d->max.Ival */
  if (d->max.Rval>d->max.Ival) d->max.Ival=d->max.Rval;
  else d->max.Rval=d->max.Ival;

  /* Set data flag */
  d->max.data=TRUE;

#ifdef DEBUG
  if (trace == 0) {
    gettimeofday(&tp, NULL);
    t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
    fprintf(stdout,"  Took %f secs:\n",t2-t1);
    fprintf(stdout,"  d->max.Mval = %f\n",d->max.Mval);
    fprintf(stdout,"  d->max.Rval = %f\n",d->max.Rval);
    fprintf(stdout,"  d->max.Ival = %f\n",d->max.Ival);
    fprintf(stdout,"  d->max.np = %d\n",d->max.np);
    fprintf(stdout,"  d->max.nv = %d\n",d->max.nv);
    fprintf(stdout,"  d->max.nv2 = %d\n",d->max.nv2);
    fflush(stdout);
  }
#endif

}

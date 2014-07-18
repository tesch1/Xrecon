/* dprocEPI.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dprocEPI.c: EPI data processing routines                                  */
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

#define SOURCEFILE "reconEPI/dprocEPI.c"

void addscaledEPIref(struct data *d,struct data *ref1,struct data *ref2)
{
  int dim1,dim2,dim3,nr;
  int nseg;
  int i,j,k,l,ix;
  double re,re2,im,im2,M;
  double level,M2,avM2;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"addscaledEPIref"); /* Set function name */
#endif

  if (ref1->datamode == NONE) return;
  if (ref2->datamode == NONE) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  nseg=(int)*val("nseg",&d->p);

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  getnoise2D(ref2,STD);
  level=*val("nlevel",&ref2->p);
  avM2=level*ref2->noise.avM2;

  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim2;k++) {
        ix = k*dim1;
        for (l=0;l<dim1;l++) {
          re=d->data[i][j][ix+l][0];
          im=d->data[i][j][ix+l][1];
          re2=ref2->data[i][j][ix+l][0];
          im2=ref2->data[i][j][ix+l][1];
          M2=re2*re2+im2*im2;
          /* Average with the scaled reference if pixel is above noise level */
          if (M2>avM2) {
            M=sqrt(re*re+im*im)/sqrt(M2);
            d->data[i][j][ix+l][0] += M*ref1->data[i][j][ix+l][0];
            d->data[i][j][ix+l][1] += M*ref1->data[i][j][ix+l][1];
            d->data[i][j][ix+l][0] /=2;
            d->data[i][j][ix+l][1] /=2;
          }
        }
      }
    }
  }

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Combining with reference data: took %f secs\n",t2-t1);
  fflush(stdout);
#endif

}

void addEPIref(struct data *d,struct data *ref)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l,ix;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"addEPIref"); /* Set function name */
#endif

  if (ref->datamode == NONE) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim2;k++) {
        ix = k*dim1;
        for (l=0;l<dim1;l++) {
          d->data[i][j][ix+l][0] += ref->data[i][j][ix+l][0];
          d->data[i][j][ix+l][1] += ref->data[i][j][ix+l][1];
          d->data[i][j][ix+l][0] /=2;
          d->data[i][j][ix+l][1] /=2;
        }
      }
    }
  }

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Combining with reference data: took %f secs\n",t2-t1);
  fflush(stdout);
#endif

}

void ftnvEPI(struct data *d)
{
  fftw_complex *data;
  fftw_plan p;
  int dim1,dim2,dim3,nr;
  int nphase,fn,nvshft;
  int i,j,k,l,ix;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"ftnvEPI"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Starting 1D FT in phase ...\n");
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Fill data to nphase */
  nvfillEPI(d);

  /* Shift for FT */
  nphase=(int)*val("nphase",&d->p);
  nvshft=nphase/2-1;

  shift2Ddata(d,0,nvshft);

  /* Phase ramp the data to correct for phase encode offset ppe */
  phaserampEPI(d,PHASE);

  d->datamode=FID;               /* Flag that data is FID */
  weightdata2D(d,PHASE);         /* Weight data in read using standard VnmrJ parameters */
  d->datamode=IMAGE;             /* Flag that data is IMAGE */

  /* Zero fill in phase */
  fn=d->fn;
  d->fn = 0;          /* Turn off zero fill in read */
  zerofill2D(d,STD);  /* zero fill */
  d->fn=fn;           /* Reset zero fill in read */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Allocate memory for a phase encode trace */
  if ((data = (fftw_complex *)fftw_malloc(dim2*sizeof(fftw_complex))) == NULL) nomem();

  /* Measured plan for ft */
  p=fftw_plan_dft_1d(dim2,data,data,FFTW_FORWARD,FFTW_MEASURE);

  /* ft  profiles ... */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (l=0;l<dim1;l++) {
        for (k=0;k<dim2;k++) {
          ix = k*dim1;
          data[k][0]=d->data[i][j][ix+l][0];
          data[k][1]=d->data[i][j][ix+l][1];
        }
        fftw_execute(p); /* ft */
        for (k=0;k<dim2;k++) {
          ix = k*dim1;
          d->data[i][j][ix+l][0]=data[k][0];
          d->data[i][j][ix+l][1]=data[k][1];
        }
      }
    }
  }
  /* ... and tidy up */
  fftw_destroy_plan(p);
  fftw_free(data);

  /* Shift to get profiles */
  nvshft=dim2/2;
  shift2Ddata(d,0,nvshft);

  /* Flag as IMAGE data */
  d->datamode=IMAGE;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  1D FT of %d x %d data: took %f secs\n",dim1,dim2,t2-t1);
  fflush(stdout);
#endif

}

void nvfillEPI(struct data *d)
{
  fftw_complex *data;
  int dim1,dim2,dim3,nr;
  int nphase,nfill;
  int i,j,k,l,ix1,ix2;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"nvfillEPI"); /* Set function name */
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Return if nothing to do */
  nphase=(int)*val("nphase",&d->p);
  if (dim2 >= nphase) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  if ((data = (fftw_complex *)fftw_malloc(nphase*dim1*sizeof(fftw_complex))) == NULL) nomem();

  nfill=nphase-dim2;
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      ix1=0;
      for (k=0;k<nfill;k++) {
        for (l=0;l<dim1;l++) {
          data[ix1][0]=0.0;
          data[ix1][1]=0.0;
          ix1++;
        }
      }
      ix2=0;
      for (k=0;k<dim2;k++) {
        for (l=0;l<dim1;l++) {
          data[ix1][0]=d->data[i][j][ix2][0];
          data[ix1][1]=d->data[i][j][ix2][1];
          ix1++; ix2++;
        }
      }
      fftw_free(d->data[i][j]);
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(nphase*dim1*sizeof(fftw_complex))) == NULL) nomem();
      ix1=0;
      for (k=0;k<nphase;k++) {
        for (l=0;l<dim1;l++) {
          d->data[i][j][ix1][0]=data[ix1][0];
          d->data[i][j][ix1][1]=data[ix1][1];
          ix1++;
        }
      }
    }
  }

  fftw_free(data);

  /* Set d->nv for the data set */
  d->nv=nphase;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Phase zerofilled for kzero: took %f secs\n",t2-t1);
  fflush(stdout);
#endif

}

void ftnpEPI(struct data *d)
{
  fftw_complex *data;
  fftw_plan p;
  int dim1,dim2,dim3,nr;
  int image,fn1,npshft;
  int i,j,k,l,ix;
  double oversample;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"ftnpEPI"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Starting 1D FT in read ...\n");
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Shift for FT */
  npshft=d->np/4;
  shift2Ddata(d,npshft,0);
  d->shift=TRUE; /* Flag that data is shifted */

  /* Phase ramp the data to correct for readout offset pro */
  phaserampEPI(d,READ);

  /* Assume parameter "tep" has been set properly to eliminate zipper */
  /* Then we can weight data in read as per usual */
  /* Process data as directed by image parameter */
  image=(int)getelem(d,"image",d->vol);
  switch (image) {
    case 0:   /* Reference, no phase-encode */
      weightdata2D(d,REFREAD);          /* Weight data in read using reference VnmrJ parameters */
      break;
    case -2:  /* Inverted Readout Reference, no phase-encode */
      weightdata2D(d,REFREAD);          /* Weight data in read using reference VnmrJ parameters */
      break;
    default:  /* Phase encoded data */
      weightdata2D(d,READ);             /* Weight data in read using standard VnmrJ parameters */
      break;
  }

  /* Zero fill in read */
  oversample=*val("oversample",&d->p);
  d->fn *= oversample;           /* Correct d->fn for oversample */
  fn1=d->fn1;                    /* Take a copy of d->fn1 */
  d->fn1 = 0;                    /* Disable zero filling in phase dimension */
  zerofill2D(d,STD);             /* zero fill */
  d->fn1=fn1;                    /* Reset d->fn1 */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Allocate memory for a trace */
  if ((data = (fftw_complex *)fftw_malloc(dim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Measured plan for ft */
  p=fftw_plan_dft_1d(dim1,data,data,FFTW_FORWARD,FFTW_MEASURE);

  /* ft  profiles ... */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim2;k++) {
        ix = k*dim1;
        for (l=0;l<dim1;l++) {
          data[l][0]=d->data[i][j][ix+l][0]; data[l][1]=d->data[i][j][ix+l][1];
        }
        fftw_execute(p); /* ft */
        for (l=0;l<dim1;l++) {
          d->data[i][j][ix+l][0]=data[l][0]; d->data[i][j][ix+l][1]=data[l][1];
        }
      }
    }
  }
  /* ... and tidy up */
  fftw_destroy_plan(p);
  fftw_free(data);

  /* Shift to get profiles */
  npshft=d->np/4;
  shift2Ddata(d,npshft,0);

  /* Flag as IMAGE data */
  d->datamode=IMAGE;

  /* If oversampled, zoom to get the requested FOV */
  if (oversample>1) zoomEPI(d);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  1D FT of %d x %d data: took %f secs\n",dim1,dim2,t2-t1);
  fflush(stdout);
#endif

}

void zoomEPI(struct data *d)
{
  int dim1,dim2,nread,fn;

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv;

  /* Requested data dimensions */
  nread=(int)*val("nread",&d->p)/2;
  fn=(int)*val("fn",&d->p)/2;

  /* Zoom in readout dimension */
  if (fn>0) zoomdata2D(d,(dim1-fn)/2,fn,0,dim2);
  else zoomdata2D(d,(dim1-nread)/2,nread,0,dim2);

}

void revreadEPI(struct data *d)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l,ix;
  fftw_complex *data;

#ifdef DEBUG
  char function[20];
  strcpy(function,"revreadEPI"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Reversing data in readout dimension ... ");
  fflush(stdout);
#endif

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Make data the new size */
  if ((data = (fftw_complex *)fftw_malloc(dim1*sizeof(fftw_complex))) == NULL) nomem();

  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim2;k++) {
        ix = k*dim1+dim1; /* not ix = k*dim1+dim1-1 because central pixel is offset ?? */
        for (l=1;l<dim1;l++) {
          data[l][0] = d->data[i][j][ix-l][0];
          data[l][1] = d->data[i][j][ix-l][1];
        }
        ix = k*dim1;
        for (l=0;l<dim1;l++) {
          d->data[i][j][ix+l][0] = data[l][0];
          d->data[i][j][ix+l][1] = data[l][1];
        }
      }
    }
  }

#ifdef DEBUG
  fprintf(stdout,"done\n");
  fflush(stdout);
#endif

}

void prepEPIref(struct data *ref1,struct data *ref2)
{
  int dim1,dim2,dim3,nr;

#ifdef DEBUG
  char function[20];
  strcpy(function,"prepEPIref"); /* Set function name */
#endif

  if (ref1->datamode == NONE) return;

  /* Initial data dimensions */
  dim1=ref1->np/2; dim2=ref1->nv; dim3=ref1->endpos-ref1->startpos; nr=ref1->nr;

  ftnpEPI(ref1);                /* FT along readout dimension */
  phaseEPI(ref1,ref2);          /* Phase correct with the reference */
  navcorrEPI(ref1);             /* Phase correct with the navigator */
  stripEPInav(ref1);            /* Strip the navigators */

  /* Initial data dimensions */
  dim1=ref1->np/2; dim2=ref1->nv; dim3=ref1->endpos-ref1->startpos; nr=ref1->nr;

  ftnvEPI(ref1);                /* FT along phase encode dimension */

}

void phaseEPIref(struct data *ref1,struct data *ref2,struct data *ref3)
{

  if (ref1->datamode != EPIREF) return;
  if (ref2->datamode != EPIREF) return;

  copy2Ddata(ref1,ref3);

  phaseEPI(ref3,ref2);

  navcorrEPI(ref3);

}

void phaseEPI(struct data *d,struct data *ref)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l,ix;
  double re,im,M,refphase,phase;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"phaseEPI"); /* Set function name */
#endif

  if (ref->datamode != EPIREF) return;
  if (spar(d,"refpc","n")) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim2;k++) {
        ix = k*dim1;
        for (l=0;l<dim1;l++) {
          re=ref->data[i][j][ix+l][0];
          im=ref->data[i][j][ix+l][1];
          refphase=atan2(im,re);
          re=d->data[i][j][ix+l][0];
          im=d->data[i][j][ix+l][1];
          phase=atan2(im,re);
          M = sqrt(re*re + im*im);
          d->data[i][j][ix+l][0]=M*cos(phase-refphase);
          d->data[i][j][ix+l][1]=M*sin(phase-refphase);
        }
      }
    }
  }

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Readout phase correction (%d traces): took %f secs\n",dim2,t2-t1);
  fflush(stdout);
#endif

}

void phaserampEPI(struct data *d,int mode)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l,ix;
  double re,im,M,theta,factor;
  double offset,fov,oversample,readstretch;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"phaserampEPI"); /* Set function name */
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  switch (mode) {
    case READ:
      /* Return if there is no phase ramp to apply */
      offset=*val("pro",&d->p);
      if (offset == 0.0) return;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif
      fov=*val("lro",&d->p);
      oversample=*val("oversample",&d->p);
      readstretch=*val("readstretch",&d->p);
      /* Set phase ramp factor to correct for offset */
      factor=2*M_PI*offset/(fov*oversample);
      /* Adjust phase ramp factor to correct for scaled FOV */
      factor*=(dim1+oversample*readstretch)/dim1;
      /* Adjust phase ramp factor to correct for Inverted Readout Reference scans */
      switch ((int)getelem(d,"image",d->vol)) {
        case -1:   /* Inverted Readout Reference, phase-encode */
          factor *=-1;
          break;
        case -2:  /* Inverted Readout Reference, no phase-encode */
          factor *=-1;
          break;
      }
      /* Apply phase ramp to generate frequency shift */
      for (i=0;i<nr;i++) {
        for (j=0;j<dim3;j++) {
          for (k=0;k<dim2;k++) {
            ix = k*dim1;
            for (l=0;l<dim1/2;l++) {
              re=d->data[i][j][ix+l][0];
              im=d->data[i][j][ix+l][1];
              M=sqrt(re*re+im*im);
              theta = atan2(im,re) + factor*(l);
              d->data[i][j][ix+l][0]=M*cos(theta);
              d->data[i][j][ix+l][1]=M*sin(theta);
            }
            for (l=dim1/2;l<dim1;l++) {
              re=d->data[i][j][ix+l][0];
              im=d->data[i][j][ix+l][1];
              M=sqrt(re*re+im*im);
              theta = atan2(im,re) + factor*(l-dim1);
              d->data[i][j][ix+l][0]=M*cos(theta);
              d->data[i][j][ix+l][1]=M*sin(theta);
            }
          }
        }
      }
#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Readout phase ramp (%d traces): took %f secs\n",dim2,t2-t1);
  fflush(stdout);
#endif
      break;
    case PHASE:
      phaseramp2D(d,PHASE);
      break;
  }

}

void navcorrEPI(struct data *d)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l,n,ix,ix2;
  int nnav,etl,nseg,shotnv;
  double re,im,M,phase,*navphase;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"navcorrEPI"); /* Set function name */
#endif

  /* Return if navigator correction not required */
  if (spar(d,"nav_type","off")) return;
  nnav=(int)*val("nnav",&d->p);
  if (nnav == 0) return;

  etl=(int)*val("etl",&d->p);
  nseg=(int)*val("nseg",&d->p);
  shotnv=nnav+etl;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Correct each shot according to last navigator of the shot */
  if ((navphase = (double *)malloc(dim1*sizeof(double))) == NULL) nomem();
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (n=0;n<nseg;n++) {
        ix = (n*(shotnv)+nnav-1)*dim1;
        for (l=0;l<dim1;l++) {
          re=d->data[i][j][ix+l][0];
          im=d->data[i][j][ix+l][1];
          navphase[l]=atan2(im,re);
        }
        ix = n*shotnv*dim1;
        for (k=0;k<shotnv;k++) {
          ix2=ix+k*dim1;
          for (l=0;l<dim1;l++) {
            re=d->data[i][j][ix2+l][0];
            im=d->data[i][j][ix2+l][1];
            phase=atan2(im,re);
            M = sqrt(re*re + im*im);
            d->data[i][j][ix2+l][0]=M*cos(phase-navphase[l]);
            d->data[i][j][ix2+l][1]=M*sin(phase-navphase[l]);
          }
        }
      }
    }
  }
  free(navphase);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Navigator phase correction (%d traces): took %f secs\n",dim2,t2-t1);
  fflush(stdout);
#endif

}

void stripEPInav(struct data *d)
{
  int i,j,k,l,n;
  int ix1,ix2;
  int dim1,dim2,dim3,nr;
  int ndim2,shotnv,nnav,etl,nseg,nphase,kzero,segzero;
  fftw_complex *data;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"stripEPInav"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Figure which segment kzero is in */
  nphase=(int)*val("nphase",&d->p);
  nseg=(int)*val("nseg",&d->p);
  segzero = (nphase/2-1)%nseg;
  segzero++;

  /* Calculate number of traces per shot */
  nnav=(int)*val("nnav",&d->p);
  etl=(int)*val("etl",&d->p);
  shotnv=nnav+etl;

  /* Calculate new data dimension */
  kzero=(int)*val("kzero",&d->p);
  ndim2 = nseg*(kzero-1) + segzero + nphase/2;

  /* Make data the new size */
  if ((data = (fftw_complex *)fftw_malloc(ndim2*dim1*sizeof(fftw_complex))) == NULL) nomem();

  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (n=0;n<nseg;n++) {
        ix2 = (n*(shotnv)+nnav)*dim1;
        for(k=0;k<etl;k++) {
          ix1 = k*nseg*dim1+n*dim1;
          if (ix1 >= ndim2*dim1 ) break;
          for (l=0;l<dim1;l++) {
            data[ix1][0]=d->data[i][j][ix2][0];
            data[ix1][1]=d->data[i][j][ix2][1];
            ix1++;
            ix2++;
          }
        }
      }

      /* free d->data[i][j], reallocate and copy data back in its place */
      fftw_free(d->data[i][j]);
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(ndim2*dim1*sizeof(fftw_complex))) == NULL) nomem();
      ix1 = 0;
      for(k=0;k<ndim2;k++) {
        for (l=0;l<dim1;l++) {
          d->data[i][j][ix1][0]=data[ix1][0];
          d->data[i][j][ix1][1]=data[ix1][1];
          ix1++;
        }
      }

    }
  }

  d->nv=ndim2;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Stripping to give %d x %d data: took %f secs\n",dim1,ndim2,t2-t1);
  fflush(stdout);
#endif

}

void weightnavs(struct data *d,double gf)
{
  int dim1,dim2,dim3,nr;
  int shotnv,nnav,etl,nseg;
  int i,j,k,l,n,ix;
  double f;
  double *weight;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"weightnavs"); /* Set function name */
#endif

  /* Return if not FID data and not shifted */
  if ((d->datamode != IMAGE) || !d->shift) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  nnav=(int)*val("nnav",&d->p);
  if (nnav == 0) return;
  etl=(int)*val("etl",&d->p);
  nseg=(int)*val("nseg",&d->p);
  shotnv=nnav+etl;

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Calculate weighting */
  if ((weight = (double *)malloc(dim1*sizeof(double))) == NULL) nomem();
  for (l=0;l<dim1;l++) {
    weight[l]=1.0;
  }

  /* Gaussian line broadening */
  for (l=0;l<dim1/2;l++) {
    /* Gaussian broadening as fraction of FOV */
    f=l*M_PI*gf;
    f=exp(-f*f);
    weight[l] *=f;
    weight[dim1-1-l] *=f;
  }

  /* Weight the navs */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (n=0;n<nseg;n++) {
        ix = (n*(shotnv))*dim1;
        for(k=0;k<nnav;k++) {
          for (l=0;l<dim1;l++) {
            d->data[i][j][ix][0] *=weight[l];
            d->data[i][j][ix][1] *=weight[l];
            ix++;
          }
        }
      }
    }
  }

  free(weight);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  gf nav weighting = %f\n",gf);
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fflush(stdout);
#endif

}

void segscale(struct data *d,struct segscale *scale)
{
  int nseg,nnav,etl,nv;
  int dim1,dim2,dim3,nr;
  int i,j,k,l,n,ix;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"segscale"); /* Set function name */
#endif

  /* Return if there is no scale data */
  if (scale->data==FALSE) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  nseg=(int)*val("nseg",&d->p);
  nnav=(int)*val("nnav",&d->p);             /* # navigators */
  etl=(int)*val("etl",&d->p);               /* echo train length */
  nv=nnav+etl;

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Scale the data */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (n=0;n<nseg;n++) {
        for (k=0;k<nv;k++) {
          ix=n*nv*dim1+k*dim1;
          for (l=0;l<dim1;l++) {
            d->data[i][j][ix+l][0] *=scale->value[j][n][l];
            d->data[i][j][ix+l][1] *=scale->value[j][n][l];
          }
        }
      }
    }
  }

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Scaling compressed segements: Took %f secs\n",t2-t1);
  /* Print scale values
  for (j=0;j<dim3;j++) {
    for (l=0;l<dim1;l++) {
      fprintf(stdout,"  slice %d, np %d:",j,l);
      for (k=0;k<nseg;k++) fprintf(stdout," %f",scale->value[j][k][l]);
      fprintf(stdout,"\n");
    }
  }
  */
  fflush(stdout);
#endif

}

void setsegscale(struct data *d,struct segscale *scale)
{
  int cseg,nseg,scaleseg,nnav,etl,nv;
  int dim1,dim2,dim3,nr;
  int i,j,k,l,n,ix;
  double re,im,max;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"setsegscale"); /* Set function name */
#endif

  /* Return if segments are not compressed and scaling not selected */
  cseg=spar(d,"cseg","y");
  nseg=(int)*val("nseg",&d->p);
  scaleseg=(int)*val("scaleseg",&d->p);
  if (!cseg || nseg<2 || !scaleseg) {
    scale->data=FALSE;
    return;
  }

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  nnav=(int)*val("nnav",&d->p);             /* # navigators */
  etl=(int)*val("etl",&d->p);               /* echo train length */
  nv=nnav+etl;                              /* # echoes per shot */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  if (d->datamode!=IMAGE) ftnpEPI(d); /* FT data if required */

  /* Allocate for scale values, [slice][segment][np] */
  if ((scale->value = (double ***)fftw_malloc(dim3*sizeof(double **))) == NULL) nomem();
  for (j=0;j<dim3;j++) {
    if ((scale->value[j] = (double **)fftw_malloc(nseg*sizeof(double *))) == NULL) nomem();
    for (n=0;n<nseg;n++)
      if ((scale->value[j][n] = (double *)fftw_malloc(dim1*sizeof(double))) == NULL) nomem();
  }

  /* Initialise scale values */
  for (j=0;j<dim3;j++) {
    for (n=0;n<nseg;n++) {
      for (l=0;l<dim1;l++) scale->value[j][n][l]=0.0;
    }
  }

  /* Fill scale values with magnitudes */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (n=0;n<nseg;n++) {
        for (k=0;k<nv;k++) {
          ix=n*nv*dim1+k*dim1;
          for (l=0;l<dim1;l++) {
            re=d->data[i][j][ix+l][0];
            im=d->data[i][j][ix+l][1];
            scale->value[j][n][l]+=sqrt(re*re+im*im);
          }
        }
      }
    }
  }

  /* Figure scale values */
  for (j=0;j<dim3;j++) {
    for (l=0;l<dim1;l++) {
      max=scale->value[j][0][l];
      for (n=0;n<nseg;n++)
        if (scale->value[j][n][l]>max) max=scale->value[j][n][l];
      for (n=0;n<nseg;n++) scale->value[j][n][l]=max/scale->value[j][n][l];
    }
  }

  /* Flag there is scale data */
  scale->data=TRUE;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Calculating scaling for compressed segements: Took %f secs\n",t2-t1);
  fflush(stdout);
#endif
}

void analyseEPInav(struct data *d)
{
  int dim1,dim2,dim3,nr;
  int shotnv,nnav,etl,nseg;
  int i,j,k,l,n,ix,nval;
  double re,im,M2,noisefrac,*noiselvl,*navphase,*navmag;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"analyseEPInav"); /* Set function name */
#endif

  nnav=(int)*val("nnav",&d->p);
  if (nnav == 0) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  etl=(int)*val("etl",&d->p);
  nseg=(int)*val("nseg",&d->p);
  shotnv=nnav+etl;

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Estimate noise so we can scale magnitude to it */
  noisefrac=0.05; nval=0;
  if ((noiselvl = (double *)malloc(nr*sizeof(double))) == NULL) nomem();
  for (i=0;i<nr;i++) {
    noiselvl[i]=0; nval=0;
    for (j=0;j<dim3;j++) {
      for (n=0;n<nseg;n++) {
        nval=0;
        ix = (n*(shotnv))*dim1;
        for(k=0;k<nnav;k++) {
          for (l=0;l<dim1*noisefrac/2;l++) {
            re=d->data[i][j][ix+l][0];
            im=d->data[i][j][ix+l][1];
            noiselvl[i]+=re*re+im*im;
            nval++;
          }
          for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
            re=d->data[i][j][ix+l][0];
            im=d->data[i][j][ix+l][1];
            noiselvl[i]+=re*re+im*im;
            nval++;
          }
        }
      }
    }
    noiselvl[i] /= nval;
  }

#ifdef DEBUG
  fprintf(stdout,"  %d pixels in each mavigator sampled for noise (noisefrac = %f)\n",nval/nnav,noisefrac);
  fflush(stdout);
#endif

  if ((navphase = (double *)malloc(nnav*nseg*dim1*sizeof(double))) == NULL) nomem();
  if ((navmag = (double *)malloc(nnav*nseg*dim1*sizeof(double))) == NULL) nomem();
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      nval=0;
      for (n=0;n<nseg;n++) {
        ix = (n*(shotnv))*dim1;
        for(k=0;k<nnav;k++) {
          for (l=0;l<dim1;l++) {
            re=d->data[i][j][ix][0];
            im=d->data[i][j][ix][1];
            M2=re*re+im*im;
            navphase[nval]=atan2(im,re);
            navmag[nval]=sqrt(M2);
            nval++;
            ix++;
          }
        }
      }
    }
    for (j=0;j<dim3;j++) {
      for (n=0;n<nseg;n++) {
        fprintf(stdout,"Slice %d seg %d:\n",j+1,n+1);
        for (l=0;l<dim1;l++) {
          fprintf(stdout,"[%5.1f] ",navmag[l]/sqrt(noiselvl[i]));
          for(k=0;k<nnav;k++) {
            ix=(n*nnav+k)*dim1+l;
            fprintf(stdout,"%6.1f ",navphase[ix]*180/M_PI);
          }
          fprintf(stdout,"\n");
        }
      }
    }
  }
  free(navphase);
  free(navmag);
  free(noiselvl);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fflush(stdout);
#endif
}

void getblockEPI(struct data *d,int volindex,int DCCflag)
{
  /* Check d->nvols has been set */
  if (d->nvols == 0) setnvolsEPI(d); /* Set the number of data volumes */

  /* Set start and end position of block */
  setblock(d,d->ns);

  /* Get the block */
  getblock(d,volindex,DCCflag);

}

void setblockEPI(struct data *d)
{
  int i,j,k,l,m,n;
  int ix1,ix2;
  int dim1,dim2,dim3,nr;
  int ndim1,ndim2,nnav,etl,nv,nseg,gradnp,rampnp,platnp,oversample,skipnp=0,offset,echo;
  int prescan=FALSE,grid=FALSE;
  int ns,slice,slix,altread,GE,cseg=FALSE,rseg;
  fftw_complex *data;
  fftw_complex **cdata;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"setblockEPI"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Check for prescan, ramp sampling, linear sampling and gridding */
  if (spar(d,"recon","prescanEPI")) prescan=TRUE;
  if (spar(d,"grid","y")) grid=TRUE;
  if (prescan) grid=FALSE; /* For prescan don't grid */

  /* New data dimensions */
  platnp=(int)(*val("platnp",&d->p));          /* # points on the readout plateau */
  rampnp=(int)(*val("rampnp",&d->p));          /* # points on a readout ramp */
  oversample=(int)(*val("oversample",&d->p));  /* oversample factor */
  platnp*=oversample;                          /* correct for oversample */
  rampnp*=oversample;                          /* correct for oversample */
  gradnp=platnp+2*rampnp;                      /* # points in each readout gradient */
  ndim1=gradnp;                                /* new dim1 */
  nnav=(int)*val("nnav",&d->p);                /* # navigators */
  etl=(int)*val("etl",&d->p);                  /* echo train length */
  nseg=(int)*val("nseg",&d->p);                /* # segments */
  ndim2=nseg*(nnav+etl);                       /* new dim2 */

  nv=nnav+etl;                                 /* # echoes per shot */
  ns=nvals("pss",&d->p);                       /* # slices */
  altread=spar(d,"altread","y");               /* alternating read gradient for alternate segments ... */
  if (nseg<2) altread=FALSE;                   /* ... only if nseg>1 ... */
  if (nseg%2==0) altread=FALSE;                /* ... and only if nseg is odd */

  GE=spar(d,"spinecho","n");                   /* gradient echo */
  if (GE) {
    cseg=spar(d,"cseg","y");                   /* compressed segments only for gradient echo */
    rseg=spar(d,"rseg","y");                   /* compressed reference segments reversed */
  }

  if (cseg) { /* compressed segments */

    /* Make data the new size */
    if ((cdata = (fftw_complex **)fftw_malloc(dim3*sizeof(fftw_complex))) == NULL) nomem();
    for (j=0;j<dim3;j++)
      if ((cdata[j] = (fftw_complex *)fftw_malloc(ndim2*ndim1*sizeof(fftw_complex))) == NULL) nomem();

    for (i=0;i<nr;i++) {
      for (j=0;j<dim3;j++) {
        slix=sliceindex(d,j);
        for (n=0;n<nseg;n++) {

          if (rseg) {   /* Reference segments reversed */
            switch ((int)getelem(d,"image",d->vol)) {
              case -1:  /* Inverted Readout Reference, with phase-encode */
                m=nseg-1-n;
                break;
              case -2:  /* Inverted Readout Reference, no phase-encode */
                m=nseg-1-n;
                break;
              default:
                m=n;
                break;
            }
          } else m=n;

          slice=(slix*nseg+n)%ns;
          if (d->pssorder[0]>=0) slice=d->pssorder[slice];
          offset=(slix*nseg+n)/ns;
          offset*=dim1;

          echo=0;
          if (altread && n%2) echo=1;
          for(k=0;k<nnav;k++) {
            if (echo%2==0) {
              ix1=(m*nv+k)*ndim1;
              ix2=k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                cdata[j][ix1][0]=d->data[i][slice][ix2][0];
                cdata[j][ix1][1]=d->data[i][slice][ix2][1];
                ix1++;
                ix2++;
              }
            } else {
              ix1=(m*nv+k)*ndim1;
              if (grid) ix2=k*ndim1+offset; else ix2=ndim1-1+k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                cdata[j][ix1][0]=d->data[i][slice][ix2][0];
                cdata[j][ix1][1]=d->data[i][slice][ix2][1];
                ix1++;
                if (grid) ix2++; else ix2--;
              }
            }
            offset += 2*skipnp;
            echo++;
          }

          /* Skip the extra readout gradient */
          offset += gradnp;
          echo++;

          for(k=nnav;k<nv;k++) {
            if (echo%2==0) {
              ix1=(m*nv+k)*ndim1;
              ix2=k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                cdata[j][ix1][0]=d->data[i][slice][ix2][0];
                cdata[j][ix1][1]=d->data[i][slice][ix2][1];
                ix1++;
                ix2++;
              }
            } else {
              ix1=(m*nv+k)*ndim1;
              if (grid) ix2=k*ndim1+offset; else ix2=ndim1-1+k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                cdata[j][ix1][0]=d->data[i][slice][ix2][0];
                cdata[j][ix1][1]=d->data[i][slice][ix2][1];
                ix1++;
                if (grid) ix2++; else ix2--;
              }
            }
            offset += 2*skipnp;
            echo++;
          }

        } /* end nseg loop */
      }
      for (j=0;j<dim3;j++) {
        /* free d->data[i][j], reallocate and copy data back in its place */
        fftw_free(d->data[i][j]);
        if ((d->data[i][j] = (fftw_complex *)fftw_malloc(ndim2*ndim1*sizeof(fftw_complex))) == NULL) nomem();
        for(k=0;k<ndim2;k++) {
          for (l=0;l<ndim1;l++) {
            ix1=k*ndim1+l;
            d->data[i][j][ix1][0]=cdata[j][ix1][0];
            d->data[i][j][ix1][1]=cdata[j][ix1][1];
          }
        }
      }
    }

    for (j=0;j<dim3;j++) free(cdata[j]);
    free(cdata);

  } else { /* Regular segments */

    /* Make data the new size */
    if ((data = (fftw_complex *)fftw_malloc(ndim2*ndim1*sizeof(fftw_complex))) == NULL) nomem();

    for (i=0;i<nr;i++) {
      for (j=0;j<dim3;j++) {

        for (n=0;n<nseg;n++) {

          offset=skipnp;
          echo=0;
          if (altread && n%2) echo=1;
          for(k=0;k<nnav;k++) {
            if (echo%2==0) {
              ix1=(n*nv+k)*ndim1;
              ix2=n*dim1+k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                data[ix1][0]=d->data[i][j][ix2][0];
                data[ix1][1]=d->data[i][j][ix2][1];
                ix1++;
                ix2++;
              }
            } else {
              ix1=(n*nv+k)*ndim1;
              if (grid) ix2=n*dim1+k*ndim1+offset; else ix2=ndim1-1+n*dim1+k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                data[ix1][0]=d->data[i][j][ix2][0];
                data[ix1][1]=d->data[i][j][ix2][1];
                ix1++;
                if (grid) ix2++; else ix2--;
              }
            }
            offset += 2*skipnp;
            echo++;
          }

          /* Skip the extra readout gradient */
          offset += gradnp;
          echo++;

          for(k=nnav;k<nv;k++) {
            if (echo%2==0) {
              ix1=(n*nv+k)*ndim1;
              ix2=n*dim1+k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                data[ix1][0]=d->data[i][j][ix2][0];
                data[ix1][1]=d->data[i][j][ix2][1];
                ix1++;
                ix2++;
              }
            } else {
              ix1=(n*nv+k)*ndim1;
              if (grid) ix2=n*dim1+k*ndim1+offset; else ix2=ndim1-1+n*dim1+k*ndim1+offset;
              for (l=0;l<ndim1;l++) {
                data[ix1][0]=d->data[i][j][ix2][0];
                data[ix1][1]=d->data[i][j][ix2][1];
                ix1++;
                if (grid) ix2++; else ix2--;
              }
            }
            offset += 2*skipnp;
            echo++;
          }

        } /* end nseg loop */

        /* free d->data[i][j], reallocate and copy data back in its place */
        fftw_free(d->data[i][j]);
        if ((d->data[i][j] = (fftw_complex *)fftw_malloc(ndim2*ndim1*sizeof(fftw_complex))) == NULL) nomem();
        for(k=0;k<ndim2;k++) {
          for (l=0;l<ndim1;l++) {
            ix1=k*ndim1+l;
            d->data[i][j][ix1][0]=data[ix1][0];
            d->data[i][j][ix1][1]=data[ix1][1];
          }
        }

      } /* end dim3 (slice) loop */
    } /* end nr loop */

    free(data);

  } /* end regular segment */

  d->np=ndim1*2; /* set d->np for the data */
  d->nv=ndim2;   /* set d->nv for the data */

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Reordering to give %d x %d data: took %f secs\n",ndim1,ndim2,t2-t1);
  fflush(stdout);
#endif

}

void setoutvolEPI(struct data *d)
{
  switch ((int)getelem(d,"image",d->vol)) {
    case 1:
      d->outvol++;
      if (d->outvol==1) return; /* We want to process the first actual scan */
      /* Advance to requested start volume if we are not already there */
      while (outvolEPI(d) && (d->outvol<d->startvol)) {
        d->vol++;
        d->outvol++;
      }
      break;
    case -3:
      if (spar(d,"imSGE","y")) d->outvol++;
      break;
    case -4:
      if (spar(d,"imSGE","y")) d->outvol++;
      break;
  }
}

int outvolEPI(struct data *d)
{
  switch ((int)getelem(d,"image",d->vol)) {
    case 1: return(1); break;
    default: return(0); break;
  }
}

void setnvolsEPI(struct data *d)
{
  int dim1,dim2,dim3,nr;
#ifdef DEBUG
  char function[20];
  strcpy(function,"setnvolsEPI"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fflush(stdout);
#endif
  /* Set data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->ns; nr=d->nr;

  /* Set nuber of "volumes" */
  d->nvols=d->fh.nblocks/nr;

#ifdef DEBUG
  fprintf(stdout,"  Number of volumes = %d\n",d->nvols);
  fflush(stdout);
#endif

  if (spar(d,"allvolumes","n")) {  /* Don't process all volumes */
    d->startvol=(int)*val("startvol",&d->p);
    d->endvol=(int)*val("endvol",&d->p);
    if (d->startvol > d->endvol) {
      /* Swap them */
      d->vol=d->endvol; d->endvol=d->startvol; d->startvol=d->vol; d->vol=0;
    }
    if (d->startvol < 0) d->startvol=0;
    if (d->endvol > d->nvols) d->endvol=d->nvols;
  } else { /* Process all volumes */
    d->startvol=0;
    d->endvol=d->nvols;
  }

  /* Set the sequence mode from seqcon and apptype parameters */
  setseqmode(d);
}

void setblockSGE(struct data *d)
{
  int i,j,k,l,n;
  int ix1,ix2,slice,slix,peix;
  int dim1,dim2,dim3,nr;
  int ndim1,ndim2,nnav,etl,ntraces,platnp,rampnp,oversample,nseg,ns,image;
  fftw_complex **data;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"setblockSGE"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* New data dimensions */
  platnp=(int)(*val("platnp",&d->p));          /* # points on the readout plateau */
  rampnp=(int)(*val("rampnp",&d->p));          /* # points on a readout ramp */
  oversample=(int)(*val("oversample",&d->p));  /* oversample factor */
  platnp*=oversample;                          /* correct for oversample */
  rampnp*=oversample;                          /* correct for oversample */
  ndim1=platnp+2*rampnp;                       /* # points in each readout gradient */
  ndim2=(int)(*val("nphase",&d->p));           /* # points in phase */

  nnav=(int)*val("nnav",&d->p);                /* # navigators */
  etl=(int)*val("etl",&d->p);                  /* echo train length */
  nseg=(int)*val("nseg",&d->p);                /* # segments */
  ntraces=nnav+etl+1;                          /* # traces per segment */
  ns=nvals("pss",&d->p);                       /* # slices */

  image=(int)getelem(d,"image",d->vol);

  /* Make data the new size */
  if ((data = (fftw_complex **)fftw_malloc(dim3*sizeof(fftw_complex))) == NULL) nomem();
  for (j=0;j<dim3;j++)
    if ((data[j] = (fftw_complex *)fftw_malloc(ndim2*ndim1*sizeof(fftw_complex))) == NULL) nomem();

  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      slix=sliceindex(d,j);
      for (n=0;n<nseg;n++) {
        switch (image) {
          case -3:
            for(k=0;k<ntraces;k++) {
              ix1=n*ntraces+k;
              if (ix1>=ndim2) break;
              peix=phase2index(d,ix1);
              slice=((ns*peix+slix)/ntraces)%ns;
              if (d->pssorder[0]>=0) slice=d->pssorder[slice];
              ix2=(ns*peix+slix)%ntraces +ntraces*((ns*peix+slix)/(ns*ntraces));
              ix1 *=ndim1;
              ix2 *=ndim1;
              for (l=0;l<ndim1;l++) {
                data[j][ix1][0]=d->data[i][slice][ix2][0];
                data[j][ix1][1]=d->data[i][slice][ix2][1];
                ix1++;
                ix2++;
              }
            }
            break;
          case -4:
            for(k=0;k<ntraces;k++) {
              ix1=n*ntraces+k;
              if (ix1>=ndim2) break;
              peix=phase2index(d,ix1);
              slice=(slix*nseg+peix/ntraces)%ns;
              if (d->pssorder[0]>=0) slice=d->pssorder[slice];
              ix2=(slix*nseg*ntraces+peix)%ntraces+ntraces*((slix*nseg*ntraces+peix)/(ns*ntraces));
              ix1 *=ndim1;
              ix2 *=ndim1;
              for (l=0;l<ndim1;l++) {
                data[j][ix1][0]=d->data[i][slice][ix2][0];
                data[j][ix1][1]=d->data[i][slice][ix2][1];
                ix1++;
                ix2++;
              }
            }
            break;
        }
      } /* end nseg loop */
      /* Zero fill if we don't have full k-space coverage */
      for(k=nseg*ntraces;k<ndim2;k++) {
        ix1=k*ndim1;
        for (l=0;l<ndim1;l++) {
          data[j][ix1][0]=0.0;
          data[j][ix1][1]=0.0;
          ix1++;
          ix2++;
        }
      }
    }
    for (j=0;j<dim3;j++) {
      /* free d->data[i][j], reallocate and copy data back in its place */
      fftw_free(d->data[i][j]);
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(ndim2*ndim1*sizeof(fftw_complex))) == NULL) nomem();
      for(k=0;k<ndim2;k++) {
        for (l=0;l<ndim1;l++) {
          ix1=k*ndim1+l;
          d->data[i][j][ix1][0]=data[j][ix1][0];
          d->data[i][j][ix1][1]=data[j][ix1][1];
        }
      }
    }
  }

  for (j=0;j<dim3;j++) free(data[j]);
  free(data);

  d->np=ndim1*2;
  d->nv=ndim2;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Reordering to give %d x %d data: took %f secs\n",ndim1,ndim2,t2-t1);
  fflush(stdout);
#endif
}

/* noise2D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* noise2D.c: 2D routines based on noise measurements                        */
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

#define SOURCEFILE "common2D/noise2D.c"

/*-----------------------*/
/*---- Some defaults ----*/
/*-----------------------*/
/* DCCnoisefraction: Fraction of FOV to use to sample noise for DC correction */
static double DCCnoisefraction=0.25;

/* Knoisefraction: Fraction of k-space FOV to use to sample noise */
static double Knoisefraction=0.05;

/* IMnoisefraction: Fraction of image space FOV to use to sample noise */
static double IMnoisefraction=0.05;

void dccorrect2D(struct data *d)
{
  int dim1,dim2,dim3,nr;
  double noisefrac;
  int i,j,k,l;
  int n;
  double redc,imdc;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"dccorrect2D"); /* Set function name */
#endif

  /* Assume we are applying to FID data that has been OPT shifted */
  if ((d->datamode == FID) && d->shift) {

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

    /* Data dimensions */
    dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

    noisefrac=DCCnoisefraction;

    for (i=0;i<nr;i++) {
      n=0;
      redc=0.0;
      imdc=0.0;
      for (j=0;j<dim3;j++) {
        for(k=dim2/2-dim2*noisefrac/2;k<dim2/2+dim2*noisefrac/2;k++) {
          for (l=dim1/2-dim1*noisefrac/2;l<dim1/2+dim1*noisefrac/2;l++) {
            redc+=d->data[i][j][k*dim1+l][0];
            imdc+=d->data[i][j][k*dim1+l][1];
            n++;
          }
        }
      }
      redc/=n;
      imdc/=n;
      for (j=0;j<dim3;j++) {
        for(k=0;k<dim2;k++) {
          for (l=0;l<dim1;l++) {
            d->data[i][j][k*dim1+l][0]-=redc;
            d->data[i][j][k*dim1+l][1]-=imdc;
          }
        }
      }
    }

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  DC correction using noise region in data: took %f secs\n",t2-t1);
  fflush(stdout);
#endif

  }

}

void getnoise2D(struct data *d,int mode)
{
  int dim1,dim2,dim3,nr;
  double noisefrac;
  int i,j,k,l;
  double re,im,M2;
  int n=0;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"getnoise2D"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
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

  zeronoise(d);
  for (i=0;i<nr;i++) {
    n=0;
    if (((d->datamode == FID) && d->shift) || ((d->datamode == IMAGE) && !d->shift)) {
      /* For shifted FID and not-shifted IMAGE the noise is at centre */
      for (j=0;j<dim3;j++) {
        for(k=dim2/2-dim2*noisefrac/2;k<dim2/2+dim2*noisefrac/2;k++) {
          for (l=dim1/2-dim1*noisefrac/2;l<dim1/2+dim1*noisefrac/2;l++) {
            re=d->data[i][j][k*dim1+l][0];
            im=d->data[i][j][k*dim1+l][1];
            M2=re*re+im*im;
            d->noise.M[i]+=sqrt(M2);
            d->noise.M2[i]+=M2;
            d->noise.Re[i]+=fabs(re);
            d->noise.Im[i]+=fabs(im);
            n++;
          }
        }
      }
    } else {
      /* For not shifted FID and shifted IMAGE the noise is at edges */
      for (j=0;j<dim3;j++) {
        for(k=0;k<dim2*noisefrac/2;k++) {
          for (l=0;l<dim1*noisefrac/2;l++) {
            re=d->data[i][j][k*dim1+l][0];
            im=d->data[i][j][k*dim1+l][1];
            M2=re*re+im*im;
            d->noise.M[i]+=sqrt(M2);
            d->noise.M2[i]+=M2;
            d->noise.Re[i]+=fabs(re);
            d->noise.Im[i]+=fabs(im);
            n++;
          }
          for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
            re=d->data[i][j][k*dim1+l][0];
            im=d->data[i][j][k*dim1+l][1];
            M2=re*re+im*im;
            d->noise.M[i]+=sqrt(M2);
            d->noise.M2[i]+=M2;
            d->noise.Re[i]+=fabs(re);
            d->noise.Im[i]+=fabs(im);
            n++;
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
            n++;
          }
          for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
            re=d->data[i][j][k*dim1+l][0];
            im=d->data[i][j][k*dim1+l][1];
            M2=re*re+im*im;
            d->noise.M[i]+=sqrt(M2);
            d->noise.M2[i]+=M2;
            d->noise.Re[i]+=fabs(re);
            d->noise.Im[i]+=fabs(im);
            n++;
          }
        }
      }
    }
    /* calculate the mean */
    d->noise.M[i] /=n;
    d->noise.M2[i] /=n;
    /* For Real and Imaginary we must consider console type.
       The DDR in VNMRS produces equal noise Re and Im channels - no quad images */
    if (!(strcmp(*sval("console",&d->p),"vnmrs"))) { /* VNMRS */
      d->noise.Re[i] += d->noise.Im[i];
      d->noise.Re[i] /=2.0;
      d->noise.Im[i] = d->noise.Re[i];
    }
    d->noise.Re[i] /=n;
    d->noise.Im[i] /=n;
  }

  /* Now average over all receivers */
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
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Noise data averaged over %d points: took %f secs\n",n,t2-t1);
  for (i=0;i<nr;i++) {
  fprintf(stdout,"  Receiver %d: M = %.3f, M2 = %.3f, Re = %.3f, Im = %.3f\n",
    i,d->noise.M[i],d->noise.M2[i],d->noise.Re[i],d->noise.Im[i]);
  }
  fprintf(stdout,"  Average:    M = %.3f, M2 = %.3f, Re = %.3f, Im = %.3f\n",
    d->noise.avM,d->noise.avM2,d->noise.avRe,d->noise.avIm);
  fflush(stdout);
#endif

}

void get2Dnoisematrix(struct data *d,int mode)
{
  int dim1,dim2,dim3,nr;
  double noisefrac;
  int h,i,j,k,l;
  int n=0;
  gsl_complex cx,cx1,cx2;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"get2Dnoisematrix"); /* Set function name */
#endif

  /* Get noise if it has not been measured */
  if (!d->noise.data) getnoise2D(d,mode);

  /* Equalize noise if it has not been done */
  if (!d->noise.equal) equalizenoise2D(d,mode);

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
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

  /* Allocate memory for noise matrix */   /* Get noise if it has not been measured */
  if (!d->noise.data) getnoise2D(d,mode);

  if ((d->noise.mat = (gsl_matrix_complex **)malloc(dim3*sizeof(gsl_matrix_complex *))) == NULL) nomem();
  for (j=0;j<dim3;j++)
    d->noise.mat[j]=(gsl_matrix_complex *)gsl_matrix_complex_calloc(nr,nr);

  for (h=0;h<nr;h++) {
    for (i=0;i<nr;i++) {
      if (((d->datamode == FID) && d->shift) || ((d->datamode == IMAGE) && !d->shift)) {
        /* For shifted FID and not-shifted IMAGE the noise is at centre */
        for (j=0;j<dim3;j++) {
          n=0;
          GSL_SET_COMPLEX(&cx,0.0,0.0);
          for(k=dim2/2-dim2*noisefrac/2;k<dim2/2+dim2*noisefrac/2;k++) {
            for (l=dim1/2-dim1*noisefrac/2;l<dim1/2+dim1*noisefrac/2;l++) {
              GSL_SET_REAL(&cx1,d->data[h][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx1,d->data[h][j][k*dim1+l][1]);
              GSL_SET_REAL(&cx2,d->data[i][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx2,d->data[i][j][k*dim1+l][1]);
              cx=gsl_complex_add(cx,gsl_complex_mul(cx1,gsl_complex_conjugate(cx2)));
              n++;
            }
          }
          /* Take the mean and normalize */
          GSL_SET_COMPLEX(&cx,GSL_REAL(cx)/(n*d->noise.avM2),GSL_IMAG(cx)/(n*d->noise.avM2));
          gsl_matrix_complex_set(d->noise.mat[j],h,i,cx);
        }
      } else {
        /* For not shifted FID and shifted IMAGE the noise is at edges */
        for (j=0;j<dim3;j++) {
          n=0;
          GSL_SET_COMPLEX(&cx,0.0,0.0);
           for(k=0;k<dim2*noisefrac/2;k++) {
            for (l=0;l<dim1*noisefrac/2;l++) {
              GSL_SET_REAL(&cx1,d->data[h][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx1,d->data[h][j][k*dim1+l][1]);
              GSL_SET_REAL(&cx2,d->data[i][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx2,d->data[i][j][k*dim1+l][1]);
              cx=gsl_complex_add(cx,gsl_complex_mul(cx1,gsl_complex_conjugate(cx2)));
              n++;
            }
            for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
              GSL_SET_REAL(&cx1,d->data[h][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx1,d->data[h][j][k*dim1+l][1]);
              GSL_SET_REAL(&cx2,d->data[i][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx2,d->data[i][j][k*dim1+l][1]);
              cx=gsl_complex_add(cx,gsl_complex_mul(cx1,gsl_complex_conjugate(cx2)));
              n++;
            }
          }
          for(k=dim2-dim2*noisefrac/2;k<dim2;k++) {
            for (l=0;l<dim1*noisefrac/2;l++) {
              GSL_SET_REAL(&cx1,d->data[h][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx1,d->data[h][j][k*dim1+l][1]);
              GSL_SET_REAL(&cx2,d->data[i][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx2,d->data[i][j][k*dim1+l][1]);
              cx=gsl_complex_add(cx,gsl_complex_mul(cx1,gsl_complex_conjugate(cx2)));
              n++;
            }
            for (l=dim1-dim1*noisefrac/2;l<dim1;l++) {
              GSL_SET_REAL(&cx1,d->data[h][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx1,d->data[h][j][k*dim1+l][1]);
              GSL_SET_REAL(&cx2,d->data[i][j][k*dim1+l][0]);
              GSL_SET_IMAG(&cx2,d->data[i][j][k*dim1+l][1]);
              cx=gsl_complex_add(cx,gsl_complex_mul(cx1,gsl_complex_conjugate(cx2)));
              n++;
            }
          }
          /* Take the mean and normalize */
          GSL_SET_COMPLEX(&cx,GSL_REAL(cx)/(n*d->noise.avM2),GSL_IMAG(cx)/(n*d->noise.avM2));
          gsl_matrix_complex_set(d->noise.mat[j],h,i,cx);
        }
      }
    }
  }

  /* Set noise matrix flag */
  d->noise.matrix=TRUE;

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  print2Dnoisematrix(d);
  fflush(stdout);
#else
  /* Print Noise Matrix, if requested */
  if (!(strcmp(*sval("printNM",&d->p),"y")))  print2Dnoisematrix(d);
#endif

}

void equalizenoise2D(struct data *d,int mode)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l;
  double *Rsf,*Isf;
  double maxval;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"equalizenoise2D"); /* Set function name */
#endif

  switch(mode) {
    case MK: /* Mask parameters */
      /* Return unless maskeqnoise='y' */
      if (strcmp(*sval("maskeqnoise",&d->p),"y")) return;
      break;
    case SM: /* Sensitivity map parameters */
      /* Return unless maskeqnoise='y' */
      if (strcmp(*sval("smapeqnoise",&d->p),"y")) return;
      break;
    default: /* Default parameters */
      /* Return unless eqnoise='y' */
      if (strcmp(*sval("eqnoise",&d->p),"y")) return;
      break;
  } /* end mode switch */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Return if only one receiver */
  if (nr < 2) return;

  /* Return if noise has not been measured */
  if (!d->noise.data) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Get the maximum noise level */
  maxval=d->noise.Re[0];
  for (i=1;i<nr;i++) if (d->noise.Re[i] > maxval) maxval=d->noise.Re[i];
  for (i=0;i<nr;i++) if (d->noise.Im[i] > maxval) maxval=d->noise.Im[i];

  /* Find scale factors */
  if ((Rsf = (double *)malloc(nr*sizeof(double))) == NULL) nomem();
  if ((Isf = (double *)malloc(nr*sizeof(double))) == NULL) nomem();
  for (i=0;i<nr;i++) Rsf[i] = maxval/d->noise.Re[i];
  for (i=0;i<nr;i++) Isf[i] = maxval/d->noise.Im[i];

  /* Scale to equalize noise */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for(k=0;k<dim2;k++) {
        for (l=0;l<dim1;l++) {
          d->data[i][j][k*dim1+l][0] *=Rsf[i];
          d->data[i][j][k*dim1+l][1] *=Isf[i];
        }
      }
    }
  }

  /* Update noise values */
/*  getnoise2D(d,mode); */

  /* Set equalize flag */
  d->noise.equal=TRUE;

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Noise data equlaized: took %f secs\n",t2-t1);
  fflush(stdout);
#endif

}

void scale2Ddata(struct data *d,double factor)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  char function[20];
  strcpy(function,"scale2Ddata"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Scale data */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for(k=0;k<dim2;k++) {
        for (l=0;l<dim1;l++) {
          d->data[i][j][k*dim1+l][0] *=factor;
          d->data[i][j][k*dim1+l][1] *=factor;
        }
      }
    }
  }

#ifdef DEBUG
  gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Data scaled by %f: took %f secs\n",factor,t2-t1);
  fflush(stdout);
#endif

  /* Zero noise data */
  zeronoise(d);

}

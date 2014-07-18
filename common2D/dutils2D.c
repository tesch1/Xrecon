/* dutils2D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dutils2D.c: 2D data utilities                                             */
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

#define SOURCEFILE "common2D/dutils2D.c"

int check2Dref(struct data *d,struct data *ref)
{
  char function[20];
  strcpy(function,"check2Dref"); /* Set function name */

  /* The data sets must have the same FOV orientation and slices */

  /* Check FOV */
  if (!checkequal(d,ref,"lro","RO length")) return(FALSE);
  if (!checkequal(d,ref,"lpe","PE length")) return(FALSE);
  if (!checkequal(d,ref,"pro","RO positions")) return(FALSE);
  if (!checkequal(d,ref,"ppe","PE positions")) return(FALSE);

  /* Check orientations */
  if (!checkequal(d,ref,"psi","orientations")) return(FALSE);
  if (!checkequal(d,ref,"phi","orientations")) return(FALSE);
  if (!checkequal(d,ref,"theta","orientations")) return(FALSE);

  /* Check slices */
  if (!checkequal(d,ref,"pss","slices")) return(FALSE);

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Reference %s compatible with data %s\n",ref->procpar,d->procpar);
  fflush(stdout);
#endif
  return(TRUE);

}

void copy2Ddata(struct data *d1,struct data *d2)
{
  int dim1,dim2,dim3,nr;
  int i,j,k;
  double *dp1,*dp2;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"copy2Ddata"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d1->np/2; dim2=d1->nv; dim3=d1->endpos-d1->startpos; nr=d1->nr;

  if (d2->datamode != NONE) clear2Ddata(d2);
  /*  if ((d2->datamode == FID) || (d2->datamode == IMAGE)) clear2Ddata(d2); */

  /* Allocate memory according to nr */
  if ((d2->data = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();
  for (i=0;i<nr;i++) { /* loop over receiver blocks */
    if ((d2->data[i] = (fftw_complex **)fftw_malloc(dim3*sizeof(fftw_complex *))) == NULL) nomem();
    for (j=0;j<dim3;j++) /* loop over slices */
      if ((d2->data[i][j] = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();
  }

  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      dp1 = *d1->data[i][j];
      dp2 = *d2->data[i][j];
      for(k=0;k<dim2*d1->np;k++) *dp2++ = *dp1++;
    }
  }

  /* Set dimensions and some detail as well */
  d2->np=d1->np; d2->nv=d1->nv; d2->ns=d1->ns; d2->nr=d1->nr;
  d2->fn=d1->fn; d2->fn1=d1->fn1;
  d2->datamode=d1->datamode;
  d2->shift=d1->shift;
  d2->zerofill=d1->zerofill;
  d2->vol=d1->vol;
  d2->nblocks=d1->nblocks; d2->block=d1->block;
  d2->startpos=d1->startpos; d2->endpos=d1->endpos;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Copying data: took  %f secs\n",t2-t1);
  fflush(stdout);
#endif

}

void print2Dnoisematrix(struct data *d)
{
  int dim3;
  int h,i,j;
  gsl_complex cx;
  char function[20];
  strcpy(function,"print2Dnoisematrix"); /* Set function name */

  if (d->noise.matrix) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    dim3=d->endpos-d->startpos;
    for (j=0;j<dim3;j++) {
      fprintf(stdout,"  Noise Matrix, slice %d:\n",j+d->startpos+1);
      for (h=0;h<d->nr;h++) {
        for (i=0;i<d->nr;i++) {
          cx=gsl_matrix_complex_get(d->noise.mat[j],h,i);
          fprintf(stdout,"  [%d,%d] = %f , %fi\n",h,i,GSL_REAL(cx),GSL_IMAG(cx));
        }
      }
    }
  }
}

void setdatapars2D(struct data *d)
{

#ifdef DEBUG
  char function[20];
  strcpy(function,"setdatapars2D"); /* Set function name */
#endif

  /* There can not be more blocks than slices */
  if (d->nblocks>d->ns) d->nblocks=d->ns;

  /* Set default start slice and endslice for one block */
  d->startpos=0;
  d->endpos=d->ns;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
#endif
}

void init2Ddata(struct data *d)
{
  /* Set some defaults */
  d->nvols=0;
  d->datamode=NONE;
  d->shift=FALSE;
  d->noise.data=FALSE;
  d->noise.matrix=FALSE;
  d->maskdata=FALSE;
  d->zerofill=FALSE;
  d->dimorder=FALSE;
}

void clear2Dall(struct data *d)
{
  clear2Ddata(d);
  clear2Dmask(d);
  clearnoise2D(d);
  cleardimorder(d);
  clearpars(&d->p);
  clearpars(&d->s);
  cleararray(&d->a);
}

void clear2Ddata(struct data *d)
{
  int dim3;
  int i,j;
  if (d->datamode != NONE) {
    dim3=d->endpos-d->startpos;
    for (i=0;i<d->nr;i++) {
      for (j=0;j<dim3;j++) fftw_free(d->data[i][j]);
      fftw_free(d->data[i]);
    }
    fftw_free(d->data);
    /* Set data mode */
    d->datamode=NONE;
  }
}

void zero2Dnoisematrix(struct data *d)
{
  int dim3;
  int h,i,j;
  gsl_complex cx;
  if (d->noise.matrix) {
    GSL_SET_COMPLEX(&cx,0.0,0.0);
    dim3=d->endpos-d->startpos;
    for (j=0;j<dim3;j++) {
      for (h=0;h<d->nr;h++) {
        for (i=0;i<d->nr;i++) {
          gsl_matrix_complex_set(d->noise.mat[j],h,i,cx);
        }
      }
    }
  }
  /* Set noise matrix flag */
  d->noise.matrix=FALSE;
}

void clearnoise2D(struct data *d)
{
  int dim3;
  int j;
  if (d->noise.data) {
    free(d->noise.M);
    free(d->noise.M2);
    free(d->noise.Re);
    free(d->noise.Im);
  }
  if (d->noise.matrix) {
    dim3=d->endpos-d->startpos;
    for (j=0;j<dim3;j++) gsl_matrix_complex_free(d->noise.mat[j]);
    free(d->noise.mat);
  }
  /* Set noise matrix flag */
  d->noise.matrix=FALSE;
  /* Set data flag */
  d->noise.data=FALSE;
}

void clear2Dmask(struct data *d)
{
  int dim3;
  int i;
  if (d->maskdata) {
    dim3=d->endpos-d->startpos;
    for (i=0;i<dim3;i++) free(d->mask[i]);
    free(d->mask);
    /* Set data flag */
    d->maskdata=FALSE;
  }
}

void cleardimorder(struct data *d)
{
  switch (d->dimorder) {
    case IM2D:
      free(d->dim2order);
      free(d->dim3order);
      break;
    case IM3D:
      free(d->dim2order);
      free(d->dim3order);
      free(d->pssorder);
      break;
  }
  d->dimorder=FALSE;
}

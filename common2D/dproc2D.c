/* dproc2D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dproc2D.c: 2D Data processing routines                                    */
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

#define SOURCEFILE "common2D/dproc2D.c"

void fft2D(struct data *d)
{
  fftw_complex *slice;
  fftw_plan p;
  int dim1,dim2,dim3,nr;
  int i,j,k;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"fft2D"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

/*
  fftw_import_wisdom_from_file(FILE *input_file);
*/

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Allocate memory for slice data */
  if ((slice = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Measured plan for ft2d */
  p=fftw_plan_dft_2d(dim2,dim1,slice,slice,FFTW_FORWARD,FFTW_MEASURE);

  /* ft2d slices ... */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for (k=0;k<dim1*dim2;k++) { slice[k][0]=d->data[i][j][k][0]; slice[k][1]=d->data[i][j][k][1]; }
      /* Estimated plan */
/*
      p=fftw_plan_dft_2d(dim2,dim1,d->data[i][j],d->data[i][j],FFTW_FORWARD,FFTW_ESTIMATE);
*/
      /* ft2d */
      fftw_execute(p);
      for (k=0;k<dim1*dim2;k++) { d->data[i][j][k][0]=slice[k][0]; d->data[i][j][k][1]=slice[k][1]; }
    }
  }
  /* ... and tidy up */
  fftw_destroy_plan(p);
  fftw_free(slice);

  /* Set datamode for image data */
  d->datamode=IMAGE;

  /* Update shifted flag */
  d->shift=FALSE; /* not shifted */

  /* Reset max pars to default values */
/*  zeromax(d); */

  /* Zero noise */
/*  zeronoise(d); */

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  ft2d (%d x %d): %d slice(s), %d receiver(s): took %f secs\n",
    dim2,dim1,dim3,nr,t2-t1);
  fflush(stdout);
#endif

}

int shiftdata2D(struct data *d,int mode)
{
  int dim1,dim2;
  int npshft=0,nvshft=0;
  char function[20];
  strcpy(function,"shiftdata2D"); /* Set function name */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv;

  /* Set shift points */
  switch(mode) {
    case STD:
      npshft=dim1/2;
      nvshft=dim2/2;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Standard shift: npshft = %d, nvshft = %d\n",npshft,nvshft);
  fflush(stdout);
#endif
      break;
    case OPT:
      /* getmax2D sets maximum and its coordinates */
      if (!d->max.data) getmax2D(d);
      npshft=d->max.np;
      nvshft=d->max.nv;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Optimized shift: npshft = %d, nvshft = %d\n",npshft,nvshft);
  fflush(stdout);
#endif
      break;
    default:
      /* Invalid mode */
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Invalid 2nd argument %s(*,'type')\n",function);
      return(1);
      break;
  } /* end mode switch */

  /* Shift the data */
  shift2Ddata(d,npshft,nvshft);

  /* Update shift flag */
  if (d->shift) d->shift=FALSE;  /* if 'shifted' set to 'not shifted' */
  else d->shift=TRUE;            /* else set flag to 'shifted' */

  return(0);
}

int shift2Ddata(struct data *d,int npshft,int nvshft)
{
  fftw_complex *slice;
  int dim1,dim2,dim3,nr;
  int i,j,k,l;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"shift2Ddata"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  npshft = %d, nvshft = %d: ",npshft,nvshft);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Allocate memory for slice data */
  if ((slice = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Loop over receivers and slices */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      /* Fill slice with shifted data */
      for (k=0;k<nvshft;k++) {
        for (l=0;l<npshft;l++) {
          slice[(k+dim2-nvshft)*dim1+l+dim1-npshft][0]=d->data[i][j][k*dim1+l][0];
          slice[(k+dim2-nvshft)*dim1+l+dim1-npshft][1]=d->data[i][j][k*dim1+l][1];
        }
        for (l=npshft;l<dim1;l++) {
          slice[(k+dim2-nvshft)*dim1+l-npshft][0]=d->data[i][j][k*dim1+l][0];
          slice[(k+dim2-nvshft)*dim1+l-npshft][1]=d->data[i][j][k*dim1+l][1];
        }
      }
      for (k=nvshft;k<dim2;k++) {
        for (l=0;l<npshft;l++) {
          slice[(k-nvshft)*dim1+l+dim1-npshft][0]=d->data[i][j][k*dim1+l][0];
          slice[(k-nvshft)*dim1+l+dim1-npshft][1]=d->data[i][j][k*dim1+l][1];
        }
        for (l=npshft;l<dim1;l++) {
          slice[(k-nvshft)*dim1+l-npshft][0]=d->data[i][j][k*dim1+l][0];
          slice[(k-nvshft)*dim1+l-npshft][1]=d->data[i][j][k*dim1+l][1];
        }
      }
      /* Copy shifted data back to d->data */
      for (k=0;k<dim2*dim1;k++) {
        d->data[i][j][k][0]=slice[k][0];
        d->data[i][j][k][1]=slice[k][1];
      }
    }
  }

  fftw_free(slice);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"took %f secs\n",t2-t1);
  fflush(stdout);
#endif

  return(0);
}

void dimorder2D(struct data *d)
{
  /* Fill dim2order with the nv phase encode order */
  d->dim2order=phaseorder(d,d->nv,"pelist");

  /* Fill pssorder with the slice order */
  d->pssorder=sliceorder(d,d->ns,"pss");

  /* Set dimorder flag */
  d->dimorder=IM2D;
}

int phaseorder2D(struct data *d)
{
  /* There is a problem if data was acquired using a 'tablib' file.
     The file is not stored with data.
     The plan is to write table values to pelist to give a permanent record
     of the order with saved data */
  /* For the moment just accomodate for FSE acquired with standard ordering */
  fftw_complex *slice;
  double *pelist;
  int *tab,*petab;
  int dim1,dim2,ns,nr;
  int etl,kzero;
  int i,j,k,l,m;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"phaseorder2D"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; ns=d->ns; nr=d->nr;

  if ((tab = (int *)malloc(dim2*sizeof(int))) == NULL) nomem();
  if ((petab = (int *)malloc(dim2*sizeof(int))) == NULL) nomem();

  etl=1;
  if (!(strcmp(*sval("apptype",&d->p),"im2Dfse"))) { /* 2D FSE */
    /* Calculate PE table for kzero=1 */
    /* Data will end up in petab */
    etl=(int)*val("etl",&d->p);
    kzero=(int)*val("kzero",&d->p);
    for (i=0;i<dim2/2;i++)
      tab[i]= -(((i+1)*(dim2/(2*etl))-1)%(dim2/2))+i/etl;
    for (i=dim2/2;i<dim2;i++)
      tab[i]= ((i-dim2/2)*(dim2/(2*etl))%(dim2/2))+(i-dim2/2)/etl+1;

    /* Adjust for actual kzero filling petab as we go */
    for (i=0;i<dim2/etl;i++) {
      for (j=0;j<etl;j++) {
        k=i*etl;
        if (j < kzero-1) k+=etl-j-1;
        else k+=j-kzero+1;
        petab[i*etl+j]=tab[k];
      }
    }
  } else {
    /* Check pelist */
    if (nvals("pelist",&d->p) != dim2) {
      free(tab);
      free(petab);
#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Phase reordering not required (took %f secs)\n",t2-t1);
  fflush(stdout);
#endif
      return(1);
    }
    pelist=val("pelist",&d->p);
    for (i=0;i<dim2;i++) petab[i]=(int)pelist[i];
  }

#ifdef DEBUG
/*
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Calculated PE table for reordering:\n");
  for (i=0;i<dim2/etl;i++) {
    for (j=0;j<etl;j++)
      fprintf(stdout,"%6.1d",petab[i*etl+j]);
    fprintf(stdout,"\n");
  }
  fflush(stdout);
*/
#endif

  /* Fill tab with the required order */
  for (i=0;i<dim2;i++) tab[i]=i-dim2/2+1;

#ifdef DEBUG
/*
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Required PE table:\n");
  for (i=0;i<dim2/etl;i++) {
    for (j=0;j<etl;j++)
      fprintf(stdout,"%6.1d",tab[i*etl+j]);
    fprintf(stdout,"\n");
  }
  fflush(stdout);
*/
#endif

  /* Allocate memory for slice data */
  if ((slice = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Fill slice with the required order and copy back to d->data */
  for (i=0;i<nr;i++) {
    for (j=0;j<ns;j++) {
      /* Fill slice with the required order */
      for (k=0;k<dim2;k++) {
        for (l=0;l<dim2;l++) {
          if (petab[k] == tab[l]) {
            for (m=0;m<dim1;m++) {
              slice[l*dim1+m][0]=d->data[i][j][k*dim1+m][0];
              slice[l*dim1+m][1]=d->data[i][j][k*dim1+m][1];
            }
            break;
          }
        }
      }
      /* Copy new PE ordering back to d->data */
      for (k=0;k<dim2*dim1;k++) {
        d->data[i][j][k][0]=slice[k][0];
        d->data[i][j][k][1]=slice[k][1];
      }
    }
  }

  fftw_free(slice);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  if (!(strcmp(*sval("apptype",&d->p),"im2Dfse"))) { /* 2D FSE */
    fprintf(stdout,"  Calculated PE table for reordering:\n");
    for (i=0;i<dim2/etl;i++) {
      for (j=0;j<etl;j++)
        fprintf(stdout,"%6.1d",petab[i*etl+j]);
      fprintf(stdout,"\n");
    }
  }
/*
  fprintf(stdout,"  Required PE table:\n");
  for (i=0;i<dim2/etl;i++) {
    for (j=0;j<etl;j++)
      fprintf(stdout,"%6.1d",tab[i*etl+j]);
    fprintf(stdout,"\n");
  }
*/
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fflush(stdout);
#endif

  free(tab);
  free(petab);

  return(0);
}


void sliceorder2D(struct data *d)
{
  int index;
  int ns,nr;
  double *pss;
  fftw_complex ***slice;
  int i,j,k;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"sliceorder2D"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Index of pss in pars struct, return if pss is not found */
  index=parindex("pss",&d->p);
  if (index < 0) return;

  /* Number of slices and receivers */
  ns=d->ns; nr=d->nr;

  /* Take a copy of the starting slice order */
  if ((pss = (double *)malloc(ns*sizeof(double))) == NULL) nomem();
  for (i=0;i<ns;i++) pss[i]=d->p.d[index][i];

  /* Sort slice positions into ascending order */
  qsort(d->p.d[index],ns,sizeof(d->p.d[index][0]),doublecmp);

  /* Allocate memory for some data pointers to the slices */
  if ((slice = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();
  for (i=0;i<nr;i++) { /* loop over receiver blocks */
    if ((slice[i] = (fftw_complex **)fftw_malloc(ns*sizeof(fftw_complex *))) == NULL) nomem();
    for (j=0;j<ns;j++) /* loop over slices */
      if ((slice[i][j] = (fftw_complex *)fftw_malloc(sizeof(fftw_complex))) == NULL) nomem();
  }

  /* Set new pointers to slices in ascending order */
  for (i=0;i<ns;i++) {
    for (j=0;j<ns;j++) {
      if (d->p.d[index][i] == pss[j]) {
        for (k=0;k<nr;k++)
          slice[k][i]=d->data[k][j];
        break;
      }
    }
  }

  /* Reset original data pointers in ascending order */
  for (i=0;i<nr;i++)
    for (j=0;j<ns;j++)
      d->data[i][j]=slice[i][j];

  fftw_free(slice);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fprintf(stdout,"            pss\n");
  fprintf(stdout,"  original\tnew\n");
  for (i=0;i<ns;i++) 
    fprintf(stdout,"  %f\t%f\n",pss[i],d->p.d[index][i]);
  fflush(stdout);
#endif

  free(pss);
}

void weightdata2D(struct data *d,int mode)
{
  int i,j,k,l;
  int ix1,ix2;
  double f;
  int dim1,dim2,dim3,nr;
  double lb,lb1,gf,gf1,sb,sb1;
  double *weight;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"weightdata2D"); /* Set function name */
#endif

  switch(mode) {
    case READ: /* Read only using standard parameters */
      /* Return if ftproc[1]=0 (as VnmrJ) */
      if ((*val("ftproc",&d->p) == 0.0)) return;
      lb=*val("lb",&d->p); /* Lorentzian */
      gf=*val("gf",&d->p); /* Gaussian */
      sb=*val("sb",&d->p); /* Sinebell */
      lb1=0.0; gf1=0.0; sb1=0.0;
      break;
    case PHASE: /* Phase only using standard parameters */
      /* Return if ftproc[1]=0 (as VnmrJ) */
      if ((*val("ftproc",&d->p) == 0.0)) return;
      lb1=*val("lb1",&d->p); /* Lorentzian */
      gf1=*val("gf1",&d->p); /* Gaussian */
      sb1=*val("sb1",&d->p); /* Sinebell */
      lb=0.0; gf=0.0; sb=0.0;
      break;
     case REFREAD: /* Read only using reference parameters */
      lb=*val("reflb",&d->p); /* Lorentzian */
      gf=*val("refgf",&d->p); /* Gaussian */
      sb=*val("refsb",&d->p); /* Sinebell */
      lb1=0.0; gf1=0.0; sb1=0.0;
      break;
    case MK: /* Mask parameters */
      lb=*val("masklb",&d->p);   /* Lorentzian */
      lb1=*val("masklb1",&d->p); /* Lorentzian */
      gf=*val("maskgf",&d->p);   /* Gaussian */
      gf1=*val("maskgf1",&d->p); /* Gaussian */
      sb=*val("masksb",&d->p);   /* Sinebell */
      sb1=*val("masksb1",&d->p); /* Sinebell */
      break;
    case SM: /* Sensitivity map parameters */
      lb=*val("smaplb",&d->p);   /* Lorentzian */
      lb1=*val("smaplb1",&d->p); /* Lorentzian */
      gf=*val("smapgf",&d->p);   /* Gaussian */
      gf1=*val("smapgf1",&d->p); /* Gaussian */
      sb=*val("smapsb",&d->p);   /* Sinebell */
      sb1=*val("smapsb1",&d->p); /* Sinebell */
      break;
    default: /* Standard parameters */
      /* Return if ftproc[1]=0 (as VnmrJ) */
      if ((*val("ftproc",&d->p) == 0.0)) return;
      lb=*val("lb",&d->p);   /* Lorentzian */
      lb1=*val("lb1",&d->p); /* Lorentzian */
      gf=*val("gf",&d->p);   /* Gaussian */
      gf1=*val("gf1",&d->p); /* Gaussian */
      sb=*val("sb",&d->p);   /* Sinebell */
      sb1=*val("sb1",&d->p); /* Sinebell */
/*
      at=*val("at",&d->p);
*/
  } /* end mode switch */

  /* Return if no weighting is active */
  if (!lb && !lb1 && !gf && !gf1 && !sb && !sb1) return;

  /* Return if not FID data and not shifted */
  if ((d->datamode != FID) || !d->shift) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Calculate weighting */
  if ((weight = (double *)malloc(dim2*dim1*sizeof(double))) == NULL) nomem();
  for(k=0;k<dim2;k++) {
    for (l=0;l<dim1;l++) {
      ix1=k*dim1+l;
      weight[ix1]=1.0;
    }
  }

  /* Lorentzian line broadening */
  if (lb) { /* lb is active */
    for (l=0;l<dim1/2;l++) {
      /* Lorentzian broadening as fraction of FOV */
      f=exp(-l*M_PI*lb);
      /* Lorentzian broadening in # pixels */
      /* f=exp(-l*M_PI*lb/dim1); */
      /* Lorentzian broadening in Hz (as VnmrJ) */
      /* f=exp(-l*at*M_PI*lb/dim1); */
      for(k=0;k<dim2;k++) {
        ix1=k*dim1+l;
        ix2=k*dim1+dim1-1-l;
        weight[ix1] *=f;
        weight[ix2] *=f;
      }
    }
  }
  if (lb1) { /* lb1 is active */
    for(k=0;k<dim2/2;k++) {
      /* Lorentzian broadening as fraction of FOV */
      f=exp(-k*M_PI*lb1);
      for (l=0;l<dim1;l++) {
        ix1=k*dim1+l;
        ix2=(dim2-1-k)*dim1+l;
        weight[ix1] *=f;
        weight[ix2] *=f;
      }
    }
  }

  /* Gaussian line broadening */
  if (gf) { /* gf is active */
    for (l=0;l<dim1/2;l++) {
      /* Gaussian broadening as fraction of FOV */
      f=l*M_PI*gf;
      f=exp(-f*f);
      for(k=0;k<dim2;k++) {
        ix1=k*dim1+l;
        ix2=k*dim1+dim1-1-l;
        weight[ix1] *=f;
        weight[ix2] *=f;
      }
    }
  }
  if (gf1) { /* gf1 is active */
    for(k=0;k<dim2/2;k++) {
      /* Gaussian broadening as fraction of FOV */
      f=k*M_PI*gf1;
      f=exp(-f*f);
      for (l=0;l<dim1;l++) {
        ix1=k*dim1+l;
        ix2=(dim2-1-k)*dim1+l;
        weight[ix1] *=f;
        weight[ix2] *=f;
      }
    }
  }

  /* Sinebell line broadening */
  if (sb) { /* sb is active */
    for (l=0;l<dim1/2;l++) {
      /* Sinebell broadening as fraction of FOV */
      if (l*M_PI*sb < M_PI/2)
        f=cos(l*M_PI*sb);
      else 
        f=0.0;
      for(k=0;k<dim2;k++) {
        ix1=k*dim1+l;
        ix2=k*dim1+dim1-1-l;
        weight[ix1] *=f;
        weight[ix2] *=f;
      }
    }
  }
  if (sb1) { /* sb1 is active */
    for(k=0;k<dim2/2;k++) {
      /* Sinebell broadening as fraction of FOV */
      if (k*M_PI*sb1 < M_PI/2)
        f=cos(k*M_PI*sb1);
      else 
        f=0.0;
      for (l=0;l<dim1;l++) {
        ix1=k*dim1+l;
        ix2=(dim2-1-k)*dim1+l;
        weight[ix1] *=f;
        weight[ix2] *=f;
      }
    }
  }

  /* Weight the data */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for(k=0;k<dim2;k++) {
        for (l=0;l<dim1;l++) {
          ix1=k*dim1+l;
          d->data[i][j][ix1][0] *=weight[ix1];
          d->data[i][j][ix1][1] *=weight[ix1];
        }
      }
    }
  }

  free(weight);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  if (lb) fprintf(stdout,"  lb weighting = %f\n",lb);
  if (lb1) fprintf(stdout,"  lb1 weighting = %f\n",lb1);
  if (gf) fprintf(stdout,"  gf weighting = %f\n",gf);
  if (gf1) fprintf(stdout,"  gf1 weighting = %f\n",gf1);
  if (sb) fprintf(stdout,"  sb weighting = %f\n",sb);
  if (sb1) fprintf(stdout,"  sb1 weighting = %f\n",sb1);
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fflush(stdout);
#endif

}

void zerofill2D(struct data *d,int mode)
{
  int i,j,k,l;
  int ix1,ix2;
  int dim1,dim2,dim3,nr;
  int fn,fn1;
  int ndim1,ndim2;
  fftw_complex *data;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"zerofill2D"); /* Set function name */
#endif

  switch(mode) {
    case MK: /* Mask parameters */
      fn=*val("maskfn",&d->p);
      fn1=*val("maskfn1",&d->p);
      break;
    default: /* Internally set */
      fn=d->fn;
      fn1=d->fn1;
  } /* end mode switch */

  /* Check that either fn or fn1 are active */
  /* If so, make sure they are exactly divisible by 4 */
  if (!fn && !fn1) return;
  if (!fn) fn = d->np;
  else { fn /=4; fn *=4; }
  if (!fn1) fn1 = d->nv*2;
  else { fn1 /=4; fn1 *=4; }

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Initial data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Make sure fn and fn1 are exactly divisible by 4 */
/*
  fn /=4; fn *=4;
  fn1 /=4; fn1 *=4;
*/

  /* New data dimensions */
  ndim1=fn/2;
  ndim2=fn1/2;

  /* Make data the new size and initialise to zero */
  if ((data = (fftw_complex *)fftw_malloc(ndim2*ndim1*sizeof(fftw_complex))) == NULL) nomem();
  for(k=0;k<ndim2;k++) {
    for (l=0;l<ndim1;l++) {
      ix1=k*ndim1+l;
      data[ix1][0]=0.0;
      data[ix1][1]=0.0;
    }
  }

  /* Copy from d->data[i][j] to data, zero filling as we go */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {

      if (ndim2 < dim2) {
        for(k=0;k<ndim2/2;k++) {

          if (ndim1 < dim1) {
            for (l=0;l<ndim1/2;l++) {
              ix1=k*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1-ndim1/2;l<dim1;l++) {
              ix1=k*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          } else {
            for (l=0;l<dim1/2;l++) {
              ix1=k*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1/2;l<dim1;l++) {
              ix1=k*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          }

        }
        for(k=dim2-ndim2/2;k<dim2;k++) {

          if (ndim1 < dim1) {
            for (l=0;l<ndim1/2;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1-ndim1/2;l<dim1;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          } else {
            for (l=0;l<dim1/2;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1/2;l<dim1;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          }

        }
      } else {
        for(k=0;k<dim2/2;k++) {

          if (ndim1 < dim1) {
            for (l=0;l<ndim1/2;l++) {
              ix1=k*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1-ndim1/2;l<dim1;l++) {
              ix1=k*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          } else {
            for (l=0;l<dim1/2;l++) {
              ix1=k*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1/2;l<dim1;l++) {
              ix1=k*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          }

        }
        for(k=dim2/2;k<dim2;k++) {

          if (ndim1 < dim1) {
            for (l=0;l<ndim1/2;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1-ndim1/2;l<dim1;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          } else {
            for (l=0;l<dim1/2;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l;
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
            for (l=dim1/2;l<dim1;l++) {
              ix1=(k-dim2+ndim2)*ndim1+l-(dim1-ndim1);
              ix2=k*dim1+l;
              data[ix1][0]=d->data[i][j][ix2][0];
              data[ix1][1]=d->data[i][j][ix2][1];
            }
          }

        }
      }

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
    }
  }

  /* Update parameters to reflect new data size */
  d->np=fn;
  d->nv=ndim2;

  /* Set zerofill flag */
  d->zerofill=TRUE;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Zero filling to give %d x %d data: took %f secs\n",ndim1,ndim2,t2-t1);
  fflush(stdout);
#endif

}

void phaseramp2D(struct data *d,int mode)
{
  int dim1,dim2,dim3,nr;
  int i,j,k,l,ix;
  double re,im,M,theta,factor;
  double offset,fov,oversample;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"phaseramp2D"); /* Set function name */
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
      if (oversample<1) oversample=1;
      /* Set phase ramp factor to correct for offset */
      factor=2*M_PI*offset/(fov*oversample);
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
      /* Return if there is no phase ramp to apply */
      offset=*val("ppe",&d->p);
      if (offset == 0.0) return;
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif
      fov=*val("lpe",&d->p);
      /* Set phase ramp factor to correct for offset */
      factor=-2*M_PI*offset/fov;
      /* Apply phase ramp to generate frequency shift */
      for (i=0;i<nr;i++) {
        for (j=0;j<dim3;j++) {
          for (k=0;k<dim2/2;k++) {
            ix = k*dim1;
            for (l=0;l<dim1;l++) {
              re=d->data[i][j][ix+l][0];
              im=d->data[i][j][ix+l][1];
              M=sqrt(re*re+im*im);
              theta = atan2(im,re) + factor*(k);
              d->data[i][j][ix+l][0]=M*cos(theta);
              d->data[i][j][ix+l][1]=M*sin(theta);
            }
          }
          for (k=dim2/2;k<dim2;k++) {
            ix = k*dim1;
            for (l=0;l<dim1;l++) {
              re=d->data[i][j][ix+l][0];
              im=d->data[i][j][ix+l][1];
              M=sqrt(re*re+im*im);
              theta = atan2(im,re) + factor*(k-dim2);
              d->data[i][j][ix+l][0]=M*cos(theta);
              d->data[i][j][ix+l][1]=M*sin(theta);
            }
          }
        }
      }
#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Phase encode phase ramp (%d traces): took %f secs\n",dim1,t2-t1);
  fflush(stdout);
#endif
      break;
  }

}

void phasedata2D(struct data *d,int mode)
{
  int dim1,dim2,dim3,nr;
  double rp,lp,lp1;
  double re,im,M,theta;
  int i,j,k,l;
  int ix;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
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
  rp=*val("rp",&d->p);
  lp=*val("lp",&d->p);
  lp1=*val("lp1",&d->p);

  /* Return if no phasing is required */
  if (!rp && !lp && !lp1) return;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Phase */
  /* Set lp and lp1 so that they adjust the phase by the specified
     phase in degrees accross the image */
  rp *=DEG2RAD;
  lp *=DEG2RAD/dim1;
  lp1 *=DEG2RAD/dim2;
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for(k=0;k<dim2/2;k++) {
        for (l=0;l<dim1/2;l++) {
          ix=k*dim1+l;
          re=d->data[i][j][ix][0];
          im=d->data[i][j][ix][1];
          M=sqrt(re*re+im*im);
          theta = atan2(im,re) + rp + lp*l + lp1*k;
          d->data[i][j][ix][0]=M*cos(theta);
          d->data[i][j][ix][1]=M*sin(theta);
        }
        for (l=dim1/2;l<dim1;l++) {
          ix=k*dim1+l;
          re=d->data[i][j][ix][0];
          im=d->data[i][j][ix][1];
          M=sqrt(re*re+im*im);
          theta = atan2(im,re) + rp + lp*(l-dim1) + lp1*k;
          d->data[i][j][ix][0]=M*cos(theta);
          d->data[i][j][ix][1]=M*sin(theta);
        }
      }
      for(k=dim2/2;k<dim2;k++) {
        for (l=0;l<dim1/2;l++) {
          ix=k*dim1+l;
          re=d->data[i][j][ix][0];
          im=d->data[i][j][ix][1];
          M=sqrt(re*re+im*im);
          theta = atan2(im,re) + rp + lp*l + lp1*(k-dim2);
          d->data[i][j][ix][0]=M*cos(theta);
          d->data[i][j][ix][1]=M*sin(theta);
        }
        for (l=dim1/2;l<dim1;l++) {
          ix=k*dim1+l;
          re=d->data[i][j][ix][0];
          im=d->data[i][j][ix][1];
          M=sqrt(re*re+im*im);
          theta = atan2(im,re) + rp + lp*(l-dim1) + lp1*(k-dim2);
          d->data[i][j][ix][0]=M*cos(theta);
          d->data[i][j][ix][1]=M*sin(theta);
        }
      }
/*
      for(k=0;k<dim2;k++) {
        for (l=0;l<dim1;l++) {
          ix=k*dim1+l;
          re=d->data[i][j][ix][0];
          im=d->data[i][j][ix][1];
          M=sqrt(re*re+im*im);
          theta = atan2(im,re) + rp +lp*(l-dim1/2) +lp1*(k-dim2/2);
          d->data[i][j][ix][0]=M*cos(theta);
          d->data[i][j][ix][1]=M*sin(theta);
        }
      }
*/
    }
  }

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs:\n",t2-t1);
  fflush(stdout);
#endif

}

void getmax2D(struct data *d)
{
  int dim1,dim2,dim3,nr;
  double re,im,M;
  int i,j,k,l;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"getmax2D"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Set some defaults */
  zeromax(d);

  /* Now get maximum and coordinates */
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

void zoomdata2D(struct data *d,int startdim1, int widthdim1, int startdim2, int widthdim2)
{
  int dim1,dim3,nr;
  int i,j,k,l;
  int ix1,ix2;
  fftw_complex *data;

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
  char function[20];
  strcpy(function,"zoomdata2D"); /* Set function name */
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Data dimensions */
  dim1=d->np/2;
  dim3=d->endpos-d->startpos;
  nr=d->nr;

  /* Make data the new size and initialise to zero */
  if ((data = (fftw_complex *)fftw_malloc(widthdim2*widthdim1*sizeof(fftw_complex))) == NULL) nomem();

  /* Copy from d->data[i][j] to data */
  for (i=0;i<nr;i++) {
    for (j=0;j<dim3;j++) {
      for(k=0;k<widthdim2;k++) {
        ix1=k*widthdim1;
        ix2=(k+startdim2)*dim1+startdim1;
        for (l=0;l<widthdim1;l++) {
          data[ix1][0]=d->data[i][j][ix2][0];
          data[ix1][1]=d->data[i][j][ix2][1];
          ix1++; ix2++;
        }
      }
      /* free d->data[i][j], reallocate and copy data back in its place */
      fftw_free(d->data[i][j]);
      if ((d->data[i][j] = (fftw_complex *)fftw_malloc(widthdim2*widthdim1*sizeof(fftw_complex))) == NULL) nomem();
      ix1=0;
      for(k=0;k<widthdim2;k++) {
        for (l=0;l<widthdim1;l++) {
          d->data[i][j][ix1][0]=data[ix1][0];
          d->data[i][j][ix1][1]=data[ix1][1];
          ix1++;
        }
      }
    }
  }

  /* Update parameters to reflect new data size */
  d->np=widthdim1*2;
  d->nv=widthdim2;

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Zoom to give %d x %d data: took %f secs\n",widthdim1,widthdim2,t2-t1);
  fflush(stdout);
#endif

}

/* rawIO3D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* rawIO3D.c: 3D raw IO routines                                             */
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

#define SOURCEFILE "common3D/rawIO3D.c"

void wrawbin3D(struct data *d,int dataorder,int type,int precision)
{
  char outdir[MAXPATHLEN];
  char function[20];

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
#endif

  strcpy(function,"wrawbin3D"); /* Set function name */

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Writing raw 3D data\n");
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Check that type is valid */
  if (!validtype(type) && (type != CX)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Invalid 3nd argument %s(*,*,'type',*)\n",function);
    return;
  }

  /* Set up output directory (outdir) according to type */
  switch(d->datamode) {
    case FID: strcpy(outdir,"raw"); break;
    default:
      switch (type) {
        case MG: strcpy(outdir,"recon"); break;
        case PH: strcpy(outdir,"recon"); break;
        case RE: strcpy(outdir,"recon"); break;
        case IM: strcpy(outdir,"recon"); break;
        case CX: strcpy(outdir,"recon"); break;
        default: strcpy(outdir,""); break;
      }
      break;
  } /* end d->datamode switch */
  switch(type) {
    case MG:  if (d->datamode == FID) strcat(outdir,"MG"); break; /* Magnitude */
    case PH:  strcat(outdir,"PH");    break; /* Phase */
    case RE:  strcat(outdir,"RE");    break; /* Real */
    case IM:  strcat(outdir,"IM");    break; /* Imaginary */
    case CX:  strcat(outdir,"CX");    break; /* Complex */
    case MK:  strcat(outdir,"mask");  break; /* Mask */
    case RMK: strcat(outdir,"maskR"); break; /* Reverse mask of magnitude */
    case SM:  strcat(outdir,"smap");  break; /* Sensitivity maps */
    case GF:  strcat(outdir,"gmap");  break; /* Geometry factor */
    case RS:  strcat(outdir,"Rsnr");  break; /* Relative SNR */
  } /* end type switch */
  switch(dataorder) {
    case D12: strcat(outdir,"_D12");  break;
    case D3:  strcat(outdir,"_D3");   break;
  }
  strcat(outdir,".raw");
  /* Select output */
  switch(type) {
    case MG:
      if (d->nr>1) {
        gen3Draw(d,'c',outdir,dataorder,type,precision);
        gen3Draw(d,'i',outdir,dataorder,type,precision);
      } else gen3Draw(d,'s',outdir,dataorder,type,precision);
      break;
    case MK:  gen3Draw(d,'s',outdir,dataorder,type,precision); break;
    case RMK: gen3Draw(d,'c',outdir,dataorder,type,precision); break;
    case GF:  gen3Draw(d,'g',outdir,dataorder,type,precision); break;
    case RS:  gen3Draw(d,'r',outdir,dataorder,type,precision); break;
    default:
      if (d->nr>1) gen3Draw(d,'i',outdir,dataorder,type,precision);
      else gen3Draw(d,'s',outdir,dataorder,type,precision);
      break;
  } /* end type switch */

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fflush(stdout);
#endif
}

void gen3Draw(struct data *d,int output,char *outdir,int dataorder,int type,int precision)
{
  char basename[MAXPATHLEN],dirname[MAXPATHLEN],filename[MAXPATHLEN];
  int ne,ns;
  int i;
  int slab,image,echo;
  int volindex;
  char function[20];
  strcpy(function,"gen3Draw"); /* Set function name */

  /* This function checks the output type requested and sets the output
     filename accordingly. 3D raw data is output according to the specified
     mode using functions w3Draw() and wcomb3Draw().
     These functions output data either from individual receivers (w3Draw)
     or using a combination of data from all receivers (wcomb3Draw). */

  /* Check for VnmrJ recon to figure if we are offline */
  if (vnmrj_recon) strcpy(basename,vnmrj_path);
  else {
    for (i=0;i<=strlen(d->file)-9;i++)
      basename[i]=d->file[i];
    basename[i]=0;
  }

  /* Number of echoes */
  ne=(int)*val("ne",&d->p);
  if (ne < 1) ne=1; /* Set ne to 1 if 'ne' does not exist */

  /* Number of slices (slabs) */
  ns=nvals("pss",&d->p);
  if (ns < 1) ns=1; /* Set ns to 1 if 'ns' does not exist */

  /* Allow for compressed multi-echo loop and multiple slabs */
  volindex=d->vol;
  image=volindex/(ne*ns);
  slab=(volindex/ne)%ns;
  echo=volindex%ne;

  switch(output) {
    case 'i': /* Individual output */
      for (i=0;i<d->nr;i++) {
        sprintf(dirname,"%s%s%.3d",basename,outdir,i+1);
        createdir(dirname);
        sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dirname,slab+1,image+1,echo+1);
        switch(dataorder) {
          case D3: w3DrawD3(filename,d,i); break;
          default: w3Draw(filename,d,i,type,precision); break;
        }
      }
      break;
    case 'c': /* Combined output (Magnitude only) */
      sprintf(dirname,"%s%s",basename,outdir);
      createdir(dirname);
      sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dirname,slab+1,image+1,echo+1);
      wcomb3Draw(filename,d,type,precision);
      break;
    case 's': /* Single receiver */
      sprintf(dirname,"%s%s",basename,outdir);
      createdir(dirname);
      sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dirname,slab+1,image+1,echo+1);
      switch(dataorder) {
        case D3: w3DrawD3(filename,d,0); break;
        default: w3Draw(filename,d,0,type,precision); break;
      }
      break;
    case 'g': /* Geometry Factor output */
      sprintf(dirname,"%s%s",basename,outdir);
      createdir(dirname);
      sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dirname,slab+1,image+1,echo+1);
      w3Draw(filename,d,0,type,precision);
      break;
    case 'r': /* Relative SNR */
      sprintf(dirname,"%s%s",basename,outdir);
      createdir(dirname);
      sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dirname,slab+1,image+1,echo+1);
      w3Draw(filename,d,1,type,precision);
      break;
    default:
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Invalid 2nd argument %s(*,'mode',*,*,*,*)\n",function);
      break;
  } /* end output switch */
}

void w3Draw(char *filename,struct data *d,int receiver,int type,int precision)
{
  FILE *f_out;
  float *floatdata;
  double *doubledata;
  double re,im,M;
  int dim1,dim2,dim3;
  int i,j,k;
  int ix1;
  char function[20];
  strcpy(function,"w3Draw"); /* Set function name */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos;

  /* Open file for writing if its the first processing block of a volume */
  if (d->block == 0) { /* The first processing block of a volume */
    if ((f_out = fopen(filename, "w")) == NULL) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to write to %s\n",filename);
      return;
    }
  } else { /* Open file for appending */
    if ((f_out = fopen(filename, "a")) == NULL) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to append to %s\n",filename);
      return;
    }
  }

  /* For Sensitivity Map, Geometry Factor and Relative SNR we just output
     the magnitude. */
  switch(type) {
    case SM: type = MG; break; /* Sensitivity Map */
    case GF: type = MG; break; /* Geometry Factor */
    case RS: type = MG; break; /* Relative SNR */
  }

  /* Allocate memory and write data */
  switch(precision) {
    case FLT32: /* 32 bit float */
      if ((floatdata = (float *)malloc(dim2*dim1*sizeof(float))) == NULL) nomem();
      switch(type) {
        case MG: /* Magnitude */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                re=fabs(d->data[receiver][k][ix1][0]);
                im=fabs(d->data[receiver][k][ix1][1]);
                M=sqrt(re*re+im*im);
                floatdata[ix1] = (float)M;
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
        case PH: /* Phase */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                re=d->data[receiver][k][ix1][0];
                im=d->data[receiver][k][ix1][1];
                /* atan2 returns values (in radians) between +PI and -PI */
                M=atan2(im,re);
                floatdata[ix1] = (float)M;
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
        case RE: /* Real */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                floatdata[ix1] = (float)d->data[receiver][k][ix1][0];
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
        case IM: /* Imaginary */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                floatdata[ix1] = (float)d->data[receiver][k][ix1][1];
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
        case CX: /* Complex */
          if ((floatdata = (float *)realloc(floatdata,dim2*dim1*2*sizeof(float))) == NULL) nomem();
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                floatdata[2*ix1] =   (float)d->data[receiver][k][ix1][0];
                floatdata[2*ix1+1] = (float)d->data[receiver][k][ix1][1];
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2*2,f_out);
          }
          break;
        case MK: /* Mask */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                floatdata[ix1] = (float)d->mask[k][ix1];
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
        case RMK: /* Reverse Mask of Magnitude */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                if (!d->mask[k][ix1]) {
                  re=fabs(d->data[receiver][k][ix1][0]);
                  im=fabs(d->data[receiver][k][ix1][1]);
                  M=sqrt(re*re+im*im);
                  floatdata[ix1] = (float)M;
                } else {
                  floatdata[ix1] = 0.0;
                }
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
        default:
          break;
      } /* end type switch */
      fclose(f_out);
      free(floatdata);
      break;
    case DBL64: /* 64 bit double */
      if ((doubledata = (double *)malloc(dim2*dim1*sizeof(double))) == NULL) nomem();
      switch(type) {
        case MG: /* Magnitude */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                re=fabs(d->data[receiver][k][ix1][0]);
                im=fabs(d->data[receiver][k][ix1][1]);
                M=sqrt(re*re+im*im);
                doubledata[ix1] = M;
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
        case PH: /* Phase */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                re=d->data[receiver][k][ix1][0];
                im=d->data[receiver][k][ix1][1];
                /* atan2 returns values (in radians) between +PI and -PI */
                M=atan2(im,re);
                doubledata[ix1] = M;
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
        case RE: /* Real */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                doubledata[ix1] = d->data[receiver][k][ix1][0];
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
        case IM: /* Imaginary */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                doubledata[ix1] = d->data[receiver][k][ix1][1];
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
        case CX: /* Complex */
          if ((doubledata = (double *)realloc(doubledata,dim2*dim1*2*sizeof(double))) == NULL) nomem();
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                doubledata[2*ix1] =   d->data[receiver][k][ix1][0];
                doubledata[2*ix1+1] = d->data[receiver][k][ix1][1];
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2*2,f_out);
          }
          break;
        case MK: /* Mask */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                doubledata[ix1] = (double)d->mask[k][ix1];
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
        case RMK: /* Reverse Mask of Magnitude */
          for (k=0;k<dim3;k++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                if (!d->mask[k][ix1]) {
                  re=fabs(d->data[receiver][k][ix1][0]);
                  im=fabs(d->data[receiver][k][ix1][1]);
                  M=sqrt(re*re+im*im);
                  doubledata[ix1] = M;
                } else {
                  doubledata[ix1] = 0.0;
                }
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
        default:
          break;
      } /* end type switch */
      fclose(f_out);
      free(doubledata);
      break;
    default:
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Invalid 5th argument %s(*,*,*,*,'precision')\n",function);
      break;
  } /* end precision switch */
}

void wcomb3Draw(char *filename,struct data *d,int type,int precision)
{
  FILE *f_out;
  float *floatdata;
  double *doubledata;
  double re,im,M;
  int dim1,dim2,dim3,nr;
  int i,j,k,l;
  int ix1;
  char function[20];
  strcpy(function,"wcomb3Draw"); /* Set function name */

  /* Check that type is valid */
  if ((type != MG) && (type != RMK)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Invalid 3rd argument %s(*,*,'type',*)\n",function);
    return;
  }

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr;

  /* Open file for writing if its the first processing block of a volume */
  if (d->block == 0) { /* The first processing block of a volume */
    if ((f_out = fopen(filename, "w")) == NULL) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to write to %s\n",filename);
      return;
    }
  } else { /* Open file for appending */
    if ((f_out = fopen(filename, "a")) == NULL) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to append to %s\n",filename);
      return;
    }
  }

  /* Allocate memory and write data */
  switch(precision) {
    case FLT32: /* 32 bit float */
      if ((floatdata = (float *)malloc(dim2*dim1*sizeof(float))) == NULL) nomem();
      switch(type) {
        case MG: /* Magnitude */
          for(l=0;l<dim3;l++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                M=0.0;
                for (k=0;k<nr;k++) {
                  re=fabs(d->data[k][l][ix1][0]);
                  im=fabs(d->data[k][l][ix1][1]);
                  M+=(re*re+im*im);
                }
                floatdata[ix1] = (float)sqrt(M);
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
        case RMK: /* Reverse mask of magnitude */
          for(l=0;l<dim3;l++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                if (!d->mask[l][ix1]) {
                  M=0.0;
                  for (k=0;k<nr;k++) {
                    re=fabs(d->data[k][l][ix1][0]);
                    im=fabs(d->data[k][l][ix1][1]);
                    M+=(re*re+im*im);
                  }
                  floatdata[ix1] = (float)sqrt(M);
                } else {
                  floatdata[ix1] = 0.0;
                }
              }
            }
            fwrite(floatdata,sizeof(float),dim1*dim2,f_out);
          }
          break;
      } /* end type switch */
      fclose(f_out);
      free(floatdata);
      break;
    case DBL64: /* 64 bit double */
      if ((doubledata = (double *)malloc(dim2*dim1*sizeof(double))) == NULL) nomem();
      switch(type) {
        case MG: /* Magnitude */
          for(l=0;l<dim3;l++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                M=0.0;
                for (k=0;k<nr;k++) {
                  re=fabs(d->data[k][l][ix1][0]);
                  im=fabs(d->data[k][l][ix1][1]);
                  M+=(re*re+im*im);
                }
                doubledata[ix1] = sqrt(M);
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
        case RMK: /* Reverse mask of magnitude */
          for(l=0;l<dim3;l++) {
            for(i=0;i<dim2;i++) {
              for (j=0;j<dim1;j++) {
                ix1=i*dim1+j;
                if (!d->mask[l][ix1]) {
                  M=0.0;
                  for (k=0;k<nr;k++) {
                    re=fabs(d->data[k][l][ix1][0]);
                    im=fabs(d->data[k][l][ix1][1]);
                    M+=(re*re+im*im);
                  }
                  doubledata[ix1] = sqrt(M);
                } else {
                  doubledata[ix1] = 0.0;
                }
              }
            }
            fwrite(doubledata,sizeof(double),dim1*dim2,f_out);
          }
          break;
      } /* end type switch */
      fclose(f_out);
      free(doubledata);
      break;
    default:
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Invalid 4th argument %s(*,*,*,'precision')\n",function);
      break;
  } /* end precision switch */

}

void w3DrawD3(char *filename,struct data *d,int receiver)
{
  FILE *f_out;
  double *doubledata;
  int dim1,dim2,dim3;
  int i,j;
  char function[20];
  strcpy(function,"w3DrawD3"); /* Set function name */

  /* Data dimensions */
  dim1=d->np/2; dim2=d->endpos-d->startpos; dim3=d->nv2;

  /* Generate fdf header if its the first processing block of a volume */
  if (d->block == 0) { /* The first processing block of a volume */
    /* Open file for writing */
    if ((f_out = fopen(filename, "w")) == NULL) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to write to %s\n",filename);
      return;
    }
  } else { /* Open file for appending */
    if ((f_out = fopen(filename, "a")) == NULL) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to append to %s\n",filename);
      return;
    }
  }

  /* Allocate memory and write data */
  if ((doubledata = (double *)malloc(dim3*2*sizeof(double))) == NULL) nomem();
  for(i=0;i<dim2*dim1;i++) {
    for (j=0;j<dim3;j++) {
      doubledata[2*j]   = d->data[receiver][i][j][0];
      doubledata[2*j+1] = d->data[receiver][i][j][1];
    }
    fwrite(doubledata,sizeof(double),dim3*2,f_out);
  }

  fclose(f_out);
  free(doubledata);

}

void rrawbin3D(struct data *d,int dataorder,int type,int precision)
{
  char indir[MAXPATHLEN];
  char function[20];

#ifdef DEBUG
  struct timeval tp;
  double t1,t2;
  int rtn;
#endif

  strcpy(function,"rrawbin3D"); /* Set function name */

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Reading raw 3D data\n");
  rtn=gettimeofday(&tp, NULL);
  t1=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fflush(stdout);
#endif

  /* Check that type is valid */
  if (type != CX) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Invalid 3rd argument %s(*,*,'type',*)\n",function);
    return;
  }

  /* Set up input directory (indir) according to type */
  switch(d->datamode) {
    case FID: strcpy(indir,"raw"); break;
    default:  strcpy(indir,"recon"); break;
  } /* end d->datamode switch */
  strcat(indir,"CX");
  switch(dataorder) {
    case D12: setblock(d,d->nv);  strcat(indir,"_D12");  break;
    case D3:  setblock(d,d->nv2); strcat(indir,"_D3");   break;
  }
  strcat(indir,".raw");
  /* Select input */
  if (d->nr>1) get3Draw(d,'i',indir,dataorder);
  else get3Draw(d,'s',indir,dataorder);

#ifdef DEBUG
  rtn=gettimeofday(&tp, NULL);
  t2=(double)tp.tv_sec+(1.e-6)*tp.tv_usec;
  fprintf(stdout,"  Took %f secs\n",t2-t1);
  fflush(stdout);
#endif
}

void cleardim3data(struct data *d)
{
  int dim3;
  int i,j;
  if (d->datamode != NONE) {
    dim3=d->np/2*(d->endpos-d->startpos);
    for (i=0;i<d->nr;i++) {
      for (j=0;j<dim3;j++) fftw_free(d->data[i][j]);
      fftw_free(d->data[i]);
    }
    fftw_free(d->data);
    /* Set data mode */
    d->datamode=NONE;
  }
}

void get3Draw(struct data *d,int mode,char *indir,int dataorder)
{
  char basename[MAXPATHLEN],dirname[MAXPATHLEN],filename[MAXPATHLEN];
  int dim1=0,dim2=0,dim3=0,nr=0;
  int volindex;
  int ne,ns;
  int slab,image,echo;
  int i,j;
  char function[20];
  strcpy(function,"get3Draw"); /* Set function name */

  /* This function checks the output type requested and sets the input
     filename accordingly. 3D raw data is input according to the specified
     mode using functions w3Dfdf() and wcomb3Dfdf().
     These functions output data either from individual receivers (w3Dfdf)
     or using a combination of data from all receivers (wcomb3Dfdf). */

  volindex=d->vol;

  /* If a previous block has been zero filled matrix size will be incorrect */
  /* Refresh from procpar values */
/*  if (d->zerofill) setdatapars(d); */

  /* Set data dimensions */
  switch (dataorder) {
    case D12:
      if (d->zerofill) d->nv2=(int)*val("nv2",&d->p);
      dim1=d->np/2; dim2=d->endpos-d->startpos; dim3=d->nv2; nr=d->nr; 
      break;
    case D3:  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos; nr=d->nr; break;
  }

  /* Allocate memory according to nr */
  if ((d->data = (fftw_complex ***)fftw_malloc(nr*sizeof(fftw_complex **))) == NULL) nomem();

  /* Allocate memory according to block size */
  switch (dataorder) {
    case D12: /* in this case we read "D12" processing block into "D3" processing block */
      for (i=0;i<nr;i++) { /* loop over receivers and malloc for "D3" processing block */
        if ((d->data[i] = (fftw_complex **)fftw_malloc(dim1*dim2*sizeof(fftw_complex *))) == NULL) nomem();
        for (j=0;j<dim1*dim2;j++) /* loop over dim1 and 1st phase encode block */
          if ((d->data[i][j] = (fftw_complex *)fftw_malloc(dim3*sizeof(fftw_complex))) == NULL) nomem();
      }
      break;
    case D3: /* in this case we read "D3" processing block into "D12" processing block */
      for (i=0;i<nr;i++) { /* loop over receivers and malloc for "D12" processing block */
        if ((d->data[i] = (fftw_complex **)fftw_malloc(dim3*sizeof(fftw_complex *))) == NULL) nomem();
        for (j=0;j<dim3;j++) /* loop over dim3 */
          if ((d->data[i][j] = (fftw_complex *)fftw_malloc(dim2*dim1*sizeof(fftw_complex))) == NULL) nomem();
      }
      break;
  }

  /* Check for VnmrJ recon to figure if we are offline */
  if (vnmrj_recon) strcpy(basename,vnmrj_path);
  else {
    for (i=0;i<=strlen(d->file)-9;i++)
      basename[i]=d->file[i];
    basename[i]=0;
  }

  /* Number of echoes */
  ne=(int)*val("ne",&d->p);
  if (ne < 1) ne=1; /* Set ne to 1 if 'ne' does not exist */

  /* Number of slices (slabs) */
  ns=nvals("pss",&d->p);
  if (ns < 1) ns=1; /* Set ns to 1 if 'ns' does not exist */

  /* Allow for compressed multi-echo loop and multiple slabs */
  image=volindex/(ne*ns);
  slab=(volindex/ne)%ns;
  echo=volindex%ne;

  switch(mode) {
    case 'i': /* Individual input */
      for (i=0;i<d->nr;i++) {
        sprintf(dirname,"%s%s%.3d",basename,indir,i+1);
        sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dirname,slab+1,image+1,echo+1);
        switch(dataorder) {
          case D12: r3Draw(filename,d,i);   break;
          case D3:  r3DrawD3(filename,d,i); break;
        }
      }
      break;
    case 's': /* Single receiver */
      sprintf(dirname,"%s%s",basename,indir);
      sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dirname,slab+1,image+1,echo+1);
      switch(dataorder) {
        case D12: r3Draw(filename,d,0);   break;
        case D3:  r3DrawD3(filename,d,0); break;
      }
      break;
    default:
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Invalid 2nd argument %s(*,'mode',*,*,*,*)\n",function);
      break;
  } /* end mode switch */
}

void r3Draw(char *filename,struct data *d,int receiver)
{
  FILE *fp;
  double *doubledata;  /* Pointer for 64-bit double data */
  int dim1,dim2,dim3,np,nv;
  int startpos;
  long offset;
  int j,k;
  int ix;
  char function[20];

#ifdef DEBUG
  strcpy(function,"r3Draw"); /* Set function name */
  if (receiver == 0) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Reading 64-bit double data\n");
    fflush(stdout);
  }
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->endpos-d->startpos; dim3=d->nv2;
  np=d->np; nv=d->nv;

  /* Start position */
  startpos=d->startpos;

  /* Allocate memory */
  if ((doubledata = (double *)malloc(np*dim2*sizeof(double))) == NULL) nomem();

  /* Check to see if input file can be opened */
  if ((fp=fopen(filename,"r")) == NULL) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Unable to open %s\n\n",filename);
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }

  /* loop over dim3 (slices or 2nd phase encodes) */
  for (j=0;j<dim3;j++) {
    offset=(long)((j*np*nv+startpos*np)*sizeof(double));
    if (fseek(fp,offset,SEEK_SET)) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to set data offset %ld\n\n",offset);
      fprintf(stdout,"  Aborting ...\n\n");
      exit(1);
    }
    if (fread(doubledata,sizeof(double),np*dim2,fp) < np*dim2) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to read %d data ponts\n\n",np*dim2);
      fprintf(stdout,"  Aborting ...\n\n");
      exit(1);
    }
    ix=0;
    for (k=0;k<dim1*dim2;k++) {
      d->data[receiver][k][j][0]=doubledata[ix++];
      d->data[receiver][k][j][1]=doubledata[ix++];
    }
  }

  /* Free memory */
  free(doubledata);

}

void r3DrawD3(char *filename,struct data *d,int receiver)
{
  FILE *fp;
  double *doubledata;  /* Pointer for 64-bit double data */
  int dim1,dim2,dim3,nv2,twodim3;
  int startpos;
  long offset=0;
  int j,k;
  int ix;
  char function[20];

#ifdef DEBUG
  strcpy(function,"r3DrawD3"); /* Set function name */
  if (receiver == 0) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Reading 64-bit double data\n");
    fflush(stdout);
  }
#endif

  /* Data dimensions */
  dim1=d->np/2; dim2=d->nv; dim3=d->endpos-d->startpos;
  nv2=d->nv2; twodim3=dim3*2;

  /* Start position */
  startpos=d->startpos;

  /* Allocate memory */
  if ((doubledata = (double *)malloc(twodim3*sizeof(double))) == NULL) nomem();

  /* Check to see if input file can be opened */
  if ((fp=fopen(filename,"r")) == NULL) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Unable to open %s\n\n",filename);
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }

  for (k=0;k<dim1*dim2;k++) {
    offset=(long)((k*nv2+startpos)*2*sizeof(double));
    if (fseek(fp,offset,SEEK_SET)) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to set data offset %ld\n\n",offset);
      fprintf(stdout,"  Aborting ...\n\n");
      exit(1);
    }
    if (fread(doubledata,sizeof(double),twodim3,fp) < twodim3) {
      fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
      fprintf(stdout,"  Unable to read %d data ponts\n\n",twodim3);
      fprintf(stdout,"  Aborting ...\n\n");
      exit(1);
    }
    ix=0;
    for (j=0;j<dim3;j++) {
      d->data[receiver][j][k][0]=doubledata[ix++];
      d->data[receiver][j][k][1]=doubledata[ix++];
    }
  }

  /* Free memory */
  free(doubledata);

}

void delrawbin3D(struct data *d,int dataorder,int type)
{
  char basename[MAXPATHLEN],dir[MAXPATHLEN],dirname[MAXPATHLEN];
  int i;
  char function[20];

  strcpy(function,"delrawdir3D"); /* Set function name */

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  Deleting raw 3D data\n");
  fflush(stdout);
#endif

  /* Check that type is valid */
  if (type != CX) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Invalid 3rd argument %s(*,*,'type',*)\n",function);
    return;
  }

  /* Set up directory (dir) according to type */
  switch(d->datamode) {
    case FID: strcpy(dir,"raw"); break;
    default:  strcpy(dir,"recon"); break;
  } /* end d->datamode switch */
  strcat(dir,"CX");
  switch(dataorder) {
    case D12: strcat(dir,"_D12");  break;
    case D3:  strcat(dir,"_D3");   break;
  }
  strcat(dir,".raw");

  /* Check for VnmrJ recon to figure if we are offline */
  if (vnmrj_recon) strcpy(basename,vnmrj_path);
  else {
    for (i=0;i<=strlen(d->file)-9;i++)
      basename[i]=d->file[i];
    basename[i]=0;
  }

  /* Remove raw data */
  if (d->nr>1) {
    for (i=0;i<d->nr;i++) {
      sprintf(dirname,"%s%s%.3d",basename,dir,i+1);
      del3Draw(d,dirname);
      remove(dirname);
    }
  } else {
    sprintf(dirname,"%s%s",basename,dir);
    del3Draw(d,dirname);
    remove(dirname);
  }
}

void del3Draw(struct data *d,char *dir)
{
  char filename[MAXPATHLEN];
  int volindex;
  int ne,ns;
  int slab,image,echo;
  char function[20];
  strcpy(function,"del3Draw"); /* Set function name */

  volindex=d->vol;

  /* Number of echoes */
  ne=(int)*val("ne",&d->p);
  if (ne < 1) ne=1; /* Set ne to 1 if 'ne' does not exist */

  /* Number of slices (slabs) */
  ns=nvals("pss",&d->p);
  if (ns < 1) ns=1; /* Set ns to 1 if 'ns' does not exist */

  /* Allow for compressed multi-echo loop and multiple slabs */
  image=volindex/(ne*ns);
  slab=(volindex/ne)%ns;
  echo=volindex%ne;

  sprintf(filename,"%s/slab%.3dimage%.3decho%.3d.raw",dir,slab+1,image+1,echo+1);
  remove(filename);
}

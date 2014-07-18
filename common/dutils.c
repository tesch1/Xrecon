/* dutils.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* dutils.c: Data utilities                                                  */
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

#define SOURCEFILE "common/dutils.c"

void opendata(char *datafile,struct data *d)
{
  char function[20];
  strcpy(function,"opendata"); /* Set function name */

  /* Check status of data file */
  if (stat(datafile,&d->buf) == -1) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Unable to access %s\n",datafile);
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }
  /* Check to see if input file can be opened */
  if ((d->fp=fopen(datafile,"r")) == NULL) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Unable to open %s\n\n",datafile);
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }
  /* Exit if input file size is less than the VNMR file header size */
  if (d->buf.st_size < sizeof(d->fh)) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Incomplete file header (?)\n");
    fprintf(stdout,"  Aborting ...\n\n");
    exit(1);
  }
  /* Copy filename */
  if ((d->file = (char *)malloc((strlen(datafile)+1)*sizeof(char))) == NULL) nomem();
  strcpy(d->file,datafile);

  /* Get file header */
  getdfh(d);

  /* Initialize data defaults */
  initdata(d);

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  %s contains %d bytes\n",d->file,(int)d->buf.st_size);
  fflush(stdout);
#endif
}

int im1D(struct data *d)
{
  if ((IM1D <= d->seqmode) && (d->seqmode < IM2D)) return(1);
  else return(0);
}

int im2D(struct data *d)
{
  if ((IM2D <= d->seqmode) && (d->seqmode < IM3D)) return(1);
  else return(0);
}

int im3D(struct data *d)
{
  if ((IM3D <= d->seqmode) && (d->seqmode < IM4D)) return(1);
  else return(0);
}

void setreffile(struct file *datafile,struct data *d,char *refpar)
{
  char filename[MAXPATHLEN],refname[MAXPATHLEN];
  int i,filenamectr=0,refnamectr=0;

  /* Get the reference name from the reference parameter */
  strcpy(refname,*sval(refpar,&d->p));

  /* If recon is straight after acquisition ... */
  if (vnmrj_recon && (!strcmp(*sval("file",&d->p),"exp"))) {
      /* Use the reference name as is */
      setfile(datafile,refname);
  } else {
    if (vnmrj_recon) { /* Recon from within VnmrJ */
      /* Use the path as defined in the parameter "file" */
      strcpy(filename,*sval("file",&d->p));
    } else { /* Recon outside of VnmrJ */
      /* Use the path as defined in d->file */
      strcpy(filename,d->file);
      filename[strlen(filename)-8]=0;  /* NULL terminate to remove .fid/fid */
    }
    /* Now use the path as defined in filename */
    for (i=0;i<strlen(filename);i++)
      if (filename[i] == '/') filenamectr=i;
    if (filenamectr>0) filenamectr++;
    for (i=0;i<strlen(refname);i++)
      if (refname[i] == '/') refnamectr=i;
    if (refnamectr>0) refnamectr++;
    for (i=refnamectr;i<strlen(refname);i++) {
      filename[filenamectr]=refname[i];
      filenamectr++;
    }
    filename[filenamectr]=0; /* NULL terminate */
    setfile(datafile,filename);
  }

}

void setfile(struct file *datafile,char *datadir)
{
  datafile->nfiles=1;

  if ((datafile->fid = (char **)malloc(sizeof(char *))) == NULL) nomem();
  if ((datafile->procpar = (char **)malloc(sizeof(char *))) == NULL) nomem();

  /* Force '*.fid' and '*.fid/' input names to be '*.fid/fid'
     and guess 'procpar' location */
  if ((strcmp(&datadir[strlen(datadir)-4],".fid") == 0)) {
    if ((datafile->fid[0] = (char *)malloc((strlen(datadir)+5)*sizeof(char))) == NULL) nomem();
    if ((datafile->procpar[0] = (char *)malloc((strlen(datadir)+9)*sizeof(char))) == NULL) nomem();
    strcpy(datafile->fid[0],datadir);
    strcat(datafile->fid[0],"/fid");
    strcpy(datafile->procpar[0],datadir);
    strcat(datafile->procpar[0],"/procpar");
  }
  else if ((strcmp(&datadir[strlen(datadir)-5],".fid/") == 0)) {
    if ((datafile->fid[0] = (char *)malloc((strlen(datadir)+4)*sizeof(char))) == NULL) nomem();
    if ((datafile->procpar[0] = (char *)malloc((strlen(datadir)+8)*sizeof(char))) == NULL) nomem();
    strcpy(datafile->fid[0],datadir);
    strcat(datafile->fid[0],"fid");
    strcpy(datafile->procpar[0],datadir);
    strcat(datafile->procpar[0],"procpar");
  }
  else if ((strcmp(&datadir[strlen(datadir)-8],".fid/fid") == 0)) {
    if ((datafile->fid[0] = (char *)malloc((strlen(datadir)+1)*sizeof(char))) == NULL) nomem();
    if ((datafile->procpar[0] = (char *)malloc((strlen(datadir)+5)*sizeof(char))) == NULL) nomem();
    strcpy(datafile->fid[0],datadir);
    strcpy(datafile->fid[0],datadir);
    datafile->procpar[0][strlen(datafile->fid[0])-3]=0; /* NULL terminate after '/' */
    strcat(datafile->procpar[0],"procpar");
  }
}

void setfn(struct data *d1,struct data *d2,double multiplier)
{
  int fn;
  if (d2->fn==0) fn=d2->np;
  else fn=d2->fn;
  fn=round2int(multiplier*fn);
  setval(&d1->p,"fn",fn);
  d1->fn=fn;
}

void setfn1(struct data *d1,struct data *d2,double multiplier)
{
  int fn1;
  if (d2->fn1==0) fn1=d2->nv*2;
  else fn1=d2->fn1;
  fn1=round2int(multiplier*fn1);
  setval(&d1->p,"fn1",fn1);
  d1->fn1=fn1;
}

void copynblocks(struct data *d1,struct data *d2)
{
  copypar("nblocks",&d1->p,&d2->p);
}

void copymaskpars(struct data *d1,struct data *d2)
{
  copypar("imMK",&d1->p,&d2->p);
  copypar("masknoisefrac",&d1->p,&d2->p);
  copypar("maskstartslice",&d1->p,&d2->p);
  copypar("maskendslice",&d1->p,&d2->p);
  copypar("maskrcvrs",&d1->p,&d2->p);
  copypar("maskwlb",&d1->p,&d2->p);
  copypar("maskwgf",&d1->p,&d2->p);
  copypar("maskwsb",&d1->p,&d2->p);
  copypar("masklb",&d1->p,&d2->p);
  copypar("masklb1",&d1->p,&d2->p);
  copypar("maskgf",&d1->p,&d2->p);
  copypar("maskgf1",&d1->p,&d2->p);
  copypar("masksb",&d1->p,&d2->p);
  copypar("masksb1",&d1->p,&d2->p);
  copypar("maskfn",&d1->p,&d2->p);
  copypar("maskfn1",&d1->p,&d2->p);
  copypar("masklvlmode",&d1->p,&d2->p);
  copypar("masklvlmax",&d1->p,&d2->p);
  copypar("masklvlnoise",&d1->p,&d2->p);
  copypar("masklvlnoisefrac",&d1->p,&d2->p);
  copypar("dfill",&d1->p,&d2->p);
  copypar("dfilldim",&d1->p,&d2->p);
  copypar("dfillfrac",&d1->p,&d2->p);
  copypar("dfillloops",&d1->p,&d2->p);
  copypar("maskeqnoise",&d1->p,&d2->p);
}

void copysmappars(struct data *d1,struct data *d2)
{
  copypar("imSM",&d1->p,&d2->p);
  copypar("smapref",&d1->p,&d2->p);
  copypar("vcoilref",&d1->p,&d2->p);
  copypar("smapmask",&d1->p,&d2->p);
  copypar("smapnoisefrac",&d1->p,&d2->p);
  copypar("smapwlb",&d1->p,&d2->p);
  copypar("smapwgf",&d1->p,&d2->p);
  copypar("smapwsb",&d1->p,&d2->p);
  copypar("smaplb",&d1->p,&d2->p);
  copypar("smaplb1",&d1->p,&d2->p);
  copypar("smapgf",&d1->p,&d2->p);
  copypar("smapgf1",&d1->p,&d2->p);
  copypar("smapsb",&d1->p,&d2->p);
  copypar("smapsb1",&d1->p,&d2->p);
  copypar("smapeqnoise",&d1->p,&d2->p);
}

void copysensepars(struct data *d1,struct data *d2)
{
  copypar("accelread",&d1->p,&d2->p);
  copypar("accelphase",&d1->p,&d2->p);
  copypar("rmapread",&d1->p,&d2->p);
  copypar("rmapphase",&d1->p,&d2->p);
  copypar("noisematrix",&d1->p,&d2->p);
  copypar("printNM",&d1->p,&d2->p);
}

void copydbh(struct datablockhead *dbh1,struct datablockhead *dbh2)
{
   dbh2->scale   = dbh1->scale;   /* scaling factor */
   dbh2->status  = dbh1->status;  /* status of data in block */
   dbh2->index   = dbh1->index;   /* block index */
   dbh2->mode    = dbh1->mode;    /* mode of data in block */
   dbh2->ctcount = dbh1->ctcount; /* ct value for FID */
   dbh2->lpval   = dbh1->lpval;   /* F2 left phase in phasefile */
   dbh2->rpval   = dbh1->rpval;   /* F2 right phase in phasefile */
   dbh2->lvl     = dbh1->lvl;     /* F2 level drift correction */
   dbh2->tlt     = dbh1->tlt;     /* F2 tilt drift correction */
}

int spar(struct data *d,char *par,char *str)
{
  if (!(strcmp(*sval(par,&d->p),str))) return(TRUE);
  return(FALSE);
}

int checkequal(struct data *d1,struct data *d2,char *par,char *comment)
{
  int i;
  double *val1,*val2;

#ifdef DEBUG
  char function[20];
  strcpy(function,"checkequal"); /* Set function name */
#endif

  if (nvals(par,&d1->p) != nvals(par,&d2->p)) {
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  %s and %s have different number of %s\n",d2->procpar,d1->procpar,comment);
  fprintf(stdout,"  Aborting recon of %s ...\n",d1->file);
  fflush(stdout);
#endif
    return(FALSE);
  }
  val1=val(par,&d1->p);
  val2=val(par,&d2->p);
  for (i=0;i<nvals(par,&d1->p);i++) {
    if (val1[i] != val2[i]) {
#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  %s and %s have different %s\n",d2->procpar,d1->procpar,comment);
  fprintf(stdout,"  Aborting recon of %s ...\n",d1->file);
  fflush(stdout);
#endif
      return(FALSE);
    }
  }
  return(TRUE);
}

double getelem(struct data *d,char *par,int image)
{
  double *dbl;
  int cycle,n,id;
  dbl=val(par,&d->p);
  if (arraycheck(par,&d->a)) { /* parameter is arrayed */
    cycle=getcycle(par,&d->a);
    n=nvals(par,&d->p);
    id=(image/cycle)%n;
    return(dbl[id]);
  } else { /* just one element */
    return(dbl[0]);
  }
}

int wprocpar(struct data *d,char *filename)
{
  FILE *fp1,*fp2;
  struct stat buf;
  char chardata;
  char function[20];
  strcpy(function,"wprocpar"); /* Set function name */

  /* Return if a procpar file already exists */
  if (!stat(filename,&buf)) return(0);

  if ((fp1=fopen(d->procpar,"r")) == NULL) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Unable to open procpar file %s\n",d->procpar);
    return(1);
  }
  if ((fp2=fopen(filename,"w")) == NULL) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Unable to write procpar file %s\n",filename);
    return(1);
  }

  while ((fscanf(fp1,"%c",&chardata)) != EOF) fprintf(fp2,"%c",chardata);

  fclose(fp1);
  fclose(fp2);

  return(0);
}

void setdatapars(struct data *d)
{
  int i;
  char rcvrs[MAXRCVRS];

#ifdef DEBUG
  char function[20];
  strcpy(function,"setdatapars"); /* Set function name */
#endif

  /* Number of points and views */
  d->np=(int)*val("np",&d->p);
  d->nv=(int)*val("nv",&d->p);
  d->nv2=(int)*val("nv2",&d->p);

  if (d->nv == 0) d->nv=1;
  if (d->nv2 == 0) d->nv2=1;

  /* Number of receivers */
  strcpy(rcvrs,*sval("rcvrs",&d->p));
  d->nr=0;
  for (i=0;i<strlen(rcvrs);i++)
    if (rcvrs[i] == 'y') d->nr++;

  /* Number of pss values = slices */
  d->ns=nvals("pss",&d->p);

  /* There must be at least one block per volume */
  d->nblocks=(int)*val("nblocks",&d->p);
  if (!d->nblocks) d->nblocks++;

  /* Number of points and views for resizing data */
  d->fn=(int)*val("fn",&d->p);
  d->fn1=(int)*val("fn1",&d->p);
  d->fn2=(int)*val("fn2",&d->p);

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  fprintf(stdout,"  d->np  = %d\n",d->np);
  fprintf(stdout,"  d->nv  = %d\n",d->nv);
  fprintf(stdout,"  d->nv2 = %d\n",d->nv2);
  fprintf(stdout,"  d->nr  = %d\n",d->nr);
  fprintf(stdout,"  d->ns  = %d\n",d->ns);
  fprintf(stdout,"  d->fn  = %d\n",d->fn);
  fprintf(stdout,"  d->fn1 = %d\n",d->fn1);
  fprintf(stdout,"  d->fn2 = %d\n",d->fn2);
  fflush(stdout);
#endif
}

void initdata(struct data *d)
{
  d->nvols=0;
  d->datamode=NONE;
  d->shift=FALSE;
  d->noise.data=FALSE;
  d->noise.matrix=FALSE;
  d->maskdata=FALSE;
  d->zerofill=FALSE;
  d->dimorder=FALSE;
  zeromax(d);
  initnoise(d);
}

void closedata(struct data *d)
{
  fclose(d->fp);
  free(d->file);
  free(d->procpar);
}

void initnoise(struct data *d)
{
  /* Initialise pointers */
  if ((d->noise.M = (double *)malloc(d->nr*sizeof(double))) == NULL) nomem();
  if ((d->noise.M2 = (double *)malloc(d->nr*sizeof(double))) == NULL) nomem();
  if ((d->noise.Re = (double *)malloc(d->nr*sizeof(double))) == NULL) nomem();
  if ((d->noise.Im = (double *)malloc(d->nr*sizeof(double))) == NULL) nomem();
  /* Zero data */
  zeronoise(d);
}

void zeronoise(struct data *d)
{
  int i;
  d->noise.samples = 0;
  for (i=0;i<d->nr;i++) {
    d->noise.M[i] = 0.0;
    d->noise.M2[i] = 0.0;
    d->noise.Re[i] = 0.0;
    d->noise.Im[i] = 0.0;
  }
  d->noise.avM = 0.0;
  d->noise.avM2 = 0.0;
  d->noise.avRe = 0.0;
  d->noise.avIm = 0.0;
  /* Set data flag */
  d->noise.data=FALSE;
}

void zeromax(struct data *d)
{
  d->max.Mval=0.0;
  d->max.Rval=0.0;
  d->max.Ival=0.0;
  d->max.np=-1;
  d->max.nv=-1;
  d->max.nv2=-1;
  /* Set data flag */
  d->max.data=FALSE;
}

/* options.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* options.c: Input options                                                  */
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

#define SOURCEFILE "common/options.c"

void usage()
{
  fprintf(stdout,
          "Usage:\n"
          "    Xrecon -h <...>\n"
          "        print this help message\n"
          "\n"
          "    Xrecon -v <paths to VnmrJ acquisition>...\n"
          "        called from VnmrJ on unsaved acquisition data, outputs\n"
          "        reconstruction to (TODO...?)\n"
          "\n"
          "    Xrecon <paths to saved data>...\n"
          "        command line invocation on saved acquisition data, outputs\n"
          "        reconstruction to (.../data.img?)\n"
          "\n");
  exit(-1);
}

void getoptions(int argc,char *argv[])
{
  int vnmrj=FALSE;
  //int ival=0;
  //double val=0;
  //char str[MAXPATHLEN];
  int i;
  char function[20];
  strcpy(function,"getoptions"); /* Set function name */

  //*str=0;

  /* Read arguments of form -x aaa -y bbb and of form -xyz */
  while ((--argc > 0) && ((*++argv)[0] == '-')) {
    char *s;

    /* Interpret each non-null character after the '-' (ie argv[0]) as an
       option */
    for (s=argv[0]+1;*s;s++) {
      switch (*s) {
        case 'v': /* Set vnmrj flag */
          vnmrj=TRUE; 
          break;
        case 'h':
          usage();
          break;
#if 0 // these dont actually do anything
        case 'i': /* Read an integer */
          sscanf(*++argv,"%d",&ival);
          argc--;
          break;
        case 'r': /* Read a real */
          sscanf(*++argv,"%lf",&val);
          argc--;
          break;
        case 's': /* Read a string */
          sscanf(*++argv,"%s",str);
          argc--;
          break;
#endif
        default:
          fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
          fprintf(stdout,"  Illegal switch option\n");
          break;
      }
    }
  }

  /* Print the input values */
#ifdef DEBUG
  if (vnmrj) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  vnmrj flag set\n");
  }
#endif
#if 0
  if (ival > 0) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Example of command line integer input (= %d)\n",ival);
  }
  if (val > 0) {
    fprintf(stdout,"\noptions.c: getoptions()\n");
    fprintf(stdout,"  Example of command line real input (= %f)\n",val);
  }
  if (*str > 0) {
    fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
    fprintf(stdout,"  Example of command line string input (= %s)\n",str);
  }
#endif

  /* print usage if no reconstductions requested */
  if (0 == argc) usage();

  /* Allocate for input file name(s) in argument list */
  if ((f.fid = (char **)malloc(argc*sizeof(char *))) == NULL) nomem();
  if ((f.procpar = (char **)malloc(argc*sizeof(char *))) == NULL) nomem();
  for (i=0;i<argc;i++) {
    if ((f.fid[i] = (char *)malloc(MAXPATHLEN*sizeof(char))) == NULL) nomem();
    if ((f.procpar[i] = (char *)malloc(MAXPATHLEN*sizeof(char))) == NULL) nomem();
  }

  /* Get input file name(s) from argument list */
  for (i=0;i<argc;i++) sscanf(*argv++,"%s",f.fid[i]);

  /* Force '*.fid' and '*.fid/' input names to be '*.fid/fid'
     and guess 'procpar' location */
  vnmrj_recon=FALSE;
  if (vnmrj) {
    vnmrj_recon=TRUE;
    strcpy(vnmrj_path,f.fid[0]);
    strcat(vnmrj_path,"/");
    strcpy(f.procpar[0],f.fid[0]);
    strcat(f.fid[0],"/acqfil/fid");
    strcat(f.procpar[0],"/curpar");
  } else {
    for (i=0;i<argc;i++) {
      if ((strcmp(&f.fid[i][strlen(f.fid[i])-4],".fid") == 0)) {
        strcpy(f.procpar[i],f.fid[i]);
        strcat(f.fid[i],"/fid");
        strcat(f.procpar[i],"/procpar");
      }
      else if ((strcmp(&f.fid[i][strlen(f.fid[i])-5],".fid/") == 0)) {
        strcpy(f.procpar[i],f.fid[i]);
        strcat(f.fid[i],"fid");
        strcat(f.procpar[i],"procpar");
      }
      else if ((strcmp(&f.fid[i][strlen(f.fid[i])-8],".fid/fid") == 0)) {
        strcpy(f.procpar[i],f.fid[i]);
        f.procpar[i][strlen(f.fid[i])-3]=0; /* NULL terminate after '/' */
        strcat(f.procpar[i],"procpar");
      }
    }
  }

  /* Set the number of files */
  f.nfiles=argc;

#ifdef DEBUG
  fprintf(stdout,"\n%s: %s()\n",SOURCEFILE,function);
  for (i=0;i<f.nfiles;i++) 
    fprintf(stdout,"  file %2d: fid:%s\n           procpar:%s\n",i,f.fid[i],f.procpar[i]);
  fflush(stdout);
#endif

}

/* recon2D.c */
/*---------------------------------------------------------------------------*/
/*                                                                           */
/* recon2D.c: 2D recon                                                       */
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

#define SOURCEFILE "recon2D/recon2D.c"

void recon2D(struct data *d)
{
  char profile[6];

#ifdef DEBUG
  char function[20];
  strcpy(function,"recon2D"); /* Set function name */
#endif

  /* Check for profile flag */
  strcpy(profile,*sval("profile",&d->p));
  if (profile[0] == 'y') {
    profile2D(d);
    return;
  }

  /* Check recon flag for non-standard recons ... */

  /* Otherwise perform standard 2D recon */
  else default2D(d);

}

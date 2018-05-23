/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *  2018 Jacob Finkenrath
 *
 * This file is part of tmLQCD.
 *
 * tmLQCD is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * tmLQCD is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with tmLQCD.  If not, see <http://www.gnu.org/licenses/>.
 * 
 * 
 * **********
 * This routine is used to allocate memory to save the inital momentas.
 * It will be used if tuning of the integrator structure is activated
 * 
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
 # include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "sse.h"

su3adj * mo_tmp=NULL;

int init_moment_field_tmp(const int V) {
  int i = 0;

/*   posix_memalign(void **memptr, size_t alignment, size_t size) */
  if( (int*)(mo_tmp = (su3adj*)calloc(4*V+1, sizeof(su3adj))) == NULL){ 
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(1);
  }
  if((void*)(moment_tmp = (su3adj**)calloc(V,sizeof(su3adj*))) == NULL) {
    printf ("malloc errno : %d\n",errno); 
    errno = 0;
    return(2);
  }
#if ( defined SSE || defined SSE2 || defined SSE3)
  moment_tmp[0] = (su3adj*)(((unsigned long int)(mo_tmp)+ALIGN_BASE)&~ALIGN_BASE);
#else
  moment_tmp[0] = mo_tmp; 
#endif
  
  for(i = 1; i < V; i++){
    moment_tmp[i] = moment_tmp[i-1]+4;
  }

  return(0);
}

void free_moment_field_tmp() {

  free(mo_tmp);  
}

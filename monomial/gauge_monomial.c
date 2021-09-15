/***********************************************************************
 *
 * Copyright (C) 2008 Carsten Urbach
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
 ***********************************************************************/

#ifdef HAVE_CONFIG_H
# include<tmlqcd_config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#ifdef TM_USE_OMP 
# include <omp.h>
#endif
#include "global.h"
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "sse.h"
#include "start.h"
#include "gettime.h"
#include "get_rectangle_staples.h"
#include "gamma.h"
#include "get_staples.h"
#include "read_input.h"
#include "measure_gauge_action.h"
#include "measure_rectangles.h"
#include "monomial/monomial.h"
#include "hamiltonian_field.h"
#include "gauge_monomial.h"
#include "fatal_error.h"
#ifdef TM_USE_QUDA
#include "quda_interface.h"
#endif

/* this function compares the gauge derivative calculated by Quda and tmLQCD */
void compare_derivative(monomial *mnl, su3adj **arg1, su3adj **arg2 ){
  if (g_proc_id==0){
    for(int ix = 0; ix < VOLUME; ix++){
      for(int mu=0; mu<4; mu++){
        double *temp1=&(arg1[ix][mu].d1);
        double *temp2=&(arg2[ix][mu].d1);
        for(int j=0; j<8; ++j){
          double diff=(temp1[j]-temp2[j])/temp1[j];
          if (sqrt(diff*diff)>1e-10){
            fprintf(stderr, "The derivative computed by QUDA and tmLQCD does not match for monomial: c0=%e c1=%e g_beta=%e\n",mnl->c0,
                                                                           mnl->c1,
                                                                           mnl->beta);
            fprintf(stderr, "The difference is %e\n", diff);
            exit(1);
          }
        }
      }
    }
  }
}


/* this function calculates the derivative of the momenta: equation 13 of Gottlieb */
void gauge_derivative(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  double factor = -mnl->c0 * g_beta/3.0;
  if(mnl->use_rectangles) {
    mnl->forcefactor = 1.;
  }
  if( mnl->external_library == QUDA_LIB ){
#ifdef TM_USE_QUDA
    if (g_debug_level > 3) {
     memcpy(ddummy[0], hf->derivative[0], (4*VOLUMEPLUSRAND+1)*sizeof(su3adj));
    }
    compute_gauge_derivative_quda(mnl, hf);

    if (g_debug_level > 3){
     su3adj **given=hf->derivative;
     hf->derivative=ddummy;
     mnl->external_library=NO_EXT_LIB;
     gauge_derivative(id, hf);
     compare_derivative(mnl,given, ddummy);
     mnl->external_library=QUDA_LIB;
     hf->derivative=given;
    }
#else
    fatal_error("Attempted to use QUDA_LIB in gauge monomial but tmLQCD has been compiled without QUDA support!", __func__);
#endif
  } else {
    #ifdef TM_USE_OMP
    #pragma omp parallel
      {
    #endif
    
      su3 ALIGN v, w;
      int i, mu;
      su3 *z;
      su3adj *xm;
    
    #ifdef TM_USE_OMP
    #pragma omp for
    #endif
      for(i = 0; i < VOLUME; i++) { 
        for(mu=0;mu<4;mu++) {
          z=&hf->gaugefield[i][mu];
          xm=&hf->derivative[i][mu];
          get_staples(&v,i,mu, (const su3**) hf->gaugefield); 
          _su3_times_su3d(w,*z,v);
          _trace_lambda_mul_add_assign((*xm), factor, w);
          
          if(mnl->use_rectangles) {
            get_rectangle_staples(&v, i, mu);
            _su3_times_su3d(w, *z, v);
            _trace_lambda_mul_add_assign((*xm), factor*mnl->c1/mnl->c0, w);
          }
        }
      }
    
    #ifdef TM_USE_OMP
      } /* OpenMP closing brace */
    #endif
  }
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return;
}

/* this function calculates the derivative of the momenta: equation 13 of Gottlieb */
void gauge_EMderivative(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  double factor = -mnl->c0 * g_beta/3.0;
  if(mnl->use_rectangles) {
    mnl->forcefactor = 1.;
  }
  
#ifdef TM_USE_OMP
#pragma omp parallel
  {
#endif

  su3 ALIGN v, w;
  int i, mu;
  su3 *z;
  su3adj *xm;

#ifdef TM_USE_OMP
#pragma omp for
#endif
  for(i = 0; i < VOLUME; i++) { 
    // electric part
    z=&hf->gaugefield[i][0];
    xm=&hf->derivative[i][0];
    get_staples(&v, i, 0, (const su3**) hf->gaugefield); 
    _su3_times_su3d(w,*z,v);
    _trace_lambda_mul_add_assign((*xm), (1.+mnl->glambda)*factor, w);
    // lambda only acts on the plaquette, effectively changing c0 in the spatial and temporal parts, c1 remains untouched
    if(mnl->use_rectangles) {
	    get_rectangle_staples(&v, i, 0);
	    _su3_times_su3d(w, *z, v);
	    _trace_lambda_mul_add_assign((*xm), factor*mnl->c1/mnl->c0, w);
    }
    // magnetic part
    for(mu=1;mu<4;mu++) {
      z=&hf->gaugefield[i][mu];
      xm=&hf->derivative[i][mu];

      get_spacelike_staples(&v, i, mu, (const su3**) hf->gaugefield); 
      _su3_times_su3d(w, *z, v);
      _trace_lambda_mul_add_assign((*xm), (1.-mnl->glambda)*factor, w);

      get_timelike_staples(&v, i, mu, (const su3**) hf->gaugefield); 
      _su3_times_su3d(w, *z, v);
      _trace_lambda_mul_add_assign((*xm), (1.+mnl->glambda)*factor, w);
      if(mnl->use_rectangles) {
      	get_rectangle_staples(&v, i, mu);
      	_su3_times_su3d(w, *z, v);
      	_trace_lambda_mul_add_assign((*xm), factor*mnl->c1/mnl->c0, w);
      }
    }
  }

#ifdef TM_USE_OMP
  } /* OpenMP closing brace */
#endif
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return;
}

void gauge_heatbath(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  mnl->c0 = 1.;
  if(mnl->use_rectangles) mnl->c0 = 1. - 8.*mnl->c1;
  
  mnl->energy0 = g_beta*(mnl->c0 * measure_gauge_action( (const su3**) hf->gaugefield, mnl->glambda));
  if(mnl->use_rectangles) {
    mnl->energy0 += g_beta*(mnl->c1 * measure_rectangles( (const su3**) hf->gaugefield));
  }
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called gauge_heatbath for id %d energy %f\n", id, mnl->energy0);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return;
}

double gauge_acc(const int id, hamiltonian_field_t * const hf) {
  tm_stopwatch_push(&g_timers);
  monomial * mnl = &monomial_list[id];
  mnl->energy1 = g_beta*(mnl->c0 * measure_gauge_action( (const su3**) hf->gaugefield, mnl->glambda));
  if(mnl->use_rectangles) {
    mnl->energy1 += g_beta*(mnl->c1 * measure_rectangles( (const su3**) hf->gaugefield));
  }
  if(g_proc_id == 0) {
    if(g_debug_level > 3) {
      printf("called gauge_acc for id %d dH = %1.10e\n", 
	     id, mnl->energy0 - mnl->energy1);
    }
  }
  tm_stopwatch_pop(&g_timers, 0, 1, mnl->name, __func__);
  return(mnl->energy0 - mnl->energy1);
}

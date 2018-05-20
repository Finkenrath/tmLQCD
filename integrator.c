/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
 *               2012 Carsten Urbach, 2018 Jacob Finkenrath
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
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include "global.h"
#include "monomial/monomial.h"
#include "update_momenta.h"
#include "update_momenta_fg.h"
#include "update_gauge.h"
#include "hamiltonian_field.h"
#include "integrator.h"

integrator Integrator;

static const double fg_chi = 0.0138888888888889;
static const double omf4_rho = 0.2539785108410595;
static const double omf4_theta = -0.03230286765269967;
static const double omf4_vartheta = 0.08398315262876693;
static const double omf4_lamb = 0.6822365335719091;

/* Scheme C' of PHYSICAL REVIEW E 66, 026701 */
static const double opt4_lamb  = 0.375; // 3/8
static const double opt4_theta = 0.166666666666667; // 1/6
static const double opt4_chi   = 0.00520833333333333333; //1/192

/* Scheme G6 of PHYSICAL REVIEW E 66, 026701 */
static const double opt6_rho   = 0.1097059723948682;
static const double opt6_theta = 0.4140632267310831;
static const double opt6_nu    = 0.2693315848935301;
static const double opt6_lamb  = 1.131980348651556;
static const double opt6_chi   = -0.01324638643416052;
static const double opt6_mu    = 0.0008642161339706166;

/* Scheme G6' of PHYSICAL REVIEW E 66, 026701 */
static const double omf6_theta = 1.0798524263824310; // 1/2+(675+75*sqrt(6))^(1/3)/30 + 5/( 2 * (675+75*sqrt(6))^(1/3))
static const double omf6_nu    = 0.35995080879414365; // theta/3
static const double omf6_lamb  = -0.14371472730265417; // - 5*theta/3*(theta-1)
static const double omf6_xi    = -0.013965254224238837; //- 5*theta^2/144+theta/36-1/288
static const double omf6_chi   = -0.039247029382345630; // 1/144-theta/36*(theta/2+1)

/* OMF8FG Parameters */
static const double omf8_a[12]={
   0.0,
   0.6922517172738832,
   -0.3183450347119991,
   0.6766724088765565,
   -0.7207972470858706,
   0.3580316862350045,
   -0.3756270611751488,
   0.3580316862350045,
   -0.7207972470858706,
   0.6766724088765565,
   -0.3183450347119991,
   0.6922517172738832
};
static const double omf8_b[12]={
   0.1839699354244402,
   0.7084389757230299,
   0.1981440445033534,
   -0.06409380745116974,
   -0.6887429532761409,
   0.1622838050764871,
   0.1622838050764871,
   -0.6887429532761409,
   -0.06409380745116974,
   0.1981440445033534,
   0.7084389757230299,
   0.1839699354244402
};
static const double omf8_c[12]={
   0.0,
   0.03976209968238716,
   0.02245403440322733,
   0.0009405266232181224,
   -0.07336500519635302,
   0.02225664796363730,
   0.02225664796363730,
   -0.07336500519635302,
   0.0009405266232181224,
   0.02245403440322733,
   0.03976209968238716,
   0.0
};
/* OPT8FG Parameters */

static const double opt8_a[12]={
   0.41009674738801111928784693005080,
   -0.34123345756052780489101697378499,
   0.25644714021068150492361761631743,
   0.27765273975812438394100476242641,
   -0.56926266869753773902939657321159,
   0.46629949890124853576794423820194,
   0.46629949890124853576794423820194,
   -0.56926266869753773902939657321159,
   0.27765273975812438394100476242641,
   0.25644714021068150492361761631743,
   -0.34123345756052780489101697378499,
   0.41009674738801111928784693005080
};
static const double opt8_b[12]={
   0.0,
   0.0048249309817414952912695842664785,
   0.17492394861090375603419001374207,
   0.29304366370957066164364546204288,
   0.047448940168459770284238136482511,
   -0.0015299863411743974499219652320477,
   -0.037422994259002571606842462603791,
   -0.0015299863411743974499219652320477,
   0.047448940168459770284238136482511,
   0.29304366370957066164364546204288,
   0.17492394861090375603419001374207,
   0.0048249309817414952912695842664785,
};
static const double opt8_c[12]={
   0.0,
   0.00014743936907797528364717244760736,
   0.00023288450531932545357194967600155,
   0.0061648659635535962497705619884752,
   -0.012307516860831240716732016960034,
   -0.000073296648559126385387017161643798,
   0.015295860994523744731993293847001,
   -0.000073296648559126385387017161643798,
   -0.012307516860831240716732016960034,
   0.0061648659635535962497705619884752,
   0.00023288450531932545357194967600155,
   0.00014743936907797528364717244760736,
};
/* second order minimal norm integration scheme */
void integrate_2mn(const double tau, const int S, const int halfstep, const double tau2);
/* second order minimal norm integration scheme in velocity version */
void integrate_2mnp(const double tau, const int S, const int halfstep, const double tau2);
/* fourth order force gradient integration scheme */
void integrate_2mnfg(const double tau, const int S, const int halfstep, const double tau2);
/* optimal fourth order force gradient integration scheme in velocity version */
void integrate_opt4fg(const double tau, const int S, const int halfstep, const double tau2);
/* optimal sixth order force gradient integration scheme in velocity version */
void integrate_opt6fg(const double tau, const int S, const int halfstep, const double tau2);
/* OMF sixth order force gradient integration scheme */
void integrate_omf6fg(const double tau, const int S, const int halfstep, const double tau2);
/* optimal eigth order force gradient integration scheme in velocity version */
void integrate_opt8fg(const double tau, const int S, const int halfstep, const double tau2);
/* OMF eigth order force gradient integration scheme */
void integrate_omf8fg(const double tau, const int S, const int halfstep, const double tau2);
/* Leap frog integration scheme*/
void integrate_leap_frog(const double tau, const int S, const int halfstep, const double tau2);
/* fourth order OMF scheme */
void integrate_omf4(const double tau, const int S, const int halfstep,const double tau2 );
/* half step function */
void dohalfstep(const double tau, const int S);

/* function to initialise the integrator, to be called once at the beginning */

int init_integrator() {
  int i, ts;
  Integrator.hf.gaugefield = (su3 **) NULL;
  Integrator.hf.momenta = (su3adj **) NULL;
  Integrator.hf.derivative = (su3adj **) NULL;
  for(i = 0; i < 10; i++) {
    Integrator.no_mnls_per_ts[i] = 0;
  }
  if(Integrator.type[Integrator.no_timescales-1] == MN2p) {
    for(i = 0; i < Integrator.no_timescales; i++) {
      Integrator.type[i] = MN2p;
      Integrator.integrate[i] = &integrate_2mnp;
    }
  }
  else if (Integrator.type[Integrator.no_timescales-1] == OPT4FG) {
    for(i = 0; i < Integrator.no_timescales; i++) {
      Integrator.type[i] = OPT4FG;
      Integrator.integrate[i] = &integrate_opt4fg;
    }
  }
  else if (Integrator.type[Integrator.no_timescales-1] == OPT6FG) {
    for(i = 0; i < Integrator.no_timescales; i++) {
      Integrator.type[i] = OPT6FG;
      Integrator.integrate[i] = &integrate_opt6fg;
    }
  }
   else if (Integrator.type[Integrator.no_timescales-1] == OPT8FG) {
    for(i = 0; i < Integrator.no_timescales; i++) {
      Integrator.type[i] = OPT8FG;
      Integrator.integrate[i] = &integrate_opt8fg;
    }
  }
  else {
    for(i = 0; i < Integrator.no_timescales; i++) {
      if(Integrator.type[i] == MN2 || Integrator.type[i] == MN2p) {
         Integrator.integrate[i] = &integrate_2mn;
      }
      else if(Integrator.type[i] == LEAPFROG) {
         Integrator.integrate[i] = &integrate_leap_frog;
      }
      else if(Integrator.type[i] == OMF4) {
         Integrator.integrate[i] = &integrate_omf4;
      }
      else if(Integrator.type[i] == MN2FG || Integrator.type[i] == OPT4FG) {
         Integrator.integrate[i] = &integrate_2mnfg;
      }
      else if(Integrator.type[i] == OMF6FG || Integrator.type[i] == OPT6FG) {
         Integrator.integrate[i] = &integrate_omf6fg;
      }
      else if(Integrator.type[i] == OMF8FG || Integrator.type[i] == OPT8FG) {
         Integrator.integrate[i] = &integrate_omf8fg;
      }

    }
  }

  for(i = 0; i < no_monomials; i++) {
    ts = monomial_list[i].timescale;
    if(ts < Integrator.no_timescales && ts > -1) {
      Integrator.mnls_per_ts[ ts ][ Integrator.no_mnls_per_ts[ts] ] = monomial_list[i].id;
      Integrator.no_mnls_per_ts[ ts ]++;
    }
    else {
      if(g_proc_id == 0) {
	fprintf(stderr, "Warning: monomial %d is not on a valid timescale and will not be integrated\n", i);
      }
    }
  }
  for(i = 0; i < Integrator.no_timescales; i++) {
    if(Integrator.no_mnls_per_ts[ i ] < 1) {
      fprintf(stderr, "Error, no monomial on timescale %d!\nAborting...\n", i);
      exit(-1);
    }
  }
  return(0);
}

/* function to set the gauge and momenta fields for the integration */

void integrator_set_fields(hamiltonian_field_t * hf) {
  Integrator.hf.gaugefield = hf->gaugefield;
  Integrator.hf.momenta = hf->momenta;
  Integrator.hf.derivative = hf->derivative;
  Integrator.hf.update_gauge_copy = hf->update_gauge_copy;
  return;
}

/* and unsets again (to NULL pointer ) */

void integrator_unset_fields() {
  Integrator.hf.gaugefield = (su3 **) NULL;
  Integrator.hf.momenta = (su3adj **) NULL;
  Integrator.hf.derivative = (su3adj **) NULL;
  return;
}

void integrate_omf8fg(const double tau, const int S, const int halfstep, const double tau2){
  int i,j=0;
  int k;
  integrator * itgr = &Integrator;
  double eps,eps2;
  double oneminus2theta = 0.5*(1.-2.*omf6_theta);
  double oneminus2lambdanu = (1.-2.*(omf6_lamb+omf6_nu));
  
  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      for (k=1;k<11;k++)
      {
         update_gauge(omf8_a[k]*eps, &itgr->hf);
         update_momenta_fg(itgr->mnls_per_ts[0], eps , omf8_b[k], omf8_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
      }
      update_gauge(omf8_a[11]*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 2*omf8_b[11]*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    
   for (k=1;k<11;k++)
   {
      update_gauge(omf8_a[k]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , omf8_b[k], omf8_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
   }
   update_gauge(omf8_a[11]*eps, &itgr->hf);
   
   if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], omf8_b[11]*(eps+eps2), itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      for (k=1;k<11;k++)
      {
         itgr->integrate[S-1](omf8_a[k]*eps, S-1, 0, omf8_a[k+1]*eps);
         update_momenta_fg(itgr->mnls_per_ts[S], eps , omf8_b[k], omf8_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
      }
      itgr->integrate[S-1](omf8_a[11]*eps, S-1, 0, omf8_a[1]*eps);
      update_momenta(itgr->mnls_per_ts[S], 2*omf8_b[11]*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    
   for (k=1;k<11;k++)
   {
      itgr->integrate[S-1](omf8_a[k]*eps, S-1, 0, omf8_a[k+1]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , omf8_b[k], omf8_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
   }
   
   if(S == itgr->no_timescales-1) 
    {
      itgr->integrate[S-1](omf8_a[11]*eps, S-1, 1, omf8_a[1]*eps);
    }
    else 
       itgr->integrate[S-1](omf8_a[11]*eps, S-1, halfstep, omf8_a[1]*eps2);
       
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], omf8_b[11]*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  return;
}


void integrate_omf6fg(const double tau, const int S, const int halfstep, const double tau2){
   int i,j=0;
  integrator * itgr = &Integrator;
  double eps,eps2;
  double oneminus2theta = 0.5*(1.-2.*omf6_theta);
  double oneminus2lambdanu = (1.-2.*(omf6_lamb+omf6_nu));
  
  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      update_gauge(omf6_theta*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2theta*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , oneminus2lambdanu, omf6_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2theta*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(omf6_theta*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 2*omf6_nu*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    
   update_gauge(omf6_theta*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2theta*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0], eps , oneminus2lambdanu, omf6_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2theta*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(omf6_theta*eps, &itgr->hf);

   if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], omf6_nu*(eps+eps2), itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      itgr->integrate[S-1](omf6_theta*eps, S-1, 0, oneminus2theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, oneminus2theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , oneminus2lambdanu, omf6_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, omf6_theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](omf6_theta*eps, S-1, 0, omf6_theta*eps);
      update_momenta(itgr->mnls_per_ts[S], 2*omf6_nu*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    }
   itgr->integrate[S-1](omf6_theta*eps, S-1, 0, oneminus2theta*eps);
   update_momenta_fg(itgr->mnls_per_ts[S], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[S], &itgr->hf);
   itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, oneminus2theta*eps);
   update_momenta_fg(itgr->mnls_per_ts[S], eps , oneminus2lambdanu, omf6_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
   itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, omf6_theta*eps);
   update_momenta_fg(itgr->mnls_per_ts[S], eps , omf6_lamb, omf6_xi , itgr->no_mnls_per_ts[S], &itgr->hf);
   
    if(S == itgr->no_timescales-1) 
    {
      itgr->integrate[S-1](omf6_theta*eps, S-1, 1, omf6_theta*eps);
    }
    else 
       itgr->integrate[S-1](omf6_theta*eps, S-1, halfstep, omf6_theta*eps2);
       
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], omf6_nu*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  return;
}

void integrate_opt8fg(const double tau, const int S, const int halfstep, const double tau2){
  int i,k;
  integrator * itgr = &Integrator;
  double eps  = tau/((double)itgr->n_int[S]);
  double eps2 = tau2/((double)itgr->n_int[S]); // dummy stepsize
  
  if(S == 0) {
    update_gauge(opt8_a[0]*eps, &itgr->hf);
    for(i = 1; i < itgr->n_int[0]; i++) {
       for(k=1;k<11;k++){
          update_momenta_fg(itgr->mnls_per_ts[0],eps , opt8_b[k], opt8_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
          update_gauge(opt8_a[k]*eps, &itgr->hf);
       }
      update_momenta_fg(itgr->mnls_per_ts[0],eps , opt8_b[11], opt8_c[11] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(2*opt8_a[11]*eps, &itgr->hf);
    }
    
   for(k=1;k<11;k++){
      update_momenta_fg(itgr->mnls_per_ts[0],eps , opt8_b[k], opt8_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(opt8_a[k]*eps, &itgr->hf);
   }
   update_momenta_fg(itgr->mnls_per_ts[0],eps , opt8_b[11], opt8_c[11] , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(opt8_a[11]*eps, &itgr->hf);
   
  }
  else {
    for(i = 0; i < itgr->n_int[S]; i++) {
   
      integrate_opt8fg(opt8_a[0]*eps, S-1, halfstep, opt8_a[1]*eps);
      for(k=1;k<12;k++){
          update_momenta_fg(itgr->mnls_per_ts[S],eps , opt8_b[k], opt8_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
          integrate_opt8fg(opt8_a[k]*eps, S-1, halfstep, opt8_a[(k+1)%12]*eps);
       }
    }
  }
}

void integrate_opt6fg(const double tau, const int S, const int halfstep, const double tau2){
   int i;
  integrator * itgr = &Integrator;
  double eps  = tau/((double)itgr->n_int[S]);
  double eps2 = tau2/((double)itgr->n_int[S]); // dummy stepsize
  double oneminus2thetarho = 0.5*(1.-2.*(opt6_theta+opt6_rho));
  double oneminus2lambdanu = (1.-2.*(opt6_lamb+opt6_nu));
  
  if(S == 0) {
    update_gauge(opt6_rho*eps, &itgr->hf);
    for(i = 1; i < itgr->n_int[0]; i++) {
      update_momenta_fg(itgr->mnls_per_ts[0], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(opt6_theta*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], opt6_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2thetarho*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0],eps , oneminus2lambdanu, opt6_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2thetarho*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], opt6_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(opt6_theta*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(2*opt6_rho*eps, &itgr->hf);
    }
   update_momenta_fg(itgr->mnls_per_ts[0], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(opt6_theta*eps, &itgr->hf);
   update_momenta(itgr->mnls_per_ts[0], opt6_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2thetarho*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0],eps , oneminus2lambdanu, opt6_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2thetarho*eps, &itgr->hf);
   update_momenta(itgr->mnls_per_ts[0], opt6_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(opt6_theta*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(opt6_rho*eps, &itgr->hf);
  }
  else {
    for(i = 0; i < itgr->n_int[S]; i++) {
      integrate_opt6fg(opt6_rho*eps, S-1, halfstep, opt6_theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[S], &itgr->hf);
      integrate_opt6fg(opt6_theta*eps, S-1, halfstep, oneminus2thetarho*eps);
      
      update_momenta(itgr->mnls_per_ts[S], opt6_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      integrate_opt6fg(oneminus2thetarho*eps, S-1, halfstep, oneminus2thetarho*eps);
      update_momenta_fg(itgr->mnls_per_ts[S],eps , oneminus2lambdanu, opt6_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
      integrate_opt6fg(oneminus2thetarho*eps, S-1, halfstep, opt6_theta*eps);
      update_momenta(itgr->mnls_per_ts[S], opt6_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      
      integrate_opt6fg(opt6_theta*eps, S-1, halfstep, opt6_rho*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[S], &itgr->hf);
      integrate_opt6fg(opt6_rho*eps, S-1, halfstep, opt6_rho*eps);
    }
  }
}

void integrate_opt4fg(const double tau, const int S, const int halfstep, const double tau2){
   int i;
  integrator * itgr = &Integrator;
  double eps  = tau/((double)itgr->n_int[S]);
  double eps2 = tau2/((double)itgr->n_int[S]); // dummy stepsize
  double oneminus2theta = 0.5*(1.-2.*opt4_theta);
  double oneminus2lambda = (1.-2.*opt4_lamb);
  
  if(S == 0) {
    update_gauge(opt4_theta*eps, &itgr->hf);
    for(i = 1; i < itgr->n_int[0]; i++) {
      update_momenta(itgr->mnls_per_ts[0], opt4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2theta*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0],eps , oneminus2lambda, opt4_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2theta*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], opt4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(2*opt4_theta*eps, &itgr->hf);
    }
   update_momenta(itgr->mnls_per_ts[0], opt4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2theta*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0],eps , oneminus2lambda, opt4_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2theta*eps, &itgr->hf);
   update_momenta(itgr->mnls_per_ts[0], opt4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(opt4_theta*eps, &itgr->hf);
  }
  else {
    for(i = 0; i < itgr->n_int[S]; i++) {
      integrate_opt4fg(opt4_theta*eps, S-1, halfstep, oneminus2theta*eps);
      update_momenta(itgr->mnls_per_ts[S], opt4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      
      integrate_opt4fg(oneminus2theta*eps, S-1, halfstep, oneminus2theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S],eps , oneminus2lambda, opt4_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
      integrate_opt4fg(oneminus2theta*eps, S-1, halfstep, oneminus2theta*eps);
      
      update_momenta(itgr->mnls_per_ts[S], opt4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      integrate_opt4fg(opt4_theta*eps, S-1, halfstep, oneminus2theta*eps);
    }
  }
}

void integrate_omf4(const double tau, const int S, const int halfstep, const double tau2) {
  int i,j=0;
  integrator * itgr = &Integrator;
  double eps,eps2;

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      update_gauge(omf4_rho*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], omf4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(omf4_theta*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge((1-2.*(omf4_theta+omf4_rho))*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(omf4_theta*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], omf4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(omf4_rho*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 2*omf4_vartheta*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    update_gauge(omf4_rho*eps, &itgr->hf);
    update_momenta(itgr->mnls_per_ts[0], omf4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(omf4_theta*eps, &itgr->hf);
    update_momenta(itgr->mnls_per_ts[0], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge((1-2.*(omf4_theta+omf4_rho))*eps, &itgr->hf);
    update_momenta(itgr->mnls_per_ts[0], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(omf4_theta*eps, &itgr->hf);
    update_momenta(itgr->mnls_per_ts[0], omf4_lamb*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(omf4_rho*eps, &itgr->hf);
    if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], omf4_vartheta*(eps+eps2), itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      itgr->integrate[S-1](omf4_rho*eps, S-1, 0, omf4_theta*eps);
      update_momenta(itgr->mnls_per_ts[S], omf4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](omf4_theta*eps, S-1, 0, (1-2.*(omf4_theta+omf4_rho))*eps);
      update_momenta(itgr->mnls_per_ts[S], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1]((1-2.*(omf4_theta+omf4_rho))*eps, S-1, 0, omf4_theta*eps);
      update_momenta(itgr->mnls_per_ts[S], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](omf4_theta*eps, S-1, 0, omf4_rho*eps);
      update_momenta(itgr->mnls_per_ts[S], omf4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](omf4_rho*eps, S-1, 0, omf4_rho*eps);
      update_momenta(itgr->mnls_per_ts[S], 2*omf4_vartheta*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    itgr->integrate[S-1](omf4_rho*eps, S-1, 0, omf4_theta*eps);
    update_momenta(itgr->mnls_per_ts[S], omf4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    itgr->integrate[S-1](omf4_theta*eps, S-1, 0, (1-2.*(omf4_theta+omf4_rho))*eps);
    update_momenta(itgr->mnls_per_ts[S], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    itgr->integrate[S-1]((1-2.*(omf4_theta+omf4_rho))*eps, S-1, 0, omf4_theta*eps);
    update_momenta(itgr->mnls_per_ts[S], 0.5*(1-2.*(omf4_lamb+omf4_vartheta))*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    itgr->integrate[S-1](omf4_theta*eps, S-1, 0, omf4_rho*eps);
    update_momenta(itgr->mnls_per_ts[S], omf4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    if(S == itgr->no_timescales-1) {
      itgr->integrate[S-1](omf4_rho*eps, S-1, 1, omf4_rho*eps);
    }
    else itgr->integrate[S-1](omf4_rho*eps, S-1, halfstep, omf4_rho*eps2);
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], omf4_vartheta*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  return;
}

/* the following are only needed locally */

void integrate_2mn(const double tau, const int S, const int halfstep, const double tau2) {
  int i,j=0;
  integrator * itgr = &Integrator;
  double eps,eps2, oneminus2lambda = (1.-2.*itgr->lambda[S]);

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      update_gauge(0.5*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], oneminus2lambda*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(0.5*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 2.*itgr->lambda[0]*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    update_gauge(0.5*eps, &itgr->hf);
    update_momenta(itgr->mnls_per_ts[0], oneminus2lambda*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(0.5*eps, &itgr->hf);
    if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], itgr->lambda[0]*(eps+eps2), itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
      update_momenta(itgr->mnls_per_ts[S], oneminus2lambda*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
      update_momenta(itgr->mnls_per_ts[S], 2*itgr->lambda[S]*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
    update_momenta(itgr->mnls_per_ts[S], oneminus2lambda*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    if(S == itgr->no_timescales-1) {
      itgr->integrate[S-1](eps/2., S-1, 1, eps/2.);
    }
    else itgr->integrate[S-1](eps/2., S-1, halfstep, eps2/2.);
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], itgr->lambda[S]*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
}

void integrate_2mnp(const double tau, const int S, const int halfstep, const double tau2) {
  int i;
  integrator * itgr = &Integrator;
  double eps  = tau/((double)itgr->n_int[S]);
  double eps2 = tau2/((double)itgr->n_int[S]); // dummy stepsize
  double oneminus2lambda = (1.-2.*itgr->lambda[S]);
  
  if(S == 0) {
    update_gauge(itgr->lambda[0]*eps, &itgr->hf);
    for(i = 1; i < itgr->n_int[0]; i++) {
      update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2lambda*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(2*itgr->lambda[0]*eps, &itgr->hf);
    }
    update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(oneminus2lambda*eps, &itgr->hf);
    update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(itgr->lambda[0]*eps, &itgr->hf);
  }
  else {
    for(i = 0; i < itgr->n_int[S]; i++) {
      integrate_2mnp(itgr->lambda[S]*eps, S-1, halfstep,oneminus2lambda*eps);
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S], &itgr->hf);

      integrate_2mnp(oneminus2lambda*eps, S-1, halfstep,itgr->lambda[S]*eps);
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S], &itgr->hf);

      integrate_2mnp(itgr->lambda[S]*eps, S-1, halfstep,itgr->lambda[S]*eps);
    }
  }
}

/* For 2MNFG lamda MUST be equal to 1/6 */
void integrate_2mnfg(const double tau, const int S, const int halfstep, const double tau2) {
  int i,j=0;
  integrator * itgr = &Integrator;
  double eps,eps2, oneminus2lambda = (1.-2.*itgr->lambda[S]);

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      update_gauge(0.5*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0],eps, oneminus2lambda, fg_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(0.5*eps, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], 2.*itgr->lambda[0]*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    update_gauge(0.5*eps, &itgr->hf);
    update_momenta_fg(itgr->mnls_per_ts[0],eps, oneminus2lambda, fg_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(0.5*eps, &itgr->hf);
    if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], itgr->lambda[0]*(eps+eps2), itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
      update_momenta_fg(itgr->mnls_per_ts[S],eps, oneminus2lambda, fg_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
      update_momenta(itgr->mnls_per_ts[S], 2*itgr->lambda[S]*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
    update_momenta_fg(itgr->mnls_per_ts[S],eps, oneminus2lambda, fg_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
    if(S == itgr->no_timescales-1) {
      itgr->integrate[S-1](eps/2., S-1, 1, eps/2.);
    }
    else itgr->integrate[S-1](eps/2., S-1, halfstep, eps2/2.);
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], itgr->lambda[S]*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
}

void integrate_leap_frog(const double tau, const int S, const int halfstep,const double tau2) {
  int i;
  integrator * itgr = &Integrator;
  double eps, eps0, eps2;

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }

  eps = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  if(S == 0) {
    eps0 = tau/((double)itgr->n_int[0]); //what is the meaning of this variable ??
    for(i = 1; i < itgr->n_int[0]; i++) {
      update_gauge(eps0, &itgr->hf);
      update_momenta(itgr->mnls_per_ts[0], eps0, itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    update_gauge(eps0, &itgr->hf);
    if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], 0.5*(eps0+eps2), itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      itgr->integrate[S-1](eps, S-1, 0, eps);
      update_momenta(itgr->mnls_per_ts[S], eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    if(S == itgr->no_timescales-1) {
      itgr->integrate[S-1](eps, S-1, 1, eps);
    }
    else itgr->integrate[S-1](eps, S-1, halfstep, eps2);
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], 0.5*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
}


void dohalfstep(const double tau, const int S) {
  integrator * itgr = &Integrator;
  double eps = tau/((double)itgr->n_int[S]);
  for(int i = S; i > 0; i--) {
    if(itgr->type[i] == LEAPFROG) {
      update_momenta(itgr->mnls_per_ts[i], 0.5*eps, itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1]);
    }
    else if((itgr->type[i] == MN2) || (itgr->type[i] == MN2FG)){
      update_momenta(itgr->mnls_per_ts[i], itgr->lambda[i]*eps, itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])*2;
    }
    else if(itgr->type[i] == OMF4) {
      update_momenta(itgr->mnls_per_ts[i], omf4_vartheta*eps, itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/omf4_rho;
    }
    else if(itgr->type[i] == OMF6FG) {
      update_momenta(itgr->mnls_per_ts[i], omf6_nu*eps, itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/omf6_theta;
    }
    else if(itgr->type[i] == OMF8FG) {
      update_momenta(itgr->mnls_per_ts[i], omf8_b[0]*eps, itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/omf8_a[1];
    }

  }
  if(itgr->type[0] == LEAPFROG) {
    update_momenta(itgr->mnls_per_ts[0], 0.5*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if((itgr->type[0] == MN2)||(itgr->type[0] == MN2FG)) {
    update_momenta(itgr->mnls_per_ts[0], itgr->lambda[0]*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if(itgr->type[0] == OMF4) {
    update_momenta(itgr->mnls_per_ts[0], omf4_vartheta*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if (itgr->type[0] == OMF6FG) {
     update_momenta(itgr->mnls_per_ts[0], omf6_nu*(eps), itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if (itgr->type[0] == OMF8FG) {
     update_momenta(itgr->mnls_per_ts[0], omf8_b[0]*(eps), itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  return;
}

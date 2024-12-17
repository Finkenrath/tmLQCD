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
int int_type_tmp[10];

static const double lambda_2mnfg = 0.166666666666667;

static const double fg_chi = 0.0138888888888889;
static const double omf4_rho = 0.2539785108410595;
static const double omf4_theta = -0.03230286765269967;
static const double omf4_vartheta = 0.08398315262876693;
static const double omf4_lamb = 0.6822365335719091;

/* Scheme C' of PHYSICAL REVIEW E 66, 026701 */
static const double opt4_lamb  = 0.375; // 3/8
static const double opt4_theta = 0.166666666666667; // 1/6
static const double opt4_chi   = 0.00520833333333333333; //1/192

/* Scheme C' of PHYSICAL REVIEW E 66, 026701 
static const double opt4_lamb  = 0.2470939580390842; // 3/8
static const double opt4_theta = 0.08935804763220157; // 1/6
static const double opt4_chi   = 0.006938106540706989; //1/192*/

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

/*OMF4FG */
static const double omf4fg_a[4]={
   0.0,
   0.281398061166772,
   0.437203877666456,
   0.281398061166772
};
static const double omf4fg_b[4]={
   0.087893686016807,
   0.412106313983193,
   0.412106313983193,
   0.087893686016807
};
static const double omf4fg_c[4]={
   0.0,
   0.003061810122369770,
   0.003061810122369770,
   0.0
};
/*OMF4FG2 */
static const double omf4fg2_a[5]={
   0.000000000000000000,
   0.192112527742946404,
   0.307887472257053596,
   0.307887472257053596, 
   0.192112527742946404
};
static const double omf4fg2_b[5]={
   0.058518726134556207, 
   0.285216224068709112,
   0.312530099593469335, 
   0.285216224068709112,
   0.058518726134556207
};
static const double omf4fg2_c[5]={
   0.0004339598806816256,
	0.0,
   0.002427475259663050,
	0.0,
	0.0004339598806816256
};

/*OPT4FG2 - eq. 78*/
static const double opt4fg2_a[6]={
   0.06419108866816235,
   0.1919807940455741,
   0.24382811728626352,
   0.24382811728626352, 
   0.1919807940455741,
   0.06419108866816235
};
static const double opt4fg2_b[6]={
   0.0,
   0.1518179640276466,
   0.2158369476787619,
   0.26469017658718297, 
   0.2158369476787619,
   0.1518179640276466
};
static const double opt4fg2_c[6]={
   0.0,
   0.0,
   0.0009628905212024874,
   0.0,
   0.0009628905212024874,
   0.0
};

/*OMF6FG2 */
static const double omf6fg2_a[6]={
   0.000000000000000000,
   0.166738123347649092,
   0.380038934430259601,
   -0.093554115555817496, 
   0.380038934430259601,
   0.166738123347649092,
};
static const double omf6fg2_b[6]={
   0.048001369933520957,
   0.263395706993534984, 
   0.188602923072944073,
   0.188602923072944073,
   0.263395706993534984,
   0.048001369933520957 
};
static const double omf6fg2_c[6]={
  -0.001709693171449844,
  0.004668083730519805,
  0.000000000000000000,
  0.000000000000000000,
  0.004668083730519805,
  -0.001709693171449844 
};
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
/* optimal fourth order force gradient integration scheme in position version */
void integrate_opt4fg(const double tau, const int S, const int halfstep, const double tau2);
/* OPT4FG2 fourth order force gradient integration scheme in position version */
void integrate_opt4fg2(const double tau, const int S, const int halfstep, const double tau2);
/* optimal sixth order force gradient integration scheme in velocity version */
void integrate_opt6fg(const double tau, const int S, const int halfstep, const double tau2);
/* OMF sixth order force gradient integration scheme */
void integrate_omf6fg(const double tau, const int S, const int halfstep, const double tau2);
/* OMF - P5 integration scheme using input parameters (same stage like omf6fg)*/
void integrate_omfp5(const double tau, const int S, const int halfstep, const double tau2);
/* optimal eigth order force gradient integration scheme in velocity version */
void integrate_omf4fg(const double tau, const int S, const int halfstep, const double tau2);
/* optimal eigth order force gradient integration scheme in velocity version */
void integrate_omf4fg2(const double tau, const int S, const int halfstep, const double tau2);
/* optimal eigth order force gradient integration scheme in velocity version */
void integrate_omf6fg2(const double tau, const int S, const int halfstep, const double tau2);
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
  
  if ((Integrator.type[Integrator.no_timescales-1] == MN2p)||(Integrator.type[Integrator.no_timescales-1] == OPT4FG)||(Integrator.type[Integrator.no_timescales-1] == OPT6FG)||(Integrator.type[Integrator.no_timescales-1] == OPT8FG)||(Integrator.type[Integrator.no_timescales-1] == OPT4FG2))
  {
   for(i = 0; i < Integrator.no_timescales; i++) {
      if((Integrator.type[i] == MN2p)||(Integrator.type[i] == MN2)||(Integrator.type[i] == LEAPFROG)) {
         Integrator.integrate[i] = &integrate_2mnp;
      }
      else if((Integrator.type[i] == OPT4FG)||(Integrator.type[i] == OMF4)||(Integrator.type[i] == MN2FG)||(Integrator.type[i] == OMF4FG)||(Integrator.type[i] == OMF4FG2)) {
         Integrator.integrate[i] = &integrate_opt4fg;
      }
      else if((Integrator.type[i] == OPT4FG2)) {
         Integrator.integrate[i] = &integrate_opt4fg2;
      }
      else if((Integrator.type[i] == OPT6FG)||(Integrator.type[i] == OMF6FG)||(Integrator.type[i] == OMF6FG2)) {
         Integrator.integrate[i] = &integrate_opt6fg;
      }
      else if((Integrator.type[i] == OPT8FG)||(Integrator.type[i] == OMF8FG)) {
         Integrator.integrate[i] = &integrate_opt8fg;
      }
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
      else if(Integrator.type[i] == OMF4FG) {
         Integrator.integrate[i] = &integrate_omf4fg;
      }
      else if(Integrator.type[i] == OMF4FG2) {
         Integrator.integrate[i] = &integrate_omf4fg2;
      }
      else if(Integrator.type[i] == OMF6FG || Integrator.type[i] == OPT6FG) {
         Integrator.integrate[i] = &integrate_omf6fg;
      }
      else if(Integrator.type[i] == OMFP5) {
         Integrator.integrate[i] = &integrate_omfp5;
      }
      else if(Integrator.type[i] == OMF6FG2) {
         Integrator.integrate[i] = &integrate_omf6fg2;
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

/* function to initialise the integrator, to be called once at the beginning */

int init4tune_integrator(int level) {
  int i, ts;
 
  for (i= 0 ; i<Integrator.no_timescales; i++)
	  int_type_tmp[i]=Integrator.type[i];
  
  for(i = 0; i < 10; i++) {
    Integrator.no_mnls_per_ts[i] = 0;
  }
  
  
  if ((Integrator.type[Integrator.no_timescales-1] == MN2p)||(Integrator.type[Integrator.no_timescales-1] == OPT4FG)||(Integrator.type[Integrator.no_timescales-1] == OPT6FG)||(Integrator.type[Integrator.no_timescales-1] == OPT8FG))
  {
   for(i = 0; i < Integrator.no_timescales; i++)
      {
      if((Integrator.type[i] == MN2p) || Integrator.type[i] == MN2 || Integrator.type[i] == LEAPFROG) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_opt4fg;
				Integrator.type[i] = OPT4FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_2mnp;
				Integrator.type[i] = MN2p;
			}
      }
      if((Integrator.type[i] == OPT4FG) || (Integrator.type[i] == MN2FG) || (Integrator.type[i] == OMF4)) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_opt6fg;
				Integrator.type[i] = OPT6FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_opt4fg;
				Integrator.type[i] = OPT4FG;
			}
      }
      if((Integrator.type[i] == OPT6FG) || (Integrator.type[i] == OMF6FG) || (Integrator.type[i] == OMF6FG2)) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_opt8fg;
				Integrator.type[i] = OPT8FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_opt6fg;
				Integrator.type[i] = OPT6FG;
			}
      }
      else if (Integrator.type[Integrator.no_timescales-1] == OPT8FG) {
         if(g_proc_id == 0) {
            fprintf(stderr, "Warning: for OPT8FG scheme no higher order scheme available\n");
      }
      }
      }
  }
  else {
    for(i = 0; i < Integrator.no_timescales; i++) {
      if(Integrator.type[i] == MN2 || Integrator.type[i] == MN2p) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_omf4;
				Integrator.type[i] = MN2FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_2mn;
				Integrator.type[i] = MN2;
			}
      }
      else if(Integrator.type[i] == LEAPFROG) {
         if (i<level){
            Integrator.integrate[i] = &integrate_omf4;
				Integrator.type[i] = OMF4;
			}
         else
			{
            Integrator.integrate[i] = &integrate_leap_frog;
				Integrator.type[i] = LEAPFROG;
			}
      }
      else if(Integrator.type[i] == MN2FG || Integrator.type[i] == OPT4FG) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_omf8fg;
				Integrator.type[i] = OMF8FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_2mnfg;
				Integrator.type[i] = MN2FG;
			}
      }
      else if(Integrator.type[i] == OMF4) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_omf8fg;
				Integrator.type[i] = OMF8FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_omf4;
				Integrator.type[i] = OMF4;
			}
      }
      else if(Integrator.type[i] == OMF4FG) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_omf8fg;
				Integrator.type[i] = OMF8FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_omf4fg;
				Integrator.type[i] = OMF4FG;
			}
      }
      else if(Integrator.type[i] == OMF4FG2) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_omf8fg;
				Integrator.type[i] = OMF8FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_omf4fg2;
				Integrator.type[i] = OMF4FG2;
			}
      }
      else if(Integrator.type[i] == OMF6FG || Integrator.type[i] == OPT6FG) {
         if (i<level)
			{
            Integrator.integrate[i] = &integrate_omf8fg;
				Integrator.type[i] = OMF8FG;
			}
         else
			{
            Integrator.integrate[i] = &integrate_omf6fg;
				Integrator.type[i] = OMF6FG;
			}
      }
      else if(Integrator.type[i] == OMF8FG || Integrator.type[i] == OPT8FG) {
            if(g_proc_id == 0) {
            fprintf(stderr, "Warning: for OMF8FG or OPT8FG scheme no higher order scheme available\n");
         }
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

void reset4tune_integrator(void) {
	int i;
	for (i= 0 ; i<Integrator.no_timescales; i++)
	  Integrator.type[i]=int_type_tmp[i];
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

void integrate_omf4fg(const double tau, const int S, const int halfstep, const double tau2){
  int i,j=0;
  int k;
  integrator * itgr = &Integrator;
  double eps,eps2;
  
  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      for (k=1;k<3;k++)
      {
         update_gauge(omf4fg_a[k]*eps, &itgr->hf);
         update_momenta_fg(itgr->mnls_per_ts[0], eps , omf4fg_b[k], omf4fg_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
      }
      update_gauge(omf4fg_a[3]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], 2*eps, omf4fg_b[3], omf4fg_c[3], itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    
   for (k=1;k<3;k++)
   {
      update_gauge(omf4fg_a[k]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , omf4fg_b[k], omf4fg_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
   }
   update_gauge(omf4fg_a[3]*eps, &itgr->hf);
   
   if(halfstep != 1) {
      update_momenta_fg(itgr->mnls_per_ts[0], eps+eps2, omf4fg_b[3],  omf4fg_c[3], itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      for (k=1;k<3;k++)
      {
         itgr->integrate[S-1](omf4fg_a[k]*eps, S-1, 0, omf4fg_a[k+1]*eps);
         update_momenta_fg(itgr->mnls_per_ts[S], eps , omf4fg_b[k], omf4fg_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
      }
      itgr->integrate[S-1](omf4fg_a[3]*eps, S-1, 0, omf4fg_a[1]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], 2*eps, omf4fg_b[3], omf4fg_c[3], itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    
   for (k=1;k<3;k++)
   {
      itgr->integrate[S-1](omf4fg_a[k]*eps, S-1, 0, omf4fg_a[k+1]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , omf4fg_b[k], omf4fg_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
   }
   
   if(S == itgr->no_timescales-1) 
    {
      itgr->integrate[S-1](omf4fg_a[3]*eps, S-1, 1, omf4fg_a[1]*eps);
    }
    else 
       itgr->integrate[S-1](omf4fg_a[3]*eps, S-1, halfstep, omf4fg_a[1]*eps2);
       
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta_fg(itgr->mnls_per_ts[S], (eps+eps2), omf4fg_b[3],  omf4fg_c[3], itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  return;
}

void integrate_omf4fg2(const double tau, const int S, const int halfstep, const double tau2){
  int i,j=0;
  int k;
  integrator * itgr = &Integrator;
  double eps,eps2;
  
  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      for (k=1;k<4;k++)
      {
         update_gauge(omf4fg2_a[k]*eps, &itgr->hf);
         update_momenta_fg(itgr->mnls_per_ts[0], eps , omf4fg2_b[k], omf4fg2_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
      }
      update_gauge(omf4fg2_a[4]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps, omf4fg2_b[4], omf4fg2_c[4], itgr->no_mnls_per_ts[0], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps, omf4fg2_b[0], omf4fg2_c[0], itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    
   for (k=1;k<4;k++)
   {
      update_gauge(omf4fg2_a[k]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , omf4fg2_b[k], omf4fg2_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
   }
   update_gauge(omf4fg2_a[4]*eps, &itgr->hf);
   
   if(halfstep != 1) {
      update_momenta_fg(itgr->mnls_per_ts[0], eps, omf4fg2_b[4],  omf4fg2_c[4], itgr->no_mnls_per_ts[0], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps2, omf4fg2_b[0],  omf4fg2_c[0], itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      for (k=1;k<4;k++)
      {
         itgr->integrate[S-1](omf4fg2_a[k]*eps, S-1, 0, omf4fg2_a[k+1]*eps);
         update_momenta_fg(itgr->mnls_per_ts[S], eps , omf4fg2_b[k], omf4fg2_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
      }
      itgr->integrate[S-1](omf4fg2_a[4]*eps, S-1, 0, omf4fg2_a[1]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps, omf4fg2_b[4], omf4fg2_c[4], itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], eps, omf4fg2_b[0], omf4fg2_c[0], itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    
   for (k=1;k<4;k++)
   {
      itgr->integrate[S-1](omf4fg2_a[k]*eps, S-1, 0, omf4fg2_a[k+1]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , omf4fg2_b[k], omf4fg2_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
   }
   
   if(S == itgr->no_timescales-1) 
    {
      itgr->integrate[S-1](omf4fg2_a[4]*eps, S-1, 1, omf4fg2_a[1]*eps);
    }
    else 
       itgr->integrate[S-1](omf4fg2_a[4]*eps, S-1, halfstep, omf4fg2_a[1]*eps2);
       
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta_fg(itgr->mnls_per_ts[S], (eps), omf4fg2_b[4],  omf4fg2_c[4], itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], (eps2), omf4fg2_b[0],  omf4fg2_c[0], itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  return;
}

void integrate_opt4fg2(const double tau, const int S, const int halfstep, const double tau2){
  int i,k;
  integrator * itgr = &Integrator;
  double eps  = tau/((double)itgr->n_int[S]);
  double eps2 = tau2/((double)itgr->n_int[S]); // dummy stepsize
  
  if(S == 0) {
    update_gauge(opt4fg2_a[0]*eps, &itgr->hf);
    for(i = 1; i < itgr->n_int[0]; i++) {
       for(k=1;k<6;k++){
          
	 if (opt4fg2_c[k]==0.0)
		update_momenta(itgr->mnls_per_ts[0], opt4fg2_b[k]*eps, itgr->no_mnls_per_ts[0], &itgr->hf);
	 else
		update_momenta_fg(itgr->mnls_per_ts[0],eps , opt4fg2_b[k], opt4fg2_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
			 
          update_gauge(opt4fg2_a[k]*eps, &itgr->hf);
       }
      update_momenta_fg(itgr->mnls_per_ts[0],eps , opt4fg2_b[5], opt4fg2_c[5] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(2*opt4fg2_a[5]*eps, &itgr->hf);
    }
    
   for(k=1;k<6;k++){
		
	if (opt4fg2_c[k]==0.0)
		update_momenta(itgr->mnls_per_ts[S], opt4fg2_b[k]*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
	else
		update_momenta_fg(itgr->mnls_per_ts[S],eps , opt4fg2_b[k], opt4fg2_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
		
      update_gauge(opt4fg2_a[k]*eps, &itgr->hf);
   }
  }
  else {
    for(i = 0; i < itgr->n_int[S]; i++) {
   
      itgr->integrate[S-1](opt4fg2_a[0]*eps, S-1, halfstep, opt4fg2_a[1]*eps);
      for(k=1;k<6;k++){
			
	 if (opt4fg2_c[k]==0.0)
		update_momenta(itgr->mnls_per_ts[S], opt4fg2_b[k]*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
	 else
		update_momenta_fg(itgr->mnls_per_ts[S],eps , opt4fg2_b[k], opt4fg2_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
         itgr->integrate[S-1](opt4fg2_a[k]*eps, S-1, halfstep, opt4fg2_a[(k+1)%6]*eps);

       }
    }
  }
}

void integrate_omf6fg2(const double tau, const int S, const int halfstep, const double tau2){
  int i,j=0;
  int k;
  integrator * itgr = &Integrator;
  double eps,eps2;
  
  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      for (k=1;k<5;k++)
      {
         update_gauge(omf6fg2_a[k]*eps, &itgr->hf);
         update_momenta_fg(itgr->mnls_per_ts[0], eps , omf6fg2_b[k], omf6fg2_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
      }
      update_gauge(omf6fg2_a[5]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps, omf6fg2_b[5], omf6fg2_c[5], itgr->no_mnls_per_ts[0], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps, omf6fg2_b[5], omf6fg2_c[5], itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    
   for (k=1;k<5;k++)
   {
      update_gauge(omf6fg2_a[k]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , omf6fg2_b[k], omf6fg2_c[k] , itgr->no_mnls_per_ts[0], &itgr->hf);
   }
   update_gauge(omf6fg2_a[5]*eps, &itgr->hf);
   
   if(halfstep != 1) {
      update_momenta_fg(itgr->mnls_per_ts[0], eps, omf6fg2_b[5],  omf6fg2_c[5], itgr->no_mnls_per_ts[0], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps2, omf6fg2_b[5],  omf6fg2_c[5], itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      for (k=1;k<5;k++)
      {
         itgr->integrate[S-1](omf6fg2_a[k]*eps, S-1, 0, omf6fg2_a[k+1]*eps);
         update_momenta_fg(itgr->mnls_per_ts[S], eps , omf6fg2_b[k], omf6fg2_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
      }
      itgr->integrate[S-1](omf6fg2_a[5]*eps, S-1, 0, omf6fg2_a[1]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps, omf6fg2_b[5], omf6fg2_c[5], itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], eps, omf6fg2_b[5], omf6fg2_c[5], itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    
   for (k=1;k<5;k++)
   {
      itgr->integrate[S-1](omf6fg2_a[k]*eps, S-1, 0, omf6fg2_a[k+1]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , omf6fg2_b[k], omf6fg2_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
   }
   
   if(S == itgr->no_timescales-1) 
    {
      itgr->integrate[S-1](omf6fg2_a[5]*eps, S-1, 1, omf6fg2_a[1]*eps);
    }
    else 
       itgr->integrate[S-1](omf6fg2_a[5]*eps, S-1, halfstep, omf6fg2_a[1]*eps2);
       
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta_fg(itgr->mnls_per_ts[S], (eps), omf6fg2_b[5],  omf6fg2_c[5], itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], (eps2), omf6fg2_b[5],  omf6fg2_c[5], itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  return;
}

void integrate_omf8fg(const double tau, const int S, const int halfstep, const double tau2){
  int i,j=0;
  int k;
  integrator * itgr = &Integrator;
  double eps,eps2;
  
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

// copied from omf6fg -- not sure if it works
void integrate_omfp5(const double tau, const int S, const int halfstep, const double tau2){
   int i,j=0;
  integrator * itgr = &Integrator;
  double eps,eps2;
  double oneminus2theta = 0.5*(1.-2.*itgr->omfp5_theta[S]);
  double oneminus2lambdanu = (1.-2.*(itgr->omfp5_lam[S]+itgr->omfp5_nu[S]));
  
  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  eps  = tau/((double)itgr->n_int[S]);
  eps2 = tau2/((double)itgr->n_int[S]);
  
  if(S == 0) {

    for(j = 1; j < itgr->n_int[0]; j++) {
      update_gauge(itgr->omfp5_theta[0]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_lam[0], itgr->omfp5_xi[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2theta*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , oneminus2lambdanu, itgr->omfp5_chi[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(oneminus2theta*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_lam[0], itgr->omfp5_xi[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_gauge(itgr->omfp5_theta[0]*eps, &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_nu[0], itgr->omfp5_mu[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_nu[0], itgr->omfp5_mu[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    
   update_gauge(itgr->omfp5_theta[0]*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_lam[0], itgr->omfp5_xi[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2theta*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0], eps , oneminus2lambdanu, itgr->omfp5_chi[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(oneminus2theta*eps, &itgr->hf);
   update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_lam[0], itgr->omfp5_xi[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
   update_gauge(itgr->omfp5_theta[0]*eps, &itgr->hf);

   if(halfstep != 1) {
      update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_nu[0], itgr->omfp5_mu[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[0], eps2 , itgr->omfp5_nu[0], itgr->omfp5_mu[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      itgr->integrate[S-1](itgr->omfp5_theta[S]*eps, S-1, 0, oneminus2theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , itgr->omfp5_lam[S], itgr->omfp5_xi[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, oneminus2theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , oneminus2lambdanu, itgr->omfp5_chi[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, itgr->omfp5_theta[S]*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , itgr->omfp5_lam[S], itgr->omfp5_xi[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](itgr->omfp5_theta[S]*eps, S-1, 0, itgr->omfp5_theta[S]*eps);
      //update_momenta(itgr->mnls_per_ts[S], 2*omf6_nu*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , itgr->omfp5_nu[S], itgr->omfp5_mu[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , itgr->omfp5_nu[S], itgr->omfp5_mu[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
    }
   itgr->integrate[S-1](itgr->omfp5_theta[S]*eps, S-1, 0, oneminus2theta*eps);
   update_momenta_fg(itgr->mnls_per_ts[S], eps , itgr->omfp5_lam[S], itgr->omfp5_xi[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
   itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, oneminus2theta*eps);
   update_momenta_fg(itgr->mnls_per_ts[S], eps , oneminus2lambdanu, itgr->omfp5_chi[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
   itgr->integrate[S-1](oneminus2theta*eps, S-1, 0, itgr->omfp5_theta[S]*eps);
   update_momenta_fg(itgr->mnls_per_ts[S], eps , itgr->omfp5_lam[S], itgr->omfp5_xi[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
   
    if(S == itgr->no_timescales-1) 
    {
      itgr->integrate[S-1](itgr->omfp5_theta[S]*eps, S-1, 1, itgr->omfp5_theta[S]*eps);
    }
    else 
       itgr->integrate[S-1](itgr->omfp5_theta[S]*eps, S-1, halfstep, itgr->omfp5_theta[S]*eps2);
       
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      //update_momenta(itgr->mnls_per_ts[S], omf6_nu*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , itgr->omfp5_nu[S], itgr->omfp5_mu[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[S], eps2 , itgr->omfp5_nu[S], itgr->omfp5_mu[S] , itgr->no_mnls_per_ts[S], &itgr->hf);
    }
  }

  if(S == itgr->no_timescales-1) {
    dohalfstep(tau, S);
  }
  return;
}
// copy end





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
   
      itgr->integrate[S-1](opt8_a[0]*eps, S-1, halfstep, opt8_a[1]*eps);
      for(k=1;k<12;k++){
          update_momenta_fg(itgr->mnls_per_ts[S],eps , opt8_b[k], opt8_c[k] , itgr->no_mnls_per_ts[S], &itgr->hf);
          itgr->integrate[S-1](opt8_a[k]*eps, S-1, halfstep, opt8_a[(k+1)%12]*eps);
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
      itgr->integrate[S-1](opt6_rho*eps, S-1, halfstep, opt6_theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](opt6_theta*eps, S-1, halfstep, oneminus2thetarho*eps);
      
      update_momenta(itgr->mnls_per_ts[S], opt6_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](oneminus2thetarho*eps, S-1, halfstep, oneminus2thetarho*eps);
      update_momenta_fg(itgr->mnls_per_ts[S],eps , oneminus2lambdanu, opt6_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](oneminus2thetarho*eps, S-1, halfstep, opt6_theta*eps);
      update_momenta(itgr->mnls_per_ts[S], opt6_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      
      itgr->integrate[S-1](opt6_theta*eps, S-1, halfstep, opt6_rho*eps);
      update_momenta_fg(itgr->mnls_per_ts[S], eps , opt6_nu, opt6_mu , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](opt6_rho*eps, S-1, halfstep, opt6_rho*eps);
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
      itgr->integrate[S-1](opt4_theta*eps, S-1, halfstep, oneminus2theta*eps);
      update_momenta(itgr->mnls_per_ts[S], opt4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      
      itgr->integrate[S-1](oneminus2theta*eps, S-1, halfstep, oneminus2theta*eps);
      update_momenta_fg(itgr->mnls_per_ts[S],eps , oneminus2lambda, opt4_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](oneminus2theta*eps, S-1, halfstep, oneminus2theta*eps);
      
      update_momenta(itgr->mnls_per_ts[S], opt4_lamb*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](opt4_theta*eps, S-1, halfstep, oneminus2theta*eps);
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
      itgr->integrate[S-1](itgr->lambda[S]*eps, S-1, halfstep,oneminus2lambda*eps);
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S], &itgr->hf);

      itgr->integrate[S-1](oneminus2lambda*eps, S-1, halfstep,itgr->lambda[S]*eps);
      update_momenta(itgr->mnls_per_ts[S], 0.5*eps, itgr->no_mnls_per_ts[S], &itgr->hf);

      itgr->integrate[S-1](itgr->lambda[S]*eps, S-1, halfstep,itgr->lambda[S]*eps);
    }
  }
}

/* For 2MNFG lamda MUST be equal to 1/6 */
void integrate_2mnfg(const double tau, const int S, const int halfstep, const double tau2) {
  int i,j=0;
  integrator * itgr = &Integrator;
  double eps,eps2, oneminus2lambda = (1.-2*lambda_2mnfg);

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
      update_momenta(itgr->mnls_per_ts[0], 2.*lambda_2mnfg, itgr->no_mnls_per_ts[0], &itgr->hf);
    }
    update_gauge(0.5*eps, &itgr->hf);
    update_momenta_fg(itgr->mnls_per_ts[0],eps, oneminus2lambda, fg_chi , itgr->no_mnls_per_ts[0], &itgr->hf);
    update_gauge(0.5*eps, &itgr->hf);
    if(halfstep != 1) {
      update_momenta(itgr->mnls_per_ts[0], lambda_2mnfg*(eps+eps2), itgr->no_mnls_per_ts[0], &itgr->hf);
    }
  }
  else {
    for(i = 1; i < itgr->n_int[S]; i++){
      itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
      update_momenta_fg(itgr->mnls_per_ts[S],eps, oneminus2lambda, fg_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
      itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
      update_momenta(itgr->mnls_per_ts[S], 2*lambda_2mnfg*eps, itgr->no_mnls_per_ts[S], &itgr->hf);
    }
    itgr->integrate[S-1](eps/2., S-1, 0, eps/2);
    update_momenta_fg(itgr->mnls_per_ts[S],eps, oneminus2lambda, fg_chi , itgr->no_mnls_per_ts[S], &itgr->hf);
    if(S == itgr->no_timescales-1) {
      itgr->integrate[S-1](eps/2., S-1, 1, eps/2.);
    }
    else itgr->integrate[S-1](eps/2., S-1, halfstep, eps2/2.);
    if(halfstep != 1 && S != itgr->no_timescales-1) {
      update_momenta(itgr->mnls_per_ts[S], lambda_2mnfg*(eps+eps2), itgr->no_mnls_per_ts[S], &itgr->hf);
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
    else if(itgr->type[i] == OMFP5) {
     // update_momenta(itgr->mnls_per_ts[i], omf6_nu*eps, itgr->no_mnls_per_ts[i], &itgr->hf);
      update_momenta_fg(itgr->mnls_per_ts[i], eps , itgr->omfp5_nu[i], itgr->omfp5_mu[i] , itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/itgr->omfp5_theta[i-1]; /*CHECK THIS*/
    }
    else if(itgr->type[i] == OMF8FG) {
      update_momenta(itgr->mnls_per_ts[i], omf8_b[0]*eps, itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/omf8_a[1];
    }
    else if(itgr->type[i] == OMF4FG) {
      update_momenta_fg(itgr->mnls_per_ts[i],eps, omf4fg_b[0], omf4fg_c[0], itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/omf4fg_a[1];
    }
    else if(itgr->type[i] == OMF4FG2) {
      update_momenta_fg(itgr->mnls_per_ts[i],eps, omf4fg2_b[0], omf4fg2_c[0], itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/omf4fg2_a[1];
    }
    else if(itgr->type[i] == OMF6FG2) {
      update_momenta_fg(itgr->mnls_per_ts[i],eps, omf6fg2_b[0], omf6fg2_c[0], itgr->no_mnls_per_ts[i], &itgr->hf);
      eps /= ((double)itgr->n_int[i-1])/omf6fg2_a[1];
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
  else if (itgr->type[0] == OMFP5) {
   //  update_momenta(itgr->mnls_per_ts[0], omf6_nu*(eps), itgr->no_mnls_per_ts[0], &itgr->hf);
     update_momenta_fg(itgr->mnls_per_ts[0], eps , itgr->omfp5_nu[0], itgr->omfp5_mu[0] , itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if (itgr->type[0] == OMF8FG) {
     update_momenta(itgr->mnls_per_ts[0], omf8_b[0]*(eps), itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if (itgr->type[0] == OMF4FG) {
     update_momenta_fg(itgr->mnls_per_ts[0],eps, omf4fg_b[0],omf4fg_c[0], itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if (itgr->type[0] == OMF4FG2) {
     update_momenta_fg(itgr->mnls_per_ts[0],eps, omf4fg2_b[0], omf4fg2_c[0], itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  else if (itgr->type[0] == OMF6FG2) {
     update_momenta_fg(itgr->mnls_per_ts[0],eps, omf6fg2_b[0],omf6fg2_c[0], itgr->no_mnls_per_ts[0], &itgr->hf);
  }
  return;
}

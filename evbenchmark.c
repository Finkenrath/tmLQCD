/***********************************************************************
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008 Carsten Urbach
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
/*******************************************************************************
*
* Benchmark program for the even-odd preconditioned Wilson-Dirac operator
*
*
*******************************************************************************/

#ifdef HAVE_CONFIG_H
# include<config.h>
#endif
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#if (defined BGL && !defined BGP)
#  include <rts.h>
#endif
#ifdef TM_USE_MPI
# include <mpi.h>
# ifdef HAVE_LIBLEMON
#  include <io/params.h>
#  include <io/gauge.h>
# endif
#endif
#ifdef TM_USE_OMP
# include <omp.h>
# include "init/init_openmp.h"
#endif
#include "gettime.h"
#include "su3.h"
#include "su3adj.h"
#include "ranlxd.h"
#include "geometry_eo.h"
#include "read_input.h"
#include "start.h"
#include "boundary.h"
#include "operator/Hopping_Matrix.h"
#include "operator/Hopping_Matrix_nocom.h"
#include "operator/tm_operators.h"
#include "global.h"
#include "xchange/xchange.h"
#include "init/init.h"
#include "test/check_geometry.h"
#include "operator/D_psi.h"
#include "phmc.h"
#include "mpi_init.h"
#include "primme_interface.h"
#include "measure_gauge_action.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_32.h"
#include "operator/tm_operators_nd.h"
#include "operator/tm_operators_nd_32.h"
#include "expo.h"

#ifdef PARALLELT
#  define SLICE (LX*LY*LZ/2)
#elif defined PARALLELXT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2))
#elif defined PARALLELXYT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ/2))
#elif defined PARALLELXYZT
#  define SLICE ((LX*LY*LZ/2)+(T*LY*LZ/2) + (T*LX*LZ/2) + (T*LX*LY/2))
#elif defined PARALLELX
#  define SLICE ((LY*LZ*T/2))
#elif defined PARALLELXY
#  define SLICE ((LY*LZ*T/2) + (LX*LZ*T/2))
#elif defined PARALLELXYZ
#  define SLICE ((LY*LZ*T/2) + (LX*LZ*T/2) + (LX*LY*T/2))
#endif

int check_xchange();

int main(int argc,char *argv[])
{
  int i,j,j_max,k,k_max;
#ifdef HAVE_LIBLEMON
  paramsXlfInfo *xlfInfo;
#endif
  int status = 0,id,op;
  double plaquette_energy;
  static double t1,t2,dt,sdt,dts,qdt,sqdt;
  double tn,tb;
  double antioptaway=0.0;
  int compute_evs,no_eigenvalues,even_odd_flag,nstore;
  double eigenvalue_precision;
  char conf_filename[50];
  su3 u,u1,u2,u3;
  su3adj p;
  
#ifdef TM_USE_MPI
  static double dt2;
  
  DUM_DERI = 6;
  DUM_MATRIX = DUM_DERI+8;
  NO_OF_SPINORFIELDS = DUM_MATRIX+2;

#  ifdef TM_USE_OMP
  int mpi_thread_provided;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &mpi_thread_provided);
#  else
  MPI_Init(&argc, &argv);
#  endif
  MPI_Comm_rank(MPI_COMM_WORLD, &g_proc_id);

#else
  g_proc_id = 0;
#endif

  g_rgi_C1 = 1.; 

    /* Read the input file */
  if((status = read_input("benchmark.input")) != 0) {
    fprintf(stderr, "Could not find input file: benchmark.input\nAborting...\n");
    exit(-1);
  }
  k_max = DUM_MATRIX + 18;
#ifdef TM_USE_OMP
  init_openmp();
#endif

  tmlqcd_mpi_init(argc, argv);


  
  if(g_proc_id==0) {
#ifdef SSE
    printf("# The code was compiled with SSE instructions\n");
#endif
#ifdef SSE2
    printf("# The code was compiled with SSE2 instructions\n");
#endif
#ifdef SSE3
    printf("# The code was compiled with SSE3 instructions\n");
#endif
#ifdef P4
    printf("# The code was compiled for Pentium4\n");
#endif
#ifdef OPTERON
    printf("# The code was compiled for AMD Opteron\n");
#endif
#ifdef _GAUGE_COPY
    printf("# The code was compiled with -D_GAUGE_COPY\n");
#endif
#ifdef BGL
    printf("# The code was compiled for Blue Gene/L\n");
#endif
#ifdef BGP
    printf("# The code was compiled for Blue Gene/P\n");
#endif
#ifdef _USE_HALFSPINOR
    printf("# The code was compiled with -D_USE_HALFSPINOR\n");
#endif    
#ifdef _USE_SHMEM
    printf("# The code was compiled with -D_USE_SHMEM\n");
#  ifdef _PERSISTENT
    printf("# The code was compiled for persistent MPI calls (halfspinor only)\n");
#  endif
#endif
#ifdef TM_USE_MPI
#  ifdef _NON_BLOCKING
    printf("# The code was compiled for non-blocking MPI calls (spinor and gauge)\n");
#  endif
#endif
    printf("\n");
    fflush(stdout);
  }
  
  
#ifdef _GAUGE_COPY
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 1);
#else
  init_gauge_field(VOLUMEPLUSRAND + g_dbw2rand, 0);
#endif
  init_geometry_indices(VOLUMEPLUSRAND + g_dbw2rand);

  
  
	/*sprintf(conf_filename, "conf.save");
	if( (i = read_gauge_field(conf_filename,g_gauge_field)) !=0) 
	{
		fprintf(stderr, "Error %d while reading gauge field from %s\n Aborting...\n", i, conf_filename);
      exit(-2);
	}*/
  
  
  if(even_odd_flag) {
    j = init_spinor_field(VOLUMEPLUSRAND/2, k_max);
  }
  else {
    j = init_spinor_field(VOLUMEPLUSRAND, 2*k_max);
  }

  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for spinor fields! Aborting...\n");
    exit(0);
  }
  j = init_moment_field(VOLUME, VOLUMEPLUSRAND + g_dbw2rand);
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for moment fields! Aborting...\n");
    exit(0);
  }
  
  if(g_proc_id == 0) {
    fprintf(stdout,"# The number of processes is %d \n",g_nproc);
    printf("# The lattice size is %d x %d x %d x %d\n",
	   (int)(T*g_nproc_t), (int)(LX*g_nproc_x), (int)(LY*g_nproc_y), (int)(g_nproc_z*LZ));
    printf("# The local lattice size is %d x %d x %d x %d\n", 
	   (int)(T), (int)(LX), (int)(LY),(int) LZ);
    if(even_odd_flag) {
      printf("# benchmarking the even/odd preconditioned Dirac operator\n");
    }
    else {
      printf("# benchmarking the standard Dirac operator\n");
    }
    fflush(stdout);
  }
  
  /* define the geometry */
  geometry();
  /* define the boundary conditions for the fermion fields */
  boundary(g_kappa);

#ifdef _USE_HALFSPINOR
  j = init_dirac_halfspinor();
  if ( j!= 0) {
    fprintf(stderr, "Not enough memory for halfspinor fields! Aborting...\n");
    exit(0);
  }
  if(g_sloppy_precision_flag == 1) {
    g_sloppy_precision = 1;
    j = init_dirac_halfspinor32();
    if ( j!= 0) {
      fprintf(stderr, "Not enough memory for 32-Bit halfspinor fields! Aborting...\n");
      exit(0);
    }
  }
#  if (defined _PERSISTENT)
  init_xchange_halffield();
#  endif
#endif  

  status = check_geometry();
  if (status != 0) {
    fprintf(stderr, "Checking of geometry failed. Unable to proceed.\nAborting....\n");
    exit(1);
  }
#if (defined TM_USE_MPI && !(defined _USE_SHMEM))
  check_xchange(); 
#endif

  start_ranlux(1, 56);
  random_gauge_field(reproduce_randomnumber_flag, g_gauge_field);

#ifdef TM_USE_MPI
  /*For parallelization: exchange the gaugefield */
  xchange_gauge(g_gauge_field);
#endif
  plaquette_energy = measure_plaquette( (const su3**) g_gauge_field);
  
  
  op=1;
  phmc_invmaxev = 1.;
  compute_evs=1;
  no_eigenvalues=1;
  eigenvalue_precision=1e-8;
  even_odd_flag=0;
  nstore=0;
  
  
  u.c00=1.0+I*1.0;
  u.c01=1.0+I*1.0;
  u.c02=1.0+I*1.0;
  u.c10=1.0+I*1.0;
  u.c11=1.0+I*1.0;
  u.c12=1.0+I*1.0;
  u.c20=1.0+I*1.0;
  u.c21=1.0+I*1.0;
  u.c22=1.0+I*1.0;
  
  restoresu3(&u1, &u);
  restoresu3o(&u2,&u);
  
  printf("DIF u1 -v3xv1 \n c01 %.4e+1i%.4e\n",creal(conj(u1.c21 * u1.c02 - u1.c22 * u1.c01)-u1.c10),cimag(conj(u1.c21 * u1.c02 - u1.c22 * u1.c01)-u1.c10));
  printf("c11 %.4e+1i%.4e\n",creal(conj(u1.c22 * u1.c00 - u1.c20 * u1.c02)-u1.c11),cimag(conj(u1.c22 * u1.c00 - u1.c20 * u1.c02)-u1.c11));
  printf("c21 %.4e+1i%.4e\n",creal(conj(u1.c20 * u1.c01 - u1.c21 * u1.c00)-u1.c12),cimag(conj(u1.c20 * u1.c01 - u1.c21 * u1.c00)-u1.c12));
  
  printf("DIF u2 -v3xv1 \n c01 %.4e+1i%.4e\n",creal(conj(u2.c21 * u2.c02 - u2.c22 * u2.c01)-u2.c10),cimag(conj(u2.c21 * u2.c02 - u2.c22 * u2.c01)-u2.c10));
  printf("c11 %.4e+1i%.4e\n",creal(conj(u2.c22 * u2.c00 - u2.c20 * u2.c02)-u2.c11),cimag(conj(u2.c22 * u2.c00 - u2.c20 * u2.c02)-u2.c11));
  printf("c21 %.4e+1i%.4e\n",creal(conj(u2.c20 * u2.c01 - u2.c21 * u2.c00)-u2.c12),cimag(conj(u2.c20 * u2.c01 - u2.c21 * u2.c00)-u2.c12));
  
  printf("\n DIF \n c00 %.4e+1i%4.3e \n",creal(u1.c00-u2.c00),cimag(u1.c00-u2.c00));
  printf(" c01 %.4e+1i%4.3e \n",creal(u1.c01-u2.c01),cimag(u1.c01-u2.c01));
  printf(" c02 %.4e+1i%4.3e \n",creal(u1.c02-u2.c02),cimag(u1.c02-u2.c02));
  printf(" c10 %.4e+1i%4.3e \n",creal(u1.c10-u2.c10),cimag(u1.c10-u2.c10));
  printf(" c11 %.4e+1i%4.3e \n",creal(u1.c11-u2.c11),cimag(u1.c11-u2.c11));
  printf(" c12 %.4e+1i%4.3e \n",creal(u1.c12-u2.c12),cimag(u1.c12-u2.c12));
  printf(" c20 %.4e+1i%4.3e \n",creal(u1.c20-u2.c20),cimag(u1.c20-u2.c20));
  printf(" c21 %.4e+1i%4.3e \n",creal(u1.c21-u2.c21),cimag(u1.c21-u2.c21));
  printf(" c22 %.4e+1i%4.3e \n",creal(u1.c22-u2.c22),cimag(u1.c22-u2.c22));
  
  p.d1=0.232;
  p.d2=0.25;
  p.d3=-0.42;
  p.d4=0.92;
  p.d5=0.121;
  p.d6=-0.232;
  p.d7=0.612;
  p.d8=-0.322;
  
  exposu3o(&u1,&p);
  exposu3(&u2,&p);
  exposu3_check(&u3,&p,40);
  
  printf("\n DIF \n c00 %.4e+1i%4.3e \n",creal(u1.c00-u2.c00),cimag(u1.c00-u2.c00));
  printf(" c01 %.4e+1i%4.3e \n",creal(u1.c01-u2.c01),cimag(u1.c01-u2.c01));
  printf(" c02 %.4e+1i%4.3e \n",creal(u1.c02-u2.c02),cimag(u1.c02-u2.c02));
  printf(" c10 %.4e+1i%4.3e \n",creal(u1.c10-u2.c10),cimag(u1.c10-u2.c10));
  printf(" c11 %.4e+1i%4.3e \n",creal(u1.c11-u2.c11),cimag(u1.c11-u2.c11));
  printf(" c12 %.4e+1i%4.3e \n",creal(u1.c12-u2.c12),cimag(u1.c12-u2.c12));
  printf(" c20 %.4e+1i%4.3e \n",creal(u1.c20-u2.c20),cimag(u1.c20-u2.c20));
  printf(" c21 %.4e+1i%4.3e \n",creal(u1.c21-u2.c21),cimag(u1.c21-u2.c21));
  printf(" c22 %.4e+1i%4.3e \n",creal(u1.c22-u2.c22),cimag(u1.c22-u2.c22));
  
  printf("\n DIF \n c00 %.4e+1i%4.3e \n",creal(u3.c00-u2.c00),cimag(u3.c00-u2.c00));
  printf(" c01 %.4e+1i%4.3e \n",creal(u3.c01-u2.c01),cimag(u3.c01-u2.c01));
  printf(" c02 %.4e+1i%4.3e \n",creal(u3.c02-u2.c02),cimag(u3.c02-u2.c02));
  printf(" c10 %.4e+1i%4.3e \n",creal(u3.c10-u2.c10),cimag(u3.c10-u2.c10));
  printf(" c11 %.4e+1i%4.3e \n",creal(u3.c11-u2.c11),cimag(u3.c11-u2.c11));
  printf(" c12 %.4e+1i%4.3e \n",creal(u3.c12-u2.c12),cimag(u3.c12-u2.c12));
  printf(" c20 %.4e+1i%4.3e \n",creal(u3.c20-u2.c20),cimag(u3.c20-u2.c20));
  printf(" c21 %.4e+1i%4.3e \n",creal(u3.c21-u2.c21),cimag(u3.c21-u2.c21));
  printf(" c22 %.4e+1i%4.3e \n",creal(u3.c22-u2.c22),cimag(u3.c22-u2.c22));
  
  printf("\n DIF \n c00 %.4e+1i%4.3e \n",creal(u3.c00-u1.c00),cimag(u3.c00-u1.c00));
  printf(" c01 %.4e+1i%4.3e \n",creal(u3.c01-u1.c01),cimag(u3.c01-u1.c01));
  printf(" c02 %.4e+1i%4.3e \n",creal(u3.c02-u1.c02),cimag(u3.c02-u1.c02));
  printf(" c10 %.4e+1i%4.3e \n",creal(u3.c10-u1.c10),cimag(u3.c10-u1.c10));
  printf(" c11 %.4e+1i%4.3e \n",creal(u3.c11-u1.c11),cimag(u3.c11-u1.c11));
  printf(" c12 %.4e+1i%4.3e \n",creal(u3.c12-u1.c12),cimag(u3.c12-u1.c12));
  printf(" c20 %.4e+1i%4.3e \n",creal(u3.c20-u1.c20),cimag(u3.c20-u1.c20));
  printf(" c21 %.4e+1i%4.3e \n",creal(u3.c21-u1.c21),cimag(u3.c21-u1.c21));
  printf(" c22 %.4e+1i%4.3e \n",creal(u3.c22-u1.c22),cimag(u3.c22-u1.c22));
  
  
  
  restoresu3(&u3, &u2);
  restoresu3o(&u2,&u2);
  
  printf("\n DIF XXXXXRESTORE NEW \n c00 %.4e+1i%4.3e \n",creal(u3.c00-u2.c00),cimag(u3.c00-u2.c00));
  printf(" c01 %.4e+1i%4.3e \n",creal(u3.c01-u2.c01),cimag(u3.c01-u2.c01));
  printf(" c02 %.4e+1i%4.3e \n",creal(u3.c02-u2.c02),cimag(u3.c02-u2.c02));
  printf(" c10 %.4e+1i%4.3e \n",creal(u3.c10-u2.c10),cimag(u3.c10-u2.c10));
  printf(" c11 %.4e+1i%4.3e \n",creal(u3.c11-u2.c11),cimag(u3.c11-u2.c11));
  printf(" c12 %.4e+1i%4.3e \n",creal(u3.c12-u2.c12),cimag(u3.c12-u2.c12));
  printf(" c20 %.4e+1i%4.3e \n",creal(u3.c20-u2.c20),cimag(u3.c20-u2.c20));
  printf(" c21 %.4e+1i%4.3e \n",creal(u3.c21-u2.c21),cimag(u3.c21-u2.c21));
  printf(" c22 %.4e+1i%4.3e \n",creal(u3.c22-u2.c22),cimag(u3.c22-u2.c22));
  
  g_c_sw=1.0;
  g_kappa=0.125;
  g_mu=0.001;
  g_mubar=0.1;
  g_epsbar=0.01;
  init_sw_fields(VOLUME);
  sw_term((const su3 **)g_gauge_field, g_kappa, g_c_sw); 
  sw_invert_nd(g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  if (compute_evs != 0) 
  {
		 tb=MPI_Wtime();
		 if (op==0)
			eigenvalues_bi(&no_eigenvalues, 5000, eigenvalue_precision, 0, &Qsw_pm_ndbipsi);
		 else
			eigenvalues(&no_eigenvalues, 5000, eigenvalue_precision, 0,0,nstore,1);
		 
		 tn=MPI_Wtime();
		if (g_proc_id==0)
			printf("Time for EV %.2e sec\n",tn-tb);
  }
  sw_term((const su3 **)g_gauge_field, g_kappa, g_c_sw); 
  sw_invert_nd(g_mubar*g_mubar - g_epsbar*g_epsbar);
  
  primme_tm_init(1,1,0,1e-8);
  
  tb=MPI_Wtime();
  primme_tm_ev();
  tn=MPI_Wtime();
  if (g_proc_id==0)
  {
		printf("Time for EV %.2e sec\n",tn-tb);
		printf("Finalize WITH EV\n");
  }
  primme_tm_finalize();
  
  
  /*
  for (id=1;id<14;id++)
  {
	   primme_tm_init();
		primme_tm_set_meth(id);
		tb=MPI_Wtime();
		primme_tm_ev();
		tn=MPI_Wtime();
		if (g_proc_id==0){
			printf("Time for EV %.2e sec for PRIMME method %d\n",tn-tb,id);}
		primme_tm_finalize();
  }
  
  primme_tmsvds_init(0);
  primme_tm_svds();*/
  
#ifdef TM_USE_OMP
  free_omp_accumulators();
#endif
  free_gauge_field();
  free_geometry_indices();
  free_spinor_field();
  free_moment_field();
#ifdef TM_USE_MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#endif
  return(0);
}

/***********************************************************************
 *
 * Copyright (C) 2018 Jacob Finkenrath
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
 * Interface for primme
 * 
 * 
 * PRIMME stores the data on the structure primme_params, which has the next fields:
 * 
 * 			Basic
 * 		PRIMME_INT n, matrix dimension.
 * 		void (* matrixMatvec )(...), matrix-vector product.
 * 		int numEvals, how many eigenpairs to find.
 * 		primme_target target, which eigenvalues to find.
 * 		int numTargetShifts, for targeting interior eigenpairs.
 *  		double * targetShifts
 * 		double eps, tolerance of the residual norm of converged eigenpairs.
 *
 * 
 *		 		For parallel programs
 * 		int numProcs, number of processes
 * 		int procID, rank of this process
 * 		PRIMME_INT nLocal, number of rows stored in this process
 * 		void (* globalSumReal )(...), sum reduction among processes
 * 
 * 
 * 			Accelerate the convergence
 * 		void (* applyPreconditioner )(...), preconditioner-vector product.
 * 		int initSize, initial vectors as approximate solutions.
 * 		int maxBasisSize
 * 		int minRestartSize
 * 		int maxBlockSize
 * 
 * 
 * 			User data
 * 		void * commInfo
 * 		void * matrix
 * 		void * preconditioner
 * 		void * convTest
 * 		void * monitor
 * 
 * 
 * 			Advanced options
 * 		PRIMME_INT ldevecs, leading dimension of the evecs.
 * 		int numOrthoConst, orthogonal constrains to the eigenvectors.
 * 		int dynamicMethodSwitch
 * 		int locking
 * 		PRIMME_INT maxMatvecs
 *  		PRIMME_INT maxOuterIterations
 * 		int intWorkSize
 * 		size_t realWorkSize
 * 		PRIMME_INT iseed [4]
 * 		int * intWork
 * 		void * realWork
 * 		double aNorm
 * 		int printLevel
 * 		FILE * outputFile
 * 		double * ShiftsForPreconditioner
 * 		primme_init initBasisMode
 * 		struct projection_params projectionParams
 * 		struct restarting_params restartingParams
 * 		struct correction_params correctionParams
 * 		struct primme_stats stats
 * 		void (* convTestFun )(...), custom convergence criterion.
 * 		PRIMME_INT ldOPs, leading dimension to use in matrixMatvec.
 * 		void (* monitorFun )(...), custom convergence history.
 *
 *******************************************************************************/

#include "primme_interface.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "complex.h"
#include "boundary.h"
#include "gettime.h"
#include "read_input.h"
#include "primme.h"
#include "global.h"
#include "linalg/assign.h"
#include "operator/tm_operators.h"
#include "operator/tm_operators_32.h"
#include "operator/tm_operators_nd.h"
#include "operator/tm_operators_nd_32.h"
primme_params primme_tm;
primme_svds_params primme_svds_tm;

void par_GlobalSumForDouble(void *sendBuf, void *recvBuf, int *count, primme_params *primme, int *ierr) 
{
	MPI_Comm communicator = MPI_COMM_WORLD;	
   if (sendBuf == recvBuf)
     *ierr = MPI_Allreduce(MPI_IN_PLACE, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator) != MPI_SUCCESS;
   else
     *ierr = MPI_Allreduce(sendBuf, recvBuf, *count, MPI_DOUBLE, MPI_SUM, communicator) != MPI_SUCCESS;

}

void Q_squared(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
   
	spinor *l1,*k1;
	double d1,d2;
	
	MPI_Barrier(MPI_COMM_WORLD);
	Qsw_pm_psi((bispinor *) y, (bispinor *) x);
	*err = 0;
}

void Qnd_squared(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_params *primme, int *err) {
	
	MPI_Barrier(MPI_COMM_WORLD);
	Qsw_pm_ndbipsi((bispinor *) y, (bispinor *) x);
	/*
	 *    
	spinor *l1,*l2,*k1,*k2;
	double d1,d2;
	
	k1=g_spinor_field[0];
	k2=g_spinor_field[1];
	l1=g_spinor_field[2];
	l2=g_spinor_field[3];
	/printf("bx[0] %.8e %.8e g_proc_id %d\n",creal(k1[0].s0.c0),cimag(k1[0].s0.c0),g_proc_id);
	/*assign(k1, (double* ) x, VOLUME);
	assign(k2, (double* ) x+24*VOLUME, VOLUME);
   
	Qsw_pm_ndpsi(l1,l2,k1,k2);

	assign(y,l1, VOLUME);
	assign((double* ) y+24*VOLUME,l2, VOLUME);*/
	*err = 0;
}

void primme_tm_init(int op, int nev, int val, double prec)
{
 	primme_initialize(&primme_tm);
	primme_tm.globalSumReal=par_GlobalSumForDouble;
	primme_tm.numProcs=g_nproc;
	primme_tm.procID=g_proc_id;
	
	if (op==0)
	{
		primme_tm.matrixMatvec = Qnd_squared;
		primme_tm.n = 24*VOLUME*g_nproc/2;
		primme_tm.nLocal= 24*VOLUME/2;
	}
	else
	{
		primme_tm.matrixMatvec = Q_squared;
		primme_tm.n = 24*VOLUME*g_nproc/4;
		primme_tm.nLocal= 24*VOLUME/4;
	}
	
	/* Set problem parameters */
    /* set problem dimension */
   primme_tm.numEvals = nev;   /* Number of wanted eigenpairs */
   primme_tm.eps = prec;      /* ||r|| <= eps * ||matrix|| */
   if (val==0)
		primme_tm.target = primme_smallest;
	else
		primme_tm.target = primme_largest;
                           /* Wanted the smallest eigenvalues 
									 primme_smallest
									 primme_largest
									 */									
	primme_tm.commInfo=MPI_COMM_WORLD;
	primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme_tm);
}

void primme_tm_set_meth(int id)
{
	if (id==1)
		primme_set_method(PRIMME_DYNAMIC, &primme_tm);
	else if (id==2)
		primme_set_method(PRIMME_DEFAULT_MIN_TIME, &primme_tm);
	else if (id==3)
		primme_set_method(PRIMME_DEFAULT_MIN_MATVECS, &primme_tm);
	else if (id==4)
		primme_set_method(PRIMME_Arnoldi, &primme_tm);
	else if (id==5)
		primme_set_method(PRIMME_GD, &primme_tm);
	else if (id==6)
		primme_set_method(PRIMME_GD_plusK, &primme_tm);
	else if (id==7)
		primme_set_method(PRIMME_RQI, &primme_tm);
	else if (id==8)
		primme_set_method(PRIMME_JDQR, &primme_tm);
	else if (id==9)
		primme_set_method(PRIMME_JDQMR, &primme_tm);
	else if (id==10)
		primme_set_method(PRIMME_JDQMR_ETol, &primme_tm);
	else if (id==11)
		primme_set_method(PRIMME_STEEPEST_DESCENT, &primme_tm);
	else if (id==12)
		primme_set_method(PRIMME_LOBPCG_OrthoBasis_Window, &primme_tm);
	else if (id==13)
		primme_set_method(PRIMME_LOBPCG_OrthoBasis_Window, &primme_tm);
}


void primme_tm_reset(void) {
	
    primme_free(&primme_tm);
}

void primme_tm_finalize(void) {
    
	primme_free(&primme_tm);
}


int primme_tm_ev(double *evals) {
	int i;
	double *evecs,*rnorms;
	
	/* Allocate space for converged Ritz values and residual norms 
   evals = (double*)malloc(primme_tm.numEvals*sizeof(complex double));*/
   evecs = (complex double*)malloc(primme_tm.nLocal*(primme_tm.numEvals)*sizeof(complex double));
   rnorms = (double*)malloc(primme_tm.numEvals*sizeof(complex double));

   /* Call primme  */
   zprimme(evals, evecs, rnorms, &primme_tm);
	
	if (g_proc_id==0){
		for (i=0; i < primme_tm.initSize; i++) {
			fprintf(primme_tm.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
				evals[i], rnorms[i]); 
		}
	}
	
	/*free(evals);*/
	free(evecs);
	free(rnorms);
}

// /*
// void Dnd_augmented(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, primme_svds_params *primme_svds, int *err) {
//    //primme_svds_params *primme_svds = (primme_svds_params*)primme->matrix;
// 	complex double *x0 , *x1, *y0 , *y1;
// 	printf("AUG\n");
// 	
// 	
// 	x0 = (complex double*)x;
// 	x1 = &x0[primme_svds->nLocal];
// 	y0 = (complex double*)y;
// 	y1 = &y0[primme_svds->nLocal];
// 	
// 	printf("%d\n",primme_svds->nLocal );
// 	/*x0 = (complex double*)x;
// 	x1 = (complex double*)x;
// 	y0 = (complex double*)y;
// 	y1 = (complex double*)y;
// 	*/
// 	
//    /* [y0; y1] <-  * [x0; x1] */
// 	int notrans=0, trans=1;
//    /* y0 <- A^t * x1 */
//    Dnd_svd(x1, ldx, y0, ldy, blockSize, &trans, primme_svds, err);
//    /* y1 <- A * x0 */
//    Dnd_svd(x0, ldx, y1, ldy, blockSize, &notrans, primme_svds, err);			  
// 				  
// 				  
// 	/*			  
// 	spinor *l1,*l2,*l3,*l4,*k1,*k2,*k3,*k4;
// 	double *dx,*dy;
// 	int V;
// 	
// 	V=24*VOLUME/2;
// 	k1=g_spinor_field[0];
// 	k2=g_spinor_field[1];
// 	k3=g_spinor_field[2];
// 	k4=g_spinor_field[3];
// 	
// 	l1=g_spinor_field[4];
// 	l2=g_spinor_field[5];
// 	l3=g_spinor_field[6];
// 	l4=g_spinor_field[7];
// 	
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	
// 	dx=(double *) x;
// 	decompact(k1,k2,(bispinor *) x);
// 	decompact(k3,k4,(bispinor *) (dx+V));
// 	
// 	/*if (*transpose)
// 	{
// 		Qsw_ndpsi(l1,l2,k3,k4);
// 		Qsw_dagger_ndpsi(l3,l4,k1,k2);
// 	}
// 	else
// 	{
// 		Qsw_ndpsi(l3,l4,k1,k2);
// 		Qsw_dagger_ndpsi(l1,l2,k3,k4);
// 	//}
// 	
// 	dy=(double *) y;
// 	compact((bispinor *) y,l1,l2);
// 	compact((bispinor *) (dy+V),l3,l4);
// 	
// 	*err = 0;*/
// 	
// }
// 
// 
// void Dnd_svd(void *x, PRIMME_INT *ldx, void *y, PRIMME_INT *ldy, int *blockSize, int * transpose,  primme_svds_params *primme, int *err) {
//    
// 	spinor *l1,*l2,*k1,*k2;
// 	
// 	k1=g_spinor_field[0];
// 	k2=g_spinor_field[1];
// 	
// 	l1=g_spinor_field[0];
// 	l2=g_spinor_field[1];
// 	
// 	printf("huhu %d\n",*transpose);
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	decompact(k1,k2,(bispinor *) x);
// 	
// 	if (*transpose)
// 		Qsw_dagger_ndpsi(l1,l2,k1,k2);
// 	else
// 		Qsw_ndpsi(l1,l2,k1,k2);
// 
// 	MPI_Barrier(MPI_COMM_WORLD);
// 	printf("end %d\n",*transpose);
// 	
// 	compact((bispinor *) y,l1,l2);
// 	*err = 0;
// }
// 
// 
// void primme_tmsvds_init(int aug)
// {
// 	primme_svds_initialize (&primme_svds_tm);
// 	//primme_svds_tm.matrixMatvec = Dnd_svd;
// 	
// 	if (aug)
// 	{
// 		primme_svds_tm.matrixMatvec = Dnd_svd;
// 		primme_svds_tm.primme.matrixMatvec = Dnd_augmented; 
// 	}
// 	else
// 		primme_svds_tm.matrixMatvec = Dnd_svd;
// 	
// 	
// 	primme_svds_tm.globalSumReal=par_GlobalSumForDouble;
// 	primme_svds_tm.n = 24*VOLUME*g_nproc/2;
// 	primme_svds_tm.m = 24*VOLUME*g_nproc/2;
// 	primme_svds_tm.numProcs=g_nproc;
// 	primme_svds_tm.procID=g_proc_id;
// 	primme_svds_tm.nLocal= 24*VOLUME/2;
// 	primme_svds_tm.mLocal= 24*VOLUME/2;
// 	/* Set problem parameters */
//     /* set problem dimension */
//    primme_svds_tm.numSvals = 1;   /* Number of wanted evs */
//    primme_svds_tm.eps = 1e-5;      /* ||r|| <= eps * ||matrix|| */
//    primme_svds_tm.target = primme_svds_smallest;
//                            /* Wanted the smallest eigenvalues 
// 									 primme_svds_smallest
// 									 primme_svds_largest
// 									 primme_svds_closest_abs
// 									 */
// 	primme_svds_tm.printLevel=3;
// 						
// 	primme_svds_tm.commInfo=MPI_COMM_WORLD;
// 	/* Set method to solve the problem 
//    primme_svds_set_method (primme_svds_augmented, PRIMME_DEFAULT_METHOD,
// 									 primme_svds_op_none, &primme_svds_tm);*/
// 	if (aug)
// 		primme_svds_set_method(primme_svds_augmented, PRIMME_DEFAULT_METHOD,
//                           PRIMME_DEFAULT_METHOD, &primme_svds_tm);
// }
// 
// void primme_tmsvds_set_meth(int id)
// {
// 
// 	primme_svds_tm.matrixMatvec = Dnd_augmented; /* MV product */
// 	primme_svds_tm.n = 24*VOLUME*g_nproc/2;
// 	primme_svds_tm.m = 24*VOLUME*g_nproc/2;
// 	primme_svds_tm.numSvals = 2; /* Number of singular values */
// 	primme_svds_set_method(primme_svds_augmented, PRIMME_DEFAULT_MIN_MATVECS,
//                           PRIMME_DEFAULT_METHOD, &primme_svds_tm);
// 	
// 	/*
// 	primme_svds_display_params(&primme_svds_tm);
// 	*/
// }
// 
// void primme_tmsvds_finalize(){
// 	
// 	primme_free(&primme_svds_tm);
// }
// 
// int primme_tm_svds(void ) {
// 	int i;
// 	double * evals,*evecs,*rnorms;
// 	
// 	primme_svds_display_params(primme_svds_tm);
// 	/* Allocate space for converged Ritz values and residual norms */
//    evals = (double*)malloc(primme_svds_tm.numSvals*sizeof(complex double));
//    evecs = (complex double*)malloc(16*primme_svds_tm.nLocal*(primme_svds_tm.numSvals)*sizeof(complex double));
//    rnorms = (double*)malloc(primme_svds_tm.numSvals*sizeof(complex double));
//   
// 	printf("start\n");
// 	
//    /* Call primme  */
//    zprimme_svds(evals, evecs, rnorms, &primme_svds_tm); 
// 	for (i=0; i < primme_svds_tm.initSize; i++) {
//       fprintf(primme_svds_tm.outputFile, "Eval[%d]: %-22.15E rnorm: %-22.15E\n", i+1,
//          evals[i], rnorms[i]); 
// 	}
// 
// 	free(evals);
// 	free(evecs);
// 	free(rnorms);
// }
// //#endif*/

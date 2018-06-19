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

int primme_tm_ev(double *evals1) {
	int i;
	double *evals,*rnorms;
        complex double *evecs;
	
	/* Allocate space for converged Ritz values and residual norms */
	evals = (double*)malloc(primme_tm.numEvals*sizeof(complex double));
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
	(*evals1)=creal(*evals);
	
	free(evals);
	free(evecs);
	free(rnorms);
}
// 
// General-purpose Photovoltaic Device Model gpvdm.com - a drift diffusion
// base/Shockley-Read-Hall model for 1st, 2nd and 3rd generation solarcells.
// The model can simulate OLEDs, Perovskite cells, and OFETs.
// 
// Copyright (C) 2008-2020 Roderick C. I. MacKenzie
// 
// https://www.gpvdm.com
// r.c.i.mackenzie at googlemail.com
// 
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of the GPVDM nor the
//       names of its contributors may be used to endorse or promote products
//       derived from this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL Roderick C. I. MacKenzie BE LIABLE FOR ANY
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

/** @file lib.c
	@brief Call the SuperLU solver.
*/


#include <stdio.h>
#include <stdlib.h>

#include "slu_ddefs.h"

#include <util.h>

void triplet_to_ccs(double *out_x,int *out_row, int *out_pointer, int *in_i,int *in_j, double *in_x,int dim,int nz)
{
	int i;
	int j=0;
	out_pointer[0]=0;
	int *add=malloc(sizeof(int)*(dim+1));
	for (i=0;i<dim+1;i++)
	{
		out_pointer[i]=0;
		add[i]=0;
	}


	for (i=0;i<nz;i++)
	{
		j=in_j[i];
		out_pointer[j+1]=out_pointer[j+1]+1;
	}

	//for (i=0;i<dim+1;i++)
	//{
	//	printf("a=%d\n",out_pointer[i]);
	//}

	int sum=0;
	for (i=0;i<dim+1;i++)
	{
		sum=sum+out_pointer[i];
		add[i]=sum;
	}


	for (i=0;i<dim+1;i++)
	{
		out_pointer[i]=add[i];
		//printf("a=%d\n",out_pointer[i]);
	}


	//getchar();
//
	for (i=0;i<nz;i++)
	{
		j=in_j[i];
		out_x[add[j]]=in_x[i];
		out_row[add[j]]=in_i[i];
		add[j]++;
	}

	free(add);

}

int umfpack_solver(struct simulation *sim,int col,int nz,int *Ti,int *Tj, long double *lTx,long double *lb)
{
int i;
void *Symbolic, *Numeric;
int status;
double *dtemp;
int *itemp;


if ((sim->last_col!=col)||(sim->last_nz!=nz))
{
	dtemp = realloc(sim->x,col*sizeof(double));
	if (dtemp==NULL)
	{
		ewe(sim,"realloc failed\n");
	}else
	{
		sim->x=dtemp;
	}


	dtemp = realloc(sim->b,col*sizeof(double));
	if (dtemp==NULL)
	{
		ewe(sim,"realloc failed\n");
	}else
	{
		sim->b=dtemp;
	}

	itemp = realloc(sim->Ap,(col+1)*sizeof(int));
	if (itemp==NULL)
	{
		ewe(sim,"realloc failed\n");
	}else
	{
		sim->Ap=itemp;
	}

	itemp = realloc(sim->Ai,(nz)*sizeof(int));
	if (itemp==NULL)
	{
		ewe(sim,"realloc failed\n");
	}else
	{
		sim->Ai=itemp;
	}

	dtemp  = realloc(sim->Ax,(nz)*sizeof(double));
	if (dtemp==NULL)
	{
		ewe(sim,"realloc failed\n");
	}else
	{
		sim->Ax=dtemp;
	}

	dtemp  = realloc(sim->Tx,(nz)*sizeof(double));
	if (dtemp==NULL)
	{
		ewe(sim,"realloc failed\n");
	}else
	{
		sim->Tx=dtemp;
	}


	sim->last_col=col;
	sim->last_nz=nz;
}

for (i=0;i<col;i++)
{
	sim->b[i]=(double)lb[i];
}

for (i=0;i<nz;i++)
{
	sim->Tx[i]=(double)lTx[i];
}


    SuperMatrix A, L, U, B;
    double   *a, *rhs;
    double   s, u, p, e, r, l;
    int      *asub, *xa;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      nrhs, info, m, n, nnz, permc_spec;
    superlu_options_t options;
    SuperLUStat_t stat;

    /* Initialize matrix A. */
    m = n = col;
    if ( !(a = doubleMalloc(nz)) )
	{
		ABORT("Malloc fails for a[].");
	}

    if ( !(asub = intMalloc(nz)) )
	{
		ABORT("Malloc fails for asub[].");
	}

    if ( !(xa = intMalloc(col+1)) )
	{
		ABORT("Malloc fails for xa[].");
	}

	triplet_to_ccs(a,asub, xa, Ti,Tj, sim->Tx,col,nz);

    /* Create matrix A in the format expected by SuperLU. */
    dCreate_CompCol_Matrix(&A, col, col, nz, a, asub, xa, SLU_NC, SLU_D, SLU_GE);
    
    /* Create right-hand side matrix B. */

    dCreate_Dense_Matrix(&B, col, 1, sim->b, col, SLU_DN, SLU_D, SLU_GE);

    if ( !(perm_r = intMalloc(col)) )
	{
		ABORT("Malloc fails for perm_r[].");
	}

    if ( !(perm_c = intMalloc(col)) )
	{
		ABORT("Malloc fails for perm_c[].");
	}
    /* Set the default input options. */
    set_default_options(&options);
    options.ColPerm = NATURAL;
	options.IterRefine = SLU_DOUBLE;
	//options.PivotGrowth=YES;
	//options.DiagPivotThresh=2.0;
	//options.ConditionNumber=YES;

    /* Initialize the statistics variables. */
    StatInit(&stat);

    /* Solve the linear system. */
    dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);

//printf(">>>>>>.%d\n",nz);
//getchar();    
   // dPrint_CompCol_Matrix("A", &A);


   // dPrint_CompCol_Matrix("U", &U);
//printf("1 %d\n",nz);
//getchar();
   // dPrint_SuperNode_Matrix("L", &L);
//printf("2>>>>>>>>>>>>>>>> %d\n",nz);
//getchar();
   // print_int_vec("\nperm_r", m, perm_r);

	for (i=0;i<col;i++)
	{
		lb[i]=(long double)sim->b[i];
	//	printf("%Le\n",lb[i]);
	}

    /* De-allocate storage */


//printf("ok %d\n",nz);
//getchar();
    //SUPERLU_FREE (rhs);
//	printf("boom\n");
//	getchar();
    SUPERLU_FREE (perm_r);
//	printf("boom\n");
//	getchar();   
 SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
//	printf("here>>\n");
//	getchar();
    Destroy_SuperMatrix_Store(&B);
//	printf("boom\n");
//	getchar();
    Destroy_SuperNode_Matrix(&L);
//	printf("boom\n");
//	getchar();
    Destroy_CompCol_Matrix(&U);
//	printf("boom\n");
//	getchar();
    StatFree(&stat);
//printf("finished\n");
//getchar();


//memcpy(b, x, col*sizeof(double));
//umfpack_toc(stats);


return 0;
}

void umfpack_solver_free(struct simulation *sim)
{

if (sim->x!=NULL)
{
	free(sim->x);
	sim->x=NULL;	
}

if (sim->b!=NULL)
{
	free(sim->b);
	sim->b=NULL;	
}

if (sim->Ap!=NULL)
{
	free(sim->Ap);
	sim->Ap=NULL;
}

if (sim->Ai!=NULL)
{
	free(sim->Ai);
	sim->Ai=NULL;
}

if (sim->Ax!=NULL)
{
	free(sim->Ax);
	sim->Ax=NULL;
}

if (sim->Tx!=NULL)
{
	free(sim->Tx);
	sim->Tx=NULL;
}

sim->last_col=0;
sim->last_nz=0;
}


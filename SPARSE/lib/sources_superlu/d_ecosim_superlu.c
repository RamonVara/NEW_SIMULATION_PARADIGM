/* =========================================================================
* Wrapper to superlu solver for EcosimPro Language:
*
* op_flag          = int specifies the operation:
*                      1, performs LU decomposition for the first time
*                      2, performs triangular solve
*                      3, free all the storage in the end
* n                = dimension of matrix
* nnz              = # of nonzero entries
* nrhs             = # of RHSs
* values           = double array containing the nonzero entries 
* colind           = int array containing the column indices of the entries 
* rowptr           = int array containing the row start
* b                = double array containing the rhs vector (gets overwritten
*                    with solution)
* f_factors        = pointer to LU factors. (If op_flag == 1, it is an output 
*                    and contains the pointer pointing to the structure of 
*                    the factored matrices. Otherwise, it it an input.
* info             = info flag from superlu
*                 = 0: successful exit   
*                 < 0: if info = -i, the i-th argument had an illegal value   
*                 > 0: if info = i, and i is   
*                 <= A->ncol: U(i,i) is exactly zero. The factorization has   
*                       been completed, but the factor U is exactly   
*                       singular, so the solution and error bounds   
*                       could not be computed.   
*                 = A->ncol+1: U is nonsingular, but RCOND is less than machine
*                       precision, meaning that the matrix is singular to
*                       working precision. Nevertheless, the solution and
*                       error bounds are computed because there are a number
*                       of situations where the computed solution can be more
*                       accurate than the value of RCOND would suggest.   
*                 > A->ncol+1: number of bytes allocated when memory allocation
*                       failure occurred, plus A->ncol.
*
* debug           = 1/0 for debug&performance info/no debug&performance info 
* =========================================================================
*/
#include "slu_ddefs.h"

#define HANDLE_SIZE  8

/* kind of integer to hold a pointer.  Use 64-bit. */
typedef long long int fptr;

typedef struct {
	 SuperMatrix *A;
    SuperMatrix *L;
    SuperMatrix *U;
    int *perm_r;
    int *perm_c;
} factors_t;

void
d_ecosim_superlu(int op_flag, int n, int nnz, int nrhs, 
                 double *values, int *colind, int *rowptr,
                 double *b, 
		           fptr *f_factors, /* a handle containing the address
				     pointing to the factored matrices */
		           int *info, double pivot_threshold,  int equil, int refine, int debug)
{
/* 
 * This routine can be called from EcosimPro Language
 *
 * op_flag (input) int
 *      Specifies the operation:
 *      = 1, performs LU decomposition for the first time
 *      = 2, performs triangular solve
 *      = 3, free all the storage in the end
 *
 * f_factors (input/output) fptr* 
 *      If op_flag == 1, it is an output and contains the pointer pointing to
 *                    the structure of the factored matrices.
 *      Otherwise, it it an input.
 *
 */
   SuperMatrix A, AC, B;
   SuperMatrix *L, *U;
   int *perm_r; /* row permutations from partial pivoting */
	int *perm_c; /* column permutation vector */
   int *etree;  /* column elimination tree */
   SCformat *Lstore;
   NCformat *Ustore;
   int  i, notran, panel_size, permc_spec, relax;
   trans_t  trant;
   mem_usage_t   mem_usage;
   superlu_options_t options;
   SuperLUStat_t stat;
   factors_t *LUfactors;
   GlobalLU_t Glu;   /* Not needed on return. */
   int *colind0;  /*counter 1-based indexing from EcosimPro arrays*/
	int *rowptr0;  
	
   trant = TRANS;
	notran = 0;

	if ( op_flag == 1 ) { /* LU decomposition */ 
		/* Set the default input options. */
		set_default_options(&options);
		/*options.Trans = TRANS;*/
		options.DiagPivotThresh =  pivot_threshold;		  
		/* Initialize the statistics variables. */
		StatInit(&stat);

		/* Adjust to 0-based indexing */
		if ( !(colind0 = intMalloc(nnz)) ) ABORT("Malloc fails for colind0[].");
		if ( !(rowptr0 = intMalloc(n+1)) ) ABORT("Malloc fails for rowptr0[].");
		for (i = 0; i < nnz; ++i) colind0[i] = colind[i] - 1;
		for (i = 0; i <= n; ++i) rowptr0[i] = rowptr[i] - 1;
 
		dCreate_CompCol_Matrix(&A, n, n, nnz, values, colind0, rowptr0,
			       SLU_NC, SLU_D, SLU_GE);	 

		L = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
		U = (SuperMatrix *) SUPERLU_MALLOC( sizeof(SuperMatrix) );
	
		if ( !(perm_r = intMalloc(n)) ) ABORT("Malloc fails for perm_r[].");
		if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
		if ( !(etree = intMalloc(n)) ) ABORT("Malloc fails for etree[].");
 	
		/*
		* Get column permutation vector perm_c[], according to permc_spec:
		*   permc_spec = 0: natural ordering 
		*   permc_spec = 1: minimum degree on structure of A'*A
		*   permc_spec = 2: minimum degree on structure of A'+A
		*   permc_spec = 3: approximate minimum degree for unsymmetric matrices
		*/   
		permc_spec = options.ColPerm;
		get_perm_c(permc_spec, &A, perm_c);
		
		sp_preorder(&options, &A, perm_c, etree, &AC);
	
		panel_size = sp_ienv(1);
		relax = sp_ienv(2);
		dgstrf(&options, &AC, relax, panel_size, etree,
                NULL, 0, perm_c, perm_r, L, U, &Glu, &stat, info);

		if(*info == 0) {
			Lstore = (SCformat *) L->Store;
			Ustore = (NCformat *) U->Store;
			if(debug!=0) {
				printf("\n");
				printf("SUPERLU returns INFO= %d\n", *info);
				printf("No of nonzeros in factor L = %d\n", Lstore->nnz);
				printf("No of nonzeros in factor U = %d\n", Ustore->nnz);
				printf("No of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz - n);
				printf("Fill ratio = %.1f\n", (float)(Lstore->nnz + Ustore->nnz - n)/nnz);
				printf("\n");
				dQuerySpace(L, U, &mem_usage);
				printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
					mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
			}
		} else {
			printf("\n");
			printf("d_ecosim_superlu returns INFO= %d\n", *info);
			if ( *info <= n ) { /* factorization completes */
				dQuerySpace(L, U, &mem_usage);
				printf("L\\U MB %.3f\ttotal MB needed %.3f\n",
					mem_usage.for_lu/1e6, mem_usage.total_needed/1e6);
			}
		}
		/* Save the LU factors and addittional addresses in the factors handle */
		LUfactors = (factors_t*) SUPERLU_MALLOC(sizeof(factors_t));
		LUfactors->L = L;
		LUfactors->U = U;
		LUfactors->perm_c = perm_c;
		LUfactors->perm_r = perm_r;
		*f_factors = (fptr) LUfactors;

		/* Free un-wanted storage */
		SUPERLU_FREE(etree);
		Destroy_SuperMatrix_Store(&A);
		Destroy_CompCol_Permuted(&AC);
		SUPERLU_FREE(rowptr0);
		SUPERLU_FREE(colind0);
		/*StatPrint(&stat);*/
		StatFree(&stat);
		return;
	} else if ( op_flag == 2 ) { /* Triangular solve */
		/* Initialize the statistics variables. */
		StatInit(&stat);

		/* Extract the LU factors in the factors handle */
		LUfactors = (factors_t*) *f_factors;
		L = LUfactors->L;
		U = LUfactors->U;
		perm_c = LUfactors->perm_c;
		perm_r = LUfactors->perm_r;

		dCreate_Dense_Matrix(&B, n, nrhs, b, n, SLU_DN, SLU_D, SLU_GE);
		
		/* Solve the system A*X=B, overwriting B with X. */
		dgstrs (trant, L, U, perm_c, perm_r, &B, &stat, info);
		
		Destroy_SuperMatrix_Store(&B);
		/*StatPrint(&stat);*/
		StatFree(&stat);
		return;
   } else if ( op_flag == 3 ) { /* Free storage */
		/* Free the LU factors in the factors handle */
		LUfactors = (factors_t*) *f_factors;
		SUPERLU_FREE (LUfactors->perm_r);
		SUPERLU_FREE (LUfactors->perm_c);
		Destroy_SuperNode_Matrix(LUfactors->L);
		Destroy_CompCol_Matrix(LUfactors->U);
      SUPERLU_FREE (LUfactors->L);
      SUPERLU_FREE (LUfactors->U);
		SUPERLU_FREE (LUfactors);
		return;
	} else {
		fprintf(stderr,"Invalid op_flag=%d passed to d_ecosim_superlu()\n",op_flag);
		*info = -1;
		return;
	}
}
